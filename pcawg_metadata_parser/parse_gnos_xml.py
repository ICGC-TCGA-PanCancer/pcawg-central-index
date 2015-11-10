#!/usr/bin/env python

# Author: Junjun Zhang

import sys
import os
import re
import glob
import xmltodict
import json
import yaml
import copy
import logging
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from elasticsearch import Elasticsearch
from collections import OrderedDict
import datetime
import dateutil.parser
from itertools import izip
from distutils.version import LooseVersion
import csv
import hashlib

logger = logging.getLogger('gnos parser')
# create console handler with a higher log level
ch = logging.StreamHandler()


def init_es(es_host, es_index):
    es = Elasticsearch([ es_host ])
    #es.indices.delete( es_index, ignore=[400, 404] )
    es.indices.create( es_index, ignore=400 )

    # create mappings
    es_mapping = open('pancan.donor.mapping.json')
    es.indices.put_mapping(index=es_index, doc_type='donor', body=es_mapping.read())
    es_mapping.close()

    es_mapping = open('pancan.file.mapping.json')
    es.indices.put_mapping(index=es_index, doc_type='bam_file', body=es_mapping.read())
    es_mapping.close()

    return es


def process_gnos_analysis(gnos_analysis, donors, vcf_entries, es_index, es, bam_output_fh, annotations):
  analysis_attrib = get_analysis_attrib(gnos_analysis)

  if analysis_attrib and analysis_attrib.get('variant_workflow_name'):  # variant call gnos entry
    donor_unique_id = analysis_attrib.get('dcc_project_code') + '::' + analysis_attrib.get('submitter_donor_id')

    if is_in_donor_blacklist(donor_unique_id):
        logger.warning('ignore blacklisted donor: {} GNOS entry: {}'
                         .format(donor_unique_id, gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    if gnos_analysis.get('study').lower().endswith('_test'):
        logger.warning('ignore variant calling entry with study ending with _test, donor: {} GNOS entry: {}'
                         .format(donor_unique_id, gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    #logger.info('Create variant calling file for donor: {}, from entry {}'
    #        .format(donor_unique_id, gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull')))

    vcf_file = create_vcf_entry(donor_unique_id, analysis_attrib, gnos_analysis, annotations)

    if not vcf_entries.get(donor_unique_id):
        vcf_entries[donor_unique_id] = {}
        vcf_entries[donor_unique_id]['vcf_entry_files'] = []
    
    vcf_entries.get(donor_unique_id)['vcf_entry_files'].append(copy.deepcopy(vcf_file))

  else:  # BAM entry
    if gnos_analysis.get('dcc_project_code') and gnos_analysis.get('dcc_project_code').upper() == 'TEST':
        logger.warning('ignore entry with dcc_project_code being TEST, GNOS entry: {}'
                         .format(gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    # if gnos_analysis.get('library_strategy') and gnos_analysis.get('library_strategy') == 'RNA-Seq':
    #     logger.warning('ignore entry with library_strategy being RNA-Seq for now, GNOS entry: {}'
    #                      .format(gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
    #     return

    if not gnos_analysis.get('aliquot_id'):
        logger.warning('ignore entry does not have aliquot_id, GNOS entry: {}'
                         .format(gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    if gnos_analysis.get('refassem_short_name') != 'unaligned' and gnos_analysis.get('refassem_short_name') != 'GRCh37':
        logger.warning('ignore entry that is aligned but not aligned to GRCh37: {}'
                         .format(gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return # completely ignore test gnos entries for now, this is the quickest way to avoid test interferes real data 

    if not analysis_attrib:
        logger.warning('ignore entry does not have ANALYSIS information, GNOS entry: {}'
                         .format(gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    if not analysis_attrib.get('dcc_project_code') or not analysis_attrib.get('submitter_donor_id') \
            or '/' in analysis_attrib.get('submitter_donor_id'):
        logger.warning('ignore entry does not have dcc_project_code or submitter_donor_id, or submitter_donor_id is invalid, GNOS entry: {}'
                         .format(gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    if not analysis_attrib.get('submitter_specimen_id') or not analysis_attrib.get('submitter_sample_id'):
        logger.warning('ignore entry does not have submitter_specimen_id or submitter_sample_id, GNOS entry: {}'
                         .format(gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    # added on Apr. 24, 2015 after discovering that one RNA-Seq uploaded to GNOS with TCGA barcode which was treated as a new donor
    if analysis_attrib.get('dcc_project_code').endswith('-US') and \
            analysis_attrib.get('submitter_donor_id').startswith('TCGA-'):
        logger.warning('ignore TCGA entry submitted with barcode, GNOS entry: {}'
                         .format(gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    donor_unique_id = analysis_attrib.get('dcc_project_code') + '::' + analysis_attrib.get('submitter_donor_id')

    if is_in_donor_blacklist(donor_unique_id):
        logger.warning('ignore blacklisted donor: {} GNOS entry: {}'
                         .format(donor_unique_id, gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    if is_test(analysis_attrib, gnos_analysis):
        logger.warning('ignore test entry: {}'.format(gnos_analysis.get('analysis_detail_uri')))
        return # completely ignore test gnos entries for now, this is the quickest way to avoid test interferes real data 

    if gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS'].get('TITLE') and gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS']['TITLE'].startswith('TCGA/ICGC PanCancer Specimen-Level Germline Variant Calling for Specimen'):
        logger.warning('ignore Annai germline call entry: {}'.format(gnos_analysis.get('analysis_detail_uri')))
        return

    if gnos_analysis.get('library_strategy') == 'RNA-Seq' and not analysis_attrib.get('workflow_name') in ('RNA-Seq_Alignment_SOP_STAR', 'Workflow_Bundle_TOPHAT2'):
        logger.warning('ignore RNA-Seq entry that is not STAR or TOPHAT2 aligned, entry: {}'.format(gnos_analysis.get('analysis_detail_uri')))
        return

    if (gnos_analysis.get('library_strategy') == 'WGS' and gnos_analysis.get('refassem_short_name') != 'unaligned'
              and not is_train_2_aligned(analysis_attrib, gnos_analysis)
            ):
        # TODO: we may create another ES index for obsoleted BAM entries
        # TODO: we will need more sophisticated check for handling BAMs that are flagged as aligned but
        #       treated as unaligned (this is actually the case for TCGA input BAM entries, maybe need a full
        #       TCGA spciment list from Marc?)
        logger.warning('ignore entry that is aligned but not by train 2 protocol: {}'
                             .format( gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return    

    if gnos_analysis.get('library_strategy') == 'WGS' and is_corrupted_train_2_alignment(analysis_attrib, gnos_analysis):
        logger.warning('ignore entry that is aligned by train 2 protocol but seems corrupted: {}'
                             .format( gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    # temporary hack here to skip any VALIDATION entries
    if gnos_analysis.get('library_strategy') == 'VALIDATION':
        logger.warning('ignore entry that is VALIDATION: {}'
                             .format( gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    #TODO: put things above into one function

    # temporary hack here to skip any BAM entries from odsc-tcga repo as it's supposed not contain
    # any BAM data, but it does, and those aligned BAMs it has overlap with what in CGHub hence causes problems
    if 'osdc-tcga' in gnos_analysis.get('analysis_detail_uri'):
        logger.warning('ignore BAM entry in osdc-tcga repo: {}'
                         .format( gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    if not donors.get(donor_unique_id):
        # create a new donor if not exist
        donors[ donor_unique_id ] = create_donor(donor_unique_id, analysis_attrib, gnos_analysis, annotations)

    else: # the donor this bam entry belongs to already exists
        # perform some comparison between existing donor and the info in the current bam entry
        if (donors[donor_unique_id].get('gnos_study') != gnos_analysis.get('study')):
            logger.warning( 'existing donor {} has study {}, but study in current gnos ao is {}'.
                            format( donor_unique_id,
                                    donors[donor_unique_id].get('gnos_study'),
                                    gnos_analysis.get('study') ) )
        # more such check may be added, no time for this now

    #logger.info('Create bam file for donor: {}, from entry {}'
    #        .format(donor_unique_id, gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull')))

    # now parse out gnos analysis object info to build bam_file doc
    bam_file = create_bam_file_entry(donor_unique_id, analysis_attrib, gnos_analysis, annotations)
     

    # only do the following when it is WGS
    if bam_file.get('library_strategy') == 'WGS':
        if 'normal' in bam_file.get('dcc_specimen_type').lower(): # normal
            if donors.get(donor_unique_id).get('normal_specimen'): # normal specimen exists
                if donors.get(donor_unique_id).get('normal_specimen').get('aliquot_id') == gnos_analysis.get('aliquot_id'):
                    if bam_file.get('is_aligned'):
                        if donors.get(donor_unique_id)['normal_specimen'].get('is_aligned'):
                            logger.info('more than one normal aligned bam for donor: {}, entry in use: {}, additional entry found in: {}'
                                  .format(donor_unique_id,
                                      donors.get(donor_unique_id).get('normal_specimen').get('gnos_metadata_url'),
                                      gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull')
                                  )
                            )
                            if (not donors.get(donor_unique_id).get('normal_specimen').get('gnos_metadata_url').split('/')[-1]
                                    == gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull').split('/')[-1]):
                                logger.warning('Two aligned BAM entries for the same normal specimen from donor: {} have different GNOS UUIDs: {} and {}'
                                    .format(donor_unique_id,
                                        donors.get(donor_unique_id).get('normal_specimen').get('gnos_metadata_url'),
                                        gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull')
                                    )
                                )
                            if donors.get(donor_unique_id).get('normal_specimen').get('upload_date') < bam_file.get(
                                    'upload_date'):  # the current one is newer
                                donors.get(donor_unique_id)['normal_specimen'].update(
                                    prepare_aggregated_specimen_level_info(copy.deepcopy(bam_file))
                                )
                                donors.get(donor_unique_id)['gnos_repo'] = bam_file.get('gnos_repo')
                        else:
                            donors.get(donor_unique_id)['normal_specimen'].update(
                                prepare_aggregated_specimen_level_info(copy.deepcopy(bam_file))
                            )
                            donors.get(donor_unique_id)['gnos_repo'] = bam_file.get('gnos_repo')
                else:
                    logger.warning('same donor: {} has different aliquot_id: {}, {} for normal specimen, entry in use: {}, additional entry found in {}'
                      .format(donor_unique_id,
                          donors.get(donor_unique_id).get('normal_specimen').get('aliquot_id'),
                          gnos_analysis.get('aliquot_id'),
                          donors.get(donor_unique_id).get('normal_specimen').get('gnos_metadata_url'),
                          gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull')
                      )
                    )                
            else:
                # add normal_specimen
                donors.get(donor_unique_id)['normal_specimen'].update(
                    prepare_aggregated_specimen_level_info(copy.deepcopy(bam_file))
                )
                # update donor's 'gnos_repo' field with normal aligned specimen
                donors.get(donor_unique_id)['gnos_repo'] = bam_file.get('gnos_repo')

        else: # not normal
            donors.get(donor_unique_id).get('all_tumor_specimen_aliquots').add(bam_file.get('aliquot_id'))
            donors.get(donor_unique_id).get('flags')['all_tumor_specimen_aliquot_counts'] = len(donors.get(donor_unique_id).get('all_tumor_specimen_aliquots'))
            if bam_file.get('is_aligned'):
                if donors.get(donor_unique_id).get('aligned_tumor_specimens'):
                    if donors.get(donor_unique_id).get('aligned_tumor_specimen_aliquots').intersection(
                            [ bam_file.get('aliquot_id') ]
                        ): # multiple alignments for the same tumor aliquot_id
                        logger.warning('more than one tumor aligned bam for donor: {} with aliquot_id: {}, additional entry found in: {}'
                              .format(donor_unique_id,
                                  bam_file.get('aliquot_id'),
                                  gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull')
                            )
                        )
                    else:
                        donors.get(donor_unique_id).get('aligned_tumor_specimens').append( copy.deepcopy(bam_file) )
                        donors.get(donor_unique_id).get('aligned_tumor_specimen_aliquots').add(bam_file.get('aliquot_id'))
                        donors.get(donor_unique_id).get('flags')['aligned_tumor_specimen_aliquot_counts'] = len(donors.get(donor_unique_id).get('aligned_tumor_specimen_aliquots'))
                else:  # create the first element of the list
                    donors.get(donor_unique_id)['aligned_tumor_specimens'] = [copy.deepcopy(bam_file)]
                    donors.get(donor_unique_id).get('aligned_tumor_specimen_aliquots').add(bam_file.get('aliquot_id'))  # set of aliquot_id
                    donors.get(donor_unique_id).get('flags')['aligned_tumor_specimen_aliquot_counts'] = 1
                    donors.get(donor_unique_id).get('flags')['has_aligned_tumor_specimen'] = True

    original_gnos = bam_file['gnos_repo']
    bam_file.update( donors[ donor_unique_id ] )
    bam_file['gnos_repo'] = original_gnos
    del bam_file['bam_files']
    del bam_file['normal_specimen']
    del bam_file['aligned_tumor_specimens']
    del bam_file['aligned_tumor_specimen_aliquots']
    del bam_file['all_tumor_specimen_aliquots']
    del bam_file['flags']
    del bam_file['rna_seq']

    
    donors[donor_unique_id]['bam_files'].append( copy.deepcopy(bam_file) )       

    # push to Elasticsearch
    # Let's not worry about this index type, it seems not that useful
    #es.index(index=es_index, doc_type='bam_file', id=bam_file['bam_gnos_ao_id'], body=json.loads( json.dumps(bam_file, default=set_default) ), timeout=90)
    bam_output_fh.write(json.dumps(bam_file, default=set_default) + '\n')


def choose_vcf_entry(vcf_entries, donor_unique_id, annotations):

    if not vcf_entries or not vcf_entries.get(donor_unique_id) or not vcf_entries.get(donor_unique_id).get('vcf_entry_files'):
        return
    
    for current_vcf_entry in vcf_entries.get(donor_unique_id).get('vcf_entry_files'):
        variant_workflow = current_vcf_entry.get('vcf_workflow_type')
        workflow_label = variant_workflow + '_variant_calling'

        if not vcf_entries.get(donor_unique_id).get(workflow_label):  # new vcf for workflow_type
            vcf_entries.get(donor_unique_id).update({workflow_label: current_vcf_entry})
        else:
            workflow_previous = vcf_entries.get(donor_unique_id).get(workflow_label)
            if workflow_previous.get('gnos_id') == current_vcf_entry.get('gnos_id'):
                if current_vcf_entry['gnos_repo'][0] in workflow_previous.get('gnos_repo'):
                    logger.warning( 'Same donor: {} has multiple variant calling with same GNOS ID: {} in the same repo: {}. This should never be possible.'
                                    .format(donor_unique_id, workflow_previous.get('gnos_id'), current_vcf_entry['gnos_repo'][0])) 
                
                else:
                    workflow_previous.get('gnos_repo').append(current_vcf_entry['gnos_repo'][0])
                    workflow_previous.get('gnos_last_modified').append(current_vcf_entry['gnos_last_modified'][0])
                    workflow_previous.get('effective_xml_md5sum').append(current_vcf_entry['effective_xml_md5sum'][0])
                    workflow_previous['exists_xml_md5sum_mismatch'] = False if len(set(workflow_previous.get('effective_xml_md5sum'))) == 1 else True
                    logger.info( 'Donor: {} has synchronized variant calling with GNOS ID: {} in the repos: {}'
                                    .format(donor_unique_id, workflow_previous.get('gnos_id'), '|'.join(workflow_previous.get('gnos_repo')))) 
            
            else:
                if current_vcf_entry['is_oct2015_entry']:
                    vcf_entries.get(donor_unique_id).update({workflow_label: current_vcf_entry})
                    logger.info(workflow_label+' results for donor: {}. Keep the oct2015_freeze_entry: {}, additional {}'
                        .format(donor_unique_id, current_vcf_entry['gnos_id'], workflow_previous['gnos_id']))
                elif current_vcf_entry['is_s3_transfer_scheduled']:
                    vcf_entries.get(donor_unique_id).update({workflow_label: current_vcf_entry})
                    logger.info(workflow_label+' results for donor: {}. Keep the one scheduled for S3 transfer: {}, additional {}'
                        .format(donor_unique_id, current_vcf_entry['gnos_id'], workflow_previous['gnos_id']))
                
                elif current_vcf_entry['is_aug2015_entry']:
                    vcf_entries.get(donor_unique_id).update({workflow_label: current_vcf_entry})
                    logger.info(workflow_label+' results for donor: {}. Keep the aug2015_freeze_entry: {}, additional {}'
                        .format(donor_unique_id, current_vcf_entry['gnos_id'], workflow_previous['gnos_id']))

                elif current_vcf_entry['is_santa_cruz_entry']:
                    vcf_entries.get(donor_unique_id).update({workflow_label: current_vcf_entry})
                    logger.info(workflow_label+' results for donor: {}. Keep the santa_cruz_freeze_entry: {}, additional {}'
                        .format(donor_unique_id, current_vcf_entry['gnos_id'], workflow_previous['gnos_id']))

                elif annotations.get(variant_workflow) and current_vcf_entry.get('gnos_id') in annotations.get(variant_workflow):
                    vcf_entries.get(donor_unique_id).update({workflow_label: current_vcf_entry})
                    logger.info(workflow_label+' results for donor: {}. Keep the one in whitelist: {}, additional {}'
                        .format(donor_unique_id, current_vcf_entry['gnos_id'], workflow_previous['gnos_id']))

                elif annotations.get(variant_workflow+'_vcf_in_jamboree') and \
                       annotations.get(variant_workflow+'_vcf_in_jamboree').get(donor_unique_id) and \
                         annotations.get(variant_workflow+'_vcf_in_jamboree').get(donor_unique_id) == current_vcf_entry.get('gnos_id'):
                    vcf_entries.get(donor_unique_id).update({workflow_label: current_vcf_entry})
                    logger.info(workflow_label+' results for donor: {}. Keep the one already saved in Jamboree: {}, additional {}'
                        .format(donor_unique_id, current_vcf_entry['gnos_id'], workflow_previous['gnos_id']))                         

                else:
                    workflow_version_current = current_vcf_entry.get('workflow_details').get('variant_workflow_version')
                    workflow_version_previous = workflow_previous.get('workflow_details').get('variant_workflow_version')
                    gnos_updated_current = current_vcf_entry.get('gnos_last_modified')[0]
                    gnos_updated_previous = workflow_previous.get('gnos_last_modified')[0]

                    if LooseVersion(workflow_version_current) > LooseVersion(workflow_version_previous): # current is newer version
                        logger.info('Newer {} variant calling result with version: {} for donor: {}, with GNOS entry: {} in {} replacing older GNOS entry {} in {}'
                            .format(variant_workflow.upper(), workflow_version_current, donor_unique_id, \
                                current_vcf_entry.get('gnos_id'), current_vcf_entry.get('gnos_repo')[0],\
                                workflow_previous.get('gnos_id'), '|'.join(workflow_previous.get('gnos_repo'))))
                        vcf_entries.get(donor_unique_id)[workflow_label] = current_vcf_entry
                    elif LooseVersion(workflow_version_current) == LooseVersion(workflow_version_previous) \
                         and gnos_updated_current > gnos_updated_previous: # current is newer
                        logger.info('Newer {} variant calling result with last modified date: {} for donor: {}, with GNOS entry: {} in {} replacing older GNOS entry {} in {}'
                            .format(variant_workflow.upper(), gnos_updated_current, donor_unique_id, \
                                current_vcf_entry.get('gnos_id'), current_vcf_entry.get('gnos_repo')[0],\
                                workflow_previous.get('gnos_id'), '|'.join(workflow_previous.get('gnos_repo'))))
                        vcf_entries.get(donor_unique_id)[workflow_label] = current_vcf_entry
                    else: # no need to replace
                        logger.warning('{} variant calling result already exist and is latest for donor: {}, ignoring entry {} in {}'
                            .format(variant_workflow.upper(), donor_unique_id, current_vcf_entry.get('gnos_id'), current_vcf_entry.get('gnos_repo')[0]))

    


def create_vcf_entry(donor_unique_id, analysis_attrib, gnos_analysis, annotations):
    files = []
    
    if isinstance(gnos_analysis.get('files').get('file'), dict):
        file_list = [gnos_analysis.get('files').get('file')]  
    elif isinstance(gnos_analysis.get('files').get('file'), list):
        file_list = gnos_analysis.get('files').get('file') 
    else:
        logger.warning('Variant calling result donor: {}, likely incorrectly populated the files section, in GNOS entry {}'
            .format(donor_unique_id, gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull')))

    for f in file_list:
        files.append({'file_name': f.get('filename'), 'file_size': f.get('filesize'), 'file_md5sum': f.get('checksum').get('#text')})

    vcf_entry = {
        #'analysis_attrib': analysis_attrib, # remove this later
        #'gnos_analysis': gnos_analysis, # remove this later
        "gnos_id": gnos_analysis.get('analysis_id'),
        "gnos_repo": [gnos_analysis.get('analysis_detail_uri').split('/cghub/')[0] + '/'],
        "gnos_last_modified": [dateutil.parser.parse(gnos_analysis.get('last_modified'))],
        "files": files,
        "study": gnos_analysis.get('study'),
        "effective_xml_md5sum": [gnos_analysis.get('_effective_xml_md5sum')],
        "is_santa_cruz_entry": True if gnos_analysis.get('analysis_id') in annotations.get('santa_cruz').get('gnos_id') else False,
        "is_aug2015_entry": True if gnos_analysis.get('analysis_id') in annotations.get('aug2015').get('gnos_id') else False,
        "is_oct2015_entry": True if gnos_analysis.get('analysis_id') in annotations.get('oct2015').get('gnos_id') else False,
        "is_s3_transfer_scheduled": True if gnos_analysis.get('analysis_id') in annotations.get('s3_transfer_scheduled') else False,
        "is_s3_transfer_completed": True if gnos_analysis.get('analysis_id') in annotations.get('s3_transfer_completed') else False,
        "exists_xml_md5sum_mismatch": False,
        "variant_calling_performed_at": gnos_analysis.get('analysis_xml').get('ANALYSIS_SET').get('ANALYSIS').get('@center_name'),
        "workflow_details": {
            "variant_workflow_name": analysis_attrib.get('variant_workflow_name'),
            "variant_workflow_version": analysis_attrib.get('variant_workflow_version'),
            "variant_pipeline_input_info": json.loads( analysis_attrib.get('variant_pipeline_input_info') ).get('workflow_inputs') if analysis_attrib.get('variant_pipeline_input_info') else [],
            "variant_pipeline_output_info": json.loads( analysis_attrib.get('variant_pipeline_output_info') ).get('workflow_outputs') if analysis_attrib.get('variant_pipeline_output_info') else [],
            "variant_qc_metrics": {},
            "variant_timing_metrics": {}
        }
    }

    qc = {}
    try:
        qc = json.loads( analysis_attrib.get('variant_qc_metrics') ).get('qc_metrics')
    except:
        logger.warning('variant_qc_metrics format incorrect: {}'.format(analysis_attrib.get('variant_qc_metrics')))

    if isinstance(qc, dict): vcf_entry.get('workflow_details')['variant_qc_metrics'] = qc

    # DO NOT KEEP timing metrics, it's way too verbose
    #timing = json.loads( analysis_attrib.get('variant_timing_metrics') ).get('timing_metrics') if analysis_attrib.get('variant_timing_metrics') else {}
    #if isinstance(timing, dict): vcf_entry.get('workflow_details')['variant_timing_metrics'] = timing

    #print json.dumps(vcf_entry)  # debugging only
    workflow_name = vcf_entry.get('workflow_details').get('variant_workflow_name')
    workflow_version = vcf_entry.get('workflow_details').get('variant_workflow_version')

    if workflow_name == 'SangerPancancerCgpCnIndelSnvStr' and (( workflow_version.startswith('1.0.') or workflow_version.startswith('1.1.'))
            and not workflow_version in ['1.0.0', '1.0.1']):
        vcf_entry['vcf_workflow_type'] = 'sanger'

    elif workflow_name.startswith('EMBLPancancer') and LooseVersion(workflow_version) >= LooseVersion('1.0.0'):
        vcf_entry['vcf_workflow_type'] = 'embl'

    elif workflow_name == 'DKFZPancancerCnIndelSnv' and LooseVersion(workflow_version) >= LooseVersion('1.0.0'):
        vcf_entry['vcf_workflow_type'] = 'dkfz'      

    elif workflow_name == 'EMBLDKFZPancancerStrCnIndelSNV' and LooseVersion(workflow_version) >= LooseVersion('1.0.5'):
        vcf_entry['vcf_workflow_type'] = 'dkfz_embl'      

    elif workflow_name == 'DKFZ_EMBL_Combined_HPC':
        vcf_entry['vcf_workflow_type'] = 'dkfz_embl'

    elif workflow_name == 'DKFZ_EMBL_Merged':
        vcf_entry['vcf_workflow_type'] = 'dkfz_embl'

    elif workflow_name == 'BROAD_MUSE_PIPELINE':
        vcf_entry.get('workflow_details')['workflow_file_subset'] = analysis_attrib.get('workflow_file_subset')
        vcf_entry.get('workflow_details')['related_file_subset_uuids'] = analysis_attrib.get('related_file_subset_uuids').split(',')
        if vcf_entry.get('workflow_details').get('workflow_file_subset') == 'broad':
            vcf_entry['vcf_workflow_type'] = 'broad'
        elif vcf_entry.get('workflow_details').get('workflow_file_subset') == 'muse':
            vcf_entry['vcf_workflow_type'] = 'muse'
        elif vcf_entry.get('workflow_details').get('workflow_file_subset') == 'broad_tar':
            vcf_entry['vcf_workflow_type'] = 'broad_tar'   
        else:
            vcf_entry['vcf_workflow_type'] = 'Unknown_broad'
            logger.warning('broad variant calling entry which has unknown file type {}, donor: {} GNOS entry: {}'
                     .format(vcf_entry.get('workflow_details')['workflow_file_subset'], donor_unique_id, gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
            
    else:
        vcf_entry['vcf_workflow_type'] = 'Unknown'
        logger.warning('the entry is variant calling but likely is test entry, donor: {} GNOS entry: {}'
            .format(donor_unique_id, gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))


    return vcf_entry


def set_default(obj):
    if isinstance(obj, datetime.datetime):
        return obj.isoformat()
    if isinstance(obj, set):
        return list(obj)
    raise TypeError


def prepare_aggregated_specimen_level_info(bam_file):
    specimen = copy.deepcopy(bam_file)
    # TODO: actual aggregation to be completed
    return specimen


def is_in_donor_blacklist(donor_unique_id):
    donor_blacklist = set([
            "PACA-CA::PCSI_0449",
            "PACA-CA::PCSI_0309",
            "LIHC-US::G1551",
            "LIHC-US::G15512",
            "TCGA_MUT_BENCHMARK_4::G15511",
            "TCGA_MUT_BENCHMARK_4::G15512",
            "PBCA-DE::SNV_CALLING_TEST"
        ])
    if donor_blacklist.intersection([donor_unique_id]):
        return True
    else:
        return False


def create_bam_file_entry(donor_unique_id, analysis_attrib, gnos_analysis, annotations):
    file_info = parse_bam_file_info(gnos_analysis.get('files').get('file'))
    bam_file = {
        "dcc_specimen_type": analysis_attrib.get('dcc_specimen_type'),
        "submitter_specimen_id": analysis_attrib.get('submitter_specimen_id'),
        "submitter_sample_id": analysis_attrib.get('submitter_sample_id'),
        "icgc_specimen_id": get_icgc_id(donor_unique_id, analysis_attrib['dcc_project_code'], analysis_attrib['submitter_specimen_id'], 'specimen', annotations),
        "icgc_sample_id": get_icgc_id(donor_unique_id, analysis_attrib['dcc_project_code'], analysis_attrib['submitter_sample_id'], 'sample', annotations),                       
        "aliquot_id": gnos_analysis.get('aliquot_id'),
        "use_cntl": analysis_attrib.get('use_cntl'),
        "total_lanes": analysis_attrib.get('total_lanes'),

        "effective_xml_md5sum": gnos_analysis.get('_effective_xml_md5sum'),
        "is_santa_cruz_entry": True if gnos_analysis.get('analysis_id') in annotations.get('santa_cruz').get('gnos_id') else False,
        "is_aug2015_entry": True if gnos_analysis.get('analysis_id') in annotations.get('aug2015').get('gnos_id') else False,
        "is_oct2015_entry": True if gnos_analysis.get('analysis_id') in annotations.get('oct2015').get('gnos_id') else False,
        "is_s3_transfer_scheduled": True if gnos_analysis.get('analysis_id') in annotations.get('s3_transfer_scheduled') else False,
        "is_s3_transfer_completed": True if gnos_analysis.get('analysis_id') in annotations.get('s3_transfer_completed') else False,

        "library_strategy": gnos_analysis.get('library_strategy'),
        "gnos_repo": gnos_analysis.get('analysis_detail_uri').split('/cghub/')[0] + '/',
        "gnos_metadata_url": gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull'),
        "refassem_short_name": gnos_analysis.get('refassem_short_name'),
        "bam_gnos_ao_id": gnos_analysis.get('analysis_id'),
        "upload_date": dateutil.parser.parse(gnos_analysis.get('upload_date')),
        "published_date": dateutil.parser.parse(gnos_analysis.get('published_date')),
        "last_modified": dateutil.parser.parse(gnos_analysis.get('last_modified')),

        "bam_file_name": file_info.get('file_name'),
        "bam_file_size": file_info.get('file_size'),
        "md5sum": file_info.get('md5sum'),

        "bai_file_name": file_info.get('bai_file_name'),
        "bai_file_size": file_info.get('bai_file_size'),
        "bai_file_md5sum": file_info.get('bai_file_md5sum'),

    }

    # much more TODO for bam file info and alignment details
    if bam_file.get('refassem_short_name') == 'unaligned' and \
            gnos_analysis.get('library_strategy') == 'WGS' :
        bam_file['is_aligned'] = False
        bam_file['bam_type'] = 'Unaligned BAM'
        bam_file['alignment'] = None  # or initiate as empty object {}, depending on how ES searches it
    elif (analysis_attrib.get('workflow_output_bam_contents') == 'unaligned'
            or gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS']['DESCRIPTION'].startswith('The BAM file includes unmapped reads extracted from specimen-level BAM with the reference alignment')
         ) and gnos_analysis.get('library_strategy') == 'WGS' : # this is actually BAM with unmapped reads
        bam_file['is_aligned'] = False
        bam_file['bam_type'] = 'Specimen level unmapped reads after BWA alignment'
        bam_file['alignment'] = None
    elif gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS']['DESCRIPTION'].startswith('Specimen-level BAM from the reference alignment') \
             and gnos_analysis.get('library_strategy') == 'WGS' :
        bam_file['is_aligned'] = True
        bam_file['bam_type'] = 'Specimen level aligned BAM'
        bam_file['alignment'] = get_alignment_detail(analysis_attrib, gnos_analysis)
    elif (gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS']['DESCRIPTION'].lower().startswith('star ') \
            or gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS']['DESCRIPTION'].lower().startswith('tophat2 ')) \
            and gnos_analysis.get('library_strategy') == 'RNA-Seq' :
        bam_file['is_aligned'] = True
        bam_file['bam_type'] = 'RNA-Seq aligned BAM'
        bam_file['alignment'] = get_rna_seq_alignment_detail(analysis_attrib, gnos_analysis)
    elif (bam_file.get('refassem_short_name') == 'unaligned' and gnos_analysis.get('library_strategy') == 'RNA-Seq'):
        bam_file['is_aligned'] = False
        bam_file['bam_type'] = 'RNA-Seq unaligned BAM'
        bam_file['alignment'] = None

    else:
        bam_file['is_aligned'] = False
        bam_file['bam_type'] = 'Unknown'
        bam_file['alignment'] = None

    
    return bam_file


def get_rna_seq_alignment_detail(analysis_attrib, gnos_analysis):
    alignment = {
        "workflow_name": analysis_attrib.get('workflow_name'),
        "workflow_version": analysis_attrib.get('workflow_version'),
        "workflow_bundle_url": analysis_attrib.get('workflow_bundle_url'),
        "workflow_source_url": analysis_attrib.get('workflow_source_url')
    }

    return alignment


def get_alignment_detail(analysis_attrib, gnos_analysis):
    alignment = {
        "data_train": "Train 2",
        "workflow_name": analysis_attrib.get('workflow_name'),
        "workflow_version": analysis_attrib.get('workflow_version'),
        "workflow_bundle_url": analysis_attrib.get('workflow_bundle_url'),
        "workflow_source_url": analysis_attrib.get('workflow_source_url'),

        "pipeline_input_info": json.loads( analysis_attrib.get('pipeline_input_info') ).get('pipeline_input_info') if analysis_attrib.get('pipeline_input_info') else [],
        "qc_metrics": json.loads( analysis_attrib.get('qc_metrics').replace('"not_collected"', 'null') ).get('qc_metrics') if analysis_attrib.get('qc_metrics') else [],
        "markduplicates_metrics": json.loads( analysis_attrib.get('markduplicates_metrics') ).get('markduplicates_metrics') if analysis_attrib.get('markduplicates_metrics') else [],
        "timing_metrics": json.loads( analysis_attrib.get('timing_metrics').replace('"not_collected"', 'null') ).get('timing_metrics') if analysis_attrib.get('timing_metrics') else [],
    }

    alignment['input_bam_summary'] = {} # TODO: do this in a function
    
    return alignment


def parse_bam_file_info(file_fragment):
    file_info = {}
    if (type(file_fragment) != list): file_fragment = [file_fragment]

    for f in file_fragment:
        f = dict(f)
        if f.get('filename').endswith('.bam'): # assume there is only one BAM file
            file_info['file_name'] = f.get('filename')
            file_info['file_size'] = int(f.get('filesize'))
            file_info['md5sum'] = f.get('checksum').get('#text')
        elif f.get('filename').endswith('.bai'): # assume there is only one BAI file
            file_info['bai_file_name'] = f.get('filename')
            file_info['bai_file_size'] = int(f.get('filesize'))
            file_info['bai_file_md5sum'] = f.get('checksum').get('#text')

    return file_info


def is_train_2_aligned(analysis_attrib, gnos_analysis):
    if ( gnos_analysis.get('refassem_short_name') == 'GRCh37'
           and analysis_attrib.get('workflow_version')
           and analysis_attrib.get('workflow_version').startswith('2.6.')
       ):
        return True
    else:
        return False


def is_corrupted_train_2_alignment(analysis_attrib, gnos_analysis):
    if ( is_train_2_aligned(analysis_attrib, gnos_analysis)
           and not gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS']['DESCRIPTION'].startswith('The BAM file includes unmapped reads extracted from specimen-level BAM with the reference alignment')
           and (not analysis_attrib.get('qc_metrics') or not analysis_attrib.get('markduplicates_metrics'))
       ):
        return True
    else:
        return False


def detect_and_low_case_uuid(submitter_id):
    uuid_pattern = re.compile('^[a-f0-9]{8}-?[a-f0-9]{4}-?4[a-f0-9]{3}-?[89ab][a-f0-9]{3}-?[a-f0-9]{12}\Z', re.I)
    uuid = submitter_id.lower() if uuid_pattern.search(submitter_id) else submitter_id
    return uuid


def get_icgc_id(donor_unique_id, dcc_project_code, submitter_id, subtype, annotations):
    submitter_id = detect_and_low_case_uuid(submitter_id)
    if dcc_project_code.endswith('-US'):
        if not annotations.get('uuid_to_barcode').get(submitter_id):
            logger.warning('donor: {}, the {} with uuid: {} has no mapping barcode'.format(donor_unique_id, subtype, submitter_id))
            return None
        submitter_id = annotations.get('uuid_to_barcode').get(submitter_id)

    if not annotations.get('icgc_'+subtype+'_id').get(dcc_project_code+'::'+submitter_id):
        logger.warning('donor: {}, the {} with pcawg_id: {} has no mapping icgc_id'.format(donor_unique_id, subtype, submitter_id))
        return None
    icgc_id = annotations.get('icgc_'+subtype+'_id').get(dcc_project_code+'::'+submitter_id)
    return icgc_id


def create_donor(donor_unique_id, analysis_attrib, gnos_analysis, annotations):
    donor = {
        'donor_unique_id': donor_unique_id,
        'submitter_donor_id': analysis_attrib['submitter_donor_id'],
        'dcc_project_code': analysis_attrib['dcc_project_code'],
        'icgc_donor_id': get_icgc_id(donor_unique_id, analysis_attrib['dcc_project_code'], analysis_attrib['submitter_donor_id'], 'donor', annotations), 
        'gnos_study': gnos_analysis.get('study'),
        'gnos_repo': gnos_analysis.get('analysis_detail_uri').split('/cghub/')[0] + '/', # can be better
        'flags': {
            'is_test': is_test(analysis_attrib, gnos_analysis),
            'is_cell_line': is_cell_line(analysis_attrib, gnos_analysis),
            'is_train2_donor': False,
            'is_train2_pilot': False,
            'is_santa_cruz_donor': True if donor_unique_id in annotations.get('santa_cruz').get('donor') else False,
            'is_aug2015_donor': True if donor_unique_id in annotations.get('aug2015').get('donor') else False,
            'is_oct2015_donor': True if donor_unique_id in annotations.get('oct2015').get('donor') else False,
            'is_normal_specimen_aligned': False,
            'are_all_tumor_specimens_aligned': False,
            'has_aligned_tumor_specimen': False,
            'aligned_tumor_specimen_aliquot_counts': 0,
            'all_tumor_specimen_aliquot_counts': 0,
            'is_sanger_variant_calling_performed': False,
            'is_dkfz_variant_calling_performed': False,
            'is_embl_variant_calling_performed': False,
            'is_dkfz_embl_variant_calling_performed': False,
            'is_broad_variant_calling_performed': False,
            'broad':{
                'broad_file_subset_exist': False,
                'broad_tar_file_subset_exist': False,
                'muse_file_subset_exist': False,
                'exist_file_subsets_mismatch': False
            },
            'variant_calling_performed': [],
            'vcf_in_jamboree': [],
            'is_normal_star_rna_seq_alignment_performed': False,
            'is_normal_tophat_rna_seq_alignment_performed': False,
            'is_tumor_star_rna_seq_alignment_performed': False,
            'is_tumor_tophat_rna_seq_alignment_performed': False,
            'exists_vcf_file_prefix_mismatch': False,
            'is_bam_used_by_variant_calling_missing': False,
            'qc_score': None,
            'exists_xml_md5sum_mismatch': False
        },
        'normal_specimen': {},
        'aligned_tumor_specimens': [],
        'aligned_tumor_specimen_aliquots': set(),
        'all_tumor_specimen_aliquots': set(),
        'bam_files': [],
        'rna_seq': {
            'alignment': {
                'normal': {},
                'tumor': []
            }
        }
    }
    try:
        if type(gnos_analysis.get('experiment_xml').get('EXPERIMENT_SET').get('EXPERIMENT')) == list:
            donor['sequencing_center'] = gnos_analysis.get('experiment_xml').get('EXPERIMENT_SET').get('EXPERIMENT')[0].get('@center_name')
        else:
            donor['sequencing_center'] = gnos_analysis.get('experiment_xml').get('EXPERIMENT_SET').get('EXPERIMENT').get('@center_name')
    except:
        logger.warning('analysis object has no sequencing_center information: {}'.format(gnos_analysis.get('analysis_detail_uri')))

    
    if not annotations.get('qc_donor_prioritization'):
        logger.warning('Missing qc_donor_prioritization annotation')
    elif annotations.get('qc_donor_prioritization').get(donor_unique_id) is not None:
        donor.get('flags')['qc_score'] = annotations.get('qc_donor_prioritization').get(donor_unique_id)
    else:
        logger.warning('No qc prioritization score for donor: {}'.format(donor_unique_id))

    return donor


def is_test(analysis_attrib, gnos_analysis):
    if (gnos_analysis.get('aliquot_id') == '85098796-a2c1-11e3-a743-6c6c38d06053'
          or gnos_analysis.get('study') == 'CGTEST'
          or gnos_analysis.get('study') == 'icgc_pancancer_vcf_test'
          or gnos_analysis.get('study').lower().endswith('_test')
        ):
        return True
    elif (analysis_attrib.get('dcc_project_code') == 'None-US'
          and analysis_attrib.get('submitter_donor_id') == 'None'
          and analysis_attrib.get('submitter_specimen_id') == 'None'
          and analysis_attrib.get('dcc_specimen_type') == 'unknown'
        ):
        return True
    # TODO: what's the criteria for determining *test* entries

    return False


def is_cell_line(analysis_attrib, gnos_analysis):
    is_cell_line = False
    if analysis_attrib.get('dcc_project_code') == 'TCGA_MUT_BENCHMARK_4':
        is_cell_line = True

    return is_cell_line


def get_analysis_attrib(gnos_analysis):
    analysis_attrib = {}
    if (not gnos_analysis['analysis_xml']['ANALYSIS_SET'].get('ANALYSIS')
          or not gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS'].get('ANALYSIS_ATTRIBUTES')
          or not gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS']['ANALYSIS_ATTRIBUTES'].get('ANALYSIS_ATTRIBUTE')
       ):
        return None

    analysis_attrib_fragment = gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS']['ANALYSIS_ATTRIBUTES']['ANALYSIS_ATTRIBUTE']
    if (type(analysis_attrib_fragment) != list): analysis_attrib_fragment = [analysis_attrib_fragment]

    for a in analysis_attrib_fragment:
        if not analysis_attrib.get(a['TAG']):
            analysis_attrib[a['TAG']] = a['VALUE']
        else:
            logger.warning('duplicated analysis attribute key: {}'.format(a['TAG']))
    return analysis_attrib


def get_gnos_analysis(f):
    with open (f, 'r') as x: xml_str = x.read()
    gnos_analysis = xmltodict.parse(xml_str).get('ResultSet').get('Result')
    add_effective_xml_md5sum(gnos_analysis, xml_str)
    return gnos_analysis


def add_effective_xml_md5sum(gnos_analysis, xml_str):
    xml_str = re.sub(r'<ResultSet .+?>', '<ResultSet>', xml_str)
    xml_str = re.sub(r'<last_modified>.+?</last_modified>', '<last_modified></last_modified>', xml_str)
    xml_str = re.sub(r'<upload_date>.+?</upload_date>', '<upload_date></upload_date>', xml_str)
    xml_str = re.sub(r'<published_date>.+?</published_date>', '<published_date></published_date>', xml_str)
    xml_str = re.sub(r'<center_name>.+?</center_name>', '<center_name></center_name>', xml_str)
    xml_str = re.sub(r'<analyte_code>.+?</analyte_code>', '<analyte_code></analyte_code>', xml_str)
    xml_str = re.sub(r'<reason>.+?</reason>', '<reason></reason>', xml_str)
    xml_str = re.sub(r'<study>.+?</study>', '<study></study>', xml_str)
    xml_str = re.sub(r'<sample_accession>.+?</sample_accession>', '<sample_accession></sample_accession>', xml_str)

    #xml_str = re.sub(r'<dcc_project_code>.+?</dcc_project_code>', '<dcc_project_code></dcc_project_code>', xml_str)
    #xml_str = re.sub(r'<participant_id>.+?</participant_id>', '<participant_id></participant_id>', xml_str)
    xml_str = re.sub(r'<dcc_specimen_type>.+?</dcc_specimen_type>', '<dcc_specimen_type></dcc_specimen_type>', xml_str)

    xml_str = re.sub(r'<specimen_id>.+?</specimen_id>', '<specimen_id></specimen_id>', xml_str)
    xml_str = re.sub(r'<sample_id>.+?</sample_id>', '<sample_id></sample_id>', xml_str)
    xml_str = re.sub(r'<use_cntl>.+?</use_cntl>', '<use_cntl></use_cntl>', xml_str)
    xml_str = re.sub(r'<library_strategy>.+?</library_strategy>', '<library_strategy></library_strategy>', xml_str)
    xml_str = re.sub(r'<platform>.+?</platform>', '<platform></platform>', xml_str)
    xml_str = re.sub(r'<refassem_short_name>.+?</refassem_short_name>', '<refassem_short_name></refassem_short_name>', xml_str)

    xml_str = re.sub(r'<STUDY_REF .+?/>', '<STUDY_REF/>', xml_str)
    xml_str = re.sub(r'</STUDY_REF>', '', xml_str)
    xml_str = re.sub(r'<STUDY_REF .+?>', '<STUDY_REF/>', xml_str)
    xml_str = re.sub(r'<ANALYSIS_SET .+?>', '<ANALYSIS_SET>', xml_str)
    xml_str = re.sub(r'<ANALYSIS .+?>', '<ANALYSIS>', xml_str)
    xml_str = re.sub(r'<EXPERIMENT_SET .+?>', '<EXPERIMENT_SET>', xml_str)
    xml_str = re.sub(r'<RUN_SET .+?>', '<RUN_SET>', xml_str)
    xml_str = re.sub(r'<analysis_detail_uri>.+?</analysis_detail_uri>', '<analysis_detail_uri></analysis_detail_uri>', xml_str)
    xml_str = re.sub(r'<analysis_submission_uri>.+?</analysis_submission_uri>', '<analysis_submission_uri></analysis_submission_uri>', xml_str)
    xml_str = re.sub(r'<analysis_data_uri>.+?</analysis_data_uri>', '<analysis_data_uri></analysis_data_uri>', xml_str)

    # we need to take care of xml properties in different order but effectively/semantically the same
    effective_eq_xml = json.dumps(xmltodict.parse(xml_str).get('ResultSet').get('Result'), indent=4, sort_keys=True)
    # print effective_eq_xml
    # sys.exit()

    gnos_analysis.update({'_effective_xml_md5sum': hashlib.md5(effective_eq_xml).hexdigest()})

    return gnos_analysis


def get_xml_files( metadata_dir, conf, repo ):
    xml_files = []
    #ao_seen = {}
    for r in conf.get('gnos_repos'):
        if repo and not r.get('repo_code') == repo:
            continue
        gnos_ao_list_file = metadata_dir + '/analysis_objects.' + r.get('repo_code') + '.tsv'
        if not os.path.isfile(gnos_ao_list_file):
            logger.warning('gnos analsysi object list file does not exist: {}'.format(gnos_ao_list_file))
            continue
        with open(gnos_ao_list_file, 'r') as list:
            for ao in list:
                ao_uuid, ao_state = str.split(ao, '\t')[0:2]
                if not ao_state == 'live': continue  # skip ao that is not live
                #if (ao_seen.get(ao_uuid)): continue  # skip ao if already added
                #ao_seen[ao_uuid] = 1  # include this one
                xml_files.append(r.get('repo_code') + '/' + ao.replace('\t', '__').replace('\n','') + '.xml')

    return xml_files


def process(metadata_dir, conf, es_index, es, donor_output_jsonl_file, bam_output_jsonl_file, repo, exclude_gnos_id_lists):
    donors = {}
    vcf_entries = {}

    # update the pc_annotation-sanger_vcf_in_jamboree files using the jamboree subdirectory files
    vcf_in_jamboree_dir = '../pcawg-operations/variant_calling/sanger_workflow/jamboree/'
    
    infiles = glob.glob(vcf_in_jamboree_dir+'/Sanger_jamboree_batch*.txt')
    outfile = 'pc_annotation-sanger_vcf_in_jamboree.tsv' # hard-code file name
    update_vcf_jamboree(infiles, outfile)    

    annotations = {}
    read_annotations(annotations, 'gnos_assignment', 'pc_annotation-gnos_assignment.yml')  # hard-code file name for now
    read_annotations(annotations, 'train2_pilot', 'pc_annotation-train2_pilot.tsv')  # hard-code file name for now
    read_annotations(annotations, 'donor_blacklist', 'pc_annotation-donor_blacklist.tsv')  # hard-code file name for now
    read_annotations(annotations, 'manual_qc_failed', 'pc_annotation-manual_qc_failed.tsv')  # hard-code file name for now
    read_annotations(annotations, 'sanger_vcf_in_jamboree', 'pc_annotation-sanger_vcf_in_jamboree.tsv')  # hard-code file name for now
    read_annotations(annotations, 'santa_cruz', '../pcawg-operations/data_releases/santa_cruz/santa_cruz_freeze_entry.tsv')
    read_annotations(annotations, 'aug2015', '../pcawg-operations/data_releases/aug2015/release_aug2015_entry.tsv')
    read_annotations(annotations, 'oct2015', '../pcawg-operations/data_releases/oct2015/release_oct2015_entry.tsv')
    read_annotations(annotations, 's3_transfer_scheduled', '../s3-transfer-operations/s3-transfer-jobs*/*/*.json')
    read_annotations(annotations, 's3_transfer_completed', '../s3-transfer-operations/s3-transfer-jobs*/completed-jobs/*.json')
    read_annotations(annotations, 'qc_donor_prioritization', 'qc_donor_prioritization.txt')
    read_annotations(annotations, 'uuid_to_barcode', 'pc_annotation-tcga_uuid2barcode.tsv')    
    read_annotations(annotations, 'icgc_donor_id', 'pc_annotation-icgc_donor_ids.csv')
    read_annotations(annotations, 'icgc_specimen_id', 'pc_annotation-icgc_specimen_ids.csv')
    read_annotations(annotations, 'icgc_sample_id', 'pc_annotation-icgc_sample_ids.csv')



    # hard-code the file name for now    
    train2_freeze_bams = read_train2_bams('../pcawg-operations/variant_calling/train2-lists/Data_Freeze_Train_2.0_GoogleDocs__2015_04_10_1150.tsv')

    # pre-exclude gnos entries when this option is chosen
    gnos_ids_to_be_excluded = set()
    if exclude_gnos_id_lists:
        files = glob.glob(exclude_gnos_id_lists)
        for fname in files:
            with open(fname) as f:
                for d in f: gnos_ids_to_be_excluded.add(d.rstrip())

    donor_fh = open(donor_output_jsonl_file, 'w')
    bam_fh = open(bam_output_jsonl_file, 'w')
    
    for f in get_xml_files( metadata_dir, conf, repo ):
        f = conf.get('output_dir') + '/__all_metadata_xml/' + f
        gnos_analysis = get_gnos_analysis(f)
        #print (json.dumps(gnos_analysis)) # debug
        if gnos_analysis:
            logger.info( 'processing xml file: {} ...'.format(f) )
            if gnos_analysis.get('analysis_id') and gnos_analysis.get('analysis_id') in gnos_ids_to_be_excluded:
                logger.warning( 'skipping xml file: {} with analysis_id: {}, as it\'s in the list to be excluded' \
                    .format(f, gnos_analysis.get('analysis_id')) )
                continue

            process_gnos_analysis( gnos_analysis, donors, vcf_entries, es_index, es, bam_fh, annotations)
        else:
            logger.warning( 'skipping invalid xml file: {}'.format(f) )

    for donor_id in donors.keys():
        donor = donors[donor_id]

        process_donor(donor, annotations, vcf_entries, conf, train2_freeze_bams)

        # push to Elasticsearch
        es.index(index=es_index, doc_type='donor', id=donor['donor_unique_id'], \
            body=json.loads(json.dumps(donor, default=set_default)), timeout=90 )
        del donor['bam_files']  # prune this before dumping JSON for Keiran
        donor_fh.write(json.dumps(donor, default=set_default) + '\n')

    donor_fh.close()
    bam_fh.close()


def update_vcf_jamboree(infilenames, outfilename):
    seen = set() # just for checking in case there are duplicated lines in jamboree files

    with open(outfilename, 'w') as fout:
        for f_index in infilenames:
            with open(f_index,'r') as fin:
                for line in fin:
                    if len(line.rstrip()) == 0: continue
                    if line in seen:
                        pass
                    else:
                        donor_unique_id, gnos_metadata_url, aliquot_id = str.split(line.rstrip(), '\t')
                        repo, gnos_id = str.split(gnos_metadata_url, 'cghub/metadata/analysisFull/')
                        fout. write(donor_unique_id+'\t'+gnos_id+'\n')
                        seen.add(line)


def read_train2_bams(filename):
    train2_bams = {}

    with open(filename, 'r') as r:
        for line in r:
            if line.startswith('dcc_project_code'): continue
            if len(line.rstrip()) == 0: continue
            dcc_project_code, donor_submitter_id, normal_aliquot_id, normal_aligned_bam_gnos_url,\
                num_tumor_samples, tumor_aliquot_id, tumor_aligned_bam_gnos_urls = str.split(line.rstrip(), '\t')

            normal_repo, normal_gnos_id = str.split(normal_aligned_bam_gnos_url, 'cghub/metadata/analysisFull/')

            train2_bams[dcc_project_code + "::" + donor_submitter_id] = {}
            train2_bams.get(dcc_project_code + "::" + donor_submitter_id)[normal_gnos_id] = \
                {"repo": normal_repo, "aliquot_id": normal_aliquot_id, "specimen_type": "normal"}

            tumor_aliquots = str.split(tumor_aliquot_id, ',')
            tumor_urls = str.split(tumor_aligned_bam_gnos_urls, ',')
            for tumor_aliquot_id, tumor_url in zip(tumor_aliquots, tumor_urls):
                tumor_repo, tumor_gnos_id = str.split(tumor_url, 'cghub/metadata/analysisFull/')
                train2_bams.get(dcc_project_code + "::" + donor_submitter_id)[tumor_gnos_id] = \
                    {"repo": tumor_repo, "aliquot_id": tumor_aliquot_id, "specimen_type": "tumor"}

    return train2_bams


def read_annotations(annotations, type, file_name):

    if type in ['s3_transfer_scheduled', 's3_transfer_completed']:
        annotations[type] = set()
        files = glob.glob(file_name)
        for f in files:
            fname = str.split(f, '/')[-1]
            gnos_id = str.split(fname, '.')[0]
            annotations[type].add(gnos_id)
    else:
        if not os.path.isfile(file_name):
            return
        with open(file_name, 'r') as r:
            if annotations.get(type): # reset annotation if exists
                del annotations[type]

            if type == 'gnos_assignment':
                annotations['gnos_assignment'] = {}
                assignment = yaml.safe_load(r)
                for repo, project_donors in assignment.iteritems():
                    for p_d in project_donors:
                        annotations['gnos_assignment'][p_d] = repo  # key is project or donor unique id, value is repo

            elif type == 'sanger_vcf_in_jamboree':
                annotations['sanger_vcf_in_jamboree'] = {}
                for line in r:
                    if line.startswith('#'): continue
                    if len(line.rstrip()) == 0: continue
                    donor_id, ao_id = str.split(line.rstrip(), '\t')
                    annotations[type][donor_id] = ao_id
                    
            elif type in ['train2_donors', 'train2_pilot', 'donor_blacklist', 'manual_qc_failed']:
                annotations[type] = set()
                for line in r:
                    if line.startswith('#'): continue
                    if len(line.rstrip()) == 0: continue
                    annotations[type].add(line.rstrip())

            elif type in ['santa_cruz', 'aug2015', 'oct2015']:
                annotations[type] = {
                    'donor': set(),
                    'gnos_id': set()
                }
                for line in r:
                    if line.startswith('#'): continue
                    if len(line.rstrip()) == 0: continue
                    donor_unique_id, gnos_id, entry_type = str.split(line.rstrip(), '\t') 
                    annotations[type]['donor'].add(donor_unique_id)
                    annotations[type]['gnos_id'].add(gnos_id)                 
            
            elif type == 'qc_donor_prioritization':
                annotations[type] = {}
                reader = csv.DictReader(r, delimiter='\t')
                for row in reader:
                    annotations[type][row.get('Unique DonorId')] = int(row.get('Issue Summary'))

            elif type == 'uuid_to_barcode':
                annotations[type] = {}
                for line in r:
                    if line.startswith('#'): continue
                    if len(line.rstrip()) == 0: continue
                    TCGA_project, subtype, uuid, barcode = str.split(line.rstrip(), '\t')
                    uuid = detect_and_low_case_uuid(uuid)
                    annotations[type][uuid] = barcode 


            elif type in ['icgc_donor_id', 'icgc_sample_id', 'icgc_specimen_id']:
                annotations[type] = {}
                subtype = type.split('_')[1]
                prefix = subtype[0:2]
                for line in r:
                    if line.startswith('#'): continue
                    if len(line.rstrip()) == 0: continue
                    icgc_id, id_pcawg, dcc_project_code, creation_release = str.split(line.rstrip(), ',')
                    id_pcawg = detect_and_low_case_uuid(id_pcawg)
                    annotations[type][dcc_project_code+'::'+id_pcawg] = prefix.upper()+icgc_id 

            else:
                logger.warning('unknown annotation type: {}'.format(type))


def process_donor(donor, annotations, vcf_entries, conf, train2_freeze_bams):
    logger.info( 'processing donor: {} ...'.format(donor.get('donor_unique_id')) )

    # check whether all tumor specimen(s) aligned
    if (donor.get('flags').get('aligned_tumor_specimen_aliquot_counts') 
            and donor.get('flags').get('aligned_tumor_specimen_aliquot_counts') == donor.get('flags').get('all_tumor_specimen_aliquot_counts')):
        donor.get('flags')['are_all_tumor_specimens_aligned'] = True

    # now build easy-to-use, specimen-level, gnos_repo-aware summary of bwa alignment status by iterating all collected bams
    aggregated_bam_info = bam_aggregation(donor['bam_files'])
    #print json.dumps(aggregated_bam_info, default=set_default)  # debug only
    
    # let's add this aggregated alignment information to donor object
    if aggregated_bam_info.get('WGS'):
        add_alignment_status_to_donor(donor, aggregated_bam_info.get('WGS'))
        
    #print json.dumps(donor.get('tumor_alignment_status'), default=set_default)  # debug only
    
    #print (json.dumps(aggregated_bam_info.get('RNA-Seq'), default=set_default))  # debug only
    if aggregated_bam_info.get('RNA-Seq'):
        add_rna_seq_status_to_donor(donor, aggregated_bam_info.get('RNA-Seq'))
        if donor.get('rna_seq').get('alignment').get('normal'):
            aliquot = donor.get('rna_seq').get('alignment').get('normal')
            if aliquot.get('tophat'):
                donor.get('flags')['is_normal_tophat_rna_seq_alignment_performed'] = True
            if aliquot.get('star'):
                donor.get('flags')['is_normal_star_rna_seq_alignment_performed'] = True

        if len(donor.get('rna_seq').get('alignment').get('tumor')) > 0:
            for aliquot in donor.get('rna_seq').get('alignment').get('tumor'):
                if aliquot.get('tophat'):
                    donor.get('flags')['is_tumor_tophat_rna_seq_alignment_performed'] = True
                if aliquot.get('star'):
                    donor.get('flags')['is_tumor_star_rna_seq_alignment_performed'] = True

    # # for debug
    # if donor.get('donor_unique_id') == 'OV-AU::AOCS-141':
    #     print json.dumps(aggregated_bam_info.get('RNA-Seq'), default=set_default)
    #     print json.dumps(donor.get('rna_seq').get('alignment'), default=set_default)
    #     sys.exit(0)
        

    if donor.get('normal_alignment_status') and donor.get('normal_alignment_status').get('aligned'):
        donor.get('flags')['is_normal_specimen_aligned'] = True
    
    # add gnos repos where complete alignments for the current donor are available
    add_gnos_repos_with_complete_alignment_set(donor)

    # add gnos repos where one alignment or all alignments for the current donor are available
    add_gnos_repos_with_alignment_result(donor)

    # add original gnos repo assignment, this is based on a manually maintained yaml file
    add_original_gnos_repo(donor, annotations['gnos_assignment'])
    if donor.get('flags').get('is_normal_specimen_aligned') and not donor.get('original_gnos_assignment'):
        logger.warning('donor with normal aligned but gnos_for_originally_aligned_at is empty, please update gnos assignment annotation for donor: {} with {}'
            .format(donor.get('donor_unique_id'), conf.get(donor.get('normal_alignment_status').get('aligned_bam').get('gnos_repo')[0])))
        # it should be pretty safe to assign it automatically for this freshly aligned normal specimen
        donor['original_gnos_assignment'] = conf.get(donor.get('normal_alignment_status').get('aligned_bam').get('gnos_repo')[0])
    add_train2_donor_flag(donor, train2_freeze_bams)
    add_train2_pilot_flag(donor, annotations['train2_pilot'])
    add_donor_blacklist_flag(donor, annotations['donor_blacklist'])
    add_manual_qc_failed_flag(donor, annotations['manual_qc_failed'])
    
    donor.get('flags')['is_sanger_vcf_in_jamboree'] = False
    if donor.get('donor_unique_id') in annotations.get('sanger_vcf_in_jamboree'):
        donor.get('flags')['is_sanger_vcf_in_jamboree'] = True
        donor.get('flags').get('vcf_in_jamboree').append('sanger')

    # choose vcf to vcf_entries by iterating all cached vcfs
    choose_vcf_entry(vcf_entries, donor.get('donor_unique_id'), annotations)

    # re-organize dkfz/embl variant call results, this function does two things:
    # 1. when the combined dkfz/embl call exists, remove the separate ones
    # 2. create combined dkfz/embl variant call entry stub when the real combined
    #    one does not exist yet, but the two separate call results do exist
    reorganize_dkfz_embl_calls(vcf_entries.get(donor.get('donor_unique_id')))

    add_vcf_entry(donor, vcf_entries.get(donor.get('donor_unique_id')))

    check_bwa_duplicates(donor, train2_freeze_bams)


def reorganize_dkfz_embl_calls(vcf_entries):
    if not vcf_entries: return

    variant_call_types = set()
    for key in vcf_entries:
        if not key.endswith('_variant_calling'): continue
        variant_call_types.add(key)

    if 'dkfz_embl_variant_calling' in variant_call_types:
        if vcf_entries.get('embl_variant_calling'):
            logger.warning('Combined dkfz/embl call exists with gnos_id: {}, removing embl call entry with gnos_id: {}'\
                .format(vcf_entries.get('dkfz_embl_variant_calling').get('gnos_id'), vcf_entries.get('embl_variant_calling').get('gnos_id')))
            vcf_entries.pop('embl_variant_calling')

        if vcf_entries.get('dkfz_variant_calling'):
            logger.warning('Combined dkfz/embl call exists with gnos_id: {}, removing dkfz call entry with gnos_id: {}'\
                .format(vcf_entries.get('dkfz_embl_variant_calling').get('gnos_id'), vcf_entries.get('dkfz_variant_calling').get('gnos_id')))
            vcf_entries.pop('dkfz_variant_calling')

    elif 'embl_variant_calling' in variant_call_types and 'dkfz_variant_calling' in variant_call_types:
        # now create the combined dkfz_embl_variant_calling stub
        vcf_entries.update({
                'dkfz_embl_variant_calling': {
                    'gnos_repo': vcf_entries.get('embl_variant_calling').get('gnos_repo'),  # this is needed for reporting purpose, get it from embl
                    'is_stub': True
                }
            })


def check_bwa_duplicates(donor, train2_freeze_bams):
    duplicated_bwa_alignment_summary = {
        'exists_gnos_xml_mismatch': False,
        'exists_gnos_xml_mismatch_in_normal': False,
        'exists_gnos_xml_mismatch_in_tumor': False,
        'exists_mismatch_bwa_bams': False,
        'exists_mismatch_bwa_bams_in_normal': False,
        'exists_mismatch_bwa_bams_in_tumor': False,
        'exists_gnos_id_mismatch': False,
        'exists_gnos_id_mismatch_in_normal': False,
        'exists_gnos_id_mismatch_in_tumor': False,
        'exists_md5sum_mismatch': False,
        'exists_md5sum_mismatch_in_normal': False,
        'exists_md5sum_mismatch_in_tumor': False,
        'exists_version_mismatch': False,
        'exists_version_mismatch_in_normal': False,
        'exists_version_mismatch_in_tumor': False,
        'exists_md5sum_mismatch_between_train2_marked_and_sanger_used': False,
        'exists_version_mismatch_between_train2_marked_and_sanger_used': False,
        'is_santa_cruz_freeze_bam_missing': False,
        'is_santa_cruz_freeze_normal_bam_missing': False,
        'is_santa_cruz_freeze_tumor_bam_missing': False,
        'is_train2_freeze_bam_missing': False,
        'is_train2_freeze_normal_bam_missing': False,
        'is_train2_freeze_tumor_bam_missing': False,
        'is_bam_used_by_sanger_missing': False,
        'is_normal_bam_used_by_sanger_missing': False,
        'is_tumor_bam_used_by_sanger_missing': False,
        'normal': {},
        '_tmp_tumor': {},
        'tumor': []
    }
    aliquots = {}
    duplicated_bwa = False

    for bam_file in donor.get('bam_files'):
        if not bam_file.get('is_aligned'): continue

        # not do it for RNA-Seq Bams
        if bam_file.get('library_strategy') == 'RNA-Seq': continue

        if aliquots.get(bam_file.get('aliquot_id')): # exists already
            duplicated_bwa = True
            aliquots.get(bam_file.get('aliquot_id')).append(bam_file)
        else:
            aliquots[bam_file.get('aliquot_id')] = [bam_file]

    if True or duplicated_bwa:  # Let's do this for all donors
        for aliquot in aliquots:
          for bam_file in aliquots.get(aliquot):
            if 'normal' in bam_file.get('dcc_specimen_type').lower():
                if duplicated_bwa_alignment_summary.get('normal'):
                    duplicated_bwa_alignment_summary.get('normal').get('aligned_bam').append(
                            {
                                'gnos_id': bam_file.get('bam_gnos_ao_id'),
                                'gnos_repo': bam_file.get('gnos_repo'),
                                'md5sum': bam_file.get('md5sum'),
                                'effective_xml_md5sum': bam_file.get('effective_xml_md5sum'),
                                'upload_date': bam_file.get('upload_date'),
                                'published_date': bam_file.get('published_date'),
                                'last_modified': bam_file.get('last_modified'),
                                'bwa_workflow_version': bam_file.get('alignment').get('workflow_version'),
                                'is_train2_bam': is_train2_bam(donor, train2_freeze_bams, bam_file.get('bam_gnos_ao_id'), 'normal'),
                                'is_used_in_sanger_variant_call': is_used_in_sanger_variant_call(donor,
                                        bam_file.get('bam_gnos_ao_id')),
                                'is_santa_cruz_entry': bam_file.get('is_santa_cruz_entry'),
                                'is_aug2015_entry': bam_file.get('is_aug2015_entry'),
                                'is_oct2015_entry': bam_file.get('is_oct2015_entry')
                            }
                        )
                else:
                    duplicated_bwa_alignment_summary['normal'] = {
                        'aliquot_id': aliquot,
                        'dcc_specimen_type': bam_file.get('dcc_specimen_type'),
                        'aligned_bam': [
                            {
                                'gnos_id': bam_file.get('bam_gnos_ao_id'),
                                'gnos_repo': bam_file.get('gnos_repo'),
                                'md5sum': bam_file.get('md5sum'),
                                'effective_xml_md5sum': bam_file.get('effective_xml_md5sum'),
                                'upload_date': bam_file.get('upload_date'),
                                'published_date': bam_file.get('published_date'),
                                'last_modified': bam_file.get('last_modified'),
                                'bwa_workflow_version': bam_file.get('alignment').get('workflow_version'),
                                'is_train2_bam': is_train2_bam(donor, train2_freeze_bams, bam_file.get('bam_gnos_ao_id'), 'normal'),
                                'is_used_in_sanger_variant_call': is_used_in_sanger_variant_call(donor,
                                        bam_file.get('bam_gnos_ao_id')),
                                'is_santa_cruz_entry': bam_file.get('is_santa_cruz_entry'),
                                'is_aug2015_entry': bam_file.get('is_aug2015_entry'),
                                'is_oct2015_entry': bam_file.get('is_oct2015_entry')
                            }
                        ]
                    }

            else: # tumor
                if not duplicated_bwa_alignment_summary.get('_tmp_tumor').get(aliquot):
                    duplicated_bwa_alignment_summary.get('_tmp_tumor')[aliquot] = {
                        'aliquot_id': aliquot,
                        'dcc_specimen_type': bam_file.get('dcc_specimen_type'),
                        'aligned_bam': []
                    }

                duplicated_bwa_alignment_summary.get('_tmp_tumor').get(aliquot).get('aligned_bam').append(
                        {
                            'gnos_id': bam_file.get('bam_gnos_ao_id'),
                            'gnos_repo': bam_file.get('gnos_repo'),
                            'md5sum': bam_file.get('md5sum'),
                            'effective_xml_md5sum': bam_file.get('effective_xml_md5sum'),
                            'upload_date': bam_file.get('upload_date'),
                            'published_date': bam_file.get('published_date'),
                            'last_modified': bam_file.get('last_modified'),
                            'bwa_workflow_version': bam_file.get('alignment').get('workflow_version'),
                            'is_train2_bam': is_train2_bam(donor, train2_freeze_bams, bam_file.get('bam_gnos_ao_id'), 'tumor'),
                            'is_used_in_sanger_variant_call': is_used_in_sanger_variant_call(donor,
                                    bam_file.get('bam_gnos_ao_id')),
                            'is_santa_cruz_entry': bam_file.get('is_santa_cruz_entry'),
                            'is_aug2015_entry': bam_file.get('is_aug2015_entry'),
                            'is_oct2015_entry': bam_file.get('is_oct2015_entry')
                        }
                    )

        for aliquot in duplicated_bwa_alignment_summary.get('_tmp_tumor'):
            duplicated_bwa_alignment_summary.get('tumor').append(duplicated_bwa_alignment_summary.get('_tmp_tumor').get(aliquot))

        del duplicated_bwa_alignment_summary['_tmp_tumor']

        # scan normal BAMs
        if duplicated_bwa_alignment_summary.get('normal'):
            b_gnos_id = None
            b_md5sum = None
            xml_md5sum = None
            b_version = None
            has_santa_cruz_n_bam = False
            has_train2_n_bam = False
            has_sanger_n_bam = False
            count_is_train2_not_sanger = 0
            count_not_train2_is_sanger = 0
            count_is_train2_is_sanger = 0
            duplicated_bwa_alignment_summary.get('normal')['exists_mismatch_bwa_bams'] = False
            duplicated_bwa_alignment_summary.get('normal')['exists_gnos_id_mismatch'] = False
            duplicated_bwa_alignment_summary.get('normal')['exists_gnos_xml_mismatch'] = False
            duplicated_bwa_alignment_summary.get('normal')['exists_md5sum_mismatch'] = False
            duplicated_bwa_alignment_summary.get('normal')['exists_version_mismatch'] = False

            for bam in duplicated_bwa_alignment_summary.get('normal').get('aligned_bam'):
                is_santa_cruz_n_bam = bam.get('is_santa_cruz_entry')
                if is_santa_cruz_n_bam: has_santa_cruz_n_bam = True
                is_train2_n_bam = bam.get('is_train2_bam')
                if is_train2_n_bam: has_train2_n_bam = True
                is_sanger_n_bam = bam.get('is_used_in_sanger_variant_call')
                if is_sanger_n_bam: has_sanger_n_bam = True

                if is_train2_n_bam and not is_sanger_n_bam: count_is_train2_not_sanger += 1
                if not is_train2_n_bam and is_sanger_n_bam: count_not_train2_is_sanger += 1
                if is_train2_n_bam and is_sanger_n_bam: count_is_train2_is_sanger += 1

                if not b_gnos_id: b_gnos_id = bam.get('gnos_id')
                if b_gnos_id and not b_gnos_id == bam.get('gnos_id'):
                    duplicated_bwa_alignment_summary['exists_gnos_id_mismatch'] = True
                    duplicated_bwa_alignment_summary['exists_gnos_id_mismatch_in_normal'] = True
                    duplicated_bwa_alignment_summary['exists_mismatch_bwa_bams'] = True
                    duplicated_bwa_alignment_summary['exists_mismatch_bwa_bams_in_normal'] = True
                    duplicated_bwa_alignment_summary.get('normal')['exists_mismatch_bwa_bams'] = True
                    duplicated_bwa_alignment_summary.get('normal')['exists_gnos_id_mismatch'] = True

                if not b_md5sum: b_md5sum = bam.get('md5sum')
                if b_md5sum and not b_md5sum == bam.get('md5sum'):
                    duplicated_bwa_alignment_summary['exists_md5sum_mismatch'] = True
                    duplicated_bwa_alignment_summary['exists_md5sum_mismatch_in_normal'] = True
                    duplicated_bwa_alignment_summary['exists_mismatch_bwa_bams'] = True
                    duplicated_bwa_alignment_summary['exists_mismatch_bwa_bams_in_normal'] = True
                    duplicated_bwa_alignment_summary.get('normal')['exists_mismatch_bwa_bams'] = True
                    duplicated_bwa_alignment_summary.get('normal')['exists_md5sum_mismatch'] = True

                if not xml_md5sum: xml_md5sum = bam.get('effective_xml_md5sum')
                if xml_md5sum and not xml_md5sum == bam.get('effective_xml_md5sum'):
                    duplicated_bwa_alignment_summary['exists_gnos_xml_mismatch'] = True
                    duplicated_bwa_alignment_summary['exists_gnos_xml_mismatch_in_normal'] = True
                    duplicated_bwa_alignment_summary.get('normal')['exists_gnos_xml_mismatch'] = True

                if not b_version: b_version = bam.get('bwa_workflow_version')
                if b_version and not b_version == bam.get('bwa_workflow_version'):
                    duplicated_bwa_alignment_summary['exists_version_mismatch'] = True
                    duplicated_bwa_alignment_summary['exists_version_mismatch_in_normal'] = True
                    duplicated_bwa_alignment_summary['exists_mismatch_bwa_bams'] = True
                    duplicated_bwa_alignment_summary['exists_mismatch_bwa_bams_in_normal'] = True
                    duplicated_bwa_alignment_summary.get('normal')['exists_mismatch_bwa_bams'] = True
                    duplicated_bwa_alignment_summary.get('normal')['exists_version_mismatch'] = True

            if donor.get('flags').get('is_santa_cruz_donor') and not has_santa_cruz_n_bam:
                duplicated_bwa_alignment_summary['is_santa_cruz_freeze_bam_missing'] = True
                duplicated_bwa_alignment_summary['is_santa_cruz_freeze_normal_bam_missing'] = True

            if donor.get('flags').get('is_train2_donor') and not has_train2_n_bam:
                duplicated_bwa_alignment_summary['is_train2_freeze_bam_missing'] = True
                duplicated_bwa_alignment_summary['is_train2_freeze_normal_bam_missing'] = True

            if donor.get('flags').get('is_sanger_variant_calling_performed') and not has_sanger_n_bam:
                duplicated_bwa_alignment_summary['is_bam_used_by_sanger_missing'] = True
                duplicated_bwa_alignment_summary['is_normal_bam_used_by_sanger_missing'] = True

            if donor.get('flags').get('is_train2_donor') and \
                    donor.get('flags').get('is_sanger_variant_calling_performed') and \
                    not count_is_train2_is_sanger and \
                    count_is_train2_not_sanger and count_not_train2_is_sanger:
                if duplicated_bwa_alignment_summary['exists_md5sum_mismatch']:
                    duplicated_bwa_alignment_summary['exists_md5sum_mismatch_between_train2_marked_and_sanger_used'] = True
                if duplicated_bwa_alignment_summary['exists_version_mismatch']:
                    duplicated_bwa_alignment_summary['exists_version_mismatch_between_train2_marked_and_sanger_used'] = True

        # scan tumor BAMs
        if duplicated_bwa_alignment_summary.get('tumor'):
            for aliquot in duplicated_bwa_alignment_summary.get('tumor'):
                b_gnos_id = None
                b_md5sum = None
                xml_md5sum = None
                b_version = None
                has_santa_cruz_t_bam = False
                has_train2_t_bam = False
                has_sanger_t_bam = False
                count_is_train2_not_sanger = 0
                count_not_train2_is_sanger = 0
                count_is_train2_is_sanger = 0
                aliquot['exists_mismatch_bwa_bams'] = False
                aliquot['exists_gnos_id_mismatch'] = False
                aliquot['exists_gnos_xml_mismatch'] = False
                aliquot['exists_md5sum_mismatch'] = False
                aliquot['exists_version_mismatch'] = False

                for bam in aliquot.get('aligned_bam'):
                    is_santa_cruz_t_bam = bam.get('is_santa_cruz_entry')
                    if is_santa_cruz_t_bam: has_santa_cruz_t_bam = True

                    is_train2_t_bam = bam.get('is_train2_bam')
                    if is_train2_t_bam: has_train2_t_bam = True

                    is_sanger_t_bam = bam.get('is_used_in_sanger_variant_call')
                    if is_sanger_t_bam: has_sanger_t_bam = True

                    if is_train2_t_bam and not is_sanger_t_bam: count_is_train2_not_sanger += 1
                    if not is_train2_t_bam and is_sanger_t_bam: count_not_train2_is_sanger += 1
                    if is_train2_t_bam and is_sanger_t_bam: count_is_train2_is_sanger += 1

                    if not b_gnos_id: b_gnos_id = bam.get('gnos_id')
                    if b_gnos_id and not b_gnos_id == bam.get('gnos_id'):
                        duplicated_bwa_alignment_summary['exists_gnos_id_mismatch'] = True
                        duplicated_bwa_alignment_summary['exists_gnos_id_mismatch_in_tumor'] = True
                        duplicated_bwa_alignment_summary['exists_mismatch_bwa_bams'] = True
                        duplicated_bwa_alignment_summary['exists_mismatch_bwa_bams_in_tumor'] = True
                        aliquot['exists_mismatch_bwa_bams'] = True
                        aliquot['exists_gnos_id_mismatch'] = True

                    if not b_md5sum: b_md5sum = bam.get('md5sum')
                    if b_md5sum and not b_md5sum == bam.get('md5sum'):
                        duplicated_bwa_alignment_summary['exists_md5sum_mismatch'] = True
                        duplicated_bwa_alignment_summary['exists_md5sum_mismatch_in_tumor'] = True
                        duplicated_bwa_alignment_summary['exists_mismatch_bwa_bams'] = True
                        duplicated_bwa_alignment_summary['exists_mismatch_bwa_bams_in_tumor'] = True
                        aliquot['exists_mismatch_bwa_bams'] = True
                        aliquot['exists_md5sum_mismatch'] = True

                    if not xml_md5sum: xml_md5sum = bam.get('effective_xml_md5sum')
                    if xml_md5sum and not xml_md5sum == bam.get('effective_xml_md5sum'):
                        duplicated_bwa_alignment_summary['exists_gnos_xml_mismatch'] = True
                        duplicated_bwa_alignment_summary['exists_gnos_xml_mismatch_in_tumor'] = True
                        aliquot['exists_gnos_xml_mismatch'] = True

                    if not b_version: b_version = bam.get('bwa_workflow_version')
                    if b_version and not b_version == bam.get('bwa_workflow_version'):
                        duplicated_bwa_alignment_summary['exists_version_mismatch'] = True
                        duplicated_bwa_alignment_summary['exists_version_mismatch_in_tumor'] = True
                        duplicated_bwa_alignment_summary['exists_mismatch_bwa_bams'] = True
                        duplicated_bwa_alignment_summary['exists_mismatch_bwa_bams_in_tumor'] = True
                        aliquot['exists_version_mismatch'] = True
                        aliquot['exists_mismatch_bwa_bams'] = True

                if donor.get('flags').get('is_santa_cruz_donor') and not has_santa_cruz_t_bam:
                    duplicated_bwa_alignment_summary['is_santa_cruz_freeze_bam_missing'] = True
                    duplicated_bwa_alignment_summary['is_santa_cruz_freeze_tumor_bam_missing'] = True

                if donor.get('flags').get('is_train2_donor') and not has_train2_t_bam:
                    duplicated_bwa_alignment_summary['is_train2_freeze_bam_missing'] = True
                    duplicated_bwa_alignment_summary['is_train2_freeze_tumor_bam_missing'] = True

                if donor.get('flags').get('is_sanger_variant_calling_performed') and not has_sanger_t_bam:
                    duplicated_bwa_alignment_summary['is_bam_used_by_sanger_missing'] = True
                    duplicated_bwa_alignment_summary['is_tumor_bam_used_by_sanger_missing'] = True

                if donor.get('flags').get('is_train2_donor') and \
                        donor.get('flags').get('is_sanger_variant_calling_performed') and \
                        not count_is_train2_is_sanger and \
                        count_is_train2_not_sanger and count_not_train2_is_sanger:
                    if duplicated_bwa_alignment_summary['exists_md5sum_mismatch']:
                        duplicated_bwa_alignment_summary['exists_md5sum_mismatch_between_train2_marked_and_sanger_used'] = True
                    if duplicated_bwa_alignment_summary['exists_version_mismatch']:
                        duplicated_bwa_alignment_summary['exists_version_mismatch_between_train2_marked_and_sanger_used'] = True

        donor['duplicated_bwa_alignment_summary'] = duplicated_bwa_alignment_summary


def is_used_in_sanger_variant_call(donor, gnos_id):
    if donor.get('variant_calling_results') and donor.get('variant_calling_results').get('sanger_variant_calling'):
        for input_gnos_entry in donor.get('variant_calling_results').get('sanger_variant_calling') \
                .get('workflow_details').get('variant_pipeline_input_info'):
            if gnos_id == input_gnos_entry.get('attributes').get('analysis_id'): return True

    return False


def is_train2_bam(donor, train2_freeze_bams, gnos_id, specimen_type):
    if donor.get('donor_unique_id') and train2_freeze_bams.get(donor.get('donor_unique_id')) \
            and train2_freeze_bams.get(donor.get('donor_unique_id')).get(gnos_id):
        if not specimen_type == train2_freeze_bams.get(donor.get('donor_unique_id')).get(gnos_id).get('specimen_type'):
            logger.warning('This should never happen: specimen type mismatch in train2 list in donor {}'
                    .format(donor.get('donor_unique_id')))
        return True
    return False


def add_vcf_entry(donor, vcf_entry):
    if not vcf_entry:
        return

    if not donor.get('variant_calling_results'): donor['variant_calling_results'] = {}

    donor['vcf_files'] = copy.deepcopy(vcf_entry.get('vcf_entry_files'))
    del vcf_entry['vcf_entry_files']
    donor.get('variant_calling_results').update(vcf_entry)

    # update the flags inside each vcf
    for workflow in ['sanger', 'embl', 'dkfz', 'dkfz_embl', 'broad', 'muse', 'broad_tar']:
      if donor.get('variant_calling_results').get(workflow + '_variant_calling'):

        # if this is a stub for dkfz_embl call, skip the rest
        if workflow == 'dkfz_embl' and donor.get('variant_calling_results').get(workflow + '_variant_calling').get('is_stub'): continue

        # add code to handle 'DKFZ_EMBL_Merged' workflow: use the dkfz output as the merged_workflow output
        if workflow == 'dkfz_embl' and \
            donor.get('variant_calling_results').get(workflow + '_variant_calling').get('workflow_details').get('variant_workflow_name') == 'DKFZ_EMBL_Merged':
            vcf_output_list = donor.get('variant_calling_results').get(workflow + '_variant_calling').get('workflow_details').get('variant_pipeline_output_info').get('dkfz').get('workflow_outputs')
        else:
            vcf_output_list = donor.get('variant_calling_results').get(workflow + '_variant_calling').get('workflow_details').get('variant_pipeline_output_info')

        if not donor.get('flags').get('all_tumor_specimen_aliquot_counts') + 1 == len(vcf_output_list):
            logger.warning(workflow + ' variant calling workflow may have missed tumour specimen for donor: {}'
                    .format(donor.get('donor_unique_id')))
            donor.get('variant_calling_results').get(workflow + '_variant_calling')['is_output_and_tumour_specimen_counts_mismatch'] = True
        else:
            donor.get('variant_calling_results').get(workflow + '_variant_calling')['is_output_and_tumour_specimen_counts_mismatch'] = False

        # add the flags of is_bam_used_by_{{workflow}}_missing, is_normal_bam_used_by_{{workflow}}_missing, is_tumor_bam_used_by_{{workflow}}_missing
        donor.get('variant_calling_results').get(workflow + '_variant_calling')['is_bam_used_by_' + workflow + '_missing'] = False
        donor.get('variant_calling_results').get(workflow + '_variant_calling')['is_normal_bam_used_by_' + workflow + '_missing'] = False
        donor.get('variant_calling_results').get(workflow + '_variant_calling')['is_tumor_bam_used_by_' + workflow + '_missing'] = False
        has_n_bam = False
        vcf_input_t_bam = set()
        tumor_alignment_bam = set()

        # scan all the vcf input  under "variant_pipeline_input_info"
        for vcf_input in donor.get('variant_calling_results').get(workflow + '_variant_calling').get('workflow_details').get('variant_pipeline_input_info'):
           if 'normal' in vcf_input.get('attributes').get('dcc_specimen_type').lower():
               # added more checks to avoid key not exist error
               if donor.get('normal_alignment_status') and \
                   donor.get('normal_alignment_status').get('aligned_bam') and \
                   donor.get('normal_alignment_status').get('aligned_bam').get('gnos_id') and \
                   vcf_input.get('attributes').get('analysis_id') == donor.get('normal_alignment_status').get('aligned_bam').get('gnos_id'): #check normal alignment
                   has_n_bam = True
           elif 'tumour' in vcf_input.get('attributes').get('dcc_specimen_type').lower(): # check the tumor
               vcf_input_t_bam.add((vcf_input.get('specimen'), vcf_input.get('attributes').get('analysis_id'))) 

           else:
               logger.warning('invalid specimen type: {} in donor: {} with aliquot_id: {}'
                    .format(vcf_input.get('attributes').get('dcc_specimen_type'), donor.get('donor_unique_id'), vcf_input.get('specimen'))
                )

        # scan all the bams in tumor_alignment_status
        if donor.get('tumor_alignment_status'):
            for tumor_alignment in donor.get('tumor_alignment_status'):
                if tumor_alignment.get('aligned_bam') and tumor_alignment.get('aligned_bam').get('gnos_id'):  # avoid key not exist error
                    tumor_alignment_bam.add((tumor_alignment.get('aliquot_id'), tumor_alignment.get('aligned_bam').get('gnos_id')))

        if not has_n_bam:
            donor.get('variant_calling_results').get(workflow + '_variant_calling')['is_normal_bam_used_by_' + workflow + '_missing'] = True
            donor.get('variant_calling_results').get(workflow + '_variant_calling')['is_bam_used_by_' + workflow + '_missing'] = True
            donor.get('flags')['is_bam_used_by_variant_calling_missing'] = True

        if vcf_input_t_bam != tumor_alignment_bam:
            donor.get('variant_calling_results').get(workflow + '_variant_calling')['is_tumor_bam_used_by_' + workflow + '_missing'] = True
            donor.get('variant_calling_results').get(workflow + '_variant_calling')['is_bam_used_by_' + workflow + '_missing'] = True
            donor.get('flags')['is_bam_used_by_variant_calling_missing'] = True

        # add the flags of exists_{workflow}_file_prefix_mismatch
        donor.get('variant_calling_results').get(workflow + '_variant_calling')['exists_' + workflow + '_file_prefix_mismatch'] = False

        # update the flags of exists_xml_md5sum_mismatch
        if donor.get('variant_calling_results').get(workflow + '_variant_calling').get('exists_xml_md5sum_mismatch'):
            donor.get('flags')['exists_xml_md5sum_mismatch'] = True

        # scan all the files under **_variant_calling
        file_prefix = set()
        for f in donor.get('variant_calling_results').get(workflow + '_variant_calling').get('files'):
            file_prefix.add(f.get('file_name').split('.')[0])
        if not file_prefix == donor.get('all_tumor_specimen_aliquots'):
            donor.get('variant_calling_results').get(workflow + '_variant_calling')['exists_' + workflow + '_file_prefix_mismatch'] = True
            donor.get('flags')['exists_vcf_file_prefix_mismatch'] = True
    
    # update the flags for sanger, dkfz_embl
    for workflow in ['sanger', 'dkfz_embl', 'embl', 'dkfz']:
        if donor.get('variant_calling_results').get(workflow + '_variant_calling'):
            donor.get('flags')['is_' + workflow + '_variant_calling_performed'] = True
            donor.get('flags').get('variant_calling_performed').append(workflow)

    #one combined flag for broad to indicate whether broad is performed well
    #donor.get('flags')['exists_mismatch_broad_file_subsets'] = False
    is_broad_file_subset_missing = False
    broad_file_subsets = set()
    for workflow in ['broad', 'muse', 'broad_tar']:
        if donor.get('variant_calling_results').get(workflow + '_variant_calling'):
            donor.get('flags').get('broad')[workflow+'_file_subset_exist'] = True
            vcf = donor.get('variant_calling_results').get(workflow + '_variant_calling')
        
            if not vcf.get('workflow_details') or not vcf.get('workflow_details').get('related_file_subset_uuids') or not vcf.get('gnos_id'):
                logger.warning('{} variant calling information for donor: {} is not completely populated'.format(workflow.upper(), donor.get('donor_unique_id')))
            else:
                current_broad_file_subsets = set(vcf.get('workflow_details').get('related_file_subset_uuids')) | set([vcf.get('gnos_id')])
            
                if not broad_file_subsets: broad_file_subsets = current_broad_file_subsets
                if broad_file_subsets and not current_broad_file_subsets == broad_file_subsets: 
                    donor.get('flags').get('broad')['exist_file_subsets_mismatch'] = True
        else:
            is_broad_file_subset_missing = True

    if not is_broad_file_subset_missing and not donor.get('flags').get('broad')['exist_file_subsets_mismatch']:
        donor.get('flags')['is_broad_variant_calling_performed'] = True
        donor.get('flags').get('variant_calling_performed').append('broad')


def add_original_gnos_repo(donor, annotation):
    if donor.get('gnos_repo'):
        del donor['gnos_repo']  # get rid of this rather confusing old flag

    if annotation.get(donor.get('donor_unique_id')):
        donor['original_gnos_assignment'] = annotation.get(donor.get('donor_unique_id'))
    elif annotation.get(donor.get('dcc_project_code')):
        donor['original_gnos_assignment'] = annotation.get(donor.get('dcc_project_code'))
    else:
        donor['original_gnos_assignment'] = None


def add_train2_donor_flag(donor, train2_freeze_bams):
    if train2_freeze_bams.get(donor.get('donor_unique_id')):
        donor.get('flags')['is_train2_donor'] = True
    else:
        donor.get('flags')['is_train2_donor'] = False


def add_train2_pilot_flag(donor, annotation):
    if donor.get('donor_unique_id') in annotation:
        donor.get('flags')['is_train2_pilot'] = True
    else:
        donor.get('flags')['is_train2_pilot'] = False


def add_donor_blacklist_flag(donor, annotation):
    if donor.get('donor_unique_id') in annotation:
        donor.get('flags')['is_donor_blacklisted'] = True
    else:
        donor.get('flags')['is_donor_blacklisted'] = False


def add_manual_qc_failed_flag(donor, annotation):
    if donor.get('donor_unique_id') in annotation:
        donor.get('flags')['is_manual_qc_failed'] = True
    else:
        donor.get('flags')['is_manual_qc_failed'] = False

def add_gnos_repos_with_alignment_result(donor):
    repos = set()

    if (donor.get('normal_alignment_status')
            and donor.get('normal_alignment_status').get('aligned_bam')):
        repos = set(donor.get('normal_alignment_status').get('aligned_bam').get('gnos_repo'))

    if donor.get('tumor_alignment_status'):
        for t in donor.get('tumor_alignment_status'):
            if t.get('aligned_bam'):
                repos.update(set(t.get('aligned_bam').get('gnos_repo')))

    donor['gnos_repos_with_alignment_result'] = repos


def add_gnos_repos_with_complete_alignment_set(donor):
    repos = set()

    if (donor.get('normal_alignment_status')
            and donor.get('normal_alignment_status').get('aligned_bam')):
        repos = set(donor.get('normal_alignment_status').get('aligned_bam').get('gnos_repo'))

    if repos and donor.get('tumor_alignment_status'):
        for t in donor.get('tumor_alignment_status'):
            if t.get('aligned_bam'):
                repos = set.intersection(repos, set(t.get('aligned_bam').get('gnos_repo')))
            else:
                repos = set()
    else:
        repos = set()

    donor['gnos_repos_with_complete_alignment_set'] = repos
    '''
    # this flag is not entirely accurate, disable it for now
    if repos:
        donor['is_alignment_completed'] = True
    else:
        donor['is_alignment_completed'] = False
    '''

def add_rna_seq_status_to_donor(donor, aggregated_bam_info):
    for aliquot_id in aggregated_bam_info.keys():
        alignment_status = aggregated_bam_info.get(aliquot_id)
        if (alignment_status.get('tophat') and 'normal' in alignment_status.get('tophat').get('dcc_specimen_type').lower()) or \
           (alignment_status.get('star') and 'normal' in alignment_status.get('star').get('dcc_specimen_type').lower()): # normal specimen
            if not donor.get('rna_seq').get('alignment').get('normal'): #no normal yet in RNA-Seq alignment of this donor
                donor.get('rna_seq').get('alignment')['normal'] = alignment_status
                if alignment_status.get('tophat') and alignment_status.get('tophat').get('exists_xml_md5sum_mismatch') or \
                   alignment_status.get('star') and alignment_status.get('star').get('exists_xml_md5sum_mismatch'):
                    donor.get('flags')['exists_xml_md5sum_mismatch'] = True
            else:
                logger.warning('more than one RNA-Seq normal aliquot found in donor: {}'.format(donor.get('donor_unique_id')))

        elif (alignment_status.get('tophat') and 'tumour' in alignment_status.get('tophat').get('dcc_specimen_type').lower()) or \
           (alignment_status.get('star') and 'tumour' in alignment_status.get('star').get('dcc_specimen_type').lower()): 
            if not donor.get('rna_seq').get('alignment').get('tumor'): #no tumor yet in RNA-Seq alignment of this donor
                donor.get('rna_seq').get('alignment')['tumor'] = []
            donor.get('rna_seq').get('alignment')['tumor'].append(copy.deepcopy(alignment_status))
            if alignment_status.get('tophat') and alignment_status.get('tophat').get('exists_xml_md5sum_mismatch') or \
                   alignment_status.get('star') and alignment_status.get('star').get('exists_xml_md5sum_mismatch'):
                donor.get('flags')['exists_xml_md5sum_mismatch'] = True   
        else:
            logger.warning('invalid aliquot_id: {} in donor: {} '
                    .format(aliquot_id, donor.get('donor_unique_id'))
                )


def add_alignment_status_to_donor(donor, aggregated_bam_info):
    for aliquot_id in aggregated_bam_info.keys():
        alignment_status = aggregated_bam_info.get(aliquot_id)
        if 'normal' in alignment_status.get('dcc_specimen_type').lower(): # normal specimen
            if not donor.get('normal_alignment_status'): # no normal yet in this donor, this is good
                donor['normal_alignment_status'] = reorganize_unaligned_bam_info(alignment_status)
                if alignment_status.get('exists_xml_md5sum_mismatch'):
                    donor.get('flags')['exists_xml_md5sum_mismatch'] = True
            else: # another normal with different aliquot_id! this is no good
                logger.warning('donor: {} has more than one normal, in use aliquot_id: {}, additional aliquot_id found: {}'
                        .format(donor.get('donor_unique_id'),
                                donor.get('normal_alignment_status').get('aliquot_id'),
                                aliquot_id)
                    )
        elif 'tumour' in alignment_status.get('dcc_specimen_type').lower(): # tumour specimen
            if not donor.get('tumor_alignment_status'):
                donor['tumor_alignment_status'] = []             
                _tmp_sample_id = []
            donor['tumor_alignment_status'].append(reorganize_unaligned_bam_info(alignment_status))
            if alignment_status.get('exists_xml_md5sum_mismatch'):
                donor.get('flags')['exists_xml_md5sum_mismatch'] = True

            if alignment_status.get('submitter_sample_id') not in _tmp_sample_id:               
                _tmp_sample_id.append(alignment_status.get('submitter_sample_id'))
            else:
                index = _tmp_sample_id.index(alignment_status.get('submitter_sample_id'))
                logger.warning('donor: {} has more than one aliquot_ids in tumour with the same submitter_sample_id: {}, one aliquot_id: {}, additional aliquot_id found: {}'
                        .format(donor.get('donor_unique_id'),
                                alignment_status.get('submitter_sample_id'),
                                donor.get('tumor_alignment_status')[index].get('aliquot_id'),
                                aliquot_id))
        else:
            logger.warning('invalid specimen type: {} in donor: {} with aliquot_id: {}'
                    .format(alignment_status.get('dcc_specimen_type'), donor.get('donor_unique_id'), aliquot_id)
                )


def update_lane_count_flags(alignment_status):
    if len(alignment_status.get('lane_count')) == 1:
        alignment_status['do_lane_counts_in_every_bam_entry_match'] = True
        if str(len(alignment_status.get('unaligned_bams'))) in alignment_status.get('lane_count'):
            alignment_status['do_lane_count_and_bam_count_match'] = True
    return alignment_status


def reorganize_unaligned_bam_info(alignment_status):
    unaligned_bams = []
    for gnos_id in alignment_status.get('unaligned_bams').keys():
        unaligned_bams.append(
            {
                "gnos_id": gnos_id,
                "bam_file_name": alignment_status.get('unaligned_bams').get(gnos_id).get('bam_file_name'),
                "md5sum": alignment_status.get('unaligned_bams').get(gnos_id).get('md5sum'),
                "gnos_repo": alignment_status.get('unaligned_bams').get(gnos_id).get('gnos_repo'),
            }
        )
    alignment_status['unaligned_bams'] = unaligned_bams
    update_lane_count_flags(alignment_status)
    return alignment_status

def create_aggregated_bam_info_dict(bam):
    aggregated_bam_info_dict = {
        "aliquot_id": bam['aliquot_id'],
        "submitter_specimen_id": bam['submitter_specimen_id'],
        "submitter_sample_id": bam['submitter_sample_id'],
        "dcc_specimen_type": bam['dcc_specimen_type'],
        "icgc_specimen_id": bam['icgc_specimen_id'],
        "icgc_sample_id": bam['icgc_sample_id'],         
        "aligned": True,
        "lane_count": set(),
        "do_lane_counts_in_every_bam_entry_match": False,
        "do_lane_count_and_bam_count_match": False,
        "exist_specimen_type_mismatch": False,
        "exist_aligned_bam_specimen_type_mismatch": False,
        "exist_unaligned_bam_specimen_type_mismatch": False,
        "exist_bam_with_unmappable_reads_specimen_type_mismatch": False,
        "exists_xml_md5sum_mismatch": False,
        "aligned_bam": {
            "gnos_id": bam['bam_gnos_ao_id'],
            "bam_file_name": bam['bam_file_name'],
            "bam_file_size": bam['bam_file_size'],
            "bam_file_md5sum": bam['md5sum'],
            "bai_file_name": bam['bai_file_name'],
            "bai_file_size": bam['bai_file_size'],
            "bai_file_md5sum": bam['bai_file_md5sum'],
            "effective_xml_md5sum": [bam['effective_xml_md5sum']],
            "gnos_last_modified": [bam['last_modified']],
            "gnos_repo": [bam['gnos_repo']],
            "is_santa_cruz_entry": bam['is_santa_cruz_entry'],
            "is_aug2015_entry": bam['is_aug2015_entry'],
            "is_oct2015_entry": bam['is_oct2015_entry'],
            "is_s3_transfer_scheduled": bam['is_s3_transfer_scheduled'],
            "is_s3_transfer_completed": bam['is_s3_transfer_completed']
         },
         "bam_with_unmappable_reads": {},
         "unaligned_bams": {}
    }
    
    return aggregated_bam_info_dict


def compare_specimen_type(specimen_type_A, specimen_type_B):
    if 'normal' in specimen_type_A.lower() and 'normal' in specimen_type_B.lower() or \
        'tumour' in specimen_type_A.lower() and 'tumour' in specimen_type_B.lower():
        return True
    else:
        return False



def bam_aggregation(bam_files):
    aggregated_bam_info_new = {}
    if not aggregated_bam_info_new.get('WGS'):
       aggregated_bam_info_new['WGS'] = {}

    aggregated_bam_info = {}

    for bam in bam_files:  # check aligned BAM(s) first
        if not bam['bam_type'] == 'Specimen level aligned BAM':
            continue

        if not aggregated_bam_info.get(bam['aliquot_id']): # new aliquot
            aggregated_bam_info[bam['aliquot_id']] = create_aggregated_bam_info_dict(bam)
        else:
            alignment_status = aggregated_bam_info.get(bam['aliquot_id'])
            
            if not compare_specimen_type(alignment_status.get('dcc_specimen_type'), bam['dcc_specimen_type']):
                alignment_status['exist_aligned_bam_specimen_type_mismatch'] = True
                alignment_status['exist_specimen_type_mismatch'] = True

            if alignment_status.get('aligned_bam').get('gnos_id') == bam['bam_gnos_ao_id']:
                if bam['gnos_repo'] in alignment_status.get('aligned_bam').get('gnos_repo'):
                    logger.warning( 'Same aliquot: {}, same GNOS ID: {} in the same GNOS repo: {} more than once. This should never be possible.'
                                    .format(
                                        bam['aliquot_id'],
                                        alignment_status.get('aligned_bam').get('gnos_id'),
                                        bam['gnos_repo']) 
                              )
                else:
                    alignment_status.get('aligned_bam').get('gnos_repo').append(bam['gnos_repo'])
                    alignment_status.get('aligned_bam').get('gnos_last_modified').append(bam['last_modified'])
                    alignment_status.get('aligned_bam').get('effective_xml_md5sum').append(bam['effective_xml_md5sum'])
                    alignment_status['exists_xml_md5sum_mismatch'] = False if len(set(alignment_status.get('aligned_bam').get('effective_xml_md5sum'))) == 1 else True
                    
            else:
                if bam['is_oct2015_entry']:
                    aggregated_bam_info[bam['aliquot_id']] = create_aggregated_bam_info_dict(bam)
                    logger.info( 'Same aliquot: {} from donor: {} has different aligned GNOS BWA BAM entries, keep the one in oct2015: {}, additional: {}'
                        .format(bam['aliquot_id'], bam['donor_unique_id'], bam['gnos_metadata_url'], alignment_status.get('aligned_bam').get('gnos_id')))
                elif bam['is_s3_transfer_scheduled']:
                    aggregated_bam_info[bam['aliquot_id']] = create_aggregated_bam_info_dict(bam)
                    logger.info( 'Same aliquot: {} from donor: {} has different aligned GNOS BWA BAM entries, keep the one scheduled for S3 transfer: {}, additional: {}'
                        .format(bam['aliquot_id'], bam['donor_unique_id'], bam['gnos_metadata_url'], alignment_status.get('aligned_bam').get('gnos_id')))

                elif bam['is_aug2015_entry']:
                    aggregated_bam_info[bam['aliquot_id']] = create_aggregated_bam_info_dict(bam)
                    logger.info( 'Same aliquot: {} from donor: {} has different aligned GNOS BWA BAM entries, keep the one in aug2015: {}, additional: {}'
                        .format(bam['aliquot_id'], bam['donor_unique_id'], bam['gnos_metadata_url'], alignment_status.get('aligned_bam').get('gnos_id')))


                elif bam['is_santa_cruz_entry']:
                    aggregated_bam_info[bam['aliquot_id']] = create_aggregated_bam_info_dict(bam)
                    logger.info( 'Same aliquot: {} from donor: {} has different aligned GNOS BWA BAM entries, keep the one in santa_cruz: {}, additional: {}'
                        .format(bam['aliquot_id'], bam['donor_unique_id'], bam['gnos_metadata_url'], alignment_status.get('aligned_bam').get('gnos_id')))

                else:
                    logger.warning( 'Same aliquot: {} from donor: {} has different aligned GNOS BWA BAM entries, in use: {}, additional: {}'
                                        .format(
                                            bam['aliquot_id'],
                                            bam['donor_unique_id'],
                                            alignment_status.get('aligned_bam').get('gnos_id'),
                                            bam['gnos_metadata_url'])
                                  )

    sort_repos_by_time(aggregated_bam_info)

    for bam in bam_files:  # now check BAM with unmappable reads that were derived from aligned BAM
        if not bam['bam_type'] == 'Specimen level unmapped reads after BWA alignment':
            continue

        if not aggregated_bam_info.get(bam['aliquot_id']): # new aliquot, too bad this is an orphaned unmapped read BAM the main aligned BAM is missing
            logger.warning('aliquot: {} has GNOS BAM entry for unmapped reads found: {}, however the main aligned BAM entry is missing'
                    .format(bam['aliquot_id'], bam['bam_gnos_ao_id'])
                )
        else:
            alignment_status = aggregated_bam_info.get(bam['aliquot_id'])
            
            if not compare_specimen_type(alignment_status.get('dcc_specimen_type'), bam['dcc_specimen_type']):
                alignment_status['exist_bam_with_unmappable_reads_specimen_type_mismatch'] = True
                alignment_status['exist_specimen_type_mismatch'] = True

            if not alignment_status.get('bam_with_unmappable_reads'):
                alignment_status['bam_with_unmappable_reads'] = {
                    "gnos_id": bam['bam_gnos_ao_id'],
                    "bam_file_name": bam['bam_file_name'],
                    "bam_file_size": bam['bam_file_size'],
                    "gnos_repo": set([bam['gnos_repo']])
                }
            elif alignment_status.get('bam_with_unmappable_reads').get('gnos_id') == bam['bam_gnos_ao_id']:
                alignment_status.get('bam_with_unmappable_reads').get('gnos_repo').add(bam['gnos_repo'])
            else:
                logger.warning( 'same aliquot: {} has different unmappable reads GNOS BWA BAM entries, in use: {}, additional: {}'
                                    .format(
                                        bam['aliquot_id'],
                                        alignment_status.get('bam_with_unmappable_reads').get('gnos_id'),
                                        bam['bam_gnos_ao_id']) 
                              )

    for bam in bam_files:  # last check original (submitted) unaligned BAM(s)
        if not bam['bam_type'] == 'Unaligned BAM':
            continue

        if not aggregated_bam_info.get(bam['aliquot_id']): # new aliquot with no aligned BAM yet
            aggregated_bam_info[bam['aliquot_id']] = {
                "aliquot_id": bam['aliquot_id'],
                "submitter_specimen_id": bam['submitter_specimen_id'],
                "submitter_sample_id": bam['submitter_sample_id'],
                "dcc_specimen_type": bam['dcc_specimen_type'],
                "aligned": False,
                "lane_count": set([bam['total_lanes']]),
                "do_lane_counts_in_every_bam_entry_match": False,
                "do_lane_count_and_bam_count_match": False,  
                "exist_specimen_type_mismatch": False,
                "exist_aligned_bam_specimen_type_mismatch": False,
                "exist_unaligned_bam_specimen_type_mismatch": False,
                "exist_bam_with_unmappable_reads_specimen_type_mismatch": False,              
                "aligned_bam": {},
                "bam_with_unmappable_reads": {},
                "unaligned_bams": {
                    bam['bam_gnos_ao_id']: {
                        "bam_file_name": bam['bam_file_name'],
                        "md5sum": bam['md5sum'],
                        "gnos_repo": set([bam['gnos_repo']])
                    }
                }
            }
        else: # aliquot already exists
            alignment_status = aggregated_bam_info.get(bam['aliquot_id'])
            alignment_status.get('lane_count').add(bam['total_lanes'])

            if not compare_specimen_type(alignment_status.get('dcc_specimen_type'), bam['dcc_specimen_type']):
                alignment_status['exist_unaligned_bam_specimen_type_mismatch'] = True
                alignment_status['exist_specimen_type_mismatch'] = True


            if alignment_status.get('unaligned_bams').get(bam['bam_gnos_ao_id']): # this unaligned bam was encountered before
                if alignment_status.get('unaligned_bams').get(bam['bam_gnos_ao_id']).get('md5sum') == bam['md5sum']: # this unaligned bam has the same md5sum with encountered one
                    alignment_status.get('unaligned_bams').get(bam['bam_gnos_ao_id']).get('gnos_repo').add(bam['gnos_repo'])
                else:
                    logger.warning( 'Unaligend lane-level BAMs with same gnos_id: {} have different md5sum, in use entry at gnos repo: {}, additional entry at gnos repo: {}'
                                    .format(
                                        bam['bam_gnos_ao_id'],
                                        alignment_status.get('unaligned_bams').get(bam['bam_gnos_ao_id']).get('gnos_repo')[-1],
                                        bam['gnos_repo'])
                                )

            else:
                alignment_status.get('unaligned_bams')[bam['bam_gnos_ao_id']] = {
                        "bam_file_name": bam['bam_file_name'],
                        "md5sum": bam['md5sum'],
                        "gnos_repo": set([bam['gnos_repo']])
                }
                

    aggregated_bam_info_new['WGS'] = aggregated_bam_info
    
    aggregated_bam_info = {}

    if not aggregated_bam_info_new.get('RNA-Seq'):
       aggregated_bam_info_new['RNA-Seq'] = {}

    for bam in bam_files:  #check RNA-Seq BAMs
        if not bam['bam_type'] == 'RNA-Seq aligned BAM':
            continue
        if not aggregated_bam_info.get(bam['aliquot_id']):  # new aliquot with RNA-Seq BAM
            aggregated_bam_info[bam['aliquot_id']] = {}
            aliquot_tmp = create_aggregated_rna_bam_info(bam)

            if 'tophat' in bam.get('alignment').get('workflow_name').lower(): 
                aggregated_bam_info.get(bam['aliquot_id'])['tophat'] = aliquot_tmp

            elif 'star' in bam.get('alignment').get('workflow_name').lower():
                aggregated_bam_info.get(bam['aliquot_id'])['star'] = aliquot_tmp
            
            else: # other unknown alignment workflows
                logger.warning('unknown RNA-Seq alignment workflows: {}'
                         .format(bam.get('alignment').get('workflow_name') ))
                return

        else:  #aliquot already exists
            alignment_status = aggregated_bam_info.get(bam['aliquot_id'])
            if 'tophat' in bam.get('alignment').get('workflow_name').lower():
                if not alignment_status.get('tophat'): # no tophat workflow for the aliquot
                    aliquot_tmp = create_aggregated_rna_bam_info(bam)
                    alignment_status['tophat'] = aliquot_tmp

                elif alignment_status.get('tophat').get('aligned_bam').get('gnos_id') == bam['bam_gnos_ao_id']:
                    if bam['gnos_repo'] in alignment_status.get('tophat').get('aligned_bam').get('gnos_repo'):
                        logger.warning( 'Same aliquot: {}, same workflow: {}, same GNOS ID: {} in the same GNOS repo: {} more than once. This should never be possible.'
                                        .format(
                                            bam['aliquot_id'],
                                            bam.get('alignment').get('workflow_name'),
                                            alignment_status.get('tophat').get('aligned_bam').get('gnos_id'),
                                            bam['gnos_repo']) 
                                  )
                    else:
                        alignment_status.get('tophat').get('aligned_bam').get('gnos_repo').append(bam['gnos_repo'])
                        alignment_status.get('tophat').get('aligned_bam').get('gnos_last_modified').append(bam['last_modified'])
                        alignment_status.get('tophat').get('aligned_bam').get('effective_xml_md5sum').append(bam['effective_xml_md5sum'])
                        alignment_status.get('tophat')['exists_xml_md5sum_mismatch'] = False if len(set(alignment_status.get('tophat').get('aligned_bam').get('effective_xml_md5sum'))) == 1 else True


                else:
                    if bam['is_oct2015_entry']:
                        aliquot_tmp = create_aggregated_rna_bam_info(bam)
                        alignment_status['tophat'] = aliquot_tmp
                        logger.info( 'Same aliquot: {} from donor: {} has different tophat aligned GNOS RNA_Seq BAM entries, keep the one in oct2015: {}, additional: {}'
                            .format(bam['aliquot_id'], bam['donor_unique_id'], bam['gnos_metadata_url'], alignment_status.get('tophat').get('aligned_bam').get('gnos_id')))
                    elif bam['is_s3_transfer_scheduled']:
                        aliquot_tmp = create_aggregated_rna_bam_info(bam)
                        alignment_status['tophat'] = aliquot_tmp
                        logger.info( 'Same aliquot: {} from donor: {} has different tophat aligned GNOS RNA_Seq BAM entries, keep the one scheduled for S3 transfer: {}, additional: {}'
                            .format(bam['aliquot_id'], bam['donor_unique_id'], bam['gnos_metadata_url'], alignment_status.get('tophat').get('aligned_bam').get('gnos_id')))

                    elif bam['is_aug2015_entry']:
                        aliquot_tmp = create_aggregated_rna_bam_info(bam)
                        alignment_status['tophat'] = aliquot_tmp
                        logger.info( 'Same aliquot: {} from donor: {} has different tophat aligned GNOS RNA_Seq BAM entries, keep the one in aug2015: {}, additional: {}'
                            .format(bam['aliquot_id'], bam['donor_unique_id'], bam['gnos_metadata_url'], alignment_status.get('tophat').get('aligned_bam').get('gnos_id')))


                    elif bam['is_santa_cruz_entry']:
                        aliquot_tmp = create_aggregated_rna_bam_info(bam)
                        alignment_status['tophat'] = aliquot_tmp
                        logger.info( 'Same aliquot: {} from donor: {} has different tophat aligned GNOS RNA_Seq BAM entries, keep the one in santa_cruz: {}, additional: {}'
                            .format(bam['aliquot_id'], bam['donor_unique_id'], bam['gnos_metadata_url'], alignment_status.get('tophat').get('aligned_bam').get('gnos_id')))

                    else:
                        logger.warning( 'Same aliquot: {} from donor: {} using same workflow: {} has different tophat aligned GNOS RNA_Seq BAM entries, in use: {}, additional: {}'
                                        .format(
                                            bam['aliquot_id'],
                                            bam['donor_unique_id'],
                                            bam.get('alignment').get('workflow_name'),
                                            alignment_status.get('tophat').get('aligned_bam').get('gnos_id'),
                                            bam['gnos_metadata_url'])
                                  )
            elif 'star' in bam.get('alignment').get('workflow_name').lower():
                if not alignment_status.get('star'): # no star workflow for the aliquot
                    aliquot_tmp = create_aggregated_rna_bam_info(bam)
                    alignment_status['star'] = aliquot_tmp

                elif alignment_status.get('star').get('aligned_bam').get('gnos_id') == bam['bam_gnos_ao_id']:
                    if bam['gnos_repo'] in alignment_status.get('star').get('aligned_bam').get('gnos_repo'):
                        logger.warning( 'Same aliquot: {}, same workflow: {}, same GNOS ID: {} in the same GNOS repo: {} more than once. This should never be possible.'
                                        .format(
                                            bam['aliquot_id'],
                                            bam.get('alignment').get('workflow_name'),
                                            alignment_status.get('star').get('aligned_bam').get('gnos_id'),
                                            bam['gnos_repo']) 
                                  )
                    else:
                        alignment_status.get('star').get('aligned_bam').get('gnos_repo').append(bam['gnos_repo'])
                        alignment_status.get('star').get('aligned_bam').get('gnos_last_modified').append(bam['last_modified'])
                        alignment_status.get('star').get('aligned_bam').get('effective_xml_md5sum').append(bam['effective_xml_md5sum'])
                        alignment_status.get('star')['exists_xml_md5sum_mismatch'] = False if len(set(alignment_status.get('star').get('aligned_bam').get('effective_xml_md5sum'))) == 1 else True

                else:
                    if bam['is_oct2015_entry']:
                        aliquot_tmp = create_aggregated_rna_bam_info(bam)
                        alignment_status['star'] = aliquot_tmp
                        logger.info( 'Same aliquot: {} from donor: {} has different star aligned GNOS RNA_Seq BAM entries, keep the one in oct2015: {}, additional: {}'
                            .format(bam['aliquot_id'], bam['donor_unique_id'], bam['gnos_metadata_url'], alignment_status.get('star').get('aligned_bam').get('gnos_id')))
                    elif bam['is_s3_transfer_scheduled']:
                        aliquot_tmp = create_aggregated_rna_bam_info(bam)
                        alignment_status['star'] = aliquot_tmp
                        logger.info( 'Same aliquot: {} from donor: {} has different star aligned GNOS RNA_Seq BAM entries, keep the one scheduled for transfer: {}, additional: {}'
                            .format(bam['aliquot_id'], bam['donor_unique_id'], bam['gnos_metadata_url'], alignment_status.get('star').get('aligned_bam').get('gnos_id')))

                    elif bam['is_aug2015_entry']:
                        aliquot_tmp = create_aggregated_rna_bam_info(bam)
                        alignment_status['star'] = aliquot_tmp
                        logger.info( 'Same aliquot: {} from donor: {} has different star aligned GNOS RNA_Seq BAM entries, keep the one in aug2015: {}, additional: {}'
                            .format(bam['aliquot_id'], bam['donor_unique_id'], bam['gnos_metadata_url'], alignment_status.get('star').get('aligned_bam').get('gnos_id')))


                    elif bam['is_santa_cruz_entry']:
                        aliquot_tmp = create_aggregated_rna_bam_info(bam)
                        alignment_status['star'] = aliquot_tmp
                        logger.info( 'Same aliquot: {} from donor: {} has different star aligned GNOS RNA_Seq BAM entries, keep the one in santa_cruz: {}, additional: {}'
                            .format(bam['aliquot_id'], bam['donor_unique_id'], bam['gnos_metadata_url'], alignment_status.get('star').get('aligned_bam').get('gnos_id')))


                    else:
                        logger.warning( 'Same aliquot: {} from donor: {} using same workflow: {} has different star aligned GNOS RNA_Seq BAM entries, in use: {}, additional: {}'
                                        .format(
                                            bam['aliquot_id'],
                                            bam['donor_unique_id'],
                                            bam.get('alignment').get('workflow_name'),
                                            alignment_status.get('star').get('aligned_bam').get('gnos_id'),
                                            bam['gnos_metadata_url'])
                                    )
            else: # other unknown alignment workflows
                logger.warning('unknown RNA-Seq alignment workflows: {}'
                         .format(bam.get('alignment').get('workflow_name') ))
                return    

    aggregated_bam_info_new['RNA-Seq'] = aggregated_bam_info

    return aggregated_bam_info_new


def create_aggregated_rna_bam_info(bam):
    aliquot_tmp = {
        "aliquot_id": bam['aliquot_id'],
        "submitter_specimen_id": bam['submitter_specimen_id'],
        "submitter_sample_id": bam['submitter_sample_id'],
        "icgc_specimen_id": bam['icgc_specimen_id'],
        "icgc_sample_id": bam['icgc_sample_id'],
        "dcc_specimen_type": bam['dcc_specimen_type'],
        "aligned": True,    
        "is_santa_cruz_entry": bam['is_santa_cruz_entry'],
        "is_aug2015_entry": bam['is_aug2015_entry'],
        "is_oct2015_entry": bam['is_oct2015_entry'],
        "is_s3_transfer_scheduled": bam['is_s3_transfer_scheduled'],  
        "is_s3_transfer_completed": bam['is_s3_transfer_completed'],
        "exists_xml_md5sum_mismatch": False,           
        "aligned_bam": {
            "gnos_repo": [bam['gnos_repo']],
            "gnos_id": bam['bam_gnos_ao_id'],
            "bam_file_name": bam['bam_file_name'],
            "bam_file_md5sum": bam['md5sum'],
            "bam_file_size": bam['bam_file_size'],
            "bai_file_name": bam['bai_file_name'],
            "bai_file_md5sum": bam['bai_file_md5sum'],
            "bai_file_size": bam['bai_file_size'],
            "gnos_last_modified": [bam['last_modified']],
            "effective_xml_md5sum": [bam['effective_xml_md5sum']]
            }
        }
    return aliquot_tmp



def sort_repos_by_time(aggregated_bam_info):
    for aliquot in aggregated_bam_info:
        agg_bam = aggregated_bam_info.get(aliquot)
        if not agg_bam.get('aligned_bam'):
            continue
        modified_dates = agg_bam.get('aligned_bam').get('gnos_last_modified')
        gnos_repos = agg_bam.get('aligned_bam').get('gnos_repo')
        agg_bam.get('aligned_bam')['gnos_last_modified'], agg_bam.get('aligned_bam')['gnos_repo'] = \
            izip(*sorted(izip(modified_dates, gnos_repos), key=lambda x: x[0]))


def find_latest_metadata_dir(output_dir):
    dir_pattern = re.compile(u'^[0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{2}-[0-9]{2}-[0-9]{2}_[A-Z]{3}$')
    metadata_dirs = []
    for dir in os.listdir(output_dir):
        if not os.path.isdir(output_dir + '/' + dir):
            continue
        if dir_pattern.search(dir):
            metadata_dirs.append(output_dir + '/' + dir)

    return sorted(metadata_dirs)[-1]


def main(argv=None):
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    parser = ArgumentParser(description="PCAWG GNOS Metadata Parser",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-c", "--config", dest="config",
             help="Configuration file for GNOS repositories", required=True)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=False)
    parser.add_argument("-r", "--gnos_repo", dest="repo",
             help="Specify which GNOS repo to process, process all repos if none specified", required=False)
    parser.add_argument("-x", "--exclude_gnos_id_lists", dest="exclude_gnos_id_lists", # don't use this option for daily cron job
             help="File(s) containing GNOS IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-s", "--es_index_suffix", dest="es_index_suffix", # don't use this option for daily cron job
             help="Single letter suffix for ES index name", required=False)

    args = parser.parse_args()
    metadata_dir = args.metadata_dir
    conf_file = args.config
    repo = args.repo
    exclude_gnos_id_lists = args.exclude_gnos_id_lists
    es_index_suffix = args.es_index_suffix
    if not es_index_suffix: es_index_suffix = ''

    with open(conf_file) as f:
        conf = yaml.safe_load(f)
        for r in conf.get('gnos_repos'):
            conf[r.get('base_url')] = r.get('repo_code')

    # output_dir
    output_dir = conf.get('output_dir')
    if metadata_dir:
        if not os.path.isdir(metadata_dir):  # TODO: should add more directory name check to make sure it's right
            sys.exit('Error: specified metadata directory does not exist!')
    else:
        metadata_dir = find_latest_metadata_dir(output_dir)  # sorted(glob.glob(output_dir + '/[0-9]*_*_*[A-Z]'))[-1] # find the directory for latest metadata list
    timestamp = str.split(metadata_dir, '/')[-1]

    logger.setLevel(logging.INFO)
    ch.setLevel(logging.WARN)

    log_file = metadata_dir + '.metadata_parser' + ('' if not repo else '.'+repo) + '.log'
    # delete old log first if exists
    if os.path.isfile(log_file): os.remove(log_file)

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)

    es_host = 'localhost:9200'
    es_index = 'p_' + ('' if not repo else repo+'_') + re.sub(r'\D', '', timestamp).replace('20','',1) + es_index_suffix
    es = init_es(es_host, es_index)

    logger.info('processing metadata list files in {} to build es index {}'.format(metadata_dir, es_index))
    process(metadata_dir, conf, es_index, es, metadata_dir+'/donor_'+es_index+'.jsonl', metadata_dir+'/bam_'+es_index+'.jsonl', repo, exclude_gnos_id_lists)

    # now update kibana dashboard
    # donor
    dashboard_name = ' ['+repo+']' if repo else ''
    with open('kibana-donor.json', 'r') as d:
        donor_dashboard = json.loads(d.read())
    donor_dashboard['index']['default'] = es_index + '/donor'
    title = 'PCAWG Donors' + dashboard_name + ' (beta)'
    donor_dashboard['title'] = title
    body = {
        'dashboard': json.dumps(donor_dashboard),
        'user': 'guest',
        'group': 'guest',
        'title': title
    }
    es.index(index='kibana-int', doc_type='dashboard', id='PCAWG Donors' + dashboard_name, body=body)

    # bam search, no need this for now, not very useful
    '''
    with open('kibana-bam.json', 'r') as d:
        bam_dashboard = json.loads(d.read())
    bam_dashboard['index']['default'] = es_index + '/bam_file'
    title = 'PCAWG BAMs' + dashboard_name + ' (beta)'
    bam_dashboard['title'] = title
    body = {
        'dashboard': json.dumps(bam_dashboard),
        'user': 'guest',
        'group': 'guest',
        'title': title
    }
    es.index(index='kibana-int', doc_type='dashboard', id='PCAWG BAMs' + dashboard_name, body=body)
    '''

    return 0


if __name__ == "__main__":
    sys.exit(main())


