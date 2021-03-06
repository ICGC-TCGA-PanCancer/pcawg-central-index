#!/usr/bin/env python


import sys
import os
import re
import glob
import json
import copy
import logging
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from elasticsearch1 import Elasticsearch
from collections import OrderedDict
import datetime
import csv

logger = logging.getLogger('generate PCAWG data release')
ch = logging.StreamHandler()
previous_release = 'mar2016'

es_queries = [
  # query 0: donors_sanger_vcf_without_missing_bams 
    {
     "fields": "donor_unique_id", 
 
    "filter": {
                    "bool": {
                      "must": [
                        {
                          "type": {
                            "value": "donor"
                          }
                        },
                        {
                          "terms":{
                            "flags.is_normal_specimen_aligned":[
                              "T"
                            ]
                          }
                        },
                        {
                          "terms":{
                            "flags.are_all_tumor_specimens_aligned":[
                              "T"
                            ]
                          }
                        }                    
                      ],
                      # "should": [
                      #   {
                      #      "terms":{
                      #         "flags.is_sanger_variant_calling_performed":[
                      #            "T"
                      #         ]
                      #      }
                      #   },
                      #   {
                      #      "bool":{
                      #          "must":[
                      #                 {
                      #                  "terms":{
                      #                     "flags.is_dkfz_embl_variant_calling_performed":[
                      #                        "T"
                      #                     ]
                      #                   }
                      #                 },
                      #                 {
                      #                    "term":{
                      #                       "flags.is_dkfz_variant_calling_performed":[
                      #                          "F"
                      #                       ]
                      #                    }
                      #                 },
                      #                 {
                      #                    "term":{
                      #                       "flags.is_embl_variant_calling_performed":[
                      #                          "F"
                      #                       ]
                      #                    }
                      #                 } 
                      #           ]
                      #      }
                      #   },
                      #   {
                      #      "terms":{
                      #         "flags.is_broad_variant_calling_performed":[
                      #            "T"
                      #         ]
                      #      }
                      #   }
                      # ],
                      "must_not": [
#                        {
#                          "terms": {
#                            "flags.is_bam_used_by_variant_calling_missing": [
#                              "T"
#                            ]
#                          }
#                        },
#                        {
#                          "terms": {
#                            "flags.exists_vcf_file_prefix_mismatch": [
#                              "T"
#                            ]
#                          }
#                        },
                        {
                          "terms": {
                            "duplicated_bwa_alignment_summary.exists_mismatch_bwa_bams": [
                              "T"
                            ]
                          }
                        },
                        # {
                        #    "terms":{
                        #       "flags.exists_xml_md5sum_mismatch":[
                        #          "T"
                        #       ]
                        #    }
                        # },
                        {
                          "terms": {
                            "flags.is_manual_qc_failed": [
                              "T"
                            ]
                          }
                        }
                      ]
                    }
                },
      "size": 10000
    }
]

def create_reorganized_donor(donor_unique_id, es_json, vcf, gnos_ids_to_be_excluded, gnos_ids_to_be_included, annotations):
    reorganized_donor = {
        'donor_unique_id': donor_unique_id,
        'wgs_exclusion_white_gray': 'Whitelist',
        'submitter_donor_id': es_json['submitter_donor_id'],
        'dcc_project_code': es_json['dcc_project_code'],
        'icgc_donor_id': es_json['icgc_donor_id'],
        previous_release+'_donor': True if es_json.get('flags').get('is_'+previous_release+'_donor') else False,
        'santa_cruz_pilot': True if es_json.get('flags').get('is_santa_cruz_donor') else False,
        'validation_by_deep_seq': True if es_json.get('flags').get('is_train2_pilot') else False,
        'TiN': es_json.get('flags').get('TiN'),
        'wgs': {
            'normal_specimen': {},
            'tumor_specimens': []
        },
        'rna_seq': {
             'normal_specimen': {},
             'tumor_specimens': []
        }
    }
    if donor_unique_id in annotations.get('blacklist'):
        reorganized_donor['wgs_exclusion_white_gray'] = 'Excluded'
    elif donor_unique_id in annotations.get('graylist'):
        reorganized_donor['wgs_exclusion_white_gray'] = 'Graylist'
    else:
        pass
    add_wgs_specimens(reorganized_donor, es_json, vcf, gnos_ids_to_be_excluded, gnos_ids_to_be_included, annotations)
    add_rna_seq_info(reorganized_donor, es_json, gnos_ids_to_be_excluded, gnos_ids_to_be_included, annotations)

    return reorganized_donor


def create_alignment(es_json, aliquot, data_type, gnos_ids_to_be_excluded, gnos_ids_to_be_included, annotations, bam_type):
    if not aliquot.get(bam_type): return
    gnos_id = aliquot.get(bam_type).get('gnos_id')
    if gnos_ids_to_be_included and not gnos_id in gnos_ids_to_be_included: return
    if gnos_ids_to_be_excluded and gnos_id in gnos_ids_to_be_excluded: return
    alignment = {
        'submitter_specimen_id': aliquot.get('submitter_specimen_id'),
        'icgc_specimen_id': aliquot.get('icgc_specimen_id'),
        'submitter_sample_id': aliquot.get('submitter_sample_id'),
        'icgc_sample_id': aliquot.get('icgc_sample_id'),
        'specimen_type': aliquot.get('dcc_specimen_type'),
        'aliquot_id': aliquot.get('aliquot_id'),
        'is_'+previous_release+'_entry': aliquot.get(bam_type).get('is_'+previous_release+'_entry') if 'wgs' in data_type else aliquot.get('is_'+previous_release+'_entry'),
        # 'gnos_repo': aliquot.get(bam_type).get('gnos_repo'),
        'gnos_repo': filter_repo(aliquot.get(bam_type).get('gnos_repo'), data_type, bam_type),
        'gnos_id': aliquot.get(bam_type).get('gnos_id'),
        'gnos_last_modified': aliquot.get(bam_type).get('gnos_last_modified')[-1],
        'files': []
    }
    if 'wgs_tumor' in data_type and 'aligned_bam' in bam_type:
        for qc_metric in ['oxog_score', 'ContEST', 'Stars']:
            alignment[qc_metric] = aliquot.get(qc_metric) if aliquot.get(qc_metric) else None

    # add the file info if exist
    for ftype in ['bam', 'bai']:
        if aliquot.get(bam_type).get(ftype+'_file_name'):
            f = {
                'file_name': aliquot.get(bam_type).get(ftype+'_file_name'),
                'file_md5sum': aliquot.get(bam_type).get(ftype+'_file_md5sum'),
                'file_size': aliquot.get(bam_type).get(ftype+'_file_size'),
            }
            alignment.get('files').append(f)
        else:
            logger.warning('{} alignment GNOS entry {} has no .{} file'.format(data_type, alignment.get('gnos_id'), ftype))
    return alignment


def get_key_map(key):
    key_map = {
      "sanger": "sanger_variant_calling",
      "dkfz": "dkfz_embl_variant_calling",
      "embl": "dkfz_embl_variant_calling",
      "broad": "broad_variant_calling",
      "muse": "muse_variant_calling",
      "broad_tar": "broad_tar_variant_calling",
      "aligned_bam": "bwa_alignment",
      "minibam": "minibam"
    }   

    return key_map.get(key)


def choose_variant_calling(es_json, vcf):
    variant_calling = set()
    if not es_json.get('variant_calling_results') or not vcf:
        return variant_calling

    for v in vcf:
        if get_key_map(v) in es_json.get('variant_calling_results').keys() and \
            not es_json.get('variant_calling_results').get(get_key_map(v)).get('is_stub'):
            variant_calling.add(get_key_map(v))
            if not check_vcf(es_json, v): 
                variant_calling.discard(get_key_map(v))

        else:
            logger.warning('donor: {} has no {}'.format(es_json.get('donor_unique_id'), get_key_map(v)))
    return variant_calling


def check_vcf(es_json, vcf_calling):
    if vcf_calling == 'broad' or vcf_calling == 'muse' or vcf_calling == 'broad_tar':
        if not es_json.get('flags').get('is_broad_variant_calling_performed'):
            return False
        elif not es_json.get('variant_calling_results').get(get_key_map('broad')).get('vcf_workflow_result_version') == 'v3':
            return False
        else: 
            return True
    elif vcf_calling == 'sanger':
        if not es_json.get('variant_calling_results').get(get_key_map('sanger')).get('vcf_workflow_result_version') == 'v3':
            return False
        else:
            return True
    else:
        return True

        

def create_variant_calling(es_json, aliquot, wgs_tumor_vcf_info, data_type, gnos_ids_to_be_excluded, gnos_ids_to_be_included):
    gnos_id = wgs_tumor_vcf_info.get('gnos_id')
    if gnos_ids_to_be_included and not gnos_id in gnos_ids_to_be_included: return
    if gnos_ids_to_be_excluded and gnos_id in gnos_ids_to_be_excluded: return
    variant_calling = {
        'submitter_specimen_id': aliquot.get('submitter_specimen_id'),
        'icgc_specimen_id': aliquot.get('icgc_specimen_id'),
        'submitter_sample_id': aliquot.get('submitter_sample_id'),
        'icgc_sample_id': aliquot.get('icgc_sample_id'),
        'specimen_type': aliquot.get('dcc_specimen_type'),
        'aliquot_id': aliquot.get('aliquot_id'),
        'is_'+previous_release+'_entry': wgs_tumor_vcf_info.get('is_'+previous_release+'_entry'),
        'gnos_repo': filter_repo(wgs_tumor_vcf_info.get('gnos_repo'), data_type),
        'gnos_id': wgs_tumor_vcf_info.get('gnos_id'),
        'gnos_last_modified': wgs_tumor_vcf_info.get('gnos_last_modified')[-1],
        'files':[]
    }
    vcf_files = wgs_tumor_vcf_info.get('files')
    if vcf_files:
        for f in vcf_files:
            if int(f.get('file_size')) == 0:
                logger.warning('donor: {} has file: {} file_size is 0'.format(es_json.get('donor_unique_id'), f.get('file_name')))
                continue
            if variant_calling.get('aliquot_id') in f.get('file_name'):
                variant_calling.get('files').append(f)   
    else:
        logger.warning('{} GNOS entry {} has no files'.format(data_type, variant_calling.get('gnos_id')))

    # do not find the vcf files for given aliquot_id
    if not variant_calling.get('files'):
        return None

    return variant_calling


def add_wgs_specimens(reorganized_donor, es_json, vcf, gnos_ids_to_be_excluded, gnos_ids_to_be_included, annotations):
    if es_json.get('normal_alignment_status'):
        wgs_normal_alignment_info = es_json.get('normal_alignment_status')
        for bam_type in ['aligned_bam', 'minibam']:
            data_type = 'wgs_normal_'+bam_type
            reorganized_donor.get('wgs').get('normal_specimen')[get_key_map(bam_type)] = create_alignment(es_json, wgs_normal_alignment_info, data_type, gnos_ids_to_be_excluded, gnos_ids_to_be_included, annotations, bam_type)

    tumor_wgs_specimen_count = 0
    if es_json.get('tumor_alignment_status'):
        variant_calling = choose_variant_calling(es_json, vcf)
        wgs_tumor_alignment_info = es_json.get('tumor_alignment_status')
        consensus_info = {}

        for aliquot in wgs_tumor_alignment_info:
            aliquot_info = OrderedDict()
            tumor_wgs_specimen_count += 1
            for bam_type in ['aligned_bam', 'minibam']:
                data_type = 'wgs_tumor_'+bam_type            
                aliquot_info[get_key_map(bam_type)] = create_alignment(es_json, aliquot, data_type, gnos_ids_to_be_excluded, gnos_ids_to_be_included, annotations, bam_type)        
            
            for vc in variant_calling:
                wgs_tumor_vcf_info = es_json.get('variant_calling_results').get(vc)
                data_type = 'wgs_'+vc        
                aliquot_info[vc] = create_variant_calling(es_json, aliquot, wgs_tumor_vcf_info, data_type, gnos_ids_to_be_excluded, gnos_ids_to_be_included)        
            reorganized_donor.get('wgs').get('tumor_specimens').append(copy.deepcopy(aliquot_info)) 

            # add consensus_variant_call results
            for ct in ['somatic', 'germline']:
                # ignore if there is not specific variant calls
                if not es_json.get('consensus_'+ct+'_variant_calls'): continue
                # initialize the dict
                if not reorganized_donor.get('consensus_'+ct+'_variant_calls'): reorganized_donor['consensus_'+ct+'_variant_calls'] = {}            
                if not consensus_info.get(ct): consensus_info[ct] = {}
                for cc in ['indel', 'snv_mnv', 'sv', 'cnv']:
                    wgs_tumor_consensus_info = es_json.get('consensus_'+ct+'_variant_calls').get(cc)
                    if not wgs_tumor_consensus_info: continue
                    if not consensus_info.get(ct).get(cc): consensus_info[ct][cc] = []

                    for consensus_vcf in wgs_tumor_consensus_info:
                        data_type = 'wgs_consensus_'+ct+'_'+cc
                        consensus_vcf_tmp = create_variant_calling(es_json, aliquot, consensus_vcf, data_type, gnos_ids_to_be_excluded, gnos_ids_to_be_included)
                        # consensus_vcf_tmp is empty if the aliquot was blacklisted aliquot
                        if not consensus_vcf_tmp: continue
                        consensus_info[ct][cc].append(consensus_vcf_tmp)
        
        for ct in ['somatic', 'germline']:
            if not consensus_info.get(ct): continue
            reorganized_donor['consensus_'+ct+'_variant_calls'].update(consensus_info.get(ct))

    reorganized_donor['tumor_wgs_specimen_count'] = tumor_wgs_specimen_count
    return reorganized_donor


def filter_repo(gnos_repo, data_type, bam_type=None):
    # osdc-icgc does not hold wgs aligned_bam
    if bam_type == 'aligned_bam' and 'wgs' in data_type and "https://gtrepo-osdc-icgc.annailabs.com/" in gnos_repo:
        gnos_repo.remove("https://gtrepo-osdc-icgc.annailabs.com/")
    # cghub went offline, all data should be sync-ed to osdc-tcga, hack for now since not all data have been sync-ed. 
    if 'https://cghub.ucsc.edu/' in gnos_repo:
        gnos_repo[gnos_repo.index('https://cghub.ucsc.edu/')] = 'https://gtrepo-osdc-tcga.annailabs.com/'

    return list(set(gnos_repo))  


def add_rna_seq_info(reorganized_donor, es_json, gnos_ids_to_be_excluded, gnos_ids_to_be_included, annotations):
    rna_seq_info = es_json.get('rna_seq').get('alignment')
    for specimen_type in rna_seq_info.keys():
        if not rna_seq_info.get(specimen_type): # the specimen_type has no alignment result
            continue
        
        if 'normal' in specimen_type:
            aliquot = rna_seq_info.get(specimen_type)
            alignment_info = {}
            for workflow_type in aliquot.keys():
                data_type = 'rna_seq_normal_'+workflow_type
                bam_type = 'aligned_bam'
                alignment_info[workflow_type] = create_alignment(es_json, aliquot.get(workflow_type), data_type, gnos_ids_to_be_excluded, gnos_ids_to_be_included, annotations, bam_type)

            reorganized_donor.get('rna_seq')[specimen_type + '_specimen'] = alignment_info
        else:
            tumor_rna_seq_specimen_count = 0
            for aliquot in rna_seq_info.get(specimen_type):
                tumor_rna_seq_specimen_count += 1
                alignment_info = {}
                for workflow_type in aliquot.keys():
                    data_type = 'rna_seq_tumor_'+workflow_type
                    bam_type = 'aligned_bam'
                    alignment_info[workflow_type] = create_alignment(es_json, aliquot.get(workflow_type), data_type, gnos_ids_to_be_excluded, gnos_ids_to_be_included, annotations, bam_type)
                reorganized_donor.get('rna_seq')[specimen_type + '_specimens'].append(copy.deepcopy(alignment_info)) 
            reorganized_donor['tumor_rna_seq_specimen_count'] = tumor_rna_seq_specimen_count

def get_donor_json(es, es_index, donor_unique_id):
    es_query_donor = {
        "query": {
            "term": {
                "donor_unique_id": donor_unique_id
            }
        }
    }
    response = es.search(index=es_index, body=es_query_donor)
    if not response['hits'].get('hits'):
        return None
    es_json = response['hits']['hits'][0]['_source']

    return es_json


def get_donors_list(es, es_index, es_queries):
    q_index = 0
    response = es.search(index=es_index, body=es_queries[q_index])
    donors_list = set()
    for p in response['hits']['hits']:
      donors_list.add(p.get('fields').get('donor_unique_id')[0])
    return donors_list 


def init_es(es_host, es_index):
    es = Elasticsearch([ es_host ])
    es.indices.create( es_index, ignore=400 )
    # create mappings
    es_mapping = open('pcawg_summary.mapping.json')
    es.indices.put_mapping(index=es_index, doc_type='donor', body=es_mapping.read())
    es_mapping.close()
    return es


def set_default(obj):
    if isinstance(obj, datetime.datetime):
        return obj.isoformat()
    if isinstance(obj, set):
        return list(obj)
    raise TypeError


def generate_tsv_file(reorganized_donor, vcf, annotations):
    donor_info = ['donor_unique_id', 'wgs_exclusion_white_gray', 'dcc_project_code', 'submitter_donor_id', 'icgc_donor_id', previous_release+'_donor','santa_cruz_pilot', 'validation_by_deep_seq', 'TiN']
    #specimen = ['submitter_specimen_id', 'icgc_specimen_id', 'submitter_sample_id', 'icgc_sample_id', 'aliquot_id']
    #alignment = ['alignment_gnos_repo', 'alignment_gnos_id', 'alignment_bam_file_name']
        
    pilot_tsv = OrderedDict()
    for d in donor_info:
        pilot_tsv[d] = reorganized_donor.get(d)
    #wgs normal specimen [aligned_bam, minibam]
    alignment = reorganized_donor.get('wgs').get('normal_specimen')
    for workflow_type in ['bwa_alignment', 'minibam']:
        generate_alignment_info(pilot_tsv, alignment, 'normal', 'wgs', workflow_type)  

    # wgs tumor specimen
    wgs_tumor_speciments = reorganized_donor.get('wgs').get('tumor_specimens')
    pilot_tsv['tumor_wgs_specimen_count'] = reorganized_donor.get('tumor_wgs_specimen_count')
    for workflow_type in ['bwa_alignment', 'minibam']:
        generate_alignment_info(pilot_tsv, wgs_tumor_speciments, 'tumor', 'wgs', workflow_type)
    # wgs variant calling
    generate_variant_calling_info(pilot_tsv, wgs_tumor_speciments, vcf, annotations)
    # initial the two columns with flag to show whether wgs has matched RNA-Seq
    for specimen_type in ['normal', 'tumor']:
        pilot_tsv[specimen_type+'_wgs_has_matched_rna_seq'] = []

    # rna_seq normal
    alignment = reorganized_donor.get('rna_seq').get('normal_specimen')
    for workflow in ['star', 'tophat']:
        generate_alignment_info(pilot_tsv, alignment, 'normal', 'rna_seq', workflow.lower()+'_alignment')
    # rna_seq tumor
    rna_seq_tumor = reorganized_donor.get('rna_seq').get('tumor_specimens')
    pilot_tsv['tumor_rna_seq_specimen_count'] = reorganized_donor.get('tumor_rna_seq_specimen_count')  
    for workflow in ['star', 'tophat']:
        generate_alignment_info(pilot_tsv, rna_seq_tumor, 'tumor', 'rna_seq', workflow.lower()+'_alignment')
    #populate the two columns with values
    for specimen_type in ['normal', 'tumor']:
        pilot_tsv[specimen_type+'_wgs_has_matched_rna_seq'] = [t in pilot_tsv[specimen_type+'_rna_seq_icgc_specimen_id'] for t in pilot_tsv[specimen_type+'_wgs_icgc_specimen_id']]

    return pilot_tsv


def generate_variant_calling_info(pilot_tsv, variant_calling, vcf, annotations):
    for v in vcf:       
        pilot_tsv[get_key_map(v)+'_repo'] = []
        pilot_tsv[get_key_map(v)+'_gnos_id'] = []
        pilot_tsv[get_key_map(v)+'_file_name_prefix'] = []
        pilot_tsv['is_'+previous_release+'_'+get_key_map(v)] = []
        if v in ['sanger', 'dkfz', 'broad']:
            pilot_tsv[get_key_map(v)+'_deprecated_gnos_id'] = annotations.get('deprecated_gnos_id').get(pilot_tsv.get('donor_unique_id')).get(v) \
                if annotations.get('deprecated_gnos_id').get(pilot_tsv.get('donor_unique_id')) and annotations.get('deprecated_gnos_id').get(pilot_tsv.get('donor_unique_id')).get(v) else None
        for specimen in variant_calling:
            if specimen.get(get_key_map(v)):
                pilot_tsv[get_key_map(v)+'_repo'] = [specimen.get(get_key_map(v)).get('gnos_repo')]
                pilot_tsv[get_key_map(v)+'_gnos_id'] = specimen.get(get_key_map(v)).get('gnos_id')
                pilot_tsv[get_key_map(v)+'_file_name_prefix'].append(specimen.get(get_key_map(v)).get('aliquot_id'))
                pilot_tsv['is_'+previous_release+'_'+get_key_map(v)] = specimen.get(get_key_map(v)).get('is_'+previous_release+'_entry')
    return pilot_tsv


def generate_alignment_info(pilot_tsv, alignment, specimen_type, sequence_type, workflow_type):
    aliquot_field = ['submitter_specimen_id', 'icgc_specimen_id', 'submitter_sample_id', 'icgc_sample_id', 'aliquot_id']
    gnos_field = ['gnos_repo', 'gnos_id']
    for d in aliquot_field:
        if pilot_tsv.get(specimen_type+'_'+sequence_type+'_'+d): continue
        pilot_tsv[specimen_type+'_'+sequence_type+'_'+d] = []
    if specimen_type == 'tumor' and sequence_type == 'wgs' and workflow_type == 'bwa_alignment':
        for qc_metric in ['oxog_score', 'ContEST', 'Stars']:
            pilot_tsv[specimen_type+'_'+sequence_type+'_'+qc_metric] = []

    for d in gnos_field:
        if pilot_tsv.get(specimen_type+'_'+sequence_type+'_'+workflow_type+'_'+d): continue
        pilot_tsv[specimen_type+'_'+sequence_type+'_'+workflow_type+'_'+d] = []
    if not pilot_tsv.get('is_'+previous_release+'_'+specimen_type+'_'+sequence_type+'_'+workflow_type):
        pilot_tsv['is_'+previous_release+'_'+specimen_type+'_'+sequence_type+'_'+workflow_type] = []
    if not pilot_tsv.get(specimen_type+'_'+sequence_type+'_'+workflow_type+'_bam_file_name'):
        pilot_tsv[specimen_type+'_'+sequence_type+'_'+workflow_type+'_bam_file_name'] = []

    if not alignment:
        logger.info('Donor: {}::{} has no {} {} at {} specimen'.format(pilot_tsv.get('dcc_project_code'), pilot_tsv.get('submitter_donor_id'), sequence_type, workflow_type, specimen_type))
    elif 'normal' in specimen_type and 'wgs' in sequence_type:
        if alignment.get(workflow_type):
            generate_alignment(aliquot_field, gnos_field, alignment.get(workflow_type), pilot_tsv, specimen_type, sequence_type, workflow_type)
    elif 'normal' in specimen_type and 'rna_seq' in sequence_type:
        if alignment.get(workflow_type.replace('_alignment', '')):
            generate_alignment(aliquot_field, gnos_field, alignment.get(workflow_type.replace('_alignment', '')), pilot_tsv, specimen_type, sequence_type, workflow_type)
    elif 'tumor' in specimen_type and 'wgs' in sequence_type:
        for specimen in alignment:
            if not specimen.get(workflow_type): continue
            generate_alignment(aliquot_field, gnos_field, specimen.get(workflow_type), pilot_tsv, specimen_type, sequence_type, workflow_type)
    elif 'tumor' in specimen_type and 'rna_seq' in sequence_type:
        for specimen in alignment:
            if not specimen.get(workflow_type.replace('_alignment', '')): continue
            generate_alignment(aliquot_field, gnos_field, specimen.get(workflow_type.replace('_alignment', '')), pilot_tsv, specimen_type, sequence_type, workflow_type)                
    else:
        logger.warning('This should never happen')

    return pilot_tsv


def generate_alignment(aliquot_field, gnos_field, alignment, pilot_tsv, specimen_type, sequence_type, workflow_type):
    if not alignment: return 
    for d in aliquot_field:
        # if not alignment.get(d) in pilot_tsv[specimen_type+'_'+sequence_type+'_'+d] and sequence_type == 'rna_seq' or sequence_type == 'wgs':
        #if not alignment.get(d) in pilot_tsv[specimen_type+'_'+sequence_type+'_'+d]:
        if (sequence_type == 'rna_seq' and workflow_type == 'tophat_alignment') or (sequence_type == 'wgs' and workflow_type == 'bwa_alignment'):
            pilot_tsv[specimen_type+'_'+sequence_type+'_'+d].append(alignment.get(d)) 
    if specimen_type == 'tumor' and sequence_type == 'wgs' and workflow_type == 'bwa_alignment':
        for qc_metric in ['oxog_score', 'ContEST', 'Stars']:
            pilot_tsv[specimen_type+'_'+sequence_type+'_'+qc_metric].append(alignment.get(qc_metric))   
    for d in gnos_field:
        pilot_tsv[specimen_type+'_'+sequence_type+'_'+workflow_type+'_'+d].append(alignment.get(d))
    pilot_tsv.get('is_'+previous_release+'_'+specimen_type+'_'+sequence_type+'_'+workflow_type).append(alignment.get('is_'+previous_release+'_entry'))
    for f in alignment.get('files'):
        if f.get('file_name').endswith('.bai'): continue
        pilot_tsv[specimen_type+'_'+sequence_type+'_'+workflow_type+'_bam_file_name'].append(f.get('file_name'))                    
    return pilot_tsv


def generate_id_list(id_lists):
    ids_list = set()
    if id_lists:
        files = glob.glob(id_lists)
        for fname in files:
            with open(fname) as f:
                for d in f: 
                    if d.startswith('#'): continue
                    ids_list.add(d.rstrip())

    return ids_list


def generate_simple_release_tsv(release_donor_json, simple_release_tsv, vcf):
    simple_release = OrderedDict()
    simple_release['donor_unique_id'] = release_donor_json.get('donor_unique_id')
    if release_donor_json.get('wgs') and release_donor_json.get('wgs').get('normal_specimen'):
        for bam_type in ['bwa_alignment', 'minibam']:
            if release_donor_json.get('wgs').get('normal_specimen').get(bam_type):
                entry = release_donor_json.get('wgs').get('normal_specimen').get(bam_type)
                simple_release['gnos_id'] = entry.get('gnos_id')
                simple_release['entry_type'] = 'normal_wgs_'+bam_type
                simple_release_tsv.append(copy.deepcopy(simple_release))

    if release_donor_json.get('wgs') and release_donor_json.get('wgs').get('tumor_specimens'):
        for aliquot in release_donor_json.get('wgs').get('tumor_specimens'):
            for bam_type in ['bwa_alignment', 'minibam']:
                if aliquot.get(bam_type):
                    entry = aliquot.get(bam_type)
                    simple_release['gnos_id'] = entry.get('gnos_id')
                    simple_release['entry_type'] = 'tumor_wgs_'+bam_type
                    simple_release_tsv.append(copy.deepcopy(simple_release))
            for v in vcf:
                if aliquot.get(get_key_map(v)):
                    entry = aliquot.get(get_key_map(v))
                    simple_release['gnos_id'] = entry.get('gnos_id')
                    simple_release['entry_type'] = get_key_map(v)
                    simple_release_tsv.append(copy.deepcopy(simple_release))
    
    if release_donor_json.get('rna_seq') and release_donor_json.get('rna_seq').get('normal_specimen'):
        entry = release_donor_json.get('rna_seq').get('normal_specimen')
        for k,v in entry.iteritems():        
            simple_release['gnos_id'] = v.get('gnos_id')
            simple_release['entry_type'] = 'normal_RNA_Seq_' + k.upper() + '_bam'
            simple_release_tsv.append(copy.deepcopy(simple_release))
    
    if release_donor_json.get('rna_seq') and release_donor_json.get('rna_seq').get('tumor_specimens'):
        for aliquot in release_donor_json.get('rna_seq').get('tumor_specimens'):
            for k,v in aliquot.iteritems():
                if not v: continue        
                simple_release['gnos_id'] = v.get('gnos_id')
                simple_release['entry_type'] = 'tumor_RNA_Seq_' + k.upper() + '_bam'
                simple_release_tsv.append(copy.deepcopy(simple_release))

    simple_release_tsv = remove_dup_items(simple_release_tsv)

    return simple_release_tsv

def remove_dup_items(simple_release_tsv):
    seen = set()
    new_simple_release_tsv = []
    for d in simple_release_tsv:
        t = tuple(d.items())
        if t not in seen:
            seen.add(t)
            new_simple_release_tsv.append(d)
    return new_simple_release_tsv

def read_annotations(annotations, type, file_name):
    if not os.path.isfile(file_name):
        return
    with open(file_name, 'r') as r:

        if type == 'deprecated_gnos_id':
            if not annotations.get(type): annotations[type] = {}
            reader = csv.DictReader(r, delimiter='\t')
            for row in reader:
                donor_unique_id = row.get('dcc_project_code')+'::'+row.get('submitter_donor_id')
                if not annotations.get(type).get(donor_unique_id): annotations[type][donor_unique_id] = {}
                for vcf in ['sanger', 'dkfz', 'broad']:
                    if not row.get(vcf+'_variant_calling_deprecated_gnos_id'): continue
                    annotations[type][donor_unique_id][vcf]=row.get(vcf+'_variant_calling_deprecated_gnos_id') 

        elif type in ['blacklist', 'graylist']:
            annotations[type] = set()
            for line in r:
                if line.startswith('#'): continue
                if len(line.rstrip()) == 0: continue
                annotations[type].add(line.rstrip())


        else:
            print('unknown annotation type: {}'.format(type))
    return annotations


def main(argv=None):

    parser = ArgumentParser(description="PCAWG Reorganized Json Donors Info Generator",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=True)
    parser.add_argument("-f", "--release_file_name", dest="release_name",
             help="Specify the release name", required=True)      
    parser.add_argument("-r", "--gnos_repo", dest="repo",
             help="Specify which GNOS repo to process, process all repos if none specified", required=False)
    parser.add_argument("-x", "--exclude_donor_id_lists", dest="exclude_donor_id_lists", 
             help="File(s) containing DONOR IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-i", "--include_donor_id_lists", dest="include_donor_id_lists", 
             help="File(s) containing DONOR IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-d", "--exclude_gnos_id_lists", dest="exclude_gnos_id_lists", 
             help="File(s) containing GNOS IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-c", "--include_gnos_id_lists", dest="include_gnos_id_lists", 
             help="File(s) containing GNOS IDs to be included, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-v", "--variant_calling", dest="vcf", nargs="*",
             help="List variant_calling types", required=False)


    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports
    release_name = args.release_name
    repo = args.repo
    exclude_donor_id_lists = args.exclude_donor_id_lists
    include_donor_id_lists = args.include_donor_id_lists
    exclude_gnos_id_lists = args.exclude_gnos_id_lists
    include_gnos_id_lists = args.include_gnos_id_lists

    vcf = args.vcf
    vcf = list(vcf) if vcf else []   

    if not os.path.isdir(metadata_dir):  # TODO: should add more directory name check to make sure it's right
        sys.exit('Error: specified metadata directory does not exist!')

    timestamp = str.split(metadata_dir, '/')[-1]
    es_index = 'p_' + ('' if not repo else repo+'_') + re.sub(r'\D', '', timestamp).replace('20','',1)
    es_index_summary = 'pcawg_summary'
    es_type = "donor"
    es_host = 'localhost:9200'

    es = Elasticsearch([es_host], timeout=600)
    es_summary = init_es(es_host, es_index_summary)

    logger.setLevel(logging.INFO)
    ch.setLevel(logging.WARN)

    log_file = re.sub(r'\.py$', '.log', os.path.basename(__file__))
    # delete old log first if exists
    if os.path.isfile(log_file): os.remove(log_file)

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)

    annotations = {}
    read_annotations(annotations, 'deprecated_gnos_id', '../pcawg-operations/lists/sanger_deprecated_gnos_id.160310.tsv')   
    read_annotations(annotations, 'deprecated_gnos_id', '../pcawg-operations/lists/dkfz_embl_deprecated_gnos_id.160310.tsv')
    read_annotations(annotations, 'deprecated_gnos_id', '../pcawg-operations/lists/broad_deprecated_gnos_id.160310.tsv')
    read_annotations(annotations, 'blacklist', '../pcawg-operations/lists/blacklist/pc_annotation-donor_blacklist.tsv')
    read_annotations(annotations, 'graylist', '../pcawg-operations/lists/graylist/pc_annotation-donor_graylist.tsv')

    if not os.path.exists(metadata_dir+'/reports/'): os.makedirs(metadata_dir+'/reports/')

    # exclude the entries with gnos_ids in gnos_ids_to_be_excluded when this option is chosen
    gnos_ids_to_be_excluded = generate_id_list(exclude_gnos_id_lists)

    # include the entries with gnos_ids in gnos_ids_to_be_included when this option is chosen
    gnos_ids_to_be_included = generate_id_list(include_gnos_id_lists)

    # # remove the gnos_ids_to_be_include from gnos_ids_to_be_excluded
    # gnos_ids_to_be_excluded.difference_update(gnos_ids_to_be_included) 

    # remove the gnos_ids_to_be_excluded from gnos_ids_to_be_included
    gnos_ids_to_be_included.difference_update(gnos_ids_to_be_excluded)

    # get the list of donors to be excluded from the whitelist: blacklist/graylist donors
    if not exclude_donor_id_lists:
        # blacklist_donors = generate_id_list('../pcawg-operations/lists/blacklist/pc_annotation-donor_blacklist.tsv')
        # graylist_donors = generate_id_list('../pcawg-operations/lists/graylist/pc_annotation-donor_graylist.tsv')
        donor_ids_to_be_excluded = annotations.get('blacklist').union(annotations.get('graylist'))
    else:    
        donor_ids_to_be_excluded = generate_id_list(exclude_donor_id_lists)

    # get the full list of donors to be included in the release if not specify the donor to be included
    if not include_donor_id_lists:  
        donor_ids_to_be_included = get_donors_list(es, es_index, es_queries)
    else:
        donor_ids_to_be_included = generate_id_list(include_donor_id_lists)

    

    for dtype in ['release', 'whitelist','blacklist', 'graylist']:
        if dtype == 'release':
            donor_fh = open(metadata_dir+'/reports/'+release_name+'.jsonl', 'w') 
            pilot_tsv_fh = open(metadata_dir + '/reports/'+release_name+'.tsv', 'w') 
            simple_tsv_fh = open(metadata_dir + '/reports/'+release_name+'_entry.tsv', 'w') 
            donors_list = donor_ids_to_be_included
        elif dtype in ['blacklist', 'graylist']:
            donor_fh = open(metadata_dir+'/reports/'+release_name+'.'+dtype+'ed_donors.jsonl', 'w') 
            pilot_tsv_fh = open(metadata_dir + '/reports/'+release_name+'.'+dtype+'ed_donors.tsv', 'w') 
            simple_tsv_fh = open(metadata_dir + '/reports/'+release_name+'_entry.'+dtype+'ed_donors.tsv', 'w') 
            donors_list = annotations.get(dtype)            
        else:
            donor_fh = open(metadata_dir+'/reports/'+release_name+'.whitelisted_donors.jsonl', 'w')
            pilot_tsv_fh = open(metadata_dir + '/reports/'+release_name+'.whitelisted_donors.tsv', 'w')
            simple_tsv_fh = open(metadata_dir + '/reports/'+release_name+'_entry.whitelisted_donors.tsv', 'w')
            # exclude the donors if they were specified on the exclude_donor_id_lists
            donors_list = donor_ids_to_be_included.difference(donor_ids_to_be_excluded)

        donors_list = sorted(donors_list)  
        simple_release_tsv = []        
        # get json doc for each donor and reorganize it
        header = True 
        for donor_unique_id in donors_list:    
            es_json = get_donor_json(es, es_index, donor_unique_id)

            if not es_json: continue
            
            reorganized_donor = create_reorganized_donor(donor_unique_id, es_json, vcf, gnos_ids_to_be_excluded, gnos_ids_to_be_included, annotations)

            # push to Elasticsearch
            if dtype == 'release':
                es_summary.index(index=es_index_summary, doc_type=es_type, id=reorganized_donor['donor_unique_id'], body=json.loads(json.dumps(reorganized_donor, default=set_default)), timeout=90 )

            donor_fh.write(json.dumps(reorganized_donor, default=set_default) + '\n')

            # generate simple tsv from reorganized donor
            simple_release_tsv = generate_simple_release_tsv(reorganized_donor, simple_release_tsv, vcf)

            # generate json for tsv file from reorganized donor
            pilot_tsv_json = generate_tsv_file(reorganized_donor, vcf, annotations)
            # write to the tsv file
            if header:
                pilot_tsv_fh.write('\t'.join(pilot_tsv_json.keys()) + '\n')
                header = False 
            line = []
            for p in pilot_tsv_json.keys():
                if isinstance(pilot_tsv_json.get(p), list):
                    field = []
                    for q in pilot_tsv_json.get(p):
                        if isinstance(q, list):
                            field.append('|'.join(q))
                        elif q is None:
                            field.append('')
                        else:
                            field.append(str(q))
                    line.append(','.join(field))

                elif pilot_tsv_json.get(p) is None:
                    line.append('')
                else:
                    line.append(str(pilot_tsv_json.get(p)))

            pilot_tsv_fh.write('\t'.join(line) + '\n')
            
        pilot_tsv_fh.close()

        donor_fh.close()

        header = True 
        for r in simple_release_tsv:
            if header:
                simple_tsv_fh.write('\t'.join(r.keys()) + '\n')
                header = False 
            # make the list of output from dict
            line = []
            for p in r.keys():
                if isinstance(r.get(p), list):
                    line.append('|'.join(r.get(p)))
                elif isinstance(r.get(p), set):
                    line.append('|'.join(list(r.get(p))))
                elif r.get(p) is None:
                    line.append('')
                else:
                    line.append(str(r.get(p)))
            simple_tsv_fh.write('\t'.join(line) + '\n') 
            
        simple_tsv_fh.close() 

    return 0

if __name__ == "__main__":
    sys.exit(main())
