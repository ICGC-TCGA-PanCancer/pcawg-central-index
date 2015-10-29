#!/usr/bin/env python


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

logger = logging.getLogger('generate PCAWG data release')
ch = logging.StreamHandler()

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
                      "should": [
                        {
                           "terms":{
                              "flags.is_sanger_variant_calling_performed":[
                                 "T"
                              ]
                           }
                        },
                        {
                           "bool":{
                               "must":[
                                      {
                                       "terms":{
                                          "flags.is_dkfz_embl_variant_calling_performed":[
                                             "T"
                                          ]
                                        }
                                      },
                                      {
                                         "term":{
                                            "flags.is_dkfz_variant_calling_performed":[
                                               "F"
                                            ]
                                         }
                                      },
                                      {
                                         "term":{
                                            "flags.is_embl_variant_calling_performed":[
                                               "F"
                                            ]
                                         }
                                      } 
                                ]
                           }
                        },
                        {
                           "terms":{
                              "flags.is_broad_variant_calling_performed":[
                                 "T"
                              ]
                           }
                        }
                      ],
                      "must_not": [
                        {
                          "terms": {
                            "flags.is_bam_used_by_variant_calling_missing": [
                              "T"
                            ]
                          }
                        },
                        {
                          "terms": {
                            "flags.exists_vcf_file_prefix_mismatch": [
                              "T"
                            ]
                          }
                        },
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
                        },
                        {
                          "terms": {
                            "flags.is_donor_blacklisted": [
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

def create_reorganized_donor(donor_unique_id, es_json, vcf, gnos_ids_to_be_excluded):
    reorganized_donor = {
        'donor_unique_id': donor_unique_id,
        'submitter_donor_id': es_json['submitter_donor_id'],
        'dcc_project_code': es_json['dcc_project_code'],
        'icgc_donor_id': es_json['icgc_donor_id'],
        'aug2015_donor': True if es_json.get('flags').get('is_aug2015_donor') else False,
        'santa_cruz_pilot': True if es_json.get('flags').get('is_santa_cruz_donor') else False,
        'validation_by_deep_seq': True if es_json.get('flags').get('is_train2_pilot') else False,
        'wgs': {
            'normal_specimen': {},
            'tumor_specimens': []
        },
        'rna_seq': {
             'normal_specimen': {},
             'tumor_specimens': []
        }
    }
    add_wgs_specimens(reorganized_donor, es_json, vcf, gnos_ids_to_be_excluded)
    add_rna_seq_info(reorganized_donor, es_json, gnos_ids_to_be_excluded)

    return reorganized_donor


def create_alignment(es_json, aliquot, data_type, gnos_ids_to_be_excluded):
    if not aliquot.get('aligned_bam'): return
    gnos_id = aliquot.get('aligned_bam').get('gnos_id')
    if gnos_ids_to_be_excluded and gnos_id in gnos_ids_to_be_excluded: return
    alignment = {
        'submitter_specimen_id': aliquot.get('submitter_specimen_id'),
        'icgc_specimen_id': aliquot.get('icgc_specimen_id'),
        'submitter_sample_id': aliquot.get('submitter_sample_id'),
        'icgc_sample_id': aliquot.get('icgc_sample_id'),
        'specimen_type': aliquot.get('dcc_specimen_type'),
        'aliquot_id': aliquot.get('aliquot_id'),
        'is_aug2015_entry': aliquot.get('aligned_bam').get('is_aug2015_entry') if 'wgs' in data_type else aliquot.get('is_aug2015_entry'),
        'gnos_repo': filter_liri_jp(es_json.get('dcc_project_code'), \
            aliquot.get('aligned_bam').get('gnos_repo'), \
            data_type, aliquot.get('aliquot_id')),
        'gnos_id': aliquot.get('aligned_bam').get('gnos_id'),
        'gnos_last_modified': aliquot.get('aligned_bam').get('gnos_last_modified')[-1],
        'files': []
    }
    # add the file info if exist
    for ftype in ['bam', 'bai']:
        if aliquot.get('aligned_bam').get(ftype+'_file_name'):
            f = {
                'file_name': aliquot.get('aligned_bam').get(ftype+'_file_name'),
                'file_md5sum': aliquot.get('aligned_bam').get(ftype+'_file_md5sum'),
                'file_size': aliquot.get('aligned_bam').get(ftype+'_file_size'),
            }
            alignment.get('files').append(f)
        else:
            logger.warning('{} alignment GNOS entry {} has no .{} file'.format(data_type, alignment.get('gnos_id'), ftype))
    return alignment


def get_formal_vcf_name(vcf):
    vcf_map = {
      "sanger": "sanger_variant_calling",
      "dkfz": "dkfz_embl_variant_calling",
      "embl": "dkfz_embl_variant_calling",
      "broad": "broad_variant_calling",
      "muse": "muse_variant_calling",
      "broad_tar": "broad_tar_variant_calling"
    }   

    return vcf_map.get(vcf)


def choose_variant_calling(es_json, vcf):
    variant_calling = set()
    if not es_json.get('variant_calling_results') or not vcf:
        return variant_calling

    for v in vcf:
        if get_formal_vcf_name(v) in es_json.get('variant_calling_results').keys() and \
            not es_json.get('variant_calling_results').get(get_formal_vcf_name(v)).get('is_stub'):
            variant_calling.add(get_formal_vcf_name(v))
            if not check_broad_vcf(es_json, v): variant_calling.discard(get_formal_vcf_name(v))
        else:
            logger.warning('donor: {} has no {}'.format(es_json.get('donor_unique_id'), get_formal_vcf_name(v)))
    return variant_calling


def check_broad_vcf(es_json, vcf_calling):
    if vcf_calling == 'broad' or vcf_calling == 'muse' or vcf_calling == 'broad_tar':
        if not es_json.get('flags').get('is_broad_variant_calling_performed'):
            return False
        else: 
            return True
    else:
        return True
        

def create_variant_calling(es_json, aliquot, wgs_tumor_vcf_info, data_type, gnos_ids_to_be_excluded):
    gnos_id = wgs_tumor_vcf_info.get('gnos_id')
    if gnos_ids_to_be_excluded and gnos_id in gnos_ids_to_be_excluded: return
    variant_calling = {
        'submitter_specimen_id': aliquot.get('submitter_specimen_id'),
        'icgc_specimen_id': aliquot.get('icgc_specimen_id'),
        'submitter_sample_id': aliquot.get('submitter_sample_id'),
        'icgc_sample_id': aliquot.get('icgc_sample_id'),
        'specimen_type': aliquot.get('dcc_specimen_type'),
        'aliquot_id': aliquot.get('aliquot_id'),
        'is_aug2015_entry': wgs_tumor_vcf_info.get('is_aug2015_entry'),
        'gnos_repo': wgs_tumor_vcf_info.get('gnos_repo'),
        'gnos_id': wgs_tumor_vcf_info.get('gnos_id'),
        'gnos_last_modified': wgs_tumor_vcf_info.get('gnos_last_modified')[-1],
        'files':[]
    }
    vcf_files = wgs_tumor_vcf_info.get('files')
    if vcf_files:
        for f in vcf_files:
            if variant_calling.get('aliquot_id') in f.get('file_name'):
                variant_calling.get('files').append(f)   
    else:
        logger.warning('{} GNOS entry {} has no files'.format(data_type, variant_calling.get('gnos_id')))

    return variant_calling


def add_wgs_specimens(reorganized_donor, es_json, vcf, gnos_ids_to_be_excluded):
    if es_json.get('normal_alignment_status'):
        wgs_normal_alignment_info = es_json.get('normal_alignment_status')
        data_type = 'wgs_normal_bwa'
        reorganized_donor.get('wgs').get('normal_specimen')['bwa_alignment'] = create_alignment(es_json, wgs_normal_alignment_info, data_type, gnos_ids_to_be_excluded)

    tumor_wgs_specimen_count = 0
    if es_json.get('tumor_alignment_status'):
        variant_calling = choose_variant_calling(es_json, vcf)
        wgs_tumor_alignment_info = es_json.get('tumor_alignment_status')
    
        for aliquot in wgs_tumor_alignment_info:
            aliquot_info = OrderedDict()
            tumor_wgs_specimen_count += 1
            data_type = 'wgs_tumor_bwa'
            aliquot_info['bwa_alignment'] = create_alignment(es_json, aliquot, data_type, gnos_ids_to_be_excluded)        
            
            for vc in variant_calling:
                wgs_tumor_vcf_info = es_json.get('variant_calling_results').get(vc)
                data_type = 'wgs_'+vc        
                aliquot_info[vc] = create_variant_calling(es_json, aliquot, wgs_tumor_vcf_info, data_type, gnos_ids_to_be_excluded)        
            reorganized_donor.get('wgs').get('tumor_specimens').append(copy.deepcopy(aliquot_info)) 

    reorganized_donor['tumor_wgs_specimen_count'] = tumor_wgs_specimen_count
    return reorganized_donor


def filter_liri_jp(project, gnos_repo, data_type, aliquot_id):
    if not project == 'LIRI-JP' or 'rna_seq' in data_type:
        return gnos_repo
    elif "https://gtrepo-riken.annailabs.com/" in gnos_repo:
        return ["https://gtrepo-riken.annailabs.com/"]
    else:
        print "This should never happen: alignment for LIRI-JP is not available at Riken repo. Alignment type: {}, aliquot_id: {}".format(data_type, aliquot_id)
        return [ gnos_repo[0] ]  # return the first one, not an entirely proper solution but gets us going


def add_rna_seq_info(reorganized_donor, es_json, gnos_ids_to_be_excluded):
    rna_seq_info = es_json.get('rna_seq').get('alignment')
    for specimen_type in rna_seq_info.keys():
        if not rna_seq_info.get(specimen_type): # the specimen_type has no alignment result
		    continue
        
        if 'normal' in specimen_type:
            aliquot = rna_seq_info.get(specimen_type)
            alignment_info = {}
            for workflow_type in aliquot.keys():
                data_type = 'rna_seq_normal_'+workflow_type
                alignment_info[workflow_type] = create_alignment(es_json, aliquot.get(workflow_type), data_type, gnos_ids_to_be_excluded)

            reorganized_donor.get('rna_seq')[specimen_type + '_specimen'] = alignment_info
        else:
            tumor_rna_seq_specimen_count = 0
            for aliquot in rna_seq_info.get(specimen_type):
                tumor_rna_seq_specimen_count += 1
                alignment_info = {}
                for workflow_type in aliquot.keys():
                    data_type = 'rna_seq_tumor_'+workflow_type
                    alignment_info[workflow_type] = create_alignment(es_json, aliquot.get(workflow_type), data_type, gnos_ids_to_be_excluded)
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
    es_mapping = open('pancan.reorganized.donor.mapping.json')
    es.indices.put_mapping(index=es_index, doc_type='donor', body=es_mapping.read())
    es_mapping.close()
    return es


def set_default(obj):
    if isinstance(obj, datetime.datetime):
        return obj.isoformat()
    if isinstance(obj, set):
        return list(obj)
    raise TypeError


def generate_tsv_file(reorganized_donor, vcf):
    donor_info = ['donor_unique_id','dcc_project_code', 'submitter_donor_id', 'icgc_donor_id', 'aug2015_donor','santa_cruz_pilot', 'validation_by_deep_seq']
    specimen = ['submitter_specimen_id', 'icgc_specimen_id', 'submitter_sample_id', 'icgc_sample_id', 'aliquot_id']
    alignment = ['alignment_gnos_repo', 'alignment_gnos_id', 'alignment_bam_file_name']
        
    pilot_tsv = OrderedDict()
    for d in donor_info:
        pilot_tsv[d] = reorganized_donor.get(d)
    #wgs normal specimen
    alignment = reorganized_donor.get('wgs').get('normal_specimen')
    generate_alignment_info(pilot_tsv, alignment, 'normal', 'wgs', 'alignment')
    # wgs tumor specimen
    wgs_tumor_speciments = reorganized_donor.get('wgs').get('tumor_specimens')
    pilot_tsv['tumor_wgs_specimen_count'] = reorganized_donor.get('tumor_wgs_specimen_count')
    generate_alignment_info(pilot_tsv, wgs_tumor_speciments, 'tumor', 'wgs', 'alignment')
    # wgs variant calling
    generate_variant_calling_info(pilot_tsv, wgs_tumor_speciments, vcf)
    # rna_seq normal
    for workflow in ['star', 'tophat']:
        alignment = reorganized_donor.get('rna_seq').get('normal_specimen')
        generate_alignment_info(pilot_tsv, alignment, 'normal', 'rna_seq', workflow.lower()+'_alignment')
    # rna_seq tumor
    rna_seq_tumor = reorganized_donor.get('rna_seq').get('tumor_specimens')
    pilot_tsv['tumor_rna_seq_specimen_count'] = reorganized_donor.get('tumor_rna_seq_specimen_count')  
    for workflow in ['star', 'tophat']:
        generate_alignment_info(pilot_tsv, rna_seq_tumor, 'tumor', 'rna_seq', workflow.lower()+'_alignment')
        
    return pilot_tsv


def generate_variant_calling_info(pilot_tsv, variant_calling, vcf):
    for v in vcf:       
        pilot_tsv[get_formal_vcf_name(v)+'_repo'] = []
        pilot_tsv[get_formal_vcf_name(v)+'_gnos_id'] = []
        pilot_tsv[get_formal_vcf_name(v)+'_file_name_prefix'] = []
        pilot_tsv['is_aug2015_'+get_formal_vcf_name(v)] = []
        for specimen in variant_calling:
            if specimen.get(get_formal_vcf_name(v)):
                pilot_tsv[get_formal_vcf_name(v)+'_repo'] = specimen.get(get_formal_vcf_name(v)).get('gnos_repo')
                pilot_tsv[get_formal_vcf_name(v)+'_gnos_id'] = specimen.get(get_formal_vcf_name(v)).get('gnos_id')
                pilot_tsv[get_formal_vcf_name(v)+'_file_name_prefix'].append(specimen.get(get_formal_vcf_name(v)).get('aliquot_id'))
                pilot_tsv['is_aug2015_'+get_formal_vcf_name(v)] = specimen.get(get_formal_vcf_name(v)).get('is_aug2015_entry')
    return pilot_tsv


def generate_alignment_info(pilot_tsv, alignment, specimen_type, sequence_type, workflow_type):
    aliquot_field = ['submitter_specimen_id', 'icgc_specimen_id', 'submitter_sample_id', 'icgc_sample_id', 'aliquot_id']
    gnos_field = ['gnos_repo', 'gnos_id']
    for d in aliquot_field:
        if pilot_tsv.get(specimen_type+'_'+sequence_type+'_'+d): continue
        pilot_tsv[specimen_type+'_'+sequence_type+'_'+d] = []
    for d in gnos_field:
        if pilot_tsv.get(specimen_type+'_'+sequence_type+'_'+workflow_type+'_'+d): continue
        pilot_tsv[specimen_type+'_'+sequence_type+'_'+workflow_type+'_'+d] = []
    if not pilot_tsv.get('is_aug2015_'+specimen_type+'_'+sequence_type+'_'+workflow_type):
    	pilot_tsv['is_aug2015_'+specimen_type+'_'+sequence_type+'_'+workflow_type] = []
    if not pilot_tsv.get(specimen_type+'_'+sequence_type+'_'+workflow_type+'_bam_file_name'):
        pilot_tsv[specimen_type+'_'+sequence_type+'_'+workflow_type+'_bam_file_name'] = []

    if not alignment:
        logger.info('Donor: {}::{} has no {} {} at {} specimen'.format(pilot_tsv.get('dcc_project_code'), pilot_tsv.get('submitter_donor_id'), sequence_type, workflow_type, specimen_type))
    elif 'normal' in specimen_type and 'wgs' in sequence_type:
        generate_alignment(aliquot_field, gnos_field, alignment.get('bwa_alignment'), pilot_tsv, specimen_type, sequence_type, workflow_type)
    elif 'normal' in specimen_type and 'rna_seq' in sequence_type:
        if alignment.get(workflow_type.replace('_alignment', '')):
            generate_alignment(aliquot_field, gnos_field, alignment.get(workflow_type.replace('_alignment', '')), pilot_tsv, specimen_type, sequence_type, workflow_type)
    elif 'tumor' in specimen_type and 'wgs' in sequence_type:
        for specimen in alignment:
            generate_alignment(aliquot_field, gnos_field, specimen.get('bwa_alignment'), pilot_tsv, specimen_type, sequence_type, workflow_type)
    elif 'tumor' in specimen_type and 'rna_seq' in sequence_type:
        for specimen in alignment:
            if specimen.get(workflow_type.replace('_alignment', '')):
                generate_alignment(aliquot_field, gnos_field, specimen.get(workflow_type.replace('_alignment', '')), pilot_tsv, specimen_type, sequence_type, workflow_type)                
    else:
        logger.warning('This should never happen')

    return pilot_tsv


def generate_alignment(aliquot_field, gnos_field, alignment, pilot_tsv, specimen_type, sequence_type, workflow_type):
    if not alignment: return 
    for d in aliquot_field:
        if not alignment.get(d) in pilot_tsv[specimen_type+'_'+sequence_type+'_'+d]:
            pilot_tsv[specimen_type+'_'+sequence_type+'_'+d].append(alignment.get(d)) 
    for d in gnos_field:
        pilot_tsv[specimen_type+'_'+sequence_type+'_'+workflow_type+'_'+d].append(alignment.get(d))
    pilot_tsv.get('is_aug2015_'+specimen_type+'_'+sequence_type+'_'+workflow_type).append(alignment.get('is_aug2015_entry'))
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
                for d in f: ids_list.add(d.rstrip())

    return ids_list


def generate_simple_release_tsv(release_donor_json, simple_release_tsv, vcf):
    simple_release = OrderedDict()
    simple_release['donor_unique_id'] = release_donor_json.get('donor_unique_id')
    if release_donor_json.get('wgs') and release_donor_json.get('wgs').get('normal_specimen') and \
        release_donor_json.get('wgs').get('normal_specimen').get('bwa_alignment'):
        entry = release_donor_json.get('wgs').get('normal_specimen').get('bwa_alignment')
        simple_release['gnos_id'] = entry.get('gnos_id')
        simple_release['entry_type'] = 'normal_wgs_bwa_bam'
        simple_release_tsv.append(copy.deepcopy(simple_release))

    if release_donor_json.get('wgs') and release_donor_json.get('wgs').get('tumor_specimens'):
        for aliquot in release_donor_json.get('wgs').get('tumor_specimens'):
            if aliquot.get('bwa_alignment'):
                entry = aliquot.get('bwa_alignment')
                simple_release['gnos_id'] = entry.get('gnos_id')
                simple_release['entry_type'] = 'tumor_wgs_bwa_bam'
                simple_release_tsv.append(copy.deepcopy(simple_release))
            for v in vcf:
                if aliquot.get(get_formal_vcf_name(v)):
                    entry = aliquot.get(get_formal_vcf_name(v))
                    simple_release['gnos_id'] = entry.get('gnos_id')
                    simple_release['entry_type'] = get_formal_vcf_name(v)
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
    parser.add_argument("-v", "--variant_calling", dest="vcf", nargs="*",
             help="List variant_calling types", required=False)


    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports
    release_name = args.release_name
    repo = args.repo
    exclude_donor_id_lists = args.exclude_donor_id_lists
    include_donor_id_lists = args.include_donor_id_lists
    exclude_gnos_id_lists = args.exclude_gnos_id_lists

    vcf = args.vcf

    vcf = list(vcf) if vcf else []   

    if not os.path.isdir(metadata_dir):  # TODO: should add more directory name check to make sure it's right
        sys.exit('Error: specified metadata directory does not exist!')

    timestamp = str.split(metadata_dir, '/')[-1]
    es_index = 'p_' + ('' if not repo else repo+'_') + re.sub(r'\D', '', timestamp).replace('20','',1)
    es_index_reorganize = 'r_' + ('' if not repo else repo+'_') + re.sub(r'\D', '', timestamp).replace('20','',1)
    es_type = "donor"
    es_host = 'localhost:9200'

    es = Elasticsearch([es_host])

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

    if not os.path.exists(metadata_dir+'/reports/'): os.makedirs(metadata_dir+'/reports/')

    donor_fh = open(metadata_dir+'/reports/'+release_name+'.jsonl', 'w')
    pilot_tsv_fh = open(metadata_dir + '/reports/'+release_name+'.tsv', 'w')  
    simple_tsv_fh = open(metadata_dir + '/reports/'+release_name+'_entry.tsv', 'w')

    # pre-exclude donors when this option is chosen
    donor_ids_to_be_excluded = generate_id_list(exclude_donor_id_lists)

	# get the list of donors
    # only process the gnos entries when this option is chosen
    donor_ids_to_be_included = generate_id_list(include_donor_id_lists)
    if not donor_ids_to_be_included:  
        donors_list = get_donors_list(es, es_index, es_queries)
    else:
        donors_list = donor_ids_to_be_included

    # exclude the donors if they were specified on the exclude_donor_id_lists
    donors_list.difference_update(donor_ids_to_be_excluded)

    # exclude the entries with gnos_ids in gnos_ids_to_be_excluded when this option is chosen
    gnos_ids_to_be_excluded = generate_id_list(exclude_gnos_id_lists)

    donors_list = sorted(donors_list)  
    simple_release_tsv = []

    # get json doc for each donor and reorganize it
    header = True 
    for donor_unique_id in donors_list:    
        
    	es_json = get_donor_json(es, es_index, donor_unique_id)
        
        reorganized_donor = create_reorganized_donor(donor_unique_id, es_json, vcf, gnos_ids_to_be_excluded)

        donor_fh.write(json.dumps(reorganized_donor, default=set_default) + '\n')

        # generate simple tsv from reorganized donor
        simple_release_tsv = generate_simple_release_tsv(reorganized_donor, simple_release_tsv, vcf)

        # generate json for tsv file from reorganized donor
        pilot_tsv_json = generate_tsv_file(reorganized_donor, vcf)
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
