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
                            "flags.is_sanger_variant_calling_performed":[
                              "T"
                            ]
                          }
                        },
                        {
                          "terms": {
                            "variant_calling_results.sanger_variant_calling.is_bam_used_by_sanger_missing": [
                              "F"
                            ]
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
                      "must_not": [
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

def create_reorganized_donor(donor_unique_id, es_json, vcf):
    reorganized_donor = {
        'donor_unique_id': donor_unique_id,
        'submitter_donor_id': es_json['submitter_donor_id'],
        'dcc_project_code': es_json['dcc_project_code'],
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
    add_wgs_specimens(reorganized_donor, es_json, vcf)
    add_rna_seq_info(reorganized_donor, es_json)

    return reorganized_donor


def create_alignment(es_json, aliquot, data_type):
    alignment = {
        'submitter_specimen_id': aliquot.get('submitter_specimen_id'),
        'submitter_sample_id': aliquot.get('submitter_sample_id'),
        'specimen_type': aliquot.get('dcc_specimen_type'),
        'aliquot_id': aliquot.get('aliquot_id'),
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
      "broad": "broad_variant_calling"
    }   

    return vcf_map.get(vcf)


def choose_variant_calling(es_json, vcf):
    variant_calling = []
    if not vcf:
        variant_calling = es_json.get('variant_calling_results').keys()
    else:
        for v in vcf:
            if get_formal_vcf_name(v) in es_json.get('variant_calling_results').keys() and \
                not es_json.get('variant_calling_results').get(get_formal_vcf_name(v)).get('is_stub'):
                variant_calling.append(get_formal_vcf_name(v))
            else:
                logger.warning('donor: {} has no {}'.format(es_json.get('donor_unique_id'), get_formal_vcf_name(v)))
    return variant_calling


def create_variant_calling(es_json, aliquot, wgs_tumor_vcf_info, data_type):
    variant_calling = {
        'submitter_specimen_id': aliquot.get('submitter_specimen_id'),
        'submitter_sample_id': aliquot.get('submitter_sample_id'),
        'specimen_type': aliquot.get('dcc_specimen_type'),
        'aliquot_id': aliquot.get('aliquot_id'),
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


def add_wgs_specimens(reorganized_donor, es_json, vcf):
    wgs_normal_alignment_info = es_json.get('normal_alignment_status')
    data_type = 'wgs_normal_bwa'
    reorganized_donor.get('wgs').get('normal_specimen')['bwa_alignment'] = create_alignment(es_json, wgs_normal_alignment_info, data_type)

    variant_calling = choose_variant_calling(es_json, vcf)
    wgs_tumor_alignment_info = es_json.get('tumor_alignment_status')

    tumor_wgs_specimen_count = 0
    for aliquot in wgs_tumor_alignment_info:
        aliquot_info = OrderedDict()
        tumor_wgs_specimen_count += 1
        data_type = 'wgs_tumor_bwa'
        aliquot_info['bwa_alignment'] = create_alignment(es_json, aliquot, data_type)        

        for vc in variant_calling:
            wgs_tumor_vcf_info = es_json.get('variant_calling_results').get(vc)
            data_type = 'wgs_'+vc        
            aliquot_info[vc] = create_variant_calling(es_json, aliquot, wgs_tumor_vcf_info, data_type)        
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


def add_rna_seq_info(reorganized_donor, es_json):
    rna_seq_info = es_json.get('rna_seq').get('alignment')
    for specimen_type in rna_seq_info.keys():
        if not rna_seq_info.get(specimen_type): # the specimen_type has no alignment result
		    continue
        if 'normal' in specimen_type:
            aliquot = rna_seq_info.get(specimen_type)
            alignment_info = {}
            for workflow_type in aliquot.keys():
                data_type = 'rna_seq_normal_'+workflow_type
                alignment_info[workflow_type] = create_alignment(es_json, aliquot.get(workflow_type), data_type)

            reorganized_donor.get('rna_seq')[specimen_type + '_specimen'] = alignment_info
        else:
            for aliquot in rna_seq_info.get(specimen_type):
                alignment_info = {}
                for workflow_type in aliquot.keys():
                    data_type = 'rna_seq_tumor_'+workflow_type
                    alignment_info[workflow_type] = create_alignment(es_json, aliquot.get(workflow_type), data_type)
                reorganized_donor.get('rna_seq')[specimen_type + '_specimens'].append(copy.deepcopy(alignment_info)) 


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

def generate_json_for_tsv_file(reorganized_donor):
    pilot_tsv_json = OrderedDict()
    pilot_tsv_json['dcc_project_code'] = reorganized_donor.get('dcc_project_code')
    pilot_tsv_json['submitter_donor_id'] = reorganized_donor.get('submitter_donor_id')
    pilot_tsv_json['santa_cruz_pilot'] = reorganized_donor.get('santa_cruz_pilot')
    pilot_tsv_json['validation_by_deep_seq'] = reorganized_donor.get('validation_by_deep_seq')
    # wgs normal specimen 
    pilot_tsv_json['normal_wgs_submitter_specimen_id'] = reorganized_donor.get('wgs').get('normal_specimen').get('bwa_alignment').get('submitter_specimen_id')
    pilot_tsv_json['normal_wgs_submitter_sample_id'] = reorganized_donor.get('wgs').get('normal_specimen').get('bwa_alignment').get('submitter_sample_id')
    pilot_tsv_json['normal_wgs_aliquot_id'] = reorganized_donor.get('wgs').get('normal_specimen').get('bwa_alignment').get('aliquot_id')
    pilot_tsv_json['normal_wgs_alignment_gnos_repo'] = [reorganized_donor.get('wgs').get('normal_specimen').get('bwa_alignment').get('gnos_repo')]
    pilot_tsv_json['normal_wgs_alignment_gnos_id'] = reorganized_donor.get('wgs').get('normal_specimen').get('bwa_alignment').get('gnos_id')
    for f in reorganized_donor.get('wgs').get('normal_specimen').get('bwa_alignment').get('files'):
        if 'bai' in f.get('file_name'): continue
        pilot_tsv_json['normal_wgs_alignment_bam_file_name'] = f.get('file_name')
    # wgs tumor specimen
    wgs_tumor_speciments = reorganized_donor.get('wgs').get('tumor_specimens')
    pilot_tsv_json['tumor_wgs_specimen_count'] = reorganized_donor.get('tumor_wgs_specimen_count')
    pilot_tsv_json['tumor_wgs_submitter_specimen_id'] = [] 
    pilot_tsv_json['tumor_wgs_submitter_sample_id'] = []
    pilot_tsv_json['tumor_wgs_aliquot_id'] = []
    pilot_tsv_json['tumor_wgs_alignment_gnos_repo'] = []
    pilot_tsv_json['tumor_wgs_alignment_gnos_id'] = []
    pilot_tsv_json['tumor_wgs_alignment_bam_file_name'] = []
    # wgs tumor sanger vcf
    pilot_tsv_json['sanger_variant_calling_repo'] = []
    pilot_tsv_json['sanger_variant_calling_gnos_id'] = wgs_tumor_speciments[0].get('sanger_variant_calling').get('gnos_id')
    pilot_tsv_json['sanger_variant_calling_file_name_prefix'] = []
    for specimen in wgs_tumor_speciments:
        pilot_tsv_json['tumor_wgs_submitter_specimen_id'].append(specimen.get('bwa_alignment').get('submitter_specimen_id'))
        pilot_tsv_json['tumor_wgs_submitter_sample_id'].append(specimen.get('bwa_alignment').get('submitter_sample_id'))
        pilot_tsv_json['tumor_wgs_aliquot_id'].append(specimen.get('bwa_alignment').get('aliquot_id'))
        pilot_tsv_json['tumor_wgs_alignment_gnos_repo'].append(specimen.get('bwa_alignment').get('gnos_repo'))
        pilot_tsv_json['tumor_wgs_alignment_gnos_id'].append(specimen.get('bwa_alignment').get('gnos_id'))
        for f in specimen.get('bwa_alignment').get('files'):
            if 'bai' in f.get('file_name'): continue
            pilot_tsv_json['tumor_wgs_alignment_bam_file_name'].append(f.get('file_name'))
        # wgs tumor sanger vcf
        pilot_tsv_json['sanger_variant_calling_repo'].append(specimen.get('sanger_variant_calling').get('gnos_repo'))
        pilot_tsv_json['sanger_variant_calling_file_name_prefix'].append(specimen.get('sanger_variant_calling').get('aliquot_id'))
    
    # rna_seq normal specimen
    pilot_tsv_json['normal_rna_seq_submitter_specimen_id'] = None
    pilot_tsv_json['normal_rna_seq_submitter_sample_id'] = None
    pilot_tsv_json['normal_rna_seq_aliquot_id'] = None
    pilot_tsv_json['normal_rna_seq_STAR_alignment_gnos_repo'] = None
    pilot_tsv_json['normal_rna_seq_STAR_alignment_gnos_id'] = None
    pilot_tsv_json['normal_rna_seq_STAR_alignment_bam_file_name'] = None
    pilot_tsv_json['normal_rna_seq_TOPHAT2_alignment_gnos_repo'] = None
    pilot_tsv_json['normal_rna_seq_TOPHAT2_alignment_gnos_id'] = None
    pilot_tsv_json['normal_rna_seq_TOPHAT2_alignment_bam_file_name'] = None

    rna_seq_normal = reorganized_donor.get('rna_seq').get('normal_specimen')
    if rna_seq_normal and rna_seq_normal.get('tophat'):
        pilot_tsv_json['normal_rna_seq_submitter_specimen_id'] = rna_seq_normal.get('tophat').get('submitter_specimen_id')
        pilot_tsv_json['normal_rna_seq_submitter_sample_id'] = rna_seq_normal.get('tophat').get('submitter_sample_id')
        pilot_tsv_json['normal_rna_seq_aliquot_id'] = rna_seq_normal.get('tophat').get('aliquot_id')
        pilot_tsv_json['normal_rna_seq_TOPHAT2_alignment_gnos_repo'] = [rna_seq_normal.get('tophat').get('gnos_repo')]
        pilot_tsv_json['normal_rna_seq_TOPHAT2_alignment_gnos_id'] = rna_seq_normal.get('tophat').get('gnos_id')
        for f in rna_seq_normal.get('tophat').get('files'):
            if 'bai' in f.get('file_name'): continue
            pilot_tsv_json['normal_rna_seq_TOPHAT2_alignment_bam_file_name'] = f.get('file_name')
    if rna_seq_normal and rna_seq_normal.get('star'):
        pilot_tsv_json['normal_rna_seq_submitter_specimen_id'] = rna_seq_normal.get('star').get('submitter_specimen_id')
        pilot_tsv_json['normal_rna_seq_submitter_sample_id'] = rna_seq_normal.get('star').get('submitter_sample_id')
        pilot_tsv_json['normal_rna_seq_aliquot_id'] = rna_seq_normal.get('star').get('aliquot_id')
        pilot_tsv_json['normal_rna_seq_STAR_alignment_gnos_repo'] = rna_seq_normal.get('star').get('gnos_repo')
        pilot_tsv_json['normal_rna_seq_STAR_alignment_gnos_id'] = rna_seq_normal.get('star').get('gnos_id')
        for f in rna_seq_normal.get('star').get('files'):
            if 'bai' in f.get('file_name'): continue
            pilot_tsv_json['normal_rna_seq_STAR_alignment_bam_file_name'] = f.get('file_name')
       
    # rna_seq tumor specimens
    pilot_tsv_json['tumor_rna_seq_submitter_specimen_id'] = []
    pilot_tsv_json['tumor_rna_seq_submitter_sample_id'] = []
    pilot_tsv_json['tumor_rna_seq_aliquot_id'] = []
    pilot_tsv_json['tumor_rna_seq_STAR_alignment_gnos_repo'] = []
    pilot_tsv_json['tumor_rna_seq_STAR_alignment_gnos_id'] = []
    pilot_tsv_json['tumor_rna_seq_STAR_alignment_bam_file_name'] = []
    pilot_tsv_json['tumor_rna_seq_TOPHAT2_alignment_gnos_repo'] = []
    pilot_tsv_json['tumor_rna_seq_TOPHAT2_alignment_gnos_id'] = []
    pilot_tsv_json['tumor_rna_seq_TOPHAT2_alignment_bam_file_name'] = []

    rna_seq_tumor = reorganized_donor.get('rna_seq').get('tumor_specimens')
    rna_seq_tumor_specimen_id = []
    rna_seq_tumor_sample_id = []
    rna_seq_tumor_aliquot_id = []
    if rna_seq_tumor:
        for rna_seq_tumor_specimen in rna_seq_tumor:
            if rna_seq_tumor_specimen.get('tophat'):
                rna_seq_tumor_specimen_id_tmp = rna_seq_tumor_specimen.get('tophat').get('submitter_specimen_id')
                rna_seq_tumor_sample_id_tmp = rna_seq_tumor_specimen.get('tophat').get('submitter_sample_id')
                rna_seq_tumor_aliquot_id_tmp = rna_seq_tumor_specimen.get('tophat').get('aliquot_id')
                pilot_tsv_json['tumor_rna_seq_TOPHAT2_alignment_gnos_repo'].append(rna_seq_tumor_specimen.get('tophat').get('gnos_repo'))
                pilot_tsv_json['tumor_rna_seq_TOPHAT2_alignment_gnos_id'].append(rna_seq_tumor_specimen.get('tophat').get('gnos_id'))
                for f in rna_seq_tumor_specimen.get('tophat').get('files'):
                    if 'bai' in f.get('file_name'): continue
                    pilot_tsv_json['tumor_rna_seq_TOPHAT2_alignment_bam_file_name'].append(f.get('file_name'))
            if rna_seq_tumor_specimen.get('star'):
                rna_seq_tumor_specimen_id_tmp = rna_seq_tumor_specimen.get('star').get('submitter_specimen_id')
                rna_seq_tumor_sample_id_tmp = rna_seq_tumor_specimen.get('star').get('submitter_sample_id')
                rna_seq_tumor_aliquot_id_tmp = rna_seq_tumor_specimen.get('star').get('aliquot_id')
                pilot_tsv_json['tumor_rna_seq_STAR_alignment_gnos_repo'].append(rna_seq_tumor_specimen.get('star').get('gnos_repo'))
                pilot_tsv_json['tumor_rna_seq_STAR_alignment_gnos_id'].append(rna_seq_tumor_specimen.get('star').get('gnos_id'))
                for f in rna_seq_tumor_specimen.get('star').get('files'):
                    if 'bai' in f.get('file_name'): continue
                    pilot_tsv_json['tumor_rna_seq_STAR_alignment_bam_file_name'].append(f.get('file_name'))
            rna_seq_tumor_specimen_id.append(rna_seq_tumor_specimen_id_tmp)
            rna_seq_tumor_sample_id.append(rna_seq_tumor_sample_id_tmp)
            rna_seq_tumor_aliquot_id.append(rna_seq_tumor_aliquot_id_tmp)
        pilot_tsv_json['tumor_rna_seq_submitter_specimen_id'] = rna_seq_tumor_specimen_id
        pilot_tsv_json['tumor_rna_seq_submitter_sample_id'] = rna_seq_tumor_sample_id
        pilot_tsv_json['tumor_rna_seq_aliquot_id'] = rna_seq_tumor_aliquot_id
    
    return pilot_tsv_json


def generate_id_list(id_lists):
    ids_list = set()
    if id_lists:
        files = glob.glob(id_lists)
        for fname in files:
            with open(fname) as f:
                for d in f: ids_list.add(d.rstrip())

    return ids_list

def generate_simple_release_tsv(release_donor_json, simple_release_tsv):
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
            if aliquot.get('sanger_variant_calling'):
                entry = aliquot.get('sanger_variant_calling')
                simple_release['gnos_id'] = entry.get('gnos_id')
                simple_release['entry_type'] = 'sanger_vcf'
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
                simple_release['gnos_id'] = v.get('gnos_id')
                simple_release['entry_type'] = 'tumor_RNA_Seq_' + k.upper() + '_bam'
                simple_release_tsv.append(copy.deepcopy(simple_release))

    return simple_release_tsv


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
    parser.add_argument("-v", "--variant_calling", dest="vcf", nargs="*",
             help="List variant_calling types", required=False)


    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports
    release_name = args.release_name
    repo = args.repo
    exclude_donor_id_lists = args.exclude_donor_id_lists
    include_donor_id_lists = args.include_donor_id_lists
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

    donors_list = sorted(donors_list)  
    simple_release_tsv = []

    # get json doc for each donor and reorganize it
    header = True 
    for donor_unique_id in donors_list:    
        
    	es_json = get_donor_json(es, es_index, donor_unique_id)
        
        reorganized_donor = create_reorganized_donor(donor_unique_id, es_json, vcf)

        donor_fh.write(json.dumps(reorganized_donor, default=set_default) + '\n')

        # generate simple tsv from reorganized donor
        simple_release_tsv = generate_simple_release_tsv(reorganized_donor, simple_release_tsv)

        # generate json for tsv file from reorganized donor
        pilot_tsv_json = generate_json_for_tsv_file(reorganized_donor)
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