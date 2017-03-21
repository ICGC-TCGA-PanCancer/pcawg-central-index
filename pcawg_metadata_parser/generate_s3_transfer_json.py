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
from elasticsearch1 import Elasticsearch
from collections import OrderedDict
import datetime
import dateutil.parser
from itertools import izip
from distutils.version import LooseVersion
import hashlib
import xml.dom.minidom
import shutil
import requests

id_service_token = os.environ.get('ICGC_TOKEN')

logger = logging.getLogger('s3 transfer json generator')
ch = logging.StreamHandler()

es_queries = [
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
              "terms": {
                "dcc_project_code": [
                    "LIRI-JP",
                    "PACA-CA",
                    "PRAD-CA",
                    "RECA-EU",
                    "PAEN-AU",
                    "PACA-AU",
                    "BOCA-UK",
                    "OV-AU",
                    "MELA-AU",
                    "BRCA-UK"
                    "PRAD-UK",
                    "CMDI-UK",
                    "LINC-JP",
                    "ORCA-IN",
                    "BTCA-SG",
                    "LAML-KR",
                    "LICA-FR",
                    "CLLE-ES",
                    "ESAD-UK"
                ]
              }
            },
            # {
            #   "terms":{
            #     "donor_unique_id":[
            #       "PACA-AU::ICGC_0088"
            #     ]
            #   }
            # },           
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
            },
            {
              "terms":{
                "flags.is_oct2015_donor":[
                  "T"
                ]
              }
            },
            # {
            #   "range":{
            #     "flags.all_tumor_specimen_aliquot_counts":{"gte": 2}
            #   }
            # },
            {
              "range":{
                    "qc_score":{"gte": 0, "lt": 10000}
                }
            }
          ],
          "must_not": [
            {
              "regexp": {
                "dcc_project_code": ".*-US"
              }
            },
            {
              "regexp": {
                "dcc_project_code": ".*-DE"
              }
            },
            {
              "terms": {
                "flags.is_bam_used_by_variant_calling_missing": [
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
            #   "terms": {
            #     "flags.exists_xml_md5sum_mismatch": [
            #       "T"
            #     ]
            #   }
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


def get_source_repo_index_pos (available_repos, chosen_gnos_repo=None):
    source_repo_rank = [
        "https://gtrepo-osdc-icgc.annailabs.com/",
        "https://gtrepo-ebi.annailabs.com/",
        "https://gtrepo-riken.annailabs.com/",
        "https://gtrepo-etri.annailabs.com/",
        "https://gtrepo-bsc.annailabs.com/",
        "https://gtrepo-dkfz.annailabs.com/"
    ]
    if chosen_gnos_repo and get_formal_repo_name(chosen_gnos_repo) in available_repos:
        source_repo_rank = [ get_formal_repo_name(chosen_gnos_repo) ]
    for r in source_repo_rank:
        try:
            if available_repos.index(r) >= 0: return available_repos.index(r)
        except:
            pass

    logger.warning('Source repo not allowed to be transferred to S3')
    return None


def get_formal_repo_name(repo):
    repo_url_to_repo = {
      "https://gtrepo-bsc.annailabs.com/": "bsc",
      "bsc": "https://gtrepo-bsc.annailabs.com/",
      "https://gtrepo-ebi.annailabs.com/": "ebi",
      "ebi": "https://gtrepo-ebi.annailabs.com/",
      "https://cghub.ucsc.edu/": "cghub",
      "cghub": "https://cghub.ucsc.edu/",
      "https://gtrepo-dkfz.annailabs.com/": "dkfz",
      "dkfz": "https://gtrepo-dkfz.annailabs.com/",
      "https://gtrepo-riken.annailabs.com/": "riken",
      "riken": "https://gtrepo-riken.annailabs.com/",
      "https://gtrepo-osdc-icgc.annailabs.com/": "osdc-icgc",
      "osdc-icgc": "https://gtrepo-osdc-icgc.annailabs.com/",
      "https://gtrepo-osdc-tcga.annailabs.com/": "osdc-tcga",
      "osdc-tcga": "https://gtrepo-osdc-tcga.annailabs.com/",
      "https://gtrepo-etri.annailabs.com/": "etri",
      "etri": "https://gtrepo-etri.annailabs.com/"
    }

    return repo_url_to_repo.get(repo)


def generate_md5_size(metadata_xml_file):
    with open (metadata_xml_file, 'r') as x: data = x.read()
    data = re.sub(r'<ResultSet .+?>', '<ResultSet>', data)

    with open('tmp.xml', 'w') as f: f.write(data)

    xml_md5 = hashlib.md5(data).hexdigest()
    xml_size = os.path.getsize('tmp.xml')

    return [xml_md5, xml_size]


def generate_object_id(filename, gnos_id):
    global id_service_token
    url = 'https://meta.icgc.org/entities'
    # try get request first
    r = requests.get(url + '?gnosId=' + gnos_id + '&fileName=' + filename,
                       headers={'Content-Type': 'application/json'})
    if not r or not r.ok:
        logger.warning('GET request unable to access metadata service: {}'.format(url))
        return ''
    elif r.json().get('totalElements') == 1:
        logger.info('GET request got the id')
        return r.json().get('content')[0].get('id')
    elif r.json().get('totalElements') > 1:
        logger.warning('GET request to metadata service return multiple matches for gnos_id: {} and filename: {}'
                          .format(gnos_id, filename))
        return ''
    elif id_service_token:  # no match then try post to create
        headers = {
            'Content-Type': 'application/json',
            'Authorization': 'Bearer ' + id_service_token
        }
        body = {
            "gnosId": gnos_id,
            "fileName": filename
        }
        r = requests.post(url, data=json.dumps(body), headers=headers)
        if not r or not r.ok:
            logger.warning('POST request failed')
            return ''
        return r.json().get('id')
    else:
        logger.info('No luck, generate FAKE ID')
        return ''

def add_wgs_normal_specimen(es_json, gnos_ids_to_be_included, gnos_ids_to_be_excluded, chosen_gnos_repo, jobs_dir):
    aliquot = es_json.get('normal_alignment_status')
    gnos_id = aliquot.get('aligned_bam').get('gnos_id')
    if gnos_ids_to_be_included and not gnos_id in gnos_ids_to_be_included: return
    if gnos_ids_to_be_excluded and gnos_id in gnos_ids_to_be_excluded: return

    aliquot_info = create_bwa_alignment(aliquot, es_json, chosen_gnos_repo)
    write_s3_transfer_json(jobs_dir, aliquot_info, gnos_ids_to_be_excluded)


def add_metadata_xml_info(obj, chosen_gnos_repo=None):
    repo = get_formal_repo_name(obj.get('gnos_repo')[ get_source_repo_index_pos(obj.get('gnos_repo'), chosen_gnos_repo) ])
    gnos_id = obj.get('gnos_id')
    ao_state = 'live'
    ao_updated = obj.get('gnos_last_modified')[ get_source_repo_index_pos(obj.get('gnos_repo'), chosen_gnos_repo) ].encode('utf8')
    ao_updated = str.split(ao_updated, '+')[0] + 'Z'
    metadata_xml_file = 'gnos_metadata/__all_metadata_xml/' + repo + '/' + gnos_id + '__' + ao_state + '__' + ao_updated + '.xml'
    metadata_xml_file_info = {
        'file_name': gnos_id + '.xml',
        'file_md5sum': generate_md5_size(metadata_xml_file)[0],
        'file_size': generate_md5_size(metadata_xml_file)[1],
        'object_id': generate_object_id(gnos_id+'.xml', gnos_id)
    }

    return metadata_xml_file_info


def get_available_repos(obj):
    repos = obj.get('gnos_repo')
    ret_repos = []
    for r in repos:
        metadata_xml_info = add_metadata_xml_info(obj, get_formal_repo_name(r))
        ret_repos.append({
              r:{
                  'file_md5sum': metadata_xml_info.get('file_md5sum'),
                  'file_size': metadata_xml_info.get('file_size')
                }
            })
    return ret_repos


def create_bwa_alignment(aliquot, es_json, chosen_gnos_repo):
    aliquot_info = {
        'data_type': 'WGS-BWA-Normal' if 'normal' in aliquot.get('dcc_specimen_type').lower() else 'WGS-BWA-Tumor',
        'project_code': es_json['dcc_project_code'],
        'submitter_donor_id': es_json['submitter_donor_id'],
        'is_santa_cruz': aliquot.get('aligned_bam').get('is_santa_cruz_entry'),
        'submitter_specimen_id': aliquot.get('submitter_specimen_id'),
        'submitter_sample_id': aliquot.get('submitter_sample_id'),
        'specimen_type': aliquot.get('dcc_specimen_type'),
        'aliquot_id': aliquot.get('aliquot_id'),
        'available_repos': get_available_repos(aliquot.get('aligned_bam')),
        'gnos_repo': [ aliquot.get('aligned_bam').get('gnos_repo')[ \
            get_source_repo_index_pos(aliquot.get('aligned_bam').get('gnos_repo'), chosen_gnos_repo) ] ],
        'gnos_id': aliquot.get('aligned_bam').get('gnos_id'),
        'files': [
            {
                'file_name': aliquot.get('aligned_bam').get('bam_file_name'),
                'file_md5sum': aliquot.get('aligned_bam').get('bam_file_md5sum'),
                'file_size': aliquot.get('aligned_bam').get('bam_file_size'),
                'object_id': generate_object_id(aliquot.get('aligned_bam').get('bam_file_name'), aliquot.get('aligned_bam').get('gnos_id'))                       
            }
        ]
    }

    # add the bai file info if exist
    if aliquot.get('aligned_bam').get('bai_file_name'):
        bai_file = {
            'file_name': aliquot.get('aligned_bam').get('bai_file_name'),
            'file_md5sum': aliquot.get('aligned_bam').get('bai_file_md5sum'),
            'file_size': aliquot.get('aligned_bam').get('bai_file_size'),
            'object_id': generate_object_id(aliquot.get('aligned_bam').get('bai_file_name'), aliquot.get('aligned_bam').get('gnos_id'))                        
        }
        aliquot_info.get('files').append(bai_file)
    else:
        logger.warning('BWA alignment GNOS entry {} has no .bai file'.format(aliquot_info.get('gnos_id')))

    # add the metadata_xml_file_info
    metadata_xml_file_info = add_metadata_xml_info(aliquot.get('aligned_bam'), chosen_gnos_repo)
    aliquot_info.get('files').append(metadata_xml_file_info)   

    return aliquot_info


def add_wgs_tumor_specimens(es_json, gnos_ids_to_be_included, gnos_ids_to_be_excluded, chosen_gnos_repo, jobs_dir):
    wgs_tumor_alignment_info = es_json.get('tumor_alignment_status')

    for aliquot in wgs_tumor_alignment_info:
        gnos_id = aliquot.get('aligned_bam').get('gnos_id')
        if gnos_ids_to_be_included and not gnos_id in gnos_ids_to_be_included: continue
        if gnos_ids_to_be_excluded and gnos_id in gnos_ids_to_be_excluded: continue

        aliquot_info = create_bwa_alignment(aliquot, es_json, chosen_gnos_repo)
        write_s3_transfer_json(jobs_dir, aliquot_info, gnos_ids_to_be_excluded)


def add_variant_calling(es_json, gnos_ids_to_be_included, gnos_ids_to_be_excluded, chosen_gnos_repo, jobs_dir, vcf, vcf_result_version):
    if not es_json.get('variant_calling_results'): return
    
    variant_callings = choose_variant_calling(es_json, vcf)
    for v in variant_callings:

        if not es_json.get('variant_calling_results').get(v): continue
        if not es_json.get('variant_calling_results').get(v).get('vcf_workflow_result_version') == vcf_result_version: continue

        wgs_tumor_vcf_info = es_json.get('variant_calling_results').get(v)

        gnos_id = wgs_tumor_vcf_info.get('gnos_id')
        if gnos_ids_to_be_included and not gnos_id in gnos_ids_to_be_included: continue
        if gnos_ids_to_be_excluded and gnos_id in gnos_ids_to_be_excluded: continue

        variant_calling = {
            'data_type': get_formal_vcf_name(v).capitalize()+'-VCF',
            'project_code': es_json['dcc_project_code'],
            'submitter_donor_id': es_json['submitter_donor_id'],  
            'vcf_workflow_result_version': wgs_tumor_vcf_info.get('vcf_workflow_result_version'),             
            'submitter_specimen_id': None,
            'submitter_sample_id': None,
            'specimen_type': None,
            'aliquot_id': None,
            'available_repos': get_available_repos(wgs_tumor_vcf_info),
            'gnos_repo': [ wgs_tumor_vcf_info.get('gnos_repo')[ \
                get_source_repo_index_pos(wgs_tumor_vcf_info.get('gnos_repo'), chosen_gnos_repo) ] ],
            'gnos_id': wgs_tumor_vcf_info.get('gnos_id'),
            'files': wgs_tumor_vcf_info.get('files')
        }
        
        # add the object_id for each file object
        for f in variant_calling.get('files'):
            if int(f.get('file_size')) == 0: 
                logger.warning('donor: {} has variant_calling file: {} file_size is 0'.format(es_json.get('donor_unique_id'), f.get('file_name')))
                variant_calling.get('files').remove(f)
            f.update({'file_size': None if f.get('file_size') == None else int(f.get('file_size'))})
            f.update({'object_id': generate_object_id(f.get('file_name'), variant_calling.get('gnos_id'))})

        # add the metadata_xml_file_info
        metadata_xml_file_info = add_metadata_xml_info(wgs_tumor_vcf_info, chosen_gnos_repo)

        variant_calling.get('files').append(metadata_xml_file_info) 

        write_s3_transfer_json(jobs_dir, variant_calling, gnos_ids_to_be_excluded)           


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


def get_formal_vcf_name(vcf):
    vcf_map = {
      "sanger": "sanger_variant_calling",
      "dkfz": "dkfz_embl_variant_calling",
      "embl": "dkfz_embl_variant_calling",
      "dkfz_embl": "dkfz_embl_variant_calling",
      "broad": "broad_variant_calling",
      "muse": "muse_variant_calling",
      "broad_tar": "broad_tar_variant_calling",
      "sanger_variant_calling": "sanger",
      "dkfz_embl_variant_calling": "dkfz_embl",
      "broad_variant_calling": "broad",
      "muse_variant_calling": "muse",
      "broad_tar_variant_calling": "broad_tar"
    }   

    return vcf_map.get(vcf)


def add_rna_seq_info(es_json, gnos_ids_to_be_included, gnos_ids_to_be_excluded, chosen_gnos_repo, jobs_dir):
    # to build pcawg santa cruz pilot dataset, this is a temporary walkaround to exclude the 130 RNA-Seq bad
    # entries from MALY-DE and CLLE-ES projects
    #if reorganized_donor.get('dcc_project_code') in ('MALY-DE', 'CLLE-ES'): return

    rna_seq_info = es_json.get('rna_seq').get('alignment')
    for specimen_type in rna_seq_info.keys():
        if not rna_seq_info.get(specimen_type): # the specimen_type has no alignment result
		    continue
        if 'normal' in specimen_type:
            aliquot = rna_seq_info.get(specimen_type)
            alignment_info = {}
            for workflow_type in aliquot.keys():
                gnos_id = aliquot.get(workflow_type).get('aligned_bam').get('gnos_id')
                if gnos_ids_to_be_included and not gnos_id in gnos_ids_to_be_included: continue
                if gnos_ids_to_be_excluded and gnos_id in gnos_ids_to_be_excluded: continue

                alignment_info[workflow_type] = create_rna_seq_alignment(aliquot, es_json, workflow_type, chosen_gnos_repo)
                write_s3_transfer_json(jobs_dir, alignment_info[workflow_type], gnos_ids_to_be_excluded)

        else:
            for aliquot in rna_seq_info.get(specimen_type):
                alignment_info = {}
                for workflow_type in aliquot.keys():
                    gnos_id = aliquot.get(workflow_type).get('aligned_bam').get('gnos_id')
                    if gnos_ids_to_be_included and not gnos_id in gnos_ids_to_be_included: continue
                    if gnos_ids_to_be_excluded and gnos_id in gnos_ids_to_be_excluded: continue

                    alignment_info[workflow_type] = create_rna_seq_alignment(aliquot, es_json, workflow_type, chosen_gnos_repo)
                    write_s3_transfer_json(jobs_dir, alignment_info[workflow_type], gnos_ids_to_be_excluded)


def create_rna_seq_alignment(aliquot, es_json, workflow_type, chosen_gnos_repo):
    alignment_info = {
        'data_type': 'RNA_Seq-'+workflow_type.capitalize()+'-Normal' if 'normal' in aliquot.get(workflow_type).get('dcc_specimen_type').lower() else 'RNA_Seq-'+workflow_type.capitalize()+'-Tumor',
        'project_code': es_json['dcc_project_code'],
        'submitter_donor_id': es_json['submitter_donor_id'],
        'is_santa_cruz': aliquot.get(workflow_type).get('is_santa_cruz_entry'),
        'submitter_specimen_id': aliquot.get(workflow_type).get('submitter_specimen_id'),
        'submitter_sample_id': aliquot.get(workflow_type).get('submitter_sample_id'),
        'specimen_type': aliquot.get(workflow_type).get('dcc_specimen_type'),
        'aliquot_id': aliquot.get(workflow_type).get('aliquot_id'),
        'available_repos': get_available_repos(aliquot.get(workflow_type).get('aligned_bam')),
        'gnos_repo': [ aliquot.get(workflow_type).get('aligned_bam').get('gnos_repo')[ \
            get_source_repo_index_pos(aliquot.get(workflow_type).get('aligned_bam').get('gnos_repo'), chosen_gnos_repo)] ],
        'gnos_id': aliquot.get(workflow_type).get('aligned_bam').get('gnos_id'),
        'files': [
            {
                'file_name': aliquot.get(workflow_type).get('aligned_bam').get('bam_file_name'),
                'file_md5sum': aliquot.get(workflow_type).get('aligned_bam').get('bam_file_md5sum'),
                'file_size': aliquot.get(workflow_type).get('aligned_bam').get('bam_file_size'),
                'object_id': generate_object_id(aliquot.get(workflow_type).get('aligned_bam').get('bam_file_name'), aliquot.get(workflow_type).get('aligned_bam').get('gnos_id'))                          
            }
        ]
    }

    # add the bai file info if exist
    if aliquot.get(workflow_type).get('aligned_bam').get('bai_file_name'):
        bai_file = {
            'file_name': aliquot.get(workflow_type).get('aligned_bam').get('bai_file_name'),
            'file_md5sum': aliquot.get(workflow_type).get('aligned_bam').get('bai_file_md5sum'),
            'file_size': aliquot.get(workflow_type).get('aligned_bam').get('bai_file_size'),
            'object_id': generate_object_id(aliquot.get(workflow_type).get('aligned_bam').get('bai_file_name'), aliquot.get(workflow_type).get('aligned_bam').get('gnos_id'))                        
        }
        alignment_info.get('files').append(bai_file)
    else:
        logger.warning('RNA-Seq alignment GNOS entry {} has no .bai file'.format(alignment_info.get('gnos_id')))

    # add the metadata_xml_file_info
    metadata_xml_file_info = add_metadata_xml_info(aliquot.get(workflow_type).get('aligned_bam'), chosen_gnos_repo)
    alignment_info.get('files').append(metadata_xml_file_info)

    return alignment_info


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


def set_default(obj):
    if isinstance(obj, datetime.datetime):
        return obj.isoformat()
    if isinstance(obj, set):
        return list(obj)
    raise TypeError


def write_s3_transfer_json(jobs_dir, transfer_json, gnos_ids_to_be_excluded):
    if transfer_json:
        gnos_id = transfer_json.get('gnos_id')

        project_code = transfer_json.get('project_code')
        donor_id = transfer_json.get('submitter_donor_id')
        specimen_id = '-' if transfer_json.get('data_type').endswith('-VCF') else transfer_json.get('submitter_specimen_id')
        data_type = transfer_json.get('data_type')
        
        json_name_list = [gnos_id, project_code, donor_id, specimen_id, data_type, 'json']
        sub_json_name = '.'.join(json_name_list[1:])

        if gnos_ids_to_be_excluded and sub_json_name in gnos_ids_to_be_excluded:
            logger.warning('{} is already in s3 scheduled'.format(sub_json_name)) 
            return

        json_name = '.'.join(json_name_list)
        json_name_new = json_name.replace(' ', '__')
        with open(jobs_dir + '/' + json_name_new, 'w') as w:
            w.write(json.dumps(transfer_json, indent=4, sort_keys=True))


def generate_id_list(id_lists):
    ids_list = set()
    if id_lists:
        files = glob.glob(id_lists)
        for fname in files:
            with open(fname) as f:
                for d in f: ids_list.add(d.rstrip())

    return ids_list


def main(argv=None):

    parser = ArgumentParser(description="Transfer Jobs Json Generator",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=True)
    parser.add_argument("-t", "--target_cloud", dest="target_cloud",
             help="Specify target_cloud of the job transfer", required=True)
    parser.add_argument("-r", "--specify source repo", dest="chosen_gnos_repo",
             help="Specify source gnos repo", required=False)
    parser.add_argument("-d", "--exclude_donor_id_lists", dest="exclude_donor_id_lists", 
             help="File(s) containing DONOR IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-c", "--include_donor_id_lists", dest="include_donor_id_lists", 
             help="File(s) containing DONOR IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-i", "--include_gnos_id_lists", dest="include_gnos_id_lists",
             help="Specify which GNOS IDs to process, process all gnos_ids if none specified", required=False)
    parser.add_argument("-x", "--exclude_gnos_id_lists", dest="exclude_gnos_id_lists", 
             help="File(s) containing GNOS IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-s", "--sequence_type", dest="seq", nargs="*",
             help="List sequence_type types", required=False)
    parser.add_argument("-v", "--variant_calling", dest="vcf", nargs="*",
             help="List variant_calling types", required=False) 
    parser.add_argument("-n", "--specify the vcf_workflow_result_version", dest="vcf_result_version", default="v2", type=str,
             help="Specify vcf_workflow_result_version", required=False)   



    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports
    target_cloud = args.target_cloud
    include_gnos_id_lists = args.include_gnos_id_lists
    exclude_gnos_id_lists = args.exclude_gnos_id_lists
    exclude_donor_id_lists = args.exclude_donor_id_lists
    include_donor_id_lists = args.include_donor_id_lists
    chosen_gnos_repo = args.chosen_gnos_repo
    seq = args.seq
    vcf = args.vcf
    vcf_result_version = args.vcf_result_version


    seq= list(seq) if seq else [] 
    vcf = list(vcf) if vcf else []   

    # pre-exclude gnos entries when this option is chosen
    gnos_ids_to_be_excluded = generate_id_list(exclude_gnos_id_lists)

    # read and parse git for the gnos_ids and fnames which are scheduled for s3 transfer
    if target_cloud == 'aws':
        git_s3_fnames = '../s3-transfer-operations/s3-transfer-jobs*/*/*.json'
    elif target_cloud == 'collab':
        git_s3_fnames = '../ceph_transfer_ops/ceph-transfer-jobs*/*/*.json'
    else:
        sys.exit('Error: unknown target_cloud!')
    files = glob.glob(git_s3_fnames)
    for f in files:
        fname = str.split(f, '/')[-1]
        gnos_id = str.split(fname, '.')[0]
        gnos_ids_to_be_excluded.add(gnos_id)
        sub_file_name = '.'.join(str.split(fname, '.')[1:])
        gnos_ids_to_be_excluded.add(sub_file_name)

    # only process the gnos entries when this option is chosen
    gnos_ids_to_be_included = generate_id_list(include_gnos_id_lists) 

    # # remove the gnos_ids_to_be_include from gnos_ids_to_be_excluded
    # gnos_ids_to_be_excluded.difference_update(gnos_ids_to_be_included) 

    # remove the gnos_ids_to_be_excluded from gnos_ids_to_be_include
    gnos_ids_to_be_included.difference_update(gnos_ids_to_be_excluded)

    if not os.path.isdir(metadata_dir):  # TODO: should add more directory name check to make sure it's right
        sys.exit('Error: specified metadata directory does not exist!')

    timestamp = str.split(metadata_dir, '/')[-1]
    es_index = 'p_' + re.sub(r'\D', '', timestamp).replace('20','',1)
    es_type = "donor"
    es_host = 'localhost:9200'

    es = Elasticsearch([es_host])

    # pre-exclude donors when this option is chosen
    donor_ids_to_be_excluded = generate_id_list(exclude_donor_id_lists)
    donor_ids_to_be_included = generate_id_list(include_donor_id_lists)
    if not donor_ids_to_be_included:  
        donors_list = get_donors_list(es, es_index, es_queries)
    else:
        donors_list = donor_ids_to_be_included

    # exclude the donors if they were specified on the exclude_donor_id_lists
    donors_list.difference_update(donor_ids_to_be_excluded)

    report_dir = re.sub(r'^generate_', '', os.path.basename(__file__))
    report_dir = re.sub(r'\.py$', '', report_dir)
    jobs_dir = metadata_dir + '/reports/' + report_dir

    if os.path.exists(jobs_dir): shutil.rmtree(jobs_dir, ignore_errors=True)  # empty the folder if exists
    os.makedirs(jobs_dir)

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

    
    # get json doc for each donor 
    for donor_unique_id in donors_list:     
        
    	es_json = get_donor_json(es, es_index, donor_unique_id)

        if seq and 'wgs' in seq:
            add_wgs_normal_specimen(es_json, gnos_ids_to_be_included, gnos_ids_to_be_excluded, chosen_gnos_repo, jobs_dir)

            add_wgs_tumor_specimens(es_json, gnos_ids_to_be_included, gnos_ids_to_be_excluded, chosen_gnos_repo, jobs_dir)

        if vcf:
            add_variant_calling(es_json, gnos_ids_to_be_included, gnos_ids_to_be_excluded, chosen_gnos_repo, jobs_dir, vcf, vcf_result_version)

        if seq and 'rna_seq' in seq:
            add_rna_seq_info(reorganized_donor, es_json, gnos_ids_to_be_included, gnos_ids_to_be_excluded, chosen_gnos_repo, jobs_dir)

    if os.path.isfile('tmp.xml'): os.remove('tmp.xml')

    return 0

if __name__ == "__main__":
    sys.exit(main())
