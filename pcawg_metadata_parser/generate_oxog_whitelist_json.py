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
import hashlib
import xml.dom.minidom
import shutil
import requests
import csv

# id_service_token = os.environ.get('ICGC_TOKEN')
icgc_project_code = os.environ.get('ICGC_PROJECT_CODE')

logger = logging.getLogger('OxOG filter json generator')
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
                "dcc_project_code": [icgc_project_code]
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
            },
            # {
            #   "range":{
            #         "flags.all_tumor_specimen_aliquot_counts":{"lt": 2}
            #     }
            # }, 
            {
               "terms":{
                  "flags.is_sanger_variant_calling_performed":[
                     "T"
                  ]
               }
            },
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
            },
            {
               "terms":{
                  "flags.is_broad_variant_calling_performed":[
                     "T"
                  ]
               }
            },
            {
              "terms":{
                "flags.is_mar2016_donor":[
                  "T"
                ]
              }
            },
            # {
            #   "range":{
            #         "qc_score":{"gte": 0, "lt": 10000}
            #     }
            # }
          ],
          "must_not": [
            # {
            #   "regexp": {
            #     "dcc_project_code": ".*-US"
            #   }
            # },
            # {
            #   "regexp": {
            #     "dcc_project_code": ".*-DE"
            #   }
            # },
            {
              "terms": {
                "flags.is_bam_used_by_variant_calling_missing": [
                  "T"
                ]
              }
            },
            {
               "terms":{
                  "flags.exists_vcf_file_prefix_mismatch":[
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


def get_source_repo_index_pos (available_repos, chosen_gnos_repo=None):
    source_repo_rank = [
        "https://gtrepo-bsc.annailabs.com/",
        "https://gtrepo-dkfz.annailabs.com/",
        "https://gtrepo-osdc-icgc.annailabs.com/",
        "https://gtrepo-ebi.annailabs.com/",
        "https://gtrepo-riken.annailabs.com/",
        "https://gtrepo-etri.annailabs.com/",
        "https://gtrepo-osdc-tcga.annailabs.com/",
        "https://cghub.ucsc.edu/"
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


def generate_object_id(filename, gnos_id, target_cloud):
    if target_cloud in ['tcga', 'gtdownload']: return ''
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


def add_metadata_xml_info(obj, chosen_gnos_repo, target_cloud):
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
        'object_id': generate_object_id(gnos_id+'.xml', gnos_id, target_cloud)
    }

    return metadata_xml_file_info


def get_available_repos(obj, target_cloud):
    repos = obj.get('gnos_repo')
    ret_repos = []
    for r in repos:
        metadata_xml_info = add_metadata_xml_info(obj, get_formal_repo_name(r), target_cloud)
        ret_repos.append({
              r:{
                  'file_md5sum': metadata_xml_info.get('file_md5sum'),
                  'file_size': metadata_xml_info.get('file_size')
                }
            })
    return ret_repos


def create_bwa_alignment(aliquot, es_json, chosen_gnos_repo, oxog_score, target_cloud):
    aliquot_info = {
        'data_type': 'WGS-BWA-Normal' if 'normal' in aliquot.get('dcc_specimen_type').lower() else 'WGS-BWA-Tumor',
        'oxog_score': oxog_score.get(aliquot.get('aliquot_id')) if oxog_score.get(aliquot.get('aliquot_id')) else None,
        'submitter_specimen_id': aliquot.get('submitter_specimen_id'),
        'submitter_sample_id': aliquot.get('submitter_sample_id'),
        'specimen_type': aliquot.get('dcc_specimen_type'),
        'aliquot_id': aliquot.get('aliquot_id'),
        'available_repos': get_available_repos(aliquot.get('aligned_bam'), target_cloud),
        'gnos_repo': [ aliquot.get('aligned_bam').get('gnos_repo')[ \
            get_source_repo_index_pos(aliquot.get('aligned_bam').get('gnos_repo'), chosen_gnos_repo) ] ],
        'gnos_id': aliquot.get('aligned_bam').get('gnos_id'),
        'files': [
            {
                'file_name': aliquot.get('aligned_bam').get('bam_file_name'),
                'file_md5sum': aliquot.get('aligned_bam').get('bam_file_md5sum'),
                'file_size': aliquot.get('aligned_bam').get('bam_file_size'),
                'object_id': generate_object_id(aliquot.get('aligned_bam').get('bam_file_name'), aliquot.get('aligned_bam').get('gnos_id'), target_cloud)                       
            }
        ]
    }

    # add the bai file info if exist
    if aliquot.get('aligned_bam').get('bai_file_name'):
        bai_file = {
            'file_name': aliquot.get('aligned_bam').get('bai_file_name'),
            'file_md5sum': aliquot.get('aligned_bam').get('bai_file_md5sum'),
            'file_size': aliquot.get('aligned_bam').get('bai_file_size'),
            'object_id': generate_object_id(aliquot.get('aligned_bam').get('bai_file_name'), aliquot.get('aligned_bam').get('gnos_id'), target_cloud)                        
        }
        aliquot_info.get('files').append(bai_file)
    else:
        logger.warning('BWA alignment GNOS entry {} has no .bai file'.format(aliquot_info.get('gnos_id')))

    # add the metadata_xml_file_info
    metadata_xml_file_info = add_metadata_xml_info(aliquot.get('aligned_bam'), chosen_gnos_repo, target_cloud)
    aliquot_info.get('files').append(metadata_xml_file_info)   

    return aliquot_info


def add_wgs_specimens(es_json, chosen_gnos_repo, jobs_dir, job_json, oxog_score, gnos_ids_in_cloud, target_cloud):
    if not es_json.get('normal_alignment_status') or not es_json.get('tumor_alignment_status'):
        logger.warning('The donor {} has no normal or tumor alignments.'.format(es_json.get('donor_unique_id')))
        return False

    # add tumor
    wgs_tumor_alignment_info = es_json.get('tumor_alignment_status')
    for aliquot in wgs_tumor_alignment_info:
        # if not oxog_score.get(aliquot.get('aliquot_id')):
        if not oxog_score.get(aliquot.get('aliquot_id')): # or float(oxog_score.get(aliquot.get('aliquot_id'))) >= 40:
            logger.warning('The aliquot {} of donor {} has no oxog_score.'.format(aliquot.get('aliquot_id'), es_json.get('donor_unique_id'))) 
            return False
        
        if not aliquot.get('aligned_bam'):
            logger.warning('The donor {} has no aligned_bam.'.format(es_json.get('donor_unique_id'))) 
            return False

        if target_cloud in ['aws', 'collab'] and not aliquot.get('aligned_bam').get('gnos_id') in gnos_ids_in_cloud:
            logger.warning('The donor {} has NOT transferred tumor bam.'.format(es_json.get('donor_unique_id'))) 
            return False             

        aliquot_info = create_bwa_alignment(aliquot, es_json, chosen_gnos_repo, oxog_score, target_cloud)
        job_json.get('tumors').append(copy.deepcopy(aliquot_info))

    # add normal
    aliquot = es_json.get('normal_alignment_status')
    if not aliquot.get('aligned_bam'):
        logger.warning('The donor {} has no aligned_bam.'.format(es_json.get('donor_unique_id'))) 
        return False

    if target_cloud in ['aws', 'collab'] and not aliquot.get('aligned_bam').get('gnos_id') in gnos_ids_in_cloud:
        logger.warning('The donor {} has NOT transferred normal bam.'.format(es_json.get('donor_unique_id'))) 
        return False 
    job_json['normal'] = create_bwa_alignment(aliquot, es_json, chosen_gnos_repo, oxog_score, target_cloud)

    return job_json
    

def add_variant_calling(es_json, chosen_gnos_repo, jobs_dir, job_json, gnos_ids_in_cloud, target_cloud):
    if not es_json.get('variant_calling_results'): return False
  
    for v in ['sanger', 'dkfz_embl', 'broad', 'muse']:

        if not es_json.get('variant_calling_results').get(get_formal_vcf_name(v)):
            logger.warning('donor: {} has no {}'.format(es_json.get('donor_unique_id'), get_formal_vcf_name(v)))
            return False

        wgs_tumor_vcf_info = es_json.get('variant_calling_results').get(get_formal_vcf_name(v))

        gnos_id = wgs_tumor_vcf_info.get('gnos_id')
        if target_cloud in ['aws', 'collab'] and not gnos_id in gnos_ids_in_cloud:
            logger.warning('The donor {} has NOT transferred {}_variant_calling.'.format(es_json.get('donor_unique_id'), v)) 
            return False 

        variant_calling = {
            'data_type': get_formal_vcf_name(v).capitalize()+'-VCF',          
            'available_repos': get_available_repos(wgs_tumor_vcf_info, target_cloud),
            'gnos_repo': [ wgs_tumor_vcf_info.get('gnos_repo')[ \
                get_source_repo_index_pos(wgs_tumor_vcf_info.get('gnos_repo'), chosen_gnos_repo) ] ],
            'gnos_id': wgs_tumor_vcf_info.get('gnos_id'),
            'files': [],
            'vcf_workflow_result_version': wgs_tumor_vcf_info.get('vcf_workflow_result_version')
        }
        
        # add the object_id for each file object
        vcf_files = wgs_tumor_vcf_info.get('files')
        for f in vcf_files:
            if int(f.get('file_size')) == 0: 
                logger.warning('donor: {} has variant_calling file: {} file_size is 0'.format(es_json.get('donor_unique_id'), f.get('file_name')))
                continue
            f.update({'file_size': None if f.get('file_size') == None else int(f.get('file_size'))})
            f.update({'object_id': generate_object_id(f.get('file_name'), variant_calling.get('gnos_id'), target_cloud)})
            variant_calling.get('files').append(f)

        # add the metadata_xml_file_info
        metadata_xml_file_info = add_metadata_xml_info(wgs_tumor_vcf_info, chosen_gnos_repo, target_cloud)

        variant_calling.get('files').append(metadata_xml_file_info)

        job_json.update({v: copy.deepcopy(variant_calling)})       
        
    return job_json


def choose_variant_calling(es_json, vcf):
    variant_calling = set()
    if not es_json.get('variant_calling_results') or not vcf:
        return variant_calling

    for v in vcf:
        if get_formal_vcf_name(v) in es_json.get('variant_calling_results').keys() and \
            not es_json.get('variant_calling_results').get(get_formal_vcf_name(v)).get('is_stub'):
            variant_calling.add(get_formal_vcf_name(v))
            if not check_vcf(es_json, v): 
                variant_calling.discard(get_formal_vcf_name(v))

        else:
            logger.warning('donor: {} has no {}'.format(es_json.get('donor_unique_id'), get_formal_vcf_name(v)))
    return variant_calling


def check_vcf(es_json, vcf_calling):
    if vcf_calling == 'broad' or vcf_calling == 'muse' or vcf_calling == 'broad_tar':
        if not es_json.get('flags').get('is_broad_variant_calling_performed'):
            return False
        elif not es_json.get('variant_calling_results').get(get_formal_vcf_name('broad')).get('vcf_workflow_result_version') == 'v3':
            return False
        else: 
            return True
    elif vcf_calling == 'sanger':
        if not es_json.get('variant_calling_results').get(get_formal_vcf_name('sanger')).get('vcf_workflow_result_version') == 'v3':
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


def write_json(jobs_dir, job_json):
    for data_type in ['normal', 'tumors', 'sanger', 'dkfz_embl', 'broad', 'muse']:
        if not job_json.get(data_type):
            logger.warning('donor: {} has no information of {}'.format(job_json.get('project_code')+'::'+job_json.get('submitter_donor_id'), data_type))
            return

    project_code = job_json.get('project_code')
    donor_id = job_json.get('submitter_donor_id')
    
    json_name_list = [project_code, donor_id, 'json']

    json_name = '.'.join(json_name_list)
    with open(jobs_dir + '/' + json_name, 'w') as w:
        w.write(json.dumps(job_json, indent=4, sort_keys=True))


def generate_id_list(id_lists):
    ids_list = set()
    if id_lists:
        files = glob.glob(id_lists)
        for fname in files:
            with open(fname) as f:
                for d in f: ids_list.add(d.rstrip())

    return ids_list

def get_oxog_scores(oxog_scores):
    oxog_score = {}
    if oxog_scores:
       files = glob.glob(oxog_scores)
       for fname in files:
          with open(fname) as f:
              reader = csv.DictReader(f, delimiter='\t')
              for row in reader:
                  if not row.get('aliquot_GUUID'): continue
                  oxog_score[row.get('aliquot_GUUID')] = row.get('picard_oxoQ')
    return oxog_score

def create_job_json(es_json):
    job_json = {
        'project_code': es_json['dcc_project_code'],
        'submitter_donor_id': es_json['submitter_donor_id'],
        'normal': {},
        'tumors': [],
        'sanger': {},
        'dkfz_embl': {},
        'broad': {}
    }
    return job_json


def main(argv=None):

    parser = ArgumentParser(description="OxOG Whitelist Jobs Json Generator",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=True)
    parser.add_argument("-s", "--oxog_scores", dest="oxog_scores",
             help="Specify the files containing oxog_scores", required=False)
    parser.add_argument("-t", "--target_cloud", dest="target_cloud",
             help="Specify target_cloud of the jobs", required=True)
    parser.add_argument("-r", "--specify source repo", dest="chosen_gnos_repo",
             help="Specify source gnos repo", required=False)
    parser.add_argument("-d", "--exclude_donor_id_lists", dest="exclude_donor_id_lists", 
             help="File(s) containing DONOR IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-c", "--include_donor_id_lists", dest="include_donor_id_lists", 
             help="File(s) containing DONOR IDs to be excluded, use filename pattern to specify the file(s)", required=False)
 


    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports
    target_cloud = args.target_cloud
    exclude_donor_id_lists = args.exclude_donor_id_lists
    include_donor_id_lists = args.include_donor_id_lists
    chosen_gnos_repo = args.chosen_gnos_repo
    oxog_scores = args.oxog_scores

    if not oxog_scores:
        oxog_scores = '../pcawg-operations/lists/broad_qc_metrics.tsv'

    if not os.path.isdir(metadata_dir):  # TODO: should add more directory name check to make sure it's right
        sys.exit('Error: specified metadata directory does not exist!')

    timestamp = str.split(metadata_dir, '/')[-1]
    es_index = 'p_' + re.sub(r'\D', '', timestamp).replace('20','',1)
    es_type = "donor"
    es_host = 'localhost:9200'

    es = Elasticsearch([es_host])


    # pre-exclude gnos entries when this option is chosen
    donor_ids_to_be_excluded = generate_id_list(exclude_donor_id_lists)

    # read and parse git for the donors which are scheduled for Oxog filtering
    git_fnames = '../oxog-ops/oxog-*-jobs/*/*.json'

    files = glob.glob(git_fnames)
    for f in files:
        fname = str.split(f, '/')[-1]
        donor_unique_id = str.split(fname, '.')[0]+'::'+str.split(fname, '.')[1]
        donor_ids_to_be_excluded.add(donor_unique_id)

    # only process the gnos entries when this option is chosen
    donor_ids_to_be_included = generate_id_list(include_donor_id_lists) 

    if not donor_ids_to_be_included:  
        donors_list = get_donors_list(es, es_index, es_queries)
    else:
        donors_list = donor_ids_to_be_included

    # exclude the donors if they were specified on the exclude_donor_id_lists
    donors_list.difference_update(donor_ids_to_be_excluded)

    # read oxog_scores files
    oxog_score = get_oxog_scores(oxog_scores)

    # get transferred gnos_ids from git
    gnos_ids_in_cloud = set()
    if target_cloud == 'aws':
        git_s3_fnames = '../s3-transfer-operations/s3-transfer-jobs*/completed-jobs/*.json'
    elif target_cloud == 'collab':
        git_s3_fnames = '../ceph_transfer_ops/ceph-transfer-jobs*/*completed-jobs/*.json'
    elif target_cloud in ['tcga', 'gtdownload']:
        git_s3_fnames = ''
    else:
        sys.exit('Error: unknown target_cloud!')
    files = glob.glob(git_s3_fnames)
    for f in files:
        fname = str.split(f, '/')[-1]
        gnos_id = str.split(fname, '.')[0]
        gnos_ids_in_cloud.add(gnos_id)


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
        # print donor_unique_id     
        
        es_json = get_donor_json(es, es_index, donor_unique_id)

        # ensure the sanger and broad have fixed sv and snv version of files
        if not es_json.get('variant_calling_results').get('sanger_variant_calling').get('vcf_workflow_result_version') == 'v3': 
            logger.warning('donor: {} has no sanger-v3 variant calling'.format(donor_unique_id))
            continue
        if not es_json.get('variant_calling_results').get('broad_variant_calling').get('vcf_workflow_result_version') == 'v3': 
            logger.warning('donor: {} has no broad-v3 variant calling'.format(donor_unique_id))
            continue

        job_json = create_job_json(es_json)       

        add_success = add_wgs_specimens(es_json, chosen_gnos_repo, jobs_dir, job_json, oxog_score, gnos_ids_in_cloud, target_cloud)
        if not add_success: continue

        add_success = add_variant_calling(es_json, chosen_gnos_repo, jobs_dir, job_json, gnos_ids_in_cloud, target_cloud)
        if not add_success: continue

        write_json(jobs_dir, job_json)


    if os.path.isfile('tmp.xml'): os.remove('tmp.xml')

    return 0

if __name__ == "__main__":
    sys.exit(main())
