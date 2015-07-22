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

id_service_token = os.environ.get('ICGC_TOKEN')

json_prefix_code = 'a'
json_prefix_start = 1
json_prefix_inc = 10

logger = logging.getLogger('s3 transfer json generator')
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
                "flags.is_santa_cruz_donor":[
                  "T"
                ]
              }
            },
            {
              "terms": {
                "dcc_project_code": [
                  "BOCA-UK", "BRCA-UK", "BTCA-SG", "CMDI-UK", "LAML-KR", "LINC-JP", "LIRI-JP"
                  "LUSC-KR", "MELA-AU", "ORCA-IN", "OV-AU", "PACA-CA",
                  "PACA-IT", "PEME-CA", "PRAD-CA", "SKCA-BR", "THCA-SA"
                ]
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
              "regexp": {
                "dcc_project_code": ".*-US"
              }
            },
            {
              "terms": {
                "duplicated_bwa_alignment_summary.exists_mismatch_bwa_bams": [
                  "T"
                ]
              }
            },
            {
              "terms": {
                "duplicated_bwa_alignment_summary.exists_gnos_xml_mismatch": [
                  "T"
                ]
              }
            },
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


def get_source_repo_index_pos (available_repos):
    source_repo_rank = (
        "https://gtrepo-bsc.annailabs.com/",
        "https://gtrepo-dkfz.annailabs.com/",
        "https://gtrepo-osdc-icgc.annailabs.com/",
        "https://gtrepo-ebi.annailabs.com/",
        "https://gtrepo-riken.annailabs.com/",
        "https://gtrepo-etri.annailabs.com/",
    )
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
      "bsc": "bsc",
      "https://gtrepo-ebi.annailabs.com/": "ebi",
      "ebi": "ebi",
      "https://cghub.ucsc.edu/": "cghub",
      "cghub": "cghub",
      "https://gtrepo-dkfz.annailabs.com/": "dkfz",
      "dkfz": "dkfz",
      "https://gtrepo-riken.annailabs.com/": "riken",
      "riken": "riken",
      "https://gtrepo-osdc-icgc.annailabs.com/": "osdc-icgc",
      "osdc-icgc": "osdc-icgc",
      "https://gtrepo-osdc-tcga.annailabs.com/": "osdc-tcga",
      "osdc-tcga": "osdc-tcga",
      "https://gtrepo-etri.annailabs.com/": "etri",
      "etri": "etri"
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
        return 'FAKE-ID'
    elif r.json().get('totalElements') == 1:
        logger.info('GET request got the id')
        return r.json().get('content')[0].get('id')
    elif r.json().get('totalElements') > 1:
        logger.warning('GET request to metadata service return multiple matches for gnos_id: {} and filename: {}'
                          .format(gnos_id, filename))
        return 'FAKE-ID'
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
            return 'FAKE-ID'
        return r.json().get('id')
    else:
        logger.info('No luck, generate FAKE ID')
        return 'FAKE-ID'


def create_reorganized_donor(donor_unique_id, es_json):
    reorganized_donor = {
        'donor_unique_id': donor_unique_id,
        'submitter_donor_id': es_json['submitter_donor_id'],
        'dcc_project_code': es_json['dcc_project_code'],
        'is_santa_cruz': True if es_json.get('flags').get('is_santa_cruz_donor') else False,
        'wgs': {
            'normal_specimen': {},
            'tumor_specimens': []
        },
        'sanger_variant_calling':{},
        'rna_seq': {
             'normal_specimen': {},
             'tumor_specimens': []
        }
    }

    add_wgs_normal_specimen(reorganized_donor, es_json)

    add_wgs_tumor_specimens(reorganized_donor, es_json)

    add_sanger_variant_calling(reorganized_donor, es_json)

    add_rna_seq_info(reorganized_donor, es_json)

    return reorganized_donor



def add_wgs_normal_specimen(reorganized_donor, es_json):
    aliquot = es_json.get('normal_alignment_status')
    aliquot_info = create_bwa_alignment(aliquot, es_json)
    reorganized_donor.get('wgs').get('normal_specimen').update(aliquot_info)

    return reorganized_donor

def add_metadata_xml_info(obj):
    repo = get_formal_repo_name(obj.get('gnos_repo')[ get_source_repo_index_pos(obj.get('gnos_repo')) ])
    gnos_id = obj.get('gnos_id')
    ao_state = 'live'
    ao_updated = obj.get('gnos_last_modified')[ get_source_repo_index_pos(obj.get('gnos_repo')) ].encode('utf8')
    ao_updated = str.split(ao_updated, '+')[0] + 'Z'
    metadata_xml_file = 'gnos_metadata/__all_metadata_xml/' + repo + '/' + gnos_id + '__' + ao_state + '__' + ao_updated + '.xml'
    metadata_xml_file_info = {
        'file_name': gnos_id + '.xml',
        'file_md5sum': generate_md5_size(metadata_xml_file)[0],
        'file_size': generate_md5_size(metadata_xml_file)[1],
        'object_id': generate_object_id(gnos_id+'.xml', gnos_id)
    }

    return metadata_xml_file_info

def create_bwa_alignment(aliquot, es_json):
    aliquot_info = {
        'data_type': 'bwa_alignment',
        'project_code': es_json['dcc_project_code'],
        'submitter_donor_id': es_json['submitter_donor_id'],
        'is_santa_cruz': aliquot.get('aligned_bam').get('is_santa_cruz_entry'),
        'submitter_specimen_id': aliquot.get('submitter_specimen_id'),
        'submitter_sample_id': aliquot.get('submitter_sample_id'),
        'specimen_type': aliquot.get('dcc_specimen_type'),
        'aliquot_id': aliquot.get('aliquot_id'),
        'available_repos': aliquot.get('aligned_bam').get('gnos_repo'),
        'gnos_repo': [ aliquot.get('aligned_bam').get('gnos_repo')[ \
                       get_source_repo_index_pos(aliquot.get('aligned_bam').get('gnos_repo')) ] ],
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
    metadata_xml_file_info = add_metadata_xml_info(aliquot.get('aligned_bam'))
    aliquot_info.get('files').append(metadata_xml_file_info)   

    return aliquot_info


def add_wgs_tumor_specimens(reorganized_donor, es_json):
    wgs_tumor_alignment_info = es_json.get('tumor_alignment_status')

    tumor_wgs_specimen_count = 0

    for aliquot in wgs_tumor_alignment_info:
        tumor_wgs_specimen_count += 1
        aliquot_info = create_bwa_alignment(aliquot, es_json)
        reorganized_donor.get('wgs').get('tumor_specimens').append(aliquot_info) 

    reorganized_donor['tumor_wgs_specimen_count'] = tumor_wgs_specimen_count

    return reorganized_donor


def add_sanger_variant_calling(reorganized_donor, es_json):
    wgs_tumor_sanger_vcf_info = es_json.get('variant_calling_results').get('sanger_variant_calling')
    sanger_vcf_files = wgs_tumor_sanger_vcf_info.get('files')

    sanger_variant_calling = {
        'data_type': 'sanger_vcf',
        'project_code': es_json['dcc_project_code'],
        'submitter_donor_id': es_json['submitter_donor_id'],  
        'is_santa_cruz': wgs_tumor_sanger_vcf_info.get('is_santa_cruz_entry'),             
        'submitter_specimen_id': None,
        'submitter_sample_id': None,
        'specimen_type': None,
        'aliquot_id': None,
        'available_repos': wgs_tumor_sanger_vcf_info.get('gnos_repo'),
        'gnos_repo': [ wgs_tumor_sanger_vcf_info.get('gnos_repo')[ \
            get_source_repo_index_pos(wgs_tumor_sanger_vcf_info.get('gnos_repo')) ] ],
        'gnos_id': wgs_tumor_sanger_vcf_info.get('gnos_id'),
        'files': wgs_tumor_sanger_vcf_info.get('files')
    }  
    
    # add the object_id for each file object
    for f in sanger_variant_calling.get('files'):
        f.update({'file_size': None if f.get('file_size') == None else int(f.get('file_size'))})
        f.update({'object_id': generate_object_id(f.get('file_name'), sanger_variant_calling.get('gnos_id'))})

    # add the metadata_xml_file_info
    metadata_xml_file_info = add_metadata_xml_info(wgs_tumor_sanger_vcf_info)

    sanger_variant_calling.get('files').append(metadata_xml_file_info)            

    reorganized_donor.get('sanger_variant_calling').update(sanger_variant_calling) 

    return reorganized_donor


def filter_liri_jp(project, gnos_repo):
    if not project == 'LIRI-JP':
        return gnos_repo
    elif "https://gtrepo-riken.annailabs.com/" in gnos_repo:
        return ["https://gtrepo-riken.annailabs.com/"]
    else:
        print "This should never happen: alignment for LIRI-JP is not available at Riken repo"
        sys.exit(1)


def add_rna_seq_info(reorganized_donor, es_json):
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
                alignment_info[workflow_type] = create_rna_seq_alignment(aliquot, es_json, workflow_type)
            
            reorganized_donor.get('rna_seq')[specimen_type + '_specimen'] = alignment_info
        else:
            for aliquot in rna_seq_info.get(specimen_type):
                alignment_info = {}
                for workflow_type in aliquot.keys():
                    alignment_info[workflow_type] = create_rna_seq_alignment(aliquot, es_json, workflow_type)

                reorganized_donor.get('rna_seq')[specimen_type + '_specimens'].append(alignment_info) 


def create_rna_seq_alignment(aliquot, es_json, workflow_type):
    alignment_info = {
        'data_type': 'RNA_Seq_'+workflow_type,
        'project_code': es_json['dcc_project_code'],
        'submitter_donor_id': es_json['submitter_donor_id'],
        'is_santa_cruz': aliquot.get(workflow_type).get('is_santa_cruz_entry'),
        'submitter_specimen_id': aliquot.get(workflow_type).get('submitter_specimen_id'),
        'submitter_sample_id': aliquot.get(workflow_type).get('submitter_sample_id'),
        'specimen_type': aliquot.get(workflow_type).get('dcc_specimen_type'),
        'aliquot_id': aliquot.get(workflow_type).get('aliquot_id'),
        'available_repos': aliquot.get(workflow_type).get('gnos_info').get('gnos_repo'),
        'gnos_repo': [ aliquot.get(workflow_type).get('gnos_info').get('gnos_repo')[ \
            get_source_repo_index_pos(aliquot.get(workflow_type).get('gnos_info').get('gnos_repo'))] ],
        'gnos_id': aliquot.get(workflow_type).get('gnos_info').get('gnos_id'),
        'files': [
            {
                'file_name': aliquot.get(workflow_type).get('gnos_info').get('bam_file_name'),
                'file_md5sum': aliquot.get(workflow_type).get('gnos_info').get('bam_file_md5sum'),
                'file_size': aliquot.get(workflow_type).get('gnos_info').get('bam_file_size'),
                'object_id': generate_object_id(aliquot.get(workflow_type).get('gnos_info').get('bam_file_name'), aliquot.get(workflow_type).get('gnos_info').get('gnos_id'))                          
            }
        ]
    }

    # add the bai file info if exist
    if aliquot.get(workflow_type).get('gnos_info').get('bai_file_name'):
        bai_file = {
            'file_name': aliquot.get(workflow_type).get('gnos_info').get('bai_file_name'),
            'file_md5sum': aliquot.get(workflow_type).get('gnos_info').get('bai_file_md5sum'),
            'file_size': aliquot.get(workflow_type).get('gnos_info').get('bai_file_size'),
            'object_id': generate_object_id(aliquot.get(workflow_type).get('gnos_info').get('bai_file_name'), aliquot.get(workflow_type).get('gnos_info').get('gnos_id'))                        
        }
        alignment_info.get('files').append(bai_file)
    else:
        logger.warning('RNA-Seq alignment GNOS entry {} has no .bai file'.format(alignment_info.get('gnos_id')))

    # add the metadata_xml_file_info
    metadata_xml_file_info = add_metadata_xml_info(aliquot.get(workflow_type).get('gnos_info'))
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
    
    donors_list = []
    for p in response['hits']['hits']:
    	donors_list.append(p.get('fields').get('donor_unique_id')[0])

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

def organize_s3_transfer(jobs_dir, reorganized_donor, gnos_ids_to_be_included, gnos_ids_to_be_excluded):

    if reorganized_donor.get('wgs').get('normal_specimen'):
        transfer_json = reorganized_donor.get('wgs').get('normal_specimen')
        write_s3_transfer_json(jobs_dir, transfer_json, gnos_ids_to_be_included, gnos_ids_to_be_excluded)

    if reorganized_donor.get('wgs').get('tumor_specimens'):
        for transfer_json in reorganized_donor.get('wgs').get('tumor_specimens'):
            write_s3_transfer_json(jobs_dir, transfer_json, gnos_ids_to_be_included, gnos_ids_to_be_excluded)

    if reorganized_donor.get('sanger_variant_calling'):
        transfer_json = reorganized_donor.get('sanger_variant_calling')
        write_s3_transfer_json(jobs_dir, transfer_json, gnos_ids_to_be_included, gnos_ids_to_be_excluded)

    if reorganized_donor.get('rna_seq').get('normal_specimen'):
        aliquot = reorganized_donor.get('rna_seq').get('normal_specimen')
        for workflow_type in aliquot.keys():
            transfer_json = aliquot.get(workflow_type)
            write_s3_transfer_json(jobs_dir, transfer_json, gnos_ids_to_be_included, gnos_ids_to_be_excluded)

    if reorganized_donor.get('rna_seq').get('tumor_specimens'):
        for aliquot in reorganized_donor.get('rna_seq').get('tumor_specimens'):
            for workflow_type in aliquot.keys():
                transfer_json = aliquot.get(workflow_type)
                write_s3_transfer_json(jobs_dir, transfer_json, gnos_ids_to_be_included, gnos_ids_to_be_excluded)


def write_s3_transfer_json(jobs_dir, transfer_json, gnos_ids_to_be_included, gnos_ids_to_be_excluded):
    global json_prefix_code, json_prefix_start, json_prefix_inc

    #if (json_prefix_start > 41): sys.exit()  # for debugging only to terminate earlier

    if transfer_json.get('is_santa_cruz'):
        gnos_id = transfer_json.get('gnos_id')
        # this is not optimistic, such filter could be done earlier in the process, but let's be it for now
        if gnos_ids_to_be_included and not gnos_id in gnos_ids_to_be_included: return
        if gnos_ids_to_be_excluded and gnos_id in gnos_ids_to_be_excluded: return

        prefix_for_priority = json_prefix_code + '0'*(6-len(str(json_prefix_start))) + str(json_prefix_start)
        project_code = transfer_json.get('project_code')
        data_type = transfer_json.get('data_type')
        json_name_list = [gnos_id, project_code, data_type, 'json']
        json_name = '.'.join(json_name_list)
        with open(jobs_dir + '/' + json_name, 'w') as w:
            w.write(json.dumps(transfer_json, indent=4, sort_keys=True))
            json_prefix_start = json_prefix_start + json_prefix_inc


def generate_gnos_id_list(gnos_id_lists):
    gnos_ids_list = set()
    if gnos_id_lists:
        files = glob.glob(gnos_id_lists)
        for fname in files:
            with open(fname) as f:
                for d in f: gnos_ids_list.add(d.rstrip())

    return gnos_ids_list


def main(argv=None):

    parser = ArgumentParser(description="S3 Transfer Jobs Json Generator",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=True)
    parser.add_argument("-i", "--include_gnos_id_lists", dest="include_gnos_id_lists",
             help="Specify which GNOS IDs to process, process all gnos_ids if none specified", required=False)
    parser.add_argument("-x", "--exclude_gnos_id_lists", dest="exclude_gnos_id_lists", 
             help="File(s) containing GNOS IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-c", "--json_prefix_code", dest="prefix_code",
             help="Json file prefix single letter code", required=False)
    parser.add_argument("-s", "--json_prefix_start", dest="prefix_start",
             help="Directory containing metadata manifest files", required=False)
    parser.add_argument("-n", "--json_prefix_inc", dest="prefix_inc",
             help="Directory containing metadata manifest files", required=False)

    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports
    include_gnos_id_lists = args.include_gnos_id_lists
    exclude_gnos_id_lists = args.exclude_gnos_id_lists

    # pre-exclude gnos entries when this option is chosen
    gnos_ids_to_be_excluded = generate_gnos_id_list(exclude_gnos_id_lists)

    # only process the gnos entries when this option is chosen
    gnos_ids_to_be_included = generate_gnos_id_list(include_gnos_id_lists)    

    global json_prefix_code, json_prefix_start, json_prefix_inc
    if args.prefix_code: json_prefix_code = args.prefix_code
    if args.prefix_start: json_prefix_start = int(args.prefix_start)
    if args.prefix_inc: json_prefix_inc = int(args.prefix_inc)

    if not os.path.isdir(metadata_dir):  # TODO: should add more directory name check to make sure it's right
        sys.exit('Error: specified metadata directory does not exist!')

    timestamp = str.split(metadata_dir, '/')[-1]
    es_index = 'p_' + re.sub(r'\D', '', timestamp).replace('20','',1)
    es_type = "donor"
    es_host = 'localhost:9200'

    es = Elasticsearch([es_host])

	# get the list of donors which is santa cruz
    donors_list = get_donors_list(es, es_index, es_queries)

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


    donor_fh = open(jobs_dir+'/s3_transfer_json.jsonl', 'w')
    
    # get json doc for each donor and reorganize it 
    for donor_unique_id in donors_list:     
        
    	es_json = get_donor_json(es, es_index, donor_unique_id)
        
        reorganized_donor = create_reorganized_donor(donor_unique_id, es_json)
    
        organize_s3_transfer(jobs_dir, reorganized_donor, gnos_ids_to_be_included, gnos_ids_to_be_excluded)

        donor_fh.write(json.dumps(reorganized_donor, default=set_default, sort_keys=True) + '\n')

    donor_fh.close()

    os.remove('tmp.xml')

    return 0


if __name__ == "__main__":
    sys.exit(main())
