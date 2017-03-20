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
import shutil



es_queries = [
  # query 0: PCAWG_full_list_donors 
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
                      "must_not": [
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

def create_gnos_entity_info(donor_unique_id, es_json):
    gnos_entity_info_list = []

    gnos_entity_info = OrderedDict()
    gnos_entity_info['donor_unique_id'] = donor_unique_id
    gnos_entity_info['submitter_donor_id'] = es_json['submitter_donor_id']
    gnos_entity_info['dcc_project_code'] = es_json['dcc_project_code']
    
    add_wgs_gnos_entity(gnos_entity_info_list, gnos_entity_info, es_json)

    return gnos_entity_info_list


def add_wgs_gnos_entity(gnos_entity_info_list, gnos_entity_info, es_json):

    if es_json.get('normal_alignment_status'):
        add_wgs_aliquot_gnos_entity(es_json.get('normal_alignment_status'), gnos_entity_info, gnos_entity_info_list)

    if es_json.get('tumor_alignment_status'):
        for aliquot in es_json.get('tumor_alignment_status'):
            add_wgs_aliquot_gnos_entity(aliquot, gnos_entity_info, gnos_entity_info_list)

    return gnos_entity_info_list


def add_wgs_aliquot_gnos_entity(aliquot, gnos_entity_info, gnos_entity_info_list):
    gnos_entity_info['library_strategy'] = 'WGS'
    gnos_entity_info['aliquot_id'] = aliquot.get('aliquot_id')
    gnos_entity_info['submitter_specimen_id'] = aliquot.get('submitter_specimen_id')
    gnos_entity_info['submitter_sample_id'] = aliquot.get('submitter_sample_id')
    gnos_entity_info['dcc_specimen_type'] = aliquot.get('dcc_specimen_type')       

    if aliquot.get('unaligned_bams'):
        gnos_entity_info['entity_type'] = 'unaligned_bams'
        for unaligned_bams in aliquot.get('unaligned_bams'):
            gnos_entity_info['gnos_id'] = unaligned_bams.get('gnos_id')
            for gnos_repo in unaligned_bams.get('gnos_repo'):
                gnos_entity_info['gnos_repo'] = gnos_repo
                gnos_entity_info['gnos_metadata_url'] = gnos_repo + 'cghub/metadata/analysisFull/' + gnos_entity_info['gnos_id']
                gnos_entity_info_list.append(copy.deepcopy(gnos_entity_info))  

    return gnos_entity_info_list


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


def set_default(obj):
    if isinstance(obj, datetime.datetime):
        return obj.isoformat()
    if isinstance(obj, set):
        return list(obj)
    raise TypeError

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

def init_report_dir(metadata_dir, report_name):
    report_dir = metadata_dir + '/reports/' + report_name
    if os.path.exists(report_dir): shutil.rmtree(report_dir, ignore_errors=True)  # empty the folder if exists
    os.makedirs(report_dir)

    return report_dir


def main(argv=None):

    parser = ArgumentParser(description="PCAWG Full List of GNOS entities Info Generator",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=True)

    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports

    if not os.path.isdir(metadata_dir):  # TODO: should add more directory name check to make sure it's right
        sys.exit('Error: specified metadata directory does not exist!')

    timestamp = str.split(metadata_dir, '/')[-1]
    es_index = 'p_' + re.sub(r'\D', '', timestamp).replace('20','',1)
    es_type = "donor"
    es_host = 'localhost:9200'

    es = Elasticsearch([es_host])

    # output result
    report_name = re.sub(r'^generate_', '', os.path.basename(__file__))
    report_name = re.sub(r'\.py$', '', report_name)
    report_dir = init_report_dir(metadata_dir, report_name)


	# get the list of donors in PCAWG
    donors_list = get_donors_list(es, es_index, es_queries)
    
    report_info_list_full = []
    # get json doc for each donor and reorganize it 
    for donor_unique_id in donors_list:            
    	es_json = get_donor_json(es, es_index, donor_unique_id)
        report_info_list_donor = create_gnos_entity_info(donor_unique_id, es_json)
        report_info_list_full.extend(report_info_list_donor)

                
    for repo in ['bsc', 'ebi', 'cghub', 'dkfz', 'riken', 'osdc-icgc', 'osdc-tcga', 'etri']:
        header = True
        report_tsv_fh = open(report_dir + '/' + repo + '.lane_level.analysis_id.txt', 'w')
        for gnos_entity in report_info_list_full: 
            if header:
                report_tsv_fh.write('\t'.join(gnos_entity.keys()) + '\n')
                header = False 
            # write to the tsv file
            if not get_formal_repo_name(gnos_entity.get('gnos_repo')) == repo: continue
            for p in gnos_entity.keys():
                if isinstance(gnos_entity.get(p), set):
                    report_tsv_fh.write('|'.join(list(gnos_entity.get(p))) + '\t')
                elif not gnos_entity.get(p):
                    report_tsv_fh.write('\t')
                else:
                    report_tsv_fh.write(str(gnos_entity.get(p)) + '\t')
            report_tsv_fh.write('\n')    
        report_tsv_fh.close()


    return 0


if __name__ == "__main__":
    sys.exit(main())
