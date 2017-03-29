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
import csv
import shutil
from operator import itemgetter

es_queries = [
{ 
    "name": "bam_gnos_objects_to_sync_from_others_into_bsc",
    "content":
              {
                 "fields":[
                     "donor_unique_id"
                     ],
                  "filter":{
                      "bool":{
                            "must": [
                              {
                               "type":{
                                  "value":"donor"
                                }
                              }
                            ],
                            "should":[{
                              "bool": {
                                 "must_not": {
                                    "terms": {
                                        "wgs.normal_specimen.bwa_alignment.gnos_repo":["https://gtrepo-bsc.annailabs.com/"]
                                      }
                                   }
                              }
                            },
                            {
                              "bool":{
                                "must_not": {
                                   "nested":{
                                        "path": "wgs.tumor_specimens",
                                        "filter":{
                                          "bool":{
                                            "must":{
                                              "terms": {
                                                 "wgs.tumor_specimens.bwa_alignment.gnos_repo":["https://gtrepo-bsc.annailabs.com/"]
                                              }
                                          }  
                                       } 
                                    } 
                                  }
                                }
                              }
                            }
                          ],
                          "must_not":[
                              {"regexp": { "dcc_project_code": ".*-US"}}
                            ]    
                          }   
                    },
                 "size": 10000
                }
},
{
    "name": "bam_gnos_objects_to_sync_from_bsc_into_ebi",
    "content":
            {
               "fields":[
                     "donor_unique_id"
                     ],  
               "filter":{
                      "bool":{
                            "must": [
                              {
                               "type":{
                                  "value":"donor"
                                }
                              }
                            ],
                            "should":[{
                              "bool": {
                                 "must_not": {
                                    "terms": {
                                        "wgs.normal_specimen.bwa_alignment.gnos_repo":[
                                        "https://gtrepo-ebi.annailabs.com/", 
                                        "https://gtrepo-riken.annailabs.com/", 
                                        "https://gtrepo-dkfz.annailabs.com/", 
                                        "https://gtrepo-osdc-icgc.annailabs.com/",
                                        "https://gtrepo-etri.annailabs.com/",
                                        "https://gtrepo-osdc-tcga.annailabs.com/",
                                        "https://cghub.ucsc.edu/"]
                                      }
                                   }
                              }
                            },
                            {
                              "bool":{
                                "must_not": {
                                   "nested":{
                                        "path": "wgs.tumor_specimens",
                                        "filter":{
                                          "bool":{
                                            "must":{
                                              "terms": {
                                                 "wgs.tumor_specimens.bwa_alignment.gnos_repo":[
                                                 "https://gtrepo-ebi.annailabs.com/", 
                                                 "https://gtrepo-riken.annailabs.com/", 
                                                 "https://gtrepo-dkfz.annailabs.com/", 
                                                 "https://gtrepo-osdc-icgc.annailabs.com/",
                                                 "https://gtrepo-etri.annailabs.com/",
                                                 "https://gtrepo-osdc-tcga.annailabs.com/",
                                                 "https://cghub.ucsc.edu/"]
                                              }
                                          }  
                                       } 
                                    } 
                                  }
                                }
                              }
                            }
                          ],
                          "must_not":[
                              {"regexp": { "dcc_project_code": ".*-US"}}
                            ]    
                          }   
                    },
                "size": 10000
        }
},

]


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


def get_donors_list(es, es_index, es_queries, q_index):
    response = es.search(index=es_index, body=es_queries[q_index].get('content'))
    
    donors_list = []
    for p in response['hits']['hits']:
      donors_list.append(p.get('fields').get('donor_unique_id')[0])

    return donors_list 

def create_report_info(donor_unique_id, es_json, q_index):
    report_info_list = []

    report_info = OrderedDict()
    report_info['dcc_project_code'] = es_json['dcc_project_code']
    report_info['submitter_donor_id'] = es_json['submitter_donor_id']
    report_info['donor_unique_id'] = donor_unique_id
    
    annotations = {}
    if q_index == 0:
        add_report_info_0(report_info, report_info_list, es_json, annotations)

    if q_index == 1:
        add_report_info_1(report_info, report_info_list, es_json, annotations)

    return report_info_list


def add_report_info_0(report_info, report_info_list, es_json, annotations):
    if es_json.get('wgs').get('normal_specimen'):
        analysis = es_json.get('wgs').get('normal_specimen').get('bwa_alignment')
        if not get_formal_repo_name('bsc') in analysis.get('gnos_repo'):
            add_report_info_analysis(report_info, report_info_list, analysis, annotations)
    if es_json.get('wgs').get('tumor_specimens'):
        for aliquot in es_json.get('wgs').get('tumor_specimens'):        
            analysis = aliquot.get('bwa_alignment')
            if not get_formal_repo_name('bsc') in analysis.get('gnos_repo'):
                add_report_info_analysis(report_info, report_info_list, analysis, annotations)
    return report_info_list


def add_report_info_1(report_info, report_info_list, es_json, annotations):
    if es_json.get('wgs').get('normal_specimen'):
        analysis = es_json.get('wgs').get('normal_specimen').get('bwa_alignment')
        if get_formal_repo_name('bsc') in analysis.get('gnos_repo') and len(analysis.get('gnos_repo'))==1:
            add_report_info_analysis(report_info, report_info_list, analysis, annotations)
    if es_json.get('wgs').get('tumor_specimens'):
        for aliquot in es_json.get('wgs').get('tumor_specimens'):        
            analysis = aliquot.get('bwa_alignment')
            if get_formal_repo_name('bsc') in analysis.get('gnos_repo') and len(analysis.get('gnos_repo'))==1:
                add_report_info_analysis(report_info, report_info_list, analysis, annotations)
    return report_info_list


def add_report_info_analysis(report_info, report_info_list, analysis, annotations):
    report_info['specimen_type'] = analysis.get('specimen_type')
    report_info['aliquot_id'] = analysis.get('aliquot_id')
    report_info['gnos_id'] = analysis.get('gnos_id')
    report_info['gnos_repo'] = analysis.get('gnos_repo')
    report_info_list.append(copy.deepcopy(report_info))
    return report_info_list


def init_report_dir(metadata_dir, report_name):
    report_dir = metadata_dir + '/reports/' + report_name
    if os.path.exists(report_dir): shutil.rmtree(report_dir, ignore_errors=True)  # empty the folder if exists
    os.makedirs(report_dir)

    return report_dir


def read_annotations(annotations, type, file_name):
    if not os.path.isfile(file_name):
        return
    with open(file_name, 'r') as r:
        if annotations.get(type): # reset annotation if exists
            del annotations[type]

        if type == 'esad-uk_reheader_uuid':
            annotations[type] = set()
            for line in r:
                if line.startswith('#'): continue
                if len(line.rstrip()) == 0: continue
                annotations[type].add(line.rstrip())
        else:
            print('unknown annotation type: {}'.format(type))
    return annotations

def main(argv=None):

    parser = ArgumentParser(description="Get Donor Info For Specific Query",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=True)
    parser.add_argument("-q", "--ES_query", dest="q_index",
             help="Specify which ES_query to be used", required=False)


    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports
    q_index = args.q_index

    if not os.path.isdir(metadata_dir):  # TODO: should add more directory name check to make sure it's right
        sys.exit('Error: specified metadata directory does not exist!')

    q_index = range(len(es_queries)) if not q_index else [int(q_index)] 

    timestamp = str.split(metadata_dir, '/')[-1]
    es_index = 'pcawg_summary'
    es_type = "donor"
    es_host = 'localhost:9200'

    es = Elasticsearch([es_host], timeout=600)
  
    # output result
    report_name = re.sub(r'^generate_', '', os.path.basename(__file__))
    report_name = re.sub(r'\.py$', '', report_name)
    report_dir = init_report_dir(metadata_dir, report_name)

    for q in q_index:
        report_tsv_fh = open(report_dir + '/' + es_queries[q].get('name') + '.txt', 'w')  

        # get the list of donors
        donors_list = get_donors_list(es, es_index, es_queries, q)
        donors_list = sorted(donors_list)  

        report_info_list_full = []
        for donor_unique_id in donors_list:
            # get json doc for each donor                 
            es_json = get_donor_json(es, es_index, donor_unique_id)
            
            report_info_list_donor = create_report_info(donor_unique_id, es_json, q)

            report_info_list_full.extend(report_info_list_donor)


        header = True  
        for r in report_info_list_full:
            if header:
                report_tsv_fh.write('\t'.join(r.keys()) + '\n')
                header = False 
            # make the list of output from dict
            line = []
            for p in r.keys():
                if isinstance(r.get(p), list):
                    line.append(', '.join(str(x) for x in r.get(p)))
                elif isinstance(r.get(p), set):
                    line.append(', '.join(list(r.get(p))))
                elif r.get(p) is None:
                    line.append('')
                else:
                    line.append(str(r.get(p)))
            report_tsv_fh.write('\t'.join(line) + '\n') 
        
        report_tsv_fh.close()            


    return 0


if __name__ == "__main__":
    sys.exit(main())
