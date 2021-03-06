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

def create_specimen_info(donor_unique_id, es_json, donor_ids):
    specimen_info_list = []

    specimen_info = OrderedDict()
    specimen_info['donor_unique_id'] = donor_unique_id
    specimen_info['donor_wgs_exclusion_white_gray'] = get_flag(donor_unique_id, donor_ids)
    specimen_info['submitter_donor_id'] = es_json['submitter_donor_id']
    specimen_info['icgc_donor_id'] = es_json['icgc_donor_id'] if es_json['icgc_donor_id'] else ''
    specimen_info['dcc_project_code'] = es_json['dcc_project_code']
    
    add_wgs_specimens(specimen_info_list, specimen_info, es_json)

    add_rna_seq_specimens(specimen_info_list, specimen_info, es_json)

    return specimen_info_list

def get_flag(donor_unique_id, donor_ids):
    for t in ['white', 'black', 'gray']:
        if not donor_unique_id in donor_ids.get(t+'list'): continue
        if t == 'black':
            return 'Excluded'
        else:
            return t.capitalize()+'list'
    return None
        

def add_wgs_specimens(specimen_info_list, specimen_info, es_json):

    if es_json.get('normal_alignment_status'):
        specimen_info['aliquot_id'] = es_json.get('normal_alignment_status').get('aliquot_id')
        specimen_info['submitter_specimen_id'] = es_json.get('normal_alignment_status').get('submitter_specimen_id')
        specimen_info['icgc_specimen_id'] = es_json.get('normal_alignment_status').get('icgc_specimen_id') \
                                            if es_json.get('normal_alignment_status').get('icgc_specimen_id') else ''
        specimen_info['submitter_sample_id'] = es_json.get('normal_alignment_status').get('submitter_sample_id')
        specimen_info['icgc_sample_id'] = es_json.get('normal_alignment_status').get('icgc_sample_id') \
                                          if es_json.get('normal_alignment_status').get('icgc_sample_id') else ''
        specimen_info['dcc_specimen_type'] = es_json.get('normal_alignment_status').get('dcc_specimen_type')
        specimen_info['library_strategy'] = 'WGS'
        specimen_info_list.append(copy.deepcopy(specimen_info))

    if es_json.get('tumor_alignment_status'):
        for aliquot in es_json.get('tumor_alignment_status'):
    	    specimen_info['aliquot_id'] = aliquot.get('aliquot_id')
            specimen_info['submitter_specimen_id'] = aliquot.get('submitter_specimen_id')
            specimen_info['icgc_specimen_id'] = aliquot.get('icgc_specimen_id') if aliquot.get('icgc_specimen_id') else ''
            specimen_info['submitter_sample_id'] = aliquot.get('submitter_sample_id')
            specimen_info['icgc_sample_id'] = aliquot.get('icgc_sample_id') if aliquot.get('icgc_sample_id') else ''
            specimen_info['dcc_specimen_type'] = aliquot.get('dcc_specimen_type')
            specimen_info['library_strategy'] = 'WGS'
            specimen_info_list.append(copy.deepcopy(specimen_info))

    return specimen_info_list


def add_rna_seq_specimens(specimen_info_list, specimen_info, es_json):
    # to build pcawg santa cruz pilot dataset, this is a temporary walkaround to exclude the 130 RNA-Seq bad
    # entries from MALY-DE and CLLE-ES projects
    #if es_json.get('dcc_project_code') in ('MALY-DE', 'CLLE-ES'): return

    rna_seq_info = es_json.get('rna_seq').get('alignment')
    for specimen_type in rna_seq_info.keys():
        if not rna_seq_info.get(specimen_type): # the specimen_type has no alignment result
		    continue
        if 'normal' in specimen_type:
            aliquot = rna_seq_info.get(specimen_type)
            specimen_info['aliquot_id'] = set()
            specimen_info['submitter_specimen_id'] = set()
            specimen_info['icgc_specimen_id'] = set()
            specimen_info['submitter_sample_id'] = set()
            specimen_info['icgc_sample_id'] = set()
            specimen_info['dcc_specimen_type'] = set()
            for workflow_type in aliquot.keys():
                specimen_info['aliquot_id'].add(aliquot.get(workflow_type).get('aliquot_id'))
                specimen_info['submitter_specimen_id'].add(aliquot.get(workflow_type).get('submitter_specimen_id'))
                specimen_info['icgc_specimen_id'].add(
                        aliquot.get(workflow_type).get('icgc_specimen_id') if aliquot.get(workflow_type).get('icgc_specimen_id') else ''
                    )
                specimen_info['submitter_sample_id'].add(aliquot.get(workflow_type).get('submitter_sample_id'))
                specimen_info['icgc_sample_id'].add(
                        aliquot.get(workflow_type).get('icgc_sample_id') if aliquot.get(workflow_type).get('icgc_sample_id') else ''
                    )
                specimen_info['dcc_specimen_type'].add(aliquot.get(workflow_type).get('dcc_specimen_type'))
                specimen_info['library_strategy'] = 'RNA-Seq'
            specimen_info_list.append(copy.deepcopy(specimen_info))

        else:
            for aliquot in rna_seq_info.get(specimen_type):
                specimen_info['aliquot_id'] = set()
                specimen_info['submitter_specimen_id'] = set()
                specimen_info['icgc_specimen_id'] = set()
                specimen_info['submitter_sample_id'] = set()
                specimen_info['icgc_sample_id'] = set()
                specimen_info['dcc_specimen_type'] = set()
                for workflow_type in aliquot.keys():
                    specimen_info['aliquot_id'].add(aliquot.get(workflow_type).get('aliquot_id'))
                    specimen_info['submitter_specimen_id'].add(aliquot.get(workflow_type).get('submitter_specimen_id'))
                    specimen_info['icgc_specimen_id'].add(
                            aliquot.get(workflow_type).get('icgc_specimen_id') if aliquot.get(workflow_type).get('icgc_specimen_id') else ''
                        )
                    specimen_info['submitter_sample_id'].add(aliquot.get(workflow_type).get('submitter_sample_id'))
                    specimen_info['icgc_sample_id'].add(
                            aliquot.get(workflow_type).get('icgc_sample_id') if aliquot.get(workflow_type).get('icgc_sample_id') else ''
                        )
                    specimen_info['dcc_specimen_type'].add(aliquot.get(workflow_type).get('dcc_specimen_type'))
                    specimen_info['library_strategy'] = 'RNA-Seq'
                specimen_info_list.append(copy.deepcopy(specimen_info))

    return specimen_info_list


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


def set_default(obj):
    if isinstance(obj, datetime.datetime):
        return obj.isoformat()
    if isinstance(obj, set):
        return list(obj)
    raise TypeError

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


def main(argv=None):

    parser = ArgumentParser(description="PCAWG Full List of Specimen Level Info Generator",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=True)
    parser.add_argument("-r", "--gnos_repo", dest="repo",
             help="Specify which GNOS repo to process, process all repos if none specified", required=False)

    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports
    repo = args.repo

    if not os.path.isdir(metadata_dir):  # TODO: should add more directory name check to make sure it's right
        sys.exit('Error: specified metadata directory does not exist!')

    timestamp = str.split(metadata_dir, '/')[-1]
    es_index = 'p_' + ('' if not repo else repo+'_') + re.sub(r'\D', '', timestamp).replace('20','',1)
    es_type = "donor"
    es_host = 'localhost:9200'

    es = Elasticsearch([es_host], timeout=600)

    donor_ids = {}
    # get the list of donors in blacklist/graylist donors
    for dtype in ['blacklist', 'graylist']:
        donor_ids[dtype] = generate_id_list('../pcawg-operations/lists/'+dtype+'/pc_annotation-donor_'+dtype+'.tsv')

    # get the list of donors to be included in the release 
    donor_ids['release'] = get_donors_list(es, es_index, es_queries)
    donor_ids['whitelist'] = donor_ids.get('release').difference(donor_ids.get('blacklist').union(donor_ids.get('graylist')))

    for dtype in ['release', 'whitelist', 'blacklist', 'graylist']:
        if dtype == 'release':
            PCAWG_specimen_tsv_fh = open(metadata_dir + '/reports/pcawg_sample_sheet.tsv', 'w')         
        else:
            PCAWG_specimen_tsv_fh = open(metadata_dir + '/reports/pcawg_sample_sheet.'+dtype+'ed_donors.tsv', 'w')
        
        donors_list = donor_ids.get(dtype)
        donors_list = sorted(donors_list)

        header = True
        # get json doc for each donor and reorganize it 
        for donor_unique_id in donors_list:
            es_json = get_donor_json(es, es_index, donor_unique_id)
            if not es_json: continue
            specimen_info_list = create_specimen_info(donor_unique_id, es_json, donor_ids)

            for specimen in specimen_info_list: 
                # write to the tsv file
                if header:
                    PCAWG_specimen_tsv_fh.write('\t'.join(specimen.keys()) + '\n')
                    header = False
                line = []
                for p in specimen.keys():
                    if isinstance(specimen.get(p), set):
                        line.append('|'.join(list(specimen.get(p))))
                    else:
                        line.append(str(specimen.get(p)))
                PCAWG_specimen_tsv_fh.write('\t'.join(line)+'\n')
            
        PCAWG_specimen_tsv_fh.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
