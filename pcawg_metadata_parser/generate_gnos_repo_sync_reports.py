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
import csv
import shutil
from operator import itemgetter



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
                        },
                        {
                          "regexp": { 
                             "dcc_project_code": ".*-US"
                            }
                        }   
                      ]
                    }
                },
      "size": 10000
    }
]

def create_report(es, es_index, es_queries, report_dir, seq, vcf):
    # get the full list of donors in PCAWG
    donors_list = get_donors_list(es, es_index, es_queries)

    report = OrderedDict()
    subreport = [
        'gnos_objects_to_sync_from_others_into_bsc',
        'gnos_objects_to_sync_from_bsc_into_ebi',
        'gnos_objects_to_delete_from_osdc_icgc',
        'gnos_objects_to_sync_from_osdc_icgc_into_ebi'
    ]

    for r in subreport:
        if not report.get(r): report[r] = []

    # get json doc for each donor and reorganize it 
    for donor_unique_id in donors_list:     
        
        es_json = get_donor_json(es, es_index, donor_unique_id)
        
        create_report_info(donor_unique_id, es_json, report, seq, vcf)
    
    for r in subreport:
        if not report.get(r): continue
        write_tsv_file(report.get(r), os.path.join(report_dir, r+'.txt')) 


def create_report_info(donor_unique_id, es_json, report, seq, vcf):
    gnos_entity_info = OrderedDict()
    gnos_entity_info['donor_unique_id'] = donor_unique_id
    gnos_entity_info['submitter_donor_id'] = es_json['submitter_donor_id']
    gnos_entity_info['dcc_project_code'] = es_json['dcc_project_code']

    if seq and 'wgs' in seq:
        add_wgs_gnos_entity(report, gnos_entity_info, es_json)

    if vcf:
        add_vcf_gnos_entity(report, gnos_entity_info, es_json, vcf)

    if seq and 'rna_seq' in seq:
        add_rna_seq_gnos_entity(report, gnos_entity_info, es_json)

    return report

def add_subreport_info(report, analysis, gnos_entity_info):
    if not get_formal_repo_name('bsc') in analysis.get('gnos_repo'):    
        report['gnos_objects_to_sync_from_others_into_bsc'].append(copy.deepcopy(gnos_entity_info))
    if get_formal_repo_name('bsc') in analysis.get('gnos_repo') and len(analysis.get('gnos_repo'))==1:
        report['gnos_objects_to_sync_from_bsc_into_ebi'].append(copy.deepcopy(gnos_entity_info))
    if get_formal_repo_name('osdc-icgc') in analysis.get('gnos_repo') and len(analysis.get('gnos_repo')) > 1:
        report['gnos_objects_to_delete_from_osdc_icgc'].append(copy.deepcopy(gnos_entity_info))        
    if get_formal_repo_name('osdc-icgc') in analysis.get('gnos_repo') and len(analysis.get('gnos_repo')) == 1:
        report['gnos_objects_to_sync_from_osdc_icgc_into_ebi'].append(copy.deepcopy(gnos_entity_info)) 
    return report

def add_wgs_gnos_entity(report, gnos_entity_info, es_json):

    if es_json.get('normal_alignment_status'):
        add_wgs_aliquot_gnos_entity(es_json.get('normal_alignment_status'), gnos_entity_info, report)

    if es_json.get('tumor_alignment_status'):
        for aliquot in es_json.get('tumor_alignment_status'):
            add_wgs_aliquot_gnos_entity(aliquot, gnos_entity_info, report)

    return report


def add_wgs_aliquot_gnos_entity(aliquot, gnos_entity_info, report):
    gnos_entity_info['library_strategy'] = 'WGS'
    gnos_entity_info['aliquot_id'] = aliquot.get('aliquot_id')
    gnos_entity_info['submitter_specimen_id'] = aliquot.get('submitter_specimen_id')
    gnos_entity_info['submitter_sample_id'] = aliquot.get('submitter_sample_id')
    gnos_entity_info['dcc_specimen_type'] = aliquot.get('dcc_specimen_type')

    for bam_type in ['aligned_bam', 'bam_with_unmappable_reads', 'minibam']:
        if aliquot.get(bam_type):
            gnos_entity_info['entity_type'] = bam_type
            analysis = aliquot.get(bam_type)
            gnos_entity_info['gnos_id'] = analysis.get('gnos_id')
            gnos_entity_info['gnos_repo'] = analysis.get('gnos_repo')
            add_subreport_info(report, analysis, gnos_entity_info)        

    return report

def add_vcf_gnos_entity(report, gnos_entity_info, es_json, vcf):
    if es_json.get('variant_calling_results'):
        gnos_entity_info['library_strategy'] = 'WGS'
        gnos_entity_info['aliquot_id'] = None
        gnos_entity_info['submitter_specimen_id'] = None
        gnos_entity_info['submitter_sample_id'] = None
        gnos_entity_info['dcc_specimen_type'] = None
        for vcf_type in [v+'_variant_calling' for v in vcf]:
            if es_json.get('variant_calling_results').get(vcf_type):
                gnos_entity_info['entity_type'] = vcf_type
                analysis = es_json.get('variant_calling_results').get(vcf_type)
                gnos_entity_info['gnos_id'] = analysis.get('gnos_id')
                gnos_entity_info['gnos_repo'] = analysis.get('gnos_repo')
                add_subreport_info(report, analysis, gnos_entity_info)

    return report


def add_rna_seq_gnos_entity(report, gnos_entity_info, es_json):
    rna_seq_info = es_json.get('rna_seq').get('alignment')
    for specimen_type in rna_seq_info.keys():
        if not rna_seq_info.get(specimen_type): # the specimen_type has no alignment result
		    continue
        if 'normal' in specimen_type:
            aliquot = rna_seq_info.get(specimen_type)
            add_rna_seq_aliquot_gnos_entity(aliquot, gnos_entity_info, report)

        else:
            for aliquot in rna_seq_info.get(specimen_type):
                add_rna_seq_aliquot_gnos_entity(aliquot, gnos_entity_info, report)

    return report

def add_rna_seq_aliquot_gnos_entity(aliquot, gnos_entity_info, report):
    gnos_entity_info['library_strategy'] = 'RNA-Seq'
    gnos_entity_info['aliquot_id'] = set()
    gnos_entity_info['submitter_specimen_id'] = set()
    gnos_entity_info['submitter_sample_id'] = set()
    gnos_entity_info['dcc_specimen_type'] = set()
    for workflow_type in aliquot.keys():
        gnos_entity_info['aliquot_id'].add(aliquot.get(workflow_type).get('aliquot_id'))
        gnos_entity_info['submitter_specimen_id'].add(aliquot.get(workflow_type).get('submitter_specimen_id'))
        gnos_entity_info['submitter_sample_id'].add(aliquot.get(workflow_type).get('submitter_sample_id'))
        gnos_entity_info['dcc_specimen_type'].add(aliquot.get(workflow_type).get('dcc_specimen_type'))

    for workflow_type in aliquot.keys():
        gnos_entity_info['entity_type'] = workflow_type
        analysis = aliquot.get(workflow_type).get('aligned_bam')
        gnos_entity_info['gnos_id'] = analysis.get('gnos_id')
        gnos_entity_info['gnos_repo'] = analysis.get('gnos_repo')
        add_subreport_info(report, analysis, gnos_entity_info)
          
    return report

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

def init_report_dir(metadata_dir, report_name):
    report_dir = metadata_dir + '/reports/' + report_name
    if os.path.exists(report_dir): shutil.rmtree(report_dir, ignore_errors=True)  # empty the folder if exists
    os.makedirs(report_dir)

    return report_dir

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

def write_tsv_file(report, filename):
    with open(filename, 'w') as fh:
        header = True  
        for r in report:
            if header:
                fh.write('\t'.join(r.keys()) + '\n')
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
            fh.write('\t'.join(line) + '\n') 


def main(argv=None):

    parser = ArgumentParser(description="PCAWG Full List of GNOS entities Info Generator",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=True)
    parser.add_argument("-s", "--sequence_type", dest="seq", nargs="*",
             help="List sequence_type types [wgs, rna_seq]", required=False)
    parser.add_argument("-v", "--variant_calling", dest="vcf", nargs="*",
             help="List variant_calling types [sanger, dkfz, broad, muse, broad_tar]", required=False) 

    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports
    seq = args.seq
    vcf = args.vcf

    seq= list(seq) if seq else [] 
    vcf = list(vcf) if vcf else [] 

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

    create_report(es, es_index, es_queries, report_dir, seq, vcf)

    return 0


if __name__ == "__main__":
    sys.exit(main())
