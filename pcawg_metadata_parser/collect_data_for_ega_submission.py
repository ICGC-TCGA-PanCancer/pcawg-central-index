#!/usr/bin/env python

import os
import re
import click
import csv
import copy
import xmltodict
import json
import sys
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from elasticsearch import Elasticsearch
import logging
import glob
from collections import OrderedDict

logger = logging.getLogger('Collect sample and gnos_xml data')
ch = logging.StreamHandler()


def generate_es_query(dcc_project_code):
    es_query = {
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
                    "dcc_project_code": [dcc_project_code]
                  }
                },
                # {
                # "terms": {
                #   "wgs.normal_specimen.bwa_alignment.gnos_repo":["https://gtrepo-bsc.annailabs.com/"]
                #   }
                # },
                # {
                #   "nested":{
                #         "path": "wgs.tumor_specimens",
                #         "filter":{
                #           "bool":{
                #             "must":{
                #               "terms": {
                #                  "wgs.tumor_specimens.bwa_alignment.gnos_repo":["https://gtrepo-bsc.annailabs.com/"]
                #               }
                #           }  
                #        } 
                #     } 
                #   }
                # }
              ],
              "must_not": [
                {
                  "regexp": {
                    "dcc_project_code": ".*-US"
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
                {
                  "terms": {
                    "flags.exists_xml_md5sum_mismatch": [
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
    return es_query


def get_donors_list(es, es_index, dcc_project_code):
    response = es.search(index=es_index, body=generate_es_query(dcc_project_code))
    
    donors_list = set()
    for p in response['hits']['hits']:
        donors_list.add(p.get('fields').get('donor_unique_id')[0])

    return donors_list 


def collect_sample(donors_list, sample_ids_to_be_included, sample_ids_to_be_excluded, dcc_project_code, ega_dir, pcawg_sample_sheet, seq, annotations):

    sample_sheet = []
    with open(pcawg_sample_sheet, 'r') as s:
        reader = csv.DictReader(s, delimiter='\t')
        for sample_info in reader:
            if not sample_info.get('dcc_project_code') == dcc_project_code: continue
            if not sample_info.get('donor_unique_id') in donors_list: continue
            sample_id = sample_info.get('icgc_sample_id')
            if sample_ids_to_be_included and not sample_id in sample_ids_to_be_included: continue
            if sample_ids_to_be_excluded and sample_id in sample_ids_to_be_excluded: continue
            library_strategy = sample_info.get('library_strategy').lower()
            if not library_strategy in seq: continue
        
            sample = OrderedDict()

            try:
                sample['subject_id'] = sample_info['icgc_donor_id']
                if annotations.get('gender') and annotations.get('gender').get(sample_info['icgc_donor_id']):
                    sample['gender'] = annotations.get('gender').get(sample_info['icgc_donor_id']) 
                else:
                    click.echo('Warning: missing gender informaion for donor: %s' % sample_info['icgc_donor_id'], err=True)
                    return
                # if annotations.get('phenotype') and annotations.get('phenotype').get(sample_info['icgc_donor_id']):
                #     sample['phenotype'] = annotations.get('gender').get(sample_info['icgc_donor_id']) 
                # else:
                #     click.echo('Warning: missing phenotype informaion for donor: %s' % sample_info['icgc_donor_id'], err=True)
                #     return

                sample['icgc_project_code'] = sample_info['dcc_project_code']

                for tag in ['icgc_donor_id', 'submitter_donor_id', 'icgc_specimen_id', 'submitter_specimen_id', 'icgc_sample_id', 'submitter_sample_id']:
                    if sample_info.get(tag):
                        sample[tag] = sample_info.get(tag)  
                    else: 
                        click.echo('Warning: missing  %s informaion for donor: %s' % tag, sample_info['icgc_donor_id'], err=True)
                        return
                sample['specimen_type'] = sample_info['dcc_specimen_type']
                sample['aliquot_id/sample_uuid'] = sample_info['aliquot_id']


            except KeyError, e:
                click.echo('Error: KeyError, %s' % str(e), err=True)
                return
            except IndexError, e:
                click.echo('Error: IndexError, %s' % str(e), err=True)
                return
            except Exception, e:
                click.echo('Error: %s' % str(e), err=True)
                return

            sample_sheet.append(copy.deepcopy(sample))

    if sample_sheet:    
        out_dir = os.path.join(ega_dir, dcc_project_code, 'sample')
        if not os.path.isdir(out_dir): os.makedirs(out_dir)
        out_file = os.path.join(out_dir, 'sample.'+dcc_project_code+'.'+str(len(sample_sheet))+'.tsv')
        write_tsv_file(sample_sheet, out_file)


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



def read_annotations(annotations, type, file_name):
    if not os.path.isfile(file_name):
        return
    with open(file_name, 'r') as r:
        if annotations.get(type): # reset annotation if exists
            del annotations[type]               
        
        if type == 'gender':
            annotations[type] = {}
            reader = csv.DictReader(r, delimiter='\t')
            for row in reader:
                if not row.get('study_donor_involved_in') == 'PCAWG': continue
                annotations[type][row.get('icgc_donor_id')] = row.get('donor_sex') if row.get('donor_sex') else None

        else:
            logger.warning('unknown annotation type: {}'.format(type))


def generate_exclude_list(ega_dir, project, gnos_sample_ids_to_be_excluded):
    # read and parse git for the gnos_ids and fnames which are scheduled for s3 transfer
    sample_fnames = os.path.join(ega_dir, project, 'sample', 'sample.'+project+'.*.tsv')
    files = glob.glob(sample_fnames)
    for f in files:
        with open(f, 'r') as fn:
            reader = csv.DictReader(fn, delimiter='\t')
            for r in reader:
                gnos_sample_ids_to_be_excluded.add(r.get('icgc_sample_id'))

    return gnos_sample_ids_to_be_excluded

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
    parser.add_argument("-p", "--dcc_project_code", dest="dcc_project_code", nargs="*",
             help="Specify dcc_project_code for EGA submission", required=False)
    parser.add_argument("-r", "--specify source repo", dest="chosen_gnos_repo",
             help="Specify source gnos repo", required=False)

    parser.add_argument("-c", "--include_donor_id_lists", dest="include_donor_id_lists", 
             help="File(s) containing DONOR IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-d", "--exclude_donor_id_lists", dest="exclude_donor_id_lists", 
             help="File(s) containing DONOR IDs to be excluded, use filename pattern to specify the file(s)", required=False)

    parser.add_argument("-i", "--include_gnos_sample_id_lists", dest="include_gnos_sample_id_lists",
             help="File(s) containing GNOS/SAMPLE IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-x", "--exclude_gnos_sample_id_lists", dest="exclude_gnos_sample_id_lists", 
             help="File(s) containing GNOS/SAMPLE IDs to be excluded, use filename pattern to specify the file(s)", required=False)

    parser.add_argument("-s", "--sequence_type", dest="seq", nargs="*",
             help="List sequence_type types[wgs, rna-seq]", required=False)
    parser.add_argument("-w", "--workflow", dest="workflow", nargs="*",
             help="List workflow types[bwa, sanger, dkfz, broad, rna-seq]", required=False)    



    args = parser.parse_args()
    dcc_project_code = list(args.dcc_project_code) if dcc_project_code else []
    chosen_gnos_repo = args.chosen_gnos_repo

    include_donor_id_lists = args.include_donor_id_lists
    exclude_donor_id_lists = args.exclude_donor_id_lists

    include_gnos_sample_id_lists = args.include_gnos_sample_id_lists
    exclude_gnos_sample_id_lists = args.exclude_gnos_sample_id_lists

    seq = args.seq
    workflow = args.workflow
    
    seq= list(seq) if seq else [] 
    workflow = list(workflow) if workflow else []   

    # pre-exclude gnos entries when this option is chosen
    gnos_sample_ids_to_be_excluded = generate_id_list(exclude_gnos_sample_id_lists)

    # # only process the gnos entries when this option is chosen
    gnos_sample_ids_to_be_included = generate_id_list(include_gnos_sample_id_lists) 

    # # remove the gnos_ids_to_be_include from gnos_sample_ids_to_be_excluded
    # gnos_sample_ids_to_be_excluded.difference_update(gnos_sample_ids_to_be_included) 

    # remove the gnos_sample_ids_to_be_excluded from gnos_sample_ids_to_be_included
    gnos_sample_ids_to_be_included.difference_update(gnos_sample_ids_to_be_excluded)


    es_index = 'pcawg_summary'
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

    pcawg_sample_sheet = '../pcawg-operations/lists/pc_annotation-pcawg_final_list.tsv'
    pcawg_gnos_id_sheet = '../pcawg-operations/data_releases/oct2015/release_oct2015_entry.tsv'


    ega_dir = '../pcawg-ega-submission'

    annotations = {}
    read_annotations(annotations, 'gender', ega_dir+'/annotation/donor.all_projects.release20.tsv')

    for project in dcc_project_code:
        donors_list = get_donors_list(es, es_index, project)
        gnos_sample_ids_to_be_excluded = generate_exclude_list(ega_dir, project, gnos_sample_ids_to_be_excluded)

        if seq:
            collect_sample(donors_list, gnos_sample_ids_to_be_included, gnos_sample_ids_to_be_excluded, project, ega_dir, pcawg_sample_sheet, seq, annotations)

        if workflow:
            collect_gnos_xml(donors_list, gnos_sample_ids_to_be_included, gnos_sample_ids_to_be_excluded, project, ega_dir, pcawg_gnos_id_sheet, workflow, annotations)


    return 0


if __name__ == "__main__":
    sys.exit(main())

