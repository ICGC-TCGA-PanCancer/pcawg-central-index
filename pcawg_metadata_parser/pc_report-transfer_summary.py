#!/usr/bin/env python

import sys
import os
import re
import json
from collections import OrderedDict
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from elasticsearch import Elasticsearch
import glob
import shutil
import copy

def get_project_donor_count(es, es_index, dcc_project_code):
    es_query_project = {
        "fields": ["donor_unique_id"], 
     
        "query": {
            "bool": {
              "must": [               
                {
                  "term": {
                    "dcc_project_code": dcc_project_code
                  }
                },
                {  
                  "term": {
                    "flags.is_normal_specimen_aligned": "T"
                  }
                },
                {
                  "term": {
                    "flags.are_all_tumor_specimens_aligned": "T"
                  }
                }
              ],
              "must_not": [
                {
                  "term": {
                    "flags.is_manual_qc_failed": "T"
                  }
                }
              ]
            }
        },
        "size": 10000
    }

    response = es.search(index=es_index, body=es_query_project)

    donors_list = set()
    for p in response['hits']['hits']:
        donors_list.add(p.get('fields').get('donor_unique_id')[0])
 
    return donors_list


def get_tumor_aliquot_count(es, es_index, donor_unique_id):
    es_query_donor = {
        "query": {
            "term": {
                "donor_unique_id": donor_unique_id
            }
        }
    }
    response = es.search(index=es_index, body=es_query_donor)

    tumor_aliquot_count = response['hits']['hits'][0]['_source']['flags']['all_tumor_specimen_aliquot_counts']
 
    return tumor_aliquot_count


def init_report_dir(metadata_dir, report_name, repo):
    report_dir = metadata_dir + '/reports/' + report_name if not repo else metadata_dir + '/reports/' + report_name + '/' + repo

    if os.path.exists(report_dir): shutil.rmtree(report_dir, ignore_errors=True)  # empty the folder if exists
    os.makedirs(report_dir)

    return report_dir


def generate_report(es, es_index, metadata_dir, report_name, timestamp, repo, jobs):
    # we need to run several queries to get facet counts for different type of donors
    report = OrderedDict()
    count_types = [
        "expected_to_be_transferred",
        "both_transferred",
        "normal_transferred_tumor_not",
        "tumor_transferred_normal_not",
        "both_not_transferred"
    ]

    for project, project_value in jobs.iteritems():
        if not report.get(project):
            report[project] = {}
        for q_index in range(len(count_types)):
            if not report[project].get(count_types[q_index]):
                report[project][count_types[q_index]] = {}
            if not report[project][count_types[q_index]].get('count'):
                report[project][count_types[q_index]]['count'] = 0
            if not report[project][count_types[q_index]].get('donor'):
                report[project][count_types[q_index]]['donors'] = set()

        donors_list = get_project_donor_count(es, es_index, project)

        report[project]['expected_to_be_transferred']['count'] = len(donors_list)
        report[project]['expected_to_be_transferred']['donors'] = copy.deepcopy(donors_list)
        report[project]['both_not_transferred']['count'] = len(donors_list)
        report[project]['both_not_transferred']['donors'] = donors_list


        for donor, donor_value in project_value.iteritems():
            donor_unique_id = project + "::" + donor
            tumor_aliquot_count = get_tumor_aliquot_count(es, es_index, donor_unique_id)

            if donor_value.get('WGS-BWA-Normal') and len(donor_value.get('WGS-BWA-Normal')) == 1:
                if donor_value.get('WGS-BWA-Tumor') and len(donor_value.get('WGS-BWA-Tumor')) == tumor_aliquot_count:
                    report[project]['both_transferred']['count'] += 1
                    report[project]['both_transferred']['donors'].add(donor_unique_id)
                    report[project]['both_not_transferred']['count'] -= 1 
                    report[project]['both_not_transferred']['donors'].discard(donor_unique_id)                    
                else:
                    report[project]['normal_transferred_tumor_not']['count'] += 1
                    report[project]['normal_transferred_tumor_not']['donors'].add(donor_unique_id)
                    report[project]['both_not_transferred']['count'] -= 1 
                    report[project]['both_not_transferred']['donors'].discard(donor_unique_id)   
            elif donor_value.get('WGS-BWA-Tumor') and len(donor_value.get('WGS-BWA-Tumor')) == tumor_aliquot_count:
                report[project]['tumor_transferred_normal_not']['count'] += 1
                report[project]['tumor_transferred_normal_not']['donors'].add(donor_unique_id)                
                report[project]['both_not_transferred']['count'] -= 1 
                report[project]['both_not_transferred']['donors'].discard(donor_unique_id)   

    report_dir = init_report_dir(metadata_dir, report_name, repo)
    
    summary_table = []
    for p in report.keys():
        summary = OrderedDict()
        summary['project'] = p
        summary['timestamp'] = timestamp
        for ctype in count_types:
            summary[ctype] = report.get(p).get(ctype).get('count') if report.get(p).get(ctype) else 0
            donors = report.get(p).get(ctype).get('donors') if report.get(p).get(ctype) else []

            if donors:
                with open(report_dir + '/' + p + '.' + ctype + '.donors.txt', 'w') as o:
                    o.write('# ' + ctype + '\n')
                    o.write('# dcc_project_code' + '\t' + 'submitter_donor_id' + '\n')
                    for d in donors:
                        # TODO: query ES to get JSON then retrieve BAM info: aligned/unaligned, gnos, bam file name etc
                        o.write(d.replace('::', '\t') + '\n')

        summary_table.append(summary)

    with open(report_dir + '/donor.json', 'w') as o:
        o.write(json.dumps(summary_table))


def main(argv=None):
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    parser = ArgumentParser(description="PCAWG Report Generator Using ES Backend",
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

    es = Elasticsearch([es_host])

    report_name_base = re.sub(r'^pc_report-', '', os.path.basename(__file__))
    report_name_base = re.sub(r'\.py$', '', report_name_base)


    # read and parse git for the gnos_ids and fnames which are completed for data transfer
    for transfer_target in ['s3', 'ceph']:
        if transfer_target == 's3':
            git_fnames = '../s3-transfer-operations/s3-transfer-jobs*/completed-jobs/*.json'
        else:
            git_fnames = '../ceph_transfer_ops/ceph-transfer-jobs*/completed-jobs/*.json'
        report_name = transfer_target + '_' + report_name_base
        files = glob.glob(git_fnames)
        jobs = dict()
        for f in files:
            fname = str.split(f, '/')[-1]
            gnos_id, dcc_project_code, donor_id, specimen_id, data_type, f_type = str.split(fname, '.')
            if not jobs.get(dcc_project_code):
                jobs[dcc_project_code] = {}
            if not jobs.get(dcc_project_code).get(donor_id):
                jobs.get(dcc_project_code)[donor_id] = {}
            if not jobs.get(dcc_project_code).get(donor_id).get(data_type):
                jobs.get(dcc_project_code).get(donor_id)[data_type] = set()       
            jobs.get(dcc_project_code).get(donor_id).get(data_type).add(gnos_id)            

        generate_report(es, es_index, metadata_dir, report_name, timestamp, repo, jobs)

    return 0


if __name__ == "__main__":
    sys.exit(main())
