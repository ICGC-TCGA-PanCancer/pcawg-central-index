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

es_queries_history = [
{ 
    "name": "bwa_alignment",
    "content":
              {
                "size": 10000,
                "aggs": {
                  "published_date": {
                    "date_histogram": {"field": "published_date", "interval": "day"},
                    "aggs": {
                      "repo": {
                        "terms": {"field": "gnos_repo"}
                      }
                    }
                  }
                },
                "query":{
                  "filtered":{
                     "filter":{
                      "terms":{
                          "entity_type": [
                              "aligned_bam"
                          ] 
                      }
                    }
                  }
                }          
              }
},
{
    "name": "sanger_variant_calling",
    "content": {
        "size": 10000,
                "aggs": {
                  "published_date": {
                    "date_histogram": {"field": "published_date", "interval": "day"},
                    "aggs": {
                      "repo": {
                        "terms": {"field": "gnos_repo"}
                      }
                    }
                  }
                },
                "query": {
                  "filtered": {
                    "filter":{
                      "terms":{
                          "entity_type": [
                              "sanger_variant_calling"
                          ] 
                      }
                  }
                }
            }        
        }
},
{
    "name": "dkfz_embl_variant_calling",
    "content": {
        "size": 10000,
                "aggs": {
                  "published_date": {
                    "date_histogram": {"field": "published_date", "interval": "day"},
                    "aggs": {
                      "repo": {
                        "terms": {"field": "gnos_repo"}
                      }
                    }
                  }
                },
                "query": {
                  "filtered": {
                    "filter":{
                      "terms":{
                          "entity_type": [
                              "dkfz_embl_variant_calling"
                          ] 
                      }
                  }
                }
            }        
        }
},
{
    "name": "broad_variant_calling",
    "content": {
        "size": 10000,
                "aggs": {
                  "published_date": {
                    "date_histogram": {"field": "published_date", "interval": "day"},
                    "aggs": {
                      "repo": {
                        "terms": {"field": "gnos_repo"}
                      }
                    }
                  }
                },
                "query": {
                  "filtered": {
                    "filter":{
                      "terms":{
                          "entity_type": [
                              "broad_variant_calling"
                          ] 
                      }
                  }
                }
            }        
        }
},
{
    "name": "rna_seq_tophat_alignment",
    "content": {
        "size": 10000,
                "aggs": {
                  "published_date": {
                    "date_histogram": {"field": "published_date", "interval": "day"},
                    "aggs": {
                      "repo": {
                        "terms": {"field": "gnos_repo"}
                      }
                    }
                  }
                },
                "query": {
                  "filtered": {
                    "filter":{
                      "terms":{
                          "entity_type": [
                              "tophat"
                          ] 
                      }
                  }
                }
            }        
        }
},
{
    "name": "rna_seq_star_alignment",
    "content": {
        "size": 10000,
                "aggs": {
                  "published_date": {
                    "date_histogram": {"field": "published_date", "interval": "day"},
                    "aggs": {
                      "repo": {
                        "terms": {"field": "gnos_repo"}
                      }
                    }
                  }
                },
                "query": {
                  "filtered": {
                    "filter":{
                      "terms":{
                          "entity_type": [
                              "star"
                          ] 
                      }
                  }
                }
            }        
        }
},

]


def create_gnos_entity_info(donor_unique_id, es_json):
    gnos_entity_info_list = []

    gnos_entity_info = OrderedDict()
    gnos_entity_info['donor_unique_id'] = donor_unique_id
    gnos_entity_info['submitter_donor_id'] = es_json['submitter_donor_id']
    gnos_entity_info['dcc_project_code'] = es_json['dcc_project_code']
    
    add_wgs_gnos_entity(gnos_entity_info_list, gnos_entity_info, es_json)

    add_vcf_gnos_entity(gnos_entity_info_list, gnos_entity_info, es_json)

    add_rna_seq_gnos_entity(gnos_entity_info_list, gnos_entity_info, es_json)

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

    if aliquot.get('aligned_bam'):
        gnos_entity_info['entity_type'] = 'aligned_bam'
        gnos_entity_info['gnos_id'] = aliquot.get('aligned_bam').get('gnos_id')
        gnos_entity_info['gnos_repo'] = sort_repos_by_time(aliquot.get('aligned_bam'))[1]
        gnos_entity_info['published_date'] = sort_repos_by_time(aliquot.get('aligned_bam'))[0]
        gnos_entity_info_list.append(copy.deepcopy(gnos_entity_info))

    return gnos_entity_info_list


def sort_repos_by_time(obj):    
    published_dates = obj.get('gnos_published_date')
    gnos_repos = obj.get('gnos_repo')
    published_dates_sort, gnos_repos_sort = izip(*sorted(izip(published_dates, gnos_repos), key=lambda x: x[0]))
    return [published_dates_sort[0], gnos_repos_sort[0]]


def add_vcf_gnos_entity(gnos_entity_info_list, gnos_entity_info, es_json):
    if es_json.get('variant_calling_results'):
        gnos_entity_info['library_strategy'] = 'WGS'
        gnos_entity_info['aliquot_id'] = None
        gnos_entity_info['submitter_specimen_id'] = None
        gnos_entity_info['submitter_sample_id'] = None
        gnos_entity_info['dcc_specimen_type'] = None
        for vcf_type in ['sanger_variant_calling', 'dkfz_embl_variant_calling', 'broad_variant_calling']:
            if es_json.get('flags').get('is_'+vcf_type+'_performed') and es_json.get('variant_calling_results').get(vcf_type):
                gnos_entity_info['entity_type'] = vcf_type
                gnos_entity_info['gnos_id'] = es_json.get('variant_calling_results').get(vcf_type).get('gnos_id')    
                gnos_entity_info['gnos_repo'] = sort_repos_by_time(es_json.get('variant_calling_results').get(vcf_type))[1]
                gnos_entity_info['published_date'] = sort_repos_by_time(es_json.get('variant_calling_results').get(vcf_type))[0]
                gnos_entity_info_list.append(copy.deepcopy(gnos_entity_info))

    return gnos_entity_info_list


def add_rna_seq_gnos_entity(gnos_entity_info_list, gnos_entity_info, es_json):
    rna_seq_info = es_json.get('rna_seq').get('alignment')
    for specimen_type in rna_seq_info.keys():
        if not rna_seq_info.get(specimen_type): # the specimen_type has no alignment result
		    continue
        if 'normal' in specimen_type:
            aliquot = rna_seq_info.get(specimen_type)
            add_rna_seq_aliquot_gnos_entity(aliquot, gnos_entity_info, gnos_entity_info_list)

        else:
            for aliquot in rna_seq_info.get(specimen_type):
                add_rna_seq_aliquot_gnos_entity(aliquot, gnos_entity_info, gnos_entity_info_list)

    return gnos_entity_info_list

def add_rna_seq_aliquot_gnos_entity(aliquot, gnos_entity_info, gnos_entity_info_list):
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
        gnos_entity_info['gnos_id'] = aliquot.get(workflow_type).get('aligned_bam').get('gnos_id')
        gnos_entity_info['gnos_repo'] = sort_repos_by_time(aliquot.get(workflow_type).get('aligned_bam'))[1]
        gnos_entity_info['published_date'] = sort_repos_by_time(aliquot.get(workflow_type).get('aligned_bam'))[0]
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


def init_es(es_host, es_index):
    es = Elasticsearch([ es_host ])
    es.indices.create( es_index, ignore=400 )
    # create mappings
    es_mapping = open('pcawg_history.mapping.json')
    es.indices.put_mapping(index=es_index, doc_type='gnos_entity', body=es_mapping.read())
    es_mapping.close()
    return es

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

def generate_report(es, es_index, es_queries, report_dir):
    repos = ['bsc', 'ebi', 'cghub', 'dkfz', 'riken', 'osdc-icgc', 'osdc-tcga', 'etri']

    for q_index in range(len(es_queries)):
        counts_per_day = OrderedDict()
        
        # initialized
        start_date = datetime.date(2014, 8, 1)
        current_date = datetime.date.today()
        step = datetime.timedelta(days=1)
        all_dates = []
        while start_date<current_date:
            counts_per_day[start_date.strftime("%Y-%m-%d")]={}
            counts_per_day[start_date.strftime("%Y-%m-%d")]['count']=0
            for r in repos:
                counts_per_day[start_date.strftime("%Y-%m-%d")][r]=0
            all_dates.append(start_date.strftime("%Y-%m-%d"))
            start_date += step
        
        # get counts per day
        response = es.search(index=es_index, body=es_queries[q_index].get('content'))
        for p in response['aggregations']['published_date'].get('buckets'):
            published_date = p.get('key_as_string').split('T')[0]
            counts_per_day[published_date]['count'] = p.get('doc_count')
            for d in p.get('repo').get('buckets'):
                counts_per_day[published_date][get_formal_repo_name(d.get('key'))] = d.get('doc_count')

        # get the cumulative sum
        counts_sum = OrderedDict()
        counts_sum_list = [['Date', 'count']+repos]
                
        for d in range(len(all_dates)):
            counts_list = [all_dates[d]]
            if d == 0: 
                counts_sum[all_dates[d]] = counts_per_day.get(all_dates[d])
            else:
                counts_sum[all_dates[d]] = {}
                counts_sum[all_dates[d]]['count'] = counts_sum[all_dates[d-1]]['count'] + counts_per_day[all_dates[d]]['count']
                for r in repos:
                    counts_sum[all_dates[d]][r] = counts_sum[all_dates[d-1]][r] + counts_per_day[all_dates[d]][r]
            counts_list.append(copy.deepcopy(counts_sum[all_dates[d]]['count']))
            for r in repos:
                counts_list.append(copy.deepcopy(counts_sum[all_dates[d]][r]))
            counts_sum_list.append(copy.deepcopy(counts_list)) 
            

        # with open(report_dir + '/' + es_queries[q_index].get('name') + '.counts.json', 'w') as o:
        #     o.write(json.dumps(counts_sum))
        with open(report_dir + '/' + es_queries[q_index].get('name') + '.counts.txt', 'w') as o: 
            for r in counts_sum_list:
                o.write('\t'.join(str(x) for x in r) + '\n')


def init_report_dir(metadata_dir, report_name):
    report_dir = metadata_dir + '/reports/' + report_name
    if os.path.exists(report_dir): shutil.rmtree(report_dir, ignore_errors=True)
    os.makedirs(report_dir)
    return report_dir


def main(argv=None):

    parser = ArgumentParser(description="PCAWG Full List of GNOS entities Info Generator",
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

    fh = open(metadata_dir+'/pcawg_history_entities_'+es_index+'.jsonl', 'w')
    tsv_fh = open(metadata_dir + '/pcawg_history_entities_' + es_index + '.tsv', 'w')

    es_index_history = 'pcawg_history'
    es_history = init_es(es_host, es_index_history)
    es_type_history = "gnos_entity"

	# get the full list of donors in PCAWG
    donors_list = get_donors_list(es, es_index, es_queries)
    
    header = True
    # get json doc for each donor and reorganize it 
    for donor_unique_id in donors_list:     
        
    	es_json = get_donor_json(es, es_index, donor_unique_id)
        
        gnos_entity_info_list = create_gnos_entity_info(donor_unique_id, es_json)
        
        for gnos_entity in gnos_entity_info_list: 
            fh.write(json.dumps(gnos_entity, default=set_default) + '\n')
            # push to Elasticsearch
            es_history.index(index=es_index_history, doc_type='gnos_entity', id=gnos_entity['gnos_id'], body=json.loads(json.dumps(gnos_entity, default=set_default)), timeout=90 )

            if header:
                tsv_fh.write('\t'.join(gnos_entity.keys()) + '\n')
                header = False 
            # write to the tsv file
            for p in gnos_entity.keys():
                if isinstance(gnos_entity.get(p), set):
                    tsv_fh.write('|'.join(list(gnos_entity.get(p))) + '\t')
                elif not gnos_entity.get(p):
                    tsv_fh.write('\t')
                else:
                    tsv_fh.write(str(gnos_entity.get(p)) + '\t')
            tsv_fh.write('\n')
    tsv_fh.close()
    fh.close()

    # output result
    report_name = re.sub(r'^generate_', '', os.path.basename(__file__))
    report_name = re.sub(r'\.py$', '', report_name)
    report_dir = init_report_dir(metadata_dir, report_name)

    generate_report(es, es_index_history, es_queries_history, report_dir)


    return 0


if __name__ == "__main__":
    sys.exit(main())
