#!/usr/bin/env python


import sys
import os
import re
import glob
import json
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from collections import OrderedDict
import datetime
import dateutil.parser
import shutil
import logging
from elasticsearch import Elasticsearch
import subprocess


logger = logging.getLogger('pcawg objects generator')
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
                "flags.exists_vcf_file_prefix_mismatch": [
                  "T"
                ]
              }
            },
            # {
            #    "terms":{
            #       "flags.exists_xml_md5sum_mismatch":[
            #          "T"
            #       ]
            #    }
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


def generate_id_list(id_lists):
    ids_list = set()
    if id_lists:
        files = glob.glob(id_lists)
        for fname in files:
            with open(fname) as f:
                for d in f: ids_list.add(d.rstrip())

    return ids_list


def get_mapping_transfer_object(annotations):
    # read and parse git for the gnos_ids and fnames which are scheduled for s3 transfer
    # s3_git_fnames = '../s3-transfer-operations/s3-transfer-jobs*/*/*.json'
    # s3_files = glob.glob(s3_git_fnames)
    ceph_git_fnames = '../ceph_transfer_ops/ceph-transfer-jobs*/*/*.json'
    ceph_files = glob.glob(ceph_git_fnames)
    # files = s3_files + ceph_files
    files = ceph_files
    fname_set = set()
    annotations['transfer_objects_map'] = {}
    for f in files:
    	fname = str.split(f, '/')[-1]
        if fname in fname_set: continue
        fname_set.add(fname)
        with open(f, 'r') as fh:
            json_obj = json.load(fh)
            gnos_id = json_obj.get('gnos_id')
            for obj in json_obj.get('files'):
                annotations['transfer_objects_map'][gnos_id+'::'+obj.get('file_name')] = obj.get('object_id')
                annotations['transfer_objects_map'][obj.get('object_id')] = gnos_id+'::'+obj.get('file_name')
    return annotations


def get_qc_object(cloud_bucket, annotations):
    # read and parse git for the gnos_ids and fnames which are scheduled for s3 transfer
    if cloud_bucket in ['aws', 'aws_public']:
        fnames = '../s3-data-qc-ops/bam_qc/*/*.json'
    elif cloud_bucket in ['ceph', 'ceph_public']:
    	fnames = '../ceph_transfer_ops/ceph-qc-jobs*/*/*.json'
    else:
    	return annotations

    files = glob.glob(fnames)
    annotations[cloud_bucket+'_qc'] = {}
    for f in files:
    	fname = str.split(f, '/')[-1]
    	gnos_id = str.split(fname, '.')[0]
    	qc_status = str.split(f, '/')[-2].split('-')[0]
        annotations[cloud_bucket+'_qc'][gnos_id] = qc_status

    return annotations


def get_bucket_url(cloud_bucket):
    bucket_url = {
	    'aws_public': 's3://oicr.icgc.meta/metadata/',
	    'aws': 's3://oicr.icgc/data/',
	    'ceph_public': 's3://oicr.icgc.meta/metadata/',
	    'ceph': 's3://oicr.icgc/data/'
	}
    return bucket_url.get(cloud_bucket)


def get_cloud_object_list(cloud_bucket, annotations):
    list_file = cloud_bucket + '.ls.out'
    if cloud_bucket in ['aws', 'aws_public']:		
        command = 'aws --profile amazon s3 ls ' + get_bucket_url(cloud_bucket) + ' > ' + list_file
    elif cloud_bucket in ['ceph', 'ceph_public']:
        command = 'aws --profile collabo --endpoint-url https://www.cancercollaboratory.org:9080 s3 ls ' + get_bucket_url(cloud_bucket) + ' > ' + list_file
    else:
    	sys.exit('Unknown cloud bucket: {}'.format(cloud_bucket))

    process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

    out, err = process.communicate()
    if process.returncode:
        # should not exit for just this error, improve it later
        sys.exit('Unable to download file from {}.\nError message: {}'.format(cloud_bucket, err))
    
    with open(list_file, 'r') as l:
    	annotations[cloud_bucket+'_objects']={}
    	for line in l:
            if '.meta' in line: continue
            if len(line.rstrip()) == 0: continue
            fields = str.split(line.rstrip())
            if not len(fields) == 4: continue
            annotations[cloud_bucket+'_objects'][fields[-1]] = fields[-2]

    return annotations

def get_donors_list(es, es_index, es_queries):
    q_index = 0
    response = es.search(index=es_index, body=es_queries[q_index])
    donors_list = set()
    for p in response['hits']['hits']:
        donors_list.add(p.get('fields').get('donor_unique_id')[0])
    return donors_list 


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


def get_cloud_object(donor_unique_id, es_json, object_fh, annotations, cloud_buckets, object_id_set):
    cloud_object = {
        'donor_unique_id': donor_unique_id,
        'submitter_donor_id': es_json['submitter_donor_id'],
        'dcc_project_code': es_json['dcc_project_code'],
        'gnos':{},
        'object_id': None,
        'qc':{},       
        'release':{},
    }
    object_id_set = add_cloud_object_info(cloud_object, es_json, object_fh, annotations, cloud_buckets, object_id_set)
    return object_id_set


def create_cloud_mapping(cloud_object, cloud_bucket, annotations):
    object_id = cloud_object.get('object_id')
    object_size = cloud_object.get('gnos').get('file_size')
    if annotations.get(cloud_bucket+'_objects'):
        if not object_id:
            cloud_object['is_'+cloud_bucket+'_transferred'] = False 
            cloud_object[cloud_bucket] = {
            'object_id': None,
            'object_id_match': None,
            'object_size_match': None
            }
            return cloud_object
        if annotations.get(cloud_bucket+'_objects').get(object_id):
            cloud_object['is_'+cloud_bucket+'_transferred'] = True
            cloud_object[cloud_bucket] = {
            'object_id': object_id,
            'object_id_match': True,
            'object_size_match': True if annotations.get(cloud_bucket+'_objects').get(object_id) == object_size else False
            }
        else:
            cloud_object['is_'+cloud_bucket+'_transferred'] = None
            logger.warning('Object: {} is missing on cloud bucket: {}'.format(object_id, get_bucket_url(cloud_bucket)))
    return cloud_object

	
def create_cloud_object(gnos_id, file_name, file_size, cloud_object, obj, annotations, cloud_buckets, object_fh, object_id_set):		
    cloud_object['gnos'] = {
	   'gnos_id': gnos_id,
	   'file_name': file_name,
	   'file_size': file_size,
	   'file_status': 'live'
    }
    cloud_object['object_id'] = annotations['transfer_objects_map'][gnos_id+'::'+file_name] if annotations['transfer_objects_map'].get(gnos_id+'::'+file_name) else None
    if cloud_object['object_id']: object_id_set.add(cloud_object['object_id'])

    for cloud_bucket in cloud_buckets:
        cloud_object['qc'][cloud_bucket+'_qc_status'] = annotations[cloud_bucket+'_qc'][gnos_id] if annotations[cloud_bucket+'_qc'].get(gnos_id) else None
    for release in ['oct2015']:
        cloud_object['release']['is_'+release+'_entry'] =  obj.get('is_'+release+'_entry')
    for cloud_bucket in cloud_buckets:
        cloud_object = create_cloud_mapping(cloud_object, cloud_bucket, annotations)
    
    # push to Elasticsearch
    # es.index(index=es_index_object, doc_type='object', body=json.loads(json.dumps(cloud_object, default=set_default)), timeout=90 )
    object_fh.write(json.dumps(cloud_object, default=set_default) + '\n')
    return object_id_set


def create_alignment_object(cloud_object, alignment, annotations, object_fh, cloud_buckets, object_id_set):
    for file_type in ['bai', 'bam']:
        if not alignment.get(file_type+'_file_name'): continue
        gnos_id = alignment.get('gnos_id')	    
        file_name = alignment.get(file_type+'_file_name')
        file_size = alignment.get(file_type+'_file_size')
        object_id_set = create_cloud_object(gnos_id, file_name, file_size, cloud_object, alignment, annotations, cloud_buckets, object_fh, object_id_set)
    object_id_set = create_xml_object(cloud_object, alignment, annotations, object_fh, cloud_buckets, object_id_set)
    return object_id_set

def create_vcf_object(cloud_object, vcf, annotations, object_fh, cloud_buckets, object_id_set):
    for key, value in vcf.iteritems():
        if value.get('is_stub'): continue
        if not value.get('files'): continue
        cloud_object['object_type'] = value.get('vcf_workflow_type')
        gnos_id = value.get('gnos_id')
        for f in value.get('files'):
            file_name = f.get('file_name')
            file_size = f.get('file_size')
            object_id_set = create_cloud_object(gnos_id, file_name, file_size, cloud_object, value, annotations, cloud_buckets, object_fh, object_id_set)		
        object_id_set = create_xml_object(cloud_object, value, annotations, object_fh, cloud_buckets, object_id_set)
    return object_id_set

def create_xml_object(cloud_object, obj, annotations, object_fh, cloud_buckets, object_id_set):
    gnos_id = obj.get('gnos_id')
    file_name = gnos_id+'.xml'
    file_size = None
    object_id_set = create_cloud_object(gnos_id, file_name, file_size, cloud_object, obj, annotations, cloud_buckets, object_fh, object_id_set)
    return object_id_set

def add_cloud_object_info(cloud_object, es_json, object_fh, annotations, cloud_buckets, object_id_set):
    if es_json.get('normal_alignment_status') and es_json.get('normal_alignment_status').get('aligned_bam'):
        alignment = es_json.get('normal_alignment_status').get('aligned_bam')
        cloud_object['object_type'] = 'bwa_alignment'
        object_id_set = create_alignment_object(cloud_object, alignment, annotations, object_fh, cloud_buckets, object_id_set)
        
    if es_json.get('tumor_alignment_status'):
    	for aliquot in es_json.get('tumor_alignment_status'):
    	    if aliquot.get('aligned_bam'):		
                alignment = aliquot.get('aligned_bam')
                cloud_object['object_type'] = 'bwa_alignment'
                object_id_set = create_alignment_object(cloud_object, alignment, annotations, object_fh, cloud_buckets, object_id_set)
                
    if es_json.get('variant_calling_results'):
    	vcf = es_json.get('variant_calling_results')
        object_id_set = create_vcf_object(cloud_object, vcf, annotations, object_fh, cloud_buckets, object_id_set)
    return object_id_set  

def create_append_objects(object_fh, annotations, cloud_bucket, object_id_set, object_id):
    cloud_object = {
        'donor_unique_id': None,
        'submitter_donor_id': None,
        'dcc_project_code': None,
        'object_id': None,
        'gnos':{},
        'qc':{},       
        'release':{}
    }
    
    cloud_object[cloud_bucket] = {
        'object_id': object_id,
        'object_id_match': None,
        'object_size_match': None
    }
    if annotations['transfer_objects_map'].get(object_id):
        cloud_object['object_id'] = object_id
        gnos_id, file_name = annotations['transfer_objects_map'].get(object_id).split('::')
        cloud_object['gnos']['gnos_id'] =  gnos_id
        cloud_object['gnos']['file_name'] = file_name
        cloud_object['gnos']['file_size'] =  None
        cloud_object['qc'][cloud_bucket+'_qc_status'] = annotations[cloud_bucket+'_qc'][gnos_id] if annotations[cloud_bucket+'_qc'].get(gnos_id) else None
    object_fh.write(json.dumps(cloud_object, default=set_default) + '\n')    

def set_default(obj):
    if isinstance(obj, datetime.datetime):
        return obj.isoformat()
    if isinstance(obj, set):
        return list(obj)
    raise TypeError      



def main(argv=None):

    parser = ArgumentParser(description="PCAWG Reorganized Json Donors Info Generator",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=True)   
    parser.add_argument("-r", "--gnos_repo", dest="repo",
             help="Specify which GNOS repo to process, process all repos if none specified", required=False)
    parser.add_argument("-x", "--exclude_donor_id_lists", dest="exclude_donor_id_lists", 
             help="File(s) containing DONOR IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-i", "--include_donor_id_lists", dest="include_donor_id_lists", 
             help="File(s) containing DONOR IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-v", "--variant_calling", dest="vcf", nargs="*",
             help="List variant_calling types", required=False)


    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports
    repo = args.repo
    exclude_donor_id_lists = args.exclude_donor_id_lists
    include_donor_id_lists = args.include_donor_id_lists
    vcf = args.vcf

    vcf = list(vcf) if vcf else []   

    if not os.path.isdir(metadata_dir):  # TODO: should add more directory name check to make sure it's right
        sys.exit('Error: specified metadata directory does not exist!')

    timestamp = str.split(metadata_dir, '/')[-1]
    es_index = 'p_' + ('' if not repo else repo+'_') + re.sub(r'\D', '', timestamp).replace('20','',1)
    es_type = "donor"
    es_host = 'localhost:9200'
    es = Elasticsearch([es_host])

    es_index_object = 'c_' + ('' if not repo else repo+'_') + re.sub(r'\D', '', timestamp).replace('20','',1)


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


    object_fh = open(metadata_dir+'/reports/pcawg_objects.jsonl', 'w')

    # pre-exclude donors when this option is chosen
    donor_ids_to_be_excluded = generate_id_list(exclude_donor_id_lists)

	# get the list of donors
    # only process the gnos entries when this option is chosen
    donor_ids_to_be_included = generate_id_list(include_donor_id_lists)
    if not donor_ids_to_be_included:  
        donors_list = get_donors_list(es, es_index, es_queries)
    else:
        donors_list = donor_ids_to_be_included

    # exclude the donors if they were specified on the exclude_donor_id_lists
    donors_list.difference_update(donor_ids_to_be_excluded)

    donors_list = sorted(donors_list) 
    print len(donors_list) 
    
    cloud_buckets = ['ceph', 'ceph_public']
    annotations = {}
    get_mapping_transfer_object(annotations)
    for cloud_bucket in cloud_buckets: #['aws','aws_public', 'ceph', 'ceph_public']:
        get_cloud_object_list(cloud_bucket, annotations)
        get_qc_object(cloud_bucket, annotations)

    print annotations.keys()

    object_id_set = set()
    for donor_unique_id in donors_list:
        es_json = get_donor_json(es, es_index, donor_unique_id)
        object_id_set = get_cloud_object(donor_unique_id, es_json, object_fh, annotations, cloud_buckets, object_id_set)

    for cloud_bucket in cloud_buckets:
        if not annotations.get(cloud_bucket+'_objects'): continue
        for k in annotations.get(cloud_bucket+'_objects').keys():
            if k in object_id_set: continue
            create_append_objects(object_fh, annotations, cloud_bucket, object_id_set, k)

    object_fh.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())

