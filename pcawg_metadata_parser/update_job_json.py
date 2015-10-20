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


def read_annotations(annotations, type, file_name):
    if not os.path.isfile(file_name):
        return
    with open(file_name, 'r') as r:
        if type in ['icgc_donor_id', 'icgc_sample_id', 'icgc_specimen_id']:
            annotations[type] = {}
            subtype = type.split('_')[1]
            prefix = subtype[0:2]
            for line in r:
                if line.startswith('#'): continue
                if len(line.rstrip()) == 0: continue
                icgc_id, id_pcawg, dcc_project_code, creation_release = str.split(line.rstrip(), ',')
                annotations[type][dcc_project_code+'::'+id_pcawg] = prefix.upper()+icgc_id 

        else:
            logger.warning('unknown annotation type: {}'.format(type))

def main(argv=None):

    parser = ArgumentParser(description="S3 Transfer Jobs Json Generator",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=True)
    parser.add_argument("-f", "--specify updated json folder", dest="updated_json_dir",
             help="Specify updated json folder", required=False)


    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports
    updated_json_dir = args.updated_json_dir

    if not os.path.isdir(metadata_dir):  # TODO: should add more directory name check to make sure it's right
        sys.exit('Error: specified metadata directory does not exist!')

    if not updated_json_dir:
    	updated_json_dir = metadata_dir + '/' + re.sub(r'\.py$', '', os.path.basename(__file__))
    
    if os.path.exists(updated_json_dir): shutil.rmtree(updated_json_dir, ignore_errors=True)  # empty the folder if exists
    os.makedirs(updated_json_dir)

    annotations = {}
    read_annotations(annotations, 'icgc_donor_id', 'pc_annotation-icgc_donor_ids.csv')
    read_annotations(annotations, 'icgc_specimen_id', 'pc_annotation-icgc_specimen_ids.csv')
    read_annotations(annotations, 'icgc_sample_id', 'pc_annotation-icgc_sample_ids.csv')

    bam_fh = open(updated_json_dir+'/bam_job.jsonl', 'w')
    sanger_fh = open(updated_json_dir+'/sanger_job.jsonl', 'w')

    # read and parse git for the gnos_ids and fnames which are scheduled for s3 transfer
    s3_git_fnames = '../s3-transfer-operations/s3-transfer-jobs*/completed-jobs/*.json'
    s3_files = glob.glob(s3_git_fnames)
    ceph_git_fnames = '../ceph_transfer_ops/ceph-transfer-jobs*/*/*.json'
    ceph_files = glob.glob(ceph_git_fnames)
    files = s3_files + ceph_files
    for f in files:
    	with open(f, 'r') as fh:
    		fname = str.split(f, '/')[-1]
    		json_obj = json.load(fh)
    		project_code = json_obj.get('project_code')
    		submitter_donor_id = json_obj.get('submitter_donor_id')
    		submitter_sample_id = json_obj.get('submitter_sample_id')
    		submitter_specimen_id = json_obj.get('submitter_specimen_id')
    		json_obj['icgc_donor_id'] = annotations.get('icgc_donor_id').get(project_code+'::'+submitter_donor_id) if annotations.get('icgc_donor_id').get(project_code+'::'+submitter_donor_id) else None
    		if 'Sanger' in f:
    			with open(updated_json_dir+'/'+fname, 'w') as w:
    				w.write(json.dumps(json_obj, indent=4, sort_keys=True))
    			sanger_fh.write(json.dumps(json_obj) + '\n')
    			continue
    		json_obj['icgc_specimen_id'] = annotations.get('icgc_specimen_id').get(project_code+'::'+submitter_specimen_id) if annotations.get('icgc_specimen_id').get(project_code+'::'+submitter_specimen_id) else None
    		json_obj['icgc_sample_id'] = annotations.get('icgc_sample_id').get(project_code+'::'+submitter_sample_id) if annotations.get('icgc_sample_id').get(project_code+'::'+submitter_sample_id) else None
    		with open(updated_json_dir+'/'+fname, 'w') as w:
    			w.write(json.dumps(json_obj, indent=4, sort_keys=True))
    		bam_fh.write(json.dumps(json_obj) + '\n')

    bam_fh.close()
    sanger_fh.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())