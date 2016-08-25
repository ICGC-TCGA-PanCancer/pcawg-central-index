#!/usr/bin/env python


import sys
import csv
import os
import re
import xmltodict
import requests
import logging
from collections import OrderedDict
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import time
import uuid
import copy
import simplejson as json
import glob
import hashlib
import subprocess
from random import randint
import shutil
import tarfile

logger = logging.getLogger('pcawg_final_consensus_vcfs_upload')

ch = logging.StreamHandler()



def get_files(call, work_dir, tumor_aliquot_ids):

    matched_files = []
    vcf_file_dir = os.path.join(work_dir, call)

    for aliquot in tumor_aliquot_ids:
        file_name_patterns = set([
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.consensus\.'+re.escape(call)+r'\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.consensus\.'+re.escape(call)+r'\.vcf\.gz\.tbi$'
            ])


        for f in glob.glob(os.path.join(vcf_file_dir, aliquot+'*')):
            file_name = os.path.basename(f)
            # print file_name
            matched_fp = None
            for fp in file_name_patterns:
                if re.match(fp, file_name):
                    matched_fp = fp
                    matched_files.append(copy.deepcopy(f))

            if matched_fp: file_name_patterns.remove(matched_fp)  # remove the file pattern that had a match
    
    # print matched_files
    # print len(matched_files)
    if file_name_patterns:
        for fp in file_name_patterns:
            logger.error('Missing expected consensus variant call result file with pattern: {} for aliquot {}'.format(fp, aliquot))
     
    return matched_files



def create_results_copies(row, create_results_copy, work_dir):
    
    donor_id = row.get('icgc_donor_id')
    tumor_aliquot_ids = row.get('tumor_wgs_aliquot_id').split('|')
    for dt in create_results_copy:
        call_results_dir = os.path.join(work_dir,'call_results_dir', dt, donor_id)
        if not os.path.isdir(call_results_dir): os.makedirs(call_results_dir)
        vcf_files = get_files(dt, work_dir, tumor_aliquot_ids)      
        copy_files(call_results_dir, vcf_files)
        generate_md5_files(call_results_dir, tumor_aliquot_ids)
        

def generate_md5_files(folder_name, tumor_aliquot_ids):
    for aliquot_id in tumor_aliquot_ids:
        for f in glob.glob(os.path.join(folder_name, aliquot_id+'*')):
            if f.endswith('md5'): continue
            md5_value = generate_md5(f)
            with open(f+'.md5', 'w') as fh: fh.write(md5_value)


def generate_md5(fname):
    hash = hashlib.md5()
    with open(fname) as f:
        for chunk in iter(lambda: f.read(4096), ""):
            hash.update(chunk)
    return hash.hexdigest()


def copy_files(target, source):
    for s in source:
        filename = os.path.basename(s).replace('consensus', 'consensus.20160826.somatic')
        shutil.copy(s, os.path.join(target, filename))

def generate_uuid():
    uuid_str = str(uuid.uuid4())
    return uuid_str
    

def generate_analysis_xmls(row, generate_analysis_xml, work_dir):
    donor_id = row.get('icgc_donor_id')
    aliquot_id = row.get('tumor_wgs_aliquot_id')
    normal_bam_url = row.get('normal_wgs_bwa_alignment_gnos_repo').split('|')[0] + 'cghub/metadata/analysisFull/' + row.get('normal_wgs_bwa_alignment_gnos_id')
    tumor_bam_url = row.get('tumor_wgs_bwa_alignment_gnos_repo').split('|')[0] + 'cghub/metadata/analysisFull/' + row.get('tumor_wgs_bwa_alignment_gnos_id')
    metadata_urls = normal_bam_url + ',' + tumor_bam_url
    project_code = row.get('dcc_project_code')
 
    for dt in generate_analysis_xml:
        call_results_dir = os.path.join(work_dir,'call_results_dir', dt, donor_id)
        vcf_files = glob.glob(os.path.join(call_results_dir, aliquot_id+'*.vcf.gz'))
        workflow_name = 'consensus_' + dt
        gnos_id = generate_uuid()
        output_dir = os.path.join(dt, 'osdc-tcga') if project_code.endswith('-US') else os.path.join(dt, 'osdc-icgc')
        study_ref_name = 'tcga_pancancer_vcf' if project_code.endswith('-US') else 'icgc_pancancer_vcf'

        command = generate_perl_command(dt, gnos_id, metadata_urls, vcf_files, output_dir, study_ref_name)

        process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        out, err = process.communicate()

        if not process.returncode:#success
            continue
        else:
            sys.exit('Error:{}'.format(err))



def generate_perl_command(call, gnos_id, metadata_urls, vcf_files, workflow_file_subset, related_file_subset_uuids, output_dir, study_ref_name):

    command =   'perl -I /home/ubuntu/gt-download-upload-wrapper/lib/ /home/ubuntu/vcf-uploader/gnos_upload_vcf.pl' +\
                 ' --key gnos_fake_key '+\
                 ' --metadata-urls ' + metadata_urls +\
                 ' --vcfs ' + ','.join(vcf_files) +\
                 ' --vcf-idxs ' + ','.join([vcf+'.tbi' for vcf in vcf_files]) +\
                 ' --vcf-md5sum-files ' + ','.join([vcf+'.md5' for vcf in vcf_files]) +\
                 ' --vcf-idx-md5sum-files ' + ','.join([vcf+'.tbi.md5' for vcf in vcf_files]) +\
                 ' --workflow-name ' + workflow_name +\
                 ' --study-refname-override icgc_pancancer_vcf ' +\
                 ' --workflow-version 1.0.0 ' +\
                 ' --workflow-src-url https://github.com/ucscCancer/pcawg_tools ' +\
                 ' --workflow-url https://github.com/ucscCancer/pcawg_tools ' +\
                 ' --skip-validate --skip-upload ' +\
                 ' --uuid ' + gnos_id +\
                 ' --outdir ' + output_dir +\
                 ' --study-refname-override ' + study_ref_name   

    return command


def generate_id_list(id_lists):
    ids_list = []
    if id_lists:
        with open(id_lists, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
            for row in reader: ids_list.append(row.get('Submitter_donor_ID'))
    return ids_list



def main(argv=None):
    parser = ArgumentParser(description="PCAWG final consensus vcf files upload",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--work_dir", dest="work_dir",
             help="Directory name containing fixed variant call files", required=True)
    parser.add_argument("-v", "--vcf_info_file", dest="vcf_info_file",
             help="vcf information file", required=False)
    parser.add_argument("-c", "--create_results_copy", dest="create_results_copy", nargs="*",
             help="create results copy for given variant call", required=False)
    parser.add_argument("-g", "--generate_analysis_xml", dest="generate_analysis_xml", nargs="*",
             help="generate analysis xml for given variant call", required=False)
    parser.add_argument("-i", "--include_donor_id_lists", dest="include_donor_id_lists", nargs="*",
             help="indicate DONOR IDs to be included", required=False)

    args = parser.parse_args()
    work_dir = args.work_dir
    vcf_info_file = args.vcf_info_file 
    create_results_copy = args.create_results_copy 
    generate_analysis_xml = args.generate_analysis_xml 
    include_donor_id_lists = args.include_donor_id_lists

    if not vcf_info_file: vcf_info_file = 'synapse_table_160403.tsv'
    if not os.path.exists(vcf_info_file): sys.exit('information file is missing')

    create_results_copy = list(create_results_copy) if create_results_copy else []
    generate_analysis_xml= list(generate_analysis_xml) if generate_analysis_xml else []
    include_donor_id_lists= list(include_donor_id_lists) if include_donor_id_lists else generate_id_list(vcf_info_file)   
 

    if not os.path.isdir(work_dir+'/call_results_dir'): os.makedirs(work_dir+'/call_results_dir')

    current_time = time.strftime("%Y-%m-%d_%H-%M-%S")

    logger.setLevel(logging.INFO)
    ch.setLevel(logging.ERROR)

    log_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), current_time + '.process.log')

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)

    with open(vcf_info_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        for row in reader:
            if not row.get('Submitter_donor_ID') in include_donor_id_lists: continue

            if create_results_copy: create_results_copies(row, create_results_copy, work_dir)

            if generate_analysis_xml: generate_analysis_xmls(row, generate_analysis_xml, work_dir)


if __name__ == "__main__":
    sys.exit(main())
