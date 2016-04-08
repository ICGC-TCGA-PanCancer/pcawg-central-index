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

logger = logging.getLogger('metadata_fix_and_upload')
# create console handler with a higher log level
ch = logging.StreamHandler()
#gnos_key = '/home/ubuntu/.ssh/gnos_key'
# gnos_key = '~/.ssh/gnos_key'


def get_files(donor_id, call, work_dir, aliquot_id):

    matched_files = []
    if call == 'muse':
        file_name_patterns = set([
                r'^.+\.somatic\.snv_mnv\.vcf\.gz$',
                # r'^.+\.somatic\.snv_mnv\.vcf\.gz\.md5$',
                r'^.+\.somatic\.snv_mnv\.vcf\.gz\.idx$'
                # r'^.+\.somatic\.snv_mnv\.vcf\.gz\.idx\.md5$'
            ])

        file_dir = 'Muse-calls'
        for f in glob.glob(os.path.join(work_dir, file_dir, donor_id+'*')):
            file_name = os.path.basename(f)
            # print file_name
            matched_fp = None
            for fp in file_name_patterns:
                if re.match(fp, file_name):
                    matched_fp = fp
                    matched_files.append(copy.deepcopy(f))

            if matched_fp: file_name_patterns.remove(matched_fp)  # remove the file pattern that had a match

        if file_name_patterns:
            for fp in file_name_patterns:
                logger.error('Missing expected variant call result file with pattern: {} for aliquot {}'.format(fp, aliquot)) 

    elif call == 'broad-v3':
        file_name_patterns = set([
                r'^.+\.germline\.indel\.vcf\.gz$',
                r'^.+\.germline\.indel\.vcf\.gz\.idx$',
                r'^.+\.somatic\.indel\.vcf\.gz$',
                r'^.+\.somatic\.indel\.vcf\.gz\.idx$',
                r'^.+\.broad-dRanger[^_].+\.somatic\.sv\.vcf\.gz$',
                r'^.+\.broad-dRanger[^_].+\.somatic\.sv\.vcf\.gz\.idx$',
                r'^.+\.broad-snowman.+\.somatic\.sv\.vcf\.gz$',
                r'^.+\.broad-snowman.+\.somatic\.sv\.vcf\.gz\.idx$',
                r'^.+\.broad-dRanger_snowman.+\.somatic\.sv\.vcf\.gz$',
                r'^.+\.broad-dRanger_snowman.+\.somatic\.sv\.vcf\.gz\.idx$',
                r'^.+\.germline\.sv\.vcf\.gz$',
                r'^.+\.germline\.sv\.vcf\.gz\.idx$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.snv_mnv\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.snv_mnv\.vcf\.gz\.idx$'
            ])
        file_dir = 'broad-fix-for-long-running-jobs'
        for f in glob.glob(os.path.join(work_dir, file_dir, aliquot_id+'*')):
            file_name = os.path.basename(f)
            # print file_name
            matched_fp = None
            for fp in file_name_patterns:
                if re.match(fp, file_name):
                    matched_fp = fp
                    matched_files.append(copy.deepcopy(f))

            if matched_fp: file_name_patterns.remove(matched_fp)  # remove the file pattern that had a match

        file_dir = 'Broad-calls/'+donor_id+'/links_for_gnos/tabix_*'
        for f in glob.glob(os.path.join(work_dir, file_dir, donor_id+'*')):
            file_name = os.path.basename(f)
            # print file_name
            matched_fp = None
            for fp in file_name_patterns:
                if re.match(fp, file_name):
                    matched_fp = fp
                    matched_files.append(copy.deepcopy(f))

            if matched_fp: file_name_patterns.remove(matched_fp)  # remove the file pattern that had a match

        if file_name_patterns:
            for fp in file_name_patterns:
                logger.warning('Missing expected variant call result file with pattern: {} for aliquot {}'.format(fp, aliquot_id))       


    elif call == 'broad_tar':
        file_dir = 'Broad-calls/'+donor_id+'/links_for_broad'
        matched_files = glob.glob(os.path.join(work_dir, file_dir, '*'))
        return matched_files            

    else:
        pass
     
    return matched_files


def create_results_copys(row, create_results_copy, work_dir):
    # with open(vcf_info_file, 'r') as f:
    #     reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
    #     for row in reader:
    donor_id = row.get('Submitter_donor_ID')
    aliquot_id = row.get('Tumour_WGS_aliquot_IDs')
    for dt in create_results_copy:
        call_results_dir = os.path.join(work_dir,'call_results_dir', dt, donor_id)
        if not os.path.isdir(call_results_dir): os.makedirs(call_results_dir)
        # if dt=='muse':
        vcf_files = get_files(donor_id, dt, work_dir, aliquot_id)
        print vcf_files
        #sys.exit(0)

        # else:
        #     pass        

        copy_files(call_results_dir, vcf_files, donor_id, aliquot_id, dt)
        generate_md5_files(call_results_dir, aliquot_id)
        

def generate_md5_files(folder_name, aliquot_id):
    for f in glob.glob(os.path.join(folder_name, aliquot_id+'*')):
        md5_value = generate_md5(f)
        with open(f+'.md5', 'w') as fh: fh.write(md5_value)


def generate_md5(fname):
    hash = hashlib.md5()
    with open(fname) as f:
        for chunk in iter(lambda: f.read(4096), ""):
            hash.update(chunk)
    return hash.hexdigest()


def copy_files(target, source, donor_id, aliquot_id, call):
    if call=='muse':
        for s in source:
            filename = os.path.basename(s).replace(donor_id, aliquot_id)
            shutil.copy(s, os.path.join(target, filename))

    elif call=='broad-v3':
        for s in source:
            filename = os.path.basename(s).replace(donor_id, aliquot_id).replace('DATECODE', '20160401')
            shutil.copy(s, os.path.join(target, filename))

    elif call=='broad_tar':
        filename = aliquot_id+'.broad.intermediate.tar'
        with tarfile.open(os.path.join(target, filename), "w") as tar:
            for s in source:
                tar.add(s) 

    else:
        pass
    

def generate_analysis_xmls(row, generate_analysis_xml, work_dir):
    donor_id = row.get('Submitter_donor_ID')
    aliquot_id = row.get('Tumour_WGS_aliquot_IDs')
    normal_bam_url = row.get('Normal_WGS_alignment_GNOS_repos').split('|')[0] + 'cghub/metadata/analysisFull/' + row.get('Normal_WGS_alignment_GNOS_analysis_ID')
    tumor_bam_url = row.get('Tumour_WGS_alignment_GNOS_repos').split('|')[0] + 'cghub/metadata/analysisFull/' + row.get('Tumour_WGS_alignment_GNOS_analysis_IDs')
    metadata_urls = normal_bam_url + ',' + tumor_bam_url
    project_code = row.get('Project_code')
 
    for dt in generate_analysis_xml:
        call_results_dir = os.path.join(work_dir,'call_results_dir', dt, donor_id)
        vcf_files = glob.glob(os.path.join(call_results_dir, aliquot_id+'*.vcf.gz')) if dt in ['muse', 'broad-v3'] else glob.glob(os.path.join(call_results_dir, aliquot_id+'.broad.intermediate.tar'))
        workflow_file_subset = dt
        gnos_id = row.get(id_mapping(dt))
        related_file_subset_uuids = [row.get('Muse_VCF_UUID'), row.get('Broad_VCF_UUID'), row.get('Broad_TAR_UUID')]
        related_file_subset_uuids.remove(row.get(id_mapping(dt)))
        output_dir = os.path.join(dt, 'osdc-tcga') if project_code.endswith('-US') else os.path.join(dt, 'osdc-icgc') 

        command = generate_perl_command(gnos_id, metadata_urls, vcf_files, workflow_file_subset, related_file_subset_uuids, output_dir)

#        print command
#        sys.exit(0)

        process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        out, err = process.communicate()
 #       print err
 #       sys.exit(0)
        if not process.returncode:#success
            continue
        else:
            sys.exit('Error:{}'.format(err))


def id_mapping(vcf):
    vcf_map = {
      "broad-v3": "Broad_VCF_UUID",
      "muse": "Muse_VCF_UUID",
      "broad_tar": "Broad_TAR_UUID"
    }   

    return vcf_map.get(vcf)


def generate_perl_command(gnos_id, metadata_urls, vcf_files, workflow_file_subset, related_file_subset_uuids, output_dir):

    command =   'perl -I /home/ubuntu/gt-download-upload-wrapper/lib/ /home/ubuntu/vcf-uploader/gnos_upload_vcf.pl' +\
                 ' --key gnos_fake_key '+\
                 ' --metadata-urls ' + metadata_urls +\
                 ' --vcfs ' + ','.join(vcf_files) +\
                 ' --vcf-idxs ' + ','.join([vcf+'.idx' for vcf in vcf_files]) +\
                 ' --vcf-md5sum-files ' + ','.join([vcf+'.md5' for vcf in vcf_files]) +\
                 ' --vcf-idx-md5sum-files ' + ','.join([vcf+'.idx.md5' for vcf in vcf_files]) +\
                 ' --workflow-name BROAD_MUSE_PIPELINE ' +\
                 ' --study-refname-override icgc_pancancer_vcf ' +\
                 ' --workflow-version 1.0.0 ' +\
                 ' --workflow-src-url https://github.com/ucscCancer/pcawg_tools ' +\
                 ' --workflow-url https://github.com/ucscCancer/pcawg_tools ' +\
                 ' --skip-validate --skip-upload ' +\
                 ' --workflow-file-subset ' + workflow_file_subset +\
                 ' --related-file-subset-uuids ' + ','.join(related_file_subset_uuids) +\
                 ' --uuid ' + gnos_id +\
                 ' --outdir ' + output_dir
    return command


def generate_id_list(id_lists):
    ids_list = []
    if id_lists:
        with open(id_lists, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
            for row in reader: ids_list.append(row.get('Submitter_donor_ID'))
    return ids_list



def main(argv=None):
    parser = ArgumentParser(description="Broad snv file fix and upload",
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
    if not os.path.exists(vcf_info_file): sys.exit('Helper file is missing')

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

    # pre-exclude donors when this option is chosen
    if not include_donor_id_lists:
        donor_ids_to_be_included = generate_id_list(include_donor_id_lists)

    with open(vcf_info_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        for row in reader:
            if not row.get('Submitter_donor_ID') in include_donor_id_lists: continue

            if create_results_copy: create_results_copys(row, create_results_copy, work_dir)

            if generate_analysis_xml: generate_analysis_xmls(row, generate_analysis_xml, work_dir)


if __name__ == "__main__":
    sys.exit(main())
