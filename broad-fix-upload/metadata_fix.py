#!/usr/bin/env python


import sys
import csv
import os
import re
import xmltodict
import requests
import logging
from collections import OrderedDict
import time
import uuid
import copy
import simplejson as json
import glob
import hashlib
import subprocess
from random import randint
import shutil

logger = logging.getLogger('metadata_fix_and_merge')
# create console handler with a higher log level
ch = logging.StreamHandler()
#gnos_key = '/home/ubuntu/.ssh/gnos_key'
gnos_key = '~/.ssh/gnos_key'


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


def download_metadata_xml(gnos_id, gnos_repo, download_dir):
    metadata_xml_dir = os.path.join(download_dir, gnos_id)
    logger.info('Download metadata xml from GNOS repo: {} for analysis object: {}'.format(gnos_repo, gnos_id))
    
    url = gnos_repo + 'cghub/metadata/analysisFull/' + gnos_id
    response = None
    try:
        response = requests.get(url, stream=True, timeout=15)
    except:
        logger.error('Unable to download metadata for: {} from {}'.format(gnos_id, url))
        sys.exit('Unable to download GNOS metadata xml, please check the log for details.')

    if not response or not response.ok:
        logger.error('Unable to download metadata for: {} from {}'.format(gnos_id, url))
        sys.exit('Unable to download GNOS metadata xml, please check the log for details.')
    else:
        metadata_xml_str = response.text

        metadata_xml_file = metadata_xml_dir + '/' + gnos_id  + '.xml'
        with open(metadata_xml_file, 'w') as f:  # write to metadata xml file now
            f.write(metadata_xml_str.encode('utf8'))

def download_datafiles(gnos_id, gnos_repo, download_dir):
    url = gnos_repo + 'cghub/data/analysis/download/' + gnos_id
    # datafiles_dir = download_dir + workflow_type
    for i in range(10):
        command =   'cd {} && '.format(download_dir) + \
                    'gtdownload -c ' + gnos_key + ' ' + url

        process = subprocess.Popen(
                command,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )

        out, err = process.communicate()

        if not process.returncode:
            break
        time.sleep(randint(1,10))  # pause a few seconds before retry
    os.remove(os.path.join(download_dir, gnos_id+'.gto'))


def generate_uuid():
    uuid_str = str(uuid.uuid4())
    return uuid_str


# def update_fixed_data_block_files(data_block_files, workflow_type, upload_dir, download_dir, gnos_id):
#     lookup_dir = os.path.join(download_dir, gnos_id)
#     fname_all = []
#     data_block_files_new = []
#     for f in data_block_files:
#         fname = str(f.get('@filename'))
#         if not fname in fname_all:
#             fname_all.append(fname)
#             fname_list = str.split(fname, '.')
#             fname_list[2] = '[0-9]*[0-9]'
#             fname_search = '.'.join(fname_list) 
#             f_fixed = glob.glob(lookup_dir + '/fixed_files/' + fname_search)
#             if not f_fixed:
#                 f_old = glob.glob(lookup_dir + '/' + fname_search )
#                 if len(f_old) == 1:
#                     f_checksum = generate_md5(f_old[0]) 
#                     if not f_checksum == f.get('@checksum'):
#                         logger.warning('file: {} in the downloads folder has different checksum with the original one for analysis object: {}'.format(fname, gnos_id))
#                         return
#                 elif not f_old:
#                     logger.warning('file: {} is missing in the downloads folder for analysis object: {}'.format(fname, gnos_id))
#                     return
#                 else:
#                     logger.warning('file: {} has duplicates in the downloads folder for analysis object: {}'.format(fname, gnos_id))
#                     return
#             elif len(f_fixed) == 1:
#                 f['@filename'] = f_fixed[0].split('/')[-1]
#                 f['@checksum'] = generate_md5(f_fixed[0])
                
#             else:
#                 logger.warning('file: {} has duplicates fixed files in the downloads folder for analysis object: {}'.format(fname, gnos_id))
#                 return
#             data_block_files_new.append(copy.deepcopy(f))

#         else: # duplicates, not keep in the data_block_files_new
#             logger.warning('file: {} is duplicated and removed from the data_block_files for analysis object: {}'.format(fname, gnos_id))            
        
#     return data_block_files_new


def generate_md5(fname):
    hash = hashlib.md5()
    with open(fname) as f:
        for chunk in iter(lambda: f.read(4096), ""):
            hash.update(chunk)
    return hash.hexdigest()


def get_gnos_analysis_object(f):
    with open (f, 'r') as x: xml_str = x.read()
    analysis_xml = xmltodict.parse(xml_str).get('ResultSet').get('Result').get('analysis_xml')
    return analysis_xml


def get_fix_donor_list (list_file):
    donors_to_be_fixed = []
    with open(list_file) as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        for row in reader:
            donors_to_be_fixed.append(row)

    return donors_to_be_fixed


def validate_work_dir(work_dir, donors_to_be_fixed):
    for donor in donors_to_be_fixed:
        caller = 'broad'
        gnos_entry_dir = os.path.join(work_dir, 'downloads', donor.get(caller + '_gnos_id'))
        if not os.path.isdir(gnos_entry_dir):
            logger.error('Expected GNOS entry does not exist: {}. Please ensure all GNOS entries are downloaded.'.format(gnos_entry_dir))
            sys.exit('Validating working directory failed, please check log for details.')
        if not os.path.isdir(os.path.join(work_dir, '../fixed_files')):
            logger.error('Expected folder for BROAD fixed files does not exist: {}'.format(os.path.join(work_dir, 'fixed_files')))
            sys.exit('Validating working directory failed, please check log for details.')
        else:
            aliquot_ids = donor.get('tumor_aliquot_ids').split('|')
            for aliquot_id in aliquot_ids:
                if not os.path.exists(os.path.join(work_dir, '../fixed_files/', aliquot_id+'.oxoG.somatic.snv_mnv.vcf.gz')) or not \
                       os.path.exists(os.path.join(work_dir, '../fixed_files/', aliquot_id+'.oxoG.somatic.snv_mnv.vcf.gz.idx')):
                    logger.error('No BROAD fixed files detected in: {} for donor: {}'.format(os.path.join(work_dir, '../fixed_files'), donor.get('donor_unique_id')))
                    sys.exit('Validating working directory failed, please check log for details.')                     


def download_metadata_files(work_dir, donors_to_be_fixed):
    for donor in donors_to_be_fixed:
        caller = 'broad'
        download_datafiles(donor.get(caller + '_gnos_id'), donor.get(caller + '_gnos_repo'), os.path.join(work_dir, 'downloads'))
        download_metadata_xml(donor.get(caller + '_gnos_id'), donor.get(caller + '_gnos_repo'), os.path.join(work_dir, 'downloads'))


def metadata_fix(work_dir, donors_to_be_fixed):
    caller = 'broad'
    for donor in donors_to_be_fixed:
        gnos_analysis_objects = {}

        upload_gnos_uuid = generate_uuid()
        upload_dir = os.path.join(work_dir, 'uploads', get_formal_repo_name(donor.get(caller + '_gnos_repo')), upload_gnos_uuid)
        if os.path.isdir(upload_dir): # this should never happen, but if happen regenerate a new UUID
            upload_gnos_uuid = generate_uuid()
            upload_dir = os.path.join(work_dir, 'uploads', get_formal_repo_name(donor.get(caller + '_gnos_repo')), upload_gnos_uuid)
        os.makedirs(upload_dir)


        gnos_entry_dir = os.path.join(work_dir, 'downloads', donor.get(caller + '_gnos_id'))
        fixed_file_dir = os.path.join(work_dir, '../fixed_files')
        xml_file = os.path.join(gnos_entry_dir, donor.get(caller + '_gnos_id') + '.xml')

        gnos_analysis_object = get_gnos_analysis_object(xml_file)
        files = get_files(gnos_analysis_object, fixed_file_dir, gnos_entry_dir)
        create_symlinks(upload_dir, files)
        apply_data_block_patches(gnos_analysis_object, files)

        create_fixed_gnos_submission(upload_dir, gnos_analysis_object) # merge xml

        generate_updated_metadata(donor, upload_gnos_uuid)


def generate_updated_metadata(donor, upload_gnos_uuid):
    pass
    


def create_fixed_gnos_submission(upload_dir, gnos_analysis_object):
    fixed_analysis_object = copy.deepcopy(gnos_analysis_object) 

    attributes = fixed_analysis_object.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_ATTRIBUTES').get('ANALYSIS_ATTRIBUTE')
    for attr in attributes:
        if attr.get('TAG') == 'workflow_file_subset':
            attr['VALUE'] = 'broad-v2'

        else:
            pass  # all others, leave it unchanged

    new_analysis_xml_str = xmltodict.unparse(fixed_analysis_object, pretty=True)
    # print new_analysis_xml  # debug only
    new_analysis_xml_file = os.path.join(upload_dir, 'analysis.xml')
    with open(new_analysis_xml_file, 'w') as f:  f.write(new_analysis_xml_str.encode('utf8'))


def get_attr(gnos_attr):
    attributes = {}
    for attr in gnos_attr: attributes[attr.get('TAG')] = attr.get('VALUE')
    return attributes


def apply_data_block_patches(gnos_analysis_object, files):
    old_files = {}
    for f in gnos_analysis_object.get('ANALYSIS_SET').get('ANALYSIS').get('DATA_BLOCK').get('FILES').get('FILE'):
        old_files.update({
                f.get('@filename'): {
                    'filetype': f.get('@filetype'),
                    'checksum_method': f.get('@checksum_method'),
                    'checksum': f.get('@checksum')
                }
            })

    file_names = [os.path.basename(f) for f in files]
    new_files = []
    for i in xrange(len(file_names)):
        new_file = OrderedDict()
        filename = file_names[i]
        realfile = files[i]

        filetype = None
        checksum_method = None
        checksum = None

        if not '/../fixed_files/' in realfile: # this is one of the old files, use previous values
            filetype = old_files.get(filename).get('filetype')
            checksum_method = old_files.get(filename).get('checksum_method')
            checksum = old_files.get(filename).get('checksum')
        else:  # fixed file
            if filename.endswith('vcf.gz'):
                filetype = 'vcf'
            elif filename.endswith('gz.idx'):
                filetype = 'idx'
            else:
                logger.warning('Unrecognized file type for file: {}'.format('filename'))
                continue

            checksum_method = 'MD5'
            checksum = generate_md5(realfile)
            if not checksum:  # should never happen
                logger.warning('Unable to generate md5sum for the fixed file: {}'.format(realfile))
                continue

        new_file.update({'@filename': filename})
        new_file.update({'@filetype': filetype})
        new_file.update({'@checksum_method': checksum_method})
        new_file.update({'@checksum': checksum})

        new_files.append(new_file)

    # assign the new files list to DATA_BLOCK
    gnos_analysis_object.get('ANALYSIS_SET').get('ANALYSIS').get('DATA_BLOCK').get('FILES')['FILE'] = new_files


def create_symlinks(target, source):
    for s in source:
        os.symlink(s, os.path.join(target, os.path.basename(s)))


def get_files(gnos_analysis_object, fixed_file_dir, gnos_entry_dir):
    file_name_patterns = set([
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.germline\.indel\.vcf\.gz$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.germline\.indel\.vcf\.gz\.idx$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.indel\.vcf\.gz$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.indel\.vcf\.gz\.idx$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.broad-dRanger\..+\.somatic\.sv\.vcf\.gz$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.broad-dRanger\..+\.somatic\.sv\.vcf\.gz\.idx$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.broad-snowman\..+\.somatic\.sv\.vcf\.gz$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.broad-snowman\..+\.somatic\.sv\.vcf\.gz\.idx$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.broad-dRanger_snowman\..+\.somatic\.sv\.vcf\.gz$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.broad-dRanger_snowman\..+\.somatic\.sv\.vcf\.gz\.idx$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.germline\.sv\.vcf\.gz$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.germline\.sv\.vcf\.gz\.idx$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.snv_mnv\.vcf\.gz$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.snv_mnv\.vcf\.gz\.idx$'
        ])

    matched_files = []
    for file_dir in (fixed_file_dir, gnos_entry_dir): # match fixed_file dir first
        for file in glob.glob(os.path.join(file_dir, "*")):
            file_name = os.path.basename(file)
            matched_fp = None
            for fp in file_name_patterns:
                if re.match(fp, file_name):
                    matched_fp = fp
                    matched_files.append(file)

            if matched_fp: file_name_patterns.remove(matched_fp)  # remove the file pattern that had a match
 
    print matched_files

    # print len(matched_files) 
    for fp in file_name_patterns:
        logger.error('Missing expected variant call result file with pattern: {}'.format(fp))
        sys.exit('Missing expected variant call result file, see log file for details.')

    return matched_files

def generate_index_files(work_dir, donors_to_be_fixed):
    pass


def detect_folder(work_dir, folder_name):
    # dectect whether uploads dir exists, stop if exists
    job_dir = os.path.join(work_dir, folder_name)
    if os.path.isdir(job_dir):
        try:
            os.rmdir(job_dir)
        except OSError as ex:
            sys.exit('\nStop: none empty job directory exists: {}. Please confirm it\'s safe to remove, then manually remove it and try this script again.\n'.format(job_dir))

    os.mkdir(job_dir)


def main():
    if len(sys.argv) == 1: sys.exit('\nMust specify working directory where the variant call fixes are kept.\nPlease refer to the SOP for details how to structure the working directory.\n')
    work_dir = sys.argv[1]
    # if os.path.isdir(work_dir): shutil.rmtree(work_dir, ignore_errors=True)  # empty the folder if exists
    # os.makedirs(work_dir+'/downloads/broad')
    work_dir = os.path.abspath(work_dir)

    if work_dir == os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test'):
        print('\nUsing \'test\' folder as working directory ...')
        test = True
        donor_list_file = 'test_donor_list.txt'
    else:
        donor_list_file = 'to_be_fixed_donor_list.txt'

    if not os.path.exists(donor_list_file): sys.exit('Helper file missing!')
    donors_to_be_fixed = get_fix_donor_list (donor_list_file)


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

    # now download metadata xml
    # if not test:
    print('\nDownloading data files and metadata XML from GNOS ...')
    # download_metadata_files(work_dir, donors_to_be_fixed)

    print('\nGenerating index files for fixed_files ...')
    # generate_index_files(work_dir, donors_to_be_fixed)

    # print('\nCreate symlink for fixed snv_mnv files ...')
    # attach_fixed_files(work_dir, donors_to_be_fixed)

    # validate working direcotry first
    print('\nValidating working directory...')
    validate_work_dir(work_dir, donors_to_be_fixed)

    # dectect whether uploads dir exists, stop if exists
    detect_folder(work_dir, 'uploads')
    detect_folder(work_dir, 'updated_metafiles')
    # upload_dir = os.path.join(work_dir, 'uploads')
    # if os.path.isdir(upload_dir):
    #     try:
    #         os.rmdir(upload_dir)
    #     except OSError as ex:
    #         sys.exit('\nStop: none empty "uploads" directory exists: {}. Please confirm it\'s safe to remove, then manually remove it and try this script again.\n'.format(upload_dir))

    # os.mkdir(upload_dir)

    # now process metadata xml fix and merge
    print('\nPreparing new GNOS submissions and updated related metadata XML files...')
    metadata_fix(work_dir, donors_to_be_fixed)

    # print('\nGenerating the updated related metadata XML files...')
    # generate_updated_metadata(work_dir, donors_to_be_fixed)
    
    print('Submission folder located at: {}'.format(os.path.join(work_dir, 'uploads')))
    print('Processing log file: {}'.format(log_file))
    print('Done!\n')


if __name__ == "__main__":
    sys.exit(main())