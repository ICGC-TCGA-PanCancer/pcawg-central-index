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

logger = logging.getLogger('RNA-Seq_metadata_merge')
# create console handler with a higher log level
ch = logging.StreamHandler()

def download_metadata_xml(gnos_id, gnos_repo, workflow_type, download_dir, merged_gnos_id):
    metadata_xml_dir = download_dir + workflow_type + '/' + merged_gnos_id
    logger.info('Download metadata xml from GNOS repo: {} for analysis object: {}'.format(gnos_repo, gnos_id))
    
    url = get_formal_repo_name(gnos_repo) + 'cghub/metadata/analysisFull/' + gnos_id
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


def generate_md5(fname):
    hash = hashlib.md5()
    with open(fname) as f:
        for chunk in iter(lambda: f.read(4096), ""):
            hash.update(chunk)
    return hash.hexdigest()


def get_gnos_object(f, obj, key, new_name, unaligned_merged_gnos_id):
    with open (f, 'r') as x: xml_str = x.read()
    for k in key:
        # replace the old filename with new filename
        old_name = new_name.replace(unaligned_merged_gnos_id, k)
        xml_str = re.sub('"'+old_name+'"', '"'+new_name+'"', xml_str)
    obj_xml = xmltodict.parse(xml_str).get('ResultSet').get('Result').get(obj+'_xml')
    return obj_xml


def get_fix_donor_list (list_file):
    donors_to_be_fixed = []
    with open(list_file) as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        for row in reader:
            donors_to_be_fixed.append(row)

    return donors_to_be_fixed


def validate_work_dir(work_dir, donors_to_be_fixed):
    for donor in donors_to_be_fixed:
        if not donor.get('count_match') == 'TRUE': continue
        for caller in ('STAR', 'TopHat2'):
            gnos_entry_dir = os.path.join(work_dir, 'downloads', caller, donor.get(caller + '_merged_gnos_id'))
            if not os.path.isdir(gnos_entry_dir):
                logger.error('Expected merged RNA-Seq data does not exist: {}.'.format(gnos_entry_dir))
                sys.exit('Validating working directory failed, please check log for details.')


def download_metadata_files(work_dir, donors_to_be_fixed):
    for donor in donors_to_be_fixed:
        if not donor.get('count_match') == 'TRUE': continue
        for caller in ('STAR', 'TopHat2'):
            gnos_ids_to_merge = donor.get(caller + '_analysis_id').split('|')
            merged_gnos_id = donor.get(caller + '_merged_gnos_id')
            gnos_repo = donor.get(caller + '_gnos_server').split('|')[0]
            for gnos_id in gnos_ids_to_merge:
                download_metadata_xml(gnos_id, gnos_repo, caller, os.path.join(work_dir, 'downloads/'), merged_gnos_id)


def metadata_fix_and_merge(work_dir, donors_to_be_fixed):
    for donor in donors_to_be_fixed:
        if not donor.get('count_match') == 'TRUE': continue
        donor_aliquot_id = donor.get('aliquot_id')
        for caller in ('STAR', 'TopHat2'):
            upload_gnos_uuid = donor.get(caller+'_merged_gnos_id')
            upload_gnos_repo = 'osdc-icgc' if not donor.get('dcc_project_code').endswith('-US') else 'osdc-tcga'
            upload_dir = os.path.join(work_dir, 'uploads', caller, upload_gnos_repo, upload_gnos_uuid)
            os.makedirs(upload_dir)
            merged_filename = donor.get(caller+'_merged_file_name')
            unaligned_gnos_id = donor.get('unaligned_analysis_id').split('|')
            unaligned_merged_gnos_id = donor.get('unaligned_merged_gnos_id')
            gnos_entry_dir = os.path.join(work_dir, 'downloads', caller, upload_gnos_uuid)
            # create symlinks for bam and bai files to upload
            create_symlinks(upload_dir, set(glob.glob(gnos_entry_dir+'/*.bam') + glob.glob(gnos_entry_dir+'/*.bam.bai')))
            # loop over all the gnos objects
            for obj in ('analysis', 'experiment', 'run'):
                gnos_objects = {}
                # loop over all the lanes
                for lane_gnos_id in donor.get(caller + '_analysis_id').split('|'):                    
                    xml_file = os.path.join(gnos_entry_dir, lane_gnos_id + '.xml')
                    gnos_object = get_gnos_object(xml_file, obj, unaligned_gnos_id, merged_filename, unaligned_merged_gnos_id)
                    gnos_objects.update({lane_gnos_id: gnos_object})

                create_merged_gnos_submission(donor_aliquot_id, caller, upload_dir, gnos_objects, obj, gnos_entry_dir) # merge xml


def create_merged_gnos_submission(donor_aliquot_id, caller, upload_dir, gnos_objects, obj, gnos_entry_dir):
    # random choose the xml of one lane as starting point, the first element in gnos_objects
    lane_object = gnos_objects.itervalues().next()
    merged_object = copy.deepcopy(lane_object)

    if obj == 'analysis':
        apply_data_block_patches(merged_object, set(glob.glob(gnos_entry_dir+'/*.bam') + glob.glob(gnos_entry_dir+'/*.bam.bai')))
        RUN_LABELS = []
        NOTES = []
        ATTR = {}

        for k, v in gnos_objects.iteritems():
            # merge RUN_LABELS
            RUN_LABELS.append(v.get('ANALYSIS_SET').get('ANALYSIS')['ANALYSIS_TYPE']['REFERENCE_ALIGNMENT']['RUN_LABELS'])

            # merge PIPELINE NOTES
            NOTES.append({'NOTE': v.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_TYPE').get('REFERENCE_ALIGNMENT').get('PROCESSING').get('PIPELINE').get('PIPE_SECTION').get('NOTES')})
            
            # deal with the attributes to add all the info into set
            attributes = v.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_ATTRIBUTES').get('ANALYSIS_ATTRIBUTE')
            for attr in attributes:
                if attr.get('TAG') in ('STUDY',
                                       'workflow_name',
                                       'workflow_version',
                                       'workflow_source_url',
                                       'workflow_bundle_url',
                                       caller.upper()+'_version',
                                       'dcc_project_code',
                                       'submitter_donor_id',
                                       'submitter_sample_id',
                                       'submitter_specimen_id',
                                       'dcc_specimen_type',
                                       ):
                    if not ATTR.get(attr.get('TAG')): ATTR[attr.get('TAG')] = set()
                    ATTR[attr.get('TAG')].add(attr.get('VALUE'))
                    #check whether all of the metadata have the same info
                    if not len(ATTR.get(attr.get('TAG'))) == 1:
                        logger.warning('RNA-Seq lanes metadata of donor: {} for caller: {} have discrepancy in attributes: {} '.format(donor_aliquot_id, caller, attr.get('TAG'))) 

        merged_object.get('ANALYSIS_SET').get('ANALYSIS')['ANALYSIS_TYPE']['REFERENCE_ALIGNMENT']['RUN_LABELS'] = RUN_LABELS
        merged_object.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_TYPE').get('REFERENCE_ALIGNMENT').get('PROCESSING').get('PIPELINE').get('PIPE_SECTION')['NOTES'] = NOTES
        
    
    elif obj == 'experiment':
        EXPERIMENT = []
        for k, v in gnos_objects.iteritems():
            EXPERIMENT.append(v.get('EXPERIMENT_SET'))
        merged_object['EXPERIMENT_SET'] = EXPERIMENT
    elif obj == 'run':
        RUN = []
        for k, v in gnos_objects.iteritems():
            RUN.append(v.get('RUN_SET'))
        merged_object['RUN_SET'] = RUN
    else:
        pass

    new_xml_str = xmltodict.unparse(merged_object, pretty=True)
    # print new_analysis_xml  # debug only
    new_xml_file = os.path.join(upload_dir, obj+'.xml')
    with open(new_xml_file, 'w') as f:  f.write(new_xml_str.encode('utf8'))


def get_attr(gnos_attr):
    attributes = {}
    for attr in gnos_attr: attributes[attr.get('TAG')] = attr.get('VALUE')
    return attributes


def apply_data_block_patches(gnos_analysis_object, real_files):
    file_names = [os.path.basename(f) for f in real_files]
    new_files = []
    for i in xrange(len(file_names)):
        new_file = OrderedDict()
        filename = file_names[i]
        realfile = real_files[i]

        filetype = None
        checksum_method = None
        checksum = None

        if filename.endswith('bam'):
            filetype = 'bam'
        elif filename.endswith('bam.bai'):
            filetype = 'bai'
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


def main():
    if len(sys.argv) == 1: sys.exit('\nMust specify working directory where the variant call fixes are kept.\nPlease refer to the SOP for details how to structure the working directory.\n')
    work_dir = sys.argv[1]
    if not os.path.isdir(work_dir):
        sys.exit('Specified working directory does not exist.')
    work_dir = os.path.abspath(work_dir)

    if work_dir == os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test'):
        print('\nUsing \'test\' folder as working directory ...')
        test = True
        donor_list_file = 'test_donor_list.txt'
    else:
        test = False
        donor_list_file = 'rnaseq_to_merge_lanes_derived_from_same_aliquot.tsv'

    if not os.path.exists(donor_list_file): sys.exit('Helper file missing!')
    donors_to_be_fixed = get_fix_donor_list(donor_list_file)

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

    # validate working direcotry first
    print('\nValidating working directory...')
    validate_work_dir(work_dir, donors_to_be_fixed)

    # dectect whether uploads dir exists, stop if exists
    upload_dir = os.path.join(work_dir, 'uploads')
    if os.path.isdir(upload_dir):
        try:
            os.rmdir(upload_dir)
        except OSError as ex:
            sys.exit('\nStop: none empty "uploads" directory exists: {}. Please confirm it\'s safe to remove, then manually remove it and try this script again.\n'.format(upload_dir))

    os.mkdir(upload_dir)

    # now download metadata xml
    if not test:
        logger.info('Downloading GNOS metadata XML for unmerged RNA-Seq data...')
        download_metadata_files(work_dir, donors_to_be_fixed)

    # now process metadata xml fix and merge
    print('Preparing new GNOS submissions...')
    metadata_fix_and_merge(work_dir, donors_to_be_fixed)

    print('Submission folder located at: {}'.format(os.path.join(work_dir, 'uploads')))
    print('Processing log file: {}'.format(log_file))
    print('Done!\n')


if __name__ == "__main__":
    sys.exit(main())