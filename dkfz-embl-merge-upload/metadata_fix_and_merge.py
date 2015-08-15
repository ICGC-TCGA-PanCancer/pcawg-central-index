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

logger = logging.getLogger('metadata_fix_and_merge')
# create console handler with a higher log level
ch = logging.StreamHandler()

def download_metadata_xml(gnos_id, gnos_repo, workflow_type, download_dir):
    metadata_xml_dir = download_dir + workflow_type + '/' + gnos_id
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


def generate_uuid():
    uuid_str = str(uuid.uuid4())
    return uuid_str


def update_fixed_data_block_files(data_block_files, workflow_type, upload_dir, download_dir, gnos_id):
    lookup_dir = download_dir + workflow_type + '/' + gnos_id
    fname_all = []
    data_block_files_new = []
    for f in data_block_files:
        fname = str(f.get('@filename'))
        if not fname in fname_all:
            fname_all.append(fname)
            fname_list = str.split(fname, '.')
            fname_list[2] = '[0-9]*[0-9]'
            fname_search = '.'.join(fname_list) 
            f_fixed = glob.glob(lookup_dir + '/fixed_files/' + fname_search)
            if not f_fixed:
                f_old = glob.glob(lookup_dir + '/' + fname_search )
                if len(f_old) == 1:
                    f_checksum = generate_md5(f_old[0]) 
                    if not f_checksum == f.get('@checksum'):
                        logger.warning('file: {} in the downloads folder has different checksum with the original one for analysis object: {}'.format(fname, gnos_id))
                        return
                elif not f_old:
                    logger.warning('file: {} is missing in the downloads folder for analysis object: {}'.format(fname, gnos_id))
                    return
                else:
                    logger.warning('file: {} has duplicates in the downloads folder for analysis object: {}'.format(fname, gnos_id))
                    return
            elif len(f_fixed) == 1:
                f['@filename'] = f_fixed[0].split('/')[-1]
                f['@checksum'] = generate_md5(f_fixed[0])
                
            else:
                logger.warning('file: {} has duplicates fixed files in the downloads folder for analysis object: {}'.format(fname, gnos_id))
                return
            data_block_files_new.append(copy.deepcopy(f))

        else: # duplicates, not keep in the data_block_files_new
            logger.warning('file: {} is duplicated and removed from the data_block_files for analysis object: {}'.format(fname, gnos_id))            
        
    return data_block_files_new


def generate_md5(file):
    with open (file, 'r') as x: data = x.read()
    md5 = hashlib.md5(data).hexdigest()
    return md5


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
        for caller in ('dkfz', 'embl'):
            gnos_entry_dir = os.path.join(work_dir, 'downloads', caller, donor.get(caller + '_gnos_id'))
            if not os.path.isdir(gnos_entry_dir):
                logger.error('Expected GNOS entry does not exist: {}. Please ensure all GNOS entries are downloaded.'.format(gnos_entry_dir))
                sys.exit('Validating working directory failed, please check log for details.')
            if caller == 'dkfz':
                if not os.path.isdir(os.path.join(gnos_entry_dir, 'fixed_files')):
                    logger.error('Expected folder for DKFZ fixed files does not exist: {}'.format(os.path.join(gnos_entry_dir, 'fixed_files')))
                    sys.exit('Validating working directory failed, please check log for details.')
                elif not set(glob.glob(gnos_entry_dir+'/fixed_files/*.gz') + glob.glob(gnos_entry_dir+'/fixed_files/*.gz.tbi')):  # no fixed file
                    logger.error('No DKFZ fixed files detected in: {}'.format(os.path.join(gnos_entry_dir, 'fixed_files')))
                    sys.exit('Validating working directory failed, please check log for details.')


def download_metadata_files(work_dir, donors_to_be_fixed):
    for donor in donors_to_be_fixed:
        for caller in ('dkfz', 'embl'):
            download_metadata_xml(donor.get(caller + '_gnos_id'),\
                donor.get(caller + '_gnos_repo'), caller, os.path.join(work_dir, 'downloads/'))


def metadata_fix_and_merge(work_dir, donors_to_be_fixed):
    for donor in donors_to_be_fixed:
        gnos_analysis_objects = {}
        dkfz_files = []

        upload_gnos_uuid = generate_uuid()
        upload_dir = os.path.join(work_dir, 'uploads', upload_gnos_uuid)
        if os.path.isdir(upload_dir): # this should never happen, but if happen regenerate a new UUID
            upload_gnos_uuid = generate_uuid()
            upload_dir = os.path.join(work_dir, 'uploads', upload_gnos_uuid)
        os.mkdir(upload_dir)

        for caller in ('dkfz', 'embl'):
            gnos_entry_dir = os.path.join(work_dir, 'downloads', caller, donor.get(caller + '_gnos_id'))
            fixed_file_dir = os.path.join(gnos_entry_dir, 'fixed_files')
            xml_file = os.path.join(gnos_entry_dir, donor.get(caller + '_gnos_id') + '.xml')

            gnos_analysis_object = get_gnos_analysis_object(xml_file)
            if caller == 'dkfz':
                dkfz_files = get_dkfz_files(gnos_analysis_object, fixed_file_dir)
                create_symlinks(upload_dir, dkfz_files)
                apply_dkfz_data_block_patches(gnos_analysis_object, dkfz_files)

                #print xmltodict.unparse(gnos_analysis_object, pretty=True)  # debug only
            else: # embl
                create_symlinks(upload_dir, set(glob.glob(gnos_entry_dir+'/*.gz') + glob.glob(gnos_entry_dir+'/*.gz.tbi')))

            gnos_analysis_objects.update({caller: gnos_analysis_object})

        create_merged_gnos_submission(upload_dir, gnos_analysis_objects) # merge xml


def create_merged_gnos_submission(upload_dir, gnos_analysis_objects):
    embl_analysis_object = gnos_analysis_objects.get('embl')
    dkfz_analysis_object = gnos_analysis_objects.get('dkfz')
    merged_analysis_object = copy.deepcopy(embl_analysis_object) # EMBL is the first part of the two calling workflow, hence choosing it as starting point

    merged_analysis_object.get('ANALYSIS_SET').get('ANALYSIS')['DESCRIPTION'] = \
        '[Description EMBL]: ' + embl_analysis_object.get('ANALYSIS_SET').get('ANALYSIS')['DESCRIPTION'] + \
        ' [Description DKFZ]: ' + dkfz_analysis_object.get('ANALYSIS_SET').get('ANALYSIS')['DESCRIPTION']

    # merge PIPE_SECTION
    dkfz_analysis_object.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_TYPE').get('REFERENCE_ALIGNMENT').get('PROCESSING').get('PIPELINE').get('PIPE_SECTION')['STEP_INDEX'] = 2
    merged_analysis_object.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_TYPE').get('REFERENCE_ALIGNMENT').get('PROCESSING').get('PIPELINE')['PIPE_SECTION'] = \
        [embl_analysis_object.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_TYPE').get('REFERENCE_ALIGNMENT').get('PROCESSING').get('PIPELINE')['PIPE_SECTION'], \
            dkfz_analysis_object.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_TYPE').get('REFERENCE_ALIGNMENT').get('PROCESSING').get('PIPELINE')['PIPE_SECTION']]

    # merge FILES
    merged_analysis_object.get('ANALYSIS_SET').get('ANALYSIS').get('DATA_BLOCK').get('FILES')['FILE'] = \
        embl_analysis_object.get('ANALYSIS_SET').get('ANALYSIS').get('DATA_BLOCK').get('FILES').get('FILE') + \
            dkfz_analysis_object.get('ANALYSIS_SET').get('ANALYSIS').get('DATA_BLOCK').get('FILES').get('FILE')

    # deal with attributes now
    embl_attributes = get_attr(embl_analysis_object.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_ATTRIBUTES').get('ANALYSIS_ATTRIBUTE'))
    dkfz_attributes = get_attr(dkfz_analysis_object.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_ATTRIBUTES').get('ANALYSIS_ATTRIBUTE'))

    attributes = merged_analysis_object.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_ATTRIBUTES').get('ANALYSIS_ATTRIBUTE')
    for attr in attributes:
        if attr.get('TAG') == 'variant_workflow_name':
            attr['VALUE'] = 'DKFZ_EMBL_Merged'

        elif attr.get('TAG') in ('dcc_project_code',
                               'submitter_donor_id',
                               'PM',
                               'alignment_workflow_name',
                               'alignment_workflow_source_url',
                               'alignment_workflow_version',
                               'alignment_workflow_bundle_url',
                               'variant_pipeline_input_info',
                               ):
            continue  # use the previous value no need to update anything

        elif attr.get('TAG') in ('variant_qc_metrics', 'variant_timing_metrics'):
            attr['VALUE'] = merge_metrics_jsons(
                attr.get('TAG'),
                embl_attributes.get(attr.get('TAG')),
                dkfz_attributes.get(attr.get('TAG'))
            )

        elif attr.get('TAG') in ('variant_pipeline_output_info'):
            attr['VALUE'] = merge_pipeline_output_info_jsons(
                embl_attributes.get(attr.get('TAG')),
                dkfz_attributes.get(attr.get('TAG'))
            )

        elif attr.get('TAG') in (  # strings
                'variant_workflow_version',
                'variant_workflow_source_url',
                'variant_workflow_bundle_url',
                'vm_instance_type',
                'vm_location_code',
                ):
            if embl_attributes.get(attr.get('TAG')) == dkfz_attributes.get(attr.get('TAG')):
                attr['VALUE'] = embl_attributes.get(attr.get('TAG'))
            else:
                attr['VALUE'] = 'EMBL: ' + embl_attributes.get(attr.get('TAG')) + \
                                ' DKFZ: ' + dkfz_attributes.get(attr.get('TAG'))

        elif attr.get('TAG') in (  # numbers
                'vm_instance_cores',
                'vm_instance_mem_gb',
                ):
            embl_v = embl_attributes.get(attr.get('TAG')) if isfloat(embl_attributes.get(attr.get('TAG'))) else ''
            dkfz_v = dkfz_attributes.get(attr.get('TAG')) if isfloat(dkfz_attributes.get(attr.get('TAG'))) else ''

            if isfloat(embl_v) and isfloat(dkfz_v):
                attr['VALUE'] = embl_v if float(embl_v) > float(dkfz_v) else dkfz_v
            elif isfloat(embl_v) and not isfloat(dkfz_v):
                attr['VALUE'] = embl_v
            elif not isfloat(embl_v) and isfloat(dkfz_v):
                attr['VALUE'] = dkfz_v
            else:
                attr['VALUE'] = ''

        else:
            pass  # all others, leave it unchanged

    new_analysis_xml_str = xmltodict.unparse(merged_analysis_object, pretty=True)
    # print new_analysis_xml  # debug only
    new_analysis_xml_file = os.path.join(upload_dir, 'analysis.xml')
    with open(new_analysis_xml_file, 'w') as f:  f.write(new_analysis_xml_str.encode('utf8'))


def isfloat(value):
  try:
    float(value)
    return True
  except:
    return False


def merge_metrics_jsons(metrics, embl_jsons_str, dkfz_jsons_str):
    metrics = metrics.replace('variant_', '')

    try:
        embl_attributes = json.loads(embl_jsons_str)
    except:
        embl_attributes = {}

    try:
        dkfz_attributes = json.loads(dkfz_jsons_str)
    except:
        dkfz_attributes = {}

    merged_metrics = {
        metrics:{
            "embl": embl_attributes.get(metrics),
            "dkfz": dkfz_attributes.get(metrics),
            "_note_": "information kept here may not be accurate for DKFZ as its results have been patched"
        }
    }
    return json.dumps(merged_metrics)


def merge_pipeline_output_info_jsons(embl_jsons_str, dkfz_jsons_str):
    try:
        embl_attributes = json.loads(embl_jsons_str)
    except:
        embl_attributes = {}

    try:
        dkfz_attributes = json.loads(dkfz_jsons_str)
    except:
        dkfz_attributes = {}

    merged_output_info = {
        'workflow_outputs':{
            "embl": embl_attributes,
            "dkfz": dkfz_attributes,
            "_note_": "information kept here may not be accurate for DKFZ as its results have been patched"
        }
    }
    return json.dumps(merged_output_info)


def get_attr(gnos_attr):
    attributes = {}
    for attr in gnos_attr: attributes[attr.get('TAG')] = attr.get('VALUE')
    return attributes


def apply_dkfz_data_block_patches(gnos_analysis_object, dkfz_files):
    old_files = {}
    for f in gnos_analysis_object.get('ANALYSIS_SET').get('ANALYSIS').get('DATA_BLOCK').get('FILES').get('FILE'):
        old_files.update({
                f.get('@filename'): {
                    'filetype': f.get('@filetype'),
                    'checksum_method': f.get('@checksum_method'),
                    'checksum': f.get('@checksum')
                }
            })

    file_names = [os.path.basename(f) for f in dkfz_files]
    new_files = []
    for i in xrange(len(file_names)):
        new_file = OrderedDict()
        filename = file_names[i]
        realfile = dkfz_files[i]

        filetype = None
        checksum_method = None
        checksum = None

        if '/fixed_files/../' in realfile: # this is one of the old files, use previous values
            filetype = old_files.get(filename).get('filetype')
            checksum_method = old_files.get(filename).get('checksum_method')
            checksum = old_files.get(filename).get('checksum')
        else:  # fixed file
            if filename.endswith('vcf.gz'):
                filetype = 'vcf'
            elif filename.endswith('gz.tbi'):
                filetype = 'idx'
            elif filename.endswith('tar.gz'):
                filetype = 'tar'
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


def get_dkfz_files(gnos_analysis_object, fixed_file_dir):
    file_name_patterns = set([
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.germline\.indel\.vcf\.gz$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.germline\.indel\.vcf\.gz\.tbi$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.indel\.vcf\.gz$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.indel\.vcf\.gz\.tbi$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.indel\.tar\.gz$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.cnv\.vcf\.gz$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.cnv\.vcf\.gz\.tbi$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.cnv\.tar\.gz$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.germline\.snv_mnv\.vcf\.gz$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.germline\.snv_mnv\.vcf\.gz\.tbi$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.snv_mnv\.vcf\.gz$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.snv_mnv\.vcf\.gz\.tbi$',
            r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.snv_mnv\.tar\.gz$',

        ])

    matched_files = []
    for file_dir in (fixed_file_dir, os.path.join(fixed_file_dir, '..')): # match fixed_file dir first
        for file in glob.glob(os.path.join(file_dir, "*")):
            file_name = os.path.basename(file)
            matched_fp = None
            for fp in file_name_patterns:
                if re.match(fp, file_name):
                    matched_fp = fp
                    matched_files.append(file)

            if matched_fp: file_name_patterns.remove(matched_fp)  # remove the file pattern that had a match

    for fp in file_name_patterns:
        if fp == '^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.germline\.indel\.vcf\.gz\.tbi$':
            logger.warning('Missing file with pattern: {}'.format(fp))
        else:
            logger.error('Missing expected variant call result file with pattern: {}'.format(fp))
            sys.exit('Missing expected variant call result file, see log file for details.')

    return matched_files


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
        logger.info('Downloading GNOS metadata XML...')
        download_metadata_files(work_dir, donors_to_be_fixed)

    # now process metadata xml fix and merge
    print('Preparing new GNOS submissions...')
    metadata_fix_and_merge(work_dir, donors_to_be_fixed)

    print('Submission folder located at: {}'.format(os.path.join(work_dir, 'uploads')))
    print('Processing log file: {}'.format(log_file))
    print('Done!\n')


if __name__ == "__main__":
    sys.exit(main())