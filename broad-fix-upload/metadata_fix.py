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

logger = logging.getLogger('metadata_fix_and_upload')
# create console handler with a higher log level
ch = logging.StreamHandler()
#gnos_key = '/home/ubuntu/.ssh/gnos_key'
# gnos_key = '~/.ssh/gnos_key'


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


def download_metadata_xml(gnos_id, gnos_repo, download_dir=None):
    logger.info('Download metadata xml from GNOS repo: {} for analysis object: {}'.format(gnos_repo, gnos_id))
    
    url = gnos_repo + 'cghub/metadata/analysisFull/' + gnos_id

    if download_dir:
        job_dir = os.path.join(download_dir, gnos_id)
    else:
        job_dir = os.path.dirname(os.path.realpath(__file__))

    command =   'cd {} && '.format(job_dir) + \
                'wget ' + url + ' -O ' + gnos_id+'.xml'
    process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
    out, err = process.communicate()

    if process.returncode:
        # should not exit for just this error, improve it later
        logger.error('Unable to download metadata for: {} from {}'.format(gnos_id, url))
        return False

    if not download_dir:
        with open(os.path.join(job_dir, gnos_id+'.xml'), 'r') as f: metadata_xml_str = f.read()
        os.remove(os.path.join(job_dir, gnos_id+'.xml'))
        return metadata_xml_str
    return True

def get_gnos_key(gnos_repo):
    if 'osdc-tcga' in gnos_repo:
        gnos_key = '~/.ssh/gnos_key-tcga'
    else:
        gnos_key = '~/.ssh/gnos_key'
    
    return gnos_key


def download_datafiles(gnos_id, gnos_repo, download_dir):
    url = gnos_repo + 'cghub/data/analysis/download/' + gnos_id
    # datafiles_dir = download_dir + workflow_type
    gnos_key = get_gnos_key(gnos_repo)
    for i in range(5):
        command =   'cd {} && '.format(download_dir) + \
                    'gtdownload -c ' + gnos_key + ' -k 4 ' + url

        process = subprocess.Popen(
                command,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )

        out, err = process.communicate()

        if not process.returncode:
            #os.remove(os.path.join(download_dir, gnos_id+'.gto'))
            return True
        time.sleep(randint(1,10))  # pause a few seconds before retry
    logger.error('Unable to download datafiles for: {} from {}'.format(gnos_id, url))
    return False


def generate_uuid():
    uuid_str = str(uuid.uuid4())
    return uuid_str


def generate_md5(fname):
    hash = hashlib.md5()
    with open(fname) as f:
        for chunk in iter(lambda: f.read(4096), ""):
            hash.update(chunk)
    return hash.hexdigest()


def get_gnos_analysis_object(f):
    with open (f, 'r') as x: xml_str = x.read()
    if xmltodict.parse(xml_str).get('ResultSet') and xmltodict.parse(xml_str).get('ResultSet').get('Result') and xmltodict.parse(xml_str).get('ResultSet').get('Result').get('analysis_xml'):
        analysis_xml = xmltodict.parse(xml_str).get('ResultSet').get('Result').get('analysis_xml')
    else:
        logger.error('Could not parse the file: {}'.format(f))
        return 
    return analysis_xml



def validate_work_dir(work_dir, donors_to_be_fixed, fixed_file_dir, caller):
    for donor in donors_to_be_fixed:
        # caller = 'broad'
        gnos_entry_dir = os.path.join(work_dir, 'downloads', donor.get(caller + '_gnos_id'))
        if not os.path.isdir(gnos_entry_dir):
            logger.error('Expected GNOS entry does not exist: {}. Please ensure all GNOS entries are downloaded.'.format(gnos_entry_dir))
            sys.exit('Validating working directory failed, please check log for details.')
        if not os.path.isdir(fixed_file_dir):
            logger.error('Expected folder for {} fixed files does not exist: {}'.format(caller, fixed_file_dir))
            sys.exit('Validating working directory failed, please check log for details.')
        else:
            if caller == 'broad':
                flist = ['broad-mutect-v3.*.gz', 'broad-mutect-v3.*.gz.idx']
            elif caller == 'oxog':
                flist = ['broad-mutect-v3.*.oxoG.v2.vcf.gz', 'broad-mutect-v3.*.oxoG.v2.vcf.gz.tbi', 'dkfz-snvCalling*.oxoG.v2.vcf.gz', 'dkfz-snvCalling*.oxoG.v2.vcf.gz.tbi', \
                         'MUSE*.oxoG.v2.vcf.gz', 'MUSE*.oxoG.v2.vcf.gz.tbi', 'svcp*.oxoG.v2.vcf.gz', 'svcp*.oxoG.v2.vcf.gz.tbi']
            else:
                sys.exit('Unrecognized workflow_name'.format(caller))
            aliquot_ids = donor.get('tumor_aliquot_ids').split('|')
            for aliquot_id in aliquot_ids:
                # if not os.path.exists(os.path.join(fixed_file_dir, aliquot_id+'.broad-mutect-v3.20160222.somatic.snv_mnv.vcf.gz')) or not \
                #        os.path.exists(os.path.join(fixed_file_dir, aliquot_id+'.broad-mutect-v3.20160222.somatic.snv_mnv.vcf.gz.idx')):
                #     logger.error('No BROAD fixed files detected in: {} for donor: {}'.format(fixed_file_dir, donor.get('donor_unique_id')))
                #     sys.exit('Validating working directory failed, please check log for details.') 
                for f in flist:
                    if not len(glob.glob(os.path.join(fixed_file_dir, aliquot_id+'.'+f))) == 1:
                        logger.error('No fixed file {} detected in: {} for donor: {}'.format(aliquot_id+'.'+f, fixed_file_dir, donor.get('donor_unique_id')))
                        sys.exit('Validating working directory failed, please check log for details.')                      


def download_metadata_files(work_dir, donors_to_be_fixed, caller):
    donors_to_be_fixed_old = copy.deepcopy(donors_to_be_fixed)
    for donor in donors_to_be_fixed_old:
        # caller = 'broad'
        gnos_entry_dir = os.path.join(work_dir, 'downloads', donor.get(caller + '_gnos_id'))
        if os.path.isdir(gnos_entry_dir): 
            logger.warning('The donor: {} has downloaded files already!'.format(donor.get('donor_unique_id')))
            continue
        success = download_datafiles(donor.get(caller + '_gnos_id'), donor.get(caller + '_gnos_repo'), os.path.join(work_dir, 'downloads'))
        if not success:
            donors_to_be_fixed.remove(donor)
            continue
        success = download_metadata_xml(donor.get(caller + '_gnos_id'), donor.get(caller + '_gnos_repo'), os.path.join(work_dir, 'downloads'))
        if not success:
            donors_to_be_fixed.remove(donor)
            continue
    return donors_to_be_fixed

def metadata_fix(work_dir, donors_to_be_fixed, fixed_file_dir, allow_partial_workflow_results, caller, workflow_version):
    # donors_to_be_fixed_old = copy.deepcopy(donors_to_be_fixed)
    fixed_donors = []
    # caller = 'broad'
    for donor in donors_to_be_fixed:
        gnos_analysis_objects = {}
        gnos_entry_dir = os.path.join(work_dir, 'downloads', donor.get(caller + '_gnos_id'))
        files = get_files(donor, fixed_file_dir, gnos_entry_dir, allow_partial_workflow_results, caller)
        if not files:
            # donors_to_be_fixed.remove(donor)
            continue

        xml_file = os.path.join(gnos_entry_dir, donor.get(caller + '_gnos_id') + '.xml')
        gnos_analysis_object = get_gnos_analysis_object(xml_file)
        if not gnos_analysis_object:
            # donors_to_be_fixed.remove(donor)
            continue  

        upload_gnos_uuid = generate_uuid()
        upload_dir = os.path.join(work_dir, 'uploads', get_formal_repo_name(donor.get(caller + '_gnos_repo')), upload_gnos_uuid)
        if os.path.isdir(upload_dir): # this should never happen, but if happen regenerate a new UUID
            upload_gnos_uuid = generate_uuid()
            upload_dir = os.path.join(work_dir, 'uploads', get_formal_repo_name(donor.get(caller + '_gnos_repo')), upload_gnos_uuid)
        os.makedirs(upload_dir)

        donor.update({caller+'-'+workflow_version+'_gnos_id': upload_gnos_uuid})

        copy_file(upload_dir, files)
        apply_data_block_patches(gnos_analysis_object, files)

        create_fixed_gnos_submission(upload_dir, gnos_analysis_object, caller, workflow_version) 

        if caller == 'broad':
            updated_metadata = generate_updated_metadata(donor, work_dir)
            if not updated_metadata:
                continue

        fixed_donors.append(donor)

    return fixed_donors

def generate_updated_metadata(donor, work_dir):
    for caller in ['muse', 'broad_tar']:
        gnos_id = donor.get(caller+'_gnos_id')
        gnos_repo = donor.get(caller+'_gnos_repo')
        update_field = 'related_file_subset_uuids'
        donor_unique_id = donor.get('donor_unique_id')
        broad_gnos_id = donor.get('broad_v2_gnos_id') if donor.get('broad_v2_gnos_id') else donor.get('broad_gnos_id')
        broad_v3_gnos_id = donor.get('broad_v3_gnos_id')
        if caller == 'muse':
            uuid_set = [broad_v3_gnos_id, donor.get('broad_tar_gnos_id')]
        else:
            uuid_set = [broad_v3_gnos_id, donor.get('muse_gnos_id')]

        xml_str = download_metadata_xml(gnos_id, gnos_repo)
        gnos_analysis = xmltodict.parse(xml_str).get('ResultSet').get('Result')
        if not gnos_analysis.get('analysis_id') == gnos_id:
            logger.warning('donor: {} has wrong gnos entry {} at repo: {}'.format(donor_unique_id, gnos_id, gnos_repo))
            return False
        analysis_xml = xmltodict.parse(xml_str).get('ResultSet').get('Result').get('analysis_xml')
        for a in analysis_xml['ANALYSIS_SET']['ANALYSIS']['ANALYSIS_ATTRIBUTES']['ANALYSIS_ATTRIBUTE']:
            if a['TAG'] == update_field:
                # uuid_set = set(a['VALUE'].split(','))
                # uuid_set.remove(broad_gnos_id)
                # uuid_set.add(broad_v3_gnos_id)
                a['VALUE'] = ','.join(uuid_set)
                logger.info('donor: {} update the {} for {}_variant_calling with gnos_id: {} at gnos_repo: {}'.format(donor_unique_id, update_field, caller, gnos_id, gnos_repo))

        # analysis_xml = {'analysis_xml': analysis_xml}
        analysis_xml_str = xmltodict.unparse(analysis_xml, pretty=True)

        # analysis_xml_str = modify_metadata(xml_str, update_field, caller, gnos_repo, gnos_id, donor.get('donor_unique_id'), donor.get('broad_v2_gnos_id'), donor.get('broad_gnos_id'))
        updated_metafiles_dir = os.path.join(work_dir, 'updated_metafiles', get_formal_repo_name(gnos_repo), gnos_id)
        if not os.path.exists(updated_metafiles_dir):
            os.makedirs(updated_metafiles_dir)
        with open(updated_metafiles_dir+'/analysis.xml', 'w') as y:
            y.write(analysis_xml_str)

    return True
    

def create_fixed_gnos_submission(upload_dir, gnos_analysis_object, caller, workflow_version):
    fixed_analysis_object = copy.deepcopy(gnos_analysis_object) 

    attributes = fixed_analysis_object.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_ATTRIBUTES').get('ANALYSIS_ATTRIBUTE')
    for attr in attributes:
        if caller == 'broad' and attr.get('TAG') == 'workflow_file_subset':
            attr['VALUE'] = caller+'-'+workflow_version
        elif caller == 'oxog' and attr.get('TAG') == 'variant_workflow_name':
            old_value = copy.deepcopy(attr.get('VALUE'))
            attr['VALUE'] = old_value+'-'+workflow_version
        else:
            pass  # all others, leave it unchanged

    new_analysis_xml_str = xmltodict.unparse(fixed_analysis_object, pretty=True)
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

        if not 'fixed_files' in realfile: # this is one of the old files, use previous values
            filetype = old_files.get(filename).get('filetype')
            checksum_method = old_files.get(filename).get('checksum_method')
            checksum = old_files.get(filename).get('checksum')
        else:  # fixed file
            if filename.endswith('vcf.gz'):
                filetype = 'vcf'
            elif filename.endswith('gz.idx') or filename.endswith('gz.tbi'):
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


def copy_file(target, source):
    for s in source:
        shutil.copy(s, target)       


def get_files(donor, fixed_file_dir, gnos_entry_dir, allow_partial_workflow_results, caller):

    if caller == 'broad':
        file_name_patterns = set([
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.germline\.indel\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.germline\.indel\.vcf\.gz\.idx$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.indel\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.indel\.vcf\.gz\.idx$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.broad-dRanger[^_].+\.somatic\.sv\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.broad-dRanger[^_].+\.somatic\.sv\.vcf\.gz\.idx$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.broad-snowman.+\.somatic\.sv\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.broad-snowman.+\.somatic\.sv\.vcf\.gz\.idx$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.broad-dRanger_snowman.+\.somatic\.sv\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.broad-dRanger_snowman.+\.somatic\.sv\.vcf\.gz\.idx$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.germline\.sv\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.germline\.sv\.vcf\.gz\.idx$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.snv_mnv\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?).+\.somatic\.snv_mnv\.vcf\.gz\.idx$'
            ])
    elif caller == 'oxog':
        file_name_patterns = set([
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_annotated_broad_indel\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_annotated_broad_indel\.vcf\.gz\.tbi$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_annotated_broad_SNV\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_annotated_broad_SNV\.vcf\.gz\.tbi$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_annotated_dkfz_embl_indel\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_annotated_dkfz_embl_indel\.vcf\.gz\.tbi$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_annotated_dkfz_embl_SNV\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_annotated_dkfz_embl_SNV\.vcf\.gz\.tbi$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_annotated_sanger_indel\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_annotated_sanger_indel\.vcf\.gz\.tbi$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_annotated_sanger_SNV\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_annotated_sanger_SNV\.vcf\.gz\.tbi$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_annotated_muse_SNV\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_annotated_muse_SNV\.vcf\.gz\.tbi$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_sanger_somatic.+\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_sanger_somatic.+\.vcf\.gz\.tbi$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_dkfz_embl_somatic.+\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_dkfz_embl_somatic.+\.vcf\.gz\.tbi$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_broad_somatic.+\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\_broad_somatic.+\.vcf\.gz\.tbi$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.broad-mutect-v3.+\.somatic\.snv_mnv.+\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.broad-mutect-v3.+\.somatic\.snv_mnv.+\.vcf\.gz\.tbi$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.dkfz-snvCalling.+\.somatic\.snv_mnv.+\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.dkfz-snvCalling.+\.somatic\.snv_mnv.+\.vcf\.gz\.tbi$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.svcp.+\.somatic\.snv_mnv.+\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.svcp.+\.somatic\.snv_mnv.+\.vcf\.gz\.tbi$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.MUSE.+\.somatic\.snv_mnv.+\.vcf\.gz$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.MUSE.+\.somatic\.snv_mnv.+\.vcf\.gz\.tbi$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.gnos_files.*\.tar$',
                r'^([a-f\d]{8}(-[a-f\d]{4}){3}-[a-f\d]{12}?)\.call_stats.*\.txt\.gz\.tar$'
            ])     
    else:
        sys.exit('Unrecognized workflow_name'.format(caller))       

    matched_files = []
    tumor_aliquot_ids = donor.get('tumor_aliquot_ids').split('|')
    for aliquot in tumor_aliquot_ids:
        for file_dir in (fixed_file_dir, gnos_entry_dir): # match fixed_file dir first
            for f in glob.glob(os.path.join(file_dir, aliquot+'*')):
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
                logger.error('Missing expected variant call result file with pattern: {} for aliquot {}'.format(fp, aliquot))
            if not allow_partial_workflow_results:
                return
                # sys.exit('Missing expected variant call result file, see log file for details.')
     
    return matched_files

def generate_index_files(work_dir, donors_to_be_fixed):
    pass

def get_fix_donor_list (fixed_file_dir, vcf_info_file, donor_ids_to_be_included, donor_ids_to_be_excluded, donor_list_file, workflow_version):
    donors_to_be_fixed = []
    aliquot_ids = set()
    for f in glob.glob(os.path.join(fixed_file_dir, "*.gz")):
        aliquot_ids.add(os.path.basename(f).split('.')[0])

    with open(vcf_info_file) as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        for row in reader:
            if row.get('vcf_workflow_result_version') and row.get('vcf_workflow_result_version') == workflow_version: continue
            if donor_ids_to_be_included and not row.get('donor_unique_id') in donor_ids_to_be_included:
                logger.warning('The donor: {} is not in the donor_ids_to_be_included, skip it!'.format(row.get('donor_unique_id')))
                continue
            if donor_ids_to_be_excluded and row.get('donor_unique_id') in donor_ids_to_be_excluded:
                logger.warning('The donor: {} is in the donor_ids_to_be_excluded, skip it!'.format(row.get('donor_unique_id')))
                continue 
            tumor_aliquot_ids = set(row.get('tumor_aliquot_ids').split('|'))
            miss = tumor_aliquot_ids.difference(aliquot_ids)
            if len(miss) == 0:
                donors_to_be_fixed.append(copy.deepcopy(row))
            elif len(miss) < len(tumor_aliquot_ids): 
                logger.warning('The donor: {} is likely a multi-tumors donor and missing fixed files for tumor aliquots: {}'.format(row.get('donor_unique_id'), '|'.join(list(miss))))
                continue
            else:
                continue       
    
    write_file(donors_to_be_fixed, donor_list_file)        

    return donors_to_be_fixed


def write_file(donors_list, filename):
    with open(filename, 'w') as fh:
        header = True  
        for r in donors_list:
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
                  


def detect_folder(work_dir, folder_name):
    # dectect whether uploads dir exists, stop if exists
    job_dir = os.path.join(work_dir, folder_name)
    if os.path.isdir(job_dir):
        try:
            os.rmdir(job_dir)
        except OSError as ex:
            sys.exit('\nStop: none empty job directory exists: {}. Please confirm it\'s safe to remove, then manually remove it and try this script again.\n'.format(job_dir))

    os.mkdir(job_dir)


def generate_id_list(id_lists):
    ids_list = set()
    if id_lists:
        files = glob.glob(id_lists)
        for fname in files:
            with open(fname) as f:
                reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
                for row in reader:
                    ids_list.add(row.get('donor_unique_id'))

    return ids_list


def main(argv=None):
    parser = ArgumentParser(description="Broad snv file fix and upload",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--work_dir", dest="work_dir",
             help="Directory name containing fixed variant call files", required=True)
    parser.add_argument("-f", "--vcf_info_file", dest="vcf_info_file",
             help="vcf information file", required=False)
    parser.add_argument("-x", "--exclude_donor_id_lists", dest="exclude_donor_id_lists", 
             help="File(s) containing DONOR IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-i", "--include_donor_id_lists", dest="include_donor_id_lists", 
             help="File(s) containing DONOR IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-s", "--allow partial workflow results", dest="allow_partial_workflow_results", default=False, type=bool,
             help="Specify whether the workflow results are complete or not", required=False) 
    parser.add_argument("-w", "--workflow name", dest="workflow_name", default='broad', type=str,
             help="Specify whether the workflow results are complete or not", required=False) 
    parser.add_argument("-v", "--workflow version", dest="workflow_version", default='v3', type=str,
             help="Specify whether the workflow results are complete or not", required=False) 

    args = parser.parse_args()
    work_dir = args.work_dir
    vcf_info_file = args.vcf_info_file 
    exclude_donor_id_lists = args.exclude_donor_id_lists
    include_donor_id_lists = args.include_donor_id_lists
    allow_partial_workflow_results = args.allow_partial_workflow_results
    workflow_name = args.workflow_name
    workflow_version = args.workflow_version

    # if len(sys.argv) == 1: sys.exit('\nMust specify working directory where the variant call fixes are kept.\nPlease refer to the SOP for details how to structure the working directory.\n')
    # work_dir = sys.argv[1]
    if not os.path.isdir(work_dir): os.makedirs(work_dir+'/downloads')
    # work_dir = os.path.abspath(work_dir)

    if 'test' in work_dir:
        print('\nUsing \'test\' folder as working directory ...')

    donor_list_file = work_dir+'_donor_list.txt'
    fixed_file_dir = workflow_name+'-'+workflow_version+'_fixed_files'

    if not os.path.isdir(fixed_file_dir): sys.exit('Fixed files are missing!')

    if not vcf_info_file:
        if workflow_name == 'broad': 
            vcf_info_file = 'broad_successful_uploads.txt'  
        elif workflow_name == 'oxog':
            vcf_info_file = 'variant_call_entries_with_broad_oxog_filter_applied.txt'
        else:
            sys.exit('Unrecognized workflow!')
    if not os.path.exists(vcf_info_file): sys.exit('Helper file is missing')

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
    donor_ids_to_be_excluded = generate_id_list(exclude_donor_id_lists)
    donor_ids_to_be_included = generate_id_list(include_donor_id_lists)

    # generate the donors_to_be_fixed list from the files in fixed_files folder
    donors_to_be_fixed = get_fix_donor_list(fixed_file_dir, vcf_info_file, donor_ids_to_be_included, donor_ids_to_be_excluded, donor_list_file, workflow_version)

    # now download data files and metadata xml
    print('\nDownloading data files and metadata XML from GNOS ...')
    donors_to_be_fixed = download_metadata_files(work_dir, donors_to_be_fixed, workflow_name)

    print('\nThe number of donors to be fixed is: {}'.format(len(donors_to_be_fixed)))
    # validate working direcotry first
    print('\nValidating working directory...')
    validate_work_dir(work_dir, donors_to_be_fixed, fixed_file_dir, workflow_name)

    # dectect whether uploads dir exists, stop if exists
    print('\nDetecting output folders...')
    detect_folder(work_dir, 'uploads')
    if workflow_name == 'broad': detect_folder(work_dir, 'updated_metafiles')

    # now process metadata xml fix and merge
    print('\nPreparing new GNOS submissions and updated related metadata XML files if needed...')
    fixed_donors = metadata_fix(work_dir, donors_to_be_fixed, fixed_file_dir, allow_partial_workflow_results, workflow_name, workflow_version)

    # write the fixed donor informaton 
    write_file(fixed_donors, 'fixed_'+donor_list_file)
    if not os.path.exists('fixed_'+donor_list_file): sys.exit('Fixed donor list file is missing!')
    
    print('Submission folder located at: {}'.format(os.path.join(work_dir, 'uploads')))
    print('Processing log file: {}'.format(log_file))
    print('Done!\n')


if __name__ == "__main__":
    sys.exit(main())