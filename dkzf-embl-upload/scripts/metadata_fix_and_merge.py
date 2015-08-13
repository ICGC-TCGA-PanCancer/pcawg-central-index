#!/usr/bin/env python


import sys
import os
import re
import subprocess
import shutil
import xmltodict
import dicttoxml
import yaml
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

logger = logging.getLogger('dkfz/embl merge and upload')
# create console handler with a higher log level
ch = logging.StreamHandler()

def download_metadata_xml(gnos_id, gnos_repo, workflow_type, download_dir):
    logger.info('download metadata xml for {} from GNOS repo: {}'.format(gnos_id, gnos_repo))

    metadata_xml_dir = download_dir + workflow_type + '/' + gnos_id
    if not os.path.exists(metadata_xml_dir):
        logger.warning('data for analysis object: {} have not been download yet'.format(gnos_id))
        return
    if workflow_type == 'dkfz':
        fixed_dir = metadata_xml_dir + '/fixed_files'
        if not os.path.exists(fixed_dir) or not os.listdir(fixed_dir):
            logger.warning('data for analysis object: {} have not been fixed yet'.format(gnos_id))
            return

    logger.info('download metadata xml from GNOS repo: {} for analysis object: {}'.format(gnos_repo, gnos_id))
    
    url = gnos_repo + 'cghub/metadata/analysisFull/' + gnos_id
    response = None
    try:
        response = requests.get(url, stream=True, timeout=5)
    except: # download failed, no need to do anything
        pass

    if not response or not response.ok:
        logger.warning('unable to download metadata for: {} from {}'.format(gnos_id, url))
        return
    else:
        metadata_xml_str = response.text
        gnos_ao = xmltodict.parse(metadata_xml_str).get('ResultSet').get('Result')
        ao_uuid = gnos_ao.get('analysis_id')

        metadata_xml_file = metadata_xml_dir + '/' + ao_uuid  + '.xml'
        with open(metadata_xml_file, 'w') as f:  # write to metadata xml file now
            f.write(metadata_xml_str.encode('utf8'))
    
        return metadata_xml_file


def generate_uuid():
    uuid_str = str(uuid.uuid4())
    return uuid_str


def process(conf, dkfz_embl_results_info):
    # download_dir
    download_dir = conf.get('download_dir')

    # upload_dir
    upload_dir = conf.get('upload_dir')

    # manifest_file
    manifest_file = conf.get('manifest_file')

    # read the manifest info
    manifest = {}
    read_manifest(manifest, manifest_file)

    with open(manifest_file, 'aw') as m:
        with open(dkfz_embl_results_info, 'r') as f:
            for line in f:
                if line.startswith('donor_unique_id'): continue
                if len(line.strip()) == 0: continue
                donor_unique_id, submitter_donor_id, dcc_project_code, embl_gnos_id, embl_gnos_repo, \
                    dkfz_gnos_id, dkfz_gnos_repo = str.split(line.rstrip(),'\t')
                if manifest.get('donor_unique_id') and donor_unique_id in manifest.get('donor_unique_id'): continue
                embl_xml_file = download_metadata_xml(embl_gnos_id, embl_gnos_repo, 'embl', download_dir)
                dkfz_xml_file = download_metadata_xml(dkfz_gnos_id, dkfz_gnos_repo, 'dkfz', download_dir)
                if embl_xml_file and os.path.isfile(embl_xml_file) and dkfz_xml_file and os.path.isfile(dkfz_xml_file):
                    merged_analysis_xml = merge_metadata_xml(embl_xml_file, dkfz_xml_file, upload_dir, download_dir, embl_gnos_id, dkfz_gnos_id)
                    if merged_analysis_xml:
                        new_gnos_id = generate_uuid()
                        upload_xml_file = upload_dir + embl_gnos_id + '.' + dkfz_gnos_id + '/' + new_gnos_id + '/analysis.xml'
                        # Create the symlinks for all the files
                        fixed_files = merged_analysis_xml.get('ANALYSIS_SET').get('ANALYSIS').get('DATA_BLOCK').get('FILES').get('FILE')
                        dst_dir = upload_dir + embl_gnos_id + '.' + dkfz_gnos_id + '/' + new_gnos_id
                        create_file_symlinks(fixed_files, download_dir, dst_dir, embl_gnos_id, dkfz_gnos_id)
                        # Write the merged_analysis_xml
                        merged_xml = xmltodict.unparse(merged_analysis_xml, pretty=True)
                        write_to_xml(upload_xml_file, merged_xml)
                        # Write to the manifest_file with new_gnos_id and show the donors which have complete the dkfz/embl merge
                        new_line = '\t'.join([donor_unique_id, submitter_donor_id, dcc_project_code, embl_gnos_id, embl_gnos_repo, dkfz_gnos_id, dkfz_gnos_repo, new_gnos_id])
                        m.write(new_line + '\n')


def merge_metadata_xml(embl_xml_file, dkfz_xml_file, upload_dir, download_dir, embl_gnos_id, dkfz_gnos_id):
    embl_analysis_xml = get_analysis_xml(embl_xml_file)
    dkfz_analysis_xml = get_analysis_xml(dkfz_xml_file)
    merged_analysis_xml = copy.deepcopy(embl_analysis_xml)

    embl_analysis = embl_analysis_xml.get('ANALYSIS_SET').get('ANALYSIS')
    dkfz_analysis = dkfz_analysis_xml.get('ANALYSIS_SET').get('ANALYSIS')
    merged_analysis = merged_analysis_xml.get('ANALYSIS_SET').get('ANALYSIS')

    # Description
    merged_analysis['DESCRIPTION'] = 'New merged dkfz/embl xml'
    # Pipe_line
    dkfz_analysis.get('ANALYSIS_TYPE').get('REFERENCE_ALIGNMENT').get('PROCESSING').get('PIPELINE').get('PIPE_SECTION')['STEP_INDEX'] = 2
    merged_analysis.get('ANALYSIS_TYPE').get('REFERENCE_ALIGNMENT').get('PROCESSING')['PIPELINE'] = \
        [embl_analysis.get('ANALYSIS_TYPE').get('REFERENCE_ALIGNMENT').get('PROCESSING').get('PIPELINE'), \
            dkfz_analysis.get('ANALYSIS_TYPE').get('REFERENCE_ALIGNMENT').get('PROCESSING').get('PIPELINE')]
    
    # Data_block
    embl_data_block_files = embl_analysis.get('DATA_BLOCK').get('FILES').get('FILE')
    embl_data_block_files = update_fixed_data_block_files(embl_data_block_files, 'embl', upload_dir, download_dir, embl_gnos_id)

    dkfz_data_block_files = dkfz_analysis.get('DATA_BLOCK').get('FILES').get('FILE')
    dkfz_data_block_files = update_fixed_data_block_files(dkfz_data_block_files, 'dkfz', upload_dir, download_dir, dkfz_gnos_id)


    if not embl_data_block_files or not dkfz_data_block_files: return
    merged_analysis.get('DATA_BLOCK').get('FILES')['FILE'] = embl_data_block_files + dkfz_data_block_files

    # Analysis_attributes
    embl_analysis_attrib = get_analysis_attrib(embl_analysis)
    dkfz_analysis_attrib = get_analysis_attrib(dkfz_analysis)
    embl_attrib_set = set(embl_analysis_attrib.items())
    dkfz_attrib_set = set(dkfz_analysis_attrib.items())
    merged_attrib_set = embl_attrib_set | dkfz_attrib_set
    analysis_attrib = generate_analysis_attrib_from_set(merged_attrib_set, embl_attrib_set, dkfz_attrib_set)
    merged_analysis.get('ANALYSIS_ATTRIBUTES')['ANALYSIS_ATTRIBUTE'] = generate_analysis_attrib_list_from_dict(analysis_attrib)

    return merged_analysis_xml


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


def create_file_symlinks(file_list, src_dir, dst_dir, embl_gnos_id, dkfz_gnos_id):
    for f in file_list:
        fname = str(f.get('@filename'))
        if 'dkfz' in fname:
            f_fixed_src = src_dir + '/dkfz/' + dkfz_gnos_id + '/fixed_files/' + f.get('@filename')
            f_src = src_dir + '/dkfz/' + dkfz_gnos_id + '/' + f.get('@filename')
        else:
            f_fixed_src = src_dir + '/embl/' + embl_gnos_id + '/fixed_files/' + f.get('@filename')
            f_src = src_dir + '/embl/' + embl_gnos_id + '/' + f.get('@filename')           
        f_dst = dst_dir + '/' + f.get('@filename')
        #print f_dst
        if not os.path.isdir(os.path.dirname(f_dst)):
            os.makedirs(os.path.dirname(f_dst))
        if os.path.isfile(f_fixed_src):
            os.symlink(f_fixed_src, f_dst)
        elif os.path.isfile(f_src):
            os.symlink(f_src, f_dst)
        else: # should not happen
            logger.warning('file: {} is missing in the downloads folder'.format(f.get('@filename')))
            return 0


def generate_md5(file):
    with open (file, 'r') as x: data = x.read()
    md5 = hashlib.md5(data).hexdigest()

    return md5



def write_to_xml(fname, xml):      
    if not os.path.exists(os.path.dirname(fname)):
        os.makedirs(os.path.dirname(fname))
    with open(fname, "w") as f:
        f.write(xml)


def generate_analysis_attrib_list_from_dict(analysis_attrib):
    analysis_attrib_list = []
    for k,v in analysis_attrib.iteritems():
        element = OrderedDict()
        element['TAG'] = k
        element['VALUE'] = v
        analysis_attrib_list.append(element)
    return analysis_attrib_list


def generate_analysis_attrib_from_set(merged_attrib_set, embl_attrib_set, dkfz_attrib_set):
    analysis_attrib = {}
    for a in merged_attrib_set:
        if not analysis_attrib.get(a[0]):
            analysis_attrib[a[0]] = a[1]  

        elif a[0] == 'variant_pipeline_output_info':
            current_output = json.loads(a[1]).get('workflow_outputs') if a[1] is not None else []
            current_files = {}
            for c_out in current_output:
                if not c_out.get('files'): continue
                current_files = c_out.get('files')

            merged_output = json.loads(analysis_attrib.get(a[0])).get('workflow_outputs') if analysis_attrib.get(a[0]) else []
            for m_out in merged_output:
                if not m_out.get('files'): continue
                merged_files = m_out.get('files')
            merged_files.update(current_files)
            analysis_attrib[a[0]] = json.dumps(merged_output)

        elif a[0] == 'variant_workflow_name':
            analysis_attrib[a[0]] = 'DKFZ_EMBL_Merged'

        else:
            if a in embl_attrib_set:
                analysis_attrib[a[0]+'_embl'] = a[1]
                analysis_attrib[a[0]+'_dkfz'] = analysis_attrib.get(a[0])
            elif a in dkfz_attrib_set:
                analysis_attrib[a[0]+'_embl'] = analysis_attrib.get(a[0]) 
                analysis_attrib[a[0]+'_dkfz'] = a[1]       
            else:
                logger.warning('unknown analysis attribute: {}'.format(a))
            del analysis_attrib[a[0]]    
    return analysis_attrib           

 
def get_analysis_attrib(analysis):
    analysis_attrib = {}
    for a in analysis['ANALYSIS_ATTRIBUTES']['ANALYSIS_ATTRIBUTE']:
        if not analysis_attrib.get(a['TAG']):
            analysis_attrib[a['TAG']] = a['VALUE']
        else:
            logger.warning('duplicated analysis attribute key: {}'.format(a['TAG']))
    return analysis_attrib


def get_analysis_xml(f):
    with open (f, 'r') as x: xml_str = x.read()
    analysis_xml = xmltodict.parse(xml_str).get('ResultSet').get('Result').get('analysis_xml')
    return analysis_xml


def read_manifest(manifest, manifest_file):
    if not os.path.isfile(manifest_file):
        with open(manifest_file, 'w') as m:
            header = ['donor_unique_id', 'submitter_donor_id', 'dcc_project_code', 'embl_gnos_id', 'embl_gnos_repo', 'dkfz_gnos_id', 'dkfz_gnos_repo']
            m.write('\t'.join(header) + '\n')

    else:
        with open(manifest_file, 'r') as r:
            manifest['donor_unique_id'] = set()
            for line in r:
                if line.startswith('donor_unique_id'): continue
                if len(line.strip()) == 0: continue
                donor_unique_id = str.split(line.rstrip(),'\t')[0]
                manifest['donor_unique_id'].add(donor_unique_id)

    return manifest


def main(argv=None):
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    current_time = time.strftime("%Y-%m-%d_%H-%M-%S_%Z")

    parser = ArgumentParser(description="PCAWG DKFZ/EMBL Variant Calling Results Merge",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-c", "--config", dest="config",
             help="Configuration file for download/upload directories", required=True)
    parser.add_argument("-l", "--dkfz_embl_results_info", dest="dkfz_embl_results_info",
             help="File list the dkfz/embl results for merge", required=True)

    args = parser.parse_args()
    dkfz_embl_results_info = args.dkfz_embl_results_info
    conf_file = args.config

    with open(conf_file) as f:
        conf = yaml.safe_load(f)


    # test_dir
    test_dir = conf.get('test_dir')

    # initialize output directory
    # mani_output_dir = initialize_output_dir(output_dir, current_time)

    logger.setLevel(logging.INFO)
    ch.setLevel(logging.WARN)

    log_file = test_dir + '/' + current_time + '.uploader.log'

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)


    process(conf, dkfz_embl_results_info)

    return 0


if __name__ == "__main__":
    sys.exit(main())