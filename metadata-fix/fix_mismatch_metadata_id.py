#!/usr/bin/env python

import re
import os
import xmltodict
import json
import hashlib
import requests
import sys
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import shutil
from operator import itemgetter
import csv
from collections import OrderedDict
import logging
import time
from random import randint
import subprocess


# logger = logging.getLogger('fix metadata')
# ch = logging.StreamHandler()

# def download_metadata_xml(gnos_repo, ao_uuid):
#     url = gnos_repo + '/cghub/metadata/analysisFull/' + ao_uuid

#     response = None

#     try:
#         response = requests.get(url, stream=True, timeout=30)
#     except: # download failed, no need to do anything
#         pass

#     if not response or not response.ok:
#         print('Unable to download {}'.format(url))
#         return None
#     else:
#         metadata_xml_str = response.text
#         return metadata_xml_str

def download_metadata_xml(gnos_repo, gnos_id, download_dir=None):
    
    url = gnos_repo + 'cghub/metadata/analysisFull/' + gnos_id

    command = 'wget ' + url + ' -O ' + gnos_id+'.xml'
    process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
    out, err = process.communicate()

    if process.returncode:
        # should not exit for just this error, improve it later
        return None
    
    with open(gnos_id+'.xml', 'r') as f: metadata_xml_str = f.read()
    os.remove(gnos_id+'.xml')
    return metadata_xml_str


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

def fix_illegal_id(xml_str, id_mapping, fix_pattern, id_types):
    for id_type in id_types:
        if not id_mapping.get(id_type): continue
        for key, value in id_mapping.get(id_type).items():
            if fix_pattern == 'whole_string':
                key = key.replace('*', '\*')
                xml_str = re.sub('>'+key+'<', '>'+value+'<', xml_str)
                xml_str = re.sub(' '+key+' ', ' '+value+' ', xml_str)
                xml_str = re.sub('"'+key+'"', '"'+value+'"', xml_str)
                xml_str = re.sub(' '+key+'<', ' '+value+'<', xml_str)
                xml_str = re.sub(' '+key+'.', ' '+value+'.', xml_str)


            elif fix_pattern == 'key_value':
                xml_dict = xmltodict.parse(xml_str)
                if xml_dict.get('ResultSet') and \
                   xml_dict.get('ResultSet').get('Result') and \
                   xml_dict.get('ResultSet').get('Result').get('analysis_xml'):
                    analysis_xml = xml_dict.get('ResultSet').get('Result').get('analysis_xml')
                    for a in analysis_xml['ANALYSIS_SET']['ANALYSIS']['ANALYSIS_ATTRIBUTES']['ANALYSIS_ATTRIBUTE']:
                        if a.get('TAG') == id_type and a.get('VALUE') == key:
                            a['VALUE'] = value
                        elif a.get('TAG') in ['variant_pipeline_input_info', 'variant_pipeline_output_info']:
                            a['VALUE'] = re.sub('"'+id_type+'":"'+key+'",', '"'+id_type+'":"'+value+'",', a.get('VALUE'))
                        else:
                            pass
                    xml_str = xmltodict.unparse(xml_dict, pretty=True)
 
                else:
                    print('Could not parse the analysis xml!')

            else:
                print('Norecognized fix_pattern!'.format(fix_pattern))

    return xml_str

def generate_metadata(xml_str, gnos_id, gnos_repo, fixed_dir, subtype):

    xml_dir = os.path.join(fixed_dir, subtype, gnos_repo, gnos_id)
    if os.path.exists(xml_dir): shutil.rmtree(xml_dir, ignore_errors=True)  # empty the folder if exists
    os.makedirs(xml_dir) 
    
    for xml_subtype in ['analysis', 'experiment', 'run']:
        if xmltodict.parse(xml_str).get('ResultSet').get('Result').get(xml_subtype+'_xml'):
            xml_subtype_str = xmltodict.unparse(
                                                xmltodict.parse(xml_str).
                                                get('ResultSet').
                                                get('Result').
                                                get(xml_subtype+'_xml'), pretty=True)
            with open(xml_dir+'/'+ xml_subtype+ '.xml', 'w') as y:
                y.write(xml_subtype_str)


def read_annotations(annotations, type, file_name, subtype):
    # file_name=os.path.join('annotation_files', file_name)
    if not os.path.isfile(file_name):
        print('The id_fixes_file does not exist!')
        return
    with open(file_name, 'r') as r:
        if type == 'ESAD-UK':
            annotations[type] = {}
            reader = csv.DictReader(r, delimiter='\t')
            for row in reader:
                if row.get('status') == 'BROKEN': continue
                url = row.get('gnos_metadata_url')
                gnos_repo = get_formal_repo_name(url.split('cghub')[0])
                gnos_id = url.split('/')[-1]
                annotations[type][gnos_id] = gnos_repo
        elif type == 'mismatch_metadata':
            # annotations[type] = {}
            reader = csv.DictReader(r, delimiter='\t')
            for row in reader:
                annotations[type][row.get('gnos_id')] = {}
                if row.get('dcc_project_code') == 'ESAD-UK' and row.get('workflow') == 'wgs_bwa_alignment':
                    if not annotations['ESAD-UK'].get(row.get('gnos_id')): continue
                    annotations[type][row.get('gnos_id')]['gnos_repo_with_good_copy'] = annotations['ESAD-UK'].get(row.get('gnos_id'))                 
                elif row.get('dcc_project_code') == 'LIRI-JP' and row.get('workflow') == 'wgs_bwa_alignment':
                    annotations[type][row.get('gnos_id')]['gnos_repo_with_good_copy'] = 'riken'
                else: # all other mismatch has no information of good copy for now
                    continue

                annotations[type][row.get('gnos_id')]['gnos_repo_all'] = [get_formal_repo_name(r) for r in row.get('gnos_repo').split('|')]
        elif type == 'id_mapping':
            if not annotations.get(type): annotations[type] = {}
            reader = csv.DictReader(r, delimiter='\t')
            for row in reader:
                donor_unique_id = row.get('donor_unique_id') if row.get('donor_unique_id') else row.get('project_code')+'::'+row.get('submitter_donor_id')
                if not annotations[type].get(donor_unique_id): annotations[type][donor_unique_id] = {}
                if not annotations[type].get(donor_unique_id).get(subtype): annotations[type][donor_unique_id][subtype] = OrderedDict()
                if not row.get('old_'+subtype) or not row.get('new_'+subtype): continue
                annotations[type][donor_unique_id][subtype][row.get('old_'+subtype)] = row.get('new_'+subtype)

        else:
            print('unknown annotation type: {}'.format(type))
    return annotations


def write_file(flist, fn):
    with open(fn, 'w') as f:
        header = True  
        for r in flist:
            if header:
                f.write('\t'.join(r.keys()) + '\n')
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
            f.write('\t'.join(line) + '\n') 


def main(argv=None):

    parser = ArgumentParser(description="Metadata fix",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--information for gnos entries", dest="fname",
             help="Specify file to update the metadata", required=True)
    parser.add_argument("-x", "--id_fixes_file", dest="id_fixes_file", 
             help="Specify file of id fixes", required=True)
    parser.add_argument("-m", "--fix_pattern", dest="fix_pattern", default="whole_string",
             help="Specify pattern of id match [whole_string, key_value]", required=False)
    parser.add_argument("-p", "--fix_type", dest="fix_type", default="id_mapping",
             help="Specify type of id fixes['id_mapping', 'mismatch_metadata']", required=False) 
    parser.add_argument("-t", "--id_types", dest="id_types", nargs="*",
             help="Specify type of id fixes['dcc_project_code', 'submitter_donor_id', 'dcc_specimen_type', 'submitter_specimen_id', 'submitter_sample_id']", required=False)                 
    parser.add_argument("-o", "--fixed_metadata_dir output folder", dest="fixed_dir", default="fixed_dir",
             help="Specify output folder for the fixed metadata", required=False)


    args = parser.parse_args()
    fname = args.fname
    id_fixes_file = args.id_fixes_file
    fix_pattern = args.fix_pattern
    fix_type = args.fix_type
    id_types = args.id_types
    fixed_dir = args.fixed_dir

    id_types = list(id_types) if id_types else ['dcc_project_code', 'submitter_donor_id', 'dcc_specimen_type', 'submitter_specimen_id', 'submitter_sample_id']

    # if not fixed_dir: fixed_dir = 'fixed_dir'
    # logger.setLevel(logging.INFO)
    # ch.setLevel(logging.WARN)

    # log_file = re.sub(r'\.py$', '.log', os.path.basename(__file__))

    # fh = logging.FileHandler(log_file)
    # fh.setLevel(logging.DEBUG)
    # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # fh.setFormatter(formatter)
    # ch.setFormatter(formatter)
    # logger.addHandler(fh)
    # logger.addHandler(ch)

    annotations = {
        'id_mapping': {},
        'mismatch_metadata':{}
    }
    for id_type in id_types:
        read_annotations(annotations, fix_type, id_fixes_file, id_type)
    # read_annotations(annotations, 'ESAD-UK', 'ESAD-UK_10broken.txt')
    # read_annotations(annotations, 'mismatch_metadata', 'specimens_with_mismatch_effective_xml_md5sum.txt')
    # read_annotations(annotations, 'id_mapping', 'ESAD-UK_id_fixes.tsv')
    # read_annotations(annotations, 'id_mapping', 'PAEN-AU_id_fixes.tsv')
    # read_annotations(annotations, 'id_mapping', 'MELA-AU_PCAWG-DCC_specimen_id_mapping.tsv', 'submitter_specimen_id')
    # read_annotations(annotations, 'id_mapping', 'OV-AU_id_fixes.tsv', 'submitter_specimen_id')
    # read_annotations(annotations, 'id_mapping', 'OV-AU_id_fixes.tsv', 'submitter_sample_id')
    # read_annotations(annotations, 'id_mapping', 'OV-AU_specimen_type_fixes.tsv', 'dcc_specimen_type')
    # read_annotations(annotations, 'id_mapping', 'PAEN-IT_project_code_fixes.tsv', 'dcc_project_code')
    # read_annotations(annotations, 'id_mapping', 'PAEN-AU_project_code_fixes.tsv', 'dcc_project_code')
    # read_annotations(annotations, 'id_mapping', 'ESAD-UK_specimen_type_fixes.tsv', 'dcc_specimen_type')
    # read_annotations(annotations, 'id_mapping', 'LIRI-JP_specimen_id_fixes.tsv', 'submitter_specimen_id')

    fixed_metadata_list = []
    with open(fname, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row.get('entity_type') and row.get('entity_type') == 'unaligned_bams': continue

            fixed_metadata = OrderedDict()
            fixed_metadata['donor_unique_id'] = row.get('donor_unique_id') if row.get('donor_unique_id') else row.get('dcc_project_code')+'::'+row.get('submitter_donor_id')
            fixed_metadata['dcc_project_code'] = row.get('dcc_project_code') if row.get('dcc_project_code') else None
            fixed_metadata['aliquot_id'] = row.get('aliquot_id') if row.get('aliquot_id') else None
            fixed_metadata['entity_type'] = row.get('entity_type') if row.get('entity_type') else None
            fixed_metadata['gnos_id'] = row.get('gnos_id') if row.get('gnos_id') else row.get('analysis_id')
            fixed_metadata['gnos_repo_original'] = get_formal_repo_name(row.get('gnos_repo')) if row.get('gnos_repo') else row.get('gnos_server')

            if not annotations.get('id_mapping').get(fixed_metadata.get('donor_unique_id')) and \
               not annotations.get('mismatch_metadata').get(fixed_metadata.get('gnos_id')): continue

            # has both id and metadata mismatch issues
            if annotations.get('id_mapping').get(fixed_metadata.get('donor_unique_id')) and annotations.get('mismatch_metadata').get(fixed_metadata.get('gnos_id')):
                fixed_metadata['gnos_repo_download'] = annotations.get('mismatch_metadata').get(fixed_metadata.get('gnos_id')).get('gnos_repo_with_good_copy')
                fixed_metadata['fixed_type'] = 'fixed_illegal_id_and_mismatch'

            elif annotations.get('id_mapping').get(fixed_metadata.get('donor_unique_id')) and not annotations.get('mismatch_metadata').get(fixed_metadata.get('gnos_id')):
                fixed_metadata['gnos_repo_download'] = fixed_metadata['gnos_repo_original']
                fixed_metadata['fixed_type'] = 'fixed_illegal_id'

            elif not annotations.get('id_mapping').get(fixed_metadata.get('donor_unique_id')) and annotations.get('mismatch_metadata').get(fixed_metadata.get('gnos_id')):
                fixed_metadata['gnos_repo_download'] = annotations.get('mismatch_metadata').get(fixed_metadata.get('gnos_id')).get('gnos_repo_with_good_copy')
                fixed_metadata['fixed_type'] = 'fixed_mismatch'
                if fixed_metadata['gnos_repo_download'] == fixed_metadata['gnos_repo_original']: continue

            else:
                print('Warning: {} not any of the above situations!!!'.format(fixed_metadata.get('gnos_id'))) 
                continue
            
            # download orignal xml 
            xml_str = download_metadata_xml(get_formal_repo_name(fixed_metadata['gnos_repo_original']), fixed_metadata['gnos_id'])
            print fixed_metadata['gnos_id']

            if not xmltodict.parse(xml_str).get('ResultSet') or not xmltodict.parse(xml_str).get('ResultSet').get('Result'): 
                print('Unable to download the xml of {} from {}'.format(fixed_metadata['gnos_id']), fixed_metadata['gnos_repo_original'])
                continue
            generate_metadata(xml_str, fixed_metadata['gnos_id'], fixed_metadata['gnos_repo_original'], fixed_dir, 'orignal')
            if fixed_metadata['fixed_type'] == 'fixed_illegal_id':
                # fix the illegal ids
                xml_str = fix_illegal_id(xml_str, annotations.get('id_mapping').get(fixed_metadata.get('donor_unique_id')), fix_pattern, id_types)
            elif fixed_metadata['fixed_type'] == 'fixed_illegal_id_and_mismatch':
                # download from the repo with metadata fixed copy
                xml_str = download_metadata_xml(get_formal_repo_name(fixed_metadata['gnos_repo_download']), fixed_metadata['gnos_id'])
                if not xmltodict.parse(xml_str).get('ResultSet') or not xmltodict.parse(xml_str).get('ResultSet').get('Result'): 
                    print('Unable to download the xml of {} from {}'.format(fixed_metadata['gnos_id']), fixed_metadata['gnos_repo_original'])
                    continue
                xml_str = fix_illegal_id(xml_str, annotations.get('id_mapping').get(fixed_metadata.get('donor_unique_id')), fix_pattern, id_types)
            elif fixed_metadata['fixed_type'] == 'fixed_mismatch':                            
                # download from the repo with metadata fixed copy
                xml_str = download_metadata_xml(get_formal_repo_name(fixed_metadata['gnos_repo_download']), fixed_metadata['gnos_id'])
                if not xmltodict.parse(xml_str).get('ResultSet') or not xmltodict.parse(xml_str).get('ResultSet').get('Result'): 
                    print('Unable to download the xml of {} from {}'.format(fixed_metadata['gnos_id']), fixed_metadata['gnos_repo_original'])
                    continue
            else:
                print('Warning: this should not happen!!!') 
                continue
            generate_metadata(xml_str, fixed_metadata['gnos_id'], fixed_metadata['gnos_repo_original'], fixed_dir, fname.split('_')[0])
            generate_metadata(xml_str, fixed_metadata['gnos_id'], fixed_metadata['gnos_repo_original'], fixed_dir, 'fixed_all')


            fixed_metadata_list.append(fixed_metadata)
    if fixed_metadata_list:
        write_file(fixed_metadata_list, fixed_dir+'/'+fname.split('_')[0]+'_fixed_metadata_list.txt')        

    return 0


if __name__ == "__main__":
    sys.exit(main()) 
