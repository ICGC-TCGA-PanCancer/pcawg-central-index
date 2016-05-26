#!/usr/bin/env python

import os
import re
import click
import csv
import copy
import xmltodict
import json
import sys
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from elasticsearch import Elasticsearch
import logging
import glob
from collections import OrderedDict
import gzip
import requests
import hashlib
import subprocess
import time
import calendar
import ftplib


logger = logging.getLogger('Collect sample and gnos_xml data')
ch = logging.StreamHandler()

ega_box_token = os.environ.get('EGA_TOKEN')

def generate_es_query(dcc_project_code):
    es_query = {
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
                  "terms": {
                    "dcc_project_code": [dcc_project_code]
                  }
                }
              ],
              "must_not": [
                {
                  "regexp": {
                    "dcc_project_code": ".*-US"
                  }
                }
              ]
            }
        },
        "size": 10000
    }
    return es_query


def get_donors_list(es, es_index, dcc_project_code):
    response = es.search(index=es_index, body=generate_es_query(dcc_project_code))
    
    donors_list = set()
    for p in response['hits']['hits']:
        donors_list.add(p.get('fields').get('donor_unique_id')[0])

    return donors_list 


def collect_sample(donors_list, sample_ids_to_be_included, sample_ids_to_be_excluded, dcc_project_code, ega_dir, pcawg_sample_sheet, seq, annotations):
    for sequence_type in seq:
        print('\nCollecting the Sample Data for sequence_type: {} of project: {}'.format(sequence_type, dcc_project_code))
        sample_sheet = []
        with open(pcawg_sample_sheet, 'r') as s:
            reader = csv.DictReader(s, delimiter='\t')
            for sample_info in reader:
                if not sample_info.get('dcc_project_code') == dcc_project_code: continue
                if not sample_info.get('donor_unique_id') in donors_list: continue
                sample_id = sample_info.get('icgc_sample_id')
                if sample_ids_to_be_included and not sample_id in sample_ids_to_be_included: continue
                if sample_ids_to_be_excluded and sample_id in sample_ids_to_be_excluded: continue
                library_strategy = sample_info.get('library_strategy').lower()
                if not library_strategy == sequence_type: continue
            
                sample = OrderedDict()

                try:
                    sample['subject_id'] = sample_info['icgc_donor_id']
                    if annotations.get('gender') and annotations.get('gender').get(sample_info.get('donor_unique_id')):
                        sample['gender'] = annotations.get('gender').get(sample_info.get('donor_unique_id')) 
                    else:
                        click.echo('Warning: missing gender informaion for donor: %s' % sample_info.get('donor_unique_id'), err=True)
                        continue
                    sample['phenotype'] = None
                    sample['icgc_project_code'] = sample_info['dcc_project_code']

                    for tag in ['icgc_donor_id', 'submitter_donor_id', 'icgc_specimen_id', 'submitter_specimen_id', 'icgc_sample_id', 'submitter_sample_id']:
                        if sample_info.get(tag):
                            sample[tag] = sample_info.get(tag)  
                        else: 
                            click.echo('Warning: missing {0} informaion for donor: {1}'.format(tag, sample_info.get('donor_unique_id')), err=True)
                            return
                    sample['specimen_type'] = sample_info['dcc_specimen_type']
                    sample['aliquot_id/sample_uuid'] = sample_info['aliquot_id']


                except KeyError, e:
                    click.echo('Error: KeyError, %s' % str(e), err=True)
                    return
                except IndexError, e:
                    click.echo('Error: IndexError, %s' % str(e), err=True)
                    return
                except Exception, e:
                    click.echo('Error: %s' % str(e), err=True)
                    return

                sample_sheet.append(copy.deepcopy(sample))     

        if sample_sheet:    
            out_dir = os.path.join(ega_dir, dcc_project_code, 'sample')
            if not os.path.isdir(out_dir): os.makedirs(out_dir)
            epoch_time = str(int(calendar.timegm(time.gmtime())))
            out_file = os.path.join(out_dir, 'sample.'+dcc_project_code+'.'+sequence_type+'_'+epoch_time+'.tsv')
            write_tsv_file(sample_sheet, out_file)

def effective_xml_md5sum(xml_str):

    xml_str = re.sub(r'<ResultSet .+?>', '<ResultSet>', xml_str)
    xml_str = re.sub(r'<last_modified>.+?</last_modified>', '<last_modified></last_modified>', xml_str)
    xml_str = re.sub(r'<upload_date>.+?</upload_date>', '<upload_date></upload_date>', xml_str)
    xml_str = re.sub(r'<published_date>.+?</published_date>', '<published_date></published_date>', xml_str)
    xml_str = re.sub(r'<center_name>.+?</center_name>', '<center_name></center_name>', xml_str)
    xml_str = re.sub(r'<analyte_code>.+?</analyte_code>', '<analyte_code></analyte_code>', xml_str)
    xml_str = re.sub(r'<reason>.+?</reason>', '<reason></reason>', xml_str)
    xml_str = re.sub(r'<study>.+?</study>', '<study></study>', xml_str)
    xml_str = re.sub(r'<sample_accession>.+?</sample_accession>', '<sample_accession></sample_accession>', xml_str)

    #xml_str = re.sub(r'<dcc_project_code>.+?</dcc_project_code>', '<dcc_project_code></dcc_project_code>', xml_str)
    #xml_str = re.sub(r'<participant_id>.+?</participant_id>', '<participant_id></participant_id>', xml_str)
    xml_str = re.sub(r'<dcc_specimen_type>.+?</dcc_specimen_type>', '<dcc_specimen_type></dcc_specimen_type>', xml_str)

    xml_str = re.sub(r'<specimen_id>.+?</specimen_id>', '<specimen_id></specimen_id>', xml_str)
    xml_str = re.sub(r'<sample_id>.+?</sample_id>', '<sample_id></sample_id>', xml_str)
    xml_str = re.sub(r'<use_cntl>.+?</use_cntl>', '<use_cntl></use_cntl>', xml_str)
    xml_str = re.sub(r'<library_strategy>.+?</library_strategy>', '<library_strategy></library_strategy>', xml_str)
    xml_str = re.sub(r'<platform>.+?</platform>', '<platform></platform>', xml_str)
    xml_str = re.sub(r'<refassem_short_name>.+?</refassem_short_name>', '<refassem_short_name></refassem_short_name>', xml_str)

    xml_str = re.sub(r'<STUDY_REF .+?/>', '<STUDY_REF/>', xml_str)
    xml_str = re.sub(r'<ANALYSIS_SET .+?>', '<ANALYSIS_SET>', xml_str)
    xml_str = re.sub(r'<ANALYSIS .+?>', '<ANALYSIS>', xml_str)
    xml_str = re.sub(r'<EXPERIMENT_SET .+?>', '<EXPERIMENT_SET>', xml_str)
    xml_str = re.sub(r'<RUN_SET .+?>', '<RUN_SET>', xml_str)
    xml_str = re.sub(r'<analysis_detail_uri>.+?</analysis_detail_uri>', '<analysis_detail_uri></analysis_detail_uri>', xml_str)
    xml_str = re.sub(r'<analysis_submission_uri>.+?</analysis_submission_uri>', '<analysis_submission_uri></analysis_submission_uri>', xml_str)
    xml_str = re.sub(r'<analysis_data_uri>.+?</analysis_data_uri>', '<analysis_data_uri></analysis_data_uri>', xml_str)


    # we need to take care of xml properties in different order but effectively/semantically the same
    effective_eq_xml = json.dumps(xmltodict.parse(xml_str).get('ResultSet').get('Result'), indent=4, sort_keys=True)

    effective_xml_md5sum = hashlib.md5(effective_eq_xml).hexdigest()
    return effective_xml_md5sum


def download_metadata_xml(gnos_id):
    
    url = 'https://gtrepo-bsc.annailabs.com/cghub/metadata/analysisFull/' + gnos_id
    response = None
    try:
        response = requests.get(url, stream=True, timeout=15)
    except:
        click.echo('Error: Unable to download metadata from %s' % url, err=True)
        sys.exit(0)

    if not response or not response.ok:
        click.echo('Error: Unable to download metadata from %s' % url, err=True)
        sys.exit(0)
    else:
        metadata_xml_str = response.text

    return metadata_xml_str

def find_cached_metadata_xml(gnos_id):

    metadata_xml_files = 'gnos_metadata/__all_metadata_xml/bsc/' + gnos_id + '__live__*.xml'
    if not glob.glob(metadata_xml_files): 
        click.echo('Warning: missing cached GNOS metadata xml in BSC for gnos_id: %s' % gnos_id, err=True)
        sys.exit(0)
    metadata_xml_file = sorted(glob.glob(metadata_xml_files))[-1]
    with open (metadata_xml_file, 'r') as x: data = x.read()
    return data


def collect_gnos_xml(donors_list, gnos_sample_ids_to_be_included, gnos_sample_ids_to_be_excluded, project, ega_dir, pcawg_gnos_id_sheet, workflow, annotations, force):
    
    for w in workflow:
        print('\nCollecting the GNOS xmls for workflow: {} of project: {}'.format(get_mapping(w), project))
        gnos_xml_dir = os.path.join(ega_dir, project, get_mapping(w), 'GNOS_xml')
        if not os.path.exists(gnos_xml_dir): os.makedirs(gnos_xml_dir)
        gnos_xml_sheet = []
        with open(pcawg_gnos_id_sheet, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if not row.get('donor_unique_id') in donors_list: continue
                gnos_id = row.get('gnos_id')
                if gnos_sample_ids_to_be_included and not gnos_id in gnos_sample_ids_to_be_included: continue
                if gnos_sample_ids_to_be_excluded and gnos_id in gnos_sample_ids_to_be_excluded: continue
                entry_type = row.get('entry_type')
                if not get_mapping(entry_type) == w: continue
                
                gnos_xml_gz_file = os.path.join(gnos_xml_dir, 'analysis.'+gnos_id+'.GNOS.xml.gz')
                if not os.path.exists(gnos_xml_gz_file+'.processed') or force:
                    latest_xml_str = download_metadata_xml(gnos_id)
                    latest_effective_xml_md5sum = effective_xml_md5sum(latest_xml_str)
                    cached_xml_str = find_cached_metadata_xml(gnos_id)
                    cached_effective_xml_md5sum = effective_xml_md5sum(cached_xml_str)
                    if not latest_effective_xml_md5sum == cached_effective_xml_md5sum:
                        click.echo('Warning: BSC gnos xml has different effective md5sum with the cached xml for gnos_id: %s' % gnos_id, err=True)
                        continue
                    
                    # gnos_xml_gz_file = os.path.join(gnos_xml_dir, 'analysis.'+gnos_id+'.GNOS.xml.gz')
                    with gzip.open(gnos_xml_gz_file, 'wb') as n:  n.write(cached_xml_str.encode('utf8'))
                xml_gz_md5sum = get_md5(gnos_xml_gz_file, False)
                gnos_xml = OrderedDict()
                gnos_xml['filename'] = gnos_id+'/analysis.'+gnos_id+'.GNOS.xml.gz.gpg'
                gnos_xml['checksum'] = annotations.get('xml_encrypted_checksum').get(gnos_xml['filename']) if annotations.get('xml_encrypted_checksum').get(gnos_xml['filename']) else None
                gnos_xml['unencrypted_checksum'] = xml_gz_md5sum
                gnos_xml_sheet.append(copy.deepcopy(gnos_xml))
        if gnos_xml_sheet:
            out_dir = os.path.join(ega_dir, 'file_info', 'GNOS_xml_file_info')
            if not os.path.isdir(out_dir): os.makedirs(out_dir)  
            staged_files = os.path.join(out_dir, project+'.'+w+'.tsv')
            write_tsv_file(gnos_xml_sheet, staged_files)


def get_md5(fname, use_shell=None):
    if use_shell:
        process = subprocess.Popen(
                'md5sum {}'.format(fname),
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )

        out, err = process.communicate()

        if process.returncode:
            # should not exit for just this error, improve it later
            print('Unable to run \'md5sum\' tool.\nError message: {}'.format(err))
            return None
        else:
            return out.split()[0]

    else:
        hash = hashlib.md5()
        with open(fname) as f:
            for chunk in iter(lambda: f.read(1024*256), ""):
                hash.update(chunk)
        return hash.hexdigest()

def get_mapping(source):
    workflow = {
       "normal_wgs_bwa_alignment": "bwa",
       "tumor_wgs_bwa_alignment": "bwa",
       "sanger_variant_calling": "sanger",
       "dkfz_embl_variant_calling": "dkfz",
       "broad_variant_calling": "broad",
       "muse_variant_calling": "muse",
       "broad_tar_variant_calling": "broad-tar",
       "normal_RNA_Seq_TOPHAT_bam": "tophat2",
       "normal_RNA_Seq_STAR_bam": "star",
       "tumor_RNA_Seq_TOPHAT_bam": "tophat2", 
       "tumor_RNA_Seq_STAR_bam": "star",
       "bwa": "analysis_alignment.PCAWG_WGS_BWA",
       "tophat2": "analysis_alignment.PCAWG_RNA-Seq_TopHat2",
       "star": "analysis_alignment.PCAWG_RNA-Seq_Star",
       "sanger": "analysis_variation.PCAWG_Variant_Sanger",
       "dkfz": "analysis_variation.PCAWG_Variant_Dkfz-embl",
       "broad": "analysis_variation.PCAWG_Variant_Broad",
       "muse": "analysis_variation.PCAWG_Variant_Muse"
    }
    return workflow.get(source)


def write_tsv_file(report, filename):
    with open(filename, 'w') as fh:
        header = True  
        for r in report:
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



def read_annotations(annotations, type, file_name):
    if not os.path.isfile(file_name):
        return
    with open(file_name, 'r') as r:
        if annotations.get(type): # reset annotation if exists
            del annotations[type]               
        
        if type == 'gender':
            annotations[type] = {}
            reader = csv.DictReader(r, delimiter='\t')
            for row in reader:
                if row.get('exist_gender_discrepancy') == 'True': continue
                annotations[type][row.get('donor_unique_id')] = row.get('gender') if row.get('gender') else None
        elif type == 'ega':
            annotations[type] = set()
            reader = csv.DictReader(r, delimiter='\t')
            for row in reader:
                if row.get('checksum') and row.get('unencrypted_checksum'):
                    annotations[type].add(row.get('filename'))
        elif type == 'project':
            annotations[type] = set()
            reader = csv.DictReader(r, delimiter='\t')
            for row in reader:
                annotations[type].add(row.get('dcc_project_code'))
        elif type == 'xml_encrypted_checksum':
            annotations[type] = {}
            reader = csv.DictReader(r, delimiter='\t')
            for row in reader:
                if row.get('path') and row.get('encrypted MD5'):
                    annotations[type][row.get('path')] = row.get('encrypted MD5')

        else:
            logger.warning('unknown annotation type: {}'.format(type))
    return annotations


def generate_exclude_list(file_pattern, gnos_sample_ids_to_be_excluded):
    files = glob.glob(file_pattern)
    for f in files:
        if f.endswith('.xml'):
            fname = str.split(f, '/')[-1]
            gnos_id = str.split(fname, '.')[1]
            gnos_sample_ids_to_be_excluded.add(gnos_id)

        elif f.endswith('.tsv'):
            with open(f, 'r') as fn:
                reader = csv.DictReader(fn, delimiter='\t')
                for r in reader:
                    if r.get('icgc_sample_id'):
                        gnos_sample_ids_to_be_excluded.add(r.get('icgc_sample_id'))
                    if r.get('filename'):
                        gnos_sample_ids_to_be_excluded.add(r.get('filename').split('/')[0])

    return gnos_sample_ids_to_be_excluded

def generate_id_list(id_lists):
    ids_list = set()
    if id_lists:
        files = glob.glob(id_lists)
        for fname in files:
            with open(fname) as f:
                for d in f:
                    if d.startswith('#'): continue 
                    ids_list.add(d.rstrip())

    return ids_list

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


def get_formal_vcf_name(vcf):
    vcf_map = {
      "sanger": "sanger_variant_calling",
      "dkfz": "dkfz_embl_variant_calling",
      "embl": "dkfz_embl_variant_calling",
      "dkfz_embl": "dkfz_embl_variant_calling",
      "broad": "broad_variant_calling",
      "muse": "muse_variant_calling",
      "broad_tar": "broad_tar_variant_calling"
    }   

    return vcf_map.get(vcf)


def generate_unstaged_files(donors_list, project, ega_dir, unstage_type, annotations, es, es_index, gnos_sample_ids_to_be_excluded, ftp, ftp_gnos_ids, file_info):
    for dt in unstage_type:
        print('\nCheck the unstaging files for data_type: {} of project: {}'.format(dt, project))
        file_pattern = os.path.join(ega_dir, project, get_mapping(dt), 'analysis','analysis.*.receipt-*.xml')
        gnos_sample_ids_to_be_excluded = generate_exclude_list(file_pattern, gnos_sample_ids_to_be_excluded)

        missing_files = set()
        for donor_unique_id in donors_list:
            es_json = get_donor_json(es, es_index, donor_unique_id)
            if dt == 'bwa':
                analysis = es_json.get('wgs').get('normal_specimen').get('bwa_alignment')
                add_files(analysis, missing_files, annotations, gnos_sample_ids_to_be_excluded, ftp, ftp_gnos_ids, file_info)
                for aliquot in es_json.get('wgs').get('tumor_specimens'):        
                    analysis = aliquot.get('bwa_alignment')
                    add_files(analysis, missing_files, annotations, gnos_sample_ids_to_be_excluded, ftp, ftp_gnos_ids, file_info)

            elif dt == 'rna_seq':
                pass

            else:
                for aliquot in es_json.get('wgs').get('tumor_specimens'):  
                    if not aliquot.get(get_formal_vcf_name(dt)):
                        break                 
                    analysis = aliquot.get(get_formal_vcf_name(vcf))
                    add_files(analysis, missing_files, annotations, gnos_sample_ids_to_be_excluded, ftp, ftp_gnos_ids, file_info)
 
        out_dir = os.path.join(ega_dir, 'file_info', 'bulk_report_of_files_missed_on_ftp_server')
        if not os.path.isdir(out_dir): os.makedirs(out_dir)
        out_file = os.path.join(out_dir, project+'.'+dt+'.tsv')
        with open(out_file, 'w') as o: o.write('\n'.join(sorted(missing_files)))        


def add_files(analysis, missing_files, annotations, gnos_sample_ids_to_be_excluded, ftp, ftp_gnos_ids, file_info):
    if gnos_sample_ids_to_be_excluded and analysis.get('gnos_id') in gnos_sample_ids_to_be_excluded: return
    # check if ftp has the folder, if not then all the files are missing for this analysis_id
    if ftp_gnos_ids and analysis.get('gnos_id') not in ftp_gnos_ids:
        filename = os.path.join(analysis.get('gnos_id'), 'analysis.'+analysis.get('gnos_id')+'.GNOS.xml.gz.gpg')
        missing_files.add(filename)
        for f in analysis.get('files'):
            filename = os.path.join(analysis.get('gnos_id'), f.get('file_name')+'.gpg')
            missing_files.add(filename)
            missing_file_info = [filename, '', f.get('file_md5sum')]
            #with open(file_info, 'a+') as fh: fh.write('\t'.join(missing_file_info) + '\n')
        return
    # if ftp has the folder, check first if the files have the information in file_info, if not check if they exist in the ftp
    filename = os.path.join(analysis.get('gnos_id'), 'analysis.'+analysis.get('gnos_id')+'.GNOS.xml.gz.gpg')
    ftp_files = None
    if not filename in annotations.get('ega'):
        ftp_files = get_ftp_files(ftp, analysis.get('gnos_id'))
        if not ftp_files: print('\nWarning: FTP did not response for analysis_id: {}'.format(analysis.get('gnos_id')))
        if ftp_files and not filename in ftp_files: missing_files.add(filename)
    for f in analysis.get('files'):
        filename = os.path.join(analysis.get('gnos_id'), f.get('file_name')+'.gpg')
        if not filename in annotations.get('ega'):
            if not ftp_files: ftp_files = get_ftp_files(ftp, analysis.get('gnos_id'))
            if not ftp_files: print('\nWarning: FTP did not response for analysis_id: {}'.format(analysis.get('gnos_id')))
            if ftp_files and not filename in ftp_files: missing_files.add(filename)
            missing_file_info = [filename, '', f.get('file_md5sum')]
            #with open(file_info, 'a+') as fh: fh.write('\t'.join(missing_file_info) + '\n')
    return missing_files

def get_ftp_files(ftp, gnos_id):
    try:
        ftp_files = set(ftp.nlst(gnos_id))
    except:
        return

    return ftp_files


def main(argv=None):

    parser = ArgumentParser(description="Transfer Jobs Json Generator",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=True)
    parser.add_argument("-p", "--dcc_project_code", dest="dcc_project_code", nargs="*",
             help="Specify dcc_project_code for EGA submission", required=False)
    parser.add_argument("-r", "--specify ega directory", dest="ega_dir",
             help="Specify ega directory", required=False)

    parser.add_argument("-c", "--include_donor_id_lists", dest="include_donor_id_lists", 
             help="File(s) containing DONOR IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-d", "--exclude_donor_id_lists", dest="exclude_donor_id_lists", 
             help="File(s) containing DONOR IDs to be excluded, use filename pattern to specify the file(s)", required=False)

    parser.add_argument("-i", "--include_gnos_sample_id_lists", dest="include_gnos_sample_id_lists",
             help="File(s) containing GNOS/SAMPLE IDs to be excluded, use filename pattern to specify the file(s)", required=False)
    parser.add_argument("-x", "--exclude_gnos_sample_id_lists", dest="exclude_gnos_sample_id_lists", 
             help="File(s) containing GNOS/SAMPLE IDs to be excluded, use filename pattern to specify the file(s)", required=False)

    parser.add_argument("-u", "--unstage_type", dest="unstage_type", nargs="*",
             help="Specify data_type for reporting the unstaging files", required=False)
    parser.add_argument("-s", "--sequence_type", dest="seq", nargs="*",
             help="List sequence_type types[wgs, rna-seq]", required=False)
    parser.add_argument("-w", "--workflow", dest="workflow", nargs="*",
             help="List workflow types[bwa, sanger, dkfz, broad, muse, tophat2, star]", required=False)

    parser.add_argument("-f", "--force", dest="force", action="store_true", 
             help="Allow to regenerate the sample/xml files even if already exists", required=False)       



    args = parser.parse_args()
    metadata_dir = args.metadata_dir
    ega_dir = args.ega_dir
    ega_dir = ega_dir if ega_dir else '../pcawg-ega-submission'

    annotations = {}
    read_annotations(annotations, 'gender', ega_dir+'/annotation/pcawg_donors_gender_info.txt')
    read_annotations(annotations, 'ega', ega_dir+'/file_info/file_info.tsv')
    read_annotations(annotations, 'project', ega_dir+'/annotation/project_info.tsv')
    read_annotations(annotations, 'xml_encrypted_checksum', ega_dir+'/annotation/pancancer_xml_encrypted_checksum_2016_02_11.tsv')


    dcc_project_code = args.dcc_project_code
    dcc_project_code = list(dcc_project_code) if dcc_project_code else list(annotations.get('project'))

    include_donor_id_lists = args.include_donor_id_lists
    exclude_donor_id_lists = args.exclude_donor_id_lists

    include_gnos_sample_id_lists = args.include_gnos_sample_id_lists
    exclude_gnos_sample_id_lists = args.exclude_gnos_sample_id_lists

    unstage_type = args.unstage_type
    seq = args.seq
    workflow = args.workflow
    force = args.force 

    unstage_type = list(unstage_type) if unstage_type else []
    seq= list(seq) if seq else [] 
    workflow = list(workflow) if workflow else []  

    donor_id_to_be_incuded = generate_id_list(include_donor_id_lists)
    donor_id_to_be_excluded = generate_id_list(exclude_donor_id_lists)

    # pre-exclude gnos entries when this option is chosen
    gnos_sample_ids_to_be_excluded = generate_id_list(exclude_gnos_sample_id_lists)

    # # only process the gnos entries when this option is chosen
    gnos_sample_ids_to_be_included = generate_id_list(include_gnos_sample_id_lists) 

    # # remove the gnos_ids_to_be_include from gnos_sample_ids_to_be_excluded
    # gnos_sample_ids_to_be_excluded.difference_update(gnos_sample_ids_to_be_included) 

    # remove the gnos_sample_ids_to_be_excluded from gnos_sample_ids_to_be_included
    gnos_sample_ids_to_be_included.difference_update(gnos_sample_ids_to_be_excluded)


    es_index = 'pcawg_summary'
    es_type = "donor"
    es_host = 'localhost:9200'

    es = Elasticsearch([es_host])

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

    pcawg_sample_sheet = '../pcawg-operations/lists/sample_sheet/pcawg_sample_sheet.tsv'
    pcawg_gnos_id_sheet = '../pcawg-operations/data_releases/may2016/release_may2016_entry.tsv'
    file_info = ega_dir+'/file_info/file_info_missing.tsv'


    for project in dcc_project_code:
        donors_list = get_donors_list(es, es_index, project)
        donors_list.difference_update(donor_id_to_be_excluded)

        if unstage_type:
            # connect with the ftp
            ftp=ftplib.FTP('ftp.ega.ebi.ac.uk', 'ega-box-520', ega_box_token)
            ftp_gnos_ids = set(ftp.nlst())
            generate_unstaged_files(donors_list, project, ega_dir, unstage_type, annotations, es, es_index, gnos_sample_ids_to_be_excluded, ftp, ftp_gnos_ids, file_info) 

        if seq:
            file_pattern = os.path.join(ega_dir, project, 'sample', 'sample.'+project+'.*.tsv')
            gnos_sample_ids_to_be_excluded = generate_exclude_list(file_pattern, gnos_sample_ids_to_be_excluded)
            collect_sample(donors_list, gnos_sample_ids_to_be_included, gnos_sample_ids_to_be_excluded, project, ega_dir, pcawg_sample_sheet, seq, annotations)

        if workflow:
            collect_gnos_xml(donors_list, gnos_sample_ids_to_be_included, gnos_sample_ids_to_be_excluded, project, ega_dir, pcawg_gnos_id_sheet, workflow, annotations, force)


    return 0


if __name__ == "__main__":
    sys.exit(main())

