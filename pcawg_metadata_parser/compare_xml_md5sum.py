#!/usr/bin/env python

#script to calculate and compare the xml md5sum for specific sections such as qc_metrics, bai/bam files
#between source_gnos_repo and target_gnos_repo

#input file has the format of: 
#reports/QC_reports/specimens_with_mismatch_effective_xml_md5sum.txt 
#-o could write to ordered xml files as well
#-i could specify the comparison results file
#-u could specify whether use the local cached xml or download the lastest xml for comparison

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


def download_metadata_xml(gnos_repo, ao_uuid):
    
    url = gnos_repo + '/cghub/metadata/analysisFull/' + ao_uuid
    response = None
    try:
        response = requests.get(url, stream=True, timeout=30)
    except: # download failed, no need to do anything
        pass

    if not response or not response.ok:
        print 'unable to download metadata'
        return 
    else:
        metadata_xml_str = response.text
        return metadata_xml_str


def find_cached_metadata_xml(metadata_dir, gnos_repo, gnos_id):
    analysis_object_file = metadata_dir+'/analysis_objects.'+get_formal_repo_name(gnos_repo)+'.tsv'
    with open(analysis_object_file, 'r') as a:
        for line in a:
            ao_id, state, updated = str.split(line.rstrip(), '\t')
            if not ao_id == gnos_id: continue
            metadata_xml_file = 'gnos_metadata/__all_metadata_xml/'+get_formal_repo_name(gnos_repo)+'/'+ gnos_id+'__'+state+'__'+updated+'.xml'
            with open (metadata_xml_file, 'r') as x: metadata_xml_str = x.read()
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

def get_analysis_attrib(gnos_analysis):
    analysis_attrib = {}
    if (not gnos_analysis['analysis_xml']['ANALYSIS_SET'].get('ANALYSIS')
          or not gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS'].get('ANALYSIS_ATTRIBUTES')
          or not gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS']['ANALYSIS_ATTRIBUTES'].get('ANALYSIS_ATTRIBUTE')
       ):
        return None
    for a in gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS']['ANALYSIS_ATTRIBUTES']['ANALYSIS_ATTRIBUTE']:
        if not analysis_attrib.get(a['TAG']):
            analysis_attrib[a['TAG']] = a['VALUE']

    return analysis_attrib

def get_file_info(file_fragment, file_type):
    file_info = {}
    if (type(file_fragment) != list): file_fragment = [file_fragment]

    for f in file_fragment:
        if not f.get('filename').endswith(file_type): continue
        file_info = f

    return file_info


def calculate_xml_md5sum(xml_str, fname, repo):
    md5sum = []
    gnos_analysis = xmltodict.parse(xml_str).get('ResultSet').get('Result')
    # take out the qc_metrics section
    analysis_attrib = get_analysis_attrib(gnos_analysis)
    qc_metrics_xml = json.dumps(json.loads(analysis_attrib.get('qc_metrics')) if analysis_attrib.get('qc_metrics') else [], indent=4, sort_keys=True)
    md5sum.append(hashlib.md5(qc_metrics_xml).hexdigest())

    # take out the different file_type section
    for file_type in ['.bam', '.bai']:
        file_xml_dict = get_file_info(gnos_analysis.get('files').get('file'), file_type)
        if file_xml_dict: 
            file_xml = json.dumps(file_xml_dict, indent=4, sort_keys=True)
            md5sum.append(hashlib.md5(file_xml).hexdigest())
        else:
            md5sum.append('missing')

    return md5sum


def main(argv=None):

    parser = ArgumentParser(description="Compare the effective xml md5sum among repos",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=True)
    parser.add_argument("-f", "--specimens_with_mismatch_effective_xml_md5sum", dest="fname",
             help="Specify file to process", required=False)
    parser.add_argument("-u", "--Whether download metadata xml", dest="download_xml",
             help="Specify whether download_metadata_xml or not", required=False)

    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports
    fname = args.fname
    download_xml = args.download_xml

    if not fname:
        fname = metadata_dir+'/reports/QC_reports/specimens_with_mismatch_effective_xml_md5sum.txt'
    
    #output file
    detail_result = fname.replace('.txt', '_details') +'.txt'
    if os.path.isfile(detail_result): os.remove(detail_result)

    with open(detail_result, 'w') as m:
        with open(fname, 'r') as f:
            for l in f:
                if l.startswith('donor_unique_id'): 
                    header = '\t'.join([l.rstrip('\n'), 'exist_effective_xml_mismatch', 'qc_metrics_md5sum', 'exist_qc_metrics_mismatch','bam_file_md5sum', \
                                                                       'exist_bam_file_mismatch', 'bai_file_md5sum', 'exist_bai_file_mismatch'])
                    m.write(header+'\n')
                    continue
                field_info = str.split(l.strip(), '\t')
                donor_unique_id = field_info[0]
                gnos_repo = field_info[6].split('|')
                gnos_id = field_info[7].split('|')[0]
                md5sum_effective = field_info[8].split('|')
                md5sum_qc_metrics = []
                md5sum_bam_file = []
                md5sum_bai_file = []
                for repo in gnos_repo:
                    if not download_xml:
                        xml_str = find_cached_metadata_xml(metadata_dir, repo, gnos_id)
                    else:
                        xml_str = download_metadata_xml(repo, gnos_id)
                    if not xml_str:
                        md5sum = ['unable_download_xml']*3
                    else:
                        md5sum = calculate_xml_md5sum(xml_str, gnos_id, repo)
                    md5sum_qc_metrics.append(md5sum[0])
                    md5sum_bam_file.append(md5sum[1])
                    md5sum_bai_file.append(md5sum[2])
                mismatch_effective = 'False' if len(set(md5sum_effective))==1 else 'True'
                mismatch_qc_metrics = 'False' if len(set(md5sum_qc_metrics))==1 else 'True'
                mismatch_bam_file = 'False' if len(set(md5sum_bam_file))==1 else 'True'
                mismatch_bai_file = 'False' if len(set(md5sum_bai_file))==1 else 'True'
                l_new = '\t'.join([l.rstrip('\n'), mismatch_effective, \
                  '|'.join(md5sum_qc_metrics), mismatch_qc_metrics, \
                  '|'.join(md5sum_bam_file), mismatch_bam_file, \
                  '|'.join(md5sum_bai_file), mismatch_bai_file])+'\n'
                m.write(l_new)

    #if os.path.isfile(fname): os.remove(fname)
    return 0


if __name__ == "__main__":
    sys.exit(main()) 
