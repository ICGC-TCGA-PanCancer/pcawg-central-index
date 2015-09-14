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
        else:
            logger.warning('duplicated analysis attribute key: {}'.format(a['TAG']))
    return analysis_attrib

def get_file_info(file_fragment, file_type):
    file_info = {}
    if (type(file_fragment) != list): file_fragment = [file_fragment]

    for f in file_fragment:
        if not f.get('filename').endswith(file_type): continue
        file_info = f

    return file_info


def calculate_xml_md5sum(xml_str, fname, repo, xml_dir):
    md5sum = []
    xml_str = re.sub(r'<ResultSet .+?>', '<ResultSet>', xml_str)
    xml_str = re.sub(r'<analysis_id>.+?</analysis_id>', '<analysis_id></analysis_id>', xml_str)
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
    #xml_str = re.sub(r'<VALUE>.*{"qc_metrics".+?</VALUE>', '<VALUE>\n{"qc_metrics"}</VALUE>', xml_str, re.DOTALL)


    xml_str = re.sub(r'<STUDY_REF .+?/>', '<STUDY_REF/>', xml_str)
    xml_str = re.sub(r'<ANALYSIS_SET .+?>', '<ANALYSIS_SET>', xml_str)
    xml_str = re.sub(r'<ANALYSIS .+?>', '<ANALYSIS>', xml_str)
    xml_str = re.sub(r'<EXPERIMENT_SET .+?>', '<EXPERIMENT_SET>', xml_str)
    xml_str = re.sub(r'<RUN_SET .+?>', '<RUN_SET>', xml_str)
    xml_str = re.sub(r'<analysis_detail_uri>.+?</analysis_detail_uri>', '<analysis_detail_uri></analysis_detail_uri>', xml_str)
    xml_str = re.sub(r'<analysis_submission_uri>.+?</analysis_submission_uri>', '<analysis_submission_uri></analysis_submission_uri>', xml_str)
    xml_str = re.sub(r'<analysis_data_uri>.+?</analysis_data_uri>', '<analysis_data_uri></analysis_data_uri>', xml_str)

    # we need to take care of xml properties in different order but effectively/semantically the same
    gnos_analysis = xmltodict.parse(xml_str).get('ResultSet').get('Result')
    effective_eq_xml = json.dumps(gnos_analysis, indent=4, sort_keys=True)
    md5sum.append(hashlib.md5(effective_eq_xml).hexdigest())

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

    if xml_dir: 
        with open(xml_dir + '/' + fname +'_ordered_'+get_formal_repo_name(repo), 'w') as y:
            y.write(effective_eq_xml)

    return md5sum


def get_gnos_repo(gnos_repo, source_gnos_repo, target_gnos_repo):
    if not source_gnos_repo or not target_gnos_repo:
        return gnos_repo

    if source_gnos_repo in gnos_repo and target_gnos_repo in gnos_repo:
        gnos_repo = [source_gnos_repo, target_gnos_repo]
    else:
        gnos_repo = []

    return gnos_repo 


def main(argv=None):

    parser = ArgumentParser(description="Compare the effective xml md5sum among repos",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=True)
    parser.add_argument("-f", "--gnos_id_file", dest="fname",
             help="Specify gnos_id file to process", required=True)
    parser.add_argument("-s", "--source_gnos_repo", dest="source_gnos_repo",
             help="Specify source_gnos_repo to do the comparation", required=False)
    parser.add_argument("-t", "--target_gnos_repo", dest="target_gnos_repo",
             help="Specify target_gnos_repo to do the comparation", required=False)
    parser.add_argument("-o", "--ordered_xml output folder", dest="xml_dir",
             help="Specify output folder for the ordered_xml", required=False)
    parser.add_argument("-i", "--Result of the comparation", dest="compare_result",
             help="Specify file for the compare result", required=False)
    parser.add_argument("-u", "--Whether download metadata xml", dest="download_xml",
             help="Specify whether download_metadata_xml or not", required=False)

    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports
    fname = args.fname
    source_gnos_repo = get_formal_repo_name(args.source_gnos_repo) if args.source_gnos_repo else None
    target_gnos_repo = get_formal_repo_name(args.target_gnos_repo) if args.source_gnos_repo else None
    xml_dir = args.xml_dir
    compare_result = args.compare_result
    download_xml = args.download_xml

    if xml_dir:
        if os.path.exists(xml_dir): shutil.rmtree(xml_dir, ignore_errors=True)  # empty the folder if exists
        os.makedirs(xml_dir)

    if not compare_result:
        compare_result = fname.replace('.txt', '') + '_' + get_formal_repo_name(source_gnos_repo) + '_' + get_formal_repo_name(target_gnos_repo) +'.txt'

    with open(compare_result, 'w') as m:
        header = ['donor_unique_id', 'gnos_repo', 'gnos_id', 'effective_xml_md5sum', 'mismatch_effective_xml', 'qc_metrics_xml_md5sum', 'mismatch_qc_metrics',\
                  'bam_file_xml_md5sum', 'mismatch_bam_file', 'bai_file_xml_md5sum', 'mismatch_bai_file']
        m.write('\t'.join(header)+'\n')
        with open(fname, 'r') as f:
            for l in f:
                if l.startswith('donor_unique_id'): continue
                field_info = str.split(l.strip(), '\t')
                donor_unique_id = field_info[0]
                gnos_repo = get_gnos_repo(field_info[6].split('|'), source_gnos_repo, target_gnos_repo)
                gnos_id = field_info[7].split('|')[0]
                md5sum_effective = []
                md5sum_qc_metrics = []
                md5sum_bam_file = []
                md5sum_bai_file = []
                if not gnos_repo: continue
                for repo in gnos_repo:
                    if not download_xml:
                        xml_str = find_cached_metadata_xml(metadata_dir, repo, gnos_id)
                    else:
                        xml_str = download_metadata_xml(repo, gnos_id)
                    if not xml_str:
                        md5sum = ['unable_download_xml']*4
                    else:
                        md5sum = calculate_xml_md5sum(xml_str, gnos_id, repo, xml_dir)
                    if not md5sum[0] in md5sum_effective: md5sum_effective.append(md5sum[0])
                    if not md5sum[1] in md5sum_qc_metrics: md5sum_qc_metrics.append(md5sum[1])
                    if not md5sum[2] in md5sum_bam_file: md5sum_bam_file.append(md5sum[2])
                    if not md5sum[3] in md5sum_bai_file: md5sum_bai_file.append(md5sum[3])
                mismatch_effective = 'False' if len(md5sum_effective)==1 else 'True'
                mismatch_qc_metrics = 'False' if len(md5sum_qc_metrics)==1 else 'True'
                mismatch_bam_file = 'False' if len(md5sum_bam_file)==1 else 'True'
                mismatch_bai_file = 'False' if len(md5sum_bai_file)==1 else 'True'
                l_new = '\t'.join([donor_unique_id, '|'.join(gnos_repo), gnos_id, \
                  '|'.join(md5sum_effective), mismatch_effective, \
                  '|'.join(md5sum_qc_metrics), mismatch_qc_metrics, \
                  '|'.join(md5sum_bam_file), mismatch_bam_file, \
                  '|'.join(md5sum_bai_file), mismatch_bai_file])+'\n'
                m.write(l_new)
    return 0



if __name__ == "__main__":
    sys.exit(main()) 
