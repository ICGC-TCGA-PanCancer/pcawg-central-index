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
from operator import itemgetter
import csv
from collections import OrderedDict


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
    file_info = []
    if (type(file_fragment) != list): file_fragment = [file_fragment]

    for f in file_fragment:
        if not f.get('filename').endswith(file_type): continue
        file_info.append(f)
    
    if file_info: file_info.sort(key=itemgetter('filename'))

    return file_info


def calculate_xml_md5sum(xml_str, workflow, xml_dir, gnos_id, gnos_repo):
    md5sum = []
    gnos_analysis = xmltodict.parse(xml_str).get('ResultSet').get('Result')
    # take out the qc_metrics section
    analysis_attrib = get_analysis_attrib(gnos_analysis)
    qc_metrics_xml = json.dumps(json.loads(analysis_attrib.get('qc_metrics')) if analysis_attrib.get('qc_metrics') else [], indent=4, sort_keys=True)
    md5sum.append(hashlib.md5(qc_metrics_xml).hexdigest())

    # take out the different file_type section
    if workflow.endswith('alignment'):
        file_types = ['.bam', '.bai']
    else:
        file_types = ['.gz', '.tbi']
    for file_type in file_types:
        file_xml = get_file_info(gnos_analysis.get('files').get('file'), file_type)
        if file_xml: 
            file_str = json.dumps(file_xml, indent=4, sort_keys=True)
            md5sum.append(hashlib.md5(file_str).hexdigest())
        else:
            md5sum.append('missing')

    xml_str = re.sub(r'<ResultSet .+?>', '<ResultSet>', xml_str)
    #xml_str = re.sub(r'<analysis_id>.+?</analysis_id>', '<analysis_id></analysis_id>', xml_str)
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
    xml_str = re.sub(r'</STUDY_REF>', '', xml_str)
    xml_str = re.sub(r'<STUDY_REF .+?>', '<STUDY_REF/>', xml_str)
    xml_str = re.sub(r'<ANALYSIS_SET .+?>', '<ANALYSIS_SET>', xml_str)
    xml_str = re.sub(r'<ANALYSIS .+?>', '<ANALYSIS>', xml_str)
    xml_str = re.sub(r'<EXPERIMENT_SET .+?>', '<EXPERIMENT_SET>', xml_str)
    xml_str = re.sub(r'<RUN_SET .+?>', '<RUN_SET>', xml_str)
    xml_str = re.sub(r'<analysis_detail_uri>.+?</analysis_detail_uri>', '<analysis_detail_uri></analysis_detail_uri>', xml_str)
    xml_str = re.sub(r'<analysis_submission_uri>.+?</analysis_submission_uri>', '<analysis_submission_uri></analysis_submission_uri>', xml_str)
    xml_str = re.sub(r'<analysis_data_uri>.+?</analysis_data_uri>', '<analysis_data_uri></analysis_data_uri>', xml_str)

    # we need to take care of xml properties in different order but effectively/semantically the same
    effective_gnos_analysis = xmltodict.parse(xml_str).get('ResultSet').get('Result')
    effective_eq_xml = json.dumps(effective_gnos_analysis, indent=4, sort_keys=True)

    if xml_dir:
        with open(xml_dir+'/'+gnos_id+'_'+get_formal_repo_name(gnos_repo), 'w') as y:
            y.write(effective_eq_xml)

    # all other parts other than the qc_metric, data_file, index_file
    xml_str = re.sub(r'<VALUE>.*{"qc_metrics".+?</VALUE>', '<VALUE>{"qc_metrics"}</VALUE>', xml_str, re.DOTALL)
    other_xml = xmltodict.parse(xml_str).get('ResultSet').get('Result')
    other_xml.pop('files')
    other_str = json.dumps(other_xml, indent=4, sort_keys=True)

    md5sum.append(hashlib.md5(other_str).hexdigest()) 

    return md5sum

def generate_subreport(fname, subreport_dir):
    for subreport in ['qc_metrics', 'data_file', 'index_file', 'other']:
        for workflow in ['wgs_bwa', 'rna_seq']:
            subreport_name = subreport_dir+'/'+workflow+'_'+subreport+'_mismatch.txt'
            subreport_list = []
            with open(fname, 'r') as s:
                reader = csv.DictReader(s, delimiter='\t')
                for row in reader:
                    if row.get('workflow').startswith(workflow) and row.get('exist_'+subreport+'_mismatch') == 'True':
                        row_order = OrderedDict()
                        for fn in reader.fieldnames:
                            row_order[fn] = row.get(fn)
                        subreport_list.append(row_order)
                if not subreport_list: continue
                write_file(subreport_list, subreport_name)


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

    parser = ArgumentParser(description="Compare the effective xml md5sum among repos",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=True)
    parser.add_argument("-f", "--specimens_with_mismatch_effective_xml_md5sum", dest="fname",
             help="Specify file to process", required=False)
    parser.add_argument("-u", "--Whether download metadata xml", dest="download_xml",
             help="Specify whether download_metadata_xml or not", required=False)
    parser.add_argument("-o", "--ordered_xml output folder", dest="xml_dir",
             help="Specify output folder for the ordered_xml if needed", required=False)

    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports
    fname = args.fname
    download_xml = args.download_xml
    xml_dir = args.xml_dir

    if not fname:
        fname = metadata_dir+'/reports/QC_reports/specimens_with_mismatch_effective_xml_md5sum.txt'

    if xml_dir:
        if os.path.exists(xml_dir): shutil.rmtree(xml_dir, ignore_errors=True)  # empty the folder if exists
        os.makedirs(xml_dir)
    
    #output file
    detail_result = fname.replace('.txt', '_details') +'.txt'
    if os.path.isfile(detail_result): os.remove(detail_result)
    # subreport output folder
    subreport_name = re.sub(r'^compare_', '', os.path.basename(__file__))
    subreport_name = re.sub(r'\.py$', '', subreport_name)
    subreport_dir = os.path.dirname(detail_result)+'/'+subreport_name
    if os.path.exists(subreport_dir): shutil.rmtree(subreport_dir, ignore_errors=True)  # empty the folder if exists
    os.makedirs(subreport_dir)


    with open(detail_result, 'w') as m:
        with open(fname, 'r') as f:
            for l in f:
                if l.startswith('donor_unique_id'): 
                    header = '\t'.join([l.rstrip('\n'), 'exist_effective_xml_mismatch', 'qc_metrics_md5sum', 'data_file_md5sum', 'index_file_md5sum', 'other_md5sum',\
                        'exist_qc_metrics_mismatch', 'exist_data_file_mismatch', 'exist_index_file_mismatch', 'exist_other_mismatch'])
                    m.write(header+'\n')
                    continue
                field_info = str.split(l.strip(), '\t')
                #donor_unique_id = l.get('donor_unique_id')
                gnos_repo = field_info[6].split('|')
                gnos_id = field_info[7]
                md5sum_effective = field_info[8].split('|')
                md5sum_qc_metrics = []
                md5sum_data_file = []
                md5sum_index_file = []
                md5sum_other = []
                for repo in gnos_repo:
                    if not download_xml:
                        xml_str = find_cached_metadata_xml(metadata_dir, repo, gnos_id)
                    else:
                        xml_str = download_metadata_xml(repo, gnos_id)
                    if not xml_str:
                        md5sum = ['unable_download_xml']*3
                    else:
                        md5sum = calculate_xml_md5sum(xml_str, field_info[5], xml_dir, gnos_id, repo)
                    md5sum_qc_metrics.append(md5sum[0])
                    md5sum_data_file.append(md5sum[1])
                    md5sum_index_file.append(md5sum[2])
                    md5sum_other.append(md5sum[3])
                mismatch_effective = 'False' if len(set(md5sum_effective))==1 else 'True'
                mismatch_qc_metrics = 'False' if len(set(md5sum_qc_metrics))==1 and not 'missing' in md5sum_qc_metrics else 'True'
                mismatch_data_file = 'False' if len(set(md5sum_data_file))==1 and not 'missing' in md5sum_data_file else 'True'
                mismatch_index_file = 'False' if len(set(md5sum_index_file))==1 and not 'missing' in md5sum_index_file else 'True'
                mismatch_other = 'False' if len(set(md5sum_other))==1 else 'True'
                l_new = '\t'.join([l.rstrip('\n'), mismatch_effective, \
                  '|'.join(md5sum_qc_metrics), '|'.join(md5sum_data_file), '|'.join(md5sum_index_file), '|'.join(md5sum_other),\
                  mismatch_qc_metrics, mismatch_data_file, mismatch_index_file, mismatch_other])+'\n'
                m.write(l_new)

    #if os.path.isfile(fname): os.remove(fname)
    # generate subreports
    generate_subreport(detail_result, subreport_dir)

    return 0


if __name__ == "__main__":
    sys.exit(main()) 
