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
import glob

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


def find_cached_metadata_xml(gnos_repo, gnos_id):

    metadata_xml_files = 'gnos_metadata/__all_metadata_xml/' + get_formal_repo_name(gnos_repo) + '/' + gnos_id + '__live__*.xml'
    metadata_xml_file = sorted(glob.glob(metadata_xml_files))[-1]
    with open (metadata_xml_file, 'r') as x: data = x.read()
    return data


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


def effective_xml_md5sum(xml_str):

    xml_str = re.sub(r'<ResultSet .+?>', '<ResultSet>', xml_str)
    #xml_str = re.sub(r'<analysis_id>.+?</analysis_id>', '<analysis_id></analysis_id>', xml_str)
    xml_str = re.sub(r'<last_modified>.+?</last_modified>', '<last_modified></last_modified>', xml_str)
    xml_str = re.sub(r'<upload_date>.+?</upload_date>', '<upload_date></upload_date>', xml_str)
    xml_str = re.sub(r'<published_date>.+?</published_date>', '<published_date></published_date>', xml_str)
    xml_str = re.sub(r'<analyte_code>.+?</analyte_code>', '<analyte_code></analyte_code>', xml_str)
    xml_str = re.sub(r'<reason>.+?</reason>', '<reason></reason>', xml_str)
    #xml_str = re.sub(r'<study>.+?</study>', '<study></study>', xml_str)
    xml_str = re.sub(r'<sample_accession>.+?</sample_accession>', '<sample_accession></sample_accession>', xml_str)
    #xml_str = re.sub(r'<dcc_project_code>.+?</dcc_project_code>', '<dcc_project_code></dcc_project_code>', xml_str)
    #xml_str = re.sub(r'<participant_id>.+?</participant_id>', '<participant_id></participant_id>', xml_str)
    #xml_str = re.sub(r'<specimen_id>.+?</specimen_id>', '<specimen_id></specimen_id>', xml_str)
    #xml_str = re.sub(r'<sample_id>.+?</sample_id>', '<sample_id></sample_id>', xml_str)
    xml_str = re.sub(r'<use_cntl>.+?</use_cntl>', '<use_cntl></use_cntl>', xml_str)
    xml_str = re.sub(r'<library_strategy>.+?</library_strategy>', '<library_strategy></library_strategy>', xml_str)
    xml_str = re.sub(r'<platform>.+?</platform>', '<platform></platform>', xml_str)
    xml_str = re.sub(r'<refassem_short_name>.+?</refassem_short_name>', '<refassem_short_name></refassem_short_name>', xml_str)

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


def generate_md5(data):
    data = re.sub(r'<ResultSet .+?>', '<ResultSet>', data)
    xml_md5 = hashlib.md5(data).hexdigest()

    return xml_md5


def main(argv=None):

    parser = ArgumentParser(description="Check the state for s3 jobs",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--jobs_folder", dest="jobs_folder",
             help="Specify jobs folder with JSONs", required=True)
    parser.add_argument("-s", "--jobs_state_file", dest="jobs_state_file",
             help="Specify output file for jobs_md5sum check state", required=True)
    parser.add_argument("-i", "--jobs_info_file", dest="jobs_info_file",
             help="Specify output file for jobs_info", required=False)

    args = parser.parse_args()

    jobs_folder = args.jobs_folder
    jobs_state_file = args.jobs_state_file
    jobs_info_file = args.jobs_info_file

    # read the info in job folder
    files = glob.glob(jobs_folder.rstrip('/') + '/*.json')
    jobs_info = set()

    if os.path.isfile(jobs_state_file): os.remove(jobs_state_file)
    
    with open(jobs_state_file, 'w') as i:
        for f in files:
            with open(f, 'r') as r:
                name_list = f.split('.')
                sub_json_name = '.'.join(name_list[3:])
                job = json.loads(r.read())
                gnos_repo = job.get('gnos_repo')[0]
                gnos_id = job.get('gnos_id')
                latest_xml_str = download_metadata_xml(gnos_repo, gnos_id)
                latest_effective_xml_md5sum = effective_xml_md5sum(latest_xml_str)
                latest_xml_md5sum = generate_md5(latest_xml_str)
                cached_xml_str = find_cached_metadata_xml(gnos_repo, gnos_id)
                cached_effective_xml_md5sum = effective_xml_md5sum(cached_xml_str)
                cached_xml_md5sum = generate_md5(cached_xml_str)
                do_effective_xml_md5sum_equal = False
                do_xml_md5sum_equal = False
                if latest_effective_xml_md5sum == cached_effective_xml_md5sum:
                    do_effective_xml_md5sum_equal = True
                if latest_xml_md5sum == cached_xml_md5sum:
                    jobs_info.add(gnos_id)
                    jobs_info.add(sub_json_name)
                    do_xml_md5sum_equal = True
                i.write('\t'.join([gnos_id+'.'+sub_json_name, gnos_repo, str(do_effective_xml_md5sum_equal), str(do_xml_md5sum_equal)])+'\n')


    # write the job info to file if specify the file name
    if jobs_info_file is not None:    
        # delete old file first if exists
        if os.path.isfile(jobs_info_file): os.remove(jobs_info_file)

        with open(jobs_info_file, 'w') as m:
            for g in jobs_info:
                m.write(g+'\n')          

    return 0


if __name__ == "__main__":
    sys.exit(main()) 
