#!/usr/bin/env python3

# Author: Junjun Zhang

import sys
import os
import xmltodict
import json
import copy
import logging
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from elasticsearch import Elasticsearch
from collections import OrderedDict
import datetime
import dateutil.parser


# create logger with 'spam_application'
logger = logging.getLogger('gnos parser')
# create console handler with a higher log level
ch = logging.StreamHandler()


def init_es(es_host, es_index):
    es = Elasticsearch([ es_host ])
    es.indices.delete( es_index, ignore=[400, 404] )
    es.indices.create( es_index, ignore=400 )

    # create mappings
    es_mapping = open('pancan.donor.mapping.json')
    es.indices.put_mapping(index=es_index, doc_type='donor', body=es_mapping.read())
    es_mapping.close()

    es_mapping = open('pancan.file.mapping.json')
    es.indices.put_mapping(index=es_index, doc_type='bam_file', body=es_mapping.read())
    es_mapping.close()

    return es


def process_gnos_analysis(gnos_analysis, donors, es_index, es, bam_output_fh):
    if ( not gnos_analysis.get('aliquot_id')
        ):
        logger.warning('ignore entry does not have aliquot_id, GNOS entry: {}'
                         .format(gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    if gnos_analysis.get('refassem_short_name') != 'unaligned' and gnos_analysis.get('refassem_short_name') != 'GRCh37':
        logger.warning('ignore entry that is aligned but not aligned to GRCh37: {}'
                         .format(gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return # completely ignore test gnos entries for now, this is the quickest way to avoid test interferes real data 

    analysis_attrib = get_analysis_attrib(gnos_analysis)
    if not analysis_attrib:
        logger.warning('ignore entry does not have ANALYSIS information, GNOS entry: {}'
                         .format(gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    if analysis_attrib.get('variant_workflow_name'): # this is test for VCF upload
        logger.warning('ignore entry that is VCF upload test, GNOS entry: {}'
                         .format(gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    if not analysis_attrib.get('dcc_project_code') or not analysis_attrib.get('submitter_donor_id'):
        logger.warning('ignore entry does not have dcc_project_code or submitter_donor_id, GNOS entry: {}'
                         .format(gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    donor_unique_id = analysis_attrib.get('dcc_project_code') + '::' + analysis_attrib.get('submitter_donor_id')
    if is_in_donor_blacklist(donor_unique_id):
        return

    if is_test(analysis_attrib, gnos_analysis):
        logger.warning('ignore test entry: {}'.format(gnos_analysis.get('analysis_detail_uri')))
        return # completely ignore test gnos entries for now, this is the quickest way to avoid test interferes real data 

    if (gnos_analysis.get('refassem_short_name') != 'unaligned'
          and not is_train_2_aligned(analysis_attrib, gnos_analysis)
        ):
        # TODO: we may create another ES index for obsoleted BAM entries
        # TODO: we will need more sophisticated check for handling BAMs that are flagged as aligned but
        #       treated as unaligned (this is actually the case for TCGA input BAM entries)
        logger.warning('ignore entry that is aligned but not by train 2 protocol: {}'
                         .format( gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull') ))
        return

    #TODO: put things above into one function

    if not donors.get(donor_unique_id):
        # create a new donor if not exist
        donors[ donor_unique_id ] = create_donor(donor_unique_id, analysis_attrib, gnos_analysis)

    else: # the donor this bam entry belongs to already exists
        # perform some comparison between existing donor and the info in the current bam entry
        if (donors[donor_unique_id].get('gnos_study') != gnos_analysis.get('study')):
            logger.warning( 'existing donor {} has study {}, but study in current gnos ao is {}'.
                            format( donor_unique_id,
                                    donors[donor_unique_id].get('gnos_study'),
                                    gnos_analysis.get('study') ) )
        # more such check may be added, no time for this now

    # now parse out gnos analysis object info to build bam_file doc
    bam_file = create_bam_file_entry(donor_unique_id, analysis_attrib, gnos_analysis)

    # this part can certainly be improved, the code is complex
    # TODO: fix here, unaligned BAM entire can also be used here to populate cooresponding section of a donor
    #       this way, we will know explicitly how many input_bam files for which specimen type in a donor

    if 'normal' in bam_file.get('dcc_specimen_type').lower(): # normal
        if donors.get(donor_unique_id).get('normal_specimen'): # normal specimen exists
            if donors.get(donor_unique_id).get('normal_specimen').get('aliquot_id') == gnos_analysis.get('aliquot_id'):
                if bam_file.get('is_aligned'):
                    if donors.get(donor_unique_id)['normal_specimen'].get('is_aligned'):
                        logger.warning('more than one normal aligned bam for donor: {}, entry in use: {}, additional entry found in: {}'
                              .format(donor_unique_id,
                                  donors.get(donor_unique_id).get('normal_specimen').get('gnos_metadata_url'),
                                  gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull')
                              )
                        )
                        if donors.get(donor_unique_id).get('normal_specimen').get('upload_date') < bam_file.get(
                                'upload_date'):
                            donors.get(donor_unique_id)['normal_specimen'].update(
                                prepare_aggregated_specimen_level_info(copy.deepcopy(bam_file))
                            )
                    else:
                        donors.get(donor_unique_id)['normal_specimen'].update(
                            prepare_aggregated_specimen_level_info(copy.deepcopy(bam_file))
                        )
            else:
                logger.warning('same donor: {} has different aliquot_id: {}, {} for normal specimen, entry in use: {}, additional entry found in {}'
                  .format(donor_unique_id,
                      donors.get(donor_unique_id).get('normal_specimen').get('aliquot_id'),
                      gnos_analysis.get('aliquot_id'),
                      donors.get(donor_unique_id).get('normal_specimen').get('gnos_metadata_url'),
                      gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull')
                  )
                )                
        else:
            donors.get(donor_unique_id)['normal_specimen'].update(
                prepare_aggregated_specimen_level_info(copy.deepcopy(bam_file))
            )

    else: # not normal
        donors.get(donor_unique_id).get('all_tumor_specimen_aliquots').add(bam_file.get('aliquot_id'))
        donors.get(donor_unique_id)['all_tumor_specimen_aliquot_counts'] = len(donors.get(donor_unique_id).get('all_tumor_specimen_aliquots'))
        if bam_file.get('is_aligned'):
            if donors.get(donor_unique_id).get('aligned_tumor_specimens'):
                if donors.get(donor_unique_id).get('aligned_tumor_specimen_aliquots').intersection(
                        bam_file.get('aliquot_id')
                    ): # multiple alignments for the same tumor aliquot_id
                    logger.warning('more than one tumor aligned bam for donor: {} with aliquot_id: {}, entry in use: {}, additional entry found in: {}'
                          .format(donor_unique_id,
                              bam_file.get('aliquot_id'),
                              donors.get(donor_unique_id).get('normal_specimen').get('gnos_metadata_url'),
                              gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull')
                        )
                    )
                else:
                    donors.get(donor_unique_id).get('aligned_tumor_specimens').append( copy.deepcopy(bam_file) )
                    donors.get(donor_unique_id).get('aligned_tumor_specimen_aliquots').add(bam_file.get('aliquot_id'))
                    donors.get(donor_unique_id)['aligned_tumor_specimen_aliquot_counts'] = len(donors.get(donor_unique_id).get('aligned_tumor_specimen_aliquots'))
            else:  # create the first element of the list
                donors.get(donor_unique_id)['aligned_tumor_specimens'] = [copy.deepcopy(bam_file)]
                donors.get(donor_unique_id).get('aligned_tumor_specimen_aliquots').add(bam_file.get('aliquot_id'))  # set of aliquot_id
                donors.get(donor_unique_id)['aligned_tumor_specimen_aliquot_counts'] = 1
                donors.get(donor_unique_id)['has_aligned_tumor_specimen'] = True


    bam_file.update( donors[ donor_unique_id ] )
    del bam_file['bam_files']
    del bam_file['normal_specimen']
    del bam_file['aligned_tumor_specimens']
    del bam_file['aligned_tumor_specimen_aliquots']
    del bam_file['aligned_tumor_specimen_aliquot_counts']
    del bam_file['has_aligned_tumor_specimen']
    del bam_file['all_tumor_specimen_aliquots']
    del bam_file['all_tumor_specimen_aliquot_counts']
    donors[donor_unique_id]['bam_files'].append( copy.deepcopy(bam_file) )

    # push to Elasticsearch
    es.index(index=es_index, doc_type='bam_file', id=bam_file['bam_gnos_ao_id'], body=json.loads( json.dumps(bam_file, default=set_default) ))
    bam_output_fh.write(json.dumps(bam_file, default=set_default) + '\n')


def set_default(obj):
    if isinstance(obj, datetime.datetime):
        return obj.isoformat()
    if isinstance(obj, set):
        return list(obj)
    raise TypeError


def prepare_aggregated_specimen_level_info(bam_file):
    specimen = copy.deepcopy(bam_file)
    # TODO: actual aggregation to be completed
    return specimen


def is_in_donor_blacklist(donor_unique_id):
    donor_blacklist = set([
            "LIHC-US::G1551",
            "LIHC-US::G15512",
            "TCGA_MUT_BENCHMARK_4::G15511",
            "TCGA_MUT_BENCHMARK_4::G15512"
        ])
    if donor_blacklist.intersection([donor_unique_id]):
        return True
    else:
        return False


def create_bam_file_entry(donor_unique_id, analysis_attrib, gnos_analysis):
    file_info = parse_bam_file_info(gnos_analysis.get('files').get('file'))
    bam_file = {
        "dcc_specimen_type": analysis_attrib.get('dcc_specimen_type'),
        "submitter_specimen_id": analysis_attrib.get('submitter_specimen_id'),
        "submitter_sample_id": analysis_attrib.get('submitter_sample_id'),
        "aliquot_id": gnos_analysis.get('aliquot_id'),
        "use_cntl": analysis_attrib.get('use_cntl'),
        "total_lanes": analysis_attrib.get('total_lanes'),

        "gnos_metadata_url": gnos_analysis.get('analysis_detail_uri').replace('analysisDetail', 'analysisFull'),
        "refassem_short_name": gnos_analysis.get('refassem_short_name'),
        "bam_gnos_ao_id": gnos_analysis.get('analysis_id'),
        "upload_date": dateutil.parser.parse(gnos_analysis.get('upload_date')),
        "published_date": dateutil.parser.parse(gnos_analysis.get('published_date')),
        "last_modified": dateutil.parser.parse(gnos_analysis.get('last_modified')),

        "bam_file_name": file_info.get('file_name'),
        "bam_file_size": file_info.get('file_size'),
        "md5sum": file_info.get('md5sum'),
    }

    # much more TODO for bam file info and alignment details
    if bam_file.get('refassem_short_name') == 'unaligned':
        bam_file['is_aligned'] = False
        bam_file['bam_type'] = 'Unaligned BAM'
        bam_file['alignment'] = None  # or initiate as empty object {}, depending on how ES searches it
    elif analysis_attrib.get('workflow_output_bam_contents') == 'unaligned': # this is actually BAM with unmapped reads
        bam_file['is_aligned'] = False
        bam_file['bam_type'] = 'Specimen level unmapped reads after BWA alignment'
        bam_file['alignment'] = None
    else:
        bam_file['is_aligned'] = True
        bam_file['bam_type'] = 'Specimen level aligned BAM'
        bam_file['alignment'] = get_alignment_detail(analysis_attrib, gnos_analysis)

    return bam_file


def get_alignment_detail(analysis_attrib, gnos_analysis):
    alignment = {
        "data_train": "Train 2",
        "workflow_name": analysis_attrib.get('workflow_name'),
        "workflow_version": analysis_attrib.get('workflow_version'),
        "workflow_bundle_url": analysis_attrib.get('workflow_bundle_url'),
        "workflow_source_url": analysis_attrib.get('workflow_source_url'),

        "pipeline_input_info": json.loads( analysis_attrib.get('pipeline_input_info') ).get('pipeline_input_info') if analysis_attrib.get('pipeline_input_info') else [],
        "qc_metrics": json.loads( analysis_attrib.get('qc_metrics').replace('"not_collected"', 'null') ).get('qc_metrics') if analysis_attrib.get('qc_metrics') else [],
        "markduplicates_metrics": json.loads( analysis_attrib.get('markduplicates_metrics') ).get('markduplicates_metrics') if analysis_attrib.get('markduplicates_metrics') else [],
        "timing_metrics": json.loads( analysis_attrib.get('timing_metrics').replace('"not_collected"', 'null') ).get('timing_metrics') if analysis_attrib.get('timing_metrics') else [],
    }

    alignment['input_bam_summary'] = {} # TODO: do this in a function
    
    return alignment


def parse_bam_file_info(file_fragment):
    file_info = {}
    if (type(file_fragment) != list): file_fragment = [file_fragment]

    for f in file_fragment:
        f = dict(f)
        if f.get('filename').endswith('.bam'): # assume there is only one BAM file
            file_info['file_name'] = f.get('filename')
            file_info['file_size'] = int(f.get('filesize'))
            file_info['md5sum'] = f.get('checksum').get('#text')

    return file_info


def is_train_2_aligned(analysis_attrib, gnos_analysis):
    if ( gnos_analysis.get('refassem_short_name') == 'GRCh37'
           and analysis_attrib.get('workflow_version')
           and analysis_attrib.get('workflow_version').startswith('2.6.')
       ):
        return True
    else:
        return False


def create_donor(donor_unique_id, analysis_attrib, gnos_analysis):
    donor = {
        'donor_unique_id': donor_unique_id,
        'submitter_donor_id': analysis_attrib['submitter_donor_id'],
        'dcc_project_code': analysis_attrib['dcc_project_code'],
        'gnos_study': gnos_analysis.get('study'),
        'gnos_repo': gnos_analysis.get('analysis_detail_uri').split('/cghub/')[0] + '/', # can be better
        'is_test': is_test(analysis_attrib, gnos_analysis),
        'is_cell_line': is_cell_line(analysis_attrib, gnos_analysis),
        'normal_specimen': {},
        'aligned_tumor_specimens': [],
        'aligned_tumor_specimen_aliquots': set(),
        'aligned_tumor_specimen_aliquot_counts': 0,
        'all_tumor_specimen_aliquots': set(),
        'all_tumor_specimen_aliquot_counts': 0,
        'has_aligned_tumor_specimen': False,
        'bam_files': []
    }
    if type(gnos_analysis.get('experiment_xml').get('EXPERIMENT_SET').get('EXPERIMENT')) == list:
        donor['sequencing_center'] = gnos_analysis.get('experiment_xml').get('EXPERIMENT_SET').get('EXPERIMENT')[0].get('@center_name')
    else:
        donor['sequencing_center'] = gnos_analysis.get('experiment_xml').get('EXPERIMENT_SET').get('EXPERIMENT').get('@center_name')

    return donor

def is_test(analysis_attrib, gnos_analysis):
    is_test = False
    if (gnos_analysis.get('aliquot_id') == '85098796-a2c1-11e3-a743-6c6c38d06053'
          or gnos_analysis.get('study') == 'CGTEST'
        ):
        is_test = True
    # TODO: what's the criteria for determining *test* entries

    return is_test


def is_cell_line(analysis_attrib, gnos_analysis):
    is_cell_line = False
    if analysis_attrib.get('dcc_project_code') == 'TCGA_MUT_BENCHMARK_4':
        is_cell_line = True

    return is_cell_line


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


def get_gnos_analysis(f):
    with open (f, 'r') as x:
        xml_str = x.read()
    return xmltodict.parse(xml_str).get('ResultSet').get('Result')


def process_dir(dir, es_index, es, donor_output_jsonl_file, bam_output_jsonl_file):
    donors = {}

    donor_fh = open(donor_output_jsonl_file, 'w')
    bam_fh = open(bam_output_jsonl_file, 'w')

    for f in os.listdir( dir ):
        if not f.endswith('.xml'): continue
        gnos_analysis = get_gnos_analysis(dir + f)
        #print (json.dumps(gnos_analysis)) # debug
        if gnos_analysis:
            logger.info( 'processing xml file: {}{} ...'.format(dir, f) )
            process_gnos_analysis( gnos_analysis, donors, es_index, es, bam_fh )
        else:
            logger.warning( 'skipping invalid xml file: {}{}'.format(dir, f) )

    for donor_id in donors.keys():
        #print(json.dumps(donors[donor_id])) # debug
        # push to Elasticsearch
        es.index(index=es_index, doc_type='donor', id=donors[donor_id]['donor_unique_id'], body=json.loads( json.dumps(donors[donor_id], default=set_default) ))
        del donors[donor_id]['bam_files']  # prune this before dumping JSON for Keiran
        donor_fh.write(json.dumps(donors[donor_id], default=set_default) + '\n')

    donor_fh.close()
    bam_fh.close()


def main(argv=None):
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    parser = ArgumentParser(description="PCAWG GNOS Metadata Parser",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--directory", dest="xml_dir",
             help="Directory where GNOS metadata XML files are included", required=True)
    parser.add_argument("-r", "--revision", dest="revision",
             help="A string for keeping revision information, will be part of the ES index name", required=True)

    args = parser.parse_args()
    xml_dir = args.xml_dir if args.xml_dir.endswith('/') else args.xml_dir + '/'
    revision = args.revision

    logger.setLevel(logging.INFO)
    ch.setLevel(logging.WARN)

    log_file = xml_dir + 'gnos_parser.log'
    # delete old log first if exists
    if os.path.isfile(log_file): os.remove(log_file)
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)

    es_host = 'localhost:9200'
    es_index = 'pancan_' + revision
    es = init_es(es_host, es_index)

    logger.info('processing xml files in {}'.format(xml_dir))
    process_dir(xml_dir, es_index, es, 'donor_'+es_index+'.jsonl', 'bam_'+es_index+'.jsonl')

    return 0


if __name__ == "__main__":
    sys.exit(main())

