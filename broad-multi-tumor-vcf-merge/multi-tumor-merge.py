#!/usr/bin/env python


import sys
import os
import re
import json
from collections import OrderedDict
from distutils.version import LooseVersion
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from copy import deepcopy
import xmltodict
import uuid
import requests
import glob


GNOS_SERVER = 'https://gtrepo-osdc-icgc.annailabs.com/'
RESULT_SETS = ['broad_vcf', 'broad_tar', 'muse_vcf']


def download_metadata(gnos_id, f):
    url = GNOS_SERVER + 'cghub/metadata/analysisFull/' + gnos_id
    response = None
    try:
        response = requests.get(url, timeout=30)
    except Exception as ex:
        sys.exit('Unable to download GNOS metadata xml for {}, error: {}'.format(gnos_id, ex))

    if not response or not response.ok:
        sys.exit('Unable to download GNOS metadata : {} from {}'.format(gnos_id, url))
    else:
        metadata_xml_str = response.content
        with open(f, 'w') as fh:  # write to metadata xml file now
            fh.write(metadata_xml_str.encode('utf8'))


def get_gnos_analysis_object(f):
    if not os.path.isfile(f):
        gnos_uuid = f.split('/')[-1].replace('.xml', '')
        download_metadata(gnos_uuid, f)

    with open (f, 'r') as x: xml_str = x.read()
    analysis_obj = xmltodict.parse(xml_str).get('ResultSet').get('Result').get('analysis_xml')
    return analysis_obj


def init_merge_obj(input_obj, input_obj_uuid, merged_gnos_object):
    merged_obj = deepcopy(input_obj)
    attrs = merged_obj.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_ATTRIBUTES').get('ANALYSIS_ATTRIBUTE')
    for attr in attrs:
        # reset to '' for these attributes
        if attr.get('TAG') in ('vm_instance_cores',
                                'vm_instance_mem_gb',
                                'vm_location_code',
                                'variant_qc_metrics',
                                'variant_timing_metrics'):
            attr['VALUE'] = ''
        elif attr.get('TAG') == 'related_file_subset_uuids':
            attr['VALUE'] = merged_gnos_object.get('related_file_subset_uuids')

    # add a new attribute to keep track of merge sources
    attrs.append(OrderedDict({
            'TAG': 'merged_from',
            'VALUE': input_obj_uuid
        }))

    return merged_obj


def merge_pipeline_info_json(info_type, merge_to, merge_from):
    json_key = ''
    if info_type == 'variant_pipeline_input_info':
        json_key = 'workflow_inputs'
    elif info_type == 'variant_pipeline_output_info':
        json_key = 'workflow_outputs'

    merge_to_obj = json.loads(merge_to).get(json_key, {})
    merge_from_obj = json.loads(merge_from).get(json_key, {})

    existing_specimens = set([ specimen.get('specimen') for specimen in merge_to_obj ])
    incoming_specimens = set([ specimen.get('specimen') for specimen in merge_from_obj ])
    new_specimens = incoming_specimens - existing_specimens

    for specimen in merge_from_obj:
        if specimen.get('specimen') in new_specimens:
            merge_to_obj.append(deepcopy(specimen))

    return json.dumps( {json_key: merge_to_obj} )


def merge_obj(merged_obj, input_obj, input_obj_uuid, merged_gnos_object):
    # merge ANALYSIS_SET.ANALYSIS.ANALYSIS_TYPE.REFERENCE_ALIGNMENT.RUN_LABELS.RUN
    input_runs = input_obj.get('ANALYSIS_SET').get('ANALYSIS').get(
        'ANALYSIS_TYPE').get('REFERENCE_ALIGNMENT').get('RUN_LABELS').get('RUN')
    if not type(input_runs) == list: input_runs = [ input_runs ]

    merged_runs = merged_obj.get('ANALYSIS_SET').get('ANALYSIS').get(
        'ANALYSIS_TYPE').get('REFERENCE_ALIGNMENT').get('RUN_LABELS').get('RUN')
    if not type(merged_runs) == list: merged_runs = [ merged_runs ]

    existing_read_group_labels = set([ run.get('@read_group_label') for run in merged_runs ])
    incoming_read_group_labels = set([ run.get('@read_group_label') for run in input_runs ])
    new_read_group_labels = incoming_read_group_labels - existing_read_group_labels

    for run in input_runs:
        if run.get('@read_group_label') in new_read_group_labels:
            merged_runs.append(deepcopy(run))


    # merge ANALYSIS_SET.ANALYSIS.ANALYSIS_TYPE.REFERENCE_ALIGNMENT.SEQ_LABELS.SEQUENCE
    input_seqs = input_obj.get('ANALYSIS_SET').get('ANALYSIS').get(
        'ANALYSIS_TYPE').get('REFERENCE_ALIGNMENT').get('SEQ_LABELS').get('SEQUENCE')
    if not type(input_seqs) == list: input_seqs = [ input_seqs ]

    merged_seqs = merged_obj.get('ANALYSIS_SET').get('ANALYSIS').get(
        'ANALYSIS_TYPE').get('REFERENCE_ALIGNMENT').get('SEQ_LABELS').get('SEQUENCE')
    if not type(merged_seqs) == list: merged_seqs = [ merged_seqs ]

    existing_seqs = set([ seq.get('@data_block_name') + '|' + seq.get('@accession') for seq in merged_seqs ])
    incoming_seqs = set([ seq.get('@data_block_name') + '|' + seq.get('@accession') for seq in input_seqs ])
    new_seqs = incoming_seqs - existing_seqs

    for seq in input_seqs:
        if seq.get('@data_block_name') + '|' + seq.get('@accession') in new_seqs:
            merged_seqs.append(deepcopy(seq))


    # merge ANALYSIS_SET.ANALYSIS.TARGETS.TARGET
    input_targets = input_obj.get('ANALYSIS_SET').get('ANALYSIS').get(
        'TARGETS').get('TARGET')
    if not type(input_targets) == list: input_targets = [ input_targets ]

    merged_targets = merged_obj.get('ANALYSIS_SET').get('ANALYSIS').get(
        'TARGETS').get('TARGET')
    if not type(merged_targets) == list: merged_targets = [ merged_targets ]

    existing_targets = set([target.get('@refname') for target in merged_targets])
    incoming_targets = set([target.get('@refname') for target in input_targets])
    new_targets = incoming_targets - existing_targets

    for target in input_targets:
        if target.get('@refname') in new_targets:
            merged_targets.append(deepcopy(target))


    # merge ANALYSIS_SET.ANALYSIS.DATA_BLOCK.FILES.FILE
    input_files = input_obj.get('ANALYSIS_SET').get('ANALYSIS').get(
        'DATA_BLOCK').get('FILES').get('FILE')
    if not type(input_files) == list: input_files = [ input_files ]

    merged_files = merged_obj.get('ANALYSIS_SET').get('ANALYSIS').get(
        'DATA_BLOCK').get('FILES').get('FILE')
    if not type(merged_files) == list: merged_files = [ merged_files ]


    existing_files = set([f.get('@filename') for f in merged_files])
    incoming_files = set([f.get('@filename') for f in input_files])
    new_files = incoming_files - existing_files

    for f in input_files:
        if f.get('@filename') in new_files:
            merged_files.append(deepcopy(f))


    # merge attributes
    input_attrs = input_obj.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_ATTRIBUTES').get('ANALYSIS_ATTRIBUTE')

    input_attrs_to_be_merged = {}
    for input_attr in input_attrs:
        if input_attr.get('TAG') in (
                                'variant_pipeline_input_info',
                                'variant_pipeline_output_info',
                                'variant_workflow_version'
                              ):
            input_attrs_to_be_merged[input_attr.get('TAG')] = \
                input_attr.get('VALUE')


    merged_attrs = merged_obj.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_ATTRIBUTES').get('ANALYSIS_ATTRIBUTE')
    for merged_attr in merged_attrs:
        if merged_attr.get('TAG') == 'merged_from':
            merged_attr['VALUE'] += ',' + input_obj_uuid
        elif merged_attr.get('TAG') == 'variant_workflow_version':
            if LooseVersion(input_attrs_to_be_merged.get('variant_workflow_version')) > \
                    LooseVersion(merged_attr.get('VALUE')):  # keep the new version
                merged_attr['VALUE'] = input_attrs_to_be_merged.get('variant_workflow_version')
        elif merged_attr.get('TAG') in (
                                'variant_pipeline_input_info',
                                'variant_pipeline_output_info'
                              ):
            merged_attr['VALUE'] = merge_pipeline_info_json(merged_attr.get('TAG'), merged_attr['VALUE'], 
                input_attrs_to_be_merged.get(merged_attr.get('TAG')))



def produce_output(donor, input_dir, result_set, merged_gnos_object):
    input_objs = os.listdir(os.path.join(input_dir, donor, result_set))
    gnos_output_dir = merged_gnos_object.get('obj_dir')

    merged_obj = {}
    for input_obj_uuid in input_objs:
        if not re.match(
                    r'[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}',
                    input_obj_uuid
                ):
            continue

        input_obj_dir = os.path.join(input_dir, donor, result_set, input_obj_uuid)
        if not os.path.isdir(input_obj_dir):
            continue

        #created symbolic links from input_obj_dir to gnos_output_dir
        for s in glob.glob(os.path.join(input_obj_dir, "*")):
            if s.endswith('.xml'): continue
            os.symlink(s, os.path.join(gnos_output_dir, os.path.basename(s)))

        input_obj = get_gnos_analysis_object(os.path.join(input_obj_dir, input_obj_uuid + '.xml'))

        if not merged_obj:
            merged_obj = init_merge_obj(input_obj, input_obj_uuid, merged_gnos_object)
        else:
            merge_obj(merged_obj, input_obj, input_obj_uuid, merged_gnos_object)

    #generate merged analysis xml file
    with open(gnos_output_dir+'/analysis.xml', 'w') as xml:
        xml.write(xmltodict.unparse(merged_obj, pretty=True))


def process_donor(donor, output_dir, input_dir):
    # create output dir if not exists
    if not os.path.isdir(output_dir): os.mkdir(output_dir)

    donor_dir = os.path.join(output_dir, donor)
    if os.path.isdir(donor_dir): os.rmdir(donor_dir)
    os.mkdir(donor_dir)

    merged_gnos_objects = {}
    uuids = [str(uuid.uuid4()), str(uuid.uuid4()), str(uuid.uuid4())]
    i = 0
    for result_set in RESULT_SETS:
        os.mkdir(os.path.join(donor_dir, result_set))

        merged_gnos_objects[result_set] = { 'obj_uuid': uuids[i] }

        if i == 0:
            merged_gnos_objects[result_set]['related_file_subset_uuids'] = \
                ','.join([uuids[1], uuids[2]])
        elif i == 1:
            merged_gnos_objects[result_set]['related_file_subset_uuids'] = \
                ','.join([uuids[0], uuids[2]])
        elif i == 2:
            merged_gnos_objects[result_set]['related_file_subset_uuids'] = \
                ','.join([uuids[0], uuids[1]])
        else:
            sys.exit('Should only have three types of result set for Broad VCF calls.')

        merged_gnos_objects[result_set]['obj_dir'] = os.path.join(donor_dir, result_set, uuids[i])
        os.mkdir(merged_gnos_objects[result_set]['obj_dir'])

        i += 1

    for result_set in merged_gnos_objects.keys():
        produce_output(donor, input_dir, result_set, merged_gnos_objects.get(result_set))



def main(argv=None):
    parser = ArgumentParser(description="Multi-tumor VCF merge",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input_dir", dest="input_dir",
             help="Directory name containing original multi-tumor VCF entries organized by donors", required=True)
    parser.add_argument("-o", "--output_dir", dest="output_dir", 
             help="Directory name containing output merged GNOS entries, default to current dir", required=False)

    args = parser.parse_args()
    input_dir = os.path.abspath(args.input_dir)
    output_dir = args.output_dir if args.output_dir else '.'

    donors = os.listdir( input_dir )
    for donor in donors:
        if not os.path.isdir(os.path.join(input_dir, donor)): continue
        process_donor(donor, output_dir, input_dir)


if __name__ == "__main__":
    sys.exit(main())