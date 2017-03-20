#!/usr/bin/env python

import yaml
import json
from elasticsearch import Elasticsearch
from parse_gnos_xml import *

repo = ''

metadata_dir = find_latest_metadata_dir('gnos_metadata')

header_printed = False

index_name = 'rna_seq_info'
index_type = 'gnos_entries'

es = Elasticsearch([ 'localhost:9200' ])
es.indices.delete( index_name, ignore=[400, 404] )
es.indices.create( index_name, ignore=400 )
es.indices.put_mapping(index=index_name, doc_type=index_type,
    body={
      "dynamic":"true",
      "dynamic_templates":[
         { 
            "template_1":{
               "mapping":{
                  "index":"not_analyzed"
               },
               "match":"*",
               "match_mapping_type":"string"
            }
         }
      ],
      "_all":{
         "enabled":False
      },
      "_source":{
         "compress":True
      },
      "properties":{
         "analysis_id":{
            "type":"string",
            "index":"not_analyzed"
         }
      }
    })


def es_indexing(doc=None):
    if not doc: return
    global es
    global index_name
    es.index(index=index_name, doc_type=index_type,
             id=doc['gnos_server']+'.'+doc['analysis_id'],
             body=doc, timeout=90 )


def get_analysis_attrib(gnos_analysis):
    analysis_attrib = {}
    if (not gnos_analysis['analysis_xml']['ANALYSIS_SET'].get('ANALYSIS')
          or not gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS'].get('ANALYSIS_ATTRIBUTES')
          or not gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS']['ANALYSIS_ATTRIBUTES'].get('ANALYSIS_ATTRIBUTE')
       ):
        return None

    analysis_attrib_fragment = gnos_analysis['analysis_xml']['ANALYSIS_SET']['ANALYSIS']['ANALYSIS_ATTRIBUTES']['ANALYSIS_ATTRIBUTE']
    if (type(analysis_attrib_fragment) != list): analysis_attrib_fragment = [analysis_attrib_fragment]

    for a in analysis_attrib_fragment:
        if not analysis_attrib.get(a['TAG']):
            analysis_attrib[a['TAG']] = a['VALUE']

    return analysis_attrib


def process_one(gnos_analysis):
    global header_printed

    analysis_attrib = get_analysis_attrib(gnos_analysis)

    if not analysis_attrib: return

    if not gnos_analysis.get('library_strategy'): return
    if not 'RNA' in gnos_analysis.get('library_strategy'): return

    info = {}

    info['analysis_id'] = gnos_analysis.get('analysis_id')
    analysis_detail_uri = gnos_analysis.get('analysis_detail_uri', '')
    info['analysis_full_uri'] = analysis_detail_uri.replace('analysisDetail', 'analysisFull')
    info['gnos_server'] = info['analysis_full_uri'].split('/')[2].split('.')[0].replace('gtrepo-','')
    info['published_date'] = gnos_analysis.get('published_date')
    info['last_modified'] = gnos_analysis.get('last_modified')
    info['aliquot_id'] = gnos_analysis.get('aliquot_id')
    info['study'] = gnos_analysis.get('study')

    info['workflow_name'] = analysis_attrib.get('workflow_name', '')
    info['workflow_version'] = analysis_attrib.get('workflow_version', '')
    info['dcc_project_code'] = analysis_attrib.get('dcc_project_code', '')
    info['submitter_donor_id'] = analysis_attrib.get('submitter_donor_id', '')
    info['submitter_specimen_id'] = analysis_attrib.get('submitter_specimen_id', '')
    info['submitter_sample_id'] = analysis_attrib.get('submitter_sample_id', '')
    info['dcc_specimen_type'] = analysis_attrib.get('dcc_specimen_type', '')

    if ( info['dcc_project_code'] == "" or
         'test' in info['study'].lower() or
         info['submitter_donor_id'] == "" or info['submitter_donor_id'] == 'N/A' or
         info['submitter_specimen_id'] == "" or
         info['submitter_sample_id'] == "" or
         info['dcc_specimen_type'] == "" ):
        return

    files = []
    if isinstance(gnos_analysis.get('files').get('file'), dict):
        file_list = [gnos_analysis.get('files').get('file')]  
    elif isinstance(gnos_analysis.get('files').get('file'), list):
        file_list = gnos_analysis.get('files').get('file') 

    file_ext = set([])
    file_md5sum = ''
    for f in file_list:
        files.append({'file_name': f.get('filename'), 'file_size': f.get('filesize'), 'file_md5sum': f.get('checksum').get('#text')})

        ext = f.get('filename').split('.')[-1]
        if ext == "gz" or ext == "bam":
            file_md5sum = f.get('checksum').get('#text')

        file_ext.add(ext)

    filenames = [ f.get('file_name') for f in files ]

    info['file_md5sum'] = file_md5sum
    info['file_name'] = ','.join(filenames)
    info['file_count'] = str(len(files))
    info['file_ext'] = ','.join(sorted(file_ext))

    #print json.dumps(info)

    field_names = sorted(info.keys())
    if not header_printed:
        print '\t'.join(field_names)
        header_printed = True

    #fields = [info.get(f,'') for f in field_names]  # for somereason this does not work
    fields = []
    for f in field_names:
        fields.append(info.get(f) if info.get(f) else '')

    print '\t'.join(fields)


    info['files'] = files

    es_indexing(doc=info)


with open('settings.yml') as f:
    conf = yaml.safe_load(f)
    for r in conf.get('gnos_repos'):
        conf[r.get('base_url')] = r.get('repo_code')


for f in get_xml_files( metadata_dir, conf, repo ):
    f = conf.get('output_dir') + '/__all_metadata_xml/' + f
    gnos_analysis = get_gnos_analysis(f)
    #print (json.dumps(gnos_analysis)) # debug
    if gnos_analysis:
        #print 'processing xml file: {} ...'.format(f)

	process_one( gnos_analysis )
