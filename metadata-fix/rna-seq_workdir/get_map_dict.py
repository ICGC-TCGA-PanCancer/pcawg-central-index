#!/usr/bin/env python

import sys
import csv
import os
import re
import xmltodict
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
import subprocess
from random import randint
import shutil
import util

#def util.write_file(flist, fn):
#    with open(fn, 'w') as f:
#        header = True  
#        for r in flist:
#            if header:
#                f.write('\t'.join(r.keys()) + '\n')
#                header = False 
#            # make the list of output from dict
#            line = []
#            for p in r.keys():
#                if isinstance(r.get(p), list):
#                    line.append('|'.join(r.get(p)))
#                elif isinstance(r.get(p), set):
#                    line.append('|'.join(list(r.get(p))))
#                elif r.get(p) is None:
#                    line.append('')
#                else:
#                    line.append(str(r.get(p)))
#            f.write('\t'.join(line) + '\n') 


fname = os.path.join(util.entry_dir, 'gnos_rnaseq_entry_uniq_mapping.txt')

aliquot_dict = {}
with open(fname, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
    for r in reader:
    	aliquot_id = r.get('aliquot_id')
    	if not aliquot_dict.get(aliquot_id): aliquot_dict[aliquot_id] = OrderedDict()
        submitter_donor_id = r.get('submitter_donor_id')
        if not aliquot_dict.get(aliquot_id).get(submitter_donor_id): aliquot_dict.get(aliquot_id)[submitter_donor_id] = OrderedDict()
        submitter_sample_id = r.get('submitter_sample_id')
        if not aliquot_dict.get(aliquot_id).get(submitter_donor_id).get(submitter_sample_id): aliquot_dict.get(aliquot_id).get(submitter_donor_id)[submitter_sample_id] =  r.get('submitter_specimen_id') if r.get('submitter_specimen_id') else None

            
fname = os.path.join(util.entry_dir, 'gnos_rnaseq_aliquot_to_fix.txt')
id_fixes_list = []
with open(fname, 'r') as f:
     reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
     for r in reader:
         if not r.get('dcc_project_code') == 'OV-AU': continue
         id_fix = OrderedDict()
         for t in ['dcc_project_code', 'aliquot_id', 'submitter_donor_id']:
             id_fix[t] = r.get(t)
         for t in ['dcc_specimen_type', 'submitter_sample_id', 'submitter_specimen_id']:
             for l in ['old', 'new']:
                 id_fix[l+'_'+t] = None
         if '|' in r.get('dcc_specimen_type'):
             types = r.get('dcc_specimen_type').split('|')
             for t in types:
                 if 'Primary' in t:       
                     id_fix['new_dcc_specimen_type'] = t
                 else:
                     id_fix['old_dcc_specimen_type'] = t              
         samples = r.get('submitter_sample_id').split('|')
         for s in samples:
             if 'AOCS-' in s: 
                     id_fix['new_submitter_sample_id']=s
                     id_fix['new_submitter_specimen_id']=aliquot_dict[r.get('aliquot_id')][r.get('submitter_donor_id')][s]
             else:
                 id_fix['old_submitter_sample_id']=s
                 id_fix['old_submitter_specimen_id']=aliquot_dict[r.get('aliquot_id')][r.get('submitter_donor_id')][s]
         id_fixes_list.append(copy.deepcopy(id_fix))       
 
util.write_file(id_fixes_list, util.id_fixes_dir+'/OV-AU_id_fixes.txt')

id_fixes_list = []
with open(fname, 'r') as f:
     reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
     for r in reader:
         if not r.get('dcc_project_code') == 'PACA-AU': continue
         id_fix = OrderedDict()
         for t in ['dcc_project_code', 'aliquot_id', 'dcc_specimen_type']:
             id_fix[t] = r.get(t)
         for t in ['submitter_donor_id', 'submitter_sample_id', 'submitter_specimen_id']:
             for l in ['old', 'new']:
                 id_fix[l+'_'+t] = None
         donors = r.get('submitter_donor_id').split('|')
         samples = r.get('submitter_sample_id').split('|')
         for s in donors:
             if 'ICGC' in s:
                 id_fix['new_submitter_donor_id']=s
                 samples.remove(r.get('aliquot_id'))
                 id_fix['new_submitter_sample_id']=samples[0]
                 id_fix['new_submitter_specimen_id']=aliquot_dict[r.get('aliquot_id')][s][id_fix.get('new_submitter_sample_id')] if aliquot_dict[r.get('aliquot_id')][s].get(id_fix.get('new_submitter_sample_id')) else 'missing' 
             else:
                 id_fix['old_submitter_donor_id']=s
                 id_fix['old_submitter_sample_id']=r.get('aliquot_id')
                 id_fix['old_submitter_specimen_id']='missing'
         id_fixes_list.append(copy.deepcopy(id_fix))

util.write_file(id_fixes_list, util.id_fixes_dir+'/PACA-AU_id_fixes.txt')


id_fixes_list = []
with open(fname, 'r') as f:
     reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
     for r in reader:
         if not r.get('dcc_project_code') == 'CLLE-ES': continue
         id_fix = OrderedDict()
         for t in ['dcc_project_code', 'aliquot_id', 'submitter_donor_id', 'submitter_sample_id', 'submitter_specimen_id']:
             id_fix[t] = r.get(t)
         for t in ['dcc_specimen_type']:
             for l in ['old', 'new']:
                 id_fix[l+'_'+t] = None
         specimen_types = r.get('dcc_specimen_type').split('|')
         for s in specimen_types:
             if 'Primary' in s:
                 id_fix['new_dcc_specimen_type']=s
             else:
                 id_fix['old_dcc_specimen_type']=s
         id_fixes_list.append(copy.deepcopy(id_fix))

util.write_file(id_fixes_list, util.id_fixes_dir+'/CLLE-ES_id_fixes.txt') 


id_fixes_list = []
with open(fname, 'r') as f:
       reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
       for r in reader:
           if not r.get('dcc_project_code') == 'LIRI-JP': continue
           id_fix = OrderedDict()
           for t in ['dcc_project_code', 'aliquot_id', 'submitter_donor_id', 'dcc_specimen_type', 'submitter_specimen_id']:
               id_fix[t] = r.get(t)
           for t in ['submitter_sample_id']:
               for l in ['old', 'new']:
                   id_fix[l+'_'+t] = None
           ids = r.get('submitter_sample_id').split('|')
           for s in ids:
               if 'RK' in s:
                   id_fix['new_submitter_sample_id']=s
               else:
                   id_fix['old_submitter_sample_id']=s
           id_fixes_list.append(copy.deepcopy(id_fix))
   
util.write_file(id_fixes_list, util.id_fixes_dir+'/LIRI-JP_id_fixes.txt')


id_fixes_list = []
with open(fname, 'r') as f:
       reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
       for r in reader:
           if not r.get('dcc_project_code') == 'CESC-US': continue
           id_fix = OrderedDict()
           for t in ['dcc_project_code', 'aliquot_id', 'submitter_donor_id', 'submitter_specimen_id']:
               id_fix[t] = r.get(t)
           for t in ['study', 'dcc_specimen_type', 'submitter_sample_id']:
               for l in ['old', 'new']:
                   id_fix[l+'_'+t] = None
           id_fix['old_study'] = 'icgc_pancancer'
           id_fix['new_study'] = 'PCAWG2.0'
           id_fix['old_dcc_specimen_type'] = 'Primary Tumour - solid tissue'
           id_fix['new_dcc_specimen_type'] = 'Primary tumour - solid tissue'
           id_fix['old_submitter_sample_id'] = 'b3b3a27c-ee9a-42af-a6d1-9af5970a98b9'
           id_fix['new_submitter_sample_id'] = '1e176d9d-dba9-4d41-946e-05b7f35eba64'
           id_fixes_list.append(copy.deepcopy(id_fix))
   
util.write_file(id_fixes_list, util.id_fixes_dir+'/CESC-US_id_fixes.txt')



fname = os.path.join(util.entry_dir, 'gnos_rnaseq_aliquot.txt')
id_fixes_list = []
with open(fname, 'r') as f:
     reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
     for r in reader:
         if not r.get('dcc_specimen_type') == 'Metastatic tumour - metastatsis to distant location': continue
         id_fix = OrderedDict()
         for t in ['dcc_project_code', 'aliquot_id', 'submitter_donor_id', 'submitter_specimen_id', 'submitter_sample_id']:
             id_fix[t] = r.get(t)
         id_fix['old_dcc_specimen_type'] = 'Metastatic tumour - metastatsis to distant location'
         id_fix['new_dcc_specimen_type'] = 'Metastatic tumour - metastasis to distant location'
         id_fixes_list.append(copy.deepcopy(id_fix))
 
util.write_file(id_fixes_list, util.id_fixes_dir+'/SKCM-US_id_fixes.txt')
