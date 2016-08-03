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


fname = util.entry_dir+'/gnos_rnaseq_entry.txt'
aliquot_dict = {}
with open(fname, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
    for r in reader:
    	aliquot_id = r.get('aliquot_id')
    	if not aliquot_dict.get(aliquot_id): aliquot_dict[aliquot_id] = OrderedDict()
    	for t in ['study', 'dcc_project_code', 'dcc_specimen_type', 'submitter_donor_id', 'submitter_sample_id', 'submitter_specimen_id' ]:
    	    if not aliquot_dict[aliquot_id].get(t): aliquot_dict[aliquot_id][t] = set()
    	    aliquot_dict[aliquot_id][t].add(r.get(t)) if r.get(t) else aliquot_dict[aliquot_id][t].add('missing')

aliquot_list = []
aliquot_to_fix_list = []
aliquot_id_to_fix_list = set()
for key, value in aliquot_dict.iteritems():
        aliquot = OrderedDict()
	aliquot['aliquot_id'] =  key
	aliquot.update(value)
	need_fix = False
	for k, v in value.iteritems():
		if not len(v) == 1: need_fix = True
	if need_fix: 
            aliquot_to_fix_list.append(copy.deepcopy(aliquot))
            aliquot_id_to_fix_list.add(key)
        aliquot_list.append(copy.deepcopy(aliquot))

util.write_file(aliquot_list, util.entry_dir+'/gnos_rnaseq_aliquot.txt')

with open(util.entry_dir+'/gnos_rnaseq_aliquot_id_to_fix.txt', 'w') as f:
    f.write('\n'.join(aliquot_id_to_fix_list))

util.write_file(aliquot_to_fix_list, util.entry_dir+'/gnos_rnaseq_aliquot_to_fix.txt')

gnos_rnaseq_entry_to_fix = []
with open(fname, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)         
    for r in reader:
        if not r.get('aliquot_id') in aliquot_id_to_fix_list: continue
        gnos_rnaseq_entry_to_fix.append(r)

util.write_file(gnos_rnaseq_entry_to_fix, util.entry_dir+'/gnos_rnaseq_entry_to_fix.txt')

gnos_rnaseq_entry_mapping = []
for g in gnos_rnaseq_entry_to_fix:
    gnew = {}
    for t in ['aliquot_id', 'study', 'dcc_project_code', 'dcc_specimen_type', 'submitter_donor_id', 'submitter_sample_id', 'submitter_specimen_id' ]:
        gnew[t]=g.get(t)
    gnos_rnaseq_entry_mapping.append(gnew)

gnos_rnaseq_entry_uniq_mapping=[dict(t) for t in set([tuple(sorted(d.items())) for d in gnos_rnaseq_entry_mapping])]

util.write_file(gnos_rnaseq_entry_uniq_mapping, util.entry_dir+'/gnos_rnaseq_entry_uniq_mapping.txt')


