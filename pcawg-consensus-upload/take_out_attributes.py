#!/usr/bin/env python


import sys
import os
import xmltodict
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import copy
import glob


def main(argv=None):
    parser = ArgumentParser(description="Take out the useless ANALYSIS_ATTRIBUTE",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--work_dir", dest="work_dir",
             help="Directory name containing fixed variant call files", required=True)

    args = parser.parse_args()
    work_dir = args.work_dir

    for repo in ('osdc-icgc', 'osdc-tcga'):
        for f in glob.glob(os.path.join(work_dir, repo, '*/analysis.xml')):
            with open(f, 'r') as x: analysis_str = x.read()
            analysis_obj = xmltodict.parse(analysis_str)
            attributes = analysis_obj.get('ANALYSIS_SET').get('ANALYSIS').get('ANALYSIS_ATTRIBUTES').get('ANALYSIS_ATTRIBUTE')
            attributes_copy = copy.deepcopy(attributes)
            for attr in attributes_copy:
                if attr.get('TAG').startswith('ega') or attr.get('TAG').startswith('icgc'):
                    attributes.remove(attr)

                else:
                    pass  # all others, leave it unchanged
            
            analysis_str_new = xmltodict.unparse(analysis_obj, pretty=True)
            with open(f, 'w') as x: x.write(analysis_str_new.encode('utf8'))

if __name__ == "__main__":
    sys.exit(main())
