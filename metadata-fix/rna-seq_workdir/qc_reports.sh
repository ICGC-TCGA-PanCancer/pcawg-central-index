#!/bin/bash

infile=$1

echo qc_reports all the rna-seq entries...
mkdir -p qc_reports

less $infile|head -1 > qc_reports/duplicated_uploadings.txt && for f in `less $infile |awk -F"\t" '{print $8, $3}'|sort -u |awk '{print $1}'|sort |uniq -c|sort -rn|head -5|awk '{print $2}'`; do less $infile |grep $f >> qc_reports/duplicated_uploadings.txt; done
less $infile|head -1 > qc_reports/orphan_unaligned_entry.txt && less $infile |grep 2bdc2df0-2a92-4a82-9833-4023e745ed21 >> qc_reports/orphan_unaligned_entry.txt
less $infile |head -1 > qc_reports/multi_laned_entry.txt
for f in `less $infile |awk -F"\t" '{print $1, $8, $9, $13}'|grep unaligned|sort -u|awk '{print $1}'|sort |uniq -c|sort -rn|awk '$1 > 1 '|awk '{print $2}'`; do less $infile |grep unaligned| grep $f >> qc_reports/multi_laned_entry.txt; done 
for f in `less $infile |awk -F"\t" '{print $1, $8, $9, $13}'|grep TopHat|sort -u|awk '{print $1}'|sort |uniq -c|sort -rn|awk '$1 > 1 '|awk '{print $2}'`; do less $infile |grep TopHat| grep $f >> qc_reports/multi_laned_entry.txt; done 
for f in `less $infile |awk -F"\t" '{print $1, $8, $9, $13}'|grep STAR|sort -u|awk '{print $1}'|sort |uniq -c|sort -rn|awk '$1 > 1 '|awk '{print $2}'`; do less $infile |grep STAR| grep $f >> qc_reports/multi_laned_entry.txt; done 
less $infile |awk -F"\t" '{print $4}'|grep -v dcc_project_code|sort |uniq -c |sort -rn > qc_reports/breakdown_by_project.txt
less $infile |awk -F"\t" '{print $5}'|grep -v dcc_specimen_type|sort |uniq -c |sort -rn > qc_reports/breakdown_by_specimen_type.txt
less $infile |awk -F"\t" '{print $7, $13, $18}'|sort|uniq -c |sort -n > qc_reports/breakdown_by_datatype.txt
less $infile|head -1 > qc_reports/has_non_CV_term_entry.txt && less $infile|grep 'peripheral blood'|grep -v tumour >> qc_reports/has_non_CV_term_entry.txt
less $infile |head -1 > qc_reports/with_TEST_in_filename_entry.txt
less $infile |grep TEST >> qc_reports/with_TEST_in_filename_entry.txt 
