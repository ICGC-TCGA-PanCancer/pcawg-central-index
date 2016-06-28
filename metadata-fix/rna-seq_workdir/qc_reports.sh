#!/bin/bash

infile=$1

echo qc_reports all the rna-seq entries...
mkdir -p qc_reports
echo duplicated_uploading
less $infile|head -1 > qc_reports/duplicated_uploadings.txt && for f in `less $infile |awk -F"\t" '{print $8, $3}'|sort -u |awk '{print $1}'|sort |uniq -c|sort -rn|head -5|awk '{print $2}'`; do less $infile |grep $f >> qc_reports/duplicated_uploadings.txt; done
echo orphan_unaligned
less $infile|head -1 > qc_reports/orphan_unaligned_entry.txt && less $infile |grep 2bdc2df0-2a92-4a82-9833-4023e745ed21 >> qc_reports/orphan_unaligned_entry.txt
echo multiple_lane_entry
less $infile |head -1 > qc_reports/multi_laned_entry.txt
for f in `less $infile |awk -F"\t" '{print $1, $8, $9, $13}'|grep unaligned|sort -u|awk '{print $1}'|sort |uniq -c|sort -rn|awk '$1 > 1 '|awk '{print $2}'`; do less $infile |grep unaligned| grep $f >> qc_reports/multi_laned_entry.txt; done 
for f in `less $infile |awk -F"\t" '{print $1, $8, $9, $13}'|grep TopHat|sort -u|awk '{print $1}'|sort |uniq -c|sort -rn|awk '$1 > 1 '|awk '{print $2}'`; do less $infile |grep TopHat| grep $f >> qc_reports/multi_laned_entry.txt; done 
for f in `less $infile |awk -F"\t" '{print $1, $8, $9, $13}'|grep STAR|sort -u|awk '{print $1}'|sort |uniq -c|sort -rn|awk '$1 > 1 '|awk '{print $2}'`; do less $infile |grep STAR| grep $f >> qc_reports/multi_laned_entry.txt; done 
echo same_sample_with_different_aliquot_ids
less $infile|head -1 > qc_reports/same_sample_id_has_different_aliquot_ids.txt && for f in `less $infile |awk -F"\t" '{print $16, $1}'|sort -u |awk '{print $1}'|sort |uniq -c|sort -rn|awk '$1 > 1'|awk '{print $2}'`; do less $infile |grep -w $f|sort >> qc_reports/same_sample_id_has_different_aliquot_ids.txt; done

echo donors_have_multiple_specimen_ids
less $infile|head -1 > qc_reports/donors_have_multiple_specimen_ids.txt
for f in `less $infile |grep tumour|awk -F"\t" '{print $4, $15, $17}'|sort -u|awk '{print $1,$2}'|sort |uniq -c|sort -rn|awk '$1 > 1'|awk '{print $3}'`; do less $infile|awk -v a=$f -F"\t" '$15 == a'|sort >> qc_reports/donors_have_multiple_specimen_ids.txt; done

echo same_specimen_with_different_sample_ids
less $infile|head -1 > qc_reports/same_specimen_id_has_different_sample_ids.txt && for f in `less $infile |awk -F"\t" '{print $17, $16}'|sort -u |awk '{print $1}'|sort |uniq -c|sort -rn|awk '$1 > 1'|awk '{print $2}'`; do less $infile |awk -v a=$f -F"\t" '$17 == a'|sort >> qc_reports/same_specimen_id_has_different_sample_ids.txt; done

echo breakdown datasets by project, specimen_type, data_type
less $infile |awk -F"\t" '{print $4}'|grep -v dcc_project_code|sort |uniq -c |sort -rn > qc_reports/breakdown_by_project.txt

less $infile |awk -F"\t" '{print $5}'|grep -v dcc_specimen_type|sort |uniq -c |sort -rn > qc_reports/breakdown_by_specimen_type.txt
less $infile |awk -F"\t" '{print $7, $13, $18}'|sort|uniq -c |sort -n > qc_reports/breakdown_by_datatype.txt
less $infile|head -1 > qc_reports/has_non_CV_term_entry.txt && less $infile|grep 'peripheral blood'|grep -v tumour >> qc_reports/has_non_CV_term_entry.txt
less $infile |head -1 > qc_reports/with_TEST_in_filename_entry.txt
less $infile |grep TEST >> qc_reports/with_TEST_in_filename_entry.txt 
