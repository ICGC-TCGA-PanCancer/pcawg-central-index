#!/bin/bash

ICGC_TOKEN_CODE=$1
CLOUD=$2
PROJECT_FILE=$3
MAX_TIME=$4


DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

COUNTER=0
while [  $COUNTER -lt $MAX_TIME ]; do

	for f in `cat $PROJECT_FILE`; do
		echo
		echo --------------------------------------
	    echo generate $CLOUD transfer jsons from $f
	    echo Script location: $DIR
		cd $DIR

		M=`find gnos_metadata -maxdepth 1 -type d -regex 'gnos_metadata/20[0-9][0-9]-[0-9][0-9].*[0-9][0-9]_[A-Z][A-Z][A-Z]' | sort | tail -1`

		echo update s3-transfer-operations git submodule
		cd ../s3-transfer-operations/
		git checkout master
		git pull
		cd $DIR

		echo update ceph-transfer-operations git submodule
		cd ../ceph_transfer_ops/
		git checkout master
		git pull
		cd $DIR

		echo update oxog-ops git submodule
		cd ../oxog-ops/
		git checkout master
		git pull
		cd $DIR

		echo update pcawg-operations git submodule where white lists are maintained
		cd ../pcawg-operations/
		git pull
		cd $DIR

		echo generating the transfer jsons from $f...
		ICGC_TOKEN=$ICGC_TOKEN_CODE ICGC_PROJECT_CODE=$f ./generate_s3_transfer_json.py -m $M -t $CLOUD -s wgs -v sanger dkfz broad muse

		JOB_NUM=`ls -l $M/reports/s3_transfer_json_$CLOUD/ |grep json|wc -l`
		Z=0

		if [ "$JOB_NUM" -gt "$Z" ]; then
			echo move job to job folders
			if [ "$CLOUD" == "collab" ]; then
				cp $M/reports/s3_transfer_json_$CLOUD/*.WGS-BWA*.json ../ceph_transfer_ops/ceph-transfer-jobs-bwa/backlog-jobs/.
				cp $M/reports/s3_transfer_json_$CLOUD/*.Dkfz_embl-VCF.json ../ceph_transfer_ops/ceph-transfer-jobs-vcf-v1/backlog-jobs/.
				cp $M/reports/s3_transfer_json_$CLOUD/*.Muse-VCF.json ../ceph_transfer_ops/ceph-transfer-jobs-vcf-v1/backlog-jobs/.
				cp $M/reports/s3_transfer_json_$CLOUD/*.Sanger-VCF.json ../ceph_transfer_ops/ceph-transfer-jobs-vcf-v3/backlog-jobs/.
				cp $M/reports/s3_transfer_json_$CLOUD/*.Broad-VCF.json ../ceph_transfer_ops/ceph-transfer-jobs-vcf-v3/backlog-jobs/.
			else
				cp $M/reports/s3_transfer_json_$CLOUD/*.WGS-BWA*.json ../s3-transfer-operations/s3-transfer-jobs-bwa/backlog-jobs/.
				cp $M/reports/s3_transfer_json_$CLOUD/*.Dkfz_embl-VCF.json ../s3-transfer-operations/s3-transfer-jobs-vcf-v1/backlog-jobs/.
				cp $M/reports/s3_transfer_json_$CLOUD/*.Muse-VCF.json ../s3-transfer-operations/s3-transfer-jobs-vcf-v1/backlog-jobs/.
				cp $M/reports/s3_transfer_json_$CLOUD/*.Sanger-VCF.json ../s3-transfer-operations/s3-transfer-jobs-vcf-v3/backlog-jobs/.
				cp $M/reports/s3_transfer_json_$CLOUD/*.Broad-VCF.json ../s3-transfer-operations/s3-transfer-jobs-vcf-v3/backlog-jobs/.
			fi
			rm -rf ${M}/reports/s3_transfer_json_$CLOUD
	        
	        echo check the job into git
	        if [ "$CLOUD" == "collab" ]; then
		        cd ../ceph_transfer_ops
		    else
		    	cd ../s3-transfer-operations
		    fi
	        git checkout master
	        git reset --hard origin/master
	        git pull
	        git add .
	        git commit -m "add $JOB_NUM new transfer jobs for project: $f"        
	        git push     
	    else
	    	echo no new job could be generated from $f
		fi
		echo
	done

	let COUNTER=COUNTER+1 
	echo The counter is $COUNTER
	# echo sleeping an hour
	# sleep 1h
    
done
