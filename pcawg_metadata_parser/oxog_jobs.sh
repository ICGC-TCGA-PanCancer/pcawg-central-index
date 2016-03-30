#!/bin/bash

CLOUD=$1
PROJECT_FILE=$2
MAX_TIME=$3

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

COUNTER=0
while [  $COUNTER -lt $MAX_TIME ]; do

	for f in `cat $PROJECT_FILE`; do
		echo
		echo --------------------------------------
	    echo generate oxog jsons from $f
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

		echo generating the jsons from $f...
		ICGC_PROJECT_CODE=$f ./generate_oxog_whitelist_json.py -m $M -t $CLOUD

		JOB_NUM=`ls -l $M/reports/oxog_whitelist_json/ |grep json|wc -l`
		Z=0

		if [ "$JOB_NUM" -gt "$Z" ]; then
			echo move job to oxog-ops
			cp $M/reports/oxog_whitelist_json/*.json ../oxog-ops/oxog-$CLOUD-jobs/backlog-jobs/.
			rm -rf ${M}/reports/oxog_whitelist_json
	        
	        echo check the job into git
	        cd ../oxog-ops/oxog-$CLOUD-jobs/backlog-jobs/
	        git checkout master
	        git reset --hard origin/master
	        git pull
	        git add .
	        git commit -m "add $JOB_NUM new jobs for project: $f"        
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
