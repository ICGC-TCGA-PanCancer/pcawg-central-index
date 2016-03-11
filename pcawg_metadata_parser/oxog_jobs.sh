#!/bin/bash

P="tcga"
D="0311"
batch="backlog-001"

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

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


# linda's token
./generate_oxog_whitelist_json.py -m $M -t $P

echo rename job folder
mv $M/reports/oxog_whitelist_json $M/reports/oxog_whitelist_json_$P
chmod 700 ${M}/reports/oxog_whitelist_json_$P

echo create the symlink
cd oxog_jobs_dir

ln -s ../${M}/reports/oxog_whitelist_json_${P}/ ${D}.${P}.jsons