#!/bin/bash
infile=$1
info=($(less $infile |head -1|awk -F"\t" '{for (i=1; i<=NF; i++) print $i}'))

echo report all the entries with mismatch info...

mkdir -p mismatch_info
ll=( 4 5 13 14 15 16)
for i in "${ll[@]}"; do
    num=`less $infile |awk -v k=$i -F"\t" '{print $1, $k}'|sort |uniq -c|awk '{print $2}'|sort |uniq -d |wc -l`
    z=0
    if [ "$num" -gt "$z" ]; then
       j=$(( $i - 1 ))
       echo ${info[j]} 
       less $infile |awk -v k=$i -F"\t" '{print $1, $k}'|sort |uniq -c|awk '{print $2}'|sort |uniq -d > tmp.txt
       ln=`grep -f tmp.txt $infile |wc -l` 
       less $infile |head -1 > mismatch_info/mismatch_"$ln"_"${info[j]}".txt
       grep -f tmp.txt $infile | sort  >> mismatch_info/mismatch_"$ln"_"${info[j]}".txt 
    fi

rm -f tmp.txt
done
