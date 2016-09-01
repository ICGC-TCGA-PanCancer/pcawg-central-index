#!/bin/bash

infile=$1
info=($(less $infile |head -1|awk -F"\t" '{for (i=1; i<=NF; i++) print $i}'))

echo report all the entries having missing info...
mkdir -p missing_info

for (( i=0; i<${#info[@]}; i++));
do
#   echo ${info[i]}
   let j=i+1
#   echo $j
#   echo $"$j"
   num=`awk -v k=$j -F"\t" '$k==""' $infile |wc -l`
#   echo $num
   z=0
   if [ "$num" -gt "$z" ]; then
      less $infile |head -1 > missing_info/missing_"$num"_${info[i]}.txt
      awk -v k=$j -F"\t" '$k==""' $infile >> missing_info/missing_"$num"_${info[i]}.txt
   fi
done

