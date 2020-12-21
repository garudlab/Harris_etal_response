#!/bin/bash

# constant Ne 10^6 model. Generate neutral simulations with this. 

model=$1
file=$2
outFile=$3
rho_in=$4
adaptive_theta=$5
selection=$6
locusLength=$7

if [ "$selection" == "True" ]; then
    MS_flag=False   
else
    MS_flag=True
fi


# read through the lengths of the short intron file. 

while read -r line;

do
    locusLength=$line
    echo $locusLength

    python ~/Jensen_response/scripts/generate_MS_commands.py $model $rho_in $adaptive_theta $selection $locusLength $MS_flag > ${file}_var

    command=`cat ${file}_var | head -1`
    eval $command

    segsites=`cat $file | grep segsites | cut -f2 -d' '`

    lineNo=`cat $file | grep -n segsites | cut -f1 -d ':'`
    (( lineNo = lineNo + 146 ))
    cat $file | head -${lineNo} | tail -146 > ${file}_cut

    #convert to MSMS format

    if [ $segsites == 0 ]
    then
	python ~/Jensen_response/scripts/convertMS_noSegSites.py 145 ${file}_MS
    else
	python ~/Jensen_response/scripts/convertMS.py ${file}_cut ${file}_MS
    fi


# compute Pi and S
# 11000 basepairs were simulated to avoid edge effects

python ~/Jensen_response/scripts/Pi_MS_TajD_perBp.py ${file}_MS ${file}_pi_s 145 $locusLength

cat ${file}_pi_s >> ${outFile}_Pi_S_TajD_perBp.txt


rm ${file}
rm ${file}_var
rm ${file}_cut
rm ${file}_MS
rm ${file}_pi_s

done < /u/home/n/ngarud/Jensen_response/data/DGRP_shortIntron/short_intron_lengths.txt

