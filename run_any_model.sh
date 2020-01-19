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

for i in `seq 1 1000`; do

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


    # cluster 

    if [ "$locusLength" -eq "10000" ]; then # test windows of length 10kb

	python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_bps -w 401 -j 50 -d 0 -s 0.5 -b 10000 -l $locusLength -m $segsites

	cat ${file}_cluster_bps >> ${outFile}_bps.txt

    elif [ "$locusLength" -eq "100000" ]; then # test windows of 400 SNP in longer simulated chromosomes. 

	python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 0.5

	cat ${file}_cluster_snps >> ${outFile}_snps.txt

    elif [ "$locusLength" -eq "350000" ]; then # simulations with selection done with longer chromosomes

	python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 0.5

	cat ${file}_cluster_snps >> ${outFile}_snps.txt

    fi


    # compute Pi and S
    # 11000 basepairs were simulated to avoid edge effects

    if [ "$locusLength" -eq "11000" ]; then
	python ~/Jensen_response/scripts/Pi_MS_TajD_perBp.py ${file}_MS ${file}_pi_s 145 $locusLength

	cat ${file}_pi_s >> ${outFile}_Pi_S_TajD_perBp.txt

    fi
    
     # compute long range LD

    if [ "$locusLength" -eq "11000" ]; then

	python ~/Jensen_response/scripts/LD_MAF_simulations.py ${file}_MS ${file}_LD_0.05_0.95 145 X 0.05 0.95 $locusLength 50

	cat ${file}_LD_0.05_0.95 >> ${outFile}_LD_0.05_0.95.txt

    fi

rm ${file}
rm ${file}_var
rm ${file}_cut
rm ${file}_MS
rm ${file}_cluster_snps
rm ${file}_cluster_bps
rm ${file}_pi_s
rm ${file}_LD_0.05_0.95
done

