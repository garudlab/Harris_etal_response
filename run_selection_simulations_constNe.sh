#!/bin/bash
# This script takes in an MS file  with 1000 simulations, and outputs a concatenated processed file. 
file=$1
rho_in=$2
adaptive_theta=$3
selection=$4
locusLength=$5
s=$6
age_rand=$7
partial_frequency=$8

outFile=~/Jensen_response/analysis/selection_simulations_constNe_rho_${rho_in}_theta_${adaptive_theta}_selection_${s}_age_${age_rand}_partialfreq_${partial_frequency}_MS

for i in `seq 1 10`; do
echo $i

python ~/Jensen_response/scripts/selection_simulations_constNe.py $rho_in $adaptive_theta $selection $locusLength $s $age_rand $partial_frequency > ${file}_var

command=`cat ${file}_var | head -1` 

eval $command

segsites=`cat $file | grep segsites | cut -f2 -d' '`

# proceed with analyzing the sweep                                   
# cut the file

lineNo=`cat $file | grep -n segsites | cut -f1 -d ':'`
(( lineNo = lineNo + 146 ))
cat $file | head -${lineNo} | tail -146 > ${file}_cut


numLines=`cat ${file}_cut | egrep '(2|3|4|5|6|7|8|9|a|b|c|d|e|f|g|h|i|j|k|l|m|n|o|p|q|r|s|t|u|v|w|x|y|z)' | wc -l`

if [ $numLines == 1 ]
then
HS=1
else
HS=0
fi


#convert to MSMS format

if [ $segsites == 0 ]
then
python ~/Jensen_response/scripts/convertMS_noSegSites.py 145 ${file}_MS
else
python ~/Jensen_response/scripts/convertMS.py ${file}_cut ${file}_MS
python ~/Jensen_response/scripts/removeSingletons.py ${file}_MS ${file}_noSingletons
fi



# cluster 


	python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 0.5

	cat ${file}_cluster_snps >> ${outFile}_snps.txt


# compute Pi and S

    if [ "$locusLength" -eq "11000" ]; then
	python ~/Jensen_response/scripts/Pi_MS_TajD_perBp.py ${file}_MS ${file}_pi_s 145 $locusLength

	cat ${file}_pi_s >> ${outFile}_Pi_S_TajD_perBp.txt

	python ~/Jensen_response/scripts/Pi_MS_TajD_perBp.py ${file}_noSingletons ${file}_pi_s_noSingletons 145 $locusLength

	cat ${file}_pi_s_noSingletons >> ${outFile}_Pi_S_TajD_perBp_noSingletons.txt


    fi




    if [ "$locusLength" -eq "11000" ]; then

     # compute long range LD
	python ~/Jensen_response/scripts/LD_MAF_simulations.py ${file}_MS ${file}_LD_0.05_0.95 145 X 0.05 0.95 $locusLength 50

	cat ${file}_LD_0.05_0.95 >> ${outFile}_LD_0.05_0.95.txt

    fi

done

rm ${file}
rm ${file}_var
rm ${file}_cut
rm ${file}_MS
rm ${file}_cluster_snps
#rm ${file}_cluster_bps
rm ${file}_pi_s
rm ${file}_pi_s_noSingletons
rm ${file}_noSingletons
rm ${file}_LD_0.05_0.95

done





