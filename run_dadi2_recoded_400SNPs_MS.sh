#!/bin/bash

# constant dadi2 model. Generate neutral simulations iwth this. 

file=$1

#Ne=1964759.10  ## this is for the Nb=0.4 model
Ne=2317075
mu=`echo | awk '{print 10^-9}'`
locusLength=100000
theta=`echo | awk '{print 4*'$Ne'*'$mu'*'$locusLength'}'`
rho=`echo | awk '{print (4*'$Ne'*'$locusLength')*(5*10^(-9))}'`

outFile=~/Jensen_response/analysis/dadi2_neutrality_recoded_MS

for i in `seq 1 10000`; do

~/./msdir/ms 145 1 -t $theta -r $rho $locusLength -eN 0.1178203 0.5891016 -eN 0.1590574 1.472754 > $file

segsites=`cat $file | grep segsites | cut -f2 -d' '`

# proceed with analyzing the sweep                                   

# cut the file

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

python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 0.5

cat ${file}_cluster_snps >> ${outFile}_snps.txt

#python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_bps -w 401 -j 50 -d 0 -s 0.5 -b 10000 -l $locusLength -m $segsites

#cat ${file}_cluster_bps >> ${outFile}_bps.txt

# compute Pi and S
python ~/Jensen_response/scripts/Pi_MS.py ${file}_MS ${file}_pi_s

cat ${file}_pi_s >> ${outFile}_Pi_S.txt


rm ${file}
rm ${file}_var
rm ${file}_cut
rm ${file}_MS
rm ${file}_cluster_snps
#rm ${file}_cluster_bps
rm ${file}_pi_s

done

