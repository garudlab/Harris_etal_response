#!/bin/bash

# constant dadi1 model. Generate neutral simulations iwth this. 

file=$1

#Ne=1965135.93 ## this is for the Nb=0.002 model 
Ne=2317520
mu=`echo | awk '{print 10^-9}'`
locusLength=10000
theta=`echo | awk '{print 4*'$Ne'*'$mu'*'$locusLength'}'`
rho=`echo | awk '{print (4*'$Ne'*'$locusLength')*(5*10^(-9))}'`

outFile=~/Jensen_response/analysis/dadi1_recoded_neutrality

for i in `seq 1 10000`; do

./msms/bin/msms 145 1 -N $Ne -t $theta -r $rho $locusLength -eN 0.2745424 0.003327787 -eN 0.2747088 1.663894 > $file

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

python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_bps -w 401 -j 50 -d 0 -s 0.5 -b 10000 -l $locusLength -m $segsites

cat ${file}_cluster_bps >> ${outFile}_bps.txt

# compute Pi and S
python ~/Jensen_response/scripts/Pi_MS_TajD.py ${file}_MS ${file}_pi_s 145

cat ${file}_pi_s >> ${outFile}_Pi_S_TajD.txt

rm ${file}
rm ${file}_var
rm ${file}_cut
rm ${file}_MS
#rm ${file}_cluster_snps
rm ${file}_cluster_bps
rm ${file}_pi_s

done

