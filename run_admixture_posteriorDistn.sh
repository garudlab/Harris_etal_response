#!/bin/bash
# This script takes in an MS file  with 1000 simulations, and outputs a concatenated processed file. 
file=$1
rho_in=$2
adaptive_theta=$3
selection=$4
locusLength=$5

#rho_in=5e-9
#adaptive_theta=10
#selection=True
#file=tmp.txt

outFile=~/Jensen_response/analysis/Admixture_rho_${rho_in}_theta_${adaptive_theta}_selection_${selection}_posteriorDistn

for i in `seq 1 10000`; do

echo $i
python ~/Jensen_response/scripts/admixture_parameters_posteriorDistn.py $rho_in $adaptive_theta $selection $locusLength > ${file}_var

command=`cat ${file}_var | head -1` 
Nac=`cat ${file}_var | head -3 | tail -1` 
s=`cat ${file}_var | head -4 | tail -1` 
age=`cat ${file}_var | head -5 | tail -1` 

eval $command

PF=`cat $file | grep -B 1 segsites | head -1 | cut -f7`
segsites=`cat $file | grep segsites | cut -f2 -d' '`
echo 'PF is' $PF
echo 'segsites is' $segsites
echo 'Nac is' $Nac
echo 's is' $s
echo 'age is' $age

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
fi


# cluster 

if [ $selection == True ]
then
#    python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 0.5 -p $PF -n $Nac -e $s -a $age

    python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_bps -w 401 -j 50 -d 0 -s 0.5 -b 10000 -l $locusLength -p $PF -n $Nac -e $s -a $age -m $segsites

else
    python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 0.5

#    python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_bps -w 401 -j 50 -d 0 -s 0.5 -b 10000 -l $locusLength -m $segsites

fi

cat ${file}_cluster_snps >> ${outFile}_snps.txt
H12=`cat ${file}_cluster_snps | cut -f9`

echo -e "$H12\t$command" >> ${outFile}_H12_command.txt

#cat ${file}_cluster_bps >> ${outFile}_bps.txt


rm ${file}
rm ${file}_var
rm ${file}_cut
rm ${file}_MS
rm ${file}_cluster_snps
#rm ${file}_cluster_bps


done

