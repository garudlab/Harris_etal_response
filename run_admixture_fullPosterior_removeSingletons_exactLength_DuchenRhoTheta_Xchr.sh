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

sample_size=37
(( sample_size_1 = sample_size + 1 ))

outFile=~/Jensen_response/analysis/Admixture_fullPosterior_rho_${rho_in}_theta_${adaptive_theta}_selection_${selection}_DuchenRhoTheta_exactLength_Xchr

for rep in `seq 1 100`; do

cat ~/Jensen_response/scripts/loc_AER_new_noHeader.txt | while read line; do
    locusLength=`echo $line| cut -d' ' -f5`
    mu=`echo $line| cut -d' ' -f6`
    rho_in=`echo $line| cut -d' ' -f7`
    
python ~/Jensen_response/scripts/admixture_parameters_posteriorDistn_DuchenRhoTheta_Xchr.py $rho_in $adaptive_theta $selection $locusLength $mu > ${file}_var

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
(( lineNo = lineNo + sample_size_1 ))
cat $file | head -${lineNo} | tail -${sample_size_1} > ${file}_cut



#convert to MSMS format

if [ $segsites == 0 ]
then
python ~/Jensen_response/scripts/convertMS_noSegSites.py $sample_size ${file}_MS
else
python ~/Jensen_response/scripts/convertMS.py ${file}_cut ${file}_MS
fi

python ~/Jensen_response/scripts/removeSingletons.py ${file}_MS ${file}_noSingletons


# cluster 

#if [ $selection == True ]
#then
#    python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 0.5 -p $PF -n $Nac -e $s -a $age

#    python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_bps -w 401 -j 50 -d 0 -s 0.5 -b 10000 -l $locusLength -p $PF -n $Nac -e $s -a $age -m $segsites

#else
#    python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 0.5

#    python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_bps -w 401 -j 50 -d 0 -s 0.5 -b 10000 -l $locusLength -m $segsites

#fi

#cat ${file}_cluster_snps >> ${outFile}_snps.txt
#cat ${file}_cluster_bps >> ${outFile}_bps.txt

python ~/Jensen_response/scripts/Pi_MS_TajD.py ${file}_MS ${file}_pi_s_singletons $sample_size

cat ${file}_pi_s_singletons >> ${outFile}_Pi_S_TajD_Singletons_DuchenRhoTheta.txt

python ~/Jensen_response/scripts/Pi_MS_TajD.py ${file}_noSingletons ${file}_pi_s_noSingletons $sample_size

cat ${file}_pi_s_noSingletons >> ${outFile}_Pi_S_TajD_noSingletons_DuchenRhoTheta.txt

rm ${file}
rm ${file}_var
rm ${file}_cut
rm ${file}_MS
#rm ${file}_cluster_snps
rm ${file}_cluster_bps
rm ${file}_noSingletons

done 

done
