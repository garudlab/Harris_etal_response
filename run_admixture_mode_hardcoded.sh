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

outFile=~/Jensen_response/analysis/Admixture_mode_corrected_hardcoded_rho_${rho_in}_theta_${adaptive_theta}_selection_${selection}



Nac=4975360  # current population size of africa
mu=10^-9
theta=`echo | awk '{print 4*'$Nac'*'$mu'*'$locusLength'}'`
recombinationRate=`echo | awk '{print (4*'$Nac'*'$locusLength'*$rho_in)}'`


totalSampleSize=145
numberOfSimulations=1
sampleSizeAfrica=0
sampleSizeEurope=0
sampleSizeAmerica=145
scaledNeEurope=0.6276
#growthRateEurope=-2.674174e-05 #changed
growthRateEurope=532.197
scaledNeAmerica=3.2127
#growthRateAmerica=-0.006059296 #changed
growthRateAmerica=120588.6
scaledTimeAdmixture=7.263e-05
proportionAdmixture=0.85
scaledTimeSplitAfricaEurope=0.009798
scaledTimeAfricaExpansion=0.119150674524 # 0.1192009 
scaledNeAfricaBottleneck=0.0001246
scaledTimeAfricaCrash=0.1192009 # 0.1192512
scaledNeAfricaAncestral=1.049994




for i in `seq 1 10000`; do

echo $i

~/./msms/bin/msms  $totalSampleSize $numberOfSimulations -N $Nac -t $theta -I 3 $sampleSizeAfrica $sampleSizeEurope $sampleSizeAmerica -n 2 $scaledNeEurope -g 2 $growthRateEurope -n 3 $scaledNeAmerica -g 3 $growthRateAmerica -es $scaledTimeAdmixture 3 $proportionAdmixture -ej $scaledTimeAdmixture 3 2 -ej $scaledTimeAdmixture 4 1 -ej $scaledTimeSplitAfricaEurope 2 1 -en $scaledTimeAfricaExpansion 1 $scaledNeAfricaBottleneck -en $scaledTimeAfricaCrash 1 $scaledNeAfricaAncestral -r $recombinationRate $locusLength > $file


# proceed with analyzing the sweep                                   

# cut the file

lineNo=`cat $file | grep -n segsites | cut -f1 -d ':'`
(( lineNo = lineNo + 146 ))
cat $file | head -${lineNo} | tail -146 > ${file}_cut

PF=`cat $file | grep -B 1 segsites | head -1 | cut -f7`
segsites=`cat $file | grep segsites | cut -f2 -d' '`


#convert to MSMS format

if [ $segsites == 0 ]
then
python ~/Jensen_response/scripts/convertMS_noSegSites.py 145 ${file}_MS
else
python ~/Jensen_response/scripts/convertMS.py ${file}_cut ${file}_MS
python ~/Jensen_response/scripts/removeSingletons.py ${file}_MS ${file}_noSingletons
fi



# cluster 

if [ $selection == True ]
then
#    python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 0.5 -p $PF -n $Nac -e $s -a $age

    python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_bps -w 401 -j 50 -d 0 -s 0.5 -b 10000 -l $locusLength -p $PF -n $Nac -e $s -a $age -m $segsites

else
#    python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 0.5

    python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_bps -w 401 -j 50 -d 0 -s 0.5 -b 10000 -l $locusLength -m $segsites

fi

#cat ${file}_cluster_snps >> ${outFile}_snps.txt
cat ${file}_cluster_bps >> ${outFile}_bps.txt

# compute Pi and S
python ~/Jensen_response/scripts/Pi_MS_TajD.py ${file}_MS ${file}_pi_s

cat ${file}_pi_s >> ${outFile}_Pi_S_TajD.txt

# compute Pi and S without singletons
python ~/Jensen_response/scripts/Pi_MS_TajD.py ${file}_noSingletons ${file}_pi_s_noSingletons

cat ${file}_pi_s_noSingletons >> ${outFile}_Pi_S_TajD_noSingletons.txt 

rm ${file}
rm ${file}_var
rm ${file}_cut
rm ${file}_MS
#rm ${file}_cluster_snps
rm ${file}_cluster_bps
rm ${file}_pi_s
rm ${file}_pi_s_noSingletons


done

