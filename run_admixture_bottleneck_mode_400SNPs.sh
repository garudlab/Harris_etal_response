#!/bin/bash

# admixture model using mode. Generate neutral simulations iwth this. 

file=$1

Nac=3100520
mu=`echo | awk '{print 10^-9}'`
locusLength=100000
theta=`echo | awk '{print 4*'$Nac'*'$mu'*'$locusLength'}'`
recombinationRate=`echo | awk '{print (4*'$Nac'*'$locusLength')*(5*10^(-9))}'`

outFile=~/Jensen_response/analysis/admixture_bottleneck_mode_corrected_neutrality.txt

totalSampleSize=145
numberOfSimulations=1
sampleSizeAfrica=0
sampleSizeEurope=0
sampleSizeAmerica=145
scaledNeEurope=0.7318321
#growthRateEurope=-6.857889e-05
growthRateEurope=850.52
scaledNeAmerica=2.968357
scaledTimeAdmixture=3.757037e-05
proportionAdmixture=0.871794
scaledNeEuropeAdmixtureBn=0.004446651
scaledTimeSplitAfricaEurope=0.006037894
scaledTimeAfricaExpansion=0.03233073
scaledNeAfricaBottleneck=0.1599473
scaledTimeAfricaCrash=0.03241136
scaledNeAfricaAncestral=1.0401




for i in `seq 1 10000`; do

~/msms/bin/msms $totalSampleSize $numberOfSimulations -N $Nac -t $theta -I 3 $sampleSizeAfrica $sampleSizeEurope $sampleSizeAmerica -n 2 $scaledNeEurope -g 2 $growthRateEurope -n 3 $scaledNeAmerica -es $scaledTimeAdmixture 3 $proportionAdmixture -en $scaledTimeAdmixture 3 $scaledNeEuropeAdmixtureBn -ej $scaledTimeAdmixture 3 2 -ej $scaledTimeAdmixture 4 1 -ej $scaledTimeSplitAfricaEurope 2 1 -en $scaledTimeAfricaExpansion 1 $scaledNeAfricaBottleneck -en $scaledTimeAfricaCrash 1 $scaledNeAfricaAncestral -r $recombinationRate $locusLength > $file

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
python ~/Jensen_response/scripts/removeSingletons.py ${file}_MS ${file}_noSingletons
fi


# cluster 

python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 0.5

#python ~/Jensen_response/scripts/H12_H2H1_simulations.py ${file}_MS 145 -o ${file}_cluster_bps -w 401 -j 50 -d 0 -s 0.5 -b 10000 -l $locusLength -m $segsites

cat ${file}_cluster_snps >> ${outFile}_snps.txt

# compute Pi and S
#python ~/Jensen_response/scripts/Pi_MS_TajD.py ${file}_MS ${file}_pi_s 145

#cat ${file}_pi_s >> ${outFile}_Pi_S_TajD.txt

# compute Pi and S without Singletons
#python ~/Jensen_response/scripts/Pi_MS_TajD.py ${file}_noSingletons ${file}_pi_s_noSingletons 145

#cat ${file}_pi_s_noSingletons >> ${outFile}_Pi_S_TajD_noSingletons.txt 


rm ${file}
rm ${file}_var
rm ${file}_cut
rm ${file}_MS
#rm ${file}_cluster_snps
rm ${file}_cluster_bps
rm ${file}_pi_s

done

