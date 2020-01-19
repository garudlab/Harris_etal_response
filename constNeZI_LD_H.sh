#!/bin/bash
# This script takes in an MS file  with 1000 simulations, and outputs a concatenated processed file. 
file=$1
outFile1=$2
outFile2=$3
sampleSize=$4

Ne=4789108
mu=`echo | awk '{print 10^-9}'`
locusLength=110000
theta=`echo | awk '{print 4*'$Ne'*'$mu'*'$locusLength'}'`
rho=`echo | awk '{print (4*'$Ne'*'$locusLength')*(10^(-8))}'` 

for w in $(seq 1 1000); do

./msms/bin/msms $sampleSize 1 -N $Ne -t $theta -r $rho $locusLength > $file


        # cut the file                                                                                               
        
        lineNo=`cat $file | grep -n segsites | cut -f1 -d ':'`
	sampleSizeMod=$sampleSize
	(( sampleSizeMod = sampleSizeMod + 1 ))
        (( lineNo = lineNo + sampleSizeMod ))
        cat $file | head -${lineNo} | tail -$sampleSizeMod > ${file}_cut

        # convert to MSMS format                                                                                    
        python convertMS_noSweep2.py ${file}_cut ${file}_MS $locusLength

        # calculate LD for different MAF classes                            	
	python LD_MAF_simulations.py ${file}_MS ${file}_LD_0.5_0.6 $sampleSize X 0.5 0.6 $locusLength 50

	python LD_MAF_simulations.py ${file}_MS ${file}_LD_0.6_0.7 $sampleSize X 0.6 0.7 $locusLength 50

	python LD_MAF_simulations.py ${file}_MS ${file}_LD_0.7_0.8 $sampleSize X 0.7 0.8 $locusLength 50

	python LD_MAF_simulations.py ${file}_MS ${file}_LD_0.8_0.9 $sampleSize X 0.8 0.9 $locusLength 50

	python LD_MAF_simulations.py ${file}_MS ${file}_LD_0.9_1 $sampleSize X 0.9 1.0 $locusLength 50

	# calculate Homozygosity 
	python filterHaplotypes32_forSimulations.py ${file}_MS ${file}_H 200 401 2 2 0 $sampleSize

	cat ${file}_H >> $outFile1
	cat ${file}_LD_0.5_0.6 >> ${outFile2}_0.5_0.6
	cat ${file}_LD_0.6_0.7 >> ${outFile2}_0.6_0.7
	cat ${file}_LD_0.7_0.8 >> ${outFile2}_0.7_0.8
	cat ${file}_LD_0.8_0.9 >> ${outFile2}_0.8_0.9
	cat ${file}_LD_0.9_1 >> ${outFile2}_0.9_1

rm ${file}
rm ${file}_cut
rm ${file}_MS
rm ${file}_LD
rm ${file}_H

done

