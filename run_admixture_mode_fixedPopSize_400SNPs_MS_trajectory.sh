#!/bin/bash
# This script takes in an MS file  with 1000 simulations, and outputs a concatenated processed file. 

rho_in=5e-9
selection=True
locusLength=100
EUr=500000


for adaptive_theta in 0.01 10; do
for NAm in 61659 1110000; do


for i in `seq 1 100`; do

outFile=~/Jensen_response/analysis/Admixture_mode_fixedPopSize_NAm${NAm}_EUr${EUr}_rho_${rho_in}_theta_${adaptive_theta}_selection_${selection}_MS_trajectory_randomstart_${i}.txt

echo $i

file=test

python ~/Jensen_response/scripts/admixture_parameters_mode_fixedPopSize_MS_trajectoryTest.py $rho_in $adaptive_theta $selection $locusLength $NAm $EUr > ${file}_var
#python ~/Jensen_response/scripts/admixture_parameters_mode_fixedPopSize_MS.py $rho_in $adaptive_theta $selection $locusLength $NAm $EUr > ${file}_var

command=`cat ${file}_var | head -1` 

eval $command

segsites_line=`cat $file | grep -n segsites | cut -f1 -d':'`

(( segsites_line_head = segsites_line - 1 ))
(( segsites_line_tail = segsites_line - 6 ))

cat $file | head -${segsites_line_head} | tail -${segsites_line_tail} > $outFile

rm ${file}
rm ${file}_var

done
done
done




