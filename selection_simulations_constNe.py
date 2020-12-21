import numpy
import sys
from scipy.stats import rv_discrete
import os.path

rho_in=float(sys.argv[1])
adaptive_theta=float(sys.argv[2])
selection=sys.argv[3] # boolean indicating if we are similating with selection
locusLength=int(sys.argv[4])
s_rand=float(sys.argv[5]) # NEW
age_rand=float(sys.argv[6]) # NEW
partial_frequency=float(sys.argv[7]) # NEW  
 
# input parameters:
#rho_in=10**-8

# this code generates the MS commands for the admixture models inferred in this paper, specifically for the fixed population size models that were inferred to fit the data. 

numberOfSimulations=1

for i in range(0, numberOfSimulations):
    

    # fixed parameters:
    Ne=2657111
    totalSampleSize=145
    smu=0 # test for 0, 0.01, and 10 for neutrality, HS, and SS respectively
    mu=10**-9 
    
    #theta:
    theta=4*Ne*mu*locusLength

    #recombination rate:
    recombinationRate=4*Ne*locusLength*rho_in


    # other hyper parameters:
    #s_rand=numpy.random.uniform(0,1)
    #age_rand=numpy.random.uniform(0,1)    

    s=2*Ne*s_rand
    s2=2*s
    age=age_rand


    # print these parameters to a file 


    if selection=='False':
        command='~/./msdir/ms ' + str(totalSampleSize) + ' 1 -t ' + str(theta) + ' -r ' + str(recombinationRate) + ' ' + str(locusLength) + ' -seed $RANDOM $RANDOM $RANDOM > $file\n'
    elif selection=='True':
        #command='~/./software/msms/bin/msms ' + str(totalSampleSize) + ' 1 -N ' + str(Ne) + ' -t ' + str(theta) + ' -r ' + str(recombinationRate) + ' ' + str(locusLength) + ' -SAA'  + ' ' + str(s2) +' -SaA ' + str(s) + ' -SF ' + str(age) + ' ' + str(partial_frequency) +  ' -SFC -Smu ' + str(adaptive_theta) + ' -Sp 0.5 -oTrace  > $file\n'

        #command='~/./software/msms/bin/msms ' + str(totalSampleSize) + ' 1 -N ' + str(Ne) + ' -t ' + str(theta) + ' -r ' + str(recombinationRate) + ' ' + str(locusLength) + ' -SAA'  + ' ' + str(s2) +' -SaA ' + str(s) + ' -SF ' + str(age) + ' 1 ' + str(partial_frequency)  + ' -Smu ' +  str(adaptive_theta) + ' -Sp 0.5 -oTrace -SFC  > $file\n'

        #command='~/./software/msms/bin/msms ' + str(totalSampleSize) + ' 1 -N ' + str(Ne) + ' -t ' + str(theta) + ' -r ' + str(recombinationRate) + ' ' + str(locusLength) + ' -SAA'  + ' ' + str(s2) +' -SaA ' + str(s)  + ' -Smu ' +  str(adaptive_theta) + ' -Sp 0.5 -oTrace -SFC  > $file\n'
        command='~/./software/msms/bin/msms ' + str(totalSampleSize) + ' 1 -N ' + str(Ne) + ' -t ' + str(theta) + ' -r ' + str(recombinationRate) + ' ' + str(locusLength) + ' -SAA'  + ' ' + str(s2) +' -SaA ' + str(s)  + ' -Smu ' +  str(adaptive_theta) + ' -Sp 0.5 -SI 0.0001 1 1.881743e-07 -oTrace -SFC  > $file\n'

    print command 
