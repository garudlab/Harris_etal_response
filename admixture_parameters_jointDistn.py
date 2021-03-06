import numpy
import sys
from scipy.stats import rv_discrete
import os.path

rho_in=float(sys.argv[1])
adaptive_theta=float(sys.argv[2])
selection=sys.argv[3] # boolean indicating if we are similating with selection
locusLength=int(sys.argv[4])

# input parameters:
#rho_in=10**-8

# read in the joint distribution of parameters:
distn=open(os.path.expanduser("~/Jensen_response/scripts/posterior_distributions_admixture/jointSamples_2.txt") , 'r')
distn.readline() #header
parameters=[]
for line in distn:
    line=line.strip().split()
    parameters.append(line)


numberOfSimulations=1

for i in range(0, numberOfSimulations):

    # randomly choose a set of parameters to simulate with
    index=numpy.random.choice(range(0,len(parameters)))
    parameter_set=parameters[index]
    
    Nac = float(parameter_set[0])
    rand_Tae = float(parameter_set[1])
    Tadm_rand = float(parameter_set[2])
    Ta = float(parameter_set[3])
    sev = float(parameter_set[4])
    Naa = float(parameter_set[5])
    Nec = float(parameter_set[6])
    Nea_log = float(parameter_set[7])
    Nnc = float(parameter_set[8])
    Nna = float(parameter_set[9])
    proportionAdmixture = float(parameter_set[10])
       
    # fixed parameters:
    totalSampleSize=145
    sampleSizeAfrica=0
    sampleSizeEurope=0
    sampleSizeAmerica=145
    smu=0 # test for 0, 0.01, and 10 for neutrality, HS, and SS respectively
    mu=10**-9 
    
    #theta:
    theta=4*Nac*mu*locusLength

    #recombination rate:
    recombinationRate=4*Nac*locusLength*rho_in

    # demographic parameters:

    # scaledNeEurope
    scaledNeEurope=Nec/Nac

    # Europe ancestral
    scaledNeEurope_ancestral=(10**Nea_log)/Nac
    
    # Time of Africa-Europe split
    Tae=(10**rand_Tae) #do not divide by 10 because we want everything in terms of generations 
    scaledTimeSplitAfricaEurope=10**rand_Tae/(4*Nac)
    
    growthRateEurope=10**(numpy.log10(scaledNeEurope_ancestral/scaledNeEurope)/Tae) -1

    # current population size of america
    scaledNeAmerica=Nnc/Nac

    # admixture time
    Tadm=10**Tadm_rand # this is in terms of generations
    scaledTimeAdmixture=(10**Tadm_rand)/(4*Nac) 
    scaledTimeAdmixture2=scaledTimeAdmixture+0.0001 # from the old code -- need to add a little extra time?
    
    
    # ancestral North America
    scaledNeAmerica_ancestral=10**Nna/Nac

    # growth rate america
    growthRateAmerica=10**(numpy.log10(scaledNeAmerica_ancestral/scaledNeAmerica)/Tadm) -1

    # time of Africa bottleneck
    scaledTimeAfricaExpansion=(Ta+ 1000)/(4*Nac) # T_a=237227 years and 1000 comes from bottleneck duration of 1000 gens 
    scaledTimeAfricaCrash=Ta/(4*Nac)

    scaledNeAfricaBottleneck=(1000/10**sev)/Nac # sev is in terms of duration/pop size, and duration is fixed at 1000 gen.
    
    # ancestral african population size
    scaledNeAfricaAncestral=Naa/Nac

    # want adaptive mutation to be 1/Ne in population:
    starting_frequency=1/(2*scaledNeAmerica*Nac)

    # other hyper parameters:
    s_rand=numpy.random.uniform(0,1)
    age_rand=numpy.random.uniform(0,1)    

    s=2*Nac*s_rand
    s2=2*s
    age=age_rand*scaledTimeAdmixture

    # print these parameters to a file 

    if selection == 'True':
        command='~/./msms/bin/msms ' + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture2) + ' 3 2 -ej ' + str(scaledTimeAdmixture2) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength) +' -SAA ' + str(s2) +' -SaA ' + str(s) +' -SI ' + str(age) +' 3 0 0 ' + str(starting_frequency) + ' -SFC -Smu ' + str(adaptive_theta) + ' -Sp 0.5 -oTrace  > $file\n'

    else:
        command='~/./msms/bin/msms ' + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture2) + ' 3 2 -ej ' + str(scaledTimeAdmixture2) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' > $file\n'

    print command + '\n' + str(Nac) + '\n' + str(s) + '\n' + str(age)
