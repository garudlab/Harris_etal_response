import numpy
import sys
from scipy.stats import rv_discrete
import os.path

rho_in=float(sys.argv[1])
adaptive_theta=float(sys.argv[2])
selection=sys.argv[3] # boolean indicating if we are similating with selection
locusLength=int(sys.argv[4])
mu=float(sys.argv[5])

# input parameters:
#rho_in=10**-8


def generate_random_number(inFile, parameter):

    # record the 95% CI for each parameter
    param_95CI={}
    param_95CI['Nac']=[2.40*10**6, 9.13*10**6]
    param_95CI['Nec']=[0.39*10**6,9.55*10**6]
    param_95CI['Nea']=[3.58,4.83]
    param_95CI['rand_Tae']=[4.69,5.86]
    param_95CI['Nnc']=[1.11*10**6,28.8*10**6]
    param_95CI['Tadm']=[2.08, 3.82]
    param_95CI['proportionAdmixture']=[0.64,0.97]
    param_95CI['Nna']=[2.20,4.79]
    param_95CI['Ta']=[0.82*10**6, 3.45*10**6] # multiplied by 10 to make in terms of generations
    param_95CI['sev']=[-0.15, 0.57]
    param_95CI['Naa']=[1.98*10**6,9.55*10**6]

    # read in the posterior distributions for each parameter
    distn = open(os.path.expanduser("~/Jensen_response/scripts/posterior_distributions_admixture/%s" %inFile) , 'r')
    distn.readline() #header
    param_val=[]
    density=[]
    for line in distn:
        line=line.strip().split()
        param_val.append(float(line[0]))
        density.append(float(line[1]))
    density=numpy.asarray(density)
    
    # make density sum to 1
    density=density/density.sum()

    #value_found=False

    # draw random value
    #while value_found == False:
    random_value=numpy.random.choice(param_val, p=density, replace=True, size=1)[0]
    # check if the random value is within the 95CI
    
    return(random_value)
    


numberOfSimulations=1

for i in range(0, numberOfSimulations):

    # fixed parameters:
    totalSampleSize=37
    sampleSizeAfrica=0
    sampleSizeEurope=0
    sampleSizeAmerica=37
    smu=0 # test for 0, 0.01, and 10 for neutrality, HS, and SS respectively
    mu=10**-9 
    
    # First set the current population size of Africa
    Nac=generate_random_number('distrNeAfr.txt','Nac')
    # the above Nac is used for scaling subsequent parameter estimations

    #theta:
    theta=3*Nac*mu*locusLength

    #recombination rate:
    recombinationRate=3*Nac*locusLength*rho_in

    # demographic parameters:

    # scaledNeEurope
    Nec=generate_random_number('distrNeEur.txt','Nec')
    scaledNeEurope=Nec/Nac

    # Europe ancestral
    Nea_log=generate_random_number('distrLogNeEurBn.txt','Nea')
    scaledNeEurope_ancestral=(10**Nea_log)/Nac
    
    # Time of Africa-Europe split
    rand_Tae=generate_random_number('distrLogTimeSplitAfrEur.txt', 'rand_Tae')  
    Tae=(10**rand_Tae) #do not divide by 10 because we want everything in terms of generations 
    scaledTimeSplitAfricaEurope=10**rand_Tae/(3*Nac)
    
    growthRateEurope=-(1 / scaledTimeSplitAfricaEurope) * numpy.log(scaledNeEurope_ancestral / scaledNeEurope)
    
    #growthRateEurope=10**(numpy.log10(scaledNeEurope_ancestral/scaledNeEurope)/Tae) -1

    # current population size of america
    Nnc=generate_random_number('distrNeAme.txt','Nnc')
    scaledNeAmerica=Nnc/Nac

    # admixture time
    Tadm_rand=generate_random_number('distrLogTimeAdm.txt','Tadm')
    Tadm=10**Tadm_rand # this is in terms of generations
    scaledTimeAdmixture=(10**Tadm_rand)/(3*Nac) 
    scaledTimeAdmixture2=scaledTimeAdmixture+0.0001 # from the old code -- need to add a little extra time?
    
    proportionAdmixture=generate_random_number('distrPropAdm.txt','proportionAdmixture')
    
    # ancestral North America
    Nna=generate_random_number('distrLogNeAmeBn.txt','Nna')
    scaledNeAmerica_ancestral=10**Nna/Nac

    # growth rate america
    growthRateAmerica=-(1 / scaledTimeAdmixture) * numpy.log(scaledNeAmerica_ancestral / scaledNeAmerica)
    
    #growthRateAmerica=10**(numpy.log10(scaledNeAmerica_ancestral/scaledNeAmerica)/Tadm) -1

    # time of Africa bottleneck
    Ta=generate_random_number('distrTimeAfrBn.txt','Ta') # this is already in terms of generations (it has been multiplied by 10)
    scaledTimeAfricaExpansion=(Ta- 1000)/(3*Nac) # T_a=237227 years and 1000 comes from bottleneck duration of 1000 gens 
    scaledTimeAfricaCrash=Ta/(3*Nac)

    sev=generate_random_number('distrSeverity.txt','sev')
    scaledNeAfricaBottleneck=(1000/10**sev)/Nac # sev is in terms of duration/pop size, and duration is fixed at 1000 gen.

    lower=1.98*10**6/Nac
    upper=9.55*10**6/Nac
    
    # ancestral african population size
    Naa=generate_random_number('distrNeAfrAnc.txt','Naa')
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
        command='~/./msms/bin/msms ' + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength) +' -SAA ' + str(s2) +' -SaA ' + str(s) +' -SI ' + str(age) +' 3 0 0 ' + str(starting_frequency) + ' -SFC -Smu ' + str(adaptive_theta) + ' -Sp 0.5 -oTrace  > $file\n'

    else:
        command='~/./msms/bin/msms ' + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' > $file\n'

    print command + '\n' + str(Nac) + '\n' + str(s) + '\n' + str(age)
