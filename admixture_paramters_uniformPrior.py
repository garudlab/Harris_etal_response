import numpy
import sys

rho_in=float(sys.argv[1])
adaptive_theta=float(sys.argv[2])
selection=sys.argv[3] # boolean indicating if we are similating with selection
locusLength=int(sys.argv[4])

# input parameters:
#rho_in=10**-8

numberOfSimulations=1

for i in range(0, numberOfSimulations):

    # fixed parameters:
    totalSampleSize=145
    sampleSizeAfrica=0
    sampleSizeEurope=0
    sampleSizeAmerica=145
    smu=0 # test for 0, 0.01, and 10 for neutrality, HS, and SS respectively
    mu=10**-9 
    
    # First set the population size of Africa
    Nac=numpy.random.uniform(2.40*10**6, 9.13*10**6)

    # the above Nac is used for all subsequent parameter estimations

    #theta:
    theta=4*Nac*mu*locusLength

    #recombination rate:
    recombinationRate=4*Nac*locusLength*rho_in

    # demographic parameters:

    # scaledNeEurope
    #[0.39*10^6, 9.55*10^6] unscaled
    lower=0.39*10**6/Nac
    upper=9.55*10**6/Nac
    scaledNeEurope=numpy.random.uniform(lower, upper) 

    rand_Nea=numpy.random.uniform(3.58,4.83)
    scaledNeEurope_ancestral=(10**rand_Nea)/Nac

    rand_Tae=numpy.random.uniform(4.69,5.86)
    Tae=(10**rand_Tae) #do not divide by 10 because we want everything in terms of generations 
    scaledTimeSplitAfricaEurope=10**rand_Tae/(4*Nac)
    
    growthRateEurope=10**(numpy.log10(scaledNeEurope_ancestral/scaledNeEurope)/Tae) -1
    
    lower=1.11*10**6/Nac
    upper=28.8*10**6/Nac
    scaledNeAmerica=numpy.random.uniform(lower, upper)
    
    Tadm_rand=numpy.random.uniform(2.08, 3.82)
    Tadm=10**Tadm_rand # this is in terms of generations
    scaledTimeAdmixture=(10**Tadm_rand)/(4*Nac) 
    scaledTimeAdmixture2=scaledTimeAdmixture+0.0001 # from the old code -- need to add a little extra time?
    
    proportionAdmixture=numpy.random.uniform(0.64,0.97)
    
    rand_Naa=numpy.random.uniform(2.20,4.79)
    scaledNeAmerica_ancestral=10**rand_Naa/Nac
    growthRateAmerica=10**(numpy.log10(scaledNeAmerica_ancestral/scaledNeAmerica)/Tadm) -1

    rand_Ta=numpy.random.uniform(0.82*10**5, 3.45*10**5)
    Ta=rand_Ta*10 # convert to generations
    scaledTimeAfricaExpansion=(Ta- 1000)/(4*Nac) # T_a=237227 years and 1000 comes from bottleneck duration of 1000 gens 

    sev_rand=numpy.random.uniform(-0.15, 0.57) # units of log 10
    scaledNeAfricaBottleneck=(1000/10**sev_rand)/Nac # sev is in terms of duration/pop size, and duration is fixed at 1000 gen.

    Ta_rand=numpy.random.uniform(0.82*10**5, 3.45*10**5)
    Ta=Ta_rand*10 #convert to generations
    scaledTimeAfricaCrash=Ta/(4*Nac)

    lower=1.98*10**6/Nac
    upper=9.55*10**6/Nac
    scaledNeAfricaAncestral=numpy.random.uniform(lower, upper)

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
