import numpy
import sys
from scipy.stats import rv_discrete
import os.path

model=str(sys.argv[1])
global rho_in
rho_in=float(sys.argv[2])
global adaptive_theta
adaptive_theta=float(sys.argv[3])
global selection
selection=sys.argv[4] # boolean indicating if we are similating with selection
global locusLength
locusLength=int(sys.argv[5])
global MS_flag
MS_flag=sys.argv[6] # boolean indicating if we want an ms vs msms command
global totalSampleSize
totalSampleSize=145
global mu
mu=10**-9
global msdir
msdir='~/./software/msdir/ms  '
global msmsdir
msmsdir='~/./software/msms/bin/msms  '

# possible models: admixture_posterior, admixture_mode, admixture_bot_mode, constNe10e6, constNe2.2e6, dadi1, dadi2



###########################


def constNe10e6():
    
    Ne=1000000
    theta=4*Ne*mu*locusLength
    recombinationRate=4*Ne*locusLength*rho_in

    if MS_flag=='True':
        command=msdir + str(totalSampleSize) + ' 1 -t ' + str(theta) + ' -r ' + str(recombinationRate) + ' ' + str(locusLength) + ' -seed $RANDOM $RANDOM $RANDOM > $file\n'
    else:
        command=msmsdir + str(totalSampleSize) + ' 1 -N ' + str(Ne) + ' -t ' + str(theta) + ' -r ' + str(recombinationRate) + ' ' + str(locusLength) + ' > $file\n'

    print command 


##########################

def constNe27e6():

    Ne=2657111
    theta=4*Ne*mu*locusLength
    recombinationRate=4*Ne*locusLength*rho_in

    if MS_flag=='True':
        command=msdir + str(totalSampleSize) + ' 1 -t ' + str(theta) + ' -r ' + str(recombinationRate) + ' ' + str(locusLength) + ' -seed $RANDOM $RANDOM $RANDOM > $file\n'
    else:
        command=msmsdir  + str(totalSampleSize) + ' 1 -N ' + str(Ne) + ' -t ' + str(theta) + ' -r ' + str(recombinationRate) + ' ' + str(locusLength) + ' > $file\n'

    print command 

#########################
def dadi1():
    # severe short bottleneck model
    Ne=2317520
    theta=4*Ne*mu*locusLength
    recombinationRate=4*Ne*locusLength*rho_in

    if MS_flag=='True':
        command=msdir + str(totalSampleSize) + ' 1 -t ' + str(theta) + ' -r ' + str(recombinationRate) + ' ' + str(locusLength) + ' -eN 0.2745424 0.003327787 -eN 0.2747088 1.663894 -seed $RANDOM $RANDOM $RANDOM > $file\n'
    else:
        command=msmsdir  + str(totalSampleSize) + ' 1 -N ' + str(Ne) + ' -t ' + str(theta) + ' -r ' + str(recombinationRate) + ' ' + str(locusLength) + ' -eN 0.2745424 0.003327787 -eN 0.2747088 1.663894 > $file\n'


    print command 

#########################
def dadi2():
    # shallow long bottleneck model. 
    Ne=2317075
    theta=4*Ne*mu*locusLength
    recombinationRate=4*Ne*locusLength*rho_in

    if MS_flag=='True':
        command=msdir + str(totalSampleSize) + ' 1 -t ' + str(theta) + ' -r ' + str(recombinationRate) + ' ' + str(locusLength) + ' -eN 0.1178203 0.5891016 -eN 0.1590574 1.472754 -seed $RANDOM $RANDOM $RANDOM > $file\n'
    else:
        command=msmsdir + str(totalSampleSize) + ' 1 -N ' + str(Ne) + ' -t ' + str(theta) + ' -r ' + str(recombinationRate) + ' ' + str(locusLength) +' -eN 0.1178203 0.5891016 -eN 0.1590574 1.472754  > $file\n'


    print command 


################

def admixture_mode_Garud2015():
    # model implemented in Garud et al. 2015
    Nac = 4975360.0
    rand_Tae = 5.29
    Tadm_rand = 3.16
    Ta = 237227.0*10 #generations
    sev = 0.21
    Naa = 5224100.0
    Nec = 3122470.0
    Nea_log = 4.23
    Nnc = 15984500.0 
    Nna = 3.4 # need to unlog
    proportionAdmixture = 0.85

    # fixed parameters:
    sampleSizeAfrica=0
    sampleSizeEurope=0
    sampleSizeAmerica=totalSampleSize
    smu=adaptive_theta # test for 0, 0.01, and 10 for neutrality, HS, and SS respectively
    
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
    
    #growthRateEurope=-(1 / scaledTimeSplitAfricaEurope) * numpy.log(scaledNeEurope_ancestral / scaledNeEurope)
    growthRateEurope=-2.674174e-05

    #growthRateEurope=10**(numpy.log10(scaledNeEurope_ancestral/scaledNeEurope)/Tae) -1

    # current population size of america
    scaledNeAmerica=Nnc/Nac

    # admixture time
    Tadm=10**Tadm_rand # this is in terms of generations
    scaledTimeAdmixture=(10**Tadm_rand)/(4*Nac) 
    #scaledTimeAdmixture2=scaledTimeAdmixture+0.0001 # from the old code -- need to add a little extra time?
    
    
    # ancestral North America
    scaledNeAmerica_ancestral=10**Nna/Nac

    # growth rate america
    #growthRateAmerica=10**(numpy.log10(scaledNeAmerica_ancestral/scaledNeAmerica)/Tadm) -1
    #growthRateAmerica=-(1 / scaledTimeAdmixture) * numpy.log(scaledNeAmerica_ancestral / scaledNeAmerica)
    growthRateAmerica=-0.006059296
        
    # time of Africa bottleneck
    scaledTimeAfricaExpansion=(Ta- 1000)/(4*Nac) # T_a=237227 years and 1000 comes from bottleneck duration of 1000 gens 
    scaledTimeAfricaCrash=Ta/(4*Nac)
    
    scaledNeAfricaBottleneck=(1000/10**sev)/Nac # sev is in terms of duration/pop size, and duration is fixed at 1000 gen.
    
    # ancestral african population size
    scaledNeAfricaAncestral=Naa/Nac
    #print scaledNeAfricaAncestral
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
        command=msmsdir + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength) +' -SAA ' + str(s2) +' -SaA ' + str(s) +' -SI ' + str(age) +' 3 0 0 ' + str(starting_frequency) + ' -SFC -Smu ' + str(adaptive_theta) + ' -Sp 0.5 -oTrace  > $file\n'

    else:
        if MS_flag=='True':
            command=msdir + str(totalSampleSize) + ' 1  -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' -seed $RANDOM $RANDOM $RANDOM > $file\n'
        else:
            command=msmsdir + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' > $file\n'
            
    print command





##########################

def admixture_mode():
    #Duchen et al. model, implemented w/ the modal parameters
    Nac = 4975360.0
    rand_Tae = 5.29
    Tadm_rand = 3.16
    Ta = 237227.0*10 #generations
    sev = 0.21
    Naa = 5224100.0
    Nec = 3122470.0
    Nea_log = 4.23
    Nnc = 15984500.0 
    Nna = 3.4 # need to unlog
    proportionAdmixture = 0.85

    # fixed parameters:
    sampleSizeAfrica=0
    sampleSizeEurope=0
    sampleSizeAmerica=totalSampleSize
    smu=adaptive_theta # test for 0, 0.01, and 10 for neutrality, HS, and SS respectively
    
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
    
    growthRateEurope=-(1 / scaledTimeSplitAfricaEurope) * numpy.log(scaledNeEurope_ancestral / scaledNeEurope)
    
    #growthRateEurope=10**(numpy.log10(scaledNeEurope_ancestral/scaledNeEurope)/Tae) -1

    # current population size of america
    scaledNeAmerica=Nnc/Nac

    # admixture time
    Tadm=10**Tadm_rand # this is in terms of generations
    scaledTimeAdmixture=(10**Tadm_rand)/(4*Nac) 
    #scaledTimeAdmixture2=scaledTimeAdmixture+0.0001 # from the old code -- need to add a little extra time?
    
    
    # ancestral North America
    scaledNeAmerica_ancestral=10**Nna/Nac

    # growth rate america
    #growthRateAmerica=10**(numpy.log10(scaledNeAmerica_ancestral/scaledNeAmerica)/Tadm) -1
    growthRateAmerica=-(1 / scaledTimeAdmixture) * numpy.log(scaledNeAmerica_ancestral / scaledNeAmerica)
        
    # time of Africa bottleneck
    scaledTimeAfricaExpansion=(Ta- 1000)/(4*Nac) # T_a=237227 years and 1000 comes from bottleneck duration of 1000 gens 
    scaledTimeAfricaCrash=Ta/(4*Nac)
    
    scaledNeAfricaBottleneck=(1000/10**sev)/Nac # sev is in terms of duration/pop size, and duration is fixed at 1000 gen.
    
    # ancestral african population size
    scaledNeAfricaAncestral=Naa/Nac
    #print scaledNeAfricaAncestral
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
        command=msmsdir  + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength) +' -SAA ' + str(s2) +' -SaA ' + str(s) +' -SI ' + str(age) +' 3 0 0 ' + str(starting_frequency) + ' -SFC -Smu ' + str(adaptive_theta) + ' -Sp 0.5 -oTrace  > $file\n'

    else:
        if MS_flag=='True':
            command=msdir + str(totalSampleSize) + ' 1  -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' -seed $RANDOM $RANDOM $RANDOM > $file\n'
        else:
            command=msmsdir + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' > $file\n'
            
    print command




##########################

def admixture_mode_Harris():
    #admixture model implemented in harris et al. with the modal parameters
    Nac = 4975360.0
    rand_Tae = 5.29
    Tadm_rand = 3.16
    Ta = 237227.0*10 #generations
    sev = 0.21
    Naa = 5224100.0
    Nec = 3122470.0
    Nea_log = 4.23
    Nnc = 15984500.0 
    Nna = 3.4 # need to unlog
    proportionAdmixture = 0.85

    # fixed parameters:
    sampleSizeAfrica=0
    sampleSizeEurope=0
    sampleSizeAmerica=totalSampleSize
    smu=adaptive_theta # test for 0, 0.01, and 10 for neutrality, HS, and SS respectively
    
    #theta:
    theta=4*Naa*mu*locusLength

    #recombination rate:
    recombinationRate=4*Naa*locusLength*rho_in

    # demographic parameters:

    # scaledNeEurope
    scaledNeEurope=Nec/Naa
    
    # Europe ancestral
    scaledNeEurope_ancestral=(10**Nea_log)/(4*Naa)

    # Time of Africa-Europe split
    Tae=(10**rand_Tae) #do not divide by 10 because we want everything in terms of generations 
    scaledTimeSplitAfricaEurope=10**rand_Tae/(4*Naa)
    
    growthRateEurope=-(1 / scaledTimeSplitAfricaEurope) * numpy.log(scaledNeEurope_ancestral / scaledNeEurope)
    
    #growthRateEurope=10**(numpy.log10(scaledNeEurope_ancestral/scaledNeEurope)/Tae) -1

    # current population size of america
    scaledNeAmerica=Nnc/Naa

    # admixture time
    Tadm=10**Tadm_rand # this is in terms of generations
    scaledTimeAdmixture=(10**Tadm_rand)/(4*Naa) 
    #scaledTimeAdmixture2=scaledTimeAdmixture+0.0001 # from the old code -- need to add a little extra time?
    
    
    # ancestral North America
    scaledNeAmerica_ancestral=10**Nna/(4*Naa)

    # growth rate america
    #growthRateAmerica=10**(numpy.log10(scaledNeAmerica_ancestral/scaledNeAmerica)/Tadm) -1
    growthRateAmerica=-(1 / scaledTimeAdmixture) * numpy.log(scaledNeAmerica_ancestral / scaledNeAmerica)
        
    # time of Africa bottleneck
    scaledTimeAfricaExpansion=(Ta- 1000)/(4*Naa) # T_a=237227 years and 1000 comes from bottleneck duration of 1000 gens 
    scaledTimeAfricaCrash=Ta/(4*Naa)
    
    scaledNeAfricaBottleneck=(1000/10**sev)/Naa # sev is in terms of duration/pop size, and duration is fixed at 1000 gen.
    
    # ancestral african population size
    scaledNeAfricaAncestral=Naa/Naa
    #print scaledNeAfricaAncestral
    # want adaptive mutation to be 1/Ne in population:
    starting_frequency=1/(2*scaledNeAmerica*Naa)

    # other hyper parameters:
    s_rand=numpy.random.uniform(0,1)
    age_rand=numpy.random.uniform(0,1)    

    s=2*Naa*s_rand
    s2=2*s
    age=age_rand*scaledTimeAdmixture

    # print these parameters to a file 

    if selection == 'True':
        command=msmsdir  + str(totalSampleSize) + ' 1 -N ' + str(Naa) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength) +' -SAA ' + str(s2) +' -SaA ' + str(s) +' -SI ' + str(age) +' 3 0 0 ' + str(starting_frequency) + ' -SFC -Smu ' + str(adaptive_theta) + ' -Sp 0.5 -oTrace -seed $RANDOM $RANDOM $RANDOM  > $file\n'

    else:
        if MS_flag=='True':
            command=msdir + str(totalSampleSize) + ' 1  -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' -seed $RANDOM $RANDOM $RANDOM > $file\n'
        else:
            command=msmsdir  + str(totalSampleSize) + ' 1 -N ' + str(Naa) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' > $file\n'
            
    print command


##################

def admixture_bot_mode():
    #Duchen et al. model of admixutre with a bottleneck
    Nac=3100520
    theta=4*Nac*mu*locusLength
    recombinationRate=4*Nac*locusLength*rho_in


    numberOfSimulations=1
    sampleSizeAfrica=0
    sampleSizeEurope=0
    sampleSizeAmerica=totalSampleSize
    scaledNeEurope=0.7318321
    growthRateEurope=850.52
    scaledNeAmerica=2.968357
    scaledTimeAdmixture=3.757037**-05
    proportionAdmixture=0.871794
    scaledNeEuropeAdmixtureBn=0.004446651
    scaledTimeSplitAfricaEurope=0.006037894
    scaledTimeAfricaExpansion=0.03233073
    scaledNeAfricaBottleneck=0.1599473
    scaledTimeAfricaCrash=0.03241136
    scaledNeAfricaAncestral=1.0401

    if MS_flag=='True':
        command=msdir + str(totalSampleSize) +' ' + str(numberOfSimulations) + ' -t ' + str(theta) + ' -I 3 ' + str(sampleSizeAfrica) + ' ' + str(sampleSizeEurope) + ' ' + str(sampleSizeAmerica) + ' -n 2 ' + str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +  ' -en ' + str(scaledTimeAdmixture) + ' 3 ' + str(scaledNeEuropeAdmixtureBn) +' -ej ' + str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) + ' 2 1 -en ' +str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' + str(recombinationRate) +' ' + str(locusLength) +  ' -seed $RANDOM $RANDOM $RANDOM  > $file\n'

    else:
        command=msmsdir + str(totalSampleSize) +' ' + str(numberOfSimulations) + ' -N ' + str(Nac) + ' -t ' + str(theta) + ' -I 3 ' + str(sampleSizeAfrica) + ' ' + str(sampleSizeEurope) + ' ' + str(sampleSizeAmerica) + ' -n 2 ' + str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +  ' -en ' + str(scaledTimeAdmixture) + ' 3 ' + str(scaledNeEuropeAdmixtureBn) +' -ej ' + str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) + ' 2 1 -en ' +str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' + str(recombinationRate) +' ' + str(locusLength) + ' > $file\n'


    print command
###################


def admixture_bot_mode_Garud2015():
    # admixture model implemented in Garud et al. 2015 wtih a bottleneck
    Nac=3100520
    theta=4*Nac*mu*locusLength
    recombinationRate=4*Nac*locusLength*rho_in


    numberOfSimulations=1
    sampleSizeAfrica=0
    sampleSizeEurope=0
    sampleSizeAmerica=totalSampleSize
    scaledNeEurope=0.7318321
    growthRateEurope=-6.857889e-05 # old code mistake
    scaledNeAmerica=2.968357
    scaledTimeAdmixture=3.757037**-05
    proportionAdmixture=0.871794
    scaledNeEuropeAdmixtureBn=0.004446651
    scaledTimeSplitAfricaEurope=0.006037894
    scaledTimeAfricaExpansion=0.03233073
    scaledNeAfricaBottleneck=0.1599473
    scaledTimeAfricaCrash=0.03241136
    scaledNeAfricaAncestral=1.0401

    if MS_flag=='True':
        command=msdir + str(totalSampleSize) +' ' + str(numberOfSimulations) + ' -t ' + str(theta) + ' -I 3 ' + str(sampleSizeAfrica) + ' ' + str(sampleSizeEurope) + ' ' + str(sampleSizeAmerica) + ' -n 2 ' + str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +  ' -en ' + str(scaledTimeAdmixture) + ' 3 ' + str(scaledNeEuropeAdmixtureBn) +' -ej ' + str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) + ' 2 1 -en ' +str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' + str(recombinationRate) +' ' + str(locusLength) +  ' -seed $RANDOM $RANDOM $RANDOM  > $file\n'

    else:
        command=msmsdir  + str(totalSampleSize) +' ' + str(numberOfSimulations) + ' -N ' + str(Nac) + ' -t ' + str(theta) + ' -I 3 ' + str(sampleSizeAfrica) + ' ' + str(sampleSizeEurope) + ' ' + str(sampleSizeAmerica) + ' -n 2 ' + str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +  ' -en ' + str(scaledTimeAdmixture) + ' 3 ' + str(scaledNeEuropeAdmixtureBn) +' -ej ' + str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) + ' 2 1 -en ' +str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' + str(recombinationRate) +' ' + str(locusLength) + ' > $file\n'


    print command
###################

def generate_random_number(inFile, parameter):
    # this function generates parameter values from the 95CI defined in Duchen et al. 

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

    value_found=False

    # draw random value
    while value_found == False:
        random_value=numpy.random.choice(param_val, p=density, replace=True, size=1)[0]
        # check if the random value is within the 95CI
        if random_value >= param_95CI[parameter][0] and random_value <= param_95CI[parameter][1]:
            value_found=True

    return(random_value)
    




#####################
def admixture_posterior():
    # Duchen et al. admixture model, drawing from the posterior. 

    # fixed parameters:
    sampleSizeAfrica=0
    sampleSizeEurope=0
    sampleSizeAmerica=totalSampleSize
    smu=adaptive_theta 
    
    # First set the current population size of Africa
    Nac=generate_random_number('distrNeAfr.txt','Nac')
    # the above Nac is used for scaling subsequent parameter estimations

    theta=4*Nac*mu*locusLength
    recombinationRate=4*Nac*locusLength*rho_in

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
    scaledTimeSplitAfricaEurope=10**rand_Tae/(4*Nac)
    
    growthRateEurope=-(1 / scaledTimeSplitAfricaEurope) * numpy.log(scaledNeEurope_ancestral / scaledNeEurope)
    
    #growthRateEurope=10**(numpy.log10(scaledNeEurope_ancestral/scaledNeEurope)/Tae) -1

    # current population size of america
    Nnc=generate_random_number('distrNeAme.txt','Nnc')
    scaledNeAmerica=Nnc/Nac

    # admixture time
    Tadm_rand=generate_random_number('distrLogTimeAdm.txt','Tadm')
    Tadm=10**Tadm_rand # this is in terms of generations
    scaledTimeAdmixture=(10**Tadm_rand)/(4*Nac) 
    
    proportionAdmixture=generate_random_number('distrPropAdm.txt','proportionAdmixture')
    
    # ancestral North America
    Nna=generate_random_number('distrLogNeAmeBn.txt','Nna')
    scaledNeAmerica_ancestral=10**Nna/Nac

    # growth rate america
    growthRateAmerica=-(1 / scaledTimeAdmixture) * numpy.log(scaledNeAmerica_ancestral / scaledNeAmerica)
    
    #growthRateAmerica=10**(numpy.log10(scaledNeAmerica_ancestral/scaledNeAmerica)/Tadm) -1

    # time of Africa bottleneck
    Ta=generate_random_number('distrTimeAfrBn.txt','Ta') # this is already in terms of generations (it has been multiplied by 10)
    scaledTimeAfricaExpansion=(Ta- 1000)/(4*Nac) # T_a=237227 years and 1000 comes from bottleneck duration of 1000 gens 
    scaledTimeAfricaCrash=Ta/(4*Nac)

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
        command=msmsdir  + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength) +' -SAA ' + str(s2) +' -SaA ' + str(s) +' -SI ' + str(age) +' 3 0 0 ' + str(starting_frequency) + ' -SFC -Smu ' + str(adaptive_theta) + ' -Sp 0.5 -oTrace  > $file\n'

    else:
        if MS_flag=='True':
            command=msdir + str(totalSampleSize) + ' 1  -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength) + ' -seed $RANDOM $RANDOM $RANDOM   > $file\n'
            
        else:
            command=msmsdir + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' > $file\n'

    print command



#####################
def admixture_posterior_Harris():
    # Harris et al. implementation of admixture model, drawing from the posterior. 

    # fixed parameters:
    sampleSizeAfrica=0
    sampleSizeEurope=0
    sampleSizeAmerica=totalSampleSize
    smu=adaptive_theta 
    

    # First, set the ANCESTRAL population size of arica
    # ancestral african population size
    Naa=generate_random_number('distrNeAfrAnc.txt','Naa')
    scaledNeAfricaAncestral=Naa/Naa


    # the above Naa is used for scaling subsequent parameter estimations

    theta=4*Naa*mu*locusLength
    recombinationRate=4*Naa*locusLength*rho_in

    # demographic parameters:

    # scaledNeEurope
    Nec=generate_random_number('distrNeEur.txt','Nec')
    scaledNeEurope=Nec/Naa

    # Europe ancestral
    Nea_log=generate_random_number('distrLogNeEurBn.txt','Nea')
    scaledNeEurope_ancestral=(10**Nea_log)/(4*Naa) # Note that Europe ancestral in Harris et al. was scaled by 4Ne rather than Ne
    
    # Time of Africa-Europe split
    rand_Tae=generate_random_number('distrLogTimeSplitAfrEur.txt', 'rand_Tae')  
    Tae=(10**rand_Tae) #do not divide by 10 because we want everything in terms of generations 
    scaledTimeSplitAfricaEurope=10**rand_Tae/(4*Naa)
    
    growthRateEurope=-(1 / scaledTimeSplitAfricaEurope) * numpy.log(scaledNeEurope_ancestral / scaledNeEurope)
    
    #growthRateEurope=10**(numpy.log10(scaledNeEurope_ancestral/scaledNeEurope)/Tae) -1

    # current population size of america
    Nnc=generate_random_number('distrNeAme.txt','Nnc')
    scaledNeAmerica=Nnc/Naa

    # admixture time
    Tadm_rand=generate_random_number('distrLogTimeAdm.txt','Tadm')
    Tadm=10**Tadm_rand # this is in terms of generations
    scaledTimeAdmixture=(10**Tadm_rand)/(4*Naa) 
    
    proportionAdmixture=generate_random_number('distrPropAdm.txt','proportionAdmixture')
    
    # ancestral North America
    Nna=generate_random_number('distrLogNeAmeBn.txt','Nna')
    scaledNeAmerica_ancestral=10**Nna/(4*Naa) # Note that America ancestral in Harris et al. was scaled by 4Ne rather than Ne 

    # growth rate america
    growthRateAmerica=-(1 / scaledTimeAdmixture) * numpy.log(scaledNeAmerica_ancestral / scaledNeAmerica)
    
    #growthRateAmerica=10**(numpy.log10(scaledNeAmerica_ancestral/scaledNeAmerica)/Tadm) -1


    # Set the current population size of Africa
    Nac=generate_random_number('distrNeAfr.txt','Nac')
    # scale it
    Nac=Nac/Naa

    # time of Africa bottleneck
    Ta=generate_random_number('distrTimeAfrBn.txt','Ta') # this is already in terms of generations (it has been multiplied by 10)
    scaledTimeAfricaExpansion=(Ta- 1000)/(4*Naa) # T_a=237227 years and 1000 comes from bottleneck duration of 1000 gens 
    scaledTimeAfricaCrash=Ta/(4*Naa)

    sev=generate_random_number('distrSeverity.txt','sev')
    scaledNeAfricaBottleneck=(1000/10**sev)/Naa # sev is in terms of duration/pop size, and duration is fixed at 1000 gen.

    lower=1.98*10**6/Naa
    upper=9.55*10**6/Naa
    

    # want adaptive mutation to be 1/Ne in population:
    # note I am not simulating selection with the Harris et al. implementation. If I do so, I need to change the starting frequency to match thiers. 
    
    starting_frequency=1/(2*scaledNeAmerica*Naa) #CHANGE

    # other hyper parameters:
    s_rand=numpy.random.uniform(0,1) #CHANGE
    age_rand=numpy.random.uniform(0,1) #CHANGE

    s=2*Nac*s_rand #CHANGE
    s2=2*s
    age=age_rand*scaledTimeAdmixture

    # print these parameters to a file 

    if selection == 'True':
        command=msmsdir + str(totalSampleSize) + ' 1 -N ' + str(Naa) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength) +' -SAA ' + str(s2) +' -SaA ' + str(s) +' -SI ' + str(age) +' 3 0 0 ' + str(starting_frequency) + ' -SFC -Smu ' + str(adaptive_theta) + ' -Sp 0.5 -oTrace  > $file\n'

    else:
        if MS_flag=='True':
            command=msdir + str(totalSampleSize) + ' 1  -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' -seed $RANDOM $RANDOM $RANDOM > $file\n'
            
        else:
            command=msmsdir + str(totalSampleSize) + ' 1 -N ' + str(Naa) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' > $file\n'

    print command




#####################
def admixture_posterior_Harris_Nac():
    # Harris et al. implementation of admixture, drawing from posterior -- this did not go into the final paper, but was something I was testing re. scaling w/ Nac. 
    # fixed parameters:
    sampleSizeAfrica=0
    sampleSizeEurope=0
    sampleSizeAmerica=totalSampleSize
    smu=adaptive_theta 
    

    # First, set the ANCESTRAL population size of arica
    # ancestral african population size
    Nac=generate_random_number('distrNeAfr.txt','Nac')

    theta=4*Nac*mu*locusLength
    recombinationRate=4*Nac*locusLength*rho_in

    # demographic parameters:

    # scaledNeEurope
    Nec=generate_random_number('distrNeEur.txt','Nec')
    scaledNeEurope=Nec/Nac

    # Europe ancestral
    Nea_log=generate_random_number('distrLogNeEurBn.txt','Nea')
    scaledNeEurope_ancestral=(10**Nea_log)/(4*Nac) # Note that Europe ancestral in Harris et al. was scaled by 4Ne rather than Ne
    
    # Time of Africa-Europe split
    rand_Tae=generate_random_number('distrLogTimeSplitAfrEur.txt', 'rand_Tae')  
    Tae=(10**rand_Tae) #do not divide by 10 because we want everything in terms of generations 
    scaledTimeSplitAfricaEurope=10**rand_Tae/(4*Nac)
    
    growthRateEurope=-(1 / scaledTimeSplitAfricaEurope) * numpy.log(scaledNeEurope_ancestral / scaledNeEurope)
    
    #growthRateEurope=10**(numpy.log10(scaledNeEurope_ancestral/scaledNeEurope)/Tae) -1

    # current population size of america
    Nnc=generate_random_number('distrNeAme.txt','Nnc')
    scaledNeAmerica=Nnc/Nac

    # admixture time
    Tadm_rand=generate_random_number('distrLogTimeAdm.txt','Tadm')
    Tadm=10**Tadm_rand # this is in terms of generations
    scaledTimeAdmixture=(10**Tadm_rand)/(4*Nac) 
    
    proportionAdmixture=generate_random_number('distrPropAdm.txt','proportionAdmixture')
    
    # ancestral North America
    Nna=generate_random_number('distrLogNeAmeBn.txt','Nna')
    scaledNeAmerica_ancestral=10**Nna/(4*Nac) # Note that America ancestral in Harris et al. was scaled by 4Ne rather than Ne 

    # growth rate america
    growthRateAmerica=-(1 / scaledTimeAdmixture) * numpy.log(scaledNeAmerica_ancestral / scaledNeAmerica)
    
    #growthRateAmerica=10**(numpy.log10(scaledNeAmerica_ancestral/scaledNeAmerica)/Tadm) -1


    # ancestral african population size
    Naa=generate_random_number('distrNeAfrAnc.txt','Naa')
    scaledNeAfricaAncestral=Naa/Nac

    # time of Africa bottleneck
    Ta=generate_random_number('distrTimeAfrBn.txt','Ta') # this is already in terms of generations (it has been multiplied by 10)
    scaledTimeAfricaExpansion=(Ta- 1000)/(4*Nac) # T_a=237227 years and 1000 comes from bottleneck duration of 1000 gens 
    scaledTimeAfricaCrash=Ta/(4*Nac)

    sev=generate_random_number('distrSeverity.txt','sev')
    scaledNeAfricaBottleneck=(1000/10**sev)/Nac # sev is in terms of duration/pop size, and duration is fixed at 1000 gen.

    lower=1.98*10**6/Nac
    upper=9.55*10**6/Nac
    

    # want adaptive mutation to be 1/Ne in population:
    # note I am not simulating selection with the Harris et al. implementation. If I do so, I need to change the starting frequency to match thiers. 
    
    starting_frequency=1/(2*scaledNeAmerica*Nac) #CHANGE

    # other hyper parameters:
    s_rand=numpy.random.uniform(0,1) #CHANGE
    age_rand=numpy.random.uniform(0,1) #CHANGE

    s=2*Nac*s_rand #CHANGE
    s2=2*s
    age=age_rand*scaledTimeAdmixture

    # print these parameters to a file 

    if selection == 'True':
        command=msmsdir + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength) +' -SAA ' + str(s2) +' -SaA ' + str(s) +' -SI ' + str(age) +' 3 0 0 ' + str(starting_frequency) + ' -SFC -Smu ' + str(adaptive_theta) + ' -Sp 0.5 -oTrace  > $file\n'

    else:
        if MS_flag=='True':
            command=msdir + str(totalSampleSize) + ' 1  -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' -seed $RANDOM $RANDOM $RANDOM  > $file\n'
            
        else:
            command=msmsdir + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) +' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + ' -en ' + str(scaledTimeAfricaCrash) + ' 1 ' + str(scaledNeAfricaAncestral) + ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' > $file\n'

    print command


#########################
def Arguello():
    
    Nac = 3980000.0 # Ne Zimbabwe present
    Tae = 20000.0 # Tsplit Netherlands
    Tadm = 179.0 #Tsplit Ithaca
    Ta = 113000 #T growth Zimbabwe
    #sev = 0.21
    Naa = 1930000.0 # Ne Africa ancestral
    Nec = 1600000.0 # Ne Netherlands present
    Nea = 37800.0 # Ne Netherlands bottleneck
    Nnc = 554000.0 #Ne Ithaca present
    Nna = 839.0 #Ne Ithaca bottleneck
    proportionAdmixture = 1-0.182 # proportion admixture from Netherlands
    #m_NA= 2*11.9 #Ithaca to Zimbabwe. Multiply by 2 because units given in 2nM whereas msdoc says 4mM
    #m_AN=2*57.9 #Zimbabwe to Ithaca
    #m_EN=2*66.7 #Netherlands to Ithaca
    #m_NE=2*64.3 #Ithaca to netherlands
    #m_AE=2*80.0 #Zimbabwe to Netherlands
    #m_EA=2*3.01 #Netherlands to Zimbabwe
    
    #m_AE_anc=2*1.74
    #m_EA_anc=2*1.78

    m_NA= 2*0.00000727553321 #Ithaca to Zimbabwe. Multiply by 2 because units given in 2nM whereas msdoc says 4mM
    m_AN=2*0.000010724282 #Zimbabwe to Ithaca
    m_EN=2*0.0000579947132 #Netherlands to Ithaca
    m_NE=2*0.0000208088427 #Ithaca to netherlands
    m_AE=2*0.000000939503141 #Zimbabwe to Netherlands
    m_EA=2*0.0000100505294 #Netherlands to Zimbabwe
    
    m_AE_anc=2*0.000000555414702 # ancestral Zimbabwe to Netherlands
    m_EA_anc=2*0.00000021912387 # ancestral Netherlands to Zimbabwe



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
    scaledNeEurope_ancestral=Nea/Nac

    # Time of Africa-Europe split
    scaledTimeSplitAfricaEurope=Tae/(4*Nac)
    
    growthRateEurope=-(1 / scaledTimeSplitAfricaEurope) * numpy.log(scaledNeEurope_ancestral / scaledNeEurope)
    
    # current population size of america
    scaledNeAmerica=Nnc/Nac

    # admixture time
    scaledTimeAdmixture=(Tadm)/(4*Nac) 
    
    # ancestral North America
    scaledNeAmerica_ancestral=Nna/Nac

    # growth rate america
    growthRateAmerica=-(1 / scaledTimeAdmixture) * numpy.log(scaledNeAmerica_ancestral / scaledNeAmerica)
    
    # time of Africa bottleneck
    scaledTimeAfricaExpansion=(Ta)/(4*Nac) 
    
    scaledNeAfricaBottleneck=(Naa)/Nac 

    # add in migration between Africa <-> NA, and NA <->Europe
    M_NA=Nac*m_NA
    M_AN=Nac*m_AN
    M_EN=Nac*m_EN
    M_NE=Nac*m_NE
    M_AE=Nac*m_AE
    M_EA=Nac*m_EA
    M_AE_anc=Nac*m_AE_anc
    M_EA_anc=Nac*m_EA_anc
    
    # parameters for selection: 

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

        command=msmsdir + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) 

        command += ' -m 1 3 ' + str(M_AN) + ' -m 3 1 ' + str(M_NA) + ' -m 2 3 ' + str(M_EN) + ' -m 3 2 ' + str(M_NE) + ' -m 1 2 ' + str(M_AE) + ' -m 2 1 ' + str(M_EA)  

        command += ' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) 

        command += ' -em ' + str(scaledTimeAdmixture) + ' 1 2 ' + str(M_AE_anc) + ' -em ' + str(scaledTimeAdmixture) + ' 2 1 ' + str(M_EA_anc)  # ancestral migration between africa and Eur before admixture

        command +=' -r ' +str(recombinationRate) + ' ' + str(locusLength) +' -SAA ' + str(s2) +' -SaA ' + str(s) +' -SI ' + str(age) +' 3 0 0 ' + str(starting_frequency) + ' -SFC -Smu ' + str(adaptive_theta) + ' -Sp 0.5 -oTrace  > $file\n'

    else:
        if MS_flag=='True':

            command=msdir  + str(totalSampleSize) + ' 1 -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) 

            command += ' -m 1 3 ' + str(M_AN) + ' -m 3 1 ' + str(M_NA) + ' -m 2 3 ' + str(M_EN) + ' -m 3 2 ' + str(M_NE) + ' -m 1 2 ' + str(M_AE) + ' -m 2 1 ' + str(M_EA)  

            command +=' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) 

            command += ' -em ' + str(scaledTimeAdmixture) + ' 1 2 ' + str(M_AE_anc) + ' -em ' + str(scaledTimeAdmixture) + ' 2 1 ' + str(M_EA_anc)  # ancestral migration between africa and Eur before admixture
        
            command += ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' -seed $RANDOM $RANDOM $RANDOM   > $file\n'

        else:

            command=msmsdir + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) 

            command += ' -m 1 3 ' + str(M_AN) + ' -m 3 1 ' + str(M_NA) + ' -m 2 3 ' + str(M_EN) + ' -m 3 2 ' + str(M_NE) + ' -m 1 2 ' + str(M_AE) + ' -m 2 1 ' + str(M_EA)  

            command +=' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) 

            command += ' -em ' + str(scaledTimeAdmixture) + ' 1 2 ' + str(M_AE_anc) + ' -em ' + str(scaledTimeAdmixture) + ' 2 1 ' + str(M_EA_anc)  # ancestral migration between africa and Eur before admixture
        
            command += ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' > $file\n'


    print command

########
# Main #
########


if model=='admixture_mode':
    admixture_mode()

if model=='constNe10e6':
    constNe10e6()

if model=='constNe2.7e6':
    constNe27e6()

if model=='dadi1':
    dadi1()

if model=='dadi2':
    dadi2()

if model=='admixture_bot_mode':
    admixture_bot_mode()

if model=='admixture_posterior':
    admixture_posterior()

if model=='admixture_posterior_Harris':
    admixture_posterior_Harris()

if model=='admixture_posterior_Harris_Nac':
    admixture_posterior_Harris_Nac()

if model=='admixture_mode_Harris':
    admixture_mode_Harris()

if model=='admixture_mode_Garud2015':
    admixture_mode_Garud2015()

if model=='admixture_bot_mode_Garud2015':
    admixture_bot_mode_Garud2015()

if model=='Arguello':
    Arguello()
