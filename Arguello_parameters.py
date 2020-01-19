import numpy
import sys
from scipy.stats import rv_discrete
import os.path

rho_in=float(sys.argv[1])
adaptive_theta=float(sys.argv[2])
selection=sys.argv[3] # boolean indicating if we are similating with selection
locusLength=int(sys.argv[4])
m_NA=float(sys.argv[5]) # NEW
m_NE=float(sys.argv[6]) # NEW
 
# input parameters:
#rho_in=10**-8


numberOfSimulations=1

for i in range(0, numberOfSimulations):
    
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
    m_NA= 2*11.9 #Ithaca to Zimbabwe. Multiply by 2 because units given in 2nM whereas msdoc says 4mM
    m_AN=2*57.9 #Zimbabwe to Ithaca
    m_EN=2*66.7 #Netherlands to Ithaca
    m_NE=2*64.3 #Ithaca to netherlands
    m_AE=2*80.0 #Zimbabwe to Netherlands
    m_EA=2*3.01 #Netherlands to Zimbabwe
    
    m_AE_anc=2*1.74
    m_EA_anc=2*1.78

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
        command='~/./msms/bin/msms ' + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) 
        command += ' -m 1 3 ' + str(M_AN) + ' -m 3 1 ' + str(M_NA) + ' -m 2 3 ' + str(M_EN) + ' -m 3 2 ' + str(M_NE) + ' -m 1 2 ' + str(M_AE) + ' -m 2 1 ' + str(M_EA)  
        command += ' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) + 

        command += ' -em ' + scaledTimeAdmixture 1 2 M_AE_anc + ' -em ' + scaledTimeAdmixture 2 1 M_EA_anc  # ancestral migration between africa and Eur before admixture

        command +=' -r ' +str(recombinationRate) + ' ' + str(locusLength) +' -SAA ' + str(s2) +' -SaA ' + str(s) +' -SI ' + str(age) +' 3 0 0 ' + str(starting_frequency) + ' -SFC -Smu ' + str(adaptive_theta) + ' -Sp 0.5 -oTrace  > $file\n'

    else:
        command='~/./msms/bin/msms ' + str(totalSampleSize) + ' 1 -N ' + str(Nac) + ' -t ' +str(theta) +' -I 3 '+ str(sampleSizeAfrica) +' ' + str(sampleSizeEurope) +' ' +str(sampleSizeAmerica) +' -n 2 '+ str(scaledNeEurope) +' -g 2 ' + str(growthRateEurope) + ' -n 3 ' + str(scaledNeAmerica) + ' -g 3 ' + str(growthRateAmerica) 
        command += ' -m 1 3 ' + str(M_AN) + ' -m 3 1 ' + str(M_NA) + ' -m 2 3 ' + str(M_EN) + ' -m 3 2 ' + str(M_NE) + ' -m 1 2 ' + str(M_AE) + ' -m 2 1 ' + str(M_EA)  
        command +=' -es ' + str(scaledTimeAdmixture) + ' 3 ' + str(proportionAdmixture) +' -ej ' +str(scaledTimeAdmixture) + ' 3 2 -ej ' + str(scaledTimeAdmixture) + ' 4 1 -ej ' + str(scaledTimeSplitAfricaEurope) +' 2 1 -en ' + str(scaledTimeAfricaExpansion) + ' 1 ' + str(scaledNeAfricaBottleneck) 

        command += ' -em ' + scaledTimeAdmixture 1 2 M_AE_anc + ' -em ' + scaledTimeAdmixture 2 1 M_EA_anc  # ancestral migration between africa and Eur before admixture
        
        command += ' -r ' +str(recombinationRate) + ' ' + str(locusLength)  + ' > $file\n'

    print command + '\n' + str(Nac) + '\n' + str(s) + '\n' + str(age)
