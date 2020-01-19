import matplotlib
matplotlib.use('Agg')
import numpy
import sys
from scipy.stats import rv_discrete
import os.path
import numpy as np

import matplotlib.pyplot as plt


rho_in='5e-9'
adaptive_theta='0'
selection='False'


def isFloat(value):
    try:
        float(value)
        return True
    except:
        return False


summary_statistics={} # store the summary statistics for all the models here. Key = [m_NA][m_NE][H12,S,Pi,avg_dist]
LD_dict={}

NAms=[2500, 61659, 1110000, 15984500, 28000000]
EUrs=[16982, 67608, 700000, 2000000, 9550000]

for NAm in NAms: 
    summary_statistics[NAm]={}
    LD_dict[NAm]={}
    for EUr in EUrs: 
        summary_statistics[NAm][EUr]={}
        summary_statistics[NAm][EUr]['H12']=[]
        summary_statistics[NAm][EUr]['H2H1']=[]
        summary_statistics[NAm][EUr]['S']=[]
        summary_statistics[NAm][EUr]['Pi']=[]
        summary_statistics[NAm][EUr]['avg_dist']=[]
        summary_statistics[NAm][EUr]['TajD']=[]
        LD_dict[NAm][EUr]={}
        LD_dict[NAm][EUr]['position']=[]
        LD_dict[NAm][EUr]['R2']=[]


for NAm in NAms: 
    for EUr in EUrs: 
        # H12
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_fixedPopSize_NAm%s_EUr%s_rho_%s_theta_%s_selection_%s_MS_snps.txt') %(NAm, EUr,rho_in, adaptive_theta, selection),'r')
        
        # iterate through the file
        # check if there are 17 lines
        # store H12 for this parameter configuration in a dictionary
        for line in inFile:
            items=line.split('\t')
            if len(items)==17:
                H12=float(items[8])
                H2H1=float(items[9])
                summary_statistics[NAm][EUr]['H12'].append(H12)
                summary_statistics[NAm][EUr]['H2H1'].append(H2H1)
        inFile.close()

        # Pi, S, TajD
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_fixedPopSize_NAm%s_EUr%s_rho_%s_theta_%s_selection_%s_MS_Pi_S_TajD_perBp.txt') %(NAm, EUr,rho_in, adaptive_theta, selection),'r')
        for line in inFile:
            items=line.strip('\n').split('\t')
            if len(items)==4:
                if isFloat(items[0]) and isFloat(items[1]) and isFloat(items[2]) and isFloat(items[3]):
                    S=float(items[0])
                    Pi=float(items[1])
                    avg_dist=float(items[2])
                    TajD=float(items[3])
                    summary_statistics[NAm][EUr]['S'].append(S)
                    summary_statistics[NAm][EUr]['Pi'].append(Pi)
                    summary_statistics[NAm][EUr]['avg_dist'].append(avg_dist)
                    summary_statistics[NAm][EUr]['TajD'].append(TajD)

        inFile.close()


        # LD
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_fixedPopSize_NAm%s_EUr%s_rho_5e-9_theta_0_selection_False_MS_LD_0.05_0.95_averaged.txt' %(NAm, EUr)),'r')

        for line in inFile:
            items=line.strip().split('\t')
            position=items[0]
            R2=items[1]
            LD_dict[NAm][EUr]['position'].append(int(position))
            LD_dict[NAm][EUr]['R2'].append(float(R2))


            
print('data is loaded...')

###########
 
# mean_H12 in the data:
data_means={}
data_means['H12']=0.01767854
data_means['H2H1']=0.8301393
data_means['TajD']=0.3835
data_means['S']=0.0578
data_means['Pi']=0.0118

# data medians
data_medians={}
data_medians['H12']=0.01545779
data_medians['H2H1']=0.8651685
data_medians['TajD']=0.3835
data_medians['S']=0.0578
data_medians['Pi']=0.0118

# load LD estimated directly from the genome
LD_dict['DGRP']={}
LD_dict['DGRP']['position']=[]
LD_dict['DGRP']['R2']=[] 

inFile=open(os.path.expanduser('~/Jensen_response/data/alldata_averaged_sorted'),'r')

for line in inFile:
    items=line.strip('\n').split('\t')
    position=items[0]
    R2=items[1]
    LD_dict['DGRP']['position'].append(int(position))
    LD_dict['DGRP']['R2'].append(float(R2))
    


##########
print('plotting...')

# make boxplots of the distributions
matplotlib.rcParams.update({'font.size': 7})

NAm_mode=15984500

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 9

for EUr_mode in EUrs:
    print('EUr %s' %EUr_mode)

    plot_no=1
    for summary_statistic in ['H12', 'H2H1', 'Pi', 'S', 'TajD']:
        plt.subplot(2, 3, plot_no)

        simulations=[]

        #bp = plt.boxplot([summary_statistics[2500][EUr_mode][summary_statistic],summary_statistics[61659][EUr_mode][summary_statistic],summary_statistics[1110000][EUr_mode][summary_statistic],summary_statistics[15984500][EUr_mode][summary_statistic], summary_statistics[28000000][EUr_mode][summary_statistic]], 0, 'k.', labels=['2500' ,'61659' ,'1110000', '15984500' ,'28000000'])
        bp = plt.boxplot([summary_statistics[61659][EUr_mode][summary_statistic],summary_statistics[1110000][EUr_mode][summary_statistic],summary_statistics[1110000][EUr_mode][summary_statistic]], 0, 'k.', labels=['61659' ,'1110000','15984500'])

        for element in ['whiskers', 'fliers', 'means', 'medians', 'caps']:
            plt.setp(bp[element], color='black')

        plt.axhline(y=data_medians[summary_statistic],color='r', label='Mean value in data')

        if plot_no==1:
            plt.ylim(0,.1)
            plt.title('Europe Ne = %s' %EUr_mode)
        elif plot_no==2:
            plt.ylim(0,1)
        elif plot_no==3:
            plt.ylim(0,0.03)
        elif plot_no==4:
            plt.ylim(0,0.1)
        elif plot_no==5:
            plt.ylim(-3,3)

        #plt.ylim(min(data_means[summary_statistic], min(summary_statistics[158][EUr_mode][summary_statistic]), min(summary_statistics[2500][EUr_mode][summary_statistic]), min(summary_statistics[61659][EUr_mode][summary_statistic]),min(summary_statistics[1110000][EUr_mode][summary_statistic]), min(summary_statistics[15984500][EUr_mode][summary_statistic]), min(summary_statistics[28000000][EUr_mode][summary_statistic]))-0.01, max(data_means[summary_statistic], max(summary_statistics[158][EUr_mode][summary_statistic]), max(summary_statistics[2500][EUr_mode][summary_statistic]), max(summary_statistics[61659][EUr_mode][summary_statistic]),max(summary_statistics[1110000][EUr_mode][summary_statistic]), max(summary_statistics[15984500][EUr_mode][summary_statistic]), max(summary_statistics[28000000][EUr_mode][summary_statistic]))+0.01)

        plt.xlabel('North American Population Size')
        plt.ylabel(summary_statistic)

        plot_no+=1

        
    ###################
    # Plot LD         # 
    ###################

    # plot LD
    plt.subplot(2, 3, 6)

    colors=['green','blue', 'red','yellow', 'purple', 'orange', 'brown']

    i=0
    for NAm in NAms:
        plt.semilogy(LD_dict[NAm][EUr_mode]['position'], LD_dict[NAm][EUr_mode]['R2'], label=NAm, color=colors[i])
        plt.xlim(0,10000)
        plt.ylim(0.005,1.05)
        plt.xlabel('Position')
        plt.ylabel('log(R^2)')
        i +=1
        
    plt.semilogy(LD_dict['DGRP']['position'], LD_dict['DGRP']['R2'], linestyle='dashed',color='black', label='DGRP')

    #plt.legend()

    plt.tight_layout()
    
    plt.savefig(os.path.expanduser('~/Jensen_response/analysis/admixture_fixedPopSize_boxplot_EUr%s.png' %EUr_mode), dpi=600)

    plt.close()








