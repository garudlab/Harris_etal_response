import matplotlib
matplotlib.use('Agg')
import numpy
import sys
from scipy.stats import rv_discrete
import os.path
import numpy as np

import matplotlib.pyplot as plt


#rho_in=float(sys.argv[1])
#adaptive_theta=float(sys.argv[2])
#selection=sys.argv[3] # boolean indicating if we are similating with selection
#locusLength=int(sys.argv[4])

rho_in='5e-9'
adaptive_theta='0'
selection='False'


summary_statistics={} # store the summary statistics for all the models here. Key = [m_NA][m_NE][H12,S,Pi,avg_dist]

proportions=[0, 0.1, 0.25, 0.4, 0.5, 0.75, 0.85, 0.9]

for prop in proportions: 
    summary_statistics[prop]={}
    summary_statistics[prop]['H12']=[]
    summary_statistics[prop]['S']=[]
    summary_statistics[prop]['Pi']=[]
    summary_statistics[prop]['avg_dist']=[]
    summary_statistics[prop]['TajD']=[]


for prop in proportions: 
        
        # H12
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_migration_diffProps%s_rho_%s_theta_%s_selection_%s_MS_snps.txt') %(prop,rho_in, adaptive_theta, selection),'r')
        
        # iterate through the file
        # check if there are 17 lines
        # store H12 for this parameter configuration in a dictionary
        for line in inFile:
            items=line.split('\t')
            if len(items)==17:
                H12=float(items[8])
                summary_statistics[prop]['H12'].append(H12)
        inFile.close()

        # Pi, S, TajD
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_migration_diffProps%s_rho_%s_theta_%s_selection_%s_MS_Pi_S_TajD.txt') %(prop,rho_in, adaptive_theta, selection),'r')
        for line in inFile:
            items=line.split('\t')
            if len(items)==4:
                S=float(items[0])/100000
                Pi=float(items[1])/100000
                avg_dist=float(items[2])/100000
                TajD=float(items[3])
                summary_statistics[prop]['S'].append(S)
                summary_statistics[prop]['Pi'].append(Pi)
                summary_statistics[prop]['avg_dist'].append(avg_dist)
                summary_statistics[prop]['TajD'].append(TajD)

        inFile.close()
            


###########
 
# mean_H12 in the data:
data_means={}
data_means['H12']=0.01767854
data_means['TajD']=0.4
data_means['S']=0.058
data_means['Pi']=0.012

##########
# make boxplots of the distributions

plot_no=1
for summary_statistic in ['H12', 'Pi', 'S', 'TajD']:
    plt.subplot(2, 2, plot_no)

    simulations=[]


    plt.boxplot([summary_statistics[0][summary_statistic],summary_statistics[0.1][summary_statistic],summary_statistics[0.25][summary_statistic],summary_statistics[0.4][summary_statistic],summary_statistics[0.5][summary_statistic], summary_statistics[0.75][summary_statistic], summary_statistics[0.85][summary_statistic],summary_statistics[0.9][summary_statistic]], labels=['0','0.1','0.25','0.4','0.5','0.75','0.85','0.9'])

    plt.axhline(y=data_means[summary_statistic],color='r', label='Mean value in data')

    plt.ylim(min(data_means[summary_statistic], min(summary_statistics[0][summary_statistic]),min(summary_statistics[0.1][summary_statistic]),min(summary_statistics[0.25][summary_statistic]),min(summary_statistics[0.4][summary_statistic]),min(summary_statistics[0.5][summary_statistic]), min(summary_statistics[0.75][summary_statistic]), min(summary_statistics[0.85][summary_statistic]),min(summary_statistics[0.9][summary_statistic]))-0.01,  max(data_means[summary_statistic], max(summary_statistics[0][summary_statistic]),max(summary_statistics[0.1][summary_statistic]),max(summary_statistics[0.25][summary_statistic]),max(summary_statistics[0.4][summary_statistic]),max(summary_statistics[0.5][summary_statistic]), max(summary_statistics[0.75][summary_statistic]), max(summary_statistics[0.85][summary_statistic]),max(summary_statistics[0.9][summary_statistic]))+0.01)

    plt.xlabel('proportion admixture')
    plt.ylabel(summary_statistic)
    plot_no+=1

plt.tight_layout()

plt.savefig(os.path.expanduser('~/Jensen_response/analysis/admixture_diffProps_boxplot.png'))




