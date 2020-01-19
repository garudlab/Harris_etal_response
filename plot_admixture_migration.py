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

for m_NA in [0, 0.1, 0.25, 0.5, 0.75]: 
    summary_statistics[m_NA]={}
    for m_NE in [0, 0.1, 0.25, 0.5, 0.75]: 
        summary_statistics[m_NA][m_NE]={}
        summary_statistics[m_NA][m_NE]['H12']=[]
        summary_statistics[m_NA][m_NE]['S']=[]
        summary_statistics[m_NA][m_NE]['Pi']=[]
        summary_statistics[m_NA][m_NE]['avg_dist']=[]
        summary_statistics[m_NA][m_NE]['TajD']=[]


for m_NA in [0, 0.1, 0.25, 0.5, 0.75]: 
    for m_NE in [0, 0.1, 0.25, 0.5, 0.75]: 
        
        # H12
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_migration_mNA%s_mNE%s_rho_%s_theta_%s_selection_%s_snps.txt') %(m_NA,m_NE,rho_in, adaptive_theta, selection),'r')
        
        # iterate through the file
        # check if there are 17 lines
        # store H12 for this parameter configuration in a dictionary
        for line in inFile:
            items=line.split('\t')
            if len(items)==17:
                H12=float(items[8])
                summary_statistics[m_NA][m_NE]['H12'].append(H12)
        inFile.close()

        # Pi, S, TajD
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_migration_mNA%s_mNE%s_rho_%s_theta_%s_selection_%s_Pi_S_TajD.txt') %(m_NA,m_NE,rho_in, adaptive_theta, selection),'r')
        for line in inFile:
            items=line.split('\t')
            if len(items)==4:
                S=float(items[0])/100000
                Pi=float(items[1])/100000
                avg_dist=float(items[2])/100000
                TajD=float(items[3])
                summary_statistics[m_NA][m_NE]['S'].append(S)
                summary_statistics[m_NA][m_NE]['Pi'].append(Pi)
                summary_statistics[m_NA][m_NE]['avg_dist'].append(avg_dist)
                summary_statistics[m_NA][m_NE]['TajD'].append(TajD)

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

    plt.boxplot([summary_statistics[0][0][summary_statistic],summary_statistics[0.1][0.1][summary_statistic],summary_statistics[0.25][0.25][summary_statistic],summary_statistics[0.5][0.5][summary_statistic],summary_statistics[0.75][0.75][summary_statistic]], labels=['0','0.1','0.25','0.5','0.75'])

    plt.axhline(y=data_means[summary_statistic],color='r', label='Mean value in data')

    plt.ylim(min(data_means[summary_statistic], min(summary_statistics[0][0][summary_statistic]), min(summary_statistics[0.1][0.1][summary_statistic]), min(summary_statistics[0.25][0.25][summary_statistic]),min(summary_statistics[0.5][0.5][summary_statistic]), min(summary_statistics[0.75][0.75][summary_statistic]))-0.01, max(data_means[summary_statistic], max(summary_statistics[0][0][summary_statistic]), max(summary_statistics[0.1][0.1][summary_statistic]), max(summary_statistics[0.25][0.25][summary_statistic]),max(summary_statistics[0.5][0.5][summary_statistic]), max(summary_statistics[0.75][0.75][summary_statistic]) )+0.01)

    plt.xlabel('m_ij')
    plt.ylabel(summary_statistic)
    plot_no+=1

plt.tight_layout()

plt.savefig(os.path.expanduser('~/Jensen_response/analysis/admixture_migration_boxplot.png'))




##############

summary_statistics_means={}
summary_statistics_lower={}
summary_statistics_upper={}
for summary_statistic in ['H12', 'Pi', 'S', 'avg_dist', 'TajD']:
    summary_statistics_means[summary_statistic]=[]
    summary_statistics_lower[summary_statistic]=[]
    summary_statistics_upper[summary_statistic]=[]

for m_ij in [0, 0.1, 0.25, 0.5, 0.75]:
    for summary_statistic in ['H12', 'Pi', 'S', 'avg_dist', 'TajD']:
        summary_statistics_means[summary_statistic].append(numpy.mean(summary_statistics[m_ij][m_ij][summary_statistic]))
        summary_statistics_lower[summary_statistic].append(numpy.percentile(summary_statistics[m_ij][m_ij][summary_statistic], 2.5))
        summary_statistics_upper[summary_statistic].append(numpy.percentile(summary_statistics[m_ij][m_ij][summary_statistic], 97.5))

plot_no=1
for summary_statistic in ['H12', 'Pi', 'S', 'TajD']:
    plt.subplot(2, 2, plot_no)
    plt.fill_between([0, 0.1, 0.25, 0.5, 0.75], summary_statistics_lower[summary_statistic], summary_statistics_upper[summary_statistic], color='b', alpha=.25)
    plt.plot([0, 0.1, 0.25, 0.5, 0.75], summary_statistics_means[summary_statistic], 'o-',label='Simulations')
    plt.axhline(y=data_means[summary_statistic],color='r', label='Mean value in data')
    plt.ylim(min(data_means[summary_statistic], min(summary_statistics_means[summary_statistic]))-0.01, max(data_means[summary_statistic],max(summary_statistics_means[summary_statistic]))+0.01)

    #plt.legend(loc='upper right')
    plt.xlabel('m_ij')
    plt.ylabel(summary_statistic)
    #plt.title('Admixture + migration, 400 SNPs')
    plot_no+=1

plt.tight_layout()

plt.savefig(os.path.expanduser('~/Jensen_response/analysis/admixture_migration_1.png'))
