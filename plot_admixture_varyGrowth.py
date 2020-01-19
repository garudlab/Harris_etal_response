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
locus_length=100000

summary_statistics={} # store the summary statistics for all the models here. Key = [m_NA][m_NE][H12,S,Pi,avg_dist]

#NAms=[158, 2500, 61659, 1110000, 15984500, 28000000]
#NAms=[158,2500, 61659, 1110000, 15984500, 28000000]
#EUrs=[3801, 16982, 67608, 39000, 3122470, 9550000]
EUr_ancs=[16982, 39000, 67608, 100000,500000, 1000000, 200000]
#EUr_curs=[67608, 100000, 200000, 1000000, 2000000, 3122470, 500000, 9550000]
EUr_curs=[1000000,2000000,3122470,9550000]

for EUr_anc in EUr_ancs: 
    summary_statistics[EUr_anc]={}
    for EUr_cur in EUr_curs:
        summary_statistics[EUr_anc][EUr_cur]={}
        summary_statistics[EUr_anc][EUr_cur]['H12']=[]
        summary_statistics[EUr_anc][EUr_cur]['S']=[]
        summary_statistics[EUr_anc][EUr_cur]['Pi']=[]
        summary_statistics[EUr_anc][EUr_cur]['avg_dist']=[]
        summary_statistics[EUr_anc][EUr_cur]['TajD']=[]


NAm=2500
#for NAm in NAms: 
for EUr_anc in EUr_ancs: 
    for EUr_cur in EUr_curs: 
        # H12

        #inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_varyGrowth_NAm_anc%s_NAm_cur%s_EUr_anc%s_EUr_cur%s_rho_%s_theta_%s_selection_%s_MS_snps.txt') %(NAm, NAm, EUr, EUr, rho_in, adaptive_theta, selection),'r')

        inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_varyGrowth_NAm_anc%s_NAm_cur%s_EUr_anc%s_EUr_cur%s_rho_%s_theta_%s_selection_%s_MS_snps.txt') %(1110000, 15984500, EUr_anc, EUr_cur, rho_in, adaptive_theta, selection),'r')
       

        # iterate through the file
        # check if there are 17 lines
        # store H12 for this parameter configuration in a dictionary
        for line in inFile:
            items=line.split('\t')
            if len(items)==17:
                H12=float(items[8])
                summary_statistics[EUr_anc][EUr_cur]['H12'].append(H12)
        inFile.close()

        # Pi, S, TajD
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_varyGrowth_NAm_anc%s_NAm_cur%s_EUr_anc%s_EUr_cur%s_rho_%s_theta_%s_selection_%s_MS_Pi_S_TajD.txt') %(1110000, 15984500, EUr_anc, EUr_cur,rho_in, adaptive_theta, selection),'r')
        for line in inFile:
            items=line.split('\t')
            if len(items)==4:
                if items[0] !='':
                    S=float(items[0])/locus_length
                    Pi=float(items[1])/locus_length
                    avg_dist=float(items[2])/locus_length
                    TajD=float(items[3])
                    summary_statistics[EUr_anc][EUr_cur]['S'].append(S)
                    summary_statistics[EUr_anc][EUr_cur]['Pi'].append(Pi)
                    summary_statistics[EUr_anc][EUr_cur]['avg_dist'].append(avg_dist)
                    summary_statistics[EUr_anc][EUr_cur]['TajD'].append(TajD)

        inFile.close()
            
print('data is loaded...')

###########
 
# mean_H12 in the data:
data_means={}
data_means['H12']=0.01767854
data_means['TajD']=0.4
data_means['S']=0.058
data_means['Pi']=0.012

##########
print('plotting...')

# make boxplots of the distributions
matplotlib.rcParams.update({'font.size': 7})


for EUr_anc in EUr_ancs:
    print('EUr_anc %s' %EUr_anc)

    plot_no=1
    for summary_statistic in ['H12', 'Pi', 'S', 'TajD']:
        plt.subplot(2, 2, plot_no)

        simulations=[]


        plt.boxplot([summary_statistics[EUr_anc][1000000][summary_statistic],summary_statistics[EUr_anc][2000000][summary_statistic],summary_statistics[EUr_anc][3122470][summary_statistic],summary_statistics[EUr_anc][9550000][summary_statistic]], labels=['1000000','2000000','3122470','9550000'])

        plt.axhline(y=data_means[summary_statistic],color='r', label='Mean value in data')

        if plot_no==1:
            plt.ylim(0,1)
            plt.title('Europe ancestral Ne = %s' %EUr_anc)
        elif plot_no==2:
            plt.ylim(0,0.2)
        elif plot_no==3:
            plt.ylim(0,0.4)
        elif plot_no==4:
            plt.ylim(-3,3)

        #plt.ylim(min(data_means[summary_statistic], min(summary_statistics[158][EUr_mode][summary_statistic]), min(summary_statistics[2500][EUr_mode][summary_statistic]), min(summary_statistics[61659][EUr_mode][summary_statistic]),min(summary_statistics[1110000][EUr_mode][summary_statistic]), min(summary_statistics[15984500][EUr_mode][summary_statistic]), min(summary_statistics[28000000][EUr_mode][summary_statistic]))-0.01, max(data_means[summary_statistic], max(summary_statistics[158][EUr_mode][summary_statistic]), max(summary_statistics[2500][EUr_mode][summary_statistic]), max(summary_statistics[61659][EUr_mode][summary_statistic]),max(summary_statistics[1110000][EUr_mode][summary_statistic]), max(summary_statistics[15984500][EUr_mode][summary_statistic]), max(summary_statistics[28000000][EUr_mode][summary_statistic]))+0.01)

        plt.xlabel('Europe current Ne')
        plt.ylabel(summary_statistic)
        plot_no+=1

        plt.tight_layout()

        plt.savefig(os.path.expanduser('~/Jensen_response/analysis/admixture_fixedPopSize_boxplot_EUr_anc%s.png' %EUr_anc))

    plt.close()
