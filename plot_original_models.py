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



def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color, linestyle='-')
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)    
    plt.setp(bp['fliers'], color=color, markersize=3)

#################


def isFloat(value):
    try:
        float(value)
        return True
    except:
        return False

#####################

def initialize_model_dict(summary_statistics,models):
    for model in models: 
        summary_statistics[model]={}
        summary_statistics[model]['H12']=[]
        summary_statistics[model]['H2H1']=[]
        summary_statistics[model]['S']=[]
        summary_statistics[model]['Pi']=[]
        summary_statistics[model]['avg_dist']=[]
        summary_statistics[model]['TajD']=[]

    return summary_statistics


summary_statistics={} # store the summary statistics for all the models here. Key = [m_NA][m_NE][H12,S,Pi,avg_dist]

models=['constNe10e6', 'constNe2.7e6', 'dadi1', 'dadi2', 'admixture_mode_Garud2015', 'admixture_bot_mode_Garud2015', 'admixture_mode', 'admixture_bot_mode', 'admixture_posterior', 'admixture_posterior_Harris', 'Arguello' ]

summary_statistics=initialize_model_dict(summary_statistics,models)



##########
for model in models: 
        
    if model !='DGRP':

        # H12
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/%s_neutrality_locusLen100000_snps.txt' %model),'r')
        
        # iterate through the file
        # check if there are 17 lines
        # store H12 for this parameter configuration in a dictionary
        for line in inFile:
            items=line.split('\t')
            if len(items)==17:
                H12=float(items[8])
                H2H1=float(items[9])
                summary_statistics[model]['H12'].append(H12)
                summary_statistics[model]['H2H1'].append(H2H1)
        inFile.close()

        # Pi, S, TajD
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/%s_neutrality_locusLen11000_Pi_S_TajD_perBp.txt' %model),'r')
        for line in inFile:
            items=line.strip('\n').split('\t')
            if len(items)==4:
                if isFloat(items[0]):
                    S=float(items[0])
                    Pi=float(items[1])
                    avg_dist=float(items[2])
                    TajD=float(items[3])
                    summary_statistics[model]['S'].append(S)
                    summary_statistics[model]['Pi'].append(Pi)
                    summary_statistics[model]['avg_dist'].append(avg_dist)
                    summary_statistics[model]['TajD'].append(TajD)

        inFile.close()
            


########
# add in the new models, which have a different file name format
EUrs=[700000]
NAms=[1110000, 15984500]


for EUr in EUrs:
    for NAm in NAms:

        model='%s_%s' %(NAm, EUr)

        summary_statistics=initialize_model_dict(summary_statistics,[model])
        models.append(model)

        inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_fixedPopSize_NAm%s_EUr%s_rho_5e-9_theta_0_selection_False_MS_snps.txt') %(NAm, EUr),'r')

        # iterate through the file
        # check if there are 17 lines
        # store H12 for this parameter configuration in a dictionary
        for line in inFile:
            items=line.split('\t')
            if len(items)==17:
                H12=float(items[8])
                H2H1=float(items[9])
                summary_statistics[model]['H12'].append(H12)
                summary_statistics[model]['H2H1'].append(H2H1)
        inFile.close()

        # Pi, S, TajD
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_fixedPopSize_NAm%s_EUr%s_rho_5e-9_theta_0_selection_False_MS_Pi_S_TajD_perBp.txt' %(NAm, EUr)),'r')
        for line in inFile:
            items=line.strip('\n').split('\t')
            if len(items)==4:
                if isFloat(items[0]):
                    S=float(items[0])
                    Pi=float(items[1])
                    avg_dist=float(items[2])
                    TajD=float(items[3])
                    summary_statistics[model]['S'].append(S)
                    summary_statistics[model]['Pi'].append(Pi)
                    summary_statistics[model]['avg_dist'].append(avg_dist)
                    summary_statistics[model]['TajD'].append(TajD)

        inFile.close()
            




############################
# load statistics for data #
############################
 
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


# data distribution of H12

DGRP_H12=[]
for chrom in ['2R','2L','3R','3L']:
    inFile=open(os.path.expanduser('~/Jensen_response/analysis/H12_scan/clusters_082712_400_50_2_2_0_chr%s_withRho_threshold.txt' %(chrom)),'r')
    for line in inFile:
        items=line.strip().split('\t')
        H12=float(items[27])
        DGRP_H12.append(H12)

DGRP_H12=np.asarray(DGRP_H12)

print len(DGRP_H12)

# read in peaks
inFile=open(os.path.expanduser('~/Jensen_response/analysis/H12_scan/peaks_082712_400_50_2_2_0_pan_509_top50.txt'),'r')

peaks=[]
for line in inFile:  
    items=line.strip().split('\t')
    H12=float(items[27])
    peaks.append(H12)

######################################
# make boxplots of the distributions #
######################################

colors=['#fed976','#feb24c', '#fd8d3c','#fc4e2a', '#e31a1c', '#b10026', '#bdd7e7','#6baed6', '#08519c','#8c510a','pink','#d9ef8b', '#66bd63']


matplotlib.rcParams.update({'font.size': 9})
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 9


plot_no=1
for summary_statistic in ['Pi', 'S', 'H12']:
    plt.subplot(2, 3, plot_no)

    if summary_statistic=='H12':

        bp = plt.boxplot(DGRP_H12, positions=[1], sym='.',patch_artist=True, labels=['DGRP'], widths=0.5)
        set_box_color(bp, 'black')

        col_i=0
        for model in models:
            bp = plt.boxplot(summary_statistics[model][summary_statistic], positions=[col_i+2], sym='.',patch_artist=True, widths=0.5, labels=[''])
            #patch.set_facecolor(colors[col_i])
            set_box_color(bp, colors[col_i]) 
            col_i+=1
        # dummy boxplot to get the colors to work for last one
        bp = plt.boxplot(summary_statistics[model][summary_statistic], positions=[col_i+10], sym='.',patch_artist=True, widths=0.5, labels=[''])
        #patch.set_facecolor('white')
        set_box_color(bp, 'white') 
        

        plt.xlim(0,len(models)+2)

        plt.plot( [1]*len(peaks),peaks, 'r.')

    else:


        col_i=0
        for model in models:
            bp = plt.boxplot(summary_statistics[model][summary_statistic], positions=[col_i+1], sym='.',patch_artist=True, widths=0.5, labels=[''])
            #patch.set_facecolor(colors[col_i])
            set_box_color(bp, colors[col_i]) 
            col_i+=1
        # dummy boxplot to get the colors to work for last one
        bp = plt.boxplot(summary_statistics[model][summary_statistic], positions=[col_i+10], sym='.',patch_artist=True, widths=0.5, labels=[''])
        #patch.set_facecolor('white')
        set_box_color(bp, 'white') 
        

        plt.xlim(0,len(models)+2)


    plt.axhline(y=data_medians[summary_statistic],linestyle='dashed',color='black', label='Median value in data')

    if summary_statistic=='H12':
        plt.ylim(0,1)
        plt.ylabel('H12 (401 SNP windows)')
    elif summary_statistic=='Pi':
        plt.ylim(0,0.03)
        plt.ylabel('Pi/bp')  
    elif summary_statistic=='S':
        plt.ylim(0,0.1)
        plt.ylabel('S/bp')   

    plot_no+=1



##########################
# Plot LD                # 
##########################

# load the data

LD_dict={}

for model in models:

    LD_dict[model]={}
    LD_dict[model]['position']=[]
    LD_dict[model]['R2']=[]

    if model =='1110000_700000':
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_fixedPopSize_NAm1110000_EUr700000_rho_5e-9_theta_0_selection_False_MS_LD_0.05_0.95_averaged.txt'))
    elif  model == '15984500_700000':
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_fixedPopSize_NAm15984500_EUr700000_rho_5e-9_theta_0_selection_False_MS_LD_0.05_0.95_averaged.txt'))
    else:
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/%s_neutrality_locusLen11000_LD_0.05_0.95_averaged.txt' %model),'r')


    for line in inFile:
        items=line.strip().split('\t')
        position=items[0]
        R2=items[1]
        LD_dict[model]['position'].append(int(position))
        LD_dict[model]['R2'].append(float(R2))



# load LD estimated directly from the genome


LD_dict['DGRP']={}
LD_dict['DGRP']['position']=[]
LD_dict['DGRP']['R2']=[] 

inFile=open(os.path.expanduser('~/Jensen_response/analysis/DGRP_LD/alldata_averaged_sorted'),'r')

for line in inFile:
    items=line.strip('\n').split('\t')
    position=items[0]
    R2=items[1]
    LD_dict['DGRP']['position'].append(int(position))
    LD_dict['DGRP']['R2'].append(float(R2))
    
#colors=['green','blue', 'red','yellow', 'purple', 'orange', 'brown','pink', 'turquoise','magenta', 'lightgray','orchid','y']

# plot LD
for xlim_range in [100,10000]:
    plt.subplot(2, 3, plot_no)
    
    i=0
    for model in models:
        
        plt.semilogy(LD_dict[model]['position'], LD_dict[model]['R2'], label=model, color=colors[i])
        plt.xlim(0,xlim_range) 
        if xlim_range == 100:
            plt.ylim(0.1,1.05)
        else:
            plt.ylim(0.005,1.05)
        plt.xlabel('Distance between SNPs (bps)')
        plt.ylabel('log($R^2$)')
        i +=1

    plt.semilogy(LD_dict['DGRP']['position'], LD_dict['DGRP']['R2'], linestyle='dashed',color='black', label='DGRP')    
       
    plot_no +=1


print(plot_no)
plt.subplot(2, 3, plot_no)
# put some text in plot number 6
plt.text(0.045, 1.0, "DGRP data",  rotation=90.)
plt.text(0.113, 1.0, "Constant $Ne$=$10^6$",  rotation=90.)
plt.text(0.18, 1.0, "Constant $Ne$=$2.7*10^6$",  rotation=90.)
plt.text(0.248, 1.0, "Severe, short bottleneck",  rotation=90.)
plt.text(0.315, 1.0, "Shallow, long bottleneck",  rotation=90.)
plt.text(0.383, 1.0, "Admixture (mode), Garud $\it{et al.}$ 2015 implementation",  rotation=90.)
plt.text(0.45, 1.0, "Admixture + bottleneck (mode), Garud $\it{et al.}$ 2015 implementation",  rotation=90.)
plt.text(0.518, 1.0, "Admixture (mode), Duchen $\it{et al.}$ 2013",  rotation=90.)
plt.text(0.585, 1.0, "Admixture + bottleneck (mode), Duchen $\it{et al.}$ 2013",  rotation=90.)
plt.text(0.653, 1.0, "Admixture (posterior distribution), Duchen $\it{et al.}$ 2013",  rotation=90.)
plt.text(0.72, 1.0, "Admixture (posterior distribution), Harris $\it{et al.}$ 2018 implementation",  rotation=90.)
plt.text(0.788, 1.0, "Arguello $\it{et al.}$ 2019",  rotation=90.)
plt.text(0.856, 1.0, "Admixture, Eur=$7*10^5$, Nac=$1.11*10^6$",  rotation=90.)
plt.text(0.925, 1.0, "Admixture, Eur=$7*10^5$, Nac=$15.9*10^6$",  rotation=90.)

plt.axis('off')


plt.tight_layout()

plt.savefig(os.path.expanduser('~/Jensen_response/analysis/original_models_LD.png'), dpi=600)

plt.close()        




#########
#########

#################################################
# make a plot of pi, s, LD                      #
#################################################

colors=['#fed976','#feb24c', '#fd8d3c','#fc4e2a', '#e31a1c', '#b10026', '#bdd7e7','#6baed6', '#08519c','#8c510a','pink','#d9ef8b', '#66bd63']


matplotlib.rcParams.update({'font.size': 9})
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 6
fig_size[1] = 6


plot_no=1
for summary_statistic in ['Pi', 'S']:
    plt.subplot(2,2, plot_no)

    col_i=0
    for model in models:
        bp = plt.boxplot(summary_statistics[model][summary_statistic], positions=[col_i], sym='.',patch_artist=True, widths=0.5, labels=[''])
        set_box_color(bp, colors[col_i])
        col_i+=1
    # dummy boxplot to get the colors to work for last one
    bp = plt.boxplot(summary_statistics[model][summary_statistic], positions=[col_i+10], sym='.',patch_artist=True, widths=0.5, labels=[''])
    set_box_color(bp, 'white')   
        
    plt.xlim(0,len(models)+2)


    plt.axhline(y=data_medians[summary_statistic],linestyle='dashed',color='black', label='Median value in data')

    if summary_statistic=='Pi':
        plt.ylim(0,0.03)
        plt.ylabel('Pi/bp')  
    elif summary_statistic=='S':
        plt.ylim(0,0.1)
        plt.ylabel('S/bp')   

    plot_no+=1





##########################
# Plot LD                # 
##########################

# load the data

LD_dict={}

for model in models:

    LD_dict[model]={}
    LD_dict[model]['position']=[]
    LD_dict[model]['R2']=[]

    if model =='1110000_700000':
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_fixedPopSize_NAm1110000_EUr700000_rho_5e-9_theta_0_selection_False_MS_LD_0.05_0.95_averaged.txt'))
    elif  model == '15984500_700000':
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_fixedPopSize_NAm15984500_EUr700000_rho_5e-9_theta_0_selection_False_MS_LD_0.05_0.95_averaged.txt'))
    else:
        inFile=open(os.path.expanduser('~/Jensen_response/analysis/%s_neutrality_locusLen11000_LD_0.05_0.95_averaged.txt' %model),'r')


    for line in inFile:
        items=line.strip().split('\t')
        position=items[0]
        R2=items[1]
        LD_dict[model]['position'].append(int(position))
        LD_dict[model]['R2'].append(float(R2))



# load LD estimated directly from the genome


LD_dict['DGRP']={}
LD_dict['DGRP']['position']=[]
LD_dict['DGRP']['R2']=[] 

inFile=open(os.path.expanduser('~/Jensen_response/analysis/DGRP_LD/alldata_averaged_sorted'),'r')

for line in inFile:
    items=line.strip('\n').split('\t')
    position=items[0]
    R2=items[1]
    LD_dict['DGRP']['position'].append(int(position))
    LD_dict['DGRP']['R2'].append(float(R2))
    
#colors=['green','blue', 'red','yellow', 'purple', 'orange', 'brown','pink', 'turquoise','magenta', 'lightgray','orchid','y']

# plot LD
for xlim_range in [100,10000]:
    plt.subplot(2, 2, plot_no)
    
    i=0
    for model in models:
        
        plt.semilogy(LD_dict[model]['position'], LD_dict[model]['R2'], label=model, color=colors[i])
        plt.xlim(0,xlim_range) 
        if xlim_range == 100:
            plt.ylim(0.1,1.05)
        else:
            plt.ylim(0.005,1.05)
        plt.xlabel('Distance between SNPs (bps)')
        plt.ylabel('log($R^2$)')
        i +=1

    plt.semilogy(LD_dict['DGRP']['position'], LD_dict['DGRP']['R2'], linestyle='dashed',color='black', label='DGRP')    
       
    plot_no +=1


#################
# print labels  #
#################
'''
print(plot_no)
plt.subplot(2,2,8)
# put some text in plot number 6
plt.text(0.045, 1.0, "DGRP data",  rotation=90.)
plt.text(0.113, 1.0, "Constant $Ne$=$10^6$",  rotation=90.)
plt.text(0.18, 1.0, "Constant $Ne$=$2.7*10^6$",  rotation=90.)
plt.text(0.248, 1.0, "Severe, short bottleneck",  rotation=90.)
plt.text(0.315, 1.0, "Shallow, long bottleneck",  rotation=90.)
plt.text(0.383, 1.0, "Admixture (mode), Garud $\it{et al.}$ 2015 implementation",  rotation=90.)
plt.text(0.45, 1.0, "Admixture + bottleneck (mode), Garud $\it{et al.}$ 2015 implementation",  rotation=90.)
plt.text(0.518, 1.0, "Admixture (mode), Duchen $\it{et al.}$ 2013",  rotation=90.)
plt.text(0.585, 1.0, "Admixture + bottleneck (mode), Duchen $\it{et al.}$ 2013",  rotation=90.)
plt.text(0.653, 1.0, "Admixture (posterior distribution), Duchen $\it{et al.}$ 2013",  rotation=90.)
plt.text(0.72, 1.0, "Admixture (posterior distribution), Harris $\it{et al.}$ 2018 implementation",  rotation=90.)
plt.text(0.788, 1.0, "Arguello $\it{et al.}$ 2019",  rotation=90.)
plt.text(0.856, 1.0, "Admixture, Eur=$7*10^5$, Nac=$1.11*10^6$",  rotation=90.)
plt.text(0.925, 1.0, "Admixture, Eur=$7*10^5$, Nac=$15.9*10^6$",  rotation=90.)
'''
#plt.axis('off')


plt.tight_layout()

plt.savefig(os.path.expanduser('~/Jensen_response/analysis/original_models_LD_2.png'), dpi=600)

plt.close()        



####################
# Plot H12         #
####################
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 5
fig_size[1] = 9

plot_no=1

# the diff. between H12_1 and H12_2 is that the second one has a zoomed in y-axis 

for summary_statistic in ['H12_1','H12_2']:
    plt.subplot(3,1, plot_no)

    if summary_statistic=='H12_1':
        summary_statistic='H12'
        bp = plt.boxplot(DGRP_H12, positions=[1], sym='.',patch_artist=True, labels=['DGRP'], widths=0.5)
        set_box_color(bp, 'black')

        col_i=0
        for model in models:
            bp = plt.boxplot(summary_statistics[model]['H12'], positions=[col_i+2], sym='.',patch_artist=True, widths=0.5, labels=[''])
            #patch.set_facecolor(colors[col_i])
            set_box_color(bp, colors[col_i]) 
            col_i+=1
        # dummy boxplot to get the colors to work for last one
        bp = plt.boxplot(summary_statistics[model]['H12'], positions=[col_i+10], sym='.',patch_artist=True, widths=0.5, labels=[''])
        #patch.set_facecolor('white')
        set_box_color(bp, 'white') 
        
        plt.xlim(0,len(models)+2)
        plt.plot( [1]*len(peaks),peaks, 'r.')

        plt.ylim(0,1)
        plt.ylabel('H12 (401 SNP windows)')


    elif summary_statistic=='H12_2':
        summary_statistic='H12'
        bp = plt.boxplot(DGRP_H12, positions=[1], sym='.',patch_artist=True, labels=['DGRP'], widths=0.5)
        set_box_color(bp, 'black')

        col_i=0
        for model in models:
            bp = plt.boxplot(summary_statistics[model]['H12'], positions=[col_i+2], sym='.',patch_artist=True, widths=0.5, labels=[''])
            #patch.set_facecolor(colors[col_i])
            set_box_color(bp, colors[col_i]) 
            col_i+=1
        # dummy boxplot to get the colors to work for last one
        bp = plt.boxplot(summary_statistics[model]['H12'], positions=[col_i+10], sym='.',patch_artist=True, widths=0.5, labels=[''])
        #patch.set_facecolor('white')
        set_box_color(bp, 'white') 
        
        plt.xlim(0,len(models)+2)
        plt.plot( [1]*len(peaks),peaks, 'r.')

        plt.ylim(0,0.1)
        plt.ylabel('H12 (401 SNP windows)')

    plt.axhline(y=data_medians['H12'],linestyle='dashed',color='black', label='Median value in data')
    plt.axhline(y=min(peaks),color='red', label='Lowest H12 peak value')


    plot_no +=1


plt.subplot(3,1,3)
# put some text in plot number 3
plt.text(-0.03, 1.0, "DGRP data")

plt.text(-0.03, 0.92, "Constant $Ne$=$10^6$")
plt.text(-0.03, 0.84, "Constant $Ne$=$2.7*10^6$")
plt.text(-0.03, 0.76, "Severe, short bottleneck")
plt.text(-0.03, 0.68, "Shallow, long bottleneck")
plt.text(-0.03, 0.60, "Admixture (mode), Garud $\it{et al.}$ 2015 implementation")
plt.text(-0.03, 0.52, "Admixture + bottleneck (mode), Garud $\it{et al.}$ 2015 implementation")
plt.text(-0.03, 0.44, "Admixture (mode), Duchen $\it{et al.}$ 2013")
plt.text(-0.03, 0.36, "Admixture + bottleneck (mode), Duchen $\it{et al.}$ 2013")
plt.text(-0.03, 0.28, "Admixture (posterior distribution), Duchen $\it{et al.}$ 2013")
plt.text(-0.03, 0.20, "Admixture (posterior distribution), Harris $\it{et al.}$ 2018 implementation")
plt.text(-0.03, 0.12, "Arguello $\it{et al.}$ 2019")
plt.text(-0.03, 0.04, "Admixture, Eur=$7*10^5$, Nac=$1.11*10^6$")
plt.text(-0.03, -0.04, "Admixture, Eur=$7*10^5$, Nac=$15.9*10^6$")

plt.axis('off')

plt.tight_layout()

plt.savefig(os.path.expanduser('~/Jensen_response/analysis/H12_plots.png'), dpi=600)

plt.close()        



















################
################
import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(4, 3)
gs.update(wspace=0.5, hspace=0.5)
#ax1 = plt.subplot(gs[0, :])
#ax2 = plt.subplot(gs[1, :-1])
#ax3 = plt.subplot(gs[1:, -1])
#ax4 = plt.subplot(gs[-1, 0])
#ax5 = plt.subplot(gs[-1, -2])

ax1 = plt.subplot(gs[0:2,0])
ax1.text(0.5, 0.5, str((2, 3, 1)),
           fontsize=18, ha='center')
ax2 = plt.subplot(gs[0:2,1])
ax3 = plt.subplot(gs[0:2,2])

ax4 = plt.subplot(gs[2:4,0])
ax5 = plt.subplot(gs[2:4,1])
ax6 = plt.subplot(gs[2:3,2])
ax7 = plt.subplot(gs[3:4,2])


col_i=0
summary_statistic='Pi'
for model in models:
    bp = ax1.boxplot(summary_statistics[model][summary_statistic], positions=[col_i+1], sym='.',patch_artist=True, widths=0.5, labels=[''])
    set_box_color(bp, colors[col_i]) 
    col_i+=1
    # dummy boxplot to get the colors to work for last one
    bp = ax1.boxplot(summary_statistics[model][summary_statistic], positions=[col_i+10], sym='.',patch_artist=True, widths=0.5, labels=[''])
    #patch.set_facecolor('white')
    set_box_color(bp, 'white') 


plt.savefig(os.path.expanduser('~/Jensen_response/analysis/tmpFig.png'), dpi=600)





'''
###############
###############




#######################################################################
# supplement: plot zoom in of H12 where y-axis ranges from 0 to 0.1   #
#######################################################################

matplotlib.rcParams.update({'font.size': 9})
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 9

plt.subplot(2,1, 1)


for summary_statistic in ['H12']:

        bp = plt.boxplot(DGRP_H12, positions=[1], sym='.',patch_artist=True, labels=['DGRP'], widths=0.5)
        set_box_color(bp, 'black')

        col_i=0
        for model in models:
            bp = plt.boxplot(summary_statistics[model][summary_statistic], positions=[col_i+2], sym='.',patch_artist=True, widths=0.5, labels=[''])
            #patch.set_facecolor(colors[col_i])
            set_box_color(bp, colors[col_i]) 
            col_i+=1
        # dummy boxplot to get the colors to work for last one
        bp = plt.boxplot(summary_statistics[model][summary_statistic], positions=[col_i+10], sym='.',patch_artist=True, widths=0.5, labels=[''])
        #patch.set_facecolor('white')
        set_box_color(bp, 'white') 
        

        plt.xlim(0,len(models)+2)

        plt.plot( [1]*len(peaks),peaks, 'r.')


        plt.axhline(y=data_medians[summary_statistic],linestyle='dashed',color='black', label='Median value in data')

        plt.ylim(0,.1)
        plt.ylabel('H12 (401 SNP windows)')



plt.subplot(2, 1, 2)
# put some text in plot number 2
plt.text(0.059, 1.1, "DGRP data",  rotation=90.)
plt.text(0.1256154, 1.1, "Constant $Ne$=$10^6$",  rotation=90.)
plt.text(0.1922308, 1.1, "Constant $Ne$=$2.7*10^6$",  rotation=90.)
plt.text(0.2588461, 1.1, "Severe, short bottleneck",  rotation=90.)
plt.text(0.3254615, 1.1, "Shallow, long bottleneck",  rotation=90.)
plt.text(0.3920769, 1.1, "Admixture (mode), Garud $\it{et al.}$ 2015 implementation",  rotation=90.)
plt.text(0.4586923, 1.1, "Admixture + bottleneck (mode), Garud $\it{et al.}$ 2015 implementation",  rotation=90.)
plt.text(0.5253077, 1.1, "Admixture (mode), Duchen $\it{et al.}$ 2013",  rotation=90.)
plt.text(0.591923, 1.1, "Admixture + bottleneck (mode), Duchen $\it{et al.}$ 2013",  rotation=90.)
plt.text(0.6585384, 1.1, "Admixture (posterior distribution), Duchen $\it{et al.}$ 2013",  rotation=90.)
plt.text(0.7251538, 1.1, "Admixture (posterior distribution), Harris $\it{et al.}$ 2018 implementation",  rotation=90.)
plt.text(0.7917692, 1.1, "Arguello $\it{et al.}$ 2019",  rotation=90.)
plt.text(0.8583846, 1.1, "Admixture, Eur=$7*10^5$, Nac=$1.11*10^6$",  rotation=90.)
plt.text(0.925, 1.1, "Admixture, Eur=$7*10^5$, Nac=$15.9*10^6$",  rotation=90.)

plt.axis('off')



plt.savefig(os.path.expanduser('~/Jensen_response/analysis/supplement_H12_zoomIn.png'), dpi=600)

plt.close()        


print('done with zoom')


########################
#                      #
## Plot simulations   ##
#                      #
########################

print 'Experiment 1: fixed population sizes'

#####################################################
# experiment 1: 
# fixed population sizes for Europe and North America
######################################################

# load the data

summary_statistics={} # store the summary statistics for all the models here. Key = [m_NA][m_NE][H12,S,Pi,avg_dist]

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

#########
# plot: # 
#########

colors=['#fed976','#feb24c', '#fd8d3c','#fc4e2a', '#e31a1c', '#b10026', '#bdd7e7','#6baed6', '#08519c','#8c510a','pink','#d9ef8b', '#66bd63']



fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 5


supplement_plot_no=3

for EUr in EUrs:

    plot_no=1
    for summary_statistic in ['Pi', 'S', 'H12']:
        plt.subplot(2, 3, plot_no)

        if summary_statistic=='H12':
            bp = plt.boxplot(DGRP_H12, positions=[1], sym='.',patch_artist=True, labels=['DGRP'], widths=0.5)
            set_box_color(bp, 'black')
            plt.plot( [1]*len(peaks),peaks, 'r.')

            col_i=0
            for NAm in NAms:
                bp = plt.boxplot(summary_statistics[NAm][EUr][summary_statistic], positions=[col_i+2], sym='.',patch_artist=True, widths=0.5, labels=[''])
                set_box_color(bp, colors[col_i]) 
                col_i+=1
            plt.labels=['1','1','1','61659' ,'1110000','15984500']

            #xtickNames = plt.setp(ax1, xticklabels=)
            plt.setp(['1','1','1','61659' ,'1110000','15984500'])

        else:
            col_i=0
            for NAm in NAms:
                bp = plt.boxplot(summary_statistics[NAm][EUr][summary_statistic], positions=[col_i+1], sym='.',patch_artist=True, widths=0.5, labels=[''])
                set_box_color(bp, colors[col_i]) 
                col_i+=1

        # dummy boxplot to get the colors to work for last one
        bp = plt.boxplot(summary_statistics[NAm][EUr][summary_statistic], positions=[col_i+10], sym='.',patch_artist=True, widths=0.5, labels=[''])
        set_box_color(bp, 'white') 
        plt.xlim(0,len(NAms)+2)

        # plot data median
        plt.axhline(y=data_medians[summary_statistic],linestyle='dashed',color='black', label='Median value in data')

        if summary_statistic=='H12':
            plt.ylim(0,0.3)
            plt.ylabel('H12 (401 SNP windows)')
        elif summary_statistic=='Pi':
            plt.ylim(0,0.03)
            plt.ylabel('Pi/bp')  
            plt.title('European $Ne$ = %s' %EUr)
        elif summary_statistic=='S':
            plt.ylim(0,0.1)
            plt.ylabel('S/bp')   

        plot_no +=1

    ##########################
    # Plot LD                # 
    ##########################

    for xlim_range in [100,10000]:
        plt.subplot(2, 3, plot_no)
    
        i=0
        for NAm in NAms:
            plt.semilogy(LD_dict[NAm][EUr]['position'], LD_dict[NAm][EUr]['R2'], color=colors[i])
            plt.xlim(0,xlim_range) 
            if xlim_range == 100:
                plt.ylim(0.1,1.05)
            else:
                plt.ylim(0.005,1.05)
            i +=1

        plt.xlabel('Distance between SNPs (bps)')
        plt.ylabel('log($R^2$)')        
        plt.semilogy(LD_dict['DGRP']['position'], LD_dict['DGRP']['R2'], linestyle='dashed',color='black', label='DGRP')    
        plot_no +=1


        
    plt.subplot(2, 3, plot_no)
    # put some text in plot number 6
    plt.text(0.128, 1.05, "DGRP data",  rotation=90.)
    plt.text(0.2704, 1.05, "2,500",  rotation=90.)
    plt.text(0.4128, 1.05, "61,659",  rotation=90.)
    plt.text(0.5552, 1.05, "1,110,000",  rotation=90.)
    plt.text(0.6976, 1.05, "15,984,500",  rotation=90.)
    plt.text(0.84, 1.05, "28,000,000",  rotation=90.)
    plt.text(0.3, 0.6, "North American $Ne$")

    plt.axis('off')

    plt.tight_layout()
    plt.savefig(os.path.expanduser('~/Jensen_response/analysis/Fig_S%s.png' %supplement_plot_no), dpi=600)
    
    supplement_plot_no +=1

    plt.close()



######################################
# experiment 2: Varying growth rates #
######################################

print 'Experiment 2: Varying growth rates'

# load the data

summary_statistics={} # store the summary statistics for all the models here. Key = [m_NA][m_NE][H12,S,Pi,avg_dist]

NAm_ancs=[2500, 61659]
NAm_curs=[1110000, 15984500, 28000000]
EUr_ancs=[16982, 67608]
EUr_curs=[700000, 2000000, 9550000]

for NAm_anc in NAm_ancs: 
    summary_statistics[NAm_anc]={}
    LD_dict[NAm_anc]={}
    for NAm_cur in NAm_curs:
        summary_statistics[NAm_anc][NAm_cur]={}
        LD_dict[NAm_anc][NAm_cur]={} 
        for EUr_anc in EUr_ancs:
            summary_statistics[NAm_anc][NAm_cur][EUr_anc]={}
            LD_dict[NAm_anc][NAm_cur][EUr_anc]={}
            for EUr_cur in EUr_curs:
                summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur]={}
                LD_dict[NAm_anc][NAm_cur][EUr_anc][EUr_cur]={}
 
                summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['H12']=[]
                summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['H2H1']=[]
                summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['S']=[]
                summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['Pi']=[]
                summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['avg_dist']=[]
                summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['TajD']=[]
                LD_dict[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['position']=[]
                LD_dict[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['R2']=[]


for NAm_anc in NAm_ancs: 
    for NAm_cur in NAm_curs:
        for EUr_anc in EUr_ancs:
            for EUr_cur in EUr_curs:
                # H12
                inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_varyGrowth_NAm_anc%s_NAm_cur%s_EUr_anc%s_EUr_cur%s_rho_%s_theta_%s_selection_%s_MS_snps.txt') %(NAm_anc, NAm_cur, EUr_anc, EUr_cur, rho_in, adaptive_theta, selection),'r')
        
                # iterate through the file
                # check if there are 17 lines
                # store H12 for this parameter configuration in a dictionary
                for line in inFile:
                    items=line.split('\t')
                    if len(items)==17:
                        H12=float(items[8])
                        H2H1=float(items[9])
                        summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['H12'].append(H12)
                        summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['H2H1'].append(H2H1)
                inFile.close()

                # Pi, S, TajD
                inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_varyGrowth_NAm_anc%s_NAm_cur%s_EUr_anc%s_EUr_cur%s_rho_%s_theta_%s_selection_%s_MS_Pi_S_TajD_perBp.txt') %(NAm_anc, NAm_cur, EUr_anc, EUr_cur,rho_in, adaptive_theta, selection),'r')
                for line in inFile:
                    items=line.strip('\n').split('\t')
                    if len(items)==4:
                        if isFloat(items[0]) and isFloat(items[1]) and isFloat(items[2]) and isFloat(items[3]):
                            S=float(items[0])
                            Pi=float(items[1])
                            avg_dist=float(items[2])
                            TajD=float(items[3])
                            summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['S'].append(S)
                            summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['Pi'].append(Pi)
                            summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['avg_dist'].append(avg_dist)
                            summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['TajD'].append(TajD)

                inFile.close()


                # LD
                inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_varyGrowth_NAm_anc%s_NAm_cur%s_EUr_anc%s_EUr_cur%s_rho_5e-9_theta_0_selection_False_MS_LD_0.05_0.95_averaged.txt' %(NAm_anc, NAm_cur, EUr_anc, EUr_cur)),'r')

                for line in inFile:
                    items=line.strip().split('\t')
                    position=items[0]
                    R2=items[1]
                    LD_dict[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['position'].append(int(position))
                    LD_dict[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['R2'].append(float(R2))

            
print('data is loaded...')

#########
# plot: # 
#########

#colors=['#fed976','#feb24c', '#fd8d3c','#fc4e2a', '#e31a1c', '#b10026', '#bdd7e7','#6baed6', '#08519c','#8c510a','pink','#d9ef8b', '#66bd63']

colors=['#fed976', '#fd8d3c', '#e31a1c', '#bdd7e7','#6baed6', '#08519c']

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 5


for EUr_anc in EUr_ancs:
    for EUr_cur in EUr_curs:
        plot_no=1
        for summary_statistic in ['Pi', 'S', 'H12']:
            plt.subplot(2, 3, plot_no)

            if summary_statistic=='H12':
                bp = plt.boxplot(DGRP_H12, positions=[1], sym='.',patch_artist=True, labels=['DGRP'], widths=0.5)
                set_box_color(bp, 'black')
                plt.plot( [1]*len(peaks),peaks, 'r.')

                col_i=0
                for NAm_anc in NAm_ancs:
                    for NAm_cur in NAm_curs:
                        bp = plt.boxplot(summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur][summary_statistic], positions=[col_i+2], sym='.',patch_artist=True, widths=0.5, labels=[''])
                        set_box_color(bp, colors[col_i]) 
                        col_i+=1

            else:
                col_i=0
                for NAm_anc in NAm_ancs:
                    for NAm_cur in NAm_curs:
                        bp = plt.boxplot(summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur][summary_statistic], positions=[col_i+1], sym='.',patch_artist=True, widths=0.5, labels=[''])
                        set_box_color(bp, colors[col_i]) 
                        col_i+=1

            # dummy boxplot to get the colors to work for last one
            bp = plt.boxplot(summary_statistics[NAm_anc][NAm_cur][EUr_anc][EUr_cur][summary_statistic], positions=[col_i+10], sym='.',patch_artist=True, widths=0.5, labels=[''])
            set_box_color(bp, 'white') 
            plt.xlim(0,len(NAm_curs)*2+2)

            # plot data median
            plt.axhline(y=data_medians[summary_statistic],linestyle='dashed',color='black', label='Median value in data')

            if summary_statistic=='H12':
                plt.ylim(0,0.3)
                plt.ylabel('H12 (401 SNP windows)')
            elif summary_statistic=='Pi':
                plt.ylim(0,0.03)
                plt.ylabel('Pi/bp')  
                plt.title('Eur_anc $Ne$ = %s, Eur_cur $Ne$ = %s' %(EUr_anc, EUr_cur) , fontsize=9)
            elif summary_statistic=='S':
                plt.ylim(0,0.1)
                plt.ylabel('S/bp')   

            plot_no +=1

        ##########################
        # Plot LD                # 
        ##########################

        for xlim_range in [100,10000]:
            plt.subplot(2, 3, plot_no)
    
            i=0
            for NAm_anc in NAm_ancs:
                for NAm_cur in NAm_curs:
                    plt.semilogy(LD_dict[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['position'], LD_dict[NAm_anc][NAm_cur][EUr_anc][EUr_cur]['R2'], color=colors[i])
                    plt.xlim(0,xlim_range) 
                    if xlim_range == 100:
                        plt.ylim(0.1,1.05)
                    else:
                        plt.ylim(0.005,1.05)
                    i +=1

            plt.xlabel('Distance between SNPs (bps)')
            plt.ylabel('log($R^2$)')        
            plt.semilogy(LD_dict['DGRP']['position'], LD_dict['DGRP']['R2'], linestyle='dashed',color='black', label='DGRP')    
            plot_no +=1



        plt.subplot(2, 3, plot_no)
        # put some text in plot number 6
        plt.text(0.11, 1.05, "DGRP data",  rotation=90.)
        plt.text(0.235, 1.05, "NAm_anc=2,500, NAm_cur=1,110,000",  rotation=90.)
        plt.text(0.36, 1.05, "NAm_anc=2,500, NAm_cur=15,984,500",  rotation=90.)
        plt.text(0.485, 1.05, "NAm_anc=2,500, NAm_cur=28,000,000",  rotation=90.)
        plt.text(0.61, 1.05, "NAm_anc=61,659, NAm_cur=1,110,000",  rotation=90.)
        plt.text(0.735, 1.05, "NAm_anc=61,659, NAm_cur=15,984,500",  rotation=90.)
        plt.text(0.86, 1.05, "NAm_anc=61,659, NAm_cur=28,000,000",  rotation=90.)
        #plt.text(0.3, 0.6, "Migration rate")

        plt.axis('off')
            
        plt.tight_layout()
        plt.savefig(os.path.expanduser('~/Jensen_response/analysis/Fig_S%s.png' %(supplement_plot_no)), dpi=600)

        supplement_plot_no+=1
        plt.close()


#####################################################
# experiment 3: 
# different proportions of admixture
######################################################

print 'Experiment 2: Different proportions of admixture'


# load the data

summary_statistics={} # store the summary statistics for all the models here. Key = [prop][H12,S,Pi,avg_dist]


proportions=[0, 0.1, 0.25, 0.4, 0.5, 0.75, 0.85, 0.9]


for prop in proportions: 
    summary_statistics[prop]={}
    LD_dict[prop]={}
    summary_statistics[prop]={}
    summary_statistics[prop]['H12']=[]
    summary_statistics[prop]['H2H1']=[]
    summary_statistics[prop]['S']=[]
    summary_statistics[prop]['Pi']=[]
    summary_statistics[prop]['avg_dist']=[]
    summary_statistics[prop]['TajD']=[]
    LD_dict[prop]={}
    LD_dict[prop]['position']=[]
    LD_dict[prop]['R2']=[]


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
            H2H1=float(items[9])
            summary_statistics[prop]['H12'].append(H12)
            summary_statistics[prop]['H2H1'].append(H2H1)
    inFile.close()

    # Pi, S, TajD
    inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_migration_diffProps%s_rho_%s_theta_%s_selection_%s_MS_Pi_S_TajD_perBp.txt') %(prop,rho_in, adaptive_theta, selection),'r')
    for line in inFile:
        items=line.strip('\n').split('\t')
        if len(items)==4:
            if isFloat(items[0]) and isFloat(items[1]) and isFloat(items[2]) and isFloat(items[3]):
                S=float(items[0])
                Pi=float(items[1])
                avg_dist=float(items[2])
                TajD=float(items[3])
                summary_statistics[prop]['S'].append(S)
                summary_statistics[prop]['Pi'].append(Pi)
                summary_statistics[prop]['avg_dist'].append(avg_dist)
                summary_statistics[prop]['TajD'].append(TajD)

    inFile.close()


    # LD
    inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_migration_diffProps%s_rho_%s_theta_%s_selection_%s_MS_LD_0.05_0.95_averaged.txt' %(prop,rho_in, adaptive_theta, selection)),'r')

    for line in inFile:
        items=line.strip().split('\t')
        position=items[0]
        R2=items[1]
        LD_dict[prop]['position'].append(int(position))
        LD_dict[prop]['R2'].append(float(R2))

            
print('data is loaded...')

#########
# plot: # 
#########

colors=['#fed976','#feb24c', '#fd8d3c','#fc4e2a', '#e31a1c', '#b10026', '#bdd7e7','#6baed6', '#08519c','#8c510a','pink','#d9ef8b', '#66bd63']



fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 5


plot_no=1
for summary_statistic in ['Pi', 'S', 'H12']:
    plt.subplot(2, 3, plot_no)

    if summary_statistic=='H12':
        bp = plt.boxplot(DGRP_H12, positions=[1], sym='.',patch_artist=True, labels=['DGRP'], widths=0.5)
        set_box_color(bp, 'black')
        plt.plot( [1]*len(peaks),peaks, 'r.')

        col_i=0
        for prop in proportions:
            bp = plt.boxplot(summary_statistics[prop][summary_statistic], positions=[col_i+2], sym='.',patch_artist=True, widths=0.5, labels=[''])
            set_box_color(bp, colors[col_i]) 
            col_i+=1

    else:
        print summary_statistic
        col_i=0
        for prop in proportions:
            bp = plt.boxplot(summary_statistics[prop][summary_statistic], positions=[col_i+1], sym='.',patch_artist=True, widths=0.5, labels=[''])
            set_box_color(bp, colors[col_i]) 
            col_i+=1

    # dummy boxplot to get the colors to work for last one
    bp = plt.boxplot(summary_statistics[prop][summary_statistic], positions=[col_i+10], sym='.',patch_artist=True, widths=0.5, labels=[''])
    set_box_color(bp, 'white') 
    plt.xlim(0,len(proportions)+2)

    # plot data median
    plt.axhline(y=data_medians[summary_statistic],linestyle='dashed',color='black', label='Median value in data')

    if summary_statistic=='H12':
        plt.ylim(0,0.3)
        plt.ylabel('H12 (401 SNP windows)')
    elif summary_statistic=='Pi':
        plt.ylim(0,0.03)
        plt.ylabel('Pi/bp')  
        plt.title('Admixture proportion')
    elif summary_statistic=='S':
        plt.ylim(0,0.1)
        plt.ylabel('S/bp')   

    plot_no +=1

##########################
# Plot LD                # 
##########################
    
for xlim_range in [100,10000]:
    plt.subplot(2, 3, plot_no)
    
    i=0
    for prop in proportions:
        plt.semilogy(LD_dict[prop]['position'], LD_dict[prop]['R2'], color=colors[i])
        plt.xlim(0,xlim_range) 
        if xlim_range == 100:
            plt.ylim(0.1,1.05)
        else:
            plt.ylim(0.005,1.05)
        i +=1

    plt.xlabel('Distance between SNPs (bps)')
    plt.ylabel('log($R^2$)')        
    plt.semilogy(LD_dict['DGRP']['position'], LD_dict['DGRP']['R2'], linestyle='dashed',color='black', label='DGRP')    
    plot_no +=1

plt.subplot(2, 3, plot_no)
# put some text in plot number 6
plt.text(0.084, 1.05, "DGRP data",  rotation=90.)
plt.text(0.184, 1.05, "0",  rotation=90.)
plt.text(0.284, 1.05, "0.1",  rotation=90.)
plt.text(0.384, 1.05, "0.25",  rotation=90.)
plt.text(0.484, 1.05, "0.4",  rotation=90.)
plt.text(0.584, 1.05, "0.5",  rotation=90.)
plt.text(0.684, 1.05, "0.75",  rotation=90.)
plt.text(0.784, 1.05, "0.85",  rotation=90.)
plt.text(0.884, 1.05, "0.9",  rotation=90.)
plt.text(0.3, 0.6, "Proportion admixture")

plt.axis('off')
    


plt.tight_layout()
plt.savefig(os.path.expanduser('~/Jensen_response/analysis/Fig_S%s.png' %(supplement_plot_no)), dpi=600)

supplement_plot_no +=1

plt.close()



#####################################################
# experiment 4: 
# Admixture + migration  
######################################################

print 'Experiment 3: Different proportions of admixture'


# load the data

summary_statistics={} # store the summary statistics for all the models here. Key = [prop][H12,S,Pi,avg_dist]


m_NAs=[0, 0.1, 0.25, 0.5, 0.75]


for m_NA in m_NAs: 
    summary_statistics[m_NA]={}
    LD_dict[m_NA]={}
    summary_statistics[m_NA]={}
    summary_statistics[m_NA]['H12']=[]
    summary_statistics[m_NA]['H2H1']=[]
    summary_statistics[m_NA]['S']=[]
    summary_statistics[m_NA]['Pi']=[]
    summary_statistics[m_NA]['avg_dist']=[]
    summary_statistics[m_NA]['TajD']=[]
    LD_dict[m_NA]={}
    LD_dict[m_NA]['position']=[]
    LD_dict[m_NA]['R2']=[]


for m_NA in m_NAs:   
    # H12
    inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_migration_mNA%s_mNE%s_rho_%s_theta_%s_selection_%s_snps.txt') %(m_NA, m_NA ,rho_in, adaptive_theta, selection),'r')
        
    # iterate through the file
    # check if there are 17 lines
    # store H12 for this parameter configuration in a dictionary
    for line in inFile:
        items=line.split('\t')
        if len(items)==17:
            H12=float(items[8])
            H2H1=float(items[9])
            summary_statistics[m_NA]['H12'].append(H12)
            summary_statistics[m_NA]['H2H1'].append(H2H1)
    inFile.close()

    # Pi, S, TajD
    inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_migration_mNA%s_mNE%s_rho_%s_theta_%s_selection_%s_Pi_S_TajD_perBp.txt') %(m_NA, m_NA,rho_in, adaptive_theta, selection),'r')
    for line in inFile:
        items=line.strip('\n').split('\t')
        if len(items)==4:
            if isFloat(items[0]) and isFloat(items[1]) and isFloat(items[2]) and isFloat(items[3]):
                S=float(items[0])
                Pi=float(items[1])
                avg_dist=float(items[2])
                TajD=float(items[3])
                summary_statistics[m_NA]['S'].append(S)
                summary_statistics[m_NA]['Pi'].append(Pi)
                summary_statistics[m_NA]['avg_dist'].append(avg_dist)
                summary_statistics[m_NA]['TajD'].append(TajD)

    inFile.close()


    # LD
    inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_migration_mNA%s_mNE%s_rho_%s_theta_%s_selection_%s_LD_0.05_0.95_averaged.txt' %(m_NA, m_NA,rho_in, adaptive_theta, selection)),'r')

    for line in inFile:
        items=line.strip().split('\t')
        position=items[0]
        R2=items[1]
        LD_dict[m_NA]['position'].append(int(position))
        LD_dict[m_NA]['R2'].append(float(R2))

            
print('data is loaded...')

#########
# plot: # 
#########

colors=['#fed976','#feb24c', '#fd8d3c','#fc4e2a', '#e31a1c', '#b10026', '#bdd7e7','#6baed6', '#08519c','#8c510a','pink','#d9ef8b', '#66bd63']



fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 5


plot_no=1
for summary_statistic in ['Pi', 'S', 'H12']:
    plt.subplot(2, 3, plot_no)

    if summary_statistic=='H12':
        bp = plt.boxplot(DGRP_H12, positions=[1], sym='.',patch_artist=True, labels=['DGRP'], widths=0.5)
        set_box_color(bp, 'black')
        plt.plot( [1]*len(peaks),peaks, 'r.')

        col_i=0
        for m_NA in m_NAs:
            bp = plt.boxplot(summary_statistics[m_NA][summary_statistic], positions=[col_i+2], sym='.',patch_artist=True, widths=0.5, labels=[''])
            set_box_color(bp, colors[col_i]) 
            col_i+=1

    else:
        print summary_statistic
        col_i=0
        for m_NA in m_NAs:
            bp = plt.boxplot(summary_statistics[m_NA][summary_statistic], positions=[col_i+1], sym='.',patch_artist=True, widths=0.5, labels=[''])
            set_box_color(bp, colors[col_i]) 
            col_i+=1

    # dummy boxplot to get the colors to work for last one
    bp = plt.boxplot(summary_statistics[m_NA][summary_statistic], positions=[col_i+10], sym='.',patch_artist=True, widths=0.5, labels=[''])
    set_box_color(bp, 'white') 
    plt.xlim(0,len(m_NAs)+2)

    # plot data median
    plt.axhline(y=data_medians[summary_statistic],linestyle='dashed',color='black', label='Median value in data')

    if summary_statistic=='H12':
        plt.ylim(0,0.3)
        plt.ylabel('H12 (401 SNP windows)')
    elif summary_statistic=='Pi':
        plt.ylim(0,0.03)
        plt.ylabel('Pi/bp')  
        plt.title('Admixture with migration')
    elif summary_statistic=='S':
        plt.ylim(0,0.1)
        plt.ylabel('S/bp')   

    plot_no +=1

##########################
# Plot LD                # 
##########################
    
for xlim_range in [100,10000]:
    plt.subplot(2, 3, plot_no)
    
    i=0
    for m_NA in m_NAs:
        plt.semilogy(LD_dict[m_NA]['position'], LD_dict[m_NA]['R2'], color=colors[i])
        plt.xlim(0,xlim_range) 
        if xlim_range == 100:
            plt.ylim(0.1,1.05)
        else:
            plt.ylim(0.005,1.05)
        i +=1

    plt.xlabel('Distance between SNPs (bps)')
    plt.ylabel('log($R^2$)')        
    plt.semilogy(LD_dict['DGRP']['position'], LD_dict['DGRP']['R2'], linestyle='dashed',color='black', label='DGRP')    
    plot_no +=1


plt.subplot(2, 3, plot_no)
# put some text in plot number 6
plt.text(0.125, 1.05, "DGRP data",  rotation=90.)
plt.text(0.268, 1.05, "0",  rotation=90.)
plt.text(0.411, 1.05, "0.1",  rotation=90.)
plt.text(0.554, 1.05, "0.25",  rotation=90.)
plt.text(0.697, 1.05, "0.5",  rotation=90.)
plt.text(0.84, 1.05, "0.75",  rotation=90.)
plt.text(0.3, 0.6, "Migration rate")

plt.axis('off')


plt.tight_layout()
plt.savefig(os.path.expanduser('~/Jensen_response/analysis/Fig_S%s.png' %(supplement_plot_no)), dpi=600)

supplement_plot_no+=1

plt.close()


'''
