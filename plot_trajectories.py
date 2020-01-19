import matplotlib
matplotlib.use('Agg')
import numpy
import sys
from scipy.stats import rv_discrete
import os.path
import numpy as np

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


rho_in='5e-9'
selection='True'

adaptive_thetas=[0.01,10]
NAms=['61659', '1110000']
EUrs=['500000']

trajectories={}

for NAm in NAms: 
    trajectories[NAm]={}
    for EUr in EUrs: 
        trajectories[NAm][EUr]={}
        for adaptive_theta in adaptive_thetas:
            trajectories[NAm][EUr][adaptive_theta]={}
            for iteration in range(1,101):
                trajectories[NAm][EUr][adaptive_theta][iteration]={'generation':[],'frequency':[]}
        
                inFile=open(os.path.expanduser('~/Jensen_response/analysis/Admixture_mode_fixedPopSize_NAm%s_EUr%s_rho_%s_theta_%s_selection_%s_MS_trajectory_randomstart_%s.txt') %(NAm, EUr,rho_in, adaptive_theta, selection, iteration),'r')
                
                inFile.readline() # skip blank
                for line in inFile:
                    
                    items=line.strip('\n').split('\t')
                    generation=float(items[0])
                    frequency=float(items[6])
                    trajectories[NAm][EUr][adaptive_theta][iteration]['generation'].append(generation)
                    trajectories[NAm][EUr][adaptive_theta][iteration]['frequency'].append(frequency)
                    

                inFile.close()



# make a vector of the onset of selection for the four scenarios

threshold=0.0001

for_boxplot={}
for NAm in NAms:
    for_boxplot[NAm]={}
    for adaptive_theta in adaptive_thetas:
        for_boxplot[NAm][adaptive_theta]=[]
        for iteration in range(1,101):
            idx=0
            found=False
            for freq in trajectories[NAm][EUr][adaptive_theta][iteration]['frequency']:
                if freq>=threshold and found==False:
                    for_boxplot[NAm][adaptive_theta].append(trajectories[NAm][EUr][adaptive_theta][iteration]['generation'][idx])
                    found=True
                idx+=1
                    



##########


print('plotting...')

# make boxplots of the distributions
matplotlib.rcParams.update({'font.size': 8})
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 6

EUr='500000'

plot_no=1
for adaptive_theta in adaptive_thetas:
    plt.subplot(1, 2, plot_no)    

    for iteration in range(1,101):
        NAm='61659'
        plt.plot(trajectories[NAm][EUr][adaptive_theta][iteration]['generation'], trajectories[NAm][EUr][adaptive_theta][iteration]['frequency'], 'r')
        NAm='1110000'
        plt.plot(trajectories[NAm][EUr][adaptive_theta][iteration]['generation'], trajectories[NAm][EUr][adaptive_theta][iteration]['frequency'],'b')
        
        red_patch = mpatches.Patch(color='red', label='NAm=61659')
        blue_patch = mpatches.Patch(color='blue', label='NAm=1110000')

        plt.legend(handles=[red_patch, blue_patch])
        plt.xlabel('generation')
        plt.ylabel('frequency')
        #plt.ylim(0,.0005)
        #plt.xlim(0.00002, 0.00004)
        plt.title('adaptive theta:' + str(adaptive_theta))
    plot_no+=1

        
plt.tight_layout()
    
plt.savefig(os.path.expanduser('~/Jensen_response/analysis/admixture_fixedPopSize_trajectories.png'), dpi=600)

plt.close()



###################

print('plotting...')

# make boxplots of the distributions
matplotlib.rcParams.update({'font.size': 8})
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 6

colors=['red','blue','red','blue']

bp=plt.boxplot([for_boxplot['61659'][0.01], for_boxplot['1110000'][0.01], for_boxplot['61659'][10], for_boxplot['1110000'][10]], labels=['NA=61659, thetaA=0.01', 'NA=1110000, thetaA=0.01', 'NA=61659, thetaA=10', 'NA=1110000, thetaA=10'],  patch_artist=True)

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

plt.title('generation when f>=%s' %threshold)
plt.ylabel('generation')
        
plt.tight_layout()
    
plt.savefig(os.path.expanduser('~/Jensen_response/analysis/admixture_fixedPopSize_trajectories_boxplots.png'), dpi=600)

plt.close()



