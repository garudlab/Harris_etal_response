# November 13, 2018
# Plot for response to Harris et al. 

#####################################################
# BF heat map                                       #
# Red and Grey scale, truncated x axis, 3 panels    #
# Panel 1: Original plot, SNP based                 #
# Panel 2: More comprehensive dem model, SNP based  #
# Panel 3: More comprehensive dem model, bp based   #
#####################################################

setwd("/Users/nanditagarud/Documents/Jensen_response/analysis")

peaks=read.table('peaks_082712_400_50_2_2_0_pan_509_top50_rhoForCalling.txt')

pdf("fig8_BFheatmap_panels.pdf",width=5.69, height=2.43, title = "Figure 2",paper="special")
#width=6.83,height=4.86
#width=4.86,height=2.43

s=as.vector(c(1,5))
split.screen(s)

change=0.06

mb1=0.1 
ml1=0.4
mt1=0 
mr1=0 

mb2=0.1 
ml2=0.35 
mt2=0 
mr2=0.05

mb3=0.1 
ml3=0.3	 
mt3=0 
mr3=0.1 

mb4=0.1 
ml4=0.25	 
mt4=0 
mr4=0.15 

mb5=0.1 
ml5=0.2	 
mt5=0 
mr5=0.2 

mb6=0.1 
ml6=0.15	 
mt6=0 
mr6=0.25 

tcl=0.2

pdjy=2.0
pdjx=-2.5
  
cex1=0.7
cex2=0.7
cex3=1

par(cex.axis=cex1) #The magnification to be used for axis annotation relative to the current setting of cex.
par(cex.lab=cex1) #The magnification to be used for x and y labels relative to the current setting of cex.
par(omi=c(0.35,0.1,0.35,0.1)) #A vector of the form c(bottom, left, top, right) giving the size of the outer margins in inches.

shift=0.02


screen(1)

par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
par(xaxs = "i")
par(yaxs = "i")

d=read.table("HS_vs_SS_likelihood_random509_admixture_sorted.txt")
d2=(d[,4])/(d[,3])
d3=matrix(d2,c(40,40))
d4=d3[,1:12]


colorMatrix=as.matrix(d4)


for (i in 1:length(d4[,1])){
	for (j in 1:length(d4[1,])){
		

		if (d3[i,j]=='NaN'){colorMatrix[i,j]="white"} else 
		{if (d3[i,j]=='Inf'){colorMatrix[i,j]="black"} else 
		{if (d3[i,j]>=0 && d3[i,j]<1){colorMatrix[i,j]="red"} else
		{if (d3[i,j]>=1 && d3[i,j]<10){colorMatrix[i,j]="gray65"} else
		{if (d3[i,j]>=10 && d3[i,j]<30){colorMatrix[i,j]="gray55"} else
		{if (d3[i,j]>=30 && d3[i,j]<100){colorMatrix[i,j]="gray35"} else
		{if (d3[i,j]>=100 && d3[i,j]<1000000){colorMatrix[i,j]="black"} else
			{colorMatrix[i,j]=3}}}}}}}
		
	}	
}


plot(c(0, 12), c(0, 40), type= "n", xaxs="i",yaxs = "i",axes=F, main='', xlim=c(0,13), ylim=c(0,41))
h=0
w=0.35
for (i in 1:length(d4[,1])){
	height=c(h-0.5, h+0.5)
	for (j in 1:length(d4[1,])){
		rect(j-0.5, height[1], j+0.5, height[2],density=-1,col=colorMatrix[i,j],lwd=0., border=colorMatrix[i,j])

	}
	h=h+1	
} 

axis(1,at=c(0,12),labels=F,tcl=0)
axis(1, at = c(0,6,12), labels=c("0","0.15","0.30"), padj=pdjx,tcl=-tcl)
axis(2,at=c(0,40),labels=F,tcl=0)
axis(2, at = c(0,10,20,30,40), labels=c("0","0.25","0.50", "0.75", "1" ),padj=pdjy,tcl=-tcl)


# overlay the observed peak's values
points(40*peaks[,28],40*peaks[,20],col='gold',pch=20, cex=0.6)
mtext(expression(H2~"/"~H1),side=2,cex=cex2,padj=-3)
mtext(expression(A),at=-0.5,padj=-0.5,cex=cex3)

# Draw a legend on the upper right:
#rect(5-0.5, 40-0.5, 5+0.5, 40+0.5,density=-1,col='black',lwd=0., border="white")
#rect(5-0.5, 38.5-0.5, 5+0.5, 38.5+0.5,density=-1,col='gray35',lwd=0., border="white")
#rect(5-0.5, 37-0.5, 5+0.5, 37+0.5,density=-1,col='gray55',lwd=0., border="white")
#rect(5-0.5, 35.5-0.5, 5+0.5, 35.5+0.5,density=-1,col='gray65',lwd=0., border="white")
#rect(5-0.5, 34-0.5, 5+0.5, 34+0.5,density=-1,col='red',lwd=0., border="white")

#mtext(expression(BF>100), at=c(8.1), line=-.5, cex=0.35)
#mtext(expression(100>~BF>=30)  , at=c(9.25), line=-0.8, cex=0.35)
#mtext(expression(30>~BF>=10)  , at=c(9), line=-1.1, cex=0.35)
#mtext(expression(10>~BF>=1)  , at=c(8.75), line=-1.4, cex=0.35)
#mtext(expression(1>~BF>=0)  , at=c(8.5), line=-1.7, cex=0.35)

#rect(4,33,13,41, density=-1,border='black', lwd=0.5)

screen(2)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
par(xaxs = "i")
par(yaxs = "i")

d=read.table("Admixture_posteriorDistn_95CI_rho_5e-9_theta_0.01_theta_10_selection_True_snps_BFs_numSim_16000.txt")
d2=(d[,4])/(d[,3])
d3=matrix(d2,c(40,40))
d4=d3[,1:12]


colorMatrix=as.matrix(d4)


for (i in 1:length(d4[,1])){
	for (j in 1:length(d4[1,])){
		

		if (d3[i,j]=='NaN'){colorMatrix[i,j]="white"} else 
		{if (d3[i,j]=='Inf'){colorMatrix[i,j]="black"} else 
		{if (d3[i,j]>=0 && d3[i,j]<1){colorMatrix[i,j]="red"} else
		{if (d3[i,j]>=1 && d3[i,j]<10){colorMatrix[i,j]="gray65"} else
		{if (d3[i,j]>=10 && d3[i,j]<30){colorMatrix[i,j]="gray55"} else
		{if (d3[i,j]>=30 && d3[i,j]<100){colorMatrix[i,j]="gray35"} else
		{if (d3[i,j]>=100 && d3[i,j]<1000000){colorMatrix[i,j]="black"} else
			{colorMatrix[i,j]=3}}}}}}}
		
	}	
}


plot(c(0, 12), c(0, 40), type= "n", xaxs="i",yaxs = "i",axes=F, main='', xlim=c(0,13))
h=0
w=0.35
for (i in 1:length(d4[,1])){
	height=c(h-0.5, h+0.5)
	for (j in 1:length(d4[1,])){
		rect(j-0.5, height[1], j+0.5, height[2],density=-1,col=colorMatrix[i,j],lwd=0., border=colorMatrix[i,j])

	}
	h=h+1	
}

axis(1,at=c(0,12),labels=F,tcl=0)
axis(1, at = c(0,6,12), labels=c("0","0.15","0.30"), padj=pdjx,tcl=-tcl)
axis(2,at=c(0,40),labels=F,tcl=0)
axis(2, at = c(0,10,20,30,40), labels=c("0","0.25","0.50", "0.75", "1" ),padj=pdjy,tcl=-tcl)

# overlay the observed peak's values
points(40*peaks[,28],40*peaks[,20],col='gold',pch=20, cex=0.6)
mtext(expression(B),at=-0.5,padj=-0.5,cex=cex3)


screen(3)

peaks=read.table('top50peaks_fixed_bp.txt')


par(mai=c(mb3,ml3,mt3,mr3)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
par(xaxs = "i")
par(yaxs = "i")

d=read.table("Admixture_posteriorDistn_95CI_rho_5e-9_theta_0.01_theta_10_selection_True_bps_BFs_numSim_16000.txt")
d2=(d[,4])/(d[,3])
d3=matrix(d2,c(40,40))
d4=d3[,1:12]


colorMatrix=as.matrix(d4)


for (i in 1:length(d4[,1])){
	for (j in 1:length(d4[1,])){
		

		if (d3[i,j]=='NaN'){colorMatrix[i,j]="white"} else 
		{if (d3[i,j]=='Inf'){colorMatrix[i,j]="black"} else 
		{if (d3[i,j]>=0 && d3[i,j]<1){colorMatrix[i,j]="red"} else
		{if (d3[i,j]>=1 && d3[i,j]<10){colorMatrix[i,j]="gray65"} else
		{if (d3[i,j]>=10 && d3[i,j]<30){colorMatrix[i,j]="gray55"} else
		{if (d3[i,j]>=30 && d3[i,j]<100){colorMatrix[i,j]="gray35"} else
		{if (d3[i,j]>=100 && d3[i,j]<1000000){colorMatrix[i,j]="black"} else
			{colorMatrix[i,j]=3}}}}}}}
		
	}	
}


plot(c(0, 12), c(0, 40), type= "n", xaxs="i",yaxs = "i",axes=F, main='', xlim=c(0,13))
h=0
w=0.35
for (i in 1:length(d4[,1])){
	height=c(h-0.5, h+0.5)
	for (j in 1:length(d4[1,])){
		rect(j-0.5, height[1], j+0.5, height[2],density=-1,col=colorMatrix[i,j],lwd=0., border=colorMatrix[i,j])

	}
	h=h+1	
}

axis(1,at=c(0,12),labels=F,tcl=0)
axis(1, at = c(0,6,12), labels=c("0","0.15","0.30"), padj=pdjx,tcl=-tcl)
axis(2,at=c(0,40),labels=F,tcl=0)
axis(2, at = c(0,10,20,30,40), labels=c("0","0.25","0.50", "0.75", "1" ),padj=pdjy,tcl=-tcl)

# overlay the observed peak's values
points(40*peaks[,9],40*peaks[,10],col='gold',pch=20, cex=0.6)

mtext(expression(H12),side=1,cex=cex2,padj=2.25, at=16)

mtext(expression(C),at=-0.5,padj=-0.5,cex=cex3)



screen(4)
peaks=read.table('peaks_082712_400_50_2_2_0_pan_509_top50_rhoForCalling.txt')

par(mai=c(mb4,ml4,mt4,mr4)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
par(xaxs = "i")
par(yaxs = "i")

d=read.table("Admixture_uniformPrior_rho_5e-9_theta_0.01_theta_10_selection_True_snps_BFs_numSim_16000.txt")
# try zeroing out entries with really low sim counts
a=d[,3]+d[,4]
#b=a<100
#d[b,3]=0
#d[b,4]=0

d2=(d[,4])/(d[,3])
d3=matrix(d2,c(40,40))
d4=d3[,1:12]


colorMatrix=as.matrix(d4)


for (i in 1:length(d4[,1])){
	for (j in 1:length(d4[1,])){
		

		if (d3[i,j]=='NaN'){colorMatrix[i,j]="white"} else 
		{if (d3[i,j]=='Inf'){colorMatrix[i,j]="black"} else 
		{if (d3[i,j]>=0 && d3[i,j]<1){colorMatrix[i,j]="red"} else
		{if (d3[i,j]>=1 && d3[i,j]<10){colorMatrix[i,j]="gray65"} else
		{if (d3[i,j]>=10 && d3[i,j]<30){colorMatrix[i,j]="gray55"} else
		{if (d3[i,j]>=30 && d3[i,j]<100){colorMatrix[i,j]="gray35"} else
		{if (d3[i,j]>=100 && d3[i,j]<1000000){colorMatrix[i,j]="black"} else
			{colorMatrix[i,j]=3}}}}}}}
		
	}	
}


plot(c(0, 12), c(0, 40), type= "n", xaxs="i",yaxs = "i",axes=F, main='', xlim=c(0,13))
h=0
w=0.35
for (i in 1:length(d4[,1])){
	height=c(h-0.5, h+0.5)
	for (j in 1:length(d4[1,])){
		rect(j-0.5, height[1], j+0.5, height[2],density=-1,col=colorMatrix[i,j],lwd=0., border=colorMatrix[i,j])

	}
	h=h+1	
}

axis(1,at=c(0,12),labels=F,tcl=0)
axis(1, at = c(0,6,12), labels=c("0","0.15","0.30"), padj=pdjx,tcl=-tcl)
axis(2,at=c(0,40),labels=F,tcl=0)
axis(2, at = c(0,10,20,30,40), labels=c("0","0.25","0.50", "0.75", "1" ),padj=pdjy,tcl=-tcl)

# overlay the observed peak's values
points(40*peaks[,28],40*peaks[,20],col='gold',pch=20, cex=0.6)
mtext(expression(D),at=-0.5,padj=-0.5,cex=cex3)



screen(5)
peaks=read.table('top50peaks_fixed_bp.txt')

par(mai=c(mb4,ml4,mt4,mr4)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
par(xaxs = "i")
par(yaxs = "i")

d=read.table("Admixture_uniformPrior_rho_5e-9_theta_0.01_theta_10_selection_True_bps_BFs_numSim_16000.txt")
# try zeroing out entries with really low sim counts
a=d[,3]+d[,4]
#b=a<100
#d[b,3]=0
#d[b,4]=0

d2=(d[,4])/(d[,3])
d3=matrix(d2,c(40,40))
d4=d3[,1:12]


colorMatrix=as.matrix(d4)


for (i in 1:length(d4[,1])){
	for (j in 1:length(d4[1,])){
		

		if (d3[i,j]=='NaN'){colorMatrix[i,j]="white"} else 
		{if (d3[i,j]=='Inf'){colorMatrix[i,j]="black"} else 
		{if (d3[i,j]>=0 && d3[i,j]<1){colorMatrix[i,j]="red"} else
		{if (d3[i,j]>=1 && d3[i,j]<10){colorMatrix[i,j]="gray65"} else
		{if (d3[i,j]>=10 && d3[i,j]<30){colorMatrix[i,j]="gray55"} else
		{if (d3[i,j]>=30 && d3[i,j]<100){colorMatrix[i,j]="gray35"} else
		{if (d3[i,j]>=100 && d3[i,j]<1000000){colorMatrix[i,j]="black"} else
			{colorMatrix[i,j]=3}}}}}}}
		
	}	
}


plot(c(0, 12), c(0, 40), type= "n", xaxs="i",yaxs = "i",axes=F, main='', xlim=c(0,13))
h=0
w=0.35
for (i in 1:length(d4[,1])){
	height=c(h-0.5, h+0.5)
	for (j in 1:length(d4[1,])){
		rect(j-0.5, height[1], j+0.5, height[2],density=-1,col=colorMatrix[i,j],lwd=0., border=colorMatrix[i,j])

	}
	h=h+1	
}

axis(1,at=c(0,12),labels=F,tcl=0)
axis(1, at = c(0,6,12), labels=c("0","0.15","0.30"), padj=pdjx,tcl=-tcl)
axis(2,at=c(0,40),labels=F,tcl=0)
axis(2, at = c(0,10,20,30,40), labels=c("0","0.25","0.50", "0.75", "1" ),padj=pdjy,tcl=-tcl)

# overlay the observed peak's values
points(40*peaks[,9],40*peaks[,10],col='gold',pch=20, cex=0.6)
mtext(expression(E),at=-0.5,padj=-0.5,cex=cex3)



close.screen(all=TRUE)

dev.off() 




############################################################################################
# Plot the distribution of H12 values observed under different demographic models and data #
############################################################################################


###############
# SNP windows #
###############


setwd("/Users/nanditagarud/Documents/Ace/neutralDemography")
a=read.table('constNe106_neutrality_400SNPs_MS_snps.txt')
#a=read.table('constantNe106_neutral.txt')
b=read.table('dadi1_recoded_neutrality_MS_snps.txt')
#b=read.table('neutralSims_dadi1_509_new.txt')
c=read.table('dadi2_neutrality_recoded_MS_snps.txt')
#c=read.table('neutralSims_dadi2_509_new.txt')
##d=read.table('admix_noSelection_509.txt')
##e=read.table('admix_bot_noSelection_5*10e-9.txt_2') # admixture + bottleneck # this has old growth rate
f=read.table('constNe2.2e6_neutrality_400SNPs_MS_snps.txt')
#f=read.table('constantNe2.7_neutral.txt') # constantNe=2.7x10^6

setwd("/Users/nanditagarud/Documents/Jensen_response/analysis/")

##g=read.table('Admixture_posteriorDistn_95CI_rho_5e-9_theta_0_selection_False_snps.txt') # add the admixture with paramteres drawn from 95% CIs
d=read.table('Admixture_mode_corrected_rho_5e-9_theta_0_selection_False_MS_snps.txt')
#d=read.table('Admixture_mode_corrected_rho_5e-9_theta_0_selection_False_snps.txt')
e=read.table('admixture_bottleneck_mode_corrected_neutrality_MS.txt_snps.txt')
#e=read.table('admixture_bottleneck_mode_corrected_neutrality.txt_snps.txt')
i=read.table('Admixture_posteriorDistn_95CI_corrected_rho_5e-9_theta_0_selection_False_MS_snps.txt')
#i=read.table('Admixture_posteriorDistn_95CI_corrected_rho_5e-9_theta_0_selection_False_snps.txt') # add the admixture with paramteres drawn from 95% CIs correcting the growth rate


setwd("/Users/nanditagarud/Documents/Ace/")

chr3R_509 = read.table('clusters_082712_400_50_2_2_0_chr3R_withRho_threshold.txt')
chr3L_509 = read.table('clusters_082712_400_50_2_2_0_chr3L_withRho_threshold.txt')
chr2L_509 = read.table('clusters_082712_400_50_2_2_0_chr2L_withRho_threshold.txt')
chr2R_509 = read.table('clusters_082712_400_50_2_2_0_chr2R_withRho_threshold.txt')
h=c(chr3R_509[,28],chr3L_509[,28],chr2R_509[,28],chr2L_509[,28])

peaks=read.table('peaks_082712_400_50_2_2_0_pan_509_top50.txt')

setwd("/Users/nanditagarud/Documents/Jensen_response/analysis")


# plot:

pdf("Figure6A.pdf",width=3.27,height=3, title = "Figure S4",paper="special")

tcl=0.2

pdjy=2.5
pdjx=-3.0
  
cex1=0.5
cex2=0.6
cex3=1

par(cex.axis=cex1) #The magnification to be used for axis annotation relative to the current setting of cex.
par(cex.lab=cex1) #The magnification to be used for x and y labels relative to the current setting of cex.
par(omi=c(0,0,0,0)) #A vector of the form c(bottom, left, top, right) giving the size of the outer margins in inches.
par(mai=c(0.5,0.5,0.5,0.5)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
shift=0.02

boxplot(h,d[1:150000,9], e[1:150000,9],i[1:150000,9],a[1:150000,26],f[1:150000,26], b[1:150000,26],c[1:150000,26], ylim=c(0,1), col=c('black','orange','purple','brown','skyblue','forestgreen','blue','red','gold'), pch=20, cex=0.4, border=c('black','orange','purple','brown','skyblue','forestgreen','blue','red','gold'), xaxs="i",yaxs = "i", xlab="",ylab="", axes=F )
points(rep(1,50), peaks[,28], col='red', pch=20, cex=0.7)

axis(1,at=c(0,8),labels=F,tcl=0)

axis(2,at=c(0,1),labels=F,tcl=0)
axis(2,at=seq(0,1,by=0.05),labels=seq(0,1,by=0.05),padj=pdjy,tcl=-tcl)
axis(3,at=c(0,8),labels=F,tcl=0)
axis(4,at=c(0,1),labels=F,tcl=0)

mtext(expression("H12"),side=2,cex=cex2,padj=-2.25)

dev.off() 

############
# Zoom in 
############

boxplot(h, d[1:150000,9], e[1:150000,9], i[1:150000,9], ylim=c(0,0.1), col=c('grey','orange','purple','brown'), ylab='H12', names=c('data','admix_mode','admix_bot_mode','admix_95CI'), main='400 SNP windows, Zoom in')
points(rep(1,50), peaks[,28], col='red', pch=20, cex=0.1.3)


###############
# 10kb windows #
###############

setwd("/Users/nanditagarud/Documents/Jensen_response/analysis")

#admix_mode=read.table('Admixture_mode_rho_5e-9_theta_0_selection_False_bps.txt')
#admix_mode_corr=read.table('Admixture_mode_corrected_rho_5e-9_theta_0_selection_False_bps.txt')
admix_mode_corr=read.table('Admixture_mode_corrected_rho_5e-9_theta_0_selection_False_MS_bps.txt')
#admix_posterior=read.table('Admixture_posteriorDistn_95CI_rho_5e-9_theta_0_selection_False_bps.txt')
#admix_posterior_corr=read.table('Admixture_posteriorDistn_95CI_corrected_rho_5e-9_theta_0_selection_False_bps.txt')
admix_posterior_corr=read.table('Admixture_posteriorDistn_95CI_corrected_rho_5e-9_theta_0_selection_False_MS_bps.txt')

#admix_uniform=read.table('Admixture_uniformPrior_rho_5e-9_theta_0_selection_False_bps.txt')
#admix_bot_mode=read.table('admixture_bottleneck_mode_neutrality.txt_bps.txt') # admixture + bottleneck
#admix_bot_mode_corr=read.table('admixture_bottleneck_mode_corrected_neutrality.txt_bps.txt') # admixture + bottleneck
admix_bot_mode_corr=read.table('admixture_bottleneck_mode_corrected_neutrality_MS.txt_bps.txt') # admixture + bottleneck

#const106=read.table('constNe106_neutrality_bps.txt') # constantNe=2.7x10^6
#const22106=read.table('constNe2.2e6_neutrality_bps.txt') # constantNe=2.7x10^6
#dadi1=read.table('dadi1_neutrality_bps.txt') # constantNe=2.7x10^6
#dadi2=read.table('dadi2_neutrality_bps.txt') # constantNe=2.7x10^6

const106=read.table('constNe106_neutrality_MS_bps.txt') # constantNe=2.7x10^6
const22106=read.table('constNe2.2e6_neutrality_MS_bps.txt') # constantNe=2.7x10^6
dadi1=read.table('dadi1_neutrality_bps.txt') # constantNe=2.7x10^6
dadi2=read.table('dadi2_neutrality_recoded_MS_bps.txt') # constantNe=2.7x10^6


# report s and pi for each of these models
i=12
mean(admix_joint[,i])/100
mean(admix_mode[,i])/100
mean(admix_posterior[,i])/100
mean(admix_uniform[,i])/100
mean(admix_bot_mode[,i])/100
mean(const106[,i])/100
mean(const22106[,i])/100
mean(dadi1[,i])/100
mean(dadi2[,i])/100


setwd("/Users/nanditagarud/Documents/Jensen_response/analysis/data_10kb_scan")

chr3R_509 = read.table('H12_scan_10kb_chr3R_rho5e-9removed.txt')
chr3L_509 = read.table('H12_scan_10kb_chr3L_rho5e-9removed.txt')
chr2L_509 = read.table('H12_scan_10kb_chr2L_rho5e-9removed.txt')
chr2R_509 = read.table('H12_scan_10kb_chr2R_rho5e-9removed.txt')
h=c(chr3R_509[,9],chr3L_509[,9],chr2R_509[,9],chr2L_509[,9])

#peaks=read.table('peaks_082712_400_50_2_2_0_pan_509_top50.txt')

setwd("/Users/nanditagarud/Documents/Jensen_response/analysis")


# plot:

pdf("Figure6B.pdf",width=3.27,height=3, title = "Figure S4",paper="special")

tcl=0.2

pdjy=2.5
pdjx=-3.0
  
cex1=0.5
cex2=0.6
cex3=1

par(cex.axis=cex1) #The magnification to be used for axis annotation relative to the current setting of cex.
par(cex.lab=cex1) #The magnification to be used for x and y labels relative to the current setting of cex.
par(omi=c(0,0,0,0)) #A vector of the form c(bottom, left, top, right) giving the size of the outer margins in inches.
par(mai=c(0.5,0.5,0.5,0.5)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
shift=0.02

boxplot(h, admix_mode_corr[1:150000,9], admix_bot_mode_corr[1:150000,9],admix_posterior_corr[1:150000,9],const106[1:150000,9],const22106[1:150000,9], dadi1[1:150000,9],dadi2[1:150000,9], ylim=c(0,1), col=c('black','orange','purple','brown','forestgreen','blue','red','gold'), pch=20, cex=0.4, border=c('black','orange','purple','brown','forestgreen','blue','red','gold'), xaxs="i",yaxs = "i", xlab="",ylab="", axes=F )
#points(rep(1,50), peaks[,28], col='red', pch=20, cex=0.7)

axis(1,at=c(0,8),labels=F,tcl=0)

axis(2,at=c(0,1),labels=F,tcl=0)
axis(2,at=seq(0,1,by=0.05),labels=seq(0,1,by=0.05),padj=pdjy,tcl=-tcl)
axis(3,at=c(0,8),labels=F,tcl=0)
axis(4,at=c(0,1),labels=F,tcl=0)

mtext(expression("H12"),side=2,cex=cex2,padj=-2.25)

dev.off() 

###############
# Zoom in     #
###############

boxplot(h, admix_mode_corr[1:150000,9], admix_bot_mode_corr[1:150000,9], admix_posterior_corr[1:150000,9], ylim=c(0,0.1), col=c('grey','orange','purple','brown'), ylab='H12', names=c('data','admix_mode','admix_bot_mode','admix_95CI'), main='10KB windows, Zoom in')


######################################
# plot the H12 scan in 10kb windows  #
######################################

setwd('/Users/nanditagarud/Documents/Jensen_response/analysis/data_10kb_scan')

# recombination >=5*10^-9
chr3R_509 = read.table('H12_scan_10kb_chr3R_rho5e-9removed.txt')
chr3L_509 = read.table('H12_scan_10kb_chr3L_rho5e-9removed.txt')
chr2L_509 = read.table('H12_scan_10kb_chr2L_rho5e-9removed.txt')
chr2R_509 = read.table('H12_scan_10kb_chr2R_rho5e-9removed.txt')

# recombination <5*10^-9
chr3R_509_low = read.table('H12_scan_10kb_chr3R_rho5e-9_only.txt')
chr3L_509_low = read.table('H12_scan_10kb_chr3L_rho5e-9_only.txt')
chr2L_509_low = read.table('H12_scan_10kb_chr2L_rho5e-9_only.txt')
chr2R_509_low = read.table('H12_scan_10kb_chr2R_rho5e-9_only.txt')

plot(chr3R_509_low[,1], chr3R_509_low[,9], col='grey', pch=20, xlab='position', ylab='H12 measured in 10kb windows', main='Chr3R')
points(chr3R_509[,1], chr3R_509[,9], pch=20 )
top_candidates=c(560, 3344, 14136, 16974)
points(chr3R_509[top_candidates,1], chr3R_509[top_candidates,9], pch=20, col='red' )

plot(chr3L_509_low[,1], chr3L_509_low[,9], col='grey', pch=20, xlab='position', ylab='H12 measured in 10kb windows', main='Chr3L')
top_candidates=c(4060, 5201, 11537)
points(chr3L_509[,1], chr3L_509[,9], pch=20 )
points(chr3L_509[top_candidates,1], chr3L_509[top_candidates,9], pch=20, col='red' )

plot(chr2R_509_low[,1], chr2R_509_low[,9], col='grey', pch=20, xlab='position', ylab='H12 measured in 10kb windows', main='Chr2R')
top_candidates=c(3188, 6236, 9991,14663)
points(chr2R_509[,1], chr2R_509[,9], pch=20 )
points(chr2R_509[top_candidates,1], chr2R_509[top_candidates,9], pch=20, col='red' )

plot(chr2L_509_low[,1], chr2L_509_low[,9], col='grey', pch=20, xlab='position', ylab='H12 measured in 10kb windows', main='Chr2L')
top_candidates=c(5227,6619,7812,13251,14834)
points(chr2L_509[,1], chr2L_509[,9], pch=20 )
points(chr2L_509[top_candidates,1], chr2L_509[top_candidates,9], pch=20, col='red' )


#######################################################
# plot H12 values vs num SNPs in 10KB window for data #
#######################################################

plot(chr3R_509[,12], chr3R_509[,9], xlab='Number of SNPs in 10Kb window', ylab='H12', main='Chr 3R')
plot(chr3L_509[,12], chr3L_509[,9], xlab='Number of SNPs in 10Kb window', ylab='H12', main='Chr 3L')
plot(chr2R_509[,12], chr2R_509[,9], xlab='Number of SNPs in 10Kb window', ylab='H12', main='Chr 2R')
plot(chr2L_509[,12], chr2L_509[,9], xlab='Number of SNPs in 10Kb window', ylab='H12', main='Chr 2L')

##############################################################
# plot H12 values vs num SNPs in 10KB window for simulations #
##############################################################

setwd("/Users/nanditagarud/Documents/Jensen_response/analysis")
admix_posterior_neutral=read.table('Admixture_posteriorDistn_95CI_rho_5e-9_theta_0_selection_False_bps.txt')
admix_posterior_HS=read.table('Admixture_posteriorDistn_95CI_rho_5e-9_theta_0.01_selection_True_bps.txt')
admix_posterior_SS=read.table('Admixture_posteriorDistn_95CI_rho_5e-9_theta_10_selection_True_bps.txt')

plot(admix_posterior_neutral[1:150000,12], admix_posterior_neutral[1:150000,9], xlab='Number of SNPs in 10Kb window', ylab='H12', main='neutrality', ylim=c(0,1))

plot(admix_posterior_HS[1:150000,12], admix_posterior_HS[1:150000,9], xlab='Number of SNPs in 10Kb window', ylab='H12', main='HS')

plot(admix_posterior_SS[1:150000,12], admix_posterior_SS[1:150000,9], xlab='Number of SNPs in 10Kb window', ylab='H12', main='SS')


##########################################################
# Compute mean pi and S for different demographic models #
##########################################################
setwd('/Users/nanditagarud/Documents/Jensen_response/analysis')
admix_posterior_Pi_S=read.table('Admixture_posteriorDistn_95CI_rho_5e-9_theta_0_selection_False_Pi_S_bps.txt')
admix_posterior_noSingletons_Pi_S=read.table('Admixture_posteriorDistn_95CI_rho_5e-9_theta_0_selection_False_Pi_S_noSingletons.txt')

admix_posterior_corrected_Pi_S=read.table('Admixture_posteriorDistn_95CI_corrected_rho_5e-9_theta_0_selection_False_Pi_S_TajD_bps.txt')
admix_posterior_corrected_noSingletons_Pi_S=read.table('Admixture_posteriorDistn_95CI_corrected_rho_5e-9_theta_0_selection_False_Pi_S_TajD_noSingletons.txt')

admix_posterior_outside_Pi_S=read.table('Admixture_posteriorDistn_outside_95CI_rho_5e-9_theta_0_selection_False_Pi_S_bps.txt')
admix_posterior_outside_noSingletons_Pi_S=read.table('Admixture_posteriorDistn_outside_95CI_rho_5e-9_theta_0_selection_False_Pi_S_noSingletons.txt')

admix_uniform_Pi_S=read.table('Admixture_uniformPrior_rho_5e-9_theta_0_selection_False_Pi_S_bps.txt') #something wrong with model
admix_uniform_noSingletons_Pi_S=read.table('Admixture_uniformPrior_rho_5e-9_theta_0_selection_False_Pi_S_noSingletons.txt')

admix_bot_mode_Pi_S=read.table('admixture_bottleneck_mode_neutrality.txt_Pi_S.txt')
admix_bot_mode_noSingletons_Pi_S=read.table('admixture_bottleneck_mode_neutrality.txt_Pi_S_noSingletons.txt') #fix

admix_bot_mode_corrected_Pi_S=read.table('admixture_bottleneck_mode_corrected_neutrality.txt_Pi_S_TajD.txt')
admix_bot_mode_corrected_no_Singletons_Pi_S=read.table('admixture_bottleneck_mode_corrected_neutrality.txt_Pi_S_TajD_noSingletons.txt')


admix_mode_Pi_S=read.table('Admixture_mode_rho_5e-9_theta_0_selection_False_Pi_S.txt')
admix_mode_noSingletons_Pi_S=read.table('Admixture_mode_rho_5e-9_theta_0_selection_False_Pi_S_noSingletons.txt')
admix_mode_DuchenRhoTheta_Pi_S=read.table('Admixture_mode_rho_5e-9_theta_0_selection_False_Pi_S_Singletons_DuchenRhoTheta.txt')
admix_mode_DuchenRhoTheta_noSingletons_Pi_S=read.table('Admixture_mode_rho_5e-9_theta_0_selection_False_Pi_S_noSingletons_DuchenRhoTheta.txt')
admix_mode_corrected_Pi_S=read.table('Admixture_mode_corrected_rho_5e-9_theta_0_selection_False_Pi_S_TajD.txt')
admix_mode_corrected_noSingletons_Pi_S=read.table('Admixture_mode_corrected_rho_5e-9_theta_0_selection_False_Pi_S_TajD_noSingletons.txt')

const106_Pi_S=read.table('constNe106_neutrality_Pi_S.txt')
const22106_Pi_S=read.table('constNe2.2e6_neutrality_Pi_S.txt')
dadi1_Pi_S=read.table('dadi1_neutrality_Pi_S.txt')
dadi2_Pi_S=read.table('dadi2_neutrality_Pi_S.txt')
dadi1_recoded_Pi_S=read.table('dadi1_recoded_neutrality_Pi_S_TajD.txt')
dadi2_recoded_Pi_S=read.table('dadi2_neutrality_recoded_Pi_S_TajD.txt')

mean(admix_posterior_Pi_S[,1])



##########################################################################################
# plot distribution of Pi, S, TajD for the 242 fragments simulated from the X chromosome #
##########################################################################################

setwd('/Users/nanditagarud/Documents/Jensen_response/analysis')

X_chr_noSingletons=read.table('Admixture_95CI_rho_5e-9_theta_0_selection_False_DuchenRhoTheta_exactLength_Xchr_Pi_S_TajD_noSingletons_DuchenRhoTheta.txt')
X_chr_Singletons=read.table('Admixture_95CI_rho_5e-9_theta_0_selection_False_DuchenRhoTheta_exactLength_Xchr_Pi_S_TajD_Singletons_DuchenRhoTheta.txt')
X_chr_fullPosterior_noSingletons=read.table('Admixture_fullPosterior_rho_5e-9_theta_0_selection_False_DuchenRhoTheta_exactLength_Xchr_Pi_S_TajD_noSingletons_DuchenRhoTheta.txt')
X_chr_fullPosterior_Singletons=read.table('Admixture_fullPosterior_rho_5e-9_theta_0_selection_False_DuchenRhoTheta_exactLength_Xchr_Pi_S_TajD_Singletons_DuchenRhoTheta.txt')


boxplot(X_chr_noSingletons, names=c('S_n','Pi_n', 'average pairwise distance','TajD'), main='X chr, 242 fragments, 95CI, no singletons')
boxplot(X_chr_fullPosterior_noSingletons, names=c('S_n','Pi_n', 'average pairwise distance','TajD'), main='X chr, 242 fragments, full posterior distribution, no singletons')
boxplot(X_chr_Singletons, names=c('S_n','Pi_n', 'average pairwise distance','TajD'), main='X chr, 242 fragments, 95CI, singletons')
boxplot(X_chr_fullPosterior_Singletons, names=c('S_n','Pi_n', 'average pairwise distance','TajD'), main='X chr, 242 fragments, full posterior distribution, singletons')


i= X_chr_Singletons[,1]<=791
boxplot(X_chr_Singletons[i,], names=c('S_n','Pi_n', 'TajD'))








#####################################################
# BF heat map                                       #
# Red and Grey scale, truncated x axis, 3 panels    #
# Panel 1: Original plot, SNP based                 #
# Panel 2: More comprehensive dem model, SNP based  #
# Panel 3: More comprehensive dem model, bp based   #
#####################################################

setwd("/Users/nanditagarud/Documents/Jensen_response/analysis")

peaks=read.table('peaks_082712_400_50_2_2_0_pan_509_top50_rhoForCalling.txt')

pdf("fig8_BFheatmap_panels.pdf",width=6.83, height=2.43, title = "Figure 2",paper="special")
#width=6.83,height=4.86
#width=4.86,height=2.43

s=as.vector(c(1,6))
split.screen(s)

change=0.06

mb1=0.1 
ml1=0.4
mt1=0 
mr1=0 

mb2=0.1 
ml2=0.35 
mt2=0 
mr2=0.05

mb3=0.1 
ml3=0.3	 
mt3=0 
mr3=0.1 

mb4=0.1 
ml4=0.25	 
mt4=0 
mr4=0.15 

mb5=0.1 
ml5=0.2	 
mt5=0 
mr5=0.2 

mb6=0.1 
ml6=0.15	 
mt6=0 
mr6=0.25 

tcl=0.2

pdjy=2.0
pdjx=-2.5
  
cex1=0.7
cex2=0.7
cex3=1

par(cex.axis=cex1) #The magnification to be used for axis annotation relative to the current setting of cex.
par(cex.lab=cex1) #The magnification to be used for x and y labels relative to the current setting of cex.
par(omi=c(0.35,0.1,0.35,0.1)) #A vector of the form c(bottom, left, top, right) giving the size of the outer margins in inches.

shift=0.02


screen(1)

par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
par(xaxs = "i")
par(yaxs = "i")

#d=read.table("Admixture_mode_fixedPopSize_NAm61659_EUr2000000_rho_5e-9_BFs_numSim_12000.txt")
d=read.table("Admixture_mode_fixedPopSize_NAm61659_EUr500000_conditionEndFreq_rho_5e-9_BFs_numSim_14000.txt")
d2=(d[,4])/(d[,3])
d3=matrix(d2,c(40,40))
d4=d3[,1:12]


colorMatrix=as.matrix(d4)


for (i in 1:length(d4[,1])){
	for (j in 1:length(d4[1,])){
		

		if (d3[i,j]=='NaN'){colorMatrix[i,j]="white"} else 
		{if (d3[i,j]=='Inf'){colorMatrix[i,j]="black"} else 
		{if (d3[i,j]>=0 && d3[i,j]<1){colorMatrix[i,j]="red"} else
		{if (d3[i,j]>=1 && d3[i,j]<10){colorMatrix[i,j]="gray65"} else
		{if (d3[i,j]>=10 && d3[i,j]<30){colorMatrix[i,j]="gray55"} else
		{if (d3[i,j]>=30 && d3[i,j]<100){colorMatrix[i,j]="gray35"} else
		{if (d3[i,j]>=100 && d3[i,j]<1000000){colorMatrix[i,j]="black"} else
			{colorMatrix[i,j]=3}}}}}}}
		
	}	
}


plot(c(0, 12), c(0, 40), type= "n", xaxs="i",yaxs = "i",axes=F, main='', xlim=c(0,13), ylim=c(0,41))
h=0
w=0.35
for (i in 1:length(d4[,1])){
	height=c(h-0.5, h+0.5)
	for (j in 1:length(d4[1,])){
		rect(j-0.5, height[1], j+0.5, height[2],density=-1,col=colorMatrix[i,j],lwd=0., border=colorMatrix[i,j])

	}
	h=h+1	
} 

axis(1,at=c(0,12),labels=F,tcl=0)
axis(1, at = c(0,6,12), labels=c("0","0.15","0.30"), padj=pdjx,tcl=-tcl)
axis(2,at=c(0,40),labels=F,tcl=0)
axis(2, at = c(0,10,20,30,40), labels=c("0","0.25","0.50", "0.75", "1" ),padj=pdjy,tcl=-tcl)


# overlay the observed peak's values
points(40*peaks[,28],40*peaks[,20],col='gold',pch=20, cex=0.6)
mtext(expression(H2~"/"~H1),side=2,cex=cex2,padj=-3)
mtext(expression(A),at=-0.5,padj=-0.5,cex=cex3)


screen(2)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
par(xaxs = "i")
par(yaxs = "i")

d=read.table("Admixture_mode_fixedPopSize_NAm61659_EUr500000_rho_5e-9_BFs_numSim_140000.txt")
d2=(d[,4])/(d[,3])
d3=matrix(d2,c(40,40))
d4=d3[,1:12]


colorMatrix=as.matrix(d4)


for (i in 1:length(d4[,1])){
	for (j in 1:length(d4[1,])){
		

		if (d3[i,j]=='NaN'){colorMatrix[i,j]="white"} else 
		{if (d3[i,j]=='Inf'){colorMatrix[i,j]="black"} else 
		{if (d3[i,j]>=0 && d3[i,j]<1){colorMatrix[i,j]="red"} else
		{if (d3[i,j]>=1 && d3[i,j]<10){colorMatrix[i,j]="gray65"} else
		{if (d3[i,j]>=10 && d3[i,j]<30){colorMatrix[i,j]="gray55"} else
		{if (d3[i,j]>=30 && d3[i,j]<100){colorMatrix[i,j]="gray35"} else
		{if (d3[i,j]>=100 && d3[i,j]<1000000){colorMatrix[i,j]="black"} else
			{colorMatrix[i,j]=3}}}}}}}
		
	}	
}


plot(c(0, 12), c(0, 40), type= "n", xaxs="i",yaxs = "i",axes=F, main='', xlim=c(0,13))
h=0
w=0.35
for (i in 1:length(d4[,1])){
	height=c(h-0.5, h+0.5)
	for (j in 1:length(d4[1,])){
		rect(j-0.5, height[1], j+0.5, height[2],density=-1,col=colorMatrix[i,j],lwd=0., border=colorMatrix[i,j])

	}
	h=h+1	
}

axis(1,at=c(0,12),labels=F,tcl=0)
axis(1, at = c(0,6,12), labels=c("0","0.15","0.30"), padj=pdjx,tcl=-tcl)
axis(2,at=c(0,40),labels=F,tcl=0)
axis(2, at = c(0,10,20,30,40), labels=c("0","0.25","0.50", "0.75", "1" ),padj=pdjy,tcl=-tcl)

# overlay the observed peak's values
points(40*peaks[,28],40*peaks[,20],col='gold',pch=20, cex=0.6)
mtext(expression(B),at=-0.5,padj=-0.5,cex=cex3)



screen(3)

par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
par(xaxs = "i")
par(yaxs = "i")

d=read.table("Admixture_mode_fixedPopSize_NAm1110000_EUr2000000_rho_5e-9_BFs_numSim_20000.txt")
d2=(d[,4])/(d[,3])
d3=matrix(d2,c(40,40))
d4=d3[,1:12]


colorMatrix=as.matrix(d4)


for (i in 1:length(d4[,1])){
	for (j in 1:length(d4[1,])){
		

		if (d3[i,j]=='NaN'){colorMatrix[i,j]="white"} else 
		{if (d3[i,j]=='Inf'){colorMatrix[i,j]="black"} else 
		{if (d3[i,j]>=0 && d3[i,j]<1){colorMatrix[i,j]="red"} else
		{if (d3[i,j]>=1 && d3[i,j]<10){colorMatrix[i,j]="gray65"} else
		{if (d3[i,j]>=10 && d3[i,j]<30){colorMatrix[i,j]="gray55"} else
		{if (d3[i,j]>=30 && d3[i,j]<100){colorMatrix[i,j]="gray35"} else
		{if (d3[i,j]>=100 && d3[i,j]<1000000){colorMatrix[i,j]="black"} else
			{colorMatrix[i,j]=3}}}}}}}
		
	}	
}


plot(c(0, 12), c(0, 40), type= "n", xaxs="i",yaxs = "i",axes=F, main='', xlim=c(0,13), ylim=c(0,41))
h=0
w=0.35
for (i in 1:length(d4[,1])){
	height=c(h-0.5, h+0.5)
	for (j in 1:length(d4[1,])){
		rect(j-0.5, height[1], j+0.5, height[2],density=-1,col=colorMatrix[i,j],lwd=0., border=colorMatrix[i,j])

	}
	h=h+1	
} 

axis(1,at=c(0,12),labels=F,tcl=0)
axis(1, at = c(0,6,12), labels=c("0","0.15","0.30"), padj=pdjx,tcl=-tcl)
axis(2,at=c(0,40),labels=F,tcl=0)
axis(2, at = c(0,10,20,30,40), labels=c("0","0.25","0.50", "0.75", "1" ),padj=pdjy,tcl=-tcl)


# overlay the observed peak's values
points(40*peaks[,28],40*peaks[,20],col='gold',pch=20, cex=0.6)
mtext(expression(H2~"/"~H1),side=2,cex=cex2,padj=-3)
mtext(expression(A),at=-0.5,padj=-0.5,cex=cex3)


screen(4)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
par(xaxs = "i")
par(yaxs = "i")

d=read.table("Admixture_mode_fixedPopSize_NAm1110000_EUr500000_rho_5e-9_BFs_numSim_20000.txt")
d2=(d[,4])/(d[,3])
d3=matrix(d2,c(40,40))
d4=d3[,1:12]


colorMatrix=as.matrix(d4)


for (i in 1:length(d4[,1])){
	for (j in 1:length(d4[1,])){
		

		if (d3[i,j]=='NaN'){colorMatrix[i,j]="white"} else 
		{if (d3[i,j]=='Inf'){colorMatrix[i,j]="black"} else 
		{if (d3[i,j]>=0 && d3[i,j]<1){colorMatrix[i,j]="red"} else
		{if (d3[i,j]>=1 && d3[i,j]<10){colorMatrix[i,j]="gray65"} else
		{if (d3[i,j]>=10 && d3[i,j]<30){colorMatrix[i,j]="gray55"} else
		{if (d3[i,j]>=30 && d3[i,j]<100){colorMatrix[i,j]="gray35"} else
		{if (d3[i,j]>=100 && d3[i,j]<1000000){colorMatrix[i,j]="black"} else
			{colorMatrix[i,j]=3}}}}}}}
		
	}	
}


plot(c(0, 12), c(0, 40), type= "n", xaxs="i",yaxs = "i",axes=F, main='', xlim=c(0,13))
h=0
w=0.35
for (i in 1:length(d4[,1])){
	height=c(h-0.5, h+0.5)
	for (j in 1:length(d4[1,])){
		rect(j-0.5, height[1], j+0.5, height[2],density=-1,col=colorMatrix[i,j],lwd=0., border=colorMatrix[i,j])

	}
	h=h+1	
}

axis(1,at=c(0,12),labels=F,tcl=0)
axis(1, at = c(0,6,12), labels=c("0","0.15","0.30"), padj=pdjx,tcl=-tcl)
axis(2,at=c(0,40),labels=F,tcl=0)
axis(2, at = c(0,10,20,30,40), labels=c("0","0.25","0.50", "0.75", "1" ),padj=pdjy,tcl=-tcl)

# overlay the observed peak's values
points(40*peaks[,28],40*peaks[,20],col='gold',pch=20, cex=0.6)
mtext(expression(B),at=-0.5,padj=-0.5,cex=cex3)





screen(5)

par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
par(xaxs = "i")
par(yaxs = "i")

d=read.table("Admixture_mode_fixedPopSize_NAm15984500_EUr2000000_rho_5e-9_BFs_numSim_12000.txt")
d2=(d[,4])/(d[,3])
d3=matrix(d2,c(40,40))
d4=d3[,1:12]


colorMatrix=as.matrix(d4)


for (i in 1:length(d4[,1])){
	for (j in 1:length(d4[1,])){
		

		if (d3[i,j]=='NaN'){colorMatrix[i,j]="white"} else 
		{if (d3[i,j]=='Inf'){colorMatrix[i,j]="black"} else 
		{if (d3[i,j]>=0 && d3[i,j]<1){colorMatrix[i,j]="red"} else
		{if (d3[i,j]>=1 && d3[i,j]<10){colorMatrix[i,j]="gray65"} else
		{if (d3[i,j]>=10 && d3[i,j]<30){colorMatrix[i,j]="gray55"} else
		{if (d3[i,j]>=30 && d3[i,j]<100){colorMatrix[i,j]="gray35"} else
		{if (d3[i,j]>=100 && d3[i,j]<1000000){colorMatrix[i,j]="black"} else
			{colorMatrix[i,j]=3}}}}}}}
		
	}	
}


plot(c(0, 12), c(0, 40), type= "n", xaxs="i",yaxs = "i",axes=F, main='', xlim=c(0,13), ylim=c(0,41))
h=0
w=0.35
for (i in 1:length(d4[,1])){
	height=c(h-0.5, h+0.5)
	for (j in 1:length(d4[1,])){
		rect(j-0.5, height[1], j+0.5, height[2],density=-1,col=colorMatrix[i,j],lwd=0., border=colorMatrix[i,j])

	}
	h=h+1	
} 

axis(1,at=c(0,12),labels=F,tcl=0)
axis(1, at = c(0,6,12), labels=c("0","0.15","0.30"), padj=pdjx,tcl=-tcl)
axis(2,at=c(0,40),labels=F,tcl=0)
axis(2, at = c(0,10,20,30,40), labels=c("0","0.25","0.50", "0.75", "1" ),padj=pdjy,tcl=-tcl)


# overlay the observed peak's values
points(40*peaks[,28],40*peaks[,20],col='gold',pch=20, cex=0.6)
mtext(expression(H2~"/"~H1),side=2,cex=cex2,padj=-3)
mtext(expression(A),at=-0.5,padj=-0.5,cex=cex3)


screen(6)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
par(xaxs = "i")
par(yaxs = "i")

d=read.table("Admixture_mode_fixedPopSize_NAm15984500_EUr500000_rho_5e-9_BFs_numSim_12000.txt")
d2=(d[,4])/(d[,3])
d3=matrix(d2,c(40,40))
d4=d3[,1:12]


colorMatrix=as.matrix(d4)


for (i in 1:length(d4[,1])){
	for (j in 1:length(d4[1,])){
		

		if (d3[i,j]=='NaN'){colorMatrix[i,j]="white"} else 
		{if (d3[i,j]=='Inf'){colorMatrix[i,j]="black"} else 
		{if (d3[i,j]>=0 && d3[i,j]<1){colorMatrix[i,j]="red"} else
		{if (d3[i,j]>=1 && d3[i,j]<10){colorMatrix[i,j]="gray65"} else
		{if (d3[i,j]>=10 && d3[i,j]<30){colorMatrix[i,j]="gray55"} else
		{if (d3[i,j]>=30 && d3[i,j]<100){colorMatrix[i,j]="gray35"} else
		{if (d3[i,j]>=100 && d3[i,j]<1000000){colorMatrix[i,j]="black"} else
			{colorMatrix[i,j]=3}}}}}}}
		
	}	
}


plot(c(0, 12), c(0, 40), type= "n", xaxs="i",yaxs = "i",axes=F, main='', xlim=c(0,13))
h=0
w=0.35
for (i in 1:length(d4[,1])){
	height=c(h-0.5, h+0.5)
	for (j in 1:length(d4[1,])){
		rect(j-0.5, height[1], j+0.5, height[2],density=-1,col=colorMatrix[i,j],lwd=0., border=colorMatrix[i,j])

	}
	h=h+1	
}

axis(1,at=c(0,12),labels=F,tcl=0)
axis(1, at = c(0,6,12), labels=c("0","0.15","0.30"), padj=pdjx,tcl=-tcl)
axis(2,at=c(0,40),labels=F,tcl=0)
axis(2, at = c(0,10,20,30,40), labels=c("0","0.25","0.50", "0.75", "1" ),padj=pdjy,tcl=-tcl)

# overlay the observed peak's values
points(40*peaks[,28],40*peaks[,20],col='gold',pch=20, cex=0.6)
mtext(expression(B),at=-0.5,padj=-0.5,cex=cex3)






close.screen(all=TRUE)

dev.off() 




#########
hard=read.table('params_theta_0.01.txt')
soft=read.table('params_theta_10.txt')

hard_filtered=c()
for (i in 1:length(hard[,1])){
	if (hard[i,1]>0 and hard[i,1]<1){
		hard_filtered=c(hard_filtered,hard[i,1])
		}	
	}
	
	
#################################################################
# Figure 2: recreation of H12 scan from PLOS Genetics paper
##################################################################

setwd('/Users/nanditagarud/Documents/Ace')

chr3R = read.table('clusters_082712_400_50_2_2_0_chr3R.txt')
chr3L = read.table('clusters_082712_400_50_2_2_0_chr3L.txt')
chr2L = read.table('clusters_082712_400_50_2_2_0_chr2L.txt')
chr2R = read.table('clusters_082712_400_50_2_2_0_chr2R.txt')

chr3R_509 = read.table('clusters_082712_400_50_2_2_0_chr3R_withRho_threshold.txt')
chr3L_509 = read.table('clusters_082712_400_50_2_2_0_chr3L_withRho_threshold.txt')
chr2L_509 = read.table('clusters_082712_400_50_2_2_0_chr2L_withRho_threshold.txt')
chr2R_509 = read.table('clusters_082712_400_50_2_2_0_chr2R_withRho_threshold.txt')

chr3R_peaks = read.table('peaks_082712_400_50_2_2_0_pan_509_top50_Chr3R.txt')
chr3L_peaks = read.table('peaks_082712_400_50_2_2_0_pan_509_top50_Chr3L.txt')
chr2L_peaks = read.table('peaks_082712_400_50_2_2_0_pan_509_top50_Chr2L.txt')
chr2R_peaks = read.table('peaks_082712_400_50_2_2_0_pan_509_top50_Chr2R.txt')

chr3R_overlap = read.table('iHS_scan_DGRPv1/overlap_H12_3R.txt')
chr3L_overlap = read.table('iHS_scan_DGRPv1/overlap_H12_3L.txt')
chr2L_overlap = read.table('iHS_scan_DGRPv1/overlap_H12_2L.txt')
chr2R_overlap = read.table('iHS_scan_DGRPv1/overlap_H12_2R.txt')


pdf("Fig_2_scan.pdf",width=6.83,height=4, title = "Figure 6a",paper="special")

s=as.vector(c(2,2))
split.screen(s)


change=0.05

mb1=0.35 
ml1=0.4 +change
mt1=0 
mr1=0.0 +change

mb2=0.35 
ml2=0.3 +change
mt2=0 
mr2=0.1 +change


tcl=0.2

pdjy=2
pdjx=-2.5
  
cex1=0.7
cex2=0.7
cex3=1

par(cex.axis=cex1) #The magnification to be used for axis annotation relative to the current setting of cex.
par(cex.lab=cex1) #The magnification to be used for x and y labels relative to the current setting of cex.
par(omi=c(0.35,0,0.35,0)) #A vector of the form c(bottom, left, top, right) giving the size of the outer margins in inches.


shift=200000

#Chr2L
screen(1)
par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size 

plot(chr2L[,1], chr2L[,26], ylim=c(0,0.27), xlim=c(0,chr2L[length(chr2L[,1]),1] +shift), col="grey",pch=20, xaxs="i",yaxs = "i",axes=F, main='',cex=0.4)
points(chr2L_509[,3],chr2L_509[,28], col='black', pch=20,cex=0.4)
points(chr2L_peaks[,3],chr2L_peaks[,28], col='red', pch=20,cex=1)
points(chr2L_overlap[,7],chr2L_overlap[,4], col='red', pch=20,cex=1)
abline(h= 0.0170749108205, col='orange', lty=1, lwd=2)


axis(1,at=c(0,chr2L[length(chr2L[,1]),1]+shift),labels=F,tcl=0)
axis(1,at=seq(0,chr2L[length(chr2L[,1]),1],5*10^6),labels=seq(0,chr2L[length(chr2L[,1]),1],5*10^6),padj=pdjx,tcl=-tcl)
axis(2,at=c(0,0.27),labels=F,tcl=0)
axis(2,at=seq(0,0.25,0.05),labels=seq(0,0.25,0.05),padj=pdjy,tcl=-tcl)
axis(3,at=c(0,chr2L[length(chr2L[,1]),1]+shift),labels=F,tcl=0)
axis(4,at=c(0,0.27),labels=F,tcl=0)

mtext("Chr2L", at=c(1.5*10^6), cex=cex2, padj=2, font=2)
mtext(expression("H12"),side=2,cex=cex2,padj=-3)
#mtext(expression(A),at=-10^3,padj=-1,cex=cex3)


#Chr2R
screen(2)
par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size 

plot(chr2R[,1], chr2R[,26], ylim=c(0,0.27), xlim=c(0,chr2R[length(chr2R[,1]),1] +shift), col="grey",pch=20, xaxs="i",yaxs = "i",axes=F, main='',cex=0.4)
points(chr2R_509[,3],chr2R_509[,28], col='black', pch=20,cex=0.4)
points(chr2R_peaks[,3],chr2R_peaks[,28], col='red', pch=20,cex=1)
points(chr2R_overlap[,7],chr2R_overlap[,4], col='red', pch=20,cex=1)
abline(h= 0.0170749108205, col='orange', lty=1, lwd=2)

axis(1,at=c(0,chr2R[length(chr2R[,1]),1]+shift),labels=F,tcl=0)
axis(1,at=seq(0,chr2R[length(chr2R[,1]),1],5*10^6),labels=seq(0,chr2R[length(chr2R[,1]),1],5*10^6),padj=pdjx,tcl=-tcl)
axis(2,at=c(0,0.27),labels=F,tcl=0)
axis(2,at=seq(0,0.25,0.05),labels=seq(0,0.25,0.05),padj=pdjy,tcl=-tcl)
axis(3,at=c(0,chr2R[length(chr2R[,1]),1]+shift),labels=F,tcl=0)
axis(4,at=c(0,0.27),labels=F,tcl=0)

mtext("Chr2R", at=c(1.5*10^6), cex=cex2, padj=2, font=2)
mtext("Cyp6g1", at=c(8074551), cex=cex2, padj=1.8, font=3)


#Chr3L
screen(3)
par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size 

plot(chr3L[,1], chr3L[,26], ylim=c(0,0.27), xlim=c(0,chr3L[length(chr3L[,1]),1] +shift), col="grey",pch=20, xaxs="i",yaxs = "i",axes=F, main='',cex=0.4)
points(chr3L_509[,3],chr3L_509[,28], col='black', pch=20,cex=0.4)
points(chr3L_peaks[,3],chr3L_peaks[,28], col='red', pch=20,cex=1)
points(chr3L_overlap[,7],chr3L_overlap[,4], col='red', pch=20,cex=1)
abline(h= 0.0170749108205, col='orange', lty=1, lwd=2)

axis(1,at=c(0,chr3L[length(chr3L[,1]),1]+shift),labels=F,tcl=0)
axis(1,at=seq(0,chr3L[length(chr3L[,1]),1],5*10^6),labels=seq(0,chr3L[length(chr3L[,1]),1],5*10^6),padj=pdjx,tcl=-tcl)
axis(2,at=c(0,0.27),labels=F,tcl=0)
axis(2,at=seq(0,0.25,0.05),labels=seq(0,0.25,0.05),padj=pdjy,tcl=-tcl)
axis(3,at=c(0,chr3L[length(chr3L[,1]),1]+shift),labels=F,tcl=0)
axis(4,at=c(0,0.27),labels=F,tcl=0)

mtext("Chr3L", at=c(1.5*10^6), cex=cex2, padj=2, font=2)
mtext(expression("H12"),side=2,cex=cex2,padj=-3)
mtext(expression("position"),side=1,cex=cex2,padj=2)



#Chr3R
screen(4)
par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size 

plot(chr3R[,1], chr3R[,26], ylim=c(0,0.27), xlim=c(0,chr3R[length(chr3R[,1]),1] +shift), col="grey",pch=20, xaxs="i",yaxs = "i",axes=F, main='',cex=0.4)
points(chr3R_509[,3],chr3R_509[,28], col='black', pch=20,cex=0.4)
points(chr3R_peaks[,3],chr3R_peaks[,28], col='red', pch=20,cex=1)
points(chr3R_overlap[,7],chr3R_overlap[,4], col='red', pch=20,cex=1)
abline(h= 0.0170749108205, col='orange', lty=1, lwd=2)

axis(1,at=c(0,chr3R[length(chr3R[,1]),1]+shift),labels=F,tcl=0)
axis(1,at=seq(0,chr3R[length(chr3R[,1]),1],5*10^6),labels=seq(0,chr3R[length(chr3R[,1]),1],5*10^6),padj=pdjx,tcl=-tcl)
axis(2,at=c(0,0.27),labels=F,tcl=0)
axis(2,at=seq(0,0.25,0.05),labels=seq(0,0.25,0.05),padj=pdjy,tcl=-tcl)
axis(3,at=c(0,chr3R[length(chr3R[,1]),1]+shift),labels=F,tcl=0)
axis(4,at=c(0,0.27),labels=F,tcl=0)

mtext("Chr3R", at=c(2*10^6), cex=cex2, padj=2, font=2)
mtext(expression("position"),side=1,cex=cex2,padj=2)
mtext("Ace", at=c(9069408), cex=cex2, padj=2, font=3)
mtext("CHKov1", at=c(21152026), cex=cex2, padj=2, font=3)



close.screen(all=TRUE)

dev.off() 




#######################################################################
# Figure 6                                                            #
# Haplotype frequency spectra for top 10 peaks   #
#######################################################################

setwd('/Users/nanditagarud/Documents/Ace')

peaks=read.table('peaks_082712_400_50_2_2_0_pan_509_top50.txt')

pdf("Figure6.pdf",width=4.25,height=4.86, title = "Figure8",paper="special")
s=as.vector(c(1,10))
split.screen(s)
par(omi=c(0.1,0.30,0.15,0.05)) #A vector of the form c(bottom, left, top, right) giving the size of the outer margins in inches.

cex2=0.7

x=1
screen(x)
# need to adjust the margins
par(mai=c(.1,0,.05,0.03))

a1 <- as.character(peaks[x,5]) 
a2 <- strsplit(a1, ",") 
a3 <- unlist(a2) 
hapDist <- as.numeric(a3) 

sampleSize=145
NumOnes=sampleSize-sum(hapDist)
wl=0; wr=50;

cols=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#bc80bd','#ccebc5','#c51b7d','#7fbc41','#f1b6da','#8c510a','#c7eae5','#35978f','#dfc27d','#d6604d')

plot(1:2,1:2,col="white",ylim=c(-1.3,sampleSize),xlim=c(-1,55),ylab="",xlab="",yaxt="n",xaxt="n",frame.plot=FALSE)

for (i in 1: NumOnes){
	height=c(i-0.35,i+0.35)
	rect(wl,height[1],wr,height[2],density=-1,col="grey",lwd=0, border='grey')}

print(x)
colNo=1
for (j in sort(hapDist)){
	height=c(height[2]+0.3,height[2]+j)
	rect(wl,height[1],wr,height[2],density=-1,col=cols[length(hapDist)-colNo+1],lwd=0, border='white')
	colNo=colNo+1
	}

# Print H12 and H2/H1 values
H12=paste(round(peaks[x,28],2))
H2H1=paste(round(peaks[x,20],2))
mtext(H12,side=1,cex=cex2,at=30, padj=-3)
mtext(H2H1,side=1,cex=cex2,at=30, padj=-1.5)
#mtext(expression(A),at=-20,padj=0,cex=0.9)
mtext("H12:",side=1,at=-30,padj=-3, cex=cex2)
mtext("H2/H1:",side=1,at=-22,padj=-1.5, cex=cex2)
mtext("Cyp6g1",at=20,padj=2, cex=cex2,font=3)

x=2
screen(x)
# need to adjust the margins
par(mai=c(.1,0,.05,0.03))

a1 <- as.character(peaks[x,5]) 
a2 <- strsplit(a1, ",") 
a3 <- unlist(a2) 
hapDist <- as.numeric(a3) 

sampleSize=145
NumOnes=sampleSize-sum(hapDist)
wl=0; wr=50;

cols=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#bc80bd','#ccebc5','#c51b7d','#7fbc41','#f1b6da','#8c510a','#c7eae5','#35978f','#dfc27d','#d6604d')

plot(1:2,1:2,col="white",ylim=c(-1.3,sampleSize),xlim=c(-1,55),ylab="",xlab="",yaxt="n",xaxt="n",frame.plot=FALSE)

for (i in 1: NumOnes){
	height=c(i-0.35,i+0.35)
	rect(wl,height[1],wr,height[2],density=-1,col="grey",lwd=0, border='grey')}

print(x)
colNo=1
for (j in sort(hapDist)){
	height=c(height[2]+0.3,height[2]+j)
	rect(wl,height[1],wr,height[2],density=-1,col=cols[length(hapDist)-colNo+1],lwd=0, border='white')
	colNo=colNo+1
	}

# Print H12 and H2/H1 values
H12=paste(round(peaks[x,28],2))
H2H1=paste(round(peaks[x,20],2))
mtext(H12,side=1,cex=cex2,padj=-3)
mtext(H2H1,side=1,cex=cex2,padj=-1.5)
mtext("CHKov1",at=27,padj=2, cex=cex2,font=3)


x=3
screen(x)
# need to adjust the margins
par(mai=c(.1,0,.05,0.03))

a1 <- as.character(peaks[x,5]) 
a2 <- strsplit(a1, ",") 
a3 <- unlist(a2) 
hapDist <- as.numeric(a3) 

sampleSize=145
NumOnes=sampleSize-sum(hapDist)
wl=0; wr=50;

cols=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#bc80bd','#ccebc5','#c51b7d','#7fbc41','#f1b6da','#8c510a','#c7eae5','#35978f','#dfc27d','#d6604d')

plot(1:2,1:2,col="white",ylim=c(-1.3,sampleSize),xlim=c(-1,55),ylab="",xlab="",yaxt="n",xaxt="n",frame.plot=FALSE)

for (i in 1: NumOnes){
	height=c(i-0.35,i+0.35)
	rect(wl,height[1],wr,height[2],density=-1,col="grey",lwd=0, border='grey')}

print(x)
colNo=1
for (j in sort(hapDist)){
	height=c(height[2]+0.3,height[2]+j)
	rect(wl,height[1],wr,height[2],density=-1,col=cols[length(hapDist)-colNo+1],lwd=0, border='white')
	colNo=colNo+1
	}

# Print H12 and H2/H1 values
H12=paste(round(peaks[x,28],2))
H2H1=paste(round(peaks[x,20],2))
mtext(H12,side=1,cex=cex2,padj=-3)
mtext(H2H1,side=1,cex=cex2,padj=-1.5)
mtext("Ace",at=25,padj=2, cex=cex2,font=3)

for (x in c(4:10)){
screen(x)
# need to adjust the margins
par(mai=c(.1,0,.05,0.03))

a1 <- as.character(peaks[x,5]) 
a2 <- strsplit(a1, ",") 
a3 <- unlist(a2) 
hapDist <- as.numeric(a3) 

sampleSize=145
NumOnes=sampleSize-sum(hapDist)
wl=0; wr=50;

cols=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#bc80bd','#ccebc5','#c51b7d','#7fbc41','#f1b6da','#8c510a','#c7eae5','#35978f','#dfc27d','#d6604d')

plot(1:2,1:2,col="white",ylim=c(-1.3,sampleSize),xlim=c(-1,55),ylab="",xlab="",yaxt="n",xaxt="n",frame.plot=FALSE)

for (i in 1: NumOnes){
	height=c(i-0.35,i+0.35)
	rect(wl,height[1],wr,height[2],density=-1,col="grey",lwd=0, border='grey')}

print(x)
colNo=1
for (j in sort(hapDist)){
	height=c(height[2]+0.3,height[2]+j)
	rect(wl,height[1],wr,height[2],density=-1,col=cols[length(hapDist)-colNo+1],lwd=0, border='white')
	colNo=colNo+1
	}

# Print H12 and H2/H1 values
H12=paste(round(peaks[x,28],2))
H2H1=paste(round(peaks[x,20],2))
mtext(H12,side=1,cex=cex2,padj=-3)
mtext(H2H1,side=1,cex=cex2,padj=-1.5)
}

close.screen(all=TRUE)
dev.off() 


#####################################################################
# Figure S16
# plot the haplotype spectra for the rest of the top 50 peaks
#####################################################################

peaks=read.table('peaks_082712_400_50_2_2_0_pan_509_top50.txt')


pdf("Figure_S16_b.pdf",width=6.83,height=4.6, title = "HapDist",paper="special")
s=as.vector(c(1,20))
split.screen(s)
par(omi=c(0.1,0.30,0.15,0.05)) #A vector of the form c(bottom, left, top, right) giving the size of the outer margins in inches.

cex2=0.6
for (x in 31:50){
#for (x in 11:30){
#screen(x-10)
screen(x-30)
# need to adjust the margins
par(mai=c(.1,0,.05,0.03))

a1 <- as.character(peaks[x,5]) 
a2 <- strsplit(a1, ",") 
a3 <- unlist(a2) 
hapDist <- as.numeric(a3) 

sampleSize=145
NumOnes=sampleSize-sum(hapDist)
wl=0; wr=50;

cols=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#bc80bd','#ccebc5','#c51b7d','#7fbc41','#f1b6da','#8c510a','#c7eae5','#35978f','#dfc27d','#d6604d')

plot(1:2,1:2,col="white",ylim=c(-1.3,sampleSize),xlim=c(-1,55),ylab="",xlab="",yaxt="n",xaxt="n",frame.plot=FALSE)

for (i in 1: NumOnes){
	height=c(i-0.35,i+0.35)
	rect(wl,height[1],wr,height[2],density=-1,col="grey",lwd=0, border='grey')}

print(x)
colNo=1
for (j in sort(hapDist)){
	height=c(height[2]+0.3,height[2]+j)
	#rect(wl,height[1],wr,height[2],density=-1,col=cols[floor(height[1])%%15+1],lwd=0, border=cols[floor(height[1])%%15+1])}
	rect(wl,height[1],wr,height[2],density=-1,col=cols[length(hapDist)-colNo+1],lwd=0, border='white')
	colNo=colNo+1
	}

# Print H12 and H2/H1 values
H12=paste(round(peaks[x,28],2))
H2H1=paste(round(peaks[x,20],2))
mtext(H12,side=1,cex=cex2,padj=-2)
mtext(H2H1,side=1,cex=cex2,padj=-0.5)

}
#mtext(expression(A),at=-1300,padj=0,cex=0.9)
mtext(expression(B),at=-1300,padj=0,cex=0.9)

mtext("H12:",side=1,at=-1305,padj=-2, cex=cex2)
mtext("H2/H1:",side=1,at=-1295,padj=-0.5, cex=cex2)

close.screen(all=TRUE)
dev.off() 



################################################################
# Plot the distribution of H12 values and fit a normal to the bulk (25th-75th percentiles) #


setwd("/Users/nanditagarud/Documents/Ace/")

chr3R_509 = read.table('clusters_082712_400_50_2_2_0_chr3R_withRho_threshold.txt')
chr3L_509 = read.table('clusters_082712_400_50_2_2_0_chr3L_withRho_threshold.txt')
chr2L_509 = read.table('clusters_082712_400_50_2_2_0_chr2L_withRho_threshold.txt')
chr2R_509 = read.table('clusters_082712_400_50_2_2_0_chr2R_withRho_threshold.txt')
h=sort(c(chr3R_509[,28],chr3L_509[,28],chr2R_509[,28],chr2L_509[,28]))

peaks=read.table('peaks_082712_400_50_2_2_0_pan_509_top50.txt')


z = hist(h) # or hist(x,plot=FALSE) to avoid the plot of the histogram
z$density = z$counts/sum(z$counts)*100
plot(z,freq=FALSE, ylim=c(0,1000))
curve(dnorm(x,mean=mean(h), sd=sd(h)),add=TRUE, yaxt="n", col='red')

lower=floor(length(h)*0.25)
upper=floor(length(h)*0.75)
bulk=h[lower:upper]

curve(dnorm(x,mean=mean(bulk), sd=sd(bulk)),add=TRUE, yaxt="n", col='green')

abline(v=peaks[1,28], col='skyblue')
abline(v=peaks[50,28], col='skyblue')

##################### 
# grab points that are within 1SD around the median
median_val=median(h)
SD=0.341*median_val 
upper=median_val+SD
lower=median_val-SD

distr=h[h<upper]
standard_dev=sd(distr)

hist(h, breaks=100, main='distribution of H12 scores', ylab='num windows', xlab='H12', col='black')
abline(v=upper, col='red')
abline(v=lower, col='red')
abline(v=median_val+standard_dev, col='skyblue')
abline(v=median_val-standard_dev, col='skyblue')
abline(v=median_val+3*standard_dev, col='blue')
abline(v=median_val-3*standard_dev, col='blue')
abline(v=median_val+5*standard_dev, col='purple')

abline(v=peaks[1,28], col='gold')
abline(v=peaks[50,28], col='green')

# what fraction of points lie within 1 sd of the median?
w=h[lower<h]
sum(h<upper)/length(h)

# what fraction of points lie within 1sd of the fitted normal?
w=h[(median_val-1*standard_dev)<h]
w2=w[w<(median_val+1*standard_dev)]
length(w2)/length(h)

#qq plot of the w2 distr.
random_norm=rnorm(length(w2), mean(w2), sd(w2))
qqplot(w2, random_norm, ylab='Random sample from normal distribution',xlab='H12 values within 1SD of fitted normal distribution', main='Quantile-quantile plot assessing fit of inferred normal distribution to the bulk of the distribution')
abline(0, 1, col = 'red')

