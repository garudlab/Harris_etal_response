#####################################################
# BF heat map                                       #
# Red and Grey scale, truncated x axis, 3 panels    #
# Panel 1: EUr 7*10^5, NAm 1.1*10^6                 #
# Panel 2: EUr 7*10^5, NAm 15.9*10^6                #
#####################################################

setwd("/Users/nanditagarud/Documents/Jensen_response/analysis")

peaks=read.table('peaks_082712_400_50_2_2_0_pan_509_top50_rhoForCalling.txt')

pdf("fig3_BFheatmap_panels.pdf",width=2.28, height=2.43, title = "Figure 3",paper="special")
#width=6.83,height=4.86
#width=4.86,height=2.43

s=as.vector(c(1,2))
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

d=read.table("Admixture_mode_fixedPopSize_NAm1110000_EUr700000_rho_5e-9_BFs_numSim_100000.txt")
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

d=read.table("Admixture_mode_fixedPopSize_NAm15984500_EUr700000_rho_5e-9_BFs_numSim_100000.txt")
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

mtext(expression(H12),side=1,cex=cex2,padj=2.25, at=-4)


close.screen(all=TRUE)

dev.off() 


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
################################################################


setwd("/Users/nanditagarud/Documents/Ace/")

chr3R_509 = read.table('clusters_082712_400_50_2_2_0_chr3R_withRho_threshold.txt')
chr3L_509 = read.table('clusters_082712_400_50_2_2_0_chr3L_withRho_threshold.txt')
chr2L_509 = read.table('clusters_082712_400_50_2_2_0_chr2L_withRho_threshold.txt')
chr2R_509 = read.table('clusters_082712_400_50_2_2_0_chr2R_withRho_threshold.txt')
h=sort(c(chr3R_509[,28],chr3L_509[,28],chr2R_509[,28],chr2L_509[,28]))

peaks=read.table('peaks_082712_400_50_2_2_0_pan_509_top50.txt')


# grab points that are within 1SD around the median
median_val=median(h)
SD=0.341*median_val 
upper=median_val+SD
lower=median_val-SD

distr=h[h<upper]
standard_dev=sd(distr)

#hist(h, breaks=100, main='distribution of H12 scores', ylab='num windows', xlab='H12', col='black')
#abline(v=upper, col='red')
#abline(v=lower, col='red')
#abline(v=median_val+standard_dev, col='skyblue')
#abline(v=median_val-standard_dev, col='skyblue')
#abline(v=median_val+3*standard_dev, col='blue')
#abline(v=median_val-3*standard_dev, col='blue')
#abline(v=median_val+5*standard_dev, col='purple')

#abline(v=peaks[1,28], col='gold')
#abline(v=peaks[50,28], col='green')

# what fraction of points lie within 1 sd of the median?
w=h[lower<h]
sum(h<upper)/length(h)

# what fraction of points lie within 1sd of the fitted normal?
w=h[(median_val-1*standard_dev)<h]
w2=w[w<(median_val+1*standard_dev)]
length(w2)/length(h)

#############################
# A) plot of 
#qq plot of the w2 distr.   #
# main='Quantile-quantile plot assessing fit of inferred normal distribution to the bulk of the distribution'
#############################
pdf("Figure_QQ_norm.pdf",width=7,height= 4, title = "Figure 2",paper="special")

s=as.vector(c(1,2))
split.screen(s)

screen(1)

random_norm=rnorm(length(w2), mean(w2), sd(w2))
qqplot(w2, random_norm, ylab='Simulated Normal',xlab='DGRP data', main='',axes=F,cex.lab=0.8)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

abline(0, 1, col = 'red')
mtext(expression(A),at=0.011,adj=0.02,cex=1.5)

#######################################
# B) Plot the distribution of H12 values for DGRP data
# demarcate where the top 50 points lie
# Also demarcate where 11SD lies
# overlay iwth simulated normal 
#################################################

screen(2)

dgrp=hist(h, breaks=seq(0, 1, length.out = 301), plot=F)
dgrp$counts=dgrp$counts/sum(dgrp$counts)
plot(dgrp, col=rgb(1,0,0,0.5),ylab='Fraction of analysis windows',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F, xlab='H12')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sample=rnorm(length(h), mean=median_val, sd=standard_dev)
sims=hist(sample, breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

for (i in seq(1,50))
{
points(peaks[i,28],0, col='red', pch=20)
}

abline(v=median_val+11*standard_dev, col='blue') # plot where 11 SD is

legend(0.08, 0.8, legend=c("DGRP data", "Simulated Normal"),
       col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch=15, cex=0.8)
mtext(expression(B),at=-0.08,padj=0,cex=1.5)

close.screen(all=TRUE)

dev.off() 

##########




dgrp=hist(h, breaks=seq(0, 1, length.out = 301), plot=F)
dgrp$counts=dgrp$counts/sum(dgrp$counts)
###############



###############
pdf("Figure_supp_tmp.pdf",width=8,height= 6, title = "Figure 2",paper="special")

s=as.vector(c(2,3))
split.screen(s)

change=0.15

mb1=0.8 
ml1=0.8
mt1=0.3 
mr1=0.1 

mb2=0.8 
ml2=0.75 
mt2=0.3 
mr2=0.15

mb3=0.8 
ml3=0.7	 
mt3=0.3 
mr3=0.1 


screen(1)

par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('constNe10e6_neutrality_locusLen100000_snps.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='Fraction of analysis windows', xlab='',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[1:length(h),9], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

for (i in seq(1,50))
{
points(peaks[i,28],0, col='red', pch=20)
}

mtext(expression(A),at=-0.08,padj=0,cex=1.5)


legend(0.08, 0.8, legend=c("DGRP data", "Simulations"),
       col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch=15, cex=0.8)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('constNe2.7e6_neutrality_locusLen100000_snps.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[1:length(h),9], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

for (i in seq(1,50))
{
points(peaks[i,28],0, col='red', pch=20)
}

mtext(expression(B),at=-0.08,padj=0,cex=1.5)


screen(3)
par(mai=c(mb3,ml3,mt3,mr3))
d=read.table('dadi1_neutrality_locusLen100000_snps.txt')

plot(dgrp, col=rgb(1,0,0,0.5), xlab='',ylab='', main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[1:length(h),9], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

for (i in seq(1,50))
{
points(peaks[i,28],0, col='red', pch=20)
}
#mtext(expression(C),at=-0.1,padj=-2,cex=0.9)
mtext(expression(C),at=-0.08,padj=0,cex=1.5)


screen(4)
par(mai=c(mb1,ml1,mt1,mr1))

d=read.table('dadi2_neutrality_locusLen100000_snps.txt')

plot(dgrp, col=rgb(1,0,0,0.5), xlab='H12',ylab='Fraction of analysis windows', main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[1:length(h),9], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T)

for (i in seq(1,50))
{
points(peaks[i,28],0, col='red', pch=20)
}


mtext(expression(D),at=-0.08,padj=0,cex=1.5)


screen(5)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('admixture_mode_Garud2015_neutrality_locusLen100000_snps.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='H12',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[1:length(h),9], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

for (i in seq(1,50))
{
points(peaks[i,28],0, col='red', pch=20)
}

mtext(expression(E),at=-0.08,padj=0,cex=1.5)




screen(6)
par(mai=c(mb3,ml3,mt3,mr3))

d=read.table('admixture_bot_mode_Garud2015_neutrality_locusLen100000_snps.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='H12',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[1:length(h),9], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

for (i in seq(1,50))
{
points(peaks[i,28],0, col='red', pch=20)
}


mtext(expression(F),at=-0.08,padj=0,cex=1.5)


close.screen(all=TRUE)

dev.off() 


##########################


###############
pdf("Figure_supp_tmp_2.pdf",width=8,height= 6, title = "Figure 2",paper="special")

s=as.vector(c(2,3))
split.screen(s)

change=0.15

mb1=0.8 
ml1=0.8
mt1=0.3 
mr1=0.1 

mb2=0.8 
ml2=0.75 
mt2=0.3 
mr2=0.15

mb3=0.8 
ml3=0.7	 
mt3=0.3 
mr3=0.1 


screen(1)

par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('admixture_mode_neutrality_locusLen100000_snps.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='Fraction of analysis windows', xlab='',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[1:length(h),9], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

for (i in seq(1,50))
{
points(peaks[i,28],0, col='red', pch=20)
}

mtext(expression(A),at=-0.08,padj=0,cex=1.5)


legend(0.08, 0.8, legend=c("DGRP data", "Simulations"),
       col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch=15, cex=0.8)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('admixture_bot_mode_neutrality_locusLen100000_snps.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[1:length(h),9], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

for (i in seq(1,50))
{
points(peaks[i,28],0, col='red', pch=20)
}

mtext(expression(B),at=-0.08,padj=0,cex=1.5)


screen(3)
par(mai=c(mb3,ml3,mt3,mr3))
d=read.table('admixture_posterior_neutrality_locusLen100000_snps.txt')

plot(dgrp, col=rgb(1,0,0,0.5), xlab='',ylab='', main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[1:length(h),9], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

for (i in seq(1,50))
{
points(peaks[i,28],0, col='red', pch=20)
}
#mtext(expression(C),at=-0.1,padj=-2,cex=0.9)
mtext(expression(C),at=-0.08,padj=0,cex=1.5)


screen(4)
par(mai=c(mb1,ml1,mt1,mr1))

d=read.table('admixture_posterior_Harris_neutrality_locusLen100000_snps.txt')

plot(dgrp, col=rgb(1,0,0,0.5), xlab='H12',ylab='Fraction of analysis windows', main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[1:length(h),9], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T)

for (i in seq(1,50))
{
points(peaks[i,28],0, col='red', pch=20)
}


mtext(expression(D),at=-0.08,padj=0,cex=1.5)


screen(5)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('Arguello_neutrality_locusLen100000_snps.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='H12',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[1:length(h),9], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

for (i in seq(1,50))
{
points(peaks[i,28],0, col='red', pch=20)
}

mtext(expression(E),at=-0.08,padj=0,cex=1.5)



close.screen(all=TRUE)

dev.off() 


#####
###############
pdf("Figure_supp_tmp_3.pdf",width=8,height= 6, title = "Figure 2",paper="special")

s=as.vector(c(2,3))
split.screen(s)


screen(1)

par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('Admixture_mode_fixedPopSize_NAm1110000_EUr700000_rho_5e-9_theta_0_selection_False_MS_snps.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='Fraction of analysis windows', xlab='',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[1:length(h),9], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

for (i in seq(1,50))
{
points(peaks[i,28],0, col='red', pch=20)
}

mtext(expression(A),at=-0.08,padj=0,cex=1.5)


legend(0.08, 0.8, legend=c("DGRP data", "Simulations"),
       col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch=15, cex=0.8)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('Admixture_mode_fixedPopSize_NAm15984500_EUr700000_rho_5e-9_theta_0_selection_False_MS_snps.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[1:length(h),9], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

for (i in seq(1,50))
{
points(peaks[i,28],0, col='red', pch=20)
}

mtext(expression(B),at=-0.08,padj=0,cex=1.5)



close.screen(all=TRUE)

dev.off() 

