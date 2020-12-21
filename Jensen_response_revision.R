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

pdf("Figure6.pdf",width=7,height=3, title = "Figure8",paper="special")
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
mtext(H12,side=1,cex=cex2,at=30, padj=-2)
mtext(H2H1,side=1,cex=cex2,at=30, padj=-0.5)
#mtext(expression(A),at=-20,padj=0,cex=0.9)
mtext("H12:",side=1,at=-14,padj=-2, cex=cex2)
mtext("H2/H1:",side=1,at=-14,padj=-0.5, cex=cex2)
mtext("Cyp6g1",at=20,padj=0, cex=cex2,font=3)

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
mtext(H12,side=1,cex=cex2,padj=-2)
mtext(H2H1,side=1,cex=cex2,padj=-0.5)
mtext("CHKov1",at=27,padj=0, cex=cex2,font=3)


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
mtext(H12,side=1,cex=cex2,padj=-2)
mtext(H2H1,side=1,cex=cex2,padj=-0.5)
mtext("Ace",at=25,padj=0, cex=cex2,font=3)

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
mtext(H12,side=1,cex=cex2,padj=-2)
mtext(H2H1,side=1,cex=cex2,padj=-0.5)
}

close.screen(all=TRUE)
dev.off() 


#####################################################################
# Figure S16
# plot the haplotype spectra for the rest of the top 50 peaks
#####################################################################

peaks=read.table('peaks_082712_400_50_2_2_0_pan_509_top50.txt')


pdf("Figure_S16_a.pdf",width=10,height=3, title = "HapDist",paper="special")
s=as.vector(c(1,20))
split.screen(s)
par(omi=c(0.1,0.30,0.15,0.05)) #A vector of the form c(bottom, left, top, right) giving the size of the outer margins in inches.

cex2=0.6
#for (x in 31:50){
for (x in 11:30){
screen(x-10)
#screen(x-30)
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
#median_val=median(h)
#SD=0.341*median_val 
#upper=median_val+SD
#lower=median_val-SD

#distr=h[h<upper]
#distr=distr[distr>lower]
#standard_dev=sd(distr)

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
#w=h[lower<h]
#sum(h<upper)/length(h)

# what fraction of points lie within 1sd of the fitted normal?
#w=h[(median_val-1*standard_dev)<h]
#w2=w[w<(median_val+1*standard_dev)]
#length(w2)/length(h)

bulk_distr <- function(h) {
   median_val=median(h)
	SD=0.341*median_val 
	upper=median_val+SD
	lower=median_val-SD
	distr=h[h<upper]
	distr=distr[distr>lower]
	# having obtained this bulk, compute the standard deve
	standard_dev=sd(distr)	
	# now grab the 1SD of points using htis new standard dev
	w=h[(median_val-1*standard_dev)<h]
	new_distr=w[w<(median_val+1*standard_dev)]
	summary=list('new_distr'=new_distr, 'median_val'=median_val, 'standard_dev'= standard_dev)
     return(summary)	
}


#############################
# A) plot of 
#qq plot of the w2 distr.   #
# main='Quantile-quantile plot assessing fit of inferred normal distribution to the bulk of the distribution'
#############################
pdf("Figure_QQ_norm.pdf",width=7,height= 7, title = "Figure 2",paper="special")

summary=bulk_distr(h)
w2=summary$new_distr
median_val=summary$median_val
standard_dev=summary$standard_dev

s=as.vector(c(2,2))
split.screen(s)

screen(1)

random_norm=rnorm(length(w2), mean(w2), sd(w2))
#random_norm=rnorm(length(w2), median_val, standard_dev)
qqplot(random_norm, w2, xlab='Simulated Normal',ylab='DGRP data (bulk)', main='',axes=F,cex.lab=0.8)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

abline(0, 1, col = 'red')
mtext(expression(A),at=0.005,adj=0.02,cex=1.5)

#######################################
# B) Plot the distribution of H12 values for DGRP data
# demarcate where the top 50 points lie
# Also demarcate where 11SD lies
# overlay iwth simulated normal 
#################################################

screen(2)

dgrp=hist(h, breaks=seq(0, 1, length.out = 301), plot=F)
dgrp$counts=dgrp$counts/sum(dgrp$counts)
plot(dgrp, col=rgb(1,0,0,0.5),ylab='Fraction of windows',main='',ylim=c(0,0.8), xlim=c(0,0.25),cex.lab=0.8,axes=F, xlab='H12')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sample=rnorm(length(h), mean(w2), sd(w2))
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



####################################################
# QQ plot for the entire DGRP distribution
#####################################################

screen(3)

random_norm=rnorm(length(h), mean(w2), sd(w2))
qqplot(random_norm, h, xlab='Simulated Normal',ylab='DGRP data', main='',axes=F,cex.lab=0.8, ylim=c(0,0.25))
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

abline(0, 1, col = 'red')
mtext(expression(C),at=0.0055,adj=0.02,cex=1.5)



###################################################
# QQ plot for a simulated neutral Const Ne distribution
###################################################

screen(4)

d=read.table('/Users/nanditagarud/Documents/Jensen_response/analysis/constNe2.7e6_neutrality_locusLen100000_snps.txt')

summary=bulk_distr(d[,9])
w2=summary$new_distr
median_val=summary$median_val
standard_dev=summary$standard_dev
random_norm=rnorm(length(d[,9]), mean(w2), sd(w2))

qqplot(random_norm, d[,9], xlab='Simulated Normal',ylab='Constant Ne model', main='',axes=F,cex.lab=0.8, ylim=c(0,0.25))
abline(0, 1, col = 'red')
mtext(expression(D),at=0.0065,adj=0.02,cex=1.5)

axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)


close.screen(all=TRUE)

dev.off() 

############################

dgrp=hist(h, breaks=seq(0, 1, length.out = 301), plot=F)
dgrp$counts=dgrp$counts/sum(dgrp$counts)
###############



###############
setwd('/Users/nanditagarud/Documents/Jensen_response/analysis')
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



#############################
# Plot S and Pi genome-wide #
#############################



setwd('/Users/nanditagarud/Documents/Jensen_response/analysis/S_and_Pi')
chr2L=read.table('Chr2L_SNPdensity.txt')
chr2R=read.table('Chr2R_SNPdensity.txt')
chr3L=read.table('Chr3L_SNPdensity.txt')
chr3R=read.table('Chr3R_SNPdensity.txt')

pdf("Fig_S_scan.pdf",width=6.83,height=4, title = "Figure 6a",paper="special")

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

plot(chr2L[,1], chr2L[,2], ylim=c(0,0.1), xlim=c(0,chr2L[length(chr2L[,1]),1] +shift), col="black",pch=20, xaxs="i",yaxs = "i",axes=F, main='',cex=0.4)


axis(1,at=c(0,chr2L[length(chr2L[,1]),1]+shift),labels=F,tcl=0)
axis(1,at=seq(0,chr2L[length(chr2L[,1]),1],5*10^6),labels=seq(0,chr2L[length(chr2L[,1]),1],5*10^6),padj=pdjx,tcl=-tcl)
axis(2,at=c(0,0.27),labels=F,tcl=0)
axis(2,at=seq(0,0.25,0.05),labels=seq(0,0.25,0.05),padj=pdjy,tcl=-tcl)
axis(3,at=c(0,chr2L[length(chr2L[,1]),1]+shift),labels=F,tcl=0)
axis(4,at=c(0,0.27),labels=F,tcl=0)

mtext("Chr2L", at=c(1.5*10^6), cex=cex2, padj=2, font=2)
mtext(expression("S/bp"),side=2,cex=cex2,padj=-2.5)
#mtext(expression(A),at=-10^3,padj=-1,cex=cex3)


#Chr2R
screen(2)
par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size 

plot(chr2R[,1], chr2R[,2], ylim=c(0,0.1), xlim=c(0,chr2R[length(chr2R[,1]),1] +shift), col="black",pch=20, xaxs="i",yaxs = "i",axes=F, main='',cex=0.4)
abline(v= 8074551, col='red',lty=2)

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

plot(chr3L[,1], chr3L[,2], ylim=c(0,0.1), xlim=c(0,chr3L[length(chr3L[,1]),1] +shift), col="black",pch=20, xaxs="i",yaxs = "i",axes=F, main='',cex=0.4)

axis(1,at=c(0,chr3L[length(chr3L[,1]),1]+shift),labels=F,tcl=0)
axis(1,at=seq(0,chr3L[length(chr3L[,1]),1],5*10^6),labels=seq(0,chr3L[length(chr3L[,1]),1],5*10^6),padj=pdjx,tcl=-tcl)
axis(2,at=c(0,0.27),labels=F,tcl=0)
axis(2,at=seq(0,0.25,0.05),labels=seq(0,0.25,0.05),padj=pdjy,tcl=-tcl)
axis(3,at=c(0,chr3L[length(chr3L[,1]),1]+shift),labels=F,tcl=0)
axis(4,at=c(0,0.27),labels=F,tcl=0)

mtext("Chr3L", at=c(1.5*10^6), cex=cex2, padj=2, font=2)
mtext(expression("S/bp"),side=2,cex=cex2,padj=-2.5)
mtext(expression("position"),side=1,cex=cex2,padj=2)



#Chr3R
screen(4)
par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size 

plot(chr3R[,1], chr3R[,2], ylim=c(0,0.1), xlim=c(0,chr3R[length(chr3R[,1]),1] +shift), col="black",pch=20, xaxs="i",yaxs = "i",axes=F, main='',cex=0.4)
abline(v= 9069408, col='red',lty=2)
abline(v= 21152026, col='red',lty=2)

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



####################
# Pi/bp genomewide #
####################

setwd('/Users/nanditagarud/Documents/Jensen_response/analysis/S_and_Pi')
chr2L=read.table('Chr2L_SNPdensity.txt')
chr2R=read.table('Chr2R_SNPdensity.txt')
chr3L=read.table('Chr3L_SNPdensity.txt')
chr3R=read.table('Chr3R_SNPdensity.txt')

pdf("Fig_Pi_scan.pdf",width=6.83,height=4, title = "Figure 6a",paper="special")

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

plot(chr2L[,1], chr2L[,3], ylim=c(0,0.015), xlim=c(0,chr2L[length(chr2L[,1]),1] +shift), col="black",pch=20, xaxs="i",yaxs = "i",axes=F, main='',cex=0.4)


axis(1,at=c(0,chr2L[length(chr2L[,1]),1]+shift),labels=F,tcl=0)
axis(1,at=seq(0,chr2L[length(chr2L[,1]),1],5*10^6),labels=seq(0,chr2L[length(chr2L[,1]),1],5*10^6),padj=pdjx,tcl=-tcl)
axis(2,at=c(0,0.015),labels=F,tcl=0)
axis(2,at=seq(0,0.015,0.01),labels=seq(0,0.015,0.01),padj=pdjy,tcl=-tcl)
axis(3,at=c(0,chr2L[length(chr2L[,1]),1]+shift),labels=F,tcl=0)
axis(4,at=c(0,0.015),labels=F,tcl=0)

mtext("Chr2L", at=c(1.5*10^6), cex=cex2, padj=2, font=2)
mtext(expression("Pi/bp"),side=2,cex=cex2,padj=-2.5)
#mtext(expression(A),at=-10^3,padj=-1,cex=cex3)


#Chr2R
screen(2)
par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size 

plot(chr2R[,1], chr2R[,3], ylim=c(0,0.015), xlim=c(0,chr2R[length(chr2R[,1]),1] +shift), col="black",pch=20, xaxs="i",yaxs = "i",axes=F, main='',cex=0.4)
abline(v= 8074551, col='red',lty=2)

axis(1,at=c(0,chr2R[length(chr2R[,1]),1]+shift),labels=F,tcl=0)
axis(1,at=seq(0,chr2R[length(chr2R[,1]),1],5*10^6),labels=seq(0,chr2R[length(chr2R[,1]),1],5*10^6),padj=pdjx,tcl=-tcl)
axis(2,at=c(0,0.015),labels=F,tcl=0)
axis(2,at=seq(0,0.015,0.01),labels=seq(0,0.015,0.01),padj=pdjy,tcl=-tcl)
axis(3,at=c(0,chr2R[length(chr2R[,1]),1]+shift),labels=F,tcl=0)
axis(4,at=c(0,0.015),labels=F,tcl=0)

mtext("Chr2R", at=c(1.5*10^6), cex=cex2, padj=2, font=2)
mtext("Cyp6g1", at=c(8074551), cex=cex2, padj=1.8, font=3)


#Chr3L
screen(3)
par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size 

plot(chr3L[,1], chr3L[,3], ylim=c(0,0.015), xlim=c(0,chr3L[length(chr3L[,1]),1] +shift), col="black",pch=20, xaxs="i",yaxs = "i",axes=F, main='',cex=0.4)

axis(1,at=c(0,chr3L[length(chr3L[,1]),1]+shift),labels=F,tcl=0)
axis(1,at=seq(0,chr3L[length(chr3L[,1]),1],5*10^6),labels=seq(0,chr3L[length(chr3L[,1]),1],5*10^6),padj=pdjx,tcl=-tcl)
axis(2,at=c(0,0.015),labels=F,tcl=0)
axis(2,at=seq(0,0.015,0.01),labels=seq(0,0.015,0.01),padj=pdjy,tcl=-tcl)
axis(3,at=c(0,chr3L[length(chr3L[,1]),1]+shift),labels=F,tcl=0)
axis(4,at=c(0,0.015),labels=F,tcl=0)

mtext("Chr3L", at=c(1.5*10^6), cex=cex2, padj=2, font=2)
mtext(expression("Pi/bp"),side=2,cex=cex2,padj=-2.5)
mtext(expression("position"),side=1,cex=cex2,padj=2)



#Chr3R
screen(4)
par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size 

plot(chr3R[,1], chr3R[,3], ylim=c(0,0.015), xlim=c(0,chr3R[length(chr3R[,1]),1] +shift), col="black",pch=20, xaxs="i",yaxs = "i",axes=F, main='',cex=0.4)
abline(v= 9069408, col='red',lty=2)
abline(v= 21152026, col='red',lty=2)

axis(1,at=c(0,chr3R[length(chr3R[,1]),1]+shift),labels=F,tcl=0)
axis(1,at=seq(0,chr3R[length(chr3R[,1]),1],5*10^6),labels=seq(0,chr3R[length(chr3R[,1]),1],5*10^6),padj=pdjx,tcl=-tcl)
axis(2,at=c(0,0.015),labels=F,tcl=0)
axis(2,at=seq(0,0.015,0.01),labels=seq(0,0.015,0.01),padj=pdjy,tcl=-tcl)
axis(3,at=c(0,chr3R[length(chr3R[,1]),1]+shift),labels=F,tcl=0)
axis(4,at=c(0,0.015),labels=F,tcl=0)

mtext("Chr3R", at=c(2*10^6), cex=cex2, padj=2, font=2)
mtext(expression("position"),side=1,cex=cex2,padj=2)
mtext("Ace", at=c(9069408), cex=cex2, padj=2, font=3)
mtext("CHKov1", at=c(21152026), cex=cex2, padj=2, font=3)


close.screen(all=TRUE)

dev.off() 


########


par(mfrow=c(2,2))
plot(chr2L[,1], chr2L[,2], xlab='position',ylab='S/bp')
plot(chr2R[,1], chr2R[,2], xlab='position',ylab='S/bp')
abline(v= 8074551, col='red')
plot(chr3L[,1], chr3L[,2], xlab='position',ylab='S/bp')
plot(chr3R[,1], chr3R[,2], xlab='position',ylab='S/bp')
abline(v= 9069408, col='red')
abline(v= 21152026, col='red')


par(mfrow=c(2,2))
plot(chr2L[,1], chr2L[,3], xlab='position',ylab='Pi/bp')
plot(chr2R[,1], chr2R[,3], xlab='position',ylab='Pi/bp')
abline(v= 8074551, col='red')
plot(chr3L[,1], chr3L[,3], xlab='position',ylab='Pi/bp')
plot(chr3R[,1], chr3R[,3], xlab='position',ylab='Pi/bp')
abline(v= 9069408, col='red')
abline(v= 21152026, col='red')


#mtext("Cyp6g1", at=c(8074551), cex=cex2, padj=1.8, font=3)
#mtext("Ace", at=c(9069408), cex=cex2, padj=2, font=3)
#mtext("CHKov1", at=c(21152026), cex=cex2, padj=2, font=3)

#################
# Plot for short introns  (perbp)

chr2L=read.table('Chr2L_SNPdensity_shortIntron.txt')
chr2R=read.table('Chr2R_SNPdensity_shortIntron.txt')
chr3L=read.table('Chr3L_SNPdensity_shortIntron.txt')
chr3R=read.table('Chr3R_SNPdensity_shortIntron.txt')

sum(c(chr2L[,2],chr2R[,2], chr3L[,2],chr3R[,2]))/870364
sum(c(chr2L[,3],chr2R[,3], chr3L[,3],chr3R[,3]))/870364

# per short intron (use this -- can actually build a distribution out of this and also the means match dadi values better)

chr2L=read.table('Chr2L_SNPdensity_perShortIntron.txt')
chr2R=read.table('Chr2R_SNPdensity_perShortIntron.txt')
chr3L=read.table('Chr3L_SNPdensity_perShortIntron.txt')
chr3R=read.table('Chr3R_SNPdensity_perShortIntron.txt')

dgrp= (c(chr2L[,3],chr2R[,3], chr3L[,3],chr3R[,3])) # s/bp
mean((c(chr2L[,3],chr2R[,3], chr3L[,3],chr3R[,3])))
mean((c(chr2L[,4],chr2R[,4], chr3L[,4],chr3R[,4])))

hist((c(chr2L[,3],chr2R[,3], chr3L[,3],chr3R[,3])))
hist((c(chr2L[,4],chr2R[,4], chr3L[,4],chr3R[,4])))
dadi_mean_s=0.0578
dadi_mean_pi=0.0118


###################################################################
# Plot distributions of S for short introns vs simulations #
###################################################################
setwd('/Users/nanditagarud/Documents/Jensen_response/analysis/S_and_Pi')
chr2L=read.table('Chr2L_SNPdensity_perShortIntron.txt')
chr2R=read.table('Chr2R_SNPdensity_perShortIntron.txt')
chr3L=read.table('Chr3L_SNPdensity_perShortIntron.txt')
chr3R=read.table('Chr3R_SNPdensity_perShortIntron.txt')

dgrp= (c(chr2L[,3],chr2R[,3], chr3L[,3],chr3R[,3])) # s/bp
dgrp=hist(dgrp, breaks=seq(0, 1, length.out = 301), plot=F)
dgrp$counts=dgrp$counts/sum(dgrp$counts)


###############
setwd('/Users/nanditagarud/Documents/Jensen_response/analysis')
pdf("Figure_supp_tmp_S_1.pdf",width=8,height= 6, title = "Figure 2",paper="special")

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

d=read.table('constNe10e6_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='Fraction of analysis windows', xlab='',main='',ylim=c(0,0.8), xlim=c(0,0.3), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,1], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(A),at=-0.08,padj=0,cex=1.5)

legend(0.08, 0.8, legend=c("DGRP data", "Simulations"),
       col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch=15, cex=0.8)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('constNe2.7e6_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,1], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(B),at=-0.08,padj=0,cex=1.5)


screen(3)
par(mai=c(mb3,ml3,mt3,mr3))
d=read.table('dadi1_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5), xlab='',ylab='', main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,1], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(C),at=-0.08,padj=0,cex=1.5)


screen(4)
par(mai=c(mb1,ml1,mt1,mr1))

d=read.table('dadi2_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5), xlab='S/bp',ylab='Fraction of analysis windows', main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,1], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T)

mtext(expression(D),at=-0.08,padj=0,cex=1.5)


screen(5)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('admixture_mode_Garud2015_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='S/bp',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,1], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(E),at=-0.08,padj=0,cex=1.5)


screen(6)
par(mai=c(mb3,ml3,mt3,mr3))

d=read.table('admixture_bot_mode_Garud2015_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='S/bp',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,1], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(F),at=-0.08,padj=0,cex=1.5)

close.screen(all=TRUE)

dev.off() 


##########################


###############
pdf("Figure_supp_tmp_S_2.pdf",width=8,height= 6, title = "Figure 2",paper="special")

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

d=read.table('admixture_mode_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='Fraction of analysis windows', xlab='',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,1], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(A),at=-0.08,padj=0,cex=1.5)

legend(0.08, 0.8, legend=c("DGRP data", "Simulations"),
       col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch=15, cex=0.8)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('admixture_bot_mode_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,1], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(B),at=-0.08,padj=0,cex=1.5)


screen(3)
par(mai=c(mb3,ml3,mt3,mr3))
d=read.table('admixture_posterior_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5), xlab='',ylab='', main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,1], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(C),at=-0.08,padj=0,cex=1.5)


screen(4)
par(mai=c(mb1,ml1,mt1,mr1))

d=read.table('admixture_posterior_Harris_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5), xlab='S/bp',ylab='Fraction of analysis windows', main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,1], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T)

mtext(expression(D),at=-0.08,padj=0,cex=1.5)


screen(5)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('Arguello_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='S/bp',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,1], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)


mtext(expression(E),at=-0.08,padj=0,cex=1.5)



close.screen(all=TRUE)

dev.off() 


#####
###############
pdf("Figure_supp_tmp_S_3.pdf",width=8,height= 6, title = "Figure 2",paper="special")

s=as.vector(c(2,3))
split.screen(s)


screen(1)

par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('Admixture_mode_fixedPopSize_NAm1110000_EUr700000_rho_5e-9_theta_0_selection_False_MS_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='Fraction of analysis windows', xlab='S/bp',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,1], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(A),at=-0.08,padj=0,cex=1.5)


legend(0.08, 0.8, legend=c("DGRP data", "Simulations"),
       col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch=15, cex=0.8)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('Admixture_mode_fixedPopSize_NAm15984500_EUr700000_rho_5e-9_theta_0_selection_False_MS_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='S/bp',main='',ylim=c(0,0.8), xlim=c(0,0.25), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,1], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(B),at=-0.08,padj=0,cex=1.5)



close.screen(all=TRUE)

dev.off() 



#########
# Pi/bp #
#########

setwd('/Users/nanditagarud/Documents/Jensen_response/analysis/S_and_Pi')
chr2L=read.table('Chr2L_SNPdensity_perShortIntron.txt')
chr2R=read.table('Chr2R_SNPdensity_perShortIntron.txt')
chr3L=read.table('Chr3L_SNPdensity_perShortIntron.txt')
chr3R=read.table('Chr3R_SNPdensity_perShortIntron.txt')

dgrp= (c(chr2L[,4],chr2R[,4], chr3L[,4],chr3R[,4])) # pi/bp
dgrp=hist(dgrp, breaks=seq(0, 1, length.out = 301), plot=F)
dgrp$counts=dgrp$counts/sum(dgrp$counts)


i=2 #column with pi
lab_pos=-0.04 #position of letter label

setwd('/Users/nanditagarud/Documents/Jensen_response/analysis')
pdf("Figure_supp_tmp_Pi_1.pdf",width=8,height= 6, title = "Figure 2",paper="special")

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

d=read.table('constNe10e6_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='Fraction of analysis windows', xlab='',main='',ylim=c(0,0.8), xlim=c(0,0.1), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,i], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(A),at=lab_pos,padj=0,cex=1.5)

legend(0.03, 0.8, legend=c("DGRP data", "Simulations"),
       col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch=15, cex=0.8)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('constNe2.7e6_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='',main='',ylim=c(0,0.8), xlim=c(0,0.1), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,i], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(B),at=lab_pos,padj=0,cex=1.5)


screen(3)
par(mai=c(mb3,ml3,mt3,mr3))
d=read.table('dadi1_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5), xlab='',ylab='', main='',ylim=c(0,0.8), xlim=c(0,0.1), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,i], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(C),at=lab_pos,padj=0,cex=1.5)


screen(4)
par(mai=c(mb1,ml1,mt1,mr1))

d=read.table('dadi2_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5), xlab='Pi/bp',ylab='Fraction of analysis windows', main='',ylim=c(0,0.8), xlim=c(0,0.1), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,i], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T)

mtext(expression(D),at=lab_pos,padj=0,cex=1.5)


screen(5)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('admixture_mode_Garud2015_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='Pi/bp',main='',ylim=c(0,0.8), xlim=c(0,0.1), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,i], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(E),at=lab_pos,padj=0,cex=1.5)


screen(6)
par(mai=c(mb3,ml3,mt3,mr3))

d=read.table('admixture_bot_mode_Garud2015_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='Pi/bp',main='',ylim=c(0,0.8), xlim=c(0,0.1), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,i], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(F),at=lab_pos,padj=0,cex=1.5)

close.screen(all=TRUE)

dev.off() 


##########################


###############
pdf("Figure_supp_tmp_Pi_2.pdf",width=8,height= 6, title = "Figure 2",paper="special")

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

d=read.table('admixture_mode_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='Fraction of analysis windows', xlab='',main='',ylim=c(0,0.8), xlim=c(0,0.1), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,i], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(A),at=lab_pos,padj=0,cex=1.5)

legend(0.03, 0.8, legend=c("DGRP data", "Simulations"),
       col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch=15, cex=0.8)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('admixture_bot_mode_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='',main='',ylim=c(0,0.8), xlim=c(0,0.1), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,i], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(B),at=lab_pos,padj=0,cex=1.5)


screen(3)
par(mai=c(mb3,ml3,mt3,mr3))
d=read.table('admixture_posterior_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5), xlab='',ylab='', main='',ylim=c(0,0.8), xlim=c(0,0.1), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,i], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(C),at=lab_pos,padj=0,cex=1.5)


screen(4)
par(mai=c(mb1,ml1,mt1,mr1))

d=read.table('admixture_posterior_Harris_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5), xlab='Pi/bp',ylab='Fraction of analysis windows', main='',ylim=c(0,0.8), xlim=c(0,0.1), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,i], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T)

mtext(expression(D),at=lab_pos,padj=0,cex=1.5)


screen(5)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('Arguello_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='Pi/bp',main='',ylim=c(0,0.8), xlim=c(0,0.1), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,i], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)


mtext(expression(E),at=lab_pos,padj=0,cex=1.5)


close.screen(all=TRUE)

dev.off() 


#####
###############
pdf("Figure_supp_tmp_Pi_3.pdf",width=8,height= 6, title = "Figure 2",paper="special")

s=as.vector(c(2,3))
split.screen(s)


screen(1)

par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('Admixture_mode_fixedPopSize_NAm1110000_EUr700000_rho_5e-9_theta_0_selection_False_MS_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='Fraction of analysis windows', xlab='Pi/bp',main='',ylim=c(0,0.8), xlim=c(0,0.1), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,i], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(A),at=lab_pos,padj=0,cex=1.5)


legend(0.03, 0.8, legend=c("DGRP data", "Simulations"),
       col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch=15, cex=0.8)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('Admixture_mode_fixedPopSize_NAm15984500_EUr700000_rho_5e-9_theta_0_selection_False_MS_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(dgrp, col=rgb(1,0,0,0.5),ylab='', xlab='Pi/bp',main='',ylim=c(0,0.8), xlim=c(0,0.1), cex.lab=0.8,axes=F)
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

sims=hist(d[,i], breaks=seq(0, 1, length.out = 301), plot=F)
sims$counts=sims$counts/sum(sims$counts)
plot(sims,col=rgb(0,0,1,0.5), add=T,cex.axes=2)

mtext(expression(B),at=lab_pos,padj=0,cex=1.5)



close.screen(all=TRUE)

dev.off() 







######################################
# test of fit of data to simulations #
######################################
setwd("/Users/nanditagarud/Documents/Ace/")
chr3R_509 = read.table('clusters_082712_400_50_2_2_0_chr3R_withRho_threshold.txt')
chr3L_509 = read.table('clusters_082712_400_50_2_2_0_chr3L_withRho_threshold.txt')
chr2L_509 = read.table('clusters_082712_400_50_2_2_0_chr2L_withRho_threshold.txt')
chr2R_509 = read.table('clusters_082712_400_50_2_2_0_chr2R_withRho_threshold.txt')
dgrp_h12=sort(c(chr3R_509[,28],chr3L_509[,28],chr2R_509[,28],chr2L_509[,28]))

setwd('/Users/nanditagarud/Documents/Jensen_response/analysis/S_and_Pi')
chr2L=read.table('Chr2L_SNPdensity_perShortIntron.txt')
chr2R=read.table('Chr2R_SNPdensity_perShortIntron.txt')
chr3L=read.table('Chr3L_SNPdensity_perShortIntron.txt')
chr3R=read.table('Chr3R_SNPdensity_perShortIntron.txt')
dgrp_s= (c(chr2L[,3],chr2R[,3], chr3L[,3],chr3R[,3])) # s/bp
dgrp_pi= (c(chr2L[,4],chr2R[,4], chr3L[,4],chr3R[,4])) # pi/bp


setwd('/Users/nanditagarud/Documents/Jensen_response/analysis')

ks_h12=c()
wilcox_h12=c()

for (name in c('constNe10e6','constNe2.7e6','dadi1','dadi2','admixture_mode_Garud2015','admixture_bot_mode_Garud2015','admixture_mode', 'admixture_bot_mode','admixture_posterior','admixture_posterior_Harris','Arguello')) {
	file=paste(name, '_neutrality_locusLen100000_snps.txt', sep = "")
	d=read.table(file)
	ks_h12 = c(ks_h12,ks.test(dgrp_h12,d[,9])$p.value)
    wilcox_h12=c(wilcox_h12, wilcox.test(dgrp_h12,d[,9])$p.value)
}

d=read.table('Admixture_mode_fixedPopSize_NAm1110000_EUr700000_rho_5e-9_theta_0_selection_False_MS_snps.txt')
ks_h12 = c(ks_h12,ks.test(dgrp_h12,d[,9])$p.value)
wilcox_h12=c(wilcox_h12, wilcox.test(dgrp_h12,d[,9])$p.value)

d=read.table('Admixture_mode_fixedPopSize_NAm15984500_EUr700000_rho_5e-9_theta_0_selection_False_MS_snps.txt')
ks_h12 = c(ks_h12,ks.test(dgrp_h12,d[,9])$p.value)
wilcox_h12=c(wilcox_h12, wilcox.test(dgrp_h12,d[,9])$p.value)


ks_pi=c()
wilcox_pi=c()
ks_s=c()
wilcox_s=c()

for (name in c('constNe10e6','constNe2.7e6','dadi1','dadi2','admixture_mode_Garud2015','admixture_bot_mode_Garud2015','admixture_mode', 'admixture_bot_mode','admixture_posterior','admixture_posterior_Harris','Arguello')) {
	file=paste(name, '_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt', sep = "")
	d=read.table(file)
	
	ks_s = c(ks_s,ks.test(dgrp_s,d[,1])$p.value)
    wilcox_s=c(wilcox_s, wilcox.test(dgrp_s,d[,1])$p.value)

    ks_pi = c(ks_pi,ks.test(dgrp_pi,d[,2])$p.value)
    wilcox_pi=c(wilcox_pi, wilcox.test(dgrp_pi,d[,2])$p.value)

}

d=read.table('Admixture_mode_fixedPopSize_NAm1110000_EUr700000_rho_5e-9_theta_0_selection_False_MS_shortIntron_Pi_S_TajD_perBp_4Fields.txt')
ks_s = c(ks_s,ks.test(dgrp_s,d[,1])$p.value)
wilcox_s=c(wilcox_s, wilcox.test(dgrp_s,d[,1])$p.value)

ks_pi = c(ks_pi,ks.test(dgrp_pi,d[,2])$p.value)
wilcox_pi=c(wilcox_pi, wilcox.test(dgrp_pi,d[,2])$p.value)

d=read.table('Admixture_mode_fixedPopSize_NAm15984500_EUr700000_rho_5e-9_theta_0_selection_False_MS_shortIntron_Pi_S_TajD_perBp_4Fields.txt')
ks_s = c(ks_s,ks.test(dgrp_s,d[,1])$p.value)
wilcox_s=c(wilcox_s, wilcox.test(dgrp_s,d[,1])$p.value)

ks_pi = c(ks_pi,ks.test(dgrp_pi,d[,2])$p.value)
wilcox_pi=c(wilcox_pi, wilcox.test(dgrp_pi,d[,2])$p.value)



################################################################
# Plot QQ-plots for fit of observed vs simulated distributions #
################################################################


# function for rmse
library(Metrics)

##########
# H12    #
##########
setwd("/Users/nanditagarud/Documents/Ace/")

chr3R_509 = read.table('clusters_082712_400_50_2_2_0_chr3R_withRho_threshold.txt')
chr3L_509 = read.table('clusters_082712_400_50_2_2_0_chr3L_withRho_threshold.txt')
chr2L_509 = read.table('clusters_082712_400_50_2_2_0_chr2L_withRho_threshold.txt')
chr2R_509 = read.table('clusters_082712_400_50_2_2_0_chr2R_withRho_threshold.txt')
h=sort(c(chr3R_509[,28],chr3L_509[,28],chr2R_509[,28],chr2L_509[,28]))


# grab points that are within 1SD around the median

bulk_distr <- function(h) {
   median_val=median(h)
	SD=0.341*median_val 
	upper=median_val+SD
	lower=median_val-SD
	distr=h[h<upper]
	distr=distr[distr>lower]
	standard_dev=sd(distr)
	w=h[(median_val-1*standard_dev)<h]
	new_distr=w[w<(median_val+1*standard_dev)]
     return(new_distr)	
}

dgrp= bulk_distr(h)

peaks=read.table('peaks_082712_400_50_2_2_0_pan_509_top50.txt')


lab_pos=-0.04 #position of letter label
i=9 # column for plotting H12

setwd('/Users/nanditagarud/Documents/Jensen_response/analysis')
pdf("Figure_supp_tmp_H12_QQ_1.pdf",width=8,height= 6, title = "Figure 2",paper="special")

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
new_distr=bulk_distr(d[,i])

w=qqplot(d[,i], h, ylab='DGRP', xlab='',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)


mtext(expression(A),at=0.004,padj=0,cex=1.5)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('constNe2.7e6_neutrality_locusLen100000_snps.txt')
new_distr=bulk_distr(d[,i])

w=qqplot(d[,i], h, ylab='', xlab='',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(B),at=0.004,padj=0,cex=1.5)


screen(3)
par(mai=c(mb3,ml3,mt3,mr3))
d=read.table('dadi1_neutrality_locusLen100000_snps.txt')
new_distr=bulk_distr(d[,i])

w=qqplot(d[,i], h,  xlab='',ylab='', main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(C),at=0.004,padj=0,cex=1.5)


screen(4)
par(mai=c(mb1,ml1,mt1,mr1))

d=read.table('dadi2_neutrality_locusLen100000_snps.txt')
new_distr=bulk_distr(d[,i])

w=qqplot(d[,i], h, xlab='Simulations',ylab='DGRP', main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(D),at=0.004,padj=0,cex=1.5)


screen(5)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('admixture_mode_Garud2015_neutrality_locusLen100000_snps.txt')
new_distr=bulk_distr(d[,i])

w=qqplot(d[,i], h, ylab='', xlab='Simulations',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)


mtext(expression(E),at=0.004,padj=0,cex=1.5)


screen(6)
par(mai=c(mb3,ml3,mt3,mr3))

d=read.table('admixture_bot_mode_Garud2015_neutrality_locusLen100000_snps.txt')
new_distr=bulk_distr(d[,i])

w=qqplot(d[,i], h, ylab='', xlab='Simulations',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)


mtext(expression(F),at=0.004,padj=0,cex=1.5)

close.screen(all=TRUE)

dev.off() 

###############

pdf("Figure_supp_tmp_H12_QQ_2.pdf",width=8,height= 6, title = "Figure 2",paper="special")

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
new_distr=bulk_distr(d[,i])

w=qqplot(d[,i], h,ylab='DGRP', xlab='',main='',cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.1,padj=0,cex=1)


mtext(expression(A),at=0,padj=0,cex=1.5)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('admixture_bot_mode_neutrality_locusLen100000_snps.txt')
new_distr=bulk_distr(d[,i])

w=qqplot(d[,i], h,ylab='', xlab='',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.1,padj=0,cex=1)


mtext(expression(B),at=-0.05,padj=0,cex=1.5)


screen(3)
par(mai=c(mb3,ml3,mt3,mr3))
d=read.table('admixture_posterior_neutrality_locusLen100000_snps.txt')
new_distr=bulk_distr(d[,i])

w=qqplot(d[,i], h, xlab='',ylab='', main='',cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.1,padj=0,cex=1)

mtext(expression(C),at=-0.15,padj=0,cex=1.5)


screen(4)
par(mai=c(mb1,ml1,mt1,mr1))

d=read.table('admixture_posterior_Harris_neutrality_locusLen100000_snps.txt')
new_distr=bulk_distr(d[,i])

w=qqplot(d[,i], h, xlab='Simulations',ylab='DGRP', main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.1,padj=0,cex=1)


mtext(expression(D),at=-0.15,padj=0,cex=1.5)


screen(5)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('Arguello_neutrality_locusLen100000_snps.txt')
new_distr=bulk_distr(d[,i])

w=qqplot(d[,i], h, ylab='', xlab='Simulations',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.05,padj=0,cex=1)


mtext(expression(E),at=0,padj=0,cex=1.5)


close.screen(all=TRUE)

dev.off() 


#####
###############
pdf("Figure_supp_tmp_H12_QQ_3.pdf",width=8,height= 6, title = "Figure 2",paper="special")

s=as.vector(c(2,3))
split.screen(s)


screen(1)

par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('Admixture_mode_fixedPopSize_NAm1110000_EUr700000_rho_5e-9_theta_0_selection_False_MS_snps.txt')
new_distr=bulk_distr(d[,i])

w=qqplot(d[,i], h, ylab='DGRP', xlab='Simulations',main='',cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(A),at= 0,padj=0.008,cex=1.5)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('Admixture_mode_fixedPopSize_NAm15984500_EUr700000_rho_5e-9_theta_0_selection_False_MS_snps.txt')
new_distr=bulk_distr(d[,i])

w=qqplot(d[,i], h, ylab='', xlab='Simulations',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)


mtext(expression(B),at= 0,padj=0.008,cex=1.5)



close.screen(all=TRUE)

dev.off() 







##########
# S/bp  #
##########

setwd('/Users/nanditagarud/Documents/Jensen_response/analysis/S_and_Pi')
chr2L=read.table('Chr2L_SNPdensity_perShortIntron.txt')
chr2R=read.table('Chr2R_SNPdensity_perShortIntron.txt')
chr3L=read.table('Chr3L_SNPdensity_perShortIntron.txt')
chr3R=read.table('Chr3R_SNPdensity_perShortIntron.txt')

dgrp= (c(chr2L[,3],chr2R[,3], chr3L[,3],chr3R[,3])) # S/bp


i=1 #column with s
lab_pos=-0.1 #position of letter label

setwd('/Users/nanditagarud/Documents/Jensen_response/analysis')
pdf("Figure_supp_tmp_S_QQ_1.pdf",width=8,height= 6, title = "Figure 2",paper="special")

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

d=read.table('constNe10e6_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i], dgrp, ylab='DGRP', xlab='',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(A),at=lab_pos,padj=0,cex=1.5)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('constNe2.7e6_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i], dgrp, ylab='', xlab='',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(B),at=lab_pos,padj=0,cex=1.5)


screen(3)
par(mai=c(mb3,ml3,mt3,mr3))
d=read.table('dadi1_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i], dgrp,  xlab='',ylab='', main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(C),at=lab_pos,padj=0,cex=1.5)


screen(4)
par(mai=c(mb1,ml1,mt1,mr1))

d=read.table('dadi2_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i], dgrp, xlab='Simulations',ylab='DGRP', main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(D),at=-0.2,padj=0,cex=1.5)


screen(5)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('admixture_mode_Garud2015_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i], dgrp, ylab='', xlab='Simulations',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(E),at=-0.25,padj=0,cex=1.5)


screen(6)
par(mai=c(mb3,ml3,mt3,mr3))

d=read.table('admixture_bot_mode_Garud2015_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i], dgrp, ylab='', xlab='Simulations',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(F),at=lab_pos,padj=0,cex=1.5)

close.screen(all=TRUE)

dev.off() 

###############

pdf("Figure_supp_tmp_S_QQ_2.pdf",width=8,height= 6, title = "Figure 2",paper="special")

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

d=read.table('admixture_mode_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i],dgrp,ylab='DGRP', xlab='',main='',cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(A),at=lab_pos,padj=0,cex=1.5)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('admixture_bot_mode_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i],dgrp,ylab='', xlab='',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(B),at=lab_pos,padj=0,cex=1.5)


screen(3)
par(mai=c(mb3,ml3,mt3,mr3))
d=read.table('admixture_posterior_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i],dgrp, xlab='',ylab='', main='',cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(C),at=lab_pos,padj=0,cex=1.5)


screen(4)
par(mai=c(mb1,ml1,mt1,mr1))

d=read.table('admixture_posterior_Harris_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i],dgrp, xlab='Simulations',ylab='DGRP', main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(D),at=lab_pos,padj=0,cex=1.5)


screen(5)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('Arguello_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i],dgrp, ylab='', xlab='Simulations',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(E),at=lab_pos,padj=0,cex=1.5)


close.screen(all=TRUE)

dev.off() 


#####
###############
pdf("Figure_supp_tmp_S_QQ_3.pdf",width=8,height= 6, title = "Figure 2",paper="special")

s=as.vector(c(2,3))
split.screen(s)


screen(1)

par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('Admixture_mode_fixedPopSize_NAm1110000_EUr700000_rho_5e-9_theta_0_selection_False_MS_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i],dgrp, ylab='DGRP', xlab='Simulations',main='',cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)


mtext(expression(A),at=lab_pos,padj=0,cex=1.5)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('Admixture_mode_fixedPopSize_NAm15984500_EUr700000_rho_5e-9_theta_0_selection_False_MS_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i],dgrp, ylab='', xlab='Simulations',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(B),at=lab_pos,padj=0,cex=1.5)



close.screen(all=TRUE)

dev.off() 




##########
# Pi/bp  #
##########

setwd('/Users/nanditagarud/Documents/Jensen_response/analysis/S_and_Pi')
chr2L=read.table('Chr2L_SNPdensity_perShortIntron.txt')
chr2R=read.table('Chr2R_SNPdensity_perShortIntron.txt')
chr3L=read.table('Chr3L_SNPdensity_perShortIntron.txt')
chr3R=read.table('Chr3R_SNPdensity_perShortIntron.txt')

dgrp= (c(chr2L[,4],chr2R[,4], chr3L[,4],chr3R[,4])) # pi/bp


i=2 #column with pi
lab_pos=-0.04 #position of letter label

setwd('/Users/nanditagarud/Documents/Jensen_response/analysis')
pdf("Figure_supp_tmp_Pi_QQ_1.pdf",width=8,height= 6, title = "Figure 2",paper="special")

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

d=read.table('constNe10e6_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i], dgrp, ylab='DGRP', xlab='',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(A),at=lab_pos,padj=0,cex=1.5)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('constNe2.7e6_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i], dgrp, ylab='', xlab='',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(B),at=lab_pos,padj=0,cex=1.5)


screen(3)
par(mai=c(mb3,ml3,mt3,mr3))
d=read.table('dadi1_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i], dgrp,  xlab='',ylab='', main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(C),at=lab_pos,padj=0,cex=1.5)


screen(4)
par(mai=c(mb1,ml1,mt1,mr1))

d=read.table('dadi2_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i], dgrp, xlab='Simulations',ylab='DGRP', main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)


mtext(expression(D),at=lab_pos,padj=0,cex=1.5)


screen(5)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('admixture_mode_Garud2015_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i], dgrp, ylab='', xlab='Simulations',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(E),at=lab_pos,padj=0,cex=1.5)


screen(6)
par(mai=c(mb3,ml3,mt3,mr3))

d=read.table('admixture_bot_mode_Garud2015_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i], dgrp, ylab='', xlab='Simulations',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)


mtext(expression(F),at=lab_pos,padj=0,cex=1.5)

close.screen(all=TRUE)

dev.off() 

###############

pdf("Figure_supp_tmp_Pi_QQ_2.pdf",width=8,height= 6, title = "Figure 2",paper="special")

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

d=read.table('admixture_mode_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i],dgrp,ylab='DGRP', xlab='',main='',cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(A),at=lab_pos,padj=0,cex=1.5)


screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('admixture_bot_mode_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i],dgrp,ylab='', xlab='',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(B),at=lab_pos,padj=0,cex=1.5)


screen(3)
par(mai=c(mb3,ml3,mt3,mr3))
d=read.table('admixture_posterior_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i],dgrp, xlab='',ylab='', main='',cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(C),at=lab_pos,padj=0,cex=1.5)


screen(4)
par(mai=c(mb1,ml1,mt1,mr1))

d=read.table('admixture_posterior_Harris_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i],dgrp, xlab='Simulations',ylab='DGRP', main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(D),at=lab_pos,padj=0,cex=1.5)


screen(5)

par(mai=c(mb2,ml2,mt2,mr2)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('Arguello_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i],dgrp, ylab='', xlab='Simulations',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)

mtext(expression(E),at=lab_pos,padj=0,cex=1.5)


close.screen(all=TRUE)

dev.off() 


#####
###############
pdf("Figure_supp_tmp_Pi_QQ_3.pdf",width=8,height= 6, title = "Figure 2",paper="special")

s=as.vector(c(2,3))
split.screen(s)


screen(1)

par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.

d=read.table('Admixture_mode_fixedPopSize_NAm1110000_EUr700000_rho_5e-9_theta_0_selection_False_MS_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i],dgrp, ylab='DGRP', xlab='Simulations',main='',cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)


mtext(expression(A),at=lab_pos,padj=0,cex=1.5)


legend(0.03, 0.8, legend=c("DGRP data", "Simulations"),
       col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch=15, cex=0.8)



screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('Admixture_mode_fixedPopSize_NAm15984500_EUr700000_rho_5e-9_theta_0_selection_False_MS_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

w=qqplot(d[,i],dgrp, ylab='', xlab='Simulations',main='', cex.lab=0.8,axes=F)
abline(a=0,b=1, col='red')
axis(1,  cex.axis=0.8)
axis(2,  cex.axis=0.8)

distance=rmse(w$x,w$y)
distance=format(round(distance, 4))
mtext(distance,at=0.01,padj=0,cex=1)


mtext(expression(B),at=lab_pos,padj=0,cex=1.5)



close.screen(all=TRUE)

dev.off() 






########################################################
# Plot H12 and H2/H1 for different selection scenarios #
########################################################

setwd('/Users/nanditagarud/Documents/Jensen_response/analysis/')

mean_h12_hard=c()
mean_h2h1_hard=c()
mean_h12_soft=c()
mean_h2h1_soft=c()

for (s in c('0','0.0001','0.001','0.005','0.01','0.05','0.1')) {
	file=paste('selection_simulations_constNe_rho_5e-9_theta_0.01_selection_', s, '_age_0_partialfreq_0.5_MS_snps.txt', sep = "")
	d=read.table(file)
	h12=mean(d[,9])
	h2h1=mean(d[,10])
	mean_h12_hard = c(mean_h12_hard, h12)
	mean_h2h1_hard = c(mean_h2h1_hard, h2h1)
	
	file=paste('selection_simulations_constNe_rho_5e-9_theta_10_selection_', s, '_age_0_partialfreq_0.5_MS_snps.txt', sep = "")	
	d=read.table(file)
	h12=mean(d[,9])
	h2h1=mean(d[,10])
	mean_h12_soft = c(mean_h12_soft, h12)
	mean_h2h1_soft = c(mean_h2h1_soft, h2h1)
	
}




pdf("Figure_H12_H2H1_sims.pdf",width=8,height= 4, title = "Figure 2",paper="special")

sel=c(0.00001,0.0001, 0.001, 0.005, 0.01, 0.05, 0.1)

s=as.vector(c(1,2))
split.screen(s)

change=0.15

mb1=0.8 
ml1=0.8
mt1=0.3 
mr1=0.4 

mb2=0.8 
ml2=0.3 
mt2=0.3 
mr2=0.8

screen(1)

par(mai=c(mb1,ml1,mt1,mr1)) #A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.


plot(sel,mean_h12_hard, type ="o",pch=19, ylab = "", xlab = "Selection Strength", col = "blue", ylim=c(0,1), yaxt='n',xaxt='n', log='x')
axis(2, col.axis = 'blue', cex.axis = 1,las = 2)
axis(1, at=c(0.00001,0.0001, 0.001, 0.01, 0.1), labels=c('0',expression(paste('10'^'-4')),expression(paste('10'^'-3')),expression(paste('10'^'-2')),expression(paste('10'^'-1'))))
mtext("H12", side = 2, line = 3, col='blue')
par(new = TRUE)
plot(sel, mean_h2h1_hard, type = "o", pch=19, xaxt = "n", yaxt = "n", ylab = "", xlab = "", col = "red", ylim=c(0,1),log='x')

#mtext(expression(A),at=-0.1,padj=0,cex=1.5)

screen(2)
par(mai=c(mb2,ml2,mt2,mr2))

d=read.table('admixture_bot_mode_neutrality_shortIntron_Pi_S_TajD_perBp_4Fields.txt')

plot(sel,mean_h12_soft, type ="o",pch=19, ylab = "", xlab = "Selection Strength", col = "blue", ylim=c(0,1), yaxt='n',xaxt='n',log='x')
axis(1, at=c(0.00001,0.0001, 0.001, 0.01, 0.1), labels=c('0',expression(paste('10'^'-4')),expression(paste('10'^'-3')),expression(paste('10'^'-2')),expression(paste('10'^'-1'))))
par(new = TRUE)
plot(sel, mean_h2h1_soft, type = "o",pch=19, xaxt = "n", yaxt = "n", ylab = "", xlab = "", col = "red", ylim=c(0,1),log='x')
#axis(side = 4)
axis(4, col.axis = 'red', cex.axis = 1,las = 2)
mtext("H2/H1", side = 4, line = 3, col='red')
legend(0.00001, 0.4, c("H12", "H2/H1"),
       col = c("blue", "red"), lty = c(1, 1))


#mtext(expression(B),at=-0.1,padj=0,cex=1.5)


close.screen(all=TRUE)

dev.off() 



#####################################################################
# Figure
# plot the haplotype spectra for simulations under selection
#####################################################################

print_sel=c(expression(paste('0'^'')),expression(paste('10'^'-4')),expression(paste('10'^'-3')),expression(paste('10'^'-2')),expression(paste('10'^'-1'))))
setwd('/Users/nanditagarud/Documents/Jensen_response/analysis/')

for (theta in c('0.01','10')){

pdf(paste("Figure_theta", theta,".pdf", sep=''),width=3.5,height=3, title = "HapDist",paper="special")
s=as.vector(c(1,5))
split.screen(s)
par(omi=c(0.1,0.30,0.15,0.05)) #A vector of the form c(bottom, left, top, right) giving the size of the outer margins in inches.

cex2=0.6

counter=1
for (s in c('0','0.0001','0.001','0.01','0.1')) {
	file=paste('selection_simulations_constNe_rho_5e-9_theta_',theta,'_selection_', s, '_age_0_partialfreq_0.5_MS_snps.txt', sep = "")
	simulation=read.table(file)


screen(counter)
# need to adjust the margins
par(mai=c(.1,0,.05,0.03))

a1 <- as.character(simulation[10,5]) # Randomly sample the 10th simulation
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

colNo=1
for (j in sort(hapDist)){
	height=c(height[2]+0.3,height[2]+j)
	#rect(wl,height[1],wr,height[2],density=-1,col=cols[floor(height[1])%%15+1],lwd=0, border=cols[floor(height[1])%%15+1])}
	rect(wl,height[1],wr,height[2],density=-1,col=cols[length(hapDist)-colNo+1],lwd=0, border='white')
	colNo=colNo+1
	}

# Print H12 and H2/H1 values
H12=paste(round(simulation[10,9],3))
H2H1=paste(round(simulation[10,10],3))
mtext(print_sel[counter], side=1, cex=cex2,padj=-1)
#mtext(H12,side=1,cex=cex2,padj=-0.5)
#mtext(H2H1,side=1,cex=cex2,padj=1)

counter = counter +1
}
#mtext(expression(A),at=-1300,padj=0,cex=0.9)
#mtext(expression(B),at=-1300,padj=0,cex=0.9)

#mtext("H12:",side=1,at=-1305,padj=-2, cex=cex2)
#mtext("H2/H1:",side=1,at=-1295,padj=-0.5, cex=cex2)

close.screen(all=TRUE)
dev.off() 
}


###############################################
# Plot QQ plots for everything and assess fit #
###############################################


