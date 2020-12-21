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



