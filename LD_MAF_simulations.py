
#! /bin/env python
# bedTileElements - tile each record in input.bed with potentially shorter elements which are written to output.bed 
import sys
from optparse import OptionParser
import copy
import csv
#import numpy
import linecache
import random
import time
import random
####################################################################################
# path psyco
####################################################################################
#sys.path.append("/home/jsp/prog/utillities/py_modules")
# to speed things up
#import psyco
#psyco.full()
#import bed



######################

def clusterHaplotypes(inFile, outFile, sampleDown, chromosome, lowerMAF, upperMAF, window, jump):
    
    # calculation of all pairs of SNPs

    countLines =  open(inFile)
    numberLines =  len(countLines.readlines())

    SNPindex_left=1
    SNPindex_right=numberLines

    if SNPindex_right != 0:
        edgeCoord1 = int(window*(float(linecache.getline(inFile,SNPindex_left).split(',')[0])))
        edgeCoord2 = int(window*(float(linecache.getline(inFile,SNPindex_right).split(',')[0])))
    else:
        edgeCoord1 = 0
        edgeCoord2 = 0

    windowIndexLeft = 0
    windowIndexRight=window
    lastSNPcoord = linecache.getline(inFile,numberLines).split(',')[0]



    # check if the left SNP has an allele frequency between 0.4 and 0.6
    SNP1_pop=linecache.getline(inFile,SNPindex_left).strip().split(',')
    MAF=calculateMAF(SNP1_pop)
    
    while MAF <lowerMAF or MAF > upperMAF:
        SNPindex_left +=1
        SNP1_pop=linecache.getline(inFile,SNPindex_left).strip().split(',')
        MAF=calculateMAF(SNP1_pop)
        edgeCoord1 = int(window*(float(linecache.getline(inFile,SNPindex_left).split(',')[0])))
        if SNPindex_left > SNPindex_right:
            SNPindex_right +=1
            edgeCoord2 = int(window*(float(linecache.getline(inFile,SNPindex_right).split(',')[0])))
            
    # calculate LD for this window
    for i in range(SNPindex_left+1, SNPindex_right):
        # get the line numbers
        SNP2_pop=linecache.getline(inFile,i).strip().split(',')
        LD_calculation(SNP1_pop, SNP2_pop, outFile)
           





#############################
        
def LD_calculation(SNP1_pop, SNP2_pop, outFile):
    
        MAF=calculateMAF(SNP1_pop[1:len(SNP1_pop)])
            # sample down
        if len(SNP1_pop)>sampleDown and len(SNP2_pop)>sampleDown: # check that the pop size is larger than sample size

            w=sorted(random.sample(xrange(1,len(SNP1_pop)), sampleDown))
      
            SNP1_down = [ SNP1_pop[z] for z in w] 
            SNP2_down = [ SNP2_pop[z] for z in w]

            # remove any Ns
            [SNP1, SNP2] = checkForNs(SNP1_down, SNP2_down)
 
            #check if the SNP is polymorphic
            if (len(set(SNP1)) > 1) and (len(set(SNP2)) > 1) and len(SNP1)>3:
                # calculate LD
                [LD_val, P_A, P_B, denom]=LD(SNP1,SNP2)
                # check if the allele frequencies are between 0.4 and 0.6:
                if P_A >=lowerMAF and P_A <=upperMAF and P_B >=lowerMAF and P_B <=upperMAF:
                    # calculate distance (bps) and which coord is bigger
                    bigCoord=''
                    smallCoord=''
                    if int(window*float(SNP1_pop[0])) > int(window*float(SNP2_pop[0])):
                        distance = int(window*float(SNP1_pop[0]))-int(window*float(SNP2_pop[0]))+1
                        bigCoord=SNP1_pop[0]
                        smallCoord=int(window*float(SNP2_pop[0]))
                    else:
                        distance = int(window*float(SNP2_pop[0]))-int(window*float(SNP1_pop[0]))
                        bigCoord=int(window*float(SNP2_pop[0]))
                        smallCoord=int(window*float(SNP1_pop[0]))
                    #write the output
                    outFile.write(chromosome + '\t' + str(int(smallCoord)-1) + '\t' + str(bigCoord) + '\t' + str(distance) + '\t' + str(format(LD_val, '.3f')) + '\t' + str(format(P_A,'.2f')) +'\t' + str(format(P_B,'.2f')) + '\t' + str(denom) + '\n')

            
##################################
def calculateMAF(SNP1):

    nucleotides1=[0,0,0,0]
    for i in range(0,len(SNP1)):
        if SNP1[i] =='A':
            nucleotides1[0] +=1
        elif SNP1[i] == 'T':
            nucleotides1[1] +=1
        elif SNP1[i] == 'G':
            nucleotides1[2] +=1
        elif SNP1[i] == 'C':
            nucleotides1[3] +=1

    # find the major allele:
    if nucleotides1.index(max(nucleotides1)) == 0:
        majorAllele1 = 'A'
    elif nucleotides1.index(max(nucleotides1)) == 1:
        majorAllele1 = 'T'
    elif nucleotides1.index(max(nucleotides1)) == 2:
        majorAllele1 = 'G'
    elif nucleotides1.index(max(nucleotides1)) == 3:
        majorAllele1 = 'C'


    #calculate frequencies
    P_A = float(max(nucleotides1))/float(sum(nucleotides1))

    return P_A




##################################
def initialize(window, SNPindex_right, SNPindex_left, inFile):

    # Add SNPs until we surpass the window
    surpassedWindow=False

    leftCoord = int(window*float(linecache.getline(inFile,SNPindex_left).split(',')[0]))
    SNPindex_right+=1
    while surpassedWindow == False:
        coord=int(window*float(linecache.getline(inFile,SNPindex_right).split(',')[0]))
        if coord <= leftCoord+window:
            SNPindex_right+=1
        else:
                surpassedWindow=True
                SNPindex_right-=1
                
    if SNPindex_right>0:
        SNPindex_left =1
    return [SNPindex_right, SNPindex_left]


##########################
def removeSNPs(jump, SNPindex_left, inFile, windowIndexLeft):
    newLeftBoundary = windowIndexLeft
    
    surpassedWindow=False
    numSNPsRemove=0

    while surpassedWindow == False:
        
        coord=int(window*float(linecache.getline(inFile,SNPindex_left).split(',')[0]))
        if coord <= newLeftBoundary:
            numSNPsRemove  += 1
            SNPindex_left +=1
        else:
            surpassedWindow=True

    return SNPindex_left


###################################
def addSNPs(jump, SNPindex_right, SNPindex_left, inFile, windowIndexRight):
    newRightBoundary = windowIndexRight
    
    surpassedWindow=False
    
    while surpassedWindow == False:
        coord=int(window*float(linecache.getline(inFile,SNPindex_right).split(',')[0]))
        
        if coord <= newRightBoundary:
            SNPindex_right +=1
        else:
            surpassedWindow=True
            SNPindex_right-=1

    if SNPindex_right >0 and SNPindex_left ==0:
        SNPindex_left =1
    return [SNPindex_right, SNPindex_left]


###################################
def sample_wr(population, k):
    # sampling with replacement
    n = len(population)
    _random, _int = random.random, int  # speed hack
    result = [None] * k
    for i in xrange(k):
        j = _int(_random() * n)
        result[i] = population[j]
    return result

################
def checkForNs(SNP1_old, SNP2_old):
    # ignore the Ns by removing them. If an N is removed from SNP1, then remove the same individual from SNP2 even if SNP2 does not have an N for that individual

    SNP1=[]
    SNP2=[]
    
    for x in range(0, len(SNP1_old)):
        if SNP1_old[x] != 'N' and SNP2_old[x] != 'N':
            SNP1.append(SNP1_old[x])
            SNP2.append(SNP2_old[x])

    return [SNP1, SNP2]

#################
def LD(SNP1, SNP2):

        
    # Frequencies of the two SNPs
    nucleotides1=[0,0,0,0]
    for i in range(0,len(SNP1)):
        if SNP1[i] =='A':
            nucleotides1[0] +=1
        elif SNP1[i] == 'T':
            nucleotides1[1] +=1
        elif SNP1[i] == 'G':
            nucleotides1[2] +=1
        elif SNP1[i] == 'C':
            nucleotides1[3] +=1

    nucleotides2=[0,0,0,0]
    for i in range(0,len(SNP2)):
        if SNP2[i] =='A':
            nucleotides2[0] +=1
        elif SNP2[i] == 'T':
            nucleotides2[1] +=1
        elif SNP2[i] == 'G':
            nucleotides2[2] +=1
        elif SNP2[i] == 'C':
            nucleotides2[3] +=1

    # find the major allele:
    if nucleotides1.index(max(nucleotides1)) == 0:
        majorAllele1 = 'A'
    elif nucleotides1.index(max(nucleotides1)) == 1:
        majorAllele1 = 'T'
    elif nucleotides1.index(max(nucleotides1)) == 2:
        majorAllele1 = 'G'
    elif nucleotides1.index(max(nucleotides1)) == 3:
        majorAllele1 = 'C'

    if nucleotides2.index(max(nucleotides2)) == 0:
        majorAllele2 = 'A'
    elif nucleotides2.index(max(nucleotides2)) == 1:
        majorAllele2 = 'T'
    elif nucleotides2.index(max(nucleotides2)) == 2:
        majorAllele2 = 'G'
    elif nucleotides2.index(max(nucleotides2)) == 3:
            majorAllele2 = 'C'

    # calculate P_AB

    AB=0
    for i in range(0,len(SNP1)):
        if SNP1[i] == majorAllele1 and SNP2[i] == majorAllele2:
            AB +=1
                
    #calculate frequencies
    P_AB = float(AB)/float(len(SNP1))
    P_A = float(max(nucleotides1))/float(len(SNP1))
    P_B = float(max(nucleotides2))/float(len(SNP2))

    # calculate LD
    LD_val = ((P_AB - P_A*P_B)**2)/(P_A*(1-P_A)*P_B*(1-P_B))

    return [LD_val, P_A, P_B, len(SNP1)]


######################
def mkOptionParser():
    """ Defines options and returns parser """
    
    usage = """%prog  <input.bed> <output.bed> <threshold>
    %prog filters out the lines that don't meet a certain threshold. """

    parser = OptionParser(usage)
   

    return parser



def main():
    """ see usage in mkOptionParser. """
    parser = mkOptionParser()
    options, args= parser.parse_args()

    if len(args) != 8:
        parser.error("Incorrect number of arguments")


#    inFN         = args[0]
    inFile       = args[0]
    outFN        = args[1]
    global sampleDown
    sampleDown = int(args[2])
    global chromosome
    chromosome = args[3]
    global lowerMAF
    lowerMAF = float(args[4])
    global upperMAF
    upperMAF =float( args[5])
    global window
    window=int(args[6])
    global jump
    jump = int(args[7])
#    if inFN == '-':
#        inFile = sys.stdin
#    else:
#        inFile      = open(inFN, 'r')

    if outFN == '-':
        outFile = sys.stdout
    else:
        outFile      = open(outFN, 'w')


    
    clusterHaplotypes(inFile, outFile, sampleDown, chromosome, lowerMAF, upperMAF, window, jump)

    
#run main
if __name__ == '__main__':
    main()
