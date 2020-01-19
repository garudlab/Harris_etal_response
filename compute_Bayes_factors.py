
#! /bin/env python
import sys
from optparse import OptionParser
import copy
import numpy
from numpy import *
import os.path


######################

def Var_Dist(inFile1, inFile2, outFile, N):

    # iterate through the HS file

    # first randomly draw N lines to study:
    entries_HS=[]
    for line in inFile1:
        entries_HS.append(line)
    entries_HS=numpy.asarray(entries_HS)
    entries_HS=numpy.random.choice(entries_HS,N,replace=False)

    entries_SS=[]
    for line in inFile2:
        entries_SS.append(line)
    entries_SS=numpy.asarray(entries_SS)
    entries_SS=numpy.random.choice(entries_SS,N,replace=False)

    # store N simulations in a vector
    
    H12_HS=[]
    H2H1_HS=[]
    
    for line in entries_HS:
        lineList=line.split('\t')
        num_fields=len(lineList)
        if len(lineList)==17:
            H12_in = float(lineList[8])
            H2H1_in = float(lineList[9]) 
            if H12_in <= 1 and H12_in >= 0 and H2H1_in <=1 and H2H1_in >=0:
                H12_HS.append(H12_in)
                H2H1_HS.append(H2H1_in)
                
    H12_HS=numpy.asarray(H12_HS)
    H2H1_HS=numpy.asarray(H2H1_HS)
    
    # iterate through the SS file
    count=0
    H12_SS=[]
    H2H1_SS=[]

    for line in entries_SS:
        lineList=line.split('\t')
        num_fields=len(lineList)
        if len(lineList)==17:
            H12_in = float(lineList[8])
            H2H1_in = float(lineList[9]) 
            if H12_in <= 1 and H12_in >= 0 and H2H1_in <=1 and H2H1_in >=0:
                H12_SS.append(H12_in)
                H2H1_SS.append(H2H1_in)
           
    H12_SS=numpy.asarray(H12_SS)
    H2H1_SS=numpy.asarray(H2H1_SS)
    
    # calculate the variances
    allH12=numpy.concatenate((H12_HS, H12_SS), axis=0)
    allH2H1=numpy.concatenate((H2H1_HS, H2H1_SS), axis=0)


    varH12=numpy.var(allH12)
    varH2H1=numpy.var(allH2H1)

    for H12 in numpy.arange(0.025, 1.025, 0.025):
        for H2H1 in numpy.arange(0.025, 1.025, 0.025):
            #print H12, H2H1
            #compute the distances from the observed H12 value
            distance_HS=((float(H12) - H12_HS)**2/varH12 + (float(H2H1) - H2H1_HS)**2/varH2H1)**0.5
            distance_SS=((float(H12) - H12_SS)**2/varH12 + (float(H2H1) - H2H1_SS)**2/varH2H1)**0.5

            # count the distances that are within 0.1
            num_HS = (distance_HS <= 0.1).sum()
            num_SS = (distance_SS <= 0.1).sum()

            #print to outfile
            outFile.write(str(H12) + '\t' + str(H2H1) + '\t' + str(num_HS) +'\t' +str(num_SS) + '\n')
    
    
###############


def mkOptionParser():
    """ Defines options and returns parser """
    
    usage = """%prog  <input1_HS> <input2_SS> <output_HS> <output_SS> <H12_obs> <H2/H1_obs> <N>
    %prog calculates the variance of H12 and H2/H1 for a set of simulations for a given peak (thetaA=0.01 for HS and thetaA=10 for SS) Also calculates the euclidean distance from the observed H12 and H2/H1 values for a given peak. Outputs N (specified by user) simulations with the distance reported. """

    parser = OptionParser(usage)
   

    return parser



def main():
    """ see usage in mkOptionParser. """
    parser = mkOptionParser()
    options, args= parser.parse_args()

    if len(args) != 4:
        parser.error("Incorrect number of arguments")


    inFN1         = os.path.expanduser(args[0])
    inFN2         = os.path.expanduser(args[1])
    outFN        = os.path.expanduser(args[2]) 
    N             = int(args[3])
    

    if inFN1 == '-':
        inFile1 = sys.stdin
    else:
        inFile1      = open(inFN1, 'r')

    if inFN2 == '-':
        inFile2 = sys.stdin
    else:
        inFile2      = open(inFN2, 'r')

    if outFN == '-':
        outFile = sys.stdout
    else:
        outFile      = open(outFN, 'w')
        
    Var_Dist(inFile1, inFile2, outFile, N)


    

#run main
if __name__ == '__main__':
    main()
