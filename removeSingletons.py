
#! /bin/env python
# bedTileElements - tile each record in input.bed with potentially shorter elements which are written to output.bed 
import sys
from optparse import OptionParser
import copy

####################################################################################
# path psyco
####################################################################################
#sys.path.append("/home/jsp/prog/utillities/py_modules")
# to speed things up
#import psyco
#psyco.full()
#import bed



######################

def clusterHaplotypes(inFile, outFile):

    for line in inFile:
        line=line.strip('\n')
        nucleotides = [0,0,0,0]
        lineList=line.split(',')
        for i in range(1, len(lineList)):
            if lineList[i] == 'A':
                nucleotides[0] +=1            
            if lineList[i] == 'T':
                nucleotides[1] +=1
            if lineList[i] == 'G':
                nucleotides[2] +=1
            if lineList[i] == 'C':
                nucleotides[3] +=1
        
        if max(nucleotides) != sum(nucleotides)-1:
            outFile.write(line+ '\n')




    
###############


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

    if len(args) != 2:
        parser.error("Incorrect number of arguments")


    inFN         = args[0]
    outFN        = args[1]
 

    if inFN == '-':
        inFile = sys.stdin
    else:
        inFile      = open(inFN, 'r')

    if outFN == '-':
        outFile = sys.stdout
    else:
        outFile      = open(outFN, 'w')



    clusterHaplotypes(inFile, outFile)


    

#run main
if __name__ == '__main__':
    main()
