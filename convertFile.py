
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
        s=''
        for i in range(1, len(line.split('\t'))):
            s+= line.split('\t')[i]+ ','
        outFile.write(s+ '\n')



    
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
