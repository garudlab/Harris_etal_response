
#! /bin/env python
# bedTileElements - tile each record in input.bed with potentially shorter elements which are written to output.bed 
import sys
from optparse import OptionParser
import copy
import csv
#import numpy
import linecache
import random

######################

def convertMS(outFile, N):
    
    outFile.write('0.5,')
    for i in range(0,N-1):
        outFile.write('A,')
    outFile.write('A')
            
###############


def mkOptionParser():
    """ Defines options and returns parser """
    
    usage = """%prog  <input file> <output file> 
    %prog converts MS output to genotypes (0/1 to A/G). """

    parser = OptionParser(usage)
   

    return parser



def main():
    """ see usage in mkOptionParser. """
    parser = mkOptionParser()
    options, args= parser.parse_args()

    if len(args) != 2:
        parser.error("Incorrect number of arguments")


    N            = int(args[0])
    outFN        = args[1]

    if outFN == '-':
        outFile = sys.stdout
    else:
        outFile      = open(outFN, 'w')

    convertMS(outFile,N)


    
#run main
if __name__ == '__main__':
    main()
