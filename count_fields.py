
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

def count_fields(inFile, outFile):
    
    for line in inFile:
        line_list = line.split()
        if len(line_list)==17:
            outFile.write(line)
            
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

    count_fields(inFile, outFile)


    
#run main
if __name__ == '__main__':
    main()
