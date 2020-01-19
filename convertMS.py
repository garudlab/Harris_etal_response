
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

def convertMS(inFile, outFile):
    
    data = {}
    dictionary = {}
    dictionary['0'] = 'A'
    dictionary['1'] = 'C'

    sweep_center = 0
    for line in inFile:
        line = line.strip()
        if line[0] == 'p':
            positions_vector_tmp = line.split(' ')
            positions_vector = positions_vector_tmp[1:len(positions_vector_tmp)]
            # Locate where 0.50000 is in the positions vector (this is where the selected allele is)
            #sweep_center = positions_vector.index('0.50000')
            # Create a vector to store each position:
            for i in range(0, len(positions_vector)):
                data[i] = []

        else:
            for x in range(0, len(line)):
                data[x].append(dictionary[line[x]])
        
    # print out
    for y in range(0, len(positions_vector)):
        s =  str(float(positions_vector[y])) + ','
        s += ','.join(data[y]) + '\n'
        outFile.write(s)
            
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

    convertMS(inFile, outFile)


    
#run main
if __name__ == '__main__':
    main()
