
#! /bin/env python
# bedTileElements - tile each record in input.bed with potentially shorter elements which are written to output.bed 
import sys
from optparse import OptionParser
import copy
import csv
#import numpy
import linecache
import random
import re
####################################################################################
# path psyco
####################################################################################
#sys.path.append("/home/jsp/prog/utillities/py_modules")
# to speed things up
#import psyco
#psyco.full()
#import bed



######################

def scrambleGenotypes(inFile, outFile):

    values={}
    #iterate1=20
    #iterate2=150
    iterate1=20
    iterate2=150

#    interval = 0+iterate1
    interval=1
    for line in inFile:
      if len(line.split('\t')) ==8 and re.match("^\d+?\.\d+?$", line.split('\t')[4]) is not None:
        LD=line.split('\t')[4]
        bp=line.strip().split('\t')[3]
        values.setdefault(interval, [])
        #values[bp].append(float(LD))
        if int(bp) <= interval and bp !=0:
            values[interval].append(float(LD))
        else:
            if interval >300:
                interval += iterate2
            else:
                interval += iterate1
            
            values.setdefault(interval, [])
            values[interval].append(float(LD))

    sortedKeys = sorted(values.keys())
    for key in sortedKeys:
        if key !=0 and len(values[key]) >0:
            outFile.write(str(key) + '\t' + str(sum(values[key])/len(values[key])) + '\n')

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

    scrambleGenotypes(inFile, outFile)


    
#run main
if __name__ == '__main__':
    main()
