
#! /bin/env python
# bedTileElements - tile each record in input.bed with potentially shorter elements which are written to output.bed 
import sys
from optparse import OptionParser
import copy
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

def clusterHaplotypes(inFile, outFile, chr):

    # read in the short intron coordinates
    inFile_SI='/u/home/n/ngarud/Jensen_response/data/DGRP_shortIntron/short_intron_sites_130_coords.txt'
    short_intron_coords= open(inFile_SI, 'r') 

    start_coords=[]
    end_coords=[]
    start=int(short_intron_coords.readline().strip().split('\t')[2])
    prev_coord=start
    start_coords.append(start)
    for line in short_intron_coords:
        items=line.strip().split('\t')
        chromosome=items[0]
        coord=int(items[2])
        if chromosome == chr: 
            if coord == prev_coord +1:
                prev_coord=coord
            else:
                end_coords.append(prev_coord)
                start_coords.append(coord)
                start=coord
                prev_coord=coord
        
    end_coords.append(coord)
        


    # read in the data

    flies = initialize()
    coord_idx=0

    for line in inFile:
        genotypes=line.strip().split(',')
        coord = int(genotypes[0])
        
        if coord < end_coords[coord_idx]:
            for i in range(1,146):
                flies[i].append(line.split(',')[i])                    

        else:
            S,Pi=calc_S_Pi(flies)
            length=end_coords[coord_idx] - start_coords[coord_idx] + 1
            #print length
            if length >10.0:
                S=S/float(length)
                Pi=Pi/float(length)
                if S >= 1.0:
                    print (str(start_coords[coord_idx]) + '\t' + str(end_coords[coord_idx]) + '\t' + str(S) + '\t' + str(Pi)  + '\n' )

                outFile.write(str(start_coords[coord_idx]) + '\t' + str(end_coords[coord_idx]) + '\t' + str(S) + '\t' + str(Pi)  + '\n' )
            coord_idx +=1
            flies=initialize()
            for i in range(1,146):
                flies[i].append(line.split(',')[i])                    



def initialize():
    sample_size=145
    nSam = sample_size +1
    flies={}
    for i in range(1,nSam):
        flies[i] = []
    return flies


def calc_S_Pi(flies):
    
    # In this definition I will calculate allele frequencies and return a vector of the frequencies in the given population

    sample_size=145
    nSam=sample_size + 1
    
    frequencies = []

    for j in range(0, len(flies[1])):
        nucleotides = [0,0,0,0]
        nucs2=[]
        for i in range(1,nSam):
            '''
            if flies[i][j] == 'A':
                nucleotides[0] +=1
            if flies[i][j] == 'T':
                nucleotides[1] +=1
            if flies[i][j] == 'G':
                nucleotides[2] +=1
            if flies[i][j] == 'C':
                nucleotides[3] +=1
            '''
            if flies[i][j] == 'A' or flies[i][j] == 'T' or flies[i][j] == 'G' or flies[i][j] == 'C':
              nucs2.append(flies[i][j])

        if len(nucs2) >=130:
            nucs3=random.sample(nucs2,130)
        
            # check whether or not this SNP has exactly two alleles           
            for i in range(0, len(nucs3)):
                if nucs3[i] == 'A':
                    nucleotides[0] +=1
                if nucs3[i] == 'T':
                    nucleotides[1] +=1
                if nucs3[i] == 'G':
                    nucleotides[2] +=1
                if nucs3[i] == 'C':
                    nucleotides[3] +=1


            counter=0
            for y in range(0, len(nucleotides)):
                if nucleotides[y]>0:
                    counter +=1
            if  counter == 2:
                frequencies.append(float(max(nucleotides))/sum(nucleotides))

        

    # Now iterate through the frequencies vector and calculate pi
    Pi=0
    for w in range(0, len(frequencies)):
        Pi += 2*frequencies[w]*(1-frequencies[w])
    Pi = Pi*(sample_size)/(float(sample_size-1))

    num_snps=len(frequencies)
    
    return num_snps,Pi



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

    if len(args) != 3:
        parser.error("Incorrect number of arguments")


    inFN         = args[0]
    outFN        = args[1]
    chr          = args[2]

    if inFN == '-':
        inFile = sys.stdin
    else:
        inFile      = open(inFN, 'r')

    if outFN == '-':
        outFile = sys.stdout
    else:
        outFile      = open(outFN, 'w')



    clusterHaplotypes(inFile, outFile, chr)


    

#run main
if __name__ == '__main__':
    main()
