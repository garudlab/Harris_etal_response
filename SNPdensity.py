
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

def clusterHaplotypes(inFile, outFile):
    # Every 1KB, count the number of SNPs and the number of SNPs with Ns in them

    lineNumber = 1
    #SNP_count = 0
    flies = initialize()
    coord_prev=0

    for line in inFile:
        genotypes=line.strip().split(',')
        coord = int(genotypes[0])
        print(str(coord) + '\t' + str(lineNumber*10000))
        

        if coord < lineNumber*10000:
            #add to the flies vector

            for i in range(1,146):
                flies[i].append(line.split(',')[i])            
            coord_prev=coord
        
        else:  # If we have surpassed the window length            
            S,Pi=calc_S_Pi(flies)
            S=S/10000.0
            Pi=Pi/10000.0
            flies = initialize()
            outFile.write(str(coord_prev) + '\t' + str(S) + '\t' + str(Pi)  + '\n' )
            print(str(coord_prev) + '\t' + str(S) + '\t' + str(Pi))
            lineNumber += 1
            #SNP_count =1

    # Print the last line

    outFile.write(str(coord_prev) + '\t' + str(S)  + '\t' + str(Pi)  + '\n' )
    


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
#    print frequencies
#    print Pi
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
