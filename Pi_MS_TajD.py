import sys
from optparse import OptionParser
import copy


######################

def clusterHaplotypes(inFile, outFile):
    #Every 10Kb, calculate Pi
    nSam=sample_size +1 
    
    lineNumber = 1
    SNP_count = 0
    flies = initialize()
    
    for line in inFile:
        coord = float(line.split(',')[0])
        
        SNP_count +=1
        #add to the flies vector
        for i in range(1,nSam):
            flies[i].append(line.split(',')[i])
        

    # calculate pi
    Pi=calcPi_2(flies)
    avgPairwiseDiff=calcPi(flies)
    TajD = calcTajD(SNP_count, avgPairwiseDiff, flies)
    outFile.write(str(SNP_count) + '\t' + str(Pi) + '\t' + str(avgPairwiseDiff) +'\t' + str(TajD) + '\n' )
    
def initialize():
    nSam = sample_size +1
    flies={}
    for i in range(1,nSam):
        flies[i] = []
    return flies

def calcPi_2(flies):
    # in this definition I will calculate allele frequencies and return a vector of the frequencies in the given population

    nSam=sample_size + 1
    
    frequencies = []

    for j in range(0, len(flies[1])):
        nucleotides = [0,0,0,0]
        for i in range(1,nSam):

            if flies[i][j] == 'A':
                nucleotides[0] +=1
            if flies[i][j] == 'T':
                nucleotides[1] +=1
            if flies[i][j] == 'G':
                nucleotides[2] +=1
            if flies[i][j] == 'C':
                nucleotides[3] +=1
        
        # check whether or not this SNP has exactly two alleles           
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

#    print frequencies
#    print Pi
    return Pi

def calcPi(flies):
    # caluclate the average pairwise difference in SNPs
    # make a nested for loop and calculate the pairwise difference between each pair of strains
    nSam = sample_size +1
    pairwiseDiffArray = []
    numComparisons = 0
    
    for i in range(1, nSam-1):
        for j in range(i+1,nSam):
            diff = pairwiseDiff(flies, i, j)
            pairwiseDiffArray.append(diff)
            numComparisons +=1

    
    avgPairwiseDiff = float(sum(pairwiseDiffArray))/float(numComparisons)
    return avgPairwiseDiff

def pairwiseDiff(flies, i, j):
    diff=0
    
    for w in range(0, len(flies[i])):
        if flies[i][w] != flies[j][w] and flies[i][w] != 'N' and flies[j][w] != 'N':
            diff+=1

    return diff


def calcTajD(numSNPs, avgPairwiseDiff, flies):
    # Need to calculate the normalizing constant

    # a1
    a1=0
    a2=0
    n=sample_size
    for i in range(1, n):
        a1 += 1/float(i)
        a2 += 1/(float(i))**2

    b1 = (n+1)/float((3*(n-1)))
    b2 = 2*(n**2 + n + 3)/float((9*n*(n-1)))
    c1 = b1 - 1/float(a1)
    c2 = b2 - float(n+2)/float(a1*n) + float(a2)/float(a1)**2
    e1 =c1/float(a1)
    e2 = c2/float(a1**2 + a2)

    if numSNPs > 0:
        TajD = (avgPairwiseDiff - numSNPs/float(a1))/(e1*numSNPs + e2*numSNPs*(numSNPs-1))**0.5
    else:
        TajD =0
        
    return TajD



    
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
    global sample_size
    sample_size         = int(args[2])

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
