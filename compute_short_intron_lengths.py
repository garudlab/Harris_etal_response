outFile=open('/u/home/n/ngarud/Jensen_response/data/DGRP_shortIntron/short_intron_lengths.txt','w')

shortIntron_lengths=[]
for chromosome in ['2R','2L','3L','3R']:
    inFN='/u/home/n/ngarud/Jensen_response/analysis/S_and_Pi/Chr' + chromosome  + '_SNPdensity_perShortIntron.txt'
    inFile=open(inFN, 'r')
    for line in inFile:
        items=line.strip().split('\t')
        start=int(items[0])
        end=int(items[1])
        length=end-start+1
        outFile.write(str(length) + '\n')
        
