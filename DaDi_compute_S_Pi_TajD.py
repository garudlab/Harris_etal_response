#! /bin/env python
import sys
from optparse import OptionParser
import copy
import os.path

# Numpy is the numerical library dadi is built upon
import scipy.optimize
import numpy
from numpy import array

import dadi



# Parse the data file to generate the data dictionary
#dd = dadi.Misc.make_data_dict('dataForDaDi_David.txt') # these include all 162 lines and is not the exact data used for any paper.    

dd = dadi.Misc.make_data_dict(os.path.expanduser('~/Jensen_response/data/short_introns/AllData_forDaDi_DGRP_IBD.txt'))   # these match the statistics reported in the PloS 2015 paper. All IBD lines removed. 

#dd = dadi.Misc.make_data_dict('ShortIntrons_all_RA_DaDi.txt')   # these were matched to the data included in Zambian paper. Very similar statistics to the above. 


# projected down to 130 samples.
fs = dadi.Spectrum.from_data_dict(dd, ['DME'], [130], polarized=False)

locusLength=738024

print(fs.Tajima_D())
print(fs.Watterson_theta()/locusLength)
print(fs.pi()/locusLength)
print(fs.S()/locusLength)

