
from __future__ import with_statement
import numpy as np
from collections import defaultdict

import generaltools
import burst
import rxte

#### READ ASCII DATA FROM FILE #############
#
# This is a useful little function that reads
# data from file and stores it in a dictionary
# with as many lists as the file had columns.
# The dictionary has the following architecture
# (e.g. for a file with three columns):
#
# {'0':[1st column data], '1':[2nd column data], '2':[3rd column data]}
#
#
# NOTE: Each element of the lists is still a *STRING*, because
# the function doesn't make an assumption about what type of data
# you're trying to read! Numbers need to be converted before using them!
#
def conversion(filename):
    f=open(filename, 'r')
    output_lists=defaultdict(list)
    for line in f:
        if not line.startswith('#'):
             line=[value for value in line.split()]
             for col, data in enumerate(line):
                 output_lists[col].append(data)
    return output_lists



def read_burst_times(filename):
    """
    Read burst start and end times from file.
    First column: start times in MET
    second column> end times in MET
    """
    data = conversion(filename)
    tstart = np.array([float(t) for t in data[0]])
    tend = np.array([float(t) for t in data[1]])
    blen = tend - tstart

    return tstart, blen



def make_bursts(datafile, bursttimefile, bary=True, fileroot="test"):

    data = rxte.RXTEData(times=None, channels=None, datafile=datafile, npcus=None, ra=None, dec=None,
                 emid = None, emiddir = "./", bary=bary)

    tstart, blen = read_burst_times(bursttimefile)

    for i,(s,l) in enumerate(zip(tstart, blen)):
        b = rxte.RXTEBurst(s, l, data.photons, data.t0, bary=bary, add_frac=0.2, fnyquist=4096.0, norm="leahy")
        b.saveburst(fileroot + "_" + str(tstart)[:7] + ".dat")
    return