
import numpy as np
from __future__ import with_statement
from collections import defaultdict

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



def read_rxte_data(filename):

    f = open(filename, 'r')
    times, barytimes, pcus, channels = [], [], [], []

    counter = 0
    for i,line in enumerate(f):
        line_split = line.split()
        if line.startswith("#") and counter == 0:
            try:
                t0 = np.float(line_split[-1])
            except ValueError:
                raise Exception("No t0 found! Aborting ...")

            counter += 1

        else:
            times.append(np.float(line_split[0]))
            barytimes.append(np.float(line_split[1]))
            pcus.append(np.float(line_split[4]))
            channels.append(np.float(line_split[5]))

    times = np.array(times)
    barytimes = np.array(barytimes)
    pcus = np.array(pcus)
    channels = np.array(channels)

    min_pcu = np.min(pcus)
    max_pcu = np.max(pcus)