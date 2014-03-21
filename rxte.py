
import glob
import numpy as np

import generaltools

class RXTEData(generaltools.Data, object):

    def __init__(self, times=None, channels=None, datafile=None, npcus=None, ra=None, dec=None, bary=True):

        if not ra is None and not dec is None:
            self.ra = ra
            self.dec = dec

        self.bary = bary

        assert datafile is not None or (times is not None and channels is not None), \
            "Either datafile or time and channels must be given!"

        if not datafile is None:
            self.readrxtedata(datafile)
        else:
            self.times = times
            self.channels = channels
            self.pcus = npcus


        return

    def readrxtedata(self, filename):

        """
        This function reads data from file.
        Assumes that data is in ascii format, with the following columns:
        TTE, barycentered TTE, 1/dt_native, bkg, N_pcus, channel value  (BARYCENTERED DATA)
        or TTE, 1/dt_native, bkg, N_pcus, channel value (NOT BARYCENTERED DATA)
        """

        f = open(filename, 'r')
        times, barytimes, pcus, channels = [], [], [], []

        if self.bary:
            baryind = 1
            pcusind= 4
            channelsind = 5

        else:
            baryind = 0
            pcusind = 3
            channelsind = 4


        counter = 0
        for i,line in enumerate(f):
            line_split = line.split()
            if line.startswith("#") and counter == 0:
                try:
                    self.t0 = np.float(line_split[-1])
                except ValueError:
                    raise Exception("No t0 found! Aborting ...")

                counter += 1

            else:
                if not line.startswith("#"):
                    times.append(np.float(line_split[0]))
                    if self.bary:
                        barytimes.append(np.float(line_split[baryind]))
                    pcus.append(np.float(line_split[pcusind]))
                    channels.append(np.float(line_split[channelsind]))
                else:
                    continue


        self.times = np.array(times)
        if self.bary:
            self.barytimes = np.array(barytimes)
        self.pcus = np.array(pcus)
        self.channels = np.array(channels)

        min_pcu = np.min(pcus)
        max_pcu = np.max(pcus)

        if not min_pcu == max_pcu:
            print("Number of pcus changes during observation!")
            self.pcus = pcus
        else:
            self.pcus = pcus[0]

        return