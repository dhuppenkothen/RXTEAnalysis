
import glob
import numpy as np

import generaltools
import lightcurve
import powerspectrum

class RXTEPhoton(generaltools.Photon, object):

    def __init__(self, time, channel, unbarycentered=None, emid=None):

        self.time = time

        ### if no energy conversion is given, assume that what was
        ### used as the channel is actually the energy
        if emid is None:
            self.channel = None
            self.energy = channel

        else:
            self.channel = int(channel)
            self.convertchannel(emid)

        if not unbarycentered is None:
            self.unbary = unbarycentered
        return

    def convertchannel(self, emid):
        self.energy = emid[self.channel]





class RXTEData(generaltools.Data, object):

    def __init__(self, times=None, channels=None, datafile=None, npcus=None, ra=None, dec=None,
                 emid = None, emiddir = None, bary=True):

        if not ra is None and not dec is None:
            self.ra = ra
            self.dec = dec

        self.bary = bary

        if emid is None:
            self.emid = self.readchannelconv(emiddir)
        else:
            self.emid = emid

        assert datafile is not None or (times is not None and channels is not None), \
            "Either datafile or time and channels must be given!"

        if not datafile is None:
            self.readrxtedata(datafile, self.emid)
        else:
            self.photons = [RXTEPhoton(t,e,emid=self.emid) for t,e in zip(times, channels)]
            self.pcus = npcus


        return

    def readrxtedata(self, filename, emid):

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


        if self.bary:
            self.photons = [RXTEPhoton(t,e,u,emid) for t,e,u in zip(barytimes, channels, times)]
        else:
            self.photons = [RXTEPhoton(t,e,emid=emid) for t,e in zip(times, channels)]

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


    def filterenergy(self, ebins):
        generaltools.Data.filterenergy(ebins[0], ebins[1])




class RXTEBurst(burst.Burst, object):

    def __init__(self, bstart, blength, photons, ttrig, bary=True, add_frac=0.2, fnyquist=4096.0, norm="leahy"):

        self.bst = bstart - add_frac*blength

        self.blen = (1.0+2*add_frac)*blength

        self.bend = self.bst + (1.0+add_frac)*blength

        self.ttrig = ttrig

        if bary:
            times = np.array([p.unbary for p in photons]) + self.ttrig
        else:
            times = np.array([p.time for p in photons]) + self.ttrig

        startind = times.searchsorted(self.bst)
        endind = times.searchsorted(self.bend)

        self.photons = photons[startind:endind]

        self.times = np.array([p.time for p in self.photons])
        self.lc = lightcurve.Lightcurve(self.times, timestep=0.5/fnyquist, tseg=self.blen)
        self.ps = powerspectrum.PowerSpectrum(self.lc, norm=norm)

        return







