from __future__ import with_statement
from collections import defaultdict

import glob
import numpy as np
import cPickle as pickle

import burst
import generaltools
import lightcurve
import powerspectrum

def conversion(filename):
    f=open(filename, 'r')
    output_lists=defaultdict(list)
    for line in f:
        if not line.startswith('#'):
             line=[value for value in line.split()]
             for col, data in enumerate(line):
                 output_lists[col].append(data)
    return output_lists

def getpickle(picklefile):
    file = open(picklefile, 'r')
    procdata = pickle.load(file)
    return procdata

def rebin_lightcurve(times, counts, n=10):

    nbins = int(len(times)/n)
    dt = times[1] - times[0]
    T = times[-1] - times[0] + dt
    bin_dt = dt*n
    bintimes = np.arange(nbins)*bin_dt + bin_dt/2.0 + times[0]

    nbins_new = int(len(counts)/n)
    counts_new = counts[:nbins_new*n]
    bincounts = np.reshape(np.array(counts_new), (nbins_new, n))
    bincounts = np.sum(bincounts, axis=1)
    bincounts = bincounts/np.float(n)

    #bincounts = np.array([np.sum(counts[i*n:i*n+n]) for i in range(nbins)])/np.float(n)
    #print("len(bintimes): " + str(len(bintimes)))
    #print("len(bincounts: " + str(len(bincounts)))
    if len(bintimes) < len(bincounts):
        bincounts = bincounts[:len(bintimes)]

    return bintimes, bincounts



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
                 emid = None, emiddir = "./", bary=True, pcus=None):

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

    def readchannelconv(self, emiddir="./"):

        filename = glob.glob(emiddir + "Energy*")
        print("filename: " + str(filename))
        assert len(filename) >= 1, "No energy conversion file present!"
        filename = filename[0]
        energydata = conversion(filename)
        energybins = np.array([float(e) for e in energydata[0]])
        return energybins

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




    def findbursts(self, dt=0.1):

        times = np.array([p.time for p in self.photons])
        lc = lightcurve.Lightcurve(times, timestep=dt)

        for t,c in zip(lc.times, lc.counts):
            print("nothing")

        return




def fitsegment(s1, s2, bin_time, bin_counts):
    print("I am here!")
    return





class RXTEBurst(burst.Burst, object):

    def __init__(self, bstart, blength, photons, ttrig, pcus=None, bary=True, add_frac=0.2, fnyquist=4096.0, norm="leahy"):

        print("bstart: " + str(bstart))
        print("blength: " + str(blength))

        self.bst = bstart - add_frac*blength

        self.blen = (1.0+2.0*add_frac)*blength

        self.bend = self.bst + (1.0+add_frac)*blength

        self.ttrig = ttrig

        print("self.bst: " + str(self.bst))
        print("self.blen: " + str(self.blen))
        print("self.bend: " + str(self.bend))
        print("calculated blen: " + str(self.bend - self.bst))

        if bary:
            times = np.array([p.unbary for p in photons]) + self.ttrig
        else:
            times = np.array([p.time for p in photons]) + self.ttrig

        #print("times[0]: " + str(times[0]))
        #print("times[-1]: " + str(times[-1]))

        startind = times.searchsorted(self.bst)
        endind = times.searchsorted(self.bend)

        if not pcus is None:
            if np.size(pcus) > 1:
                self.pcus = pcus[startind:endind]
            else:
                self.pcus = pcus

        print("startind: " + str(startind))
        print("endind: " + str(endind))
        print("len photons: " + str(len(photons)))

        self.photons = photons[startind:endind]
        print("len photons: " + str(len(self.photons)))
        if startind == endind:
            raise ZeroCountsException(0)

        self.times = np.array([p.time for p in self.photons])
        print("len times: " + str(len(self.times)))
        self.lc = lightcurve.Lightcurve(self.times, timestep=0.5/fnyquist, tseg=self.blen)
        self.ps = powerspectrum.PowerSpectrum(self.lc, norm=norm)

        return


class ZeroCountsException(Exception):

    def __init__(self, ncounts):
        self.ncounts = ncounts
        return






