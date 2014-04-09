from __future__ import with_statement
from collections import defaultdict

import glob
from operator import itemgetter
import numpy as np
import cPickle as pickle
import pyfits

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

    def __init__(self, time, channel, pcu, unbarycentered=None, emid=None):

        self.time = time

        ### if no energy conversion is given, assume that what was
        ### used as the channel is actually the energy
        if emid is None:
            self.channel = None
            self.energy = channel

        else:
            self.channel = int(channel)
            self.convertchannel(emid)

        self.pcu = pcu

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

        if emid is None and not emiddir is None:
            self.emid = self.readchannelconv(emiddir)
        elif not emid is None:
            self.emid = emid
        else:
            self.emid = np.zeros(300)

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
            #print("barytimes[0]: %f" %barytimes[0])
            #print("channels[0]: %f" %channels[0])
            #print("times[0]: %f" %times[0])
            #print("pcus[0]: %f" %pcus[0])
            self.photons = [RXTEPhoton(t,e,p,u,emid) for t,e,u,p in zip(barytimes, channels, times, pcus)]
        else:
            self.photons = [RXTEPhoton(t,e,p,emid=emid) for t,e,p in zip(times, channels, pcus)]

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

        print("times[0]: " + str(times[0]))
        print("times[-1]: " + str(times[-1]))

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

    def deadtime_correction(self, bst_unbary, bend_unbary, std1dir="./", vle_correction="max"):

        """
                Dead time corrections for RXTE data. Reads out Standard 1 files,
                then finds the right time interval and calculates the VLE count rates and source
                count rates. Returns arrays with the correction to the periodogram and
                the corrected periodogram.

                bst_unbary: start time of the burst in total *unbarycentered* seconds MET
                bend_unbary: end time of the burst in total *unbarycentered* seconds MET
                std1dir:    directory with standard 1 data, should be the directory that has
                            subdirectories ace, acs, cal, ..., pca etc.
                vle_correction: for VLE correction, compute either "mean" or "max" of the VLE count rate
        """


        ### find standard1 data filename
        if std1dir[-1] == "/":
            std1_subdir = std1dir.split("/")[-2]
            std1_all_dir = std1dir[:-(len(std1_subdir)+1)]
        else:
            std1_subdir = std1dir.split("/")[-1]
            std1_all_dir = std1dir[:-len(std1_subdir)]
        print("Subdirectory with Standard 1 Data is %s" %std1_subdir)
        print("Subdirectory with Standard 1 summary file is %s" %std1_all_dir)

        ## open list with standard 1 filenames
        std1f = open(std1_all_dir + "std1.std", "r")
        std1_files = std1f.readlines()
        std1_split = [s.split("/") for s in std1_files]

        ## std1_subdir is the directory that corresponds to the last directory in std1dir,
        ## note that Lucy's code gives these in *absolute* paths, which are on her computer, not mine
        std1_subdir_all = [s[-3] for s in std1_split]
        std1_files_all = [s[-1][:-1] for s in std1_split]

        ## loop through list to find right subdirectory and filename
        std1file = []
        for s,f in zip(std1_subdir_all, std1_files_all):
            if s == std1_subdir:
                std1file.append(f)
            else:
                continue

        if len(std1file) == 0:
            raise NoDeadTimeFileException()


        all_times, all_xe_total, all_vlecnt, all_vpcnt, all_remainingcnt = [], [], [], [], []



        for f in std1file:


            if std1dir[-1] == "/":
                std1path = std1dir + "pca/" + f
            else:
                std1path = std1dir + "/pca/" + f

            std1 = glob.glob(std1path+"*")[0]
            hdulist = pyfits.open(std1)
            data = hdulist[1].data
            time = data["Time"]
            tstart = time[0]
            tend = time[-1]+0.125*1024.0
            times = np.arange(tstart, tend, 0.125)

            all_times.extend(times)

            ### VLE events
            vlecnt = data["VLECnt"].flatten()
            all_vlecnt.extend(vlecnt)

            ### good xenon events
            xecntpcu0 = data["XeCntPcu0"].flatten()
            xecntpcu1 = data["XeCntPcu1"].flatten()
            xecntpcu2 = data["XeCntPcu2"].flatten()
            xecntpcu3 = data["XeCntPcu3"].flatten()
            xecntpcu4 = data["XeCntPcu4"].flatten()

            ### propane layer events
            vpcnt = data["VpCnt"]

            ### coincident events
            remainingcnt = data["RemainingCnt"]

            xe_total = xecntpcu0+xecntpcu1+xecntpcu2+xecntpcu3+xecntpcu4

            all_xe_total.extend(xe_total)


        ## check that times are actually sorted
        dt = np.array(all_times[1:]) - np.array(all_times[:-1])
        dt_subzero = np.where(dt < 0.0)[0]

        ### if there are time differences between bins <0, then the list isn't sorted
        ### this is bad and we need to fix it!
        if len(dt_subzero) > 0:
            allzip = zip(all_times, all_vlecnt, all_xe_total)
            allzip_sorted = sorted(allzip, key=itemgetter(0))

            all_times = np.array(allzip_sorted)[:,0]
            all_vlecnt = np.array(allzip_sorted)[:,1]
            all_xe_total = np.array(allzip_sorted)[:,2]


        ### figure out number of PCUs:
        min_pcus = np.min(self.pcus)
        max_pcus = np.max(self.pcus)

        if not min_pcus == max_pcus:
            print("Number of PCUs changes during burst. I might as well give up now. Setting to max(pcus) instead.")
            npcus = max_pcus
        else:
            npcus = max_pcus

        ts = all_times.searchsorted(bst_unbary)
        te = all_times.searchsorted(bend_unbary)+1

        burst_times = all_times[ts:te]
        if vle_correction == "max":
            burst_vle = np.max(all_vlecnt[ts:te])/float(npcus)
            burst_xe_total = np.max(all_xe_total[ts:te])/float(npcus)
        elif vle_correction == "mean":
            burst_vle = np.mean(all_vlecnt[ts:te])/float(npcus)
            burst_xe_total = np.mean(all_xe_total[ts:te])/float(npcus)


        freq = self.ps.freq
        nfreq = len(freq)

        ### dead time for VLE correction
        tau_vle = 1.7e-4

        ### VLE correction
        p_vle = 2.0*burst_vle*burst_xe_total*(tau**2.0)*(np.sin(np.pi*tau*freq)/(np.pi*tau*freq))**2.0


        ### dead time for chain correction
        tau_d = 1.0e-5
        binsize = self.lc.res
        fnyquist = np.max(freq)

        p1 = 2.0*(1.0 - )



        return


class NoDeadTimeFileException(Exception):
    def __init__(self):
        print("No Standard 1 data file found!")
        return


class ZeroCountsException(Exception):

    def __init__(self, ncounts):
        self.ncounts = ncounts
        return






