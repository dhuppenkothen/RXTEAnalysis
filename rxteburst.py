
#from __future__ import with_statement
#import numpy as np
import fnmatch
import os
import argparse
#import scipy.special
#import scipy.stats
#import scipy.optimize
#import scipy
import glob

try:
    import xbblocks
except ImportError:
    print("You need to have the xbblocks script from Peter Williams' pwpy repository in "
          "order to find bursts!")


import matplotlib
matplotlib.use("Agg")

from pylab import *
#rc("font", size=20, family="serif", serif="Computer Sans")
#rc("text", usetex=True)

#import lightcurve

try:
    import burst
    print("Burst sucessfully imported!")
except ImportError:
    print("Module Burst not found!")

import rxte
import findbursts


def search_filenames_recursively(testdir, testexpression):
    matches = []
    for root, dirnames, filenames in os.walk(testdir):
        for filename in fnmatch.filter(filenames, testexpression):
            matches.append(os.path.join(root, filename))

    return matches





def read_burst_times(filename):
    """
    Read burst start and end times from file.
    First column: start times in MET
    second column> end times in MET
    """
    data = rxte.conversion(filename)
    tstart = np.array([float(t) for t in data[0]])
    tend = np.array([float(t) for t in data[1]])
    blen = tend - tstart

    return tstart, blen



def make_bursts(datafile, bursttimefile, bary=True, fileroot="test"):

    len_datafile = len(datafile.split("/")[-1])

    data = rxte.RXTEData(times=None, channels=None, datafile=datafile, npcus=None, ra=None, dec=None,
                 emid = None, emiddir=datafile[:-len_datafile], bary=bary)

    tstart, blen = read_burst_times(bursttimefile)

    for i,(s,l) in enumerate(zip(tstart, blen)):
        #print("First photon: " + str(data.photons[0].unbary))
        #print("Last photon: " + str(data.photons[-1].unbary))
        #print("start time: " + str(s-data.t0))
        if data.photons[0].unbary <= s-data.t0 <= data.photons[-1].unbary:
            try:
                b = rxte.RXTEBurst(s, l, data.photons, data.t0, bary=bary, add_frac=0.2, fnyquist=4096.0, norm="leahy",
                                   pcus = data.pcus)
                b.save_burst(fileroot + "_" + str(s) + "_burst.dat")
            except rxte.ZeroCountsException:
                print("No counts in burst!")
                continue

        else:
            continue
    return


def all_bursts(datadir = "./", data_expression="*.asc", bursttimefile="bursts.dat"):

    filenames = search_filenames_recursively(datadir, data_expression)
    for f in filenames:
        fsplit = f.split("/")
        froot = fsplit[1]
        flen = len(fsplit[-1])
        fdir = f[:-flen]
        print("froot: " + str(froot))
        make_bursts(f, bursttimefile, fileroot=fdir+froot)

    return

def plot_lightcurves(datadir="./"):

    filenames = glob.glob(datadir + "*burst.dat")
    for f in filenames:
        b = rxte.getpickle(f)
        bintimes, bincounts = rxte.rebin_lightcurve(b.lc.time, b.lc.countrate, 10)
        plot(bintimes, bincounts/np.max(b.pcus), lw=2, color='black', linestyle='steps-mid')
        xlabel("Time since trigger [s]")
        ylabel("Count rate [cts/s]")
        fsplit = f.split("_")
        title(fsplit[0] + ", tstart = " + fsplit[1] + ", pcus = " + str(np.max(b.pcus)))
        fnew = fsplit[0] + "_" + fsplit[1] + "_lc.png"
        savefig(fnew, format="png")
        close()

    return


def bayesian_analysis(nwalker=500, niter=200, nsim=1000, fitmethod='powell', datadir="./", froot="test"):

    filenames = glob.glob(datadir + froot + "*burst.dat")

    for f in filenames:
        print("I am on burst %s" %f)
        b = rxte.getpickle(f)
        fsplit = f.split("/")
        namestr = fsplit[-1][:-10]
        ps = b.ps
        if len(ps.ps) > 4096:
            binps = ps.rebinps(1.0)
            m = int(binps.df/ps.df)
            b.ps = binps
        else:
            m = 1
        b.bayesian_analysis(namestr=namestr, nchain=nwalker, niter=niter, nsim=nsim, m=m, fitmethod=fitmethod)

        b.save_burst(namestr+"_burstfile.dat")
    return

def main():

    if extract_bursts:
        print("Running all_bursts ...")
        assert clargs.bfile, "No file with burst start times!"
        all_bursts(bursttimefile=clargs.bfile)
    if plot_lcs:
        plot_lightcurves()

    if analysis:
        bayesian_analysis(nwalker=clargs.nwalker, niter=clargs.niter, nsim=clargs.nsim, datadir=clargs.datadir,
                          froot=clargs.froot)

    print("Done!")
    return



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Make burst files out of RXTE data!')
    parser.add_argument('-b', '--bfile', action='store', dest='bfile', required=False,
                         help='File with burst start times')
    parser.add_argument("-d", "--data-dir", action="store", dest="datadir", required=False, default="./",
                        help="Directory where data is located")



    parser.add_argument("-e", "--extract-bursts", action="store_true", dest='extract',
                        help="Would you like to extract bursts from data?")
    parser.add_argument("-l", "--plot-lightcurves", action="store_true", dest='plot_lc',
                        help="Would you like to plot light curves to file?")


    parser.add_argument("-a", "--analysis", action="store_true", dest="analysis",
                        help="Would you like to run the Bayesian PSD analysis on a bunch of files?")

    parser.add_argument("-w", "--nwalker", action="store", dest="nwalker", required=False, default=500, type=int,
                        help="Number of walkers for MCMC run")
    parser.add_argument("-i", "--niter", action="store", dest="niter", required=False, default=200, type=int,
                        help="Number of iterations for MCMC run")
    parser.add_argument("-s", "--nsim", action="store", dest="nsim", required=False, default=1000, type=int,
                        help="Number of periodograms to simulate from posterior distribution")
    parser.add_argument("--fr", "--froot", action="store", dest="froot", default="", required=False,
                        help="File root that can be specified to run the analysis only on a subset of bursts")

    clargs = parser.parse_args()
    extract_bursts = clargs.extract
    plot_lcs = clargs.plot_lc
    analysis = clargs.analysis
    main()
