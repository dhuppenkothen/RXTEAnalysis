
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
import findbursts

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
    len_processed = len(datafile.split("/")[-2])


    data = rxte.RXTEData(times=None, channels=None, datafile=datafile, npcus=None, ra=None, dec=None,
                        emid = None, emiddir=None, bary=bary)
                 #emiddir=datafile[:-len_datafile], bary=bary)
    print("data unbarycentered t_start = %f" %data.photons[0].unbary)
    print("data barycentered t_start = %f" %data.photons[0].time)
    #print("bary: " + str(bary))
    #print("all photons, unbarycentered: " + str([p.unbary for p in data.photons]))

    tstart, blen = read_burst_times(bursttimefile)

    for i,(s,l) in enumerate(zip(tstart, blen)):
        #print("First photon: " + str(data.photons[0].unbary))
        #print("Last photon: " + str(data.photons[-1].unbary))
        print("start time: " + str(s-data.t0))
        if data.photons[0].unbary <= s-data.t0 <= data.photons[-1].unbary:
            try:
                b = rxte.RXTEBurst(s, l, data.photons, data.t0, bary=bary, add_frac=0.2, fnyquist=2048.0, norm="leahy",
                                   pcus = data.pcus)
                b.ps_corr = b.deadtime_correction(std1dir=datafile[:-(len_datafile+len_processed+1)])
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

def plot_lightcurves(datadir="./", plot_ps = True):

    filenames = glob.glob(datadir + "*burst.dat")
    for f in filenames:
        b = rxte.getpickle(f)
        if plot_ps:
            fig = figure(figsize=(24, 9))
            plt.subplots_adjust(top=0.9, bottom=0.1, left=0.05, right=0.95, wspace=0.2, hspace=0.2)
            ax = fig.add_subplot(121)
        else:
            fig = figure(figsize=(12,9))
            plt.subplots_adjust(top=0.9, bottom=0.1, left=0.05, right=0.95, wspace=0.2, hspace=0.2)
            subplot(111)
        #bintimes, bincounts = rxte.rebin_lightcurve(b.lc.time, b.lc.countrate, 10)
        #plot(bintimes, bincounts/np.max(b.pcus), lw=2, color='black', linestyle='steps-mid')
        plot(b.lc.time, b.lc.countrate/np.max(b.pcus), lw=2, color='black', linestyle='steps-mid')

        xlabel("Time since trigger [s]")
        ylabel("Count rate [cts/s]")
        if plot_ps:
            ax = fig.add_subplot(122)

            loglog(b.ps.freq[1:], b.ps.ps[1:], lw=2, color='black', linestyle='steps-mid',
                        label="Periodogram")
            if not b.ps_corr is None:
                loglog(b.ps_corr.freq[1:], b.ps_corr.ps[1:], lw=2, color='cyan', linestyle='steps-mid',
                        label="Dead-time corrected periodogram")

            nfreq = np.array(b.ps.freq).searchsorted(250.0)
            psmean = np.mean(b.ps.ps[nfreq:])
            psvar = np.var(b.ps.ps[nfreq:])
            npowers = len(b.ps.ps[nfreq:])

            if not b.ps_corr is None:
                nfreq_corr = np.array(b.ps_corr.freq).searchsorted(250.0)
                psmean_corr = np.mean(b.ps_corr.ps[nfreq:])
                psvar_corr = np.var(b.ps_corr.ps[nfreq:])
                npowers_corr = len(b.ps_corr.ps[nfreq:])


            theovar = 2.0

            axis([np.min(b.ps.freq[1:]), np.max(b.ps.freq[1:]), np.min(b.ps.ps[1:])/2.0, np.max(b.ps.ps[1:])*2.0])

            if not b.ps_corr is None:
                ax.text(0.05, 0.1, r"uncorrected: $\mu = %.3f$, $\sigma^2 = %.3f$" %(psmean, psvar) + " (%.2f)" %theovar +
                                           "\n" + "for %i powers \n" %npowers +
                                    r"corrected: $\mu = %.3f$, $\sigma^2 = %.3f$" %(psmean_corr, psvar_corr) + " (%.2f)" %theovar +
                                           "\n" + "for %i powers \n" %npowers_corr,
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes,
                                color='green', fontsize=20)

            else:
                ax.text(0.05, 0.1, r"$\mu = %.3f$, $\sigma^2 = %.3f$" %(psmean, psvar) + " (%.2f)" %theovar +
                                           "\n" + "for %i powers" %npowers,
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes,
                                color='green', fontsize=20)

            hlines(2.0, b.ps.freq[1], b.ps.freq[-1], lw=2, color='red', linestyle='dashed',
                   label="Theoretical Poisson level")
            legend()

            xlabel("Frequency [Hz]", fontsize=22)
            ylabel("Leahy Power", fontsize=22)
        fsplit = f.split("_")
        title(fsplit[0] + ", tstart = " + fsplit[1] + ", pcus = " + str(np.max(b.pcus)))
        fnew = fsplit[0] + "_" + fsplit[1] + "_lc.png"
        savefig(fnew, format="png")
        close()

    return


def bayesian_analysis(nwalker=500, niter=200, nsim=1000, fnyquist=2048.0, fitmethod='powell', datadir="./", froot="test"):

    logfile = findbursts.TwoPrint(datadir + "sgr1900_bayesiananalysis.dat")

    filenames = glob.glob(datadir + froot + "*burst.dat")

    for f in filenames:
        logfile("I am on burst %s" %f)
        b = rxte.getpickle(f)
        fsplit = f.split("/")
        namestr = fsplit[-1][:-10]
        if not b.ps_corr is None:
            logfile("Running on dead-time corrected periodogram")
            ps = b.ps_corr
        else:
            logfile("No dead-time corrected Periodogram for this burst. Running uncorrected periodogram instead.")
            ps = b.ps
        if len(ps.ps) > fnyquist:
            logfile("Frequency resolution > 1 Hz: rebinning periodogram to 1 Hz")
            binps = ps.rebinps(1.0)
            m = int(binps.df/ps.df)
            b.ps = binps
        else:
            m = 1
        logfile("Now running Bayesian Analysis")
        logfile("Output saved with filename root %s." %namestr)
        b.bayesian_analysis(namestr=namestr, nchain=nwalker, niter=niter, nsim=nsim, m=m, fitmethod=fitmethod)
        logfile("All done. Saving data and periodogram search results in %s"%(namestr+"_burstfile.dat"))
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
                          froot=clargs.froot, fitmethod=fitmethod)

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

    parser.add_argument("--fitmethod", action="store", dest="fitmethod", required=False, default="bfgs",
                        help="Method to use for fitting the data. Default is BFGS. For choices consult the "
                             "scipy manual.")

    clargs = parser.parse_args()
    extract_bursts = clargs.extract
    plot_lcs = clargs.plot_lc
    analysis = clargs.analysis
    fitmethod = clargs.fitmethod
    main()
