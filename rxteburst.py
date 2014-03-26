
#from __future__ import with_statement
import numpy as np
import fnmatch
import os
import argparse
import scipy.special
import scipy.stats
import scipy.optimize
import scipy
import glob

import matplotlib
matplotlib.use("Agg")

from pylab import *

rc("font", size=20, family="serif", serif="Computer Sans")
rc("text", usetex=True)

import lightcurve

try:
    import burst
    print("Burst sucessfully imported!")
except ImportError:
    print("Module Burst not found!")
import rxte


def straight(x, a,b):
    res = a*x+b
    return res

class PoissonPosterior(object):

    def __init__(self, times, counts, model):
        self.times = np.array(times)
        #print("len(times): " + str(len(times)))
        self.counts = np.array(counts)
        #print("len(counts): " + str(len(counts)))
        self.model = model
        self.delta = self.times[1] - self.times[0]
        self.countrate = self.counts/self.delta
        print('mean count rate: ' + str(np.mean(self.countrate)))
        #print("len(countrate): " + str(len(self.countrate)))

    def logprior(self, theta):
        return 1.0

    def _log_likelihood(self, lambdas, data):

        #print("len(lambdas): " + str(len(lambdas)))
        #print("len(data): " + str(len(data)))

        llike = -np.sum(lambdas) + np.sum(data*np.log(lambdas))\
                -np.sum(scipy.special.gammaln(data + 1))

        return llike

    def loglikelihood(self, theta):

        #print("len(times): " + str(len(self.times)))

        lambdas = self.model(self.times, *theta)
        #print("len(lambdas): " + str(len(lambdas)))
        return self._log_likelihood(lambdas, self.countrate)


    def logposterior(self, theta):
        return self.logprior(theta) + self.loglikelihood(theta)


    def __call__(self, theta, neg=False):

        if neg:
            return -self.logposterior(theta)
        else:
            return self.logposterior(theta)



def fit_bkg(times, counts, model, theta_init):

    lpost = PoissonPosterior(times, counts, model)


    scipy_version = scipy.__version__
    if scipy_version >= "0.11.0":
        res = scipy.optimize.minimize(lpost, theta_init, method="Nelder-Mead", args=(True,))

    if res.success:
        print("Minimization successful!")
    else:
        print("Minimization NOT successful!")
        print(res.message)

    popt = res.x
    exit_status= res.status
    success = res.success

    return popt, exit_status, success


def interpolate_bkg(s1_times, s1_counts, searchbin_time, searchbin_counts, s2_times=None, s2_counts=None,
                    plotlc = None, model=straight):

    s1_times = np.array(s1_times)
    s1_counts = np.array(s1_counts)

    if not s2_times is None and not s2_counts is None:
        s2_times = np.array(s2_times)
        s2_counts = np.array(s2_counts)

    if model == straight:
        slope_init1 = 0.0
        norm_init1 = np.mean(s1_counts)/(s1_times[1] - s1_times[0])
        theta_init1 = [slope_init1, norm_init1]

        if not s2_times is None and not s2_counts is None:
            slope_init2 = 0.0
            norm_init2 = np.mean(s2_counts)/(s2_times[1]-s2_times[0])
            theta_init2 = [slope_init2, norm_init2]

    else:
        raise Exception("Model not recognized!")


    #print("norm_init: " + str(norm_init1))
    #print("s1_times: " + str(s1_times))
    #print("s1_counts: " + str(s1_counts))

    popt1, exit_status1, success1 = fit_bkg(s1_times, s1_counts, model, theta_init1)
    #print("popt1: " + str(popt1))

    searchbin_model1 = model(searchbin_time, *popt1)


    if not plotlc is None:
        figure()
        plot(s1_times, s1_counts/(s1_times[1]-s1_times[0]), lw=2, color='black', linestyle='steps-mid')
        s1_model = model(s1_times, *popt1)
        #print("len(s1_model): " + str(len(s1_model)))
        plot(s1_times, s1_model, lw=2, color='red', linestyle='steps-mid')
        scatter(searchbin_time, searchbin_model1)
        savefig(plotlc + ".png", format="png")
        close()

    if not s2_times is None and not s2_counts is None:
        popt2, exit_status2, success2 = fit_bkg(s2_times, s2_counts, model, theta_init2)
        searchbin_model2 = model(searchbin_time, *popt2)

    else:
        searchbin_model2 = None
        popt2 = None
        exit_status2 = None
        success2 = None

    #print("searchbin_model1: " + str(searchbin_model1))
    #print("searchbin_model2: " + str(searchbin_model2))
    return searchbin_model1, popt1, exit_status1, success1, searchbin_model2, popt2, exit_status2, success2


def compare_to_poisson(searchbin_model, searchbin_countrate):

    print("searchbin_model: " + str(searchbin_model))
    print("searchbin_countrate: " + str(searchbin_countrate))


    rv = scipy.stats.poisson(searchbin_model)

    cdf = rv.cdf(searchbin_countrate)

    print("cdf: " + str(cdf))

    pval = 1.0 - cdf

    return pval


def search_data(times, dt=0.1, dt_small=0.01, tseg=5.0, bin_distance=3.0, model=straight, plotlc="test"):

    ### make sure time array starts at 0
    times = np.array(times) - times[0]

    ### make a light curve with pre-determined time step
    lc = lightcurve.Lightcurve(times, timestep=dt)
    lcsmall = lightcurve.Lightcurve(times, timestep=dt_small)

    pvals_all, popt_all1, popt_all2, es_all1, es_all2, success_all1, success_all2 = [], [], [], [], [], [], []

    for i,(t,c) in enumerate(zip(lc.time, lc.counts)):
        searchbin_time = t
        searchbin_counts = c
        searchbin_countrate = c/lc.res

        #print("searchbin_time: " + str(searchbin_time))
        #print("searchbin_counts: " + str(searchbin_counts))

        if searchbin_time < bin_distance+tseg:
            print("Only segment after bin to search taken into account")
            s1_start = lcsmall.time.searchsorted(searchbin_time+bin_distance)
            s1_end = lcsmall.time.searchsorted(searchbin_time+bin_distance+tseg)

            s1_times = lcsmall.time[s1_start:s1_end]
            s1_counts = lcsmall.counts[s1_start:s1_end]

            s1_times = [t1 for t1,c1 in zip(s1_times, s1_counts) if c1 < 1000.0]
            s1_counts = [c1 for c1 in s1_counts if c1 < 1000.0]

            #print("start time: " + str(s1_times[0]))
            #print("end time: " + str(s1_times[-1]))


            if len(s1_counts) < 10 :
                pvals_all.append(None)
                continue


            s2_times = None
            s2_counts = None

        elif bin_distance+tseg <= searchbin_time <= times[-1]-(bin_distance+tseg):
            print("Both segments taken into account")
            s1_start = lcsmall.time.searchsorted(searchbin_time-(bin_distance+tseg))
            s1_end = lcsmall.time.searchsorted(searchbin_time-(bin_distance))

            s1_times = lcsmall.time[s1_start:s1_end]
            s1_counts = lcsmall.counts[s1_start:s1_end]


            s1_times = [t1 for t1,c1 in zip(s1_times, s1_counts) if c1 < 1000.0]
            s1_counts = [c1 for c1 in s1_counts if c1 < 1000.0]

            if len(s1_counts) < 10:
                s1_times = None
                s1_counts = None

            s2_start = lcsmall.time.searchsorted(searchbin_time+bin_distance)
            s2_end = lc.time.searchsorted(searchbin_time+(bin_distance+tseg))

            s2_times = lcsmall.time[s2_start:s2_end]
            s2_counts = lcsmall.counts[s2_start:s2_end]

            s2_times = [t2 for t2,c2 in zip(s2_times, s2_counts) if c2 < 1000.0]
            s2_counts = [c2 for c2 in s2_counts if c2 < 1000.0]

            if len(s2_counts) < 10:
                s2_times = None
                s2_counts = None

            if s1_counts is None and s2_counts is None:
                pvals_all.append(None)
                continue


        else:
            print("Only prior segment taken into account")
            s1_start = lcsmall.time.searchsorted(searchbin_time-(bin_distance+tseg))
            s1_end = lcsmall.time.searchsorted(searchbin_time-(bin_distance))

            s1_times = lcsmall.time[s1_start:s1_end]
            s1_counts = lcsmall.counts[s1_start:s1_end]

            s1_times = [t1 for t1,c1 in zip(s1_times, s1_counts) if c1 < 1000.0]
            s1_counts = [c1 for c1 in s1_counts if c1 < 1000.0]

            if len(s1_counts) < 10:
                pvals_all.append(None)
                continue

            s2_times =None
            s2_counts = None

        if not plotlc is None:
            plotlc_new = plotlc + str(i)

        else:
            plotlc_new = None

        searchbin_model1, popt1, exit_status1, success1, searchbin_model2, popt2, exit_status2, success2,= \
            interpolate_bkg(s1_times, s1_counts, searchbin_time, searchbin_counts, s2_times=s2_times,
                            s2_counts=s2_counts, plotlc=plotlc_new, model=model)

        if not searchbin_model2 is None:
            print("I AM HERE")
            searchbin_model = (searchbin_model1 + searchbin_model2)/2.0
            popt_all2.append(popt2)
            es_all2.append(exit_status2)
            success_all2.append(success2)
        else:
            print("I AM HERE OTHERWISE")
            searchbin_model = searchbin_model1


        pval = compare_to_poisson(searchbin_model*lc.res, searchbin_counts)

        print("The p-value for this bin is: p = " + str(pval))

        pvals_all.append(pval)
        popt_all1.append(popt1)
        es_all1.append(exit_status1)
        success_all1.append(success1)

        burstdict = {"lc":lc, "pvals":pvals_all, "popt1":popt_all1, "popt2":popt_all2, "exit1":es_all1,
                     "exit2":es_all2, "success1":success_all1, "success2":success_all2}

    return burstdict


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


def bayesian_analysis(nwalker=500, niter=200, nsim=1000, datadir="./", froot="test"):

    filenames = glob.glob(datadir + froot + "*burst.dat")

    for f in filenames:
        b = rxte.getpickle(f)
        fsplit = f.split("/")
        namestr = fsplit[-1][:-10]
        b.bayesian_analysis(namestr=namestr, nchain=nwalker, niter=niter, nsim=nsim)

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
    parser.add_argument("-f", "--froot", action="store", dest="froot", default="", required=False,
                        help="File root that can be specified to run the analysis only on a subset of bursts")

    clargs = parser.parse_args()
    extract_bursts = clargs.extract
    plot_lcs = clargs.plot_lc
    analysis = clargs.analysis
    main()