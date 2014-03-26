
import numpy as np
import argparse
import scipy.special
import scipy.stats
import scipy.optimize
import scipy
import glob

try:
    import xbblocks
except ImportError:
    print("You need to have the xbblocks script from Peter Williams' pwpy repository in "
          "order to find bursts!")


import matplotlib
matplotlib.use("Agg")

from pylab import *

rc("font", size=20, family="serif", serif="Computer Sans")
rc("text", usetex=True)

import lightcurve



counts_cutoff = 600.0

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
        #print('mean count rate: ' + str(np.mean(self.countrate)))
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

    #print("searchbin_model: " + str(searchbin_model))
    #print("searchbin_countrate: " + str(searchbin_countrate))


    rv = scipy.stats.poisson(searchbin_model)

    cdf = rv.cdf(searchbin_countrate)

    #print("cdf: " + str(cdf))

    pval = 1.0 - cdf

    return pval


def search_data(times, dt=0.1, dt_small=0.01, tseg=5.0, bin_distance=3.0, model=straight, plotlc="test",
                sig_threshold=1.0e-8):

    ### make sure time array starts at 0
    times = np.array(times)

    ### make a light curve with pre-determined time step
    lc = lightcurve.Lightcurve(times, timestep=dt)
    lcsmall = lightcurve.Lightcurve(times, timestep=dt_small)

    pvals_all, popt_all1, popt_all2, es_all1, es_all2, success_all1, success_all2 = [], [], [], [], [], [], []

    significant_all = []

    for i,(t,c) in enumerate(zip(lc.time, lc.counts)):
        searchbin_time = t
        searchbin_counts = c
        searchbin_countrate = c/lc.res

        ## if there are no counts in bin, then don't search for a burst
        ## but put a place-holder in p-value list such that I can directly
        ## translate from the positions of the bin to the
        if c == 0.0:
            #print("Leaving out bin")
            pvals_all.append(2.0)
            continue

        #print("searchbin_time: " + str(searchbin_time))
        #print("searchbin_counts: " + str(searchbin_counts))

        if searchbin_time < bin_distance+tseg:
            #print("Only segment after bin to search taken into account")
            s1_start = lcsmall.time.searchsorted(searchbin_time+bin_distance)
            s1_end = lcsmall.time.searchsorted(searchbin_time+bin_distance+tseg)

            s1_times = lcsmall.time[s1_start:s1_end]
            s1_counts = lcsmall.counts[s1_start:s1_end]

            s1_times = [t1 for t1,c1 in zip(s1_times, s1_counts) if 0.0 < c1 < counts_cutoff]
            s1_counts = [c1 for c1 in s1_counts if 0.0 < c1 < counts_cutoff]

            #print("start time: " + str(s1_times[0]))
            #print("end time: " + str(s1_times[-1]))


            if len(s1_counts) < 50:
                pvals_all.append(2.0)
                continue


            s2_times = None
            s2_counts = None

        elif bin_distance+tseg <= searchbin_time <= times[-1]-(bin_distance+tseg):
            #print("Both segments taken into account")
            s1_start = lcsmall.time.searchsorted(searchbin_time-(bin_distance+tseg))
            s1_end = lcsmall.time.searchsorted(searchbin_time-(bin_distance))

            s1_times = lcsmall.time[s1_start:s1_end]
            s1_counts = lcsmall.counts[s1_start:s1_end]


            s1_times = [t1 for t1,c1 in zip(s1_times, s1_counts) if 0.0 < c1 < counts_cutoff]
            s1_counts = [c1 for c1 in s1_counts if 0.0 < c1 < counts_cutoff]


            s2_start = lcsmall.time.searchsorted(searchbin_time+bin_distance)
            s2_end = lcsmall.time.searchsorted(searchbin_time+(bin_distance+tseg))

            s2_times = lcsmall.time[s2_start:s2_end]
            s2_counts = lcsmall.counts[s2_start:s2_end]

            s2_times = [t2 for t2,c2 in zip(s2_times, s2_counts) if 0.0 < c2 < counts_cutoff]
            s2_counts = [c2 for c2 in s2_counts if 0.0 < c2 < counts_cutoff]

            if len(s2_counts) < 50:
                s2_times = None
                s2_counts = None

            if len(s1_counts) < 50:
                s1_times = s2_times
                s1_counts = s2_times
                s2_times = None
                s2_counts = None

            if s1_counts is None and s2_counts is None:
                pvals_all.append(2.0)
                continue


        else:
            #print("Only prior segment taken into account")
            s1_start = lcsmall.time.searchsorted(searchbin_time-(bin_distance+tseg))
            s1_end = lcsmall.time.searchsorted(searchbin_time-(bin_distance))

            s1_times = lcsmall.time[s1_start:s1_end]
            s1_counts = lcsmall.counts[s1_start:s1_end]

            s1_times = [t1 for t1,c1 in zip(s1_times, s1_counts) if 0.0 < c1 < counts_cutoff]
            s1_counts = [c1 for c1 in s1_counts if 0.0 < c1 < counts_cutoff]

            if len(s1_counts) < 50:
                pvals_all.append(2.0)
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

        #print("searchbin_model1: " + str(searchbin_model1))
        #print("searchbin_model2: " + str(searchbin_model2))


        if not searchbin_model2 is None:
            #print("I AM HERE")
            searchbin_model = (searchbin_model1 + searchbin_model2)/2.0
            popt_all2.append(popt2)
            es_all2.append(exit_status2)
            success_all2.append(success2)
        else:
            #print("I AM HERE OTHERWISE")
            searchbin_model = searchbin_model1


        pval = compare_to_poisson(searchbin_model*lc.res, searchbin_counts)

        print("The p-value for this bin is: p = " + str(pval))
        if pval < sig_threshold:
            significant_all.append({"time":searchbin_time, "counts":searchbin_counts, "countrate":searchbin_countrate,
                                    "pval":pval})

        pvals_all.append(pval)
        popt_all1.append(popt1)
        es_all1.append(exit_status1)
        success_all1.append(success1)

        burstdict = {"lc":lc, "pvals":pvals_all, "popt1":popt_all1, "popt2":popt_all2, "exit1":es_all1,
                     "exit2":es_all2, "success1":success_all1, "success2":success_all2}

    return burstdict, significant_all



def burst_parameters(times, pvals_dict, sig_threshold=1.0e-8, p0=0.05, nbootstrap=512, plotlc=None, froot="test"):

    ## read out all p-values
    pvals_all = np.array(pvals_dict["pvals"])

    ## find all p-values below the significance threshold defined by the user
    pvals_sig_ind = np.where(pvals_all <= sig_threshold)[0]

    ## find distance between indices
    ind_distance = np.array([pvals_sig_ind[i+1]-p for i,p in enumerate(pvals_sig_ind[:-1])])

    ## find indices where distance between events is 1:
    small_dist = np.where(ind_distance == 1)[0]

    if len(pvals_sig_ind) > 1:

        ### delete those indices that directly follow a previous index
        ### assume that consecutive indices are due to bursts that are longer than
        ### the time step in the light curve
        pvals_sig_new = []
        for i,p in enumerate(pvals_sig_ind[1:]):
            if not i in small_dist:
                pvals_sig_new.append(p)

    else:
        pvals_sig_new = pvals_sig_ind


    print("Found " + str(len(pvals_sig_new)) + " putative bursts. Running Bayesian Blocks ...")

    ## extract original light curve
    lc = pvals_dict["lc"]

    sigtimes = [lc.time[i] for i in pvals_sig_new]

    tseg = 2.0

    all_info = []

    for i,s in enumerate(sigtimes):
        print("I am on burst " + str(i))
        ## extract a region 2 seconds before and after each significant time bin
        si = times.searchsorted(s-tseg)
        ei = times.searchsorted(s+tseg)

        tnew = times[si:ei]

        if not plotlc is None:
            plotlc_new = plotlc + "_b" + str(i)
        else:
            plotlc_new = plotlc

        info = bayesian_blocks(tnew, p0, nbootstrap, plotlc_new)
        if len(info.bsrates) <=2 :
            print("This is not a burst!")
            continue
        else:
            all_info.append(info)

    edges = np.array([p.burst_edges for p in all_info])

    return all_info


def bayesian_blocks(times, p0=0.05, nbootstrap=512, plotlc=None):

    ## run Bayesian blocks algorithm
    info = xbblocks.bsttbblock(times, [np.min(times)], [np.max(times)], nbootstrap=nbootstrap, p0=p0)


    ### make a light curve
    lc = lightcurve.Lightcurve(times, timestep=0.01)
    ## extract edges
    edges = list(info.ledges)
    edges.append(info.redges[-1])
    edges = np.array(edges)

    burst_edges = [edges[1]-0.02, edges[-2]+0.02]
    info.burst_edges = burst_edges

    if not plotlc is None:

        min_time = edges[1]-0.1
        max_time = edges[-2]+0.1

        figure(figsize=(15,9))
        plot(lc.time, lc.countrate, lw=2, color='black', linestyle='steps-mid', label="Data light curve")
        fill_between(edges.repeat(2)[1:-1], info.bsrates.repeat(2), facecolor='white', lw=2, edgecolor='red',
                     label="Bayesian Blocks representation")
        plot(lc.time, np.zeros(len(lc.time)), lw=2, color='red', label="Bayesian Blocks representation")
        max_cr = np.max(lc.countrate)+0.1*np.max(lc.countrate)
        axis([min_time, max_time, 0.2, max_cr])
        vlines(edges[1]-0.02, 0.1, max_cr, lw=2, color='green', linestyle='dashed', label="Burst edges")
        vlines(edges[-2]+0.02, 0.1, max_cr, lw=2, color='green', linestyle='dashed')
        legend()
        savefig(plotlc + "_bblocks.png", format="png")
        close()

    return info



