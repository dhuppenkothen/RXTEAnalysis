import matplotlib
matplotlib.use("Agg")

from pylab import *

#rc("font", size=20, family="serif", serif="Computer Sans")
#rc("text", usetex=True)


import numpy as np
import scipy.special
import scipy.stats
import scipy.optimize
import scipy
import cPickle as pickle

try:
    import xbblocks
except ImportError:
    print("You need to have the xbblocks script from Peter Williams' pwpy repository in "
          "order to find bursts!")


import lightcurve


import subprocess

def git_hash(gitdir = "./"):
    return subprocess.check_output(["git", "--git-dir", gitdir+"/.git", "rev-parse", "HEAD"])

counts_cutoff = 600.0

def straight(x, a,b):
    res = a*x+b
    return res



class TwoPrint(object):

    def __init__(self,filename):
        self.file = open(filename, "w")
        self.filename = filename
        self.file.write("##\n")
        self.close()
        return

    def __call__(self, printstr):
        print(printstr)
        self.file = open(self.filename, "a")
        self.file.write(printstr + "\n")
        self.close()
        return

    def close(self):
        self.file.close()
        return



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

        #print("lambdas: " + str(lambdas))
        #print("data: " + str(data))

        llike = -np.sum(lambdas) + np.sum(data*np.log(lambdas))\
                -np.sum(scipy.special.gammaln(data + 1))

        return llike

    def loglikelihood(self, theta):

        #print("theta: " + str(theta))
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

    #if res.success:
        #print("Minimization successful!")
    if not res.success:
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

        if searchbin_time < times[0]+bin_distance+tseg:
            #print("Only segment after bin to search taken into account")
            s1_start = lcsmall.time.searchsorted(searchbin_time+bin_distance)
            s1_end = lcsmall.time.searchsorted(searchbin_time+bin_distance+tseg)

            s1_times = lcsmall.time[s1_start:s1_end]
            s1_counts = lcsmall.counts[s1_start:s1_end]

            s1_times = [t1 for t1,c1 in zip(s1_times, s1_counts) if 0.0 < c1 < counts_cutoff*lcsmall.res]
            s1_counts = [c1 for c1 in s1_counts if 0.0 < c1 < counts_cutoff*lcsmall.res]

            #print("len(s1_times): " + str(len(s1_times)))
            #print("len(s1_counts): " + str(len(s1_counts)))



            #print("start time: " + str(s1_times[0]))
            #print("end time: " + str(s1_times[-1]))

            if len(s1_counts) < 50:
                pvals_all.append(2.0)
                continue


            s2_times = None
            s2_counts = None

        elif times[0]+bin_distance+tseg <= searchbin_time <= times[-1]-(bin_distance+tseg):
            #print("Both segments taken into account")
            s1_start = lcsmall.time.searchsorted(searchbin_time-(bin_distance+tseg))
            s1_end = lcsmall.time.searchsorted(searchbin_time-(bin_distance))

            s1_times = lcsmall.time[s1_start:s1_end]
            s1_counts = lcsmall.counts[s1_start:s1_end]


            s1_times = [t1 for t1,c1 in zip(s1_times, s1_counts) if 0.0 < c1 < counts_cutoff*lcsmall.res]
            s1_counts = [c1 for c1 in s1_counts if 0.0 < c1 < counts_cutoff*lcsmall.res]

            #print("len(s1_times): " + str(len(s1_times)))
            #print("len(s1_counts): " + str(len(s1_counts)))


            s2_start = lcsmall.time.searchsorted(searchbin_time+bin_distance)
            s2_end = lcsmall.time.searchsorted(searchbin_time+(bin_distance+tseg))

            s2_times = lcsmall.time[s2_start:s2_end]
            s2_counts = lcsmall.counts[s2_start:s2_end]

            s2_times = [t2 for t2,c2 in zip(s2_times, s2_counts) if 0.0 < c2 < counts_cutoff*lcsmall.res]
            s2_counts = [c2 for c2 in s2_counts if 0.0 < c2 < counts_cutoff*lcsmall.res]

            #print("len(s1_times): " + str(len(s2_times)))
            #print("len(s1_counts): " + str(len(s2_counts)))



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

            s1_times = [t1 for t1,c1 in zip(s1_times, s1_counts) if 0.0 < c1 < counts_cutoff*lcsmall.res]
            s1_counts = [c1 for c1 in s1_counts if 0.0 < c1 < counts_cutoff*lcsmall.res]

            #print("len(s1_times): " + str(len(s1_times)))
            #print("len(s1_counts): " + str(len(s1_counts)))


            if len(s1_counts) < 50:
                pvals_all.append(2.0)
                continue

            s2_times =None
            s2_counts = None

        #plotlc = "test"

        if not plotlc is None:
            plotlc_new = plotlc + str(i)

        else:
            plotlc_new = None

        try:
            searchbin_model1, popt1, exit_status1, success1, searchbin_model2, popt2, exit_status2, success2,= \
                interpolate_bkg(s1_times, s1_counts, searchbin_time, searchbin_counts, s2_times=s2_times,
                                s2_counts=s2_counts, plotlc=plotlc_new, model=model)

        except FloatingPointError:
            print("Fitting failed! Leaving out bin.")
            pvals_all.append(2.0)
            continue

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

        if isnan(pval):
            pval = 2.0

        #print("The p-value for this bin is: p = " + str(pval))
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



def burst_parameters(times, pvals_dict, t0=None, sig_threshold=1.0e-8, p0=0.05, nbootstrap=512, plotlc=None,
                     parfile=None):

    if parfile is None:
        parfile = TwoPrint(plotlc+"_parameters.dat")
    ## read out all p-values
    pvals_all = np.array(pvals_dict["pvals"])

    #print("pvals_all: " + str(pvals_all))

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


    parfile("Found " + str(len(pvals_sig_new)) + " putative bursts. Running Bayesian Blocks ...")

    ## extract original light curve
    lc = pvals_dict["lc"]

    sigtimes = [lc.time[i] for i in pvals_sig_new]

    tseg_before = 1.0
    tseg_after = 4.0

    parfile("Time added before bin with outlier for Bayesian Blocks search: %f s.\n" % tseg_before)
    parfile("Time added after bin with outlier for Bayesian Blocks search: %f s.\n" % tseg_after)

    burst_file = open(plotlc + "_burst_times.dat", "w")

    if not t0 is None:
        burst_file.write("#Start time \t End time (both in MET)\n")

    else:
        burst_file.write("#Start time \t End time (in time since first photon)\n")


    all_info = []

    for i,s in enumerate(sigtimes):
        #parfile("I am on burst " + str(i))
        ## extract a region 2 seconds before and after each significant time bin
        si = times.searchsorted(s-tseg_before)
        ei = times.searchsorted(s+tseg_after)

        tnew = times[si:ei]
        #parfile("tseg: " + str(tnew[-1] - tnew[0]))

        if not plotlc is None:
            plotlc_new = plotlc + "_b" + str(i)
        else:
            plotlc_new = plotlc

        info = bayesian_blocks(tnew, p0, nbootstrap, plotlc_new, parfile)
        if len(info.bsrates) <=2 :
            print("This is not a burst!")
            continue
        else:

            edges = np.array(info.burst_edges)

            for j,e in enumerate(edges):
                data_file = open(plotlc_new + "_n" + str(j) + "_data.dat", "w")
                si_new = tnew.searchsorted(e[0])
                ei_new = tnew.searchsorted(e[1])
                tnew_burst = tnew[si_new:ei_new]
                if not t0 is None:
                    tnew_burst = np.array(tnew_burst) + t0
                    data_file.write("#Time in MET seconds \n")
                else:
                    data_file.write("#Time since obs start in seconds\n]")
                for t in tnew_burst:
                    data_file.write(str(t) + "\n")

                data_file.close()

            print("edges: " + str(edges))
            if not t0 is None:
                edges += t0

            for e in edges:
                burst_file.write(str(e[0]) + "\t" + str(e[1]) + "\n")

            all_info.append(info)

    burst_file.close()

    return all_info


def bayesian_blocks(times, p0=0.05, nbootstrap=512, plotlc=None, parfile=None):

    ## run Bayesian blocks algorithm
    info = xbblocks.bsttbblock(times, [np.min(times)], [np.max(times)], nbootstrap=nbootstrap, p0=p0)


    ### make a light curve
    lc = lightcurve.Lightcurve(times, timestep=0.01)
    ## extract edges
    edges = list(info.ledges)
    edges.append(info.redges[-1])
    edges = np.array(edges)


    ## bins with lowest count rate is probably background,
    ## so use those as background level
    sorted_rates = np.sort(info.bsrates)
    bkgmean = np.mean(sorted_rates[:2])

    rv = scipy.stats.norm(bkgmean, np.sqrt(bkgmean))

    pvals = [1.0 - rv.cdf(b) for b in info.bsrates]

    burst_ind = [1 if p < 0.01 else 0 for p in pvals]
    burst_diff = [burst_ind[i+1] - b  for i,b in enumerate(burst_ind[:-1])]

    burst_start, burst_end = [], []

    ## time in seconds to add on either side of the burst to make sure I catch
    ## all of the burst data
    add_time = 0.02
    if not parfile is None:
        parfile("Time added to either side of the burst: %f" % add_time)

    for i,b in enumerate(burst_diff):
        if b == 1:
            burst_start.append(info.ledges[i+1]-add_time)
        elif b == -1:
            if not len(burst_start) == 0:
                if burst_start[-1] < info.redges[i]:
                    burst_end.append(info.redges[i]+add_time)
                else:
                    continue
            else:
                continue


    #burst_edges = [edges[1]-0.02, edges[-2]+0.02]
    info.burst_edges = zip(burst_start, burst_end)

    if not plotlc is None:

        #min_time = edges[1]-0.1
        #max_time = edges[-2]+0.1
        min_time = lc.time[0]
        max_time = lc.time[-1]

        figure(figsize=(15,9))
        plot(lc.time, lc.countrate, lw=2, color='black', linestyle='steps-mid', label="Data light curve")
        fill_between(edges.repeat(2)[1:-1], info.bsrates.repeat(2), facecolor='white', lw=2, edgecolor='red',
                     label="Bayesian Blocks representation")
        plot(lc.time, np.zeros(len(lc.time)), lw=2, color='red', label="Bayesian Blocks representation")
        max_cr = np.max(lc.countrate)+0.1*np.max(lc.countrate)
        axis([min_time, max_time, 0.2, max_cr])
        for i,b in enumerate(info.burst_edges):
            vlines(b[0], 0.1, max_cr, lw=2, color='green', linestyle='dashed', label="Burst edges, burst " + str(i))
            vlines(b[1], 0.1, max_cr, lw=2, color='green', linestyle='dashed')
        legend(prop={"size":16})
        if len(info.bsrates) <= 2:
            plt.title("Not a burst")
        else:
            plt.title("It's a burst!")
        savefig(plotlc + "_bblocks.png", format="png")
        close()

    return info



def extract_bursts(times, t0=None, p0=0.05, nbootstrap=512, sig_threshold=1.0e-7, dt=0.1, dt_small=0.01,
                   tseg=5.0, bin_distance=3.0, froot="test", parfile=None):


    if parfile is None:
       parfile = TwoPrint(froot + "_parameters.dat")

    parfile("Searching data from %s to %f.\n" % (times[0], times[-1]))
    parfile("Bin size of the light curve when searching for outliers: dt = %f s\n" % dt)
    parfile("Bin size of the light curve used for fitting during outlier search: dt_small = %f s\n" % dt_small)
    parfile("Size of segment used for fitting background is: t_seg = %f \n" % tseg)
    parfile("Distance between bin searched and segments used in background fitting: bin_distance = %f \n" % bin_distance)
    parfile("Significance threshold for outlier search is: p < %.10f \n" % sig_threshold)
    parfile("Results saved in files with root %s \n" % froot)

    pval_dict, sig_all = search_data(times, dt, dt_small, tseg, bin_distance, straight, None,
                sig_threshold)


    parfile("Significance threshold for Bayesian Blocks, block acceptance: p = %f" % p0)
    parfile("Number of bootstrapping operations in Bayesian Blocks: n = %f" % nbootstrap)

    all_info = burst_parameters(times, pval_dict, t0, sig_threshold, p0, nbootstrap, froot, parfile)

    f = open(froot + "burst_info_dict.dat", "w")
    pickle.dump(all_info,f)
    f.close()

    return

