
import argparse
import findbursts
import rxte

import numpy as np


def find_bursts(datafile, bary=True, sig_threshold=1.0e-7, nbootstrap=200, froot="test"):

    len_datafile = len(datafile.split("/")[-1])
    filename = datafile.split("/")[-1]
    datadir = datafile[:-len_datafile]

    data = rxte.RXTEData(datafile=datafile, emiddir=datadir, bary=bary)
    times = np.array([p.time for p in data.photons])

    ## prepare data: sometimes, there are photons that have arrival times screwed up, so I exclude them
    dt = times[1:] - times[:-1]
    dt_neg = np.where(dt < 0.0)[0]

    if len(dt_neg) >=1:
        print("I am here!")
        times = list(times)
        for i,t in enumerate(dt_neg):
            print("bla")
            p = times.pop(t-i)

    times = np.array(times)

    findbursts.extract_bursts(times, t0=data.t0, tseg=5.0, bin_distance=1.5, nbootstrap=nbootstrap,
                              sig_threshold=sig_threshold, froot=froot)

    return


def main():

    find_bursts(datafile, True, sig_threshold, nbootstrap, froot)

    return


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Find bursts in RXTE data! Using Bayesian Blocks! Awesomeness!')


    parser.add_argument("-f", "--filename", action="store", dest="datafile",
                        help="Data file name")

    parser.add_argument("-s", "--significance", action="store", dest="sig_threshold", required=False,
                        default=1.0e-7, type=float, help="Significance threshold for burst search")

    parser.add_argument("-n", "--nbootstrap", action="store", dest="nbootstrap", required=False, default=512,
                        type=int, help="Number of bootstrap operations in Bayesian Blocks algorithm")

    parser.add_argument("-r", "--froot", action="store", dest="froot", default="test", required=False,
                        help="Root for output plots and data files")


    clargs = parser.parse_args()

    datafile = clargs.datafile
    sig_threshold = clargs.sig_threshold
    nbootstrap = clargs.nbootstrap
    froot = clargs.froot

    main()
