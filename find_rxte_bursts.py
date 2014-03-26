
import argparse
import findbursts


def find_bursts(datafile, bary=True, sig_threshold=1.0e-7, nbootstrap=200, froot="test"):

    len_datafile = len(datafile.split("/")[-1])
    filename = datafile.split("/")[-1]
    datadir = datafile[:-len_datafile]

    data = rxte.RXTEData(datafile=datafile, emiddir=datadir, bary=bary)
    times = np.array([p.time for p in data.photons])

    findbursts.extract_bursts(times, t0=data.t0, tseg=5.0, bin_distance=1.5, nbootstrap=nbootstrap,
                              sig_threshold=sig_threshold, froot=froot)

    return


def main():

    find_bursts(datafile, True, sig_threshold, nbootstrap, froot)

    return


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Find bursts in RXTE data! Using Bayesian Blocks! Awesomeness!')


    parser.add_argument("-b", "--find-bursts", action='store_true', dest='find_bursts', required=False,
                        help="Do uou want to search for bursts in a data set?")

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
