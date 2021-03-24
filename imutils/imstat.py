#!/usr/bin/env python
"""
Calculate statistical results for FITS images
"""
import os
import sys
import re
import argparse
import logging
import textwrap
import os.path
from astropy.io import fits
from astropy import stats
import numpy as np


# put parent directory into sys.path
bp = os.path.dirname(os.path.realpath(__file__)).split(os.sep)
modpath = os.sep.join(bp[:-1] + ["lib"])
sys.path.insert(0, modpath)

#  local imports
try:
    import imutils as iu
    import mutils as mu
except ImportError as e:
    logging.error("Import failed: %s", e)
    exit(1)


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """\
            Calculate statistical quantities for image")
                                    """
        ),
        epilog=textwrap.dedent(
            """\
                               """
        ),
    )
    parser.add_argument(
        "fitsfile", nargs="+", metavar="file", help="input fits file(s)"
    )
    parser.add_argument(
        "--quicklook",
        action="store_true",
        help="estimate signal, noise, counts/sec in adus",
    )
    sgroup = parser.add_argument_group(
        "stats", "select statistics and regions" " (exclusive of quicklook)"
    )
    sgroup.add_argument(
        "--region", nargs="+", metavar="reg", help='2d-slicespec: "rows,cols"'
    )
    sgroup.add_argument(
        "--datasec", action="store_true", help="perform stats on DATASEC region"
    )
    sgroup.add_argument(
        "--overscan",
        action="store_true",
        help="perform stats on serial overscan region",
    )
    sgroup.add_argument(
        "--poverscan",
        action="store_true",
        help="perform stats on parllel overscan region",
    )
    sgroup.add_argument(
        "--stats",
        nargs="+",
        metavar="stat",
        help="select from: {mean median stddev min max}",
    )
    # ---------------------------------------------------------------
    hgroup = parser.add_mutually_exclusive_group()
    hgroup.add_argument(
        "--hduname", nargs="+", metavar="idn", help="process HDU list by names"
    )
    hgroup.add_argument(
        "--hduindex", nargs="+", type=int, metavar="idx", help="process HDU list by ids"
    )
    # ---------------------------------------------------------------
    parser.add_argument(
        "--sbias",
        nargs="?",
        const="byrow",
        choices=[
            "mean",
            "median",
            "byrow",
            "byrowe2v",
            "byrowsmooth",
            "byrowsmoothe2v",
        ],
        help="perform bias estimate removal using serial overscan",
    )
    parser.add_argument(
        "--pbias",
        nargs="?",
        const="bycol",
        choices=["mean", "median", "bycol", "bycolfilter", "bycolsmooth"],
        help="perform bias estimate removal using par overscan",
    )
    parser.add_argument(
        "--rstats",
        action="store_true",
        help="use sigma_clipped_stats() for avg,med,std",
    )
    parser.add_argument(
        "--tearing",
        nargs="?",
        metavar="nrows",
        const="datasec",
        help="add tearing metric:" " nrows|divisdero|datasec(default)",
    )
    parser.add_argument(
        "--dipoles", action="store_true", help="add dipole metric to quicklook output"
    )
    parser.add_argument(
        "--threshold",
        nargs=1,
        metavar="thresh",
        type=float,
        help="count number of pixels above threshold",
    )
    parser.add_argument(
        "--info", action="store_true", help="print the info() table summarizing file"
    )
    parser.add_argument(
        "--debug", action="store_true", help="print additional debugging messages"
    )
    parser.add_argument(
        "--noheadings",
        action="store_true",
        default=False,
        help="Don't print column heads for stats",
    )
    return parser.parse_args()


def imstat():
    """main logic:"""
    optlist = parse_args()
    mu.init_logging(optlist.debug)
    mu.init_warnings()
    ncalls.counter = 0
    # begin processing -- loop over files
    for ffile in optlist.fitsfile:
        try:
            hdulist = fits.open(ffile)
        except IOError as ioerr:
            logging.error("IOError: %s", ioerr)
            exit(1)
        if optlist.info:  # just print the image info per file
            hdulist.info()
            continue
        if not optlist.noheadings:  # print filename
            print("#")
            print("# {}".format(os.path.basename(ffile)))
        # Construct a list of the HDU's to work on
        hduids = iu.get_requested_image_hduids(
            hdulist, optlist.hduname, optlist.hduindex
        )
        if optlist.quicklook:
            quicklook(optlist, hduids, hdulist)
        else:
            stats_proc(optlist, hduids, hdulist)
        ncalls.counter = 0  # reset per file, triggers headers


def stats_proc(optlist, hduids, hdulist):
    """print statistics for region according to options
    """
    # Process each HDU in the list "hduids"
    for hduid in hduids:
        hdu = hdulist[hduid]
        name = hdu.name
        if not optlist.sbias and not optlist.pbias:
            pass
        else:
            iu.subtract_bias(optlist.sbias, optlist.pbias, hdu)
        slices = []
        (datasec, soscan, poscan) = iu.get_data_oscan_slices(hdu)
        if optlist.datasec:
            slices.append(datasec)
        if optlist.overscan:
            slices.append(soscan)
        if optlist.poverscan:
            slices.append(poscan)
        if optlist.region:
            for reg in optlist.region:  # if there are regions
                logging.debug("processing %s", reg)
                slice_spec = iu.parse_region(reg)
                if slice_spec:
                    slices.append(slice_spec)
                else:
                    logging.error("skipping region %s", reg)

        if len(slices) == 0:
            stats_print(optlist, hduid, name, hdu.data, None)
        for slice_spec in slices:
            y1, y2 = slice_spec[0].start or "", slice_spec[0].stop or ""
            x1, x2 = slice_spec[1].start or "", slice_spec[1].stop or ""
            reg = "{}:{},{}:{}".format(y1, y2, x1, x2)
            stats_print(optlist, hduid, name, hdu.data[slice_spec], reg)


def stats_print(optlist, sid, name, buf, reg):
    """perform and print the given statistics quantities
    """
    if not optlist.stats:
        optlist.stats = ["mean", "median", "stddev", "min", "max"]

    if optlist.rstats:
        mean_str, median_str, stddev_str = "rmean", "rmedian", "rstddev"
    else:
        mean_str, median_str, stddev_str = "mean", "median", "stddev"

    if not optlist.noheadings and ncalls.counter == 0:
        print("#{:>3s} {:>9s}".format("id", "HDUname"), end="")
        if [stat for stat in optlist.stats if re.match(r"^mea", stat)]:
            print(" {:>9s}".format(mean_str), end="")
        if [stat for stat in optlist.stats if re.match(r"^med", stat)]:
            print(" {:>9s}".format(median_str), end="")
        if [stat for stat in optlist.stats if re.match(r"^std", stat)]:
            print(" {:>8s}".format(stddev_str), end="")
        if [stat for stat in optlist.stats if re.match(r"^min", stat)]:
            print(" {:>8s}".format("min"), end="")
        if [stat for stat in optlist.stats if re.match(r"^max", stat)]:
            print(" {:>8s}".format("max"), end="")
        if reg:
            print("  {:20s}".format("region"), end="")
        print("")  # newline)

    if not optlist.noheadings:
        print(" {:3d} {:>9s}".format(sid, name), end="")

    if optlist.rstats:
        avg, med, std = stats.sigma_clipped_stats(buf, sigma=2.7)
    else:
        avg, med, std = np.mean(buf), np.median(buf), np.std(buf)

    if [stat for stat in optlist.stats if re.match(r"^mea", stat)]:
        print(" {:>9.4g}".format(avg), end="")
    if [stat for stat in optlist.stats if re.match(r"^med", stat)]:
        print(" {:>9.4g}".format(med), end="")
    if [stat for stat in optlist.stats if re.match(r"^std", stat)]:
        print(" {:>8.3g}".format(std), end="")
    if [stat for stat in optlist.stats if re.match(r"^min", stat)]:
        print(" {:>8.4g}".format(np.min(buf)), end="")
    if [stat for stat in optlist.stats if re.match(r"^max", stat)]:
        print(" {:>8.4g}".format(np.max(buf)), end="")

    if reg:
        reg = re.sub(r"^\[*([^\]]*)\]*$", r"\1", reg)
        print("  {:20s}".format(reg), end="")
    print("")  # newline)
    ncalls()  # track call count, acts like static variable)


def quicklook(optlist, hduids, hdulist):
    """print quicklook for hdu's according to options
    """
    try:
        expt = float(hdulist[0].header["EXPTIME"])
    except KeyError as ke:
        emsg = "KeyError: {}".format(ke)
        logging.warninging(emsg)
        emsg = "adu/sec won't be available"
        logging.warninging(emsg)
        expt = 0.0

    # perform and print the given statistics quantities
    # fields are: mean, bias, signal, noise, adu/s
    quick_fields = [
        "mean",
        "bias",
        "signal",
        "noise",
        "adu/sec",
        "eper:s-cte",
        "eper:p-cte",
    ]
    if optlist.tearing:
        quick_fields.append("tearing")
    if optlist.dipoles:
        quick_fields.append("dipoles")
    if optlist.threshold:
        quick_fields.append("threshold")

    for hduid in hduids:
        #
        hdu = hdulist[hduid]
        name = hdu.name

        if not optlist.sbias and not optlist.pbias:
            pass
        else:
            iu.subtract_bias(optlist.sbias, optlist.pbias, hdu)

        # get datasec, serial overscan, parallel overscan as slices
        (datasec, soscan, poscan) = iu.get_data_oscan_slices(hdu)
        if not datasec or not soscan or not poscan:
            logging.error("Could not get DATASEC or overscan specs for %s", name)
            exit(1)

        if optlist.rstats:
            median_str, bias_str, noise_str = "rmedian", "rbias", "rnoise"
        else:
            median_str, bias_str, noise_str = "median", "bias", "noise"

        if not optlist.noheadings and ncalls.counter == 0:
            print("#{:>3s} {:>9s}".format("id", "HDUname"), end="")
            if "mean" in quick_fields:
                print(" {:>9s}".format(median_str), end="")
            if "bias" in quick_fields:
                print(" {:>9s}".format(bias_str), end="")
            if "signal" in quick_fields:
                print(" {:>9s}".format("signal"), end="")
            if "noise" in quick_fields:
                print(" {:>8s}".format(noise_str), end="")
            if "adu/sec" in quick_fields and expt > 0:
                print("{:>9s}".format("adu/sec"), end="")
            if "eper:s-cte" in quick_fields:
                print("{:>9s}".format("s-cte"), end="")
            if "eper:p-cte" in quick_fields:
                print("{:>9s}".format("p-cte"), end="")
            if "tearing" in quick_fields:
                if re.match(r"^data", optlist.tearing):
                    trows = int(datasec[0].stop - 1)
                elif re.match(r"^div", optlist.tearing):
                    trows = 100
                else:
                    trows = int(optlist.tearing)
                print("  {:s}({:>4d}r){:s}".format("tml", trows, "tmr"), end="")
            if "dipoles" in quick_fields:
                print("{:>9s}".format("%dipoles"), end="")
            if "threshold" in quick_fields:
                print("{:>9s}".format("N>thresh"), end="")
            print("")  # newline)

        if not optlist.noheadings:
            print(" {:3d} {:>9s}".format(hduid, name), end="")

        # noise evaluated in smaller region to avoid any gradient effects
        y0 = int(0.6 * datasec[0].start) + int(0.4 * datasec[0].stop)
        y1 = int(0.4 * datasec[0].start) + int(0.6 * datasec[0].stop)
        sx0 = int(0.95 * soscan[1].start) + int(0.05 * soscan[1].stop)
        if optlist.rstats:
            avg, med, std = stats.sigma_clipped_stats(hdu.data[datasec])
            sig_mean = med
            avg, med, std = stats.sigma_clipped_stats(hdu.data[soscan])
            bias_mean = med
            avg, med, std = stats.sigma_clipped_stats(hdu.data[y0:y1, sx0:])
            noise = std
        else:
            sig_mean = np.median(hdu.data[datasec])
            bias_mean = np.median(hdu.data[soscan])
            noise = np.std(hdu.data[y0:y1, sx0:])

        if "mean" in quick_fields:
            print(" {:>9.6g}".format(sig_mean), end="")
        if "bias" in quick_fields:
            print(" {:>9.6g}".format(bias_mean), end="")
        if "signal" in quick_fields:
            signal = sig_mean - bias_mean
            print(" {:>9.6g}".format(signal), end="")
        if "noise" in quick_fields:
            print(" {:>8.3f}".format(noise), end="")
        if "adu/sec" in quick_fields and expt > 0:
            print(" {:>8.3f}".format(float(signal) / expt), end="")
        if "eper:s-cte" in quick_fields:
            logging.debug("s-cte------------------")
            if signal < 5.0 * noise:
                print(" {:>8s}".format("None"), end="")
            else:
                scte = iu.eper_serial(hdu)
                if scte:
                    print(" {:>8.6f}".format(scte), end="")
                else:
                    print(" {:>8s}".format("None"), end="")
        # ---------
        if "eper:p-cte" in quick_fields:
            logging.debug("p-cte------------------")
            if signal < 5.0 * noise:
                print(" {:>8s}".format("None"), end="")
            else:
                pcte = iu.eper_parallel(hdu)
                if pcte:
                    print(" {:>8.6f}".format(pcte), end="")
                else:
                    print(" {:>8s}".format("None"), end="")
        # ---------
        if "tearing" in quick_fields:
            logging.debug("tearing check----------")
            tml, tmr = tearing_metric(hdu.data[datasec], trows)
            print(" {:>5.2f}    {:>5.2f}".format(tml, tmr), end="")
        # ---------
        if "dipoles" in quick_fields:
            logging.debug("dipoles check----------")
            ndipole = count_dipoles(hdu.data[datasec])
            print(
                "{:>9.2f}".format(
                    100.0 * float(2 * ndipole) / (np.size(hdu.data[datasec]))
                ),
                end="",
            )
        # ---------
        if "threshold" in quick_fields:
            logging.debug("threshold check----------")
            print(
                "{:>9d}".format(
                    np.count_nonzero(hdu.data[datasec] > optlist.threshold)
                ),
                end="",
            )
        # ---------
        print("")  # newline)
        ncalls()  # track call count, acts like static variable)


def eper_serial(hdu):
    """
    Given datasec and serial overscan as slices, calculate
    eper using the first ecols=3 columns of serial overscan
    """
    datasec, soscan, poscan = iu.get_data_oscan_slices(hdu)
    ecols = 3  # number of columns used for eper signal
    ncols = datasec[1].stop - datasec[1].start

    # define signal region: mid 20% in y, last 5% in x)
    y0 = int(0.6 * datasec[0].start) + int(0.4 * datasec[0].stop)
    y1 = int(0.4 * datasec[0].start) + int(0.6 * datasec[0].stop)
    sx0 = int(0.05 * datasec[1].start) + int(0.95 * datasec[1].stop)
    sx1 = datasec[1].stop - 1
    logging.debug("s_r=%s[%s:%s,%s:%s]", hdu.name, y0, y1, sx0, sx1)
    s_r = (slice(y0, y1), slice(sx0, sx1))

    # define bias region: same rows as signal, exclude 1st cols in x
    logging.debug(
        "b_r=%s[%s:%s,%s:%s]", hdu.name, y0, y1, soscan[1].start + ecols, soscan[1].stop
    )
    b_r = (slice(y0, y1), slice(soscan[1].start + ecols, soscan[1].stop))
    # define eper region: same rows as signal, 1st ecols cols in x
    logging.debug(
        "e_r=%s[%s:%s,%s:%s]",
        hdu.name,
        y0,
        y1,
        soscan[1].start,
        soscan[1].start + ecols,
    )
    e_r = (slice(y0, y1), slice(soscan[1].start, soscan[1].start + ecols))

    bias_mean = np.mean(hdu.data[b_r])
    # signal
    l_n = np.mean(hdu.data[s_r]) - bias_mean
    logging.debug("l_n=%10.6g", l_n)
    # deferred charge is avg of sum of ecols above bias
    l_nn = np.mean((hdu.data[e_r] - bias_mean).sum(axis=1))
    logging.debug("l_nn=%10.6g", l_nn)
    if l_n > 0.0:
        eper = 1 - (l_nn / (ncols * l_n))
        return eper
    else:
        return None


def eper_parallel_old(hdu):
    """
    Given datasec and parallel overscan as slices, calculate
    eper using the first erows=3 rows of parallel overscan
    """
    datasec, soscan, poscan = iu.get_data_oscan_slices(hdu)
    erows = 3  # number of rows used for eper signal
    nrows = datasec[0].stop - datasec[0].start

    # define signal region: last 1% in y
    y0 = int(0.01 * datasec[0].start) + int(0.99 * datasec[0].stop)
    y1 = datasec[0].stop - 1
    sx0 = datasec[1].start
    sx1 = datasec[1].stop
    logging.debug("s_r=%s[%s:%s,%s:%s]", hdu.name, y0, y1, sx0, sx1)
    s_r = (slice(y0, y1), slice(sx0, sx1))

    # define bias region: same cols as signal, exclude 1st rows in y
    logging.debug(
        "b_r=%s[%s:%s,%s:%s]",
        hdu.name,
        # poscan[0].start + erows,
        poscan[0].stop - erows,
        poscan[0].stop,
        sx0,
        sx1,
    )
    b_r = (slice(poscan[0].stop - erows, poscan[0].stop), slice(sx0, sx1))
    # define eper region: same rows as signal, 1st erows cols in x
    logging.debug(
        "e_r=%s[%s:%s,%s:%s]",
        hdu.name,
        poscan[0].start,
        poscan[0].start + erows,
        sx0,
        sx1,
    )
    e_r = (slice(poscan[0].start, poscan[0].start + erows), slice(sx0, sx1))

    bias_mean_row = np.mean(hdu.data[b_r], axis=0)
    # signal
    l_n = np.mean(np.median(hdu.data[s_r], axis=0) - bias_mean_row)
    logging.debug("l_n=%10.6g", l_n)
    # deferred charge is avg of sum of ecols above bias
    l_nn = np.mean((hdu.data[e_r] - bias_mean_row).sum(axis=0))
    logging.debug("l_nn=%10.6g", l_nn)
    if l_n > 0.0:
        eper = 1 - (l_nn / (nrows * l_n))
        return eper
    else:
        return None


def tearing_metric(buf, trows):
    """
    buf is one segment (w/out pre/over-scan) of an lsst ccd
    return the fraction of pixels in the first and last column, (tml, tmr),
    that are less than 1.5 stddev away from the mean of the
    nearby ~50 pixels in the same row as the pixel being evaluated
    If (tml, tmr) are > O(0.5) then tearing may be present.
    If they are well below 0.5 it is very unlikely
    """
    # left side
    arr = np.mean(buf[10:trows, 3:50], axis=1)
    astd = np.std(buf[10:trows, 3:50], axis=1)
    arr = np.abs((1.0 * buf[10:trows, 0] - arr) / astd)  # col[0] diff in std's
    tml = (1.0 * np.size(arr) - np.searchsorted(arr, 1.5)) / np.size(arr)
    # right side
    arr = np.mean(buf[10:trows, -50:-3], axis=1)
    astd = np.std(buf[10:trows, -50:-3], axis=1)
    arr = np.abs((1.0 * buf[10:trows, -1] - arr) / astd)  # col[-1] diff
    tmr = (1.0 * np.size(arr) - np.searchsorted(arr, 1.5)) / np.size(arr)
    return (tml, tmr)


def count_dipoles(buf):
    """
    buf is one segment (w/out pre/over-scan) of an lsst ccd
    count dipoles via:
    -- use upper 10% of array rows
    -- flatten in column major order
    -- scale array in units of stdev from mean
    -- find adjacent pairs where |A(n)-A(n+1)| > 5
    -- count them
    """
    (nrows, ncols) = np.shape(buf)
    logging.debug("count_dipoles():using subarray [%s:%s,:]", -int(nrows / 10), -1)
    arr = buf[-int(nrows / 10) : -1, :].flatten("F")
    avg, med, std = stats.sigma_clipped_stats(arr)
    logging.debug("clipped stats: avg:%.3g med:%s stdev:%.3g", avg, med, std)
    arr = (arr - avg) / std
    ndipole = 0
    for i in range(0, np.size(arr) - 1):
        if (np.sign(arr[i + 1] * arr[i]) == -1) and abs(arr[i + 1] - arr[i]) > 5:
            ndipole += 1
    return ndipole


def ncalls():
    """maintain a counter
    """
    ncalls.counter += 1


if __name__ == "__main__":
    imstat()
