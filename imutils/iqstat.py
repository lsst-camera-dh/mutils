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
from astropy.table import Table
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
    sys.exit(1)


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """\
            Calculate statistical quantities for image quality header info
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
        "--bias",
        action="store_true",
        help="auto bias estimate removal by CCD type (itl, e2v)",
    )
    parser.add_argument(
        "--sbias",
        nargs="?",
        const="dbloscan",
        metavar="method",
        help=textwrap.dedent(
            """\
            bias estimate removal using method in:
                {mean,median,byrow,byrowsmooth,[dbloscan],none}
            """
        ),
    )
    parser.add_argument(
        "--pbias",
        nargs="?",
        const="bycolfilter",
        choices=["mean", "median", "bycol", "bycolfilter", "bycolsmooth", "none"],
        help="perform bias estimate removal using par overscan",
    )
    parser.add_argument(
        "--rstats",
        action="store_true",
        help="use sigma_clipped_stats() for avg,med,std",
    )
    parser.add_argument(
        "--info", action="store_true", help="print the info() table summarizing file"
    )
    parser.add_argument(
        "--debug", action="store_true", help="print additional debugging messages"
    )
    parser.add_argument(
        "--human",
        action="store_true",
        default=False,
        help="Don't print column heads for stats",
    )
    return parser.parse_args()


def iqstat():
    """main logic:"""
    optlist = parse_args()
    mu.init_logging(optlist.debug)
    mu.init_warnings()
    ncalls.counter = 0
    # begin processing -- loop over files
    for ffile in optlist.fitsfile:
        try:
            hdulist = fits.open(ffile, mmap=True)
        except IOError as ioerr:
            logging.error("IOError: %s", ioerr)
            sys.exit(1)
        if optlist.info:  # just print the image info per file
            hdulist.info()
            continue
        if optlist.bias:  # auto set [sp]bias, overriding existing
            try:
                optlist.sbias, optlist.pbias = iu.auto_biastype(hdulist)
            except KeyError as kerr:
                logging.error(kerr)
                sys.exit(1)
            except ValueError as verr:
                logging.error(verr)
                sys.exit(1)

        # get basename
        # Construct a list of the HDU's to work on
        hduids = iu.get_requested_image_hduids(
            hdulist, optlist.hduname, optlist.hduindex
        )
        quicklook(optlist, hduids, hdulist)
        ncalls.counter = 0  # reset per file, triggers headers


def quicklook(optlist, hduids, hdulist):
    """print quicklook for hdu's according to options"""
    # perform and print the given statistics quantities

    rows = []
    meta = dict()
    meta["image"] = os.path.basename(hdulist.filename())
    for hduid in hduids:
        #
        hdu = hdulist[hduid]
        name = hdu.name
        names = ["id", "extname", "median", "bias", "signal", "noise"]

        # get so header key:values
        try:
            dateobs = hdulist[0].header["DATE-OBS"]
        except KeyError as ke:
            dateobs = None
            logging.debug(f"DATE-OBS header key not present {ke}")
        try:
            exptime = float(hdulist[0].header["EXPTIME"])
        except KeyError as ke:
            exptime = None
            logging.debug(f"EXPTIME header key not present {ke}")
        try:
            darktime = float(hdulist[0].header["DARKTIME"])
        except KeyError as ke:
            darktime = None
            logging.debug(f"DARKTIME header key not present {ke}")
        if not exptime and not darktime:
            logging.debug("EXPTIME|DARKTIME not in header, adu/sec won't be available")
        try:
            imgtype = hdulist[0].header["IMGTYPE"]
        except KeyError as ke:
            imgtype = None
            logging.debug(f"IMGTYPE header key not present {ke}")

        if not optlist.sbias and not optlist.pbias:
            pass
        else:
            iu.subtract_bias(optlist.sbias, optlist.pbias, hdu)

        # get datasec, serial overscan, parallel overscan as slices
        (datasec, soscan, poscan) = iu.get_data_oscan_slices(hdu)
        if not datasec or not soscan or not poscan:
            logging.error("Could not get DATASEC or overscan specs for %s", name)
            sys.exit(1)

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

        bkgrd = sig_mean - bias_mean

        row = [hduid, name, sig_mean, bias_mean, bkgrd, noise]
        if dateobs:
            meta["dateobs"] = dateobs
        if darktime:
            meta["darktime"] = darktime
            # adus/sec
            row.append(bkgrd / darktime)
            names.append("rate")
        if exptime:
            meta["exptime"] = exptime
            # cte
            eper_serial = iu.eper_serial(hdu)
            if eper_serial:
                row.append(eper_serial)
            else:
                row.append(0.0)
            names.append("scte")

            eper_parallel = iu.eper_parallel(hdu)
            if eper_parallel:
                row.append(eper_parallel)
            else:
                row.append(0.0)
            names.append("pcte")
        if imgtype:
            meta["imgtype"] = imgtype

        ncalls()  # track call count, acts like static variable)
        rows.append(row)

    # names = ["id", "name", "median", "bias", "signal",
    #           "noise", "adus/sec", "eper:s-cte", "eper:p-cte"]

    fmt = ["d", "s", ".6g", ".5g", ".6g", ".3f"]
    tbl = Table(rows=rows, names=names, meta=meta)
    tbl["median"].info.format = "9.6g"
    tbl["bias"].info.format = "9.5g"
    tbl["signal"].info.format = "9.6g"
    tbl["noise"].info.format = "8.3f"
    if darktime:
        tbl["rate"].info.format = "8.3f"
        fmt.append(".3f")
    if exptime:
        tbl["scte"].info.format = "8.6f"
        tbl["pcte"].info.format = "8.6f"
        fmt.append(".6f")
        fmt.append(".6f")

    for key in tbl.meta:
        print(f"# {key}: {tbl.meta[key]}")
    if not optlist.human:  # one item per line
        for row in rows:
            segment = row[1]
            # print("Telemetry: {", end="")
            for label, value, ffmt in zip(names[2:], row[2:], fmt[2:]):
                print(f"Telemetry: {segment}/{label}:{value:{ffmt}}")
            # print("}")
    else:  # fits table
        tbl.pprint_all()


def ncalls():
    """maintain a counter"""
    ncalls.counter += 1


if __name__ == "__main__":
    iqstat()
