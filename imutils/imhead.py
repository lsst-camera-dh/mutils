#!/usr/bin/env python
"""Print fits file headers with a variety of options
"""
import os
import sys
import argparse
import logging
from astropy.io import fits

# put parent directory into sys.path
bp = os.path.dirname(os.path.realpath(__file__)).split(os.sep)
modpath = os.sep.join(bp[:-1] + ['lib'])
sys.path.insert(0, modpath)

#  local imports
try:
    import imutils as iu
    import  mutils as mu
except ImportError as e:
    logging.error('Import failed: %s', e)
    exit(1)


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        description="Print headers (Primary+Image+[Others])")
    parser.add_argument("fitsfile", nargs='+',
                        metavar="file", help="input fits file(s)")
    parser.add_argument("--info", action='store_true',
                        help="print the info() table summarizing file")
    parser.add_argument("--hduname", nargs='+',
                        metavar='idn', help="HDU names to print")
    parser.add_argument("--hduindex", nargs='+', type=int,
                        metavar='idx', help="HDU ids to print")
    parser.add_argument("--all",
                        action='store_true', help="print all extensions")
    parser.add_argument("--debug",
                        action='store_true', help="print debugging info")
    return parser.parse_args()
    # return opts


def main():
    """print out headers"""
    optlist = parse_args()
    mu.init_logging(optlist.debug)
    mu.init_warnings()

    # processing -- loop over files
    for ffile in optlist.fitsfile:
        try:
            hdulist = fits.open(ffile)
        except IOError as ioerr:
            emsg = "IOError: {}".format(ioerr)
            logging.error(emsg)
            exit(1)
        if optlist.info:  # just print the image info and exit
            hdulist.info()
        else:
            hduids = iu.get_requested_hduids(hdulist, optlist.hduname,
                                             optlist.hduindex)
            for hduid in hduids:
                print("#--------{}:{}--------".format(
                    hduid, hdulist[hduid].name))
                print(repr(hdulist[hduid].header))
        hdulist.close()

    sys.exit(0)


if __name__ == '__main__':
    main()
