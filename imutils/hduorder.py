#!/usr/bin/env python
"""
Create new copy of MEF fits file with ImageHDU order permuted
according to command line options.
"""

import os
import sys
import argparse
from random import shuffle
from astropy.io import fits


# put parent directory into sys.path
bp = os.path.dirname(os.path.realpath(__file__)).split(os.sep)
modpath = os.sep.join(bp[:-1] + ['lib'])
sys.path.insert(0, modpath)

#  local imports
try:
    import mutils as mu
except ImportError as e:
    logging.error('Import failed: %s', e)
    exit(1)

def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        description="Copy MEF fits file with reordered ImageHDUs "
        "All other HDU types have order preserved")
    parser.add_argument("ifitsfile", help="input fits file")
    parser.add_argument("ofitsfile", nargs='?', help="output fits file")
    parser.add_argument("--info", action='store_true',
                        help="print the info() table summarizing file")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--byname", action='store_true',
                       help="write out ImageHDU's sorted by name")
    group.add_argument("--byindex", action='store_true',
                       help="write out ImageHDU's sorted by index")
    group.add_argument("--randomorder", action='store_true',
                       help="write out ImageHDU's in random order")
    group.add_argument("--aslisted", nargs='+', type=int, metavar='idx',
                       help="write out ImageHDU's as listed"
                       " (may need \"--\" to terminate list)")
    return parser.parse_args()


def hdu_writer(opts):
    """print the headers according to the options
    """
    seglist = {}      # dict to hold segment name, index as k,v
    otherlist = {}    # dict to hold non-image HDU's (except for primary)
    pindex = 0

    hdulist = fits.open(opts.ifitsfile)
    # just print the image info with indexes, names, sizes
    if opts.info:
        hdulist.info()
        return 0

    # process primary and Image headers, and other HDUs
    # build dicts etc. to facilitate processing
    for hdu in hdulist:
        index = hdulist.index_of(hdu.name)
        if isinstance(hdu, fits.ImageHDU):
            seglist[hdu.name] = index
        elif isinstance(hdu, fits.PrimaryHDU):
            pindex = index
        else:
            otherlist[hdu.name] = index

    # create the output file and append the primary hdu
    hdulisto = fits.HDUList()
    hdulisto.append(hdulist[pindex])

    # process the HDU's according to command line options
    if opts.byindex:
        # in original index order
        for hdu in hdulist:
            if isinstance(hdu, fits.ImageHDU):
                hdulisto.append(hdu)
    if opts.byname:
        # sorted by name
        for name, index in sorted(seglist.items()):
            hdulisto.append(hdulist[index])

    if opts.randomorder:
        indexlist = seglist.values()
        shuffle(indexlist)
        for idx in indexlist:
            hdulisto.append(hdulist[idx])

    if opts.aslisted:
        for idx in opts.aslisted:
            hdulisto.append(hdulist[idx])

    ids = list(otherlist.values())
    for index in sorted(ids):
        hdu = hdulist[index]
        hdulisto.append(hdu)

    hdulisto.info()
    if opts.ofitsfile:
        hdulisto.writeto(opts.ofitsfile, overwrite=True)
        hdulisto.close()
    hdulist.close()


if __name__ == '__main__':
    optlist = parse_args()
    hdu_writer(optlist)
