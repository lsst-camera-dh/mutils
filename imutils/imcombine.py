#!/usr/bin/env python
"""
Combine images via averaging, median with various options.
The main use is to generate medianed biases, darks, flats.
"""
import sys
import argparse
import logging
import datetime
import textwrap
import os.path
from astropy.io import fits

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

# module scope objects
# none yet


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """
        Image combination
                                    """
        ),
        epilog=textwrap.dedent(
            """
        Given a set of input images, produce an output image
        that is an average or median on the inputs with optional
        bias subtraction and other processing options.
        Note a "--" is often needed to indicate the end of options.
                               """
        ),
    )
    # positional args
    parser.add_argument(
        "fitsfile", nargs="+", metavar="file", help="input fits file(s)"
    )
    # processing options
    #
    parser.add_argument("--ifile", nargs=1, type=str, help="output fits_image name")
    parser.add_argument("--result", nargs=1, type=str, help="output fits_image name")
    # processing choices for image combining
    mgroup = parser.add_mutually_exclusive_group()
    mgroup.add_argument(
        "--median", action="store_true", help="use median to combine images"
    )
    mgroup.add_argument(
        "--mean", action="store_true", help="use mean to combine images"
    )
    mgroup.add_argument(
        "--sigmaclipped",
        const=3.0,
        metavar="nsig",
        nargs="?",
        help="use sigmaclipped mean to combine images",
    )
    mgroup.add_argument(
        "--rank",
        const=50,
        metavar="percentile",
        nargs="?",
        help="use percentile, eg. median=50",
    )
    # hdu selection
    hgroup = parser.add_mutually_exclusive_group()
    hgroup.add_argument(
        "--hduname", nargs="+", default=[], metavar="idn", help="list of HDU names"
    )
    hgroup.add_argument(
        "--hduindex",
        nargs="+",
        type=int,
        default=[],
        metavar="idx",
        help="list of HDU ids",
    )
    # limit output to an ROI
    parser.add_argument("--region", nargs=1, help='2d-slicespec: "rows,cols"')
    # bias subtraction controls
    parser.add_argument("--bimage", nargs=1, help="subtract a bias image")
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
    # define scaling region
    parser.add_argument(
        "--scaling", default=None, nargs="?", help='2d-slicespec: "rows,cols"'
    )
    # additional info
    parser.add_argument(
        "--verbose", action="store_true", help="print progress and process info"
    )
    parser.add_argument(
        "--debug", action="store_true", help="print additional debugging messages"
    )
    return parser.parse_args()


def main():
    """main logic:"""
    optlist = parse_args()
    mu.init_logging(optlist.debug)
    mu.init_warnings()
    verbose = optlist.verbose

    # convert iraf style ROI to slice format
    region = None
    if optlist.region:
        region = iu.parse_region(optlist.region)

    # Prepare scaling region in slice format
    scaling = None
    if optlist.scaling:
        scaling = iu.parse_region(optlist.scaling)

    # build file list
    if optlist.fitsfile:  # input files listed on cmd line
        ifiles = optlist.fitsfile
    elif optlist.ifile:  # input files listed in one or more files
        if not ifiles:
            ifiles = []
        for b in mu.file_to_tokens(optlist.ifile):
            ifiles.extend(b)
    ifiles = sorted(list(set(ifiles)))  # remove duplicates
    if verbose:
        print(f"using {len(ifiles)} input files")
    if optlist.debug:
        logging.debug("input files:")
        for ff in ifiles:
            logging.debug("  %s", ff)

    # get a list of verified images as hdulists
    # prepare input files for use (open and verify)
    if optlist.bimage:  # include in verification
        ifiles.append(optlist.bimage)
    iimages = iu.files_to_hdulists(ifiles, True)
    if optlist.bimage:
        bimage = iimages.pop()  # remove & assign last as bias image
    else:
        bimage = None

    # create output image with primary header and updates
    hdulisto = iu.create_output_hdulist(iimages[0], sys.argv)

    # get all requested hduids
    hduids_to_proc = iu.get_requested_hduids(
        iimages[0], optlist.hduname, optlist.hduindex
    )
    if hduids_to_proc is None:
        logging.error("No HDUs found or requested")
        sys.exit(1)
    # get just requested hduids with image data to be combined
    hduids_to_comb = iu.get_requested_image_hduids(
        iimages[0], optlist.hduname, optlist.hduindex
    )
    if hduids_to_comb is None:
        logging.error("No data HDUs found or requested")
        sys.exit(1)

    # choose the method to combine images
    if optlist.median:
        method = ["median"]
        if len(iimages) < 3:
            logging.warning("image count %d < 3, can only choose mean", len(iimages))
            sys.exit()
    elif optlist.mean:
        method = ["mean"]
    elif optlist.sigmaclipped:
        method = ["sigmaclipped", optlist.sigmaclipped]
        if len(iimages) < 5:
            logging.warning("image count %d < 5, can only choose mean", len(iimages))
            sys.exit()
    elif optlist.rank:
        method = ["rank", optlist.rank]
        if len(iimages) < 5:
            logging.warning(
                "image count %d < 5, can only choose median, mean", len(iimages)
            )
            sys.exit()
    else:
        method = ["median"]  # default

    # prepare the output image hdus using the first image as a template
    for hduid, hdui in enumerate(iimages[0]):
        # hdu image has data to combine
        if hduid in hduids_to_comb:
            logging.debug(f"processing hdu index {hduid}")
            hduo = iu.init_image_hdu(hdui, hdulisto, region)
            # this is the main algorithm/function
            iu.image_combine_hdu(
                iimages,
                hduid,
                method,
                region,
                bimage,
                optlist.sbias,
                optlist.ptype,
                scaling,
                hduo,
            )
            # finish up this hdu
            hduo.update_header()
            dtstr = datetime.datetime.utcnow().isoformat(timespec="milliseconds")
            hduo.add_checksum(dtstr)

        # just append if image hdu has no data (eg. tables etc.)
        if hduid in hduids_to_proc:
            if not isinstance(
                hdui, (fits.ImageHDU, fits.CompImageHDU, fits.PrimaryHDU)
            ):
                # append extensions that contain non-image data
                logging.debug(f"appending hdu index {hduid}")
                hdulisto.append(hdui)

    # write the output file
    hdulisto.writeto(optlist.result[0], overwrite=True)
    sys.exit(0)
    # ------------------------------------


if __name__ == "__main__":
    main()
