#!/usr/bin/env python3
"""
Display line plots from fits images
"""

import os
import sys
import re
import argparse
import logging
import textwrap
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# put ../lib directory into sys.path
bp = os.path.dirname(os.path.realpath(__file__)).split(os.sep)
modpath = os.sep.join(bp[:-1] + ["lib"])
sys.path.insert(0, modpath)

#  local imports
try:
    import mutils as mu
    import plotutils as pu
    import imutils as iu
except ImportError as e:
    logging.error("Import failed: %s", e)
    sys.exit(1)


def parse_args():
    """handle command line"""
    # style_list = ['default', 'classic'] + sorted(
    #      style for style in plt.style.available if style != 'classic')
    style_list = ["default"] + sorted(
        ["fast", "ggplot", "seaborn-poster", "seaborn-notebook"]
    )
    #
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """\
        Line plots of row/column regions
                                    """
        ),
        epilog=textwrap.dedent(
            """\
        For each file a plot is produced yielding a grid of plots unless the "--overlay"
        option is used which makes just one plot.  Within each file the specified
        hdu's (or all) are plotted for each region given.  Each region is collapsed
        to 1-d by taking the mean, median, clipped median, timeordered, madstd on
        the other axis from the plot.  That is for a row plot, each column of the
        region is mapped to a scalar.  The "--sbias" and "--pbias" options provide
        bias subtraction.  N.B. a "--" is often needed before the file(s) to
        indicate the end of options.  """
        ),
    )
    parser.add_argument(
        "fitsfile", nargs="+", metavar="file", help="input fits file(s)"
    )
    # row or column plot stuff affecting drawing the plots
    pgroup = parser.add_mutually_exclusive_group()
    pgroup.add_argument(
        "--row",
        nargs="*",
        metavar="2d-slicespec",
        help='row plot, fmt: "rowspec, colspec"',
    )
    pgroup.add_argument(
        "--col",
        nargs="*",
        metavar="2d-slicespec",
        help='col plot, fmt: "rowspec, colspec"',
    )
    parser.add_argument(
        "--ltype",
        default="median",
        choices=[
            "median",
            "mean",
            "clipped",
            "timeorder",
            "madstd",
            "stddev",
            "max",
            "min",
        ],
        help="line type: 2d->1d method, default: %(default)s",
    )
    parser.add_argument(
        "--offset",
        nargs="?",
        const="median",
        help='offset choices: "mean", "median", "delta", or <value> (eg 27000)',
    )
    parser.add_argument("--overlay", action="store_true", help="all lines in one plot")

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
        const="bycol",
        choices=["mean", "median", "bycol", "bycolfilter", "bycolsmooth", "none"],
        help="perform bias estimate removal using par overscan",
    )
    # x-axis matplotlib sharex exclusive
    xgroup = parser.add_mutually_exclusive_group()
    xgroup.add_argument(
        "--sharex",
        dest="sharex",
        action="store_true",
        help="Share X-axis range, this is default",
    )
    xgroup.add_argument(
        "--no-sharex",
        dest="sharex",
        action="store_false",
        help="Don't share Y-axis range",
    )
    xgroup.set_defaults(sharex=False)
    # y-axis matplotlib sharey exclusive
    ygroup = parser.add_mutually_exclusive_group()
    ygroup.add_argument(
        "--sharey",
        dest="sharey",
        action="store_true",
        help="Share Y-axis range, this is default",
    )
    ygroup.add_argument(
        "--no-sharey",
        dest="sharey",
        action="store_false",
        help="Don't share Y-axis range",
    )
    ygroup.set_defaults(sharey=True)
    # hdu name|index exclusive
    hgroup = parser.add_mutually_exclusive_group()
    hgroup.add_argument(
        "--hduname", nargs="+", metavar="idn", help="process HDU list by names"
    )
    hgroup.add_argument(
        "--hduindex", nargs="+", type=int, metavar="idx", help="process HDU list by ids"
    )
    parser.add_argument(
        "--info", action="store_true", help="print the info() table summarizing file"
    )
    # title
    parser.add_argument(
        "--title", nargs="?", metavar="Title", const="auto", help="specify Title String"
    )
    #
    pgroup = parser.add_mutually_exclusive_group()  # all to stdout
    pgroup.add_argument(
        "--nolegends",
        dest="placement",
        const="None",
        action="store_const",
        help="Don't place any legends",
    )
    pgroup.add_argument(
        "--insidelegends",
        dest="placement",
        const="inside",
        action="store_const",
        help="Force legends to be inside each plot",
    )
    pgroup.add_argument(
        "--outsidelegends",
        dest="placement",
        const="outside",
        action="store_const",
        help="Force legends to be outside each plot",
    )
    pgroup.set_defaults(placement="heuristic")
    # fluff
    parser.add_argument("--saveplot", metavar="filename.<pdf|png|..>", help="save as")
    parser.add_argument(
        "--logy",
        action="store_true",
        help="log y-axis ",
    )
    parser.add_argument(
        "--xlimits",
        nargs=2,
        type=float,
        required=False,
        help="left right",
    )
    parser.add_argument(
        "--ylimits",
        nargs=2,
        type=float,
        required=False,
        help="lower upper",
    )
    parser.add_argument(
        "--smooth",
        nargs="?",
        type=int,
        const=1,
        metavar="ksize",
        required=False,
        help="smooth lines w/Gaussian1d kernel of size [1]",
    )
    parser.add_argument(
        "--wcs",
        nargs=1,
        metavar="wcs-x",
        required=False,
        help="use wcs x transform",
    )
    parser.add_argument(
        "--layout", default="landscape", help='"landscape"|"portrait"|"nxm"'
    )
    parser.add_argument(
        "--style",
        default="ggplot",
        required=False,
        choices=style_list,
        help="default: %(default)s",
    )
    parser.add_argument(
        "--steps",
        action="store_const",
        required=False,
        const="steps-mid",
        default="default",
        help="Step style: 'steps-mid' best for short lines",
    )
    parser.add_argument("--dpi", type=int, help="set screen dpi")
    parser.add_argument(
        "--debug", action="store_true", help="print additional debugging messages"
    )
    parser.add_argument(
        "--nomemmap",
        action="store_true",
        default=False,
        help="use if memmap fails, default: %(default)s",
    )
    return parser.parse_args()


def implot():
    """main logic:"""
    optlist = parse_args()
    mu.init_logging(optlist.debug)
    mu.init_warnings()
    logging.debug("optlist: %s", optlist)

    # update/override some critical parameters
    plt.style.use(optlist.style)
    pu.update_rcparams()
    # uncomment to use latex
    #  plt.rcParams["text.usetex"] = True
    #  plt.rcParams["font.size"] = 12
    #  plt.rc("text.latex", preamble=r"\usepackage{underscore}")

    fig, axes = pu.get_fig_and_axes(
        len(optlist.fitsfile),
        optlist.layout,
        optlist.overlay,
        optlist.sharex,
        optlist.sharey,
        optlist.dpi,
    )

    fsize = fig.get_size_inches()
    logging.debug("width= %5.2f, height= %5.2f", fsize[0], fsize[1])
    logging.debug("len(axes)=%d", len(axes))
    logging.debug("axes.shape= %s", axes.shape)
    nprows, npcols = (axes.shape[0], axes.shape[1])
    logging.debug("nprows= %d, npcols= %d", nprows, npcols)

    pu.set_fig_title(optlist.title, optlist.fitsfile, fig)
    nfiles = len(optlist.fitsfile)

    # findex = 0
    # for ffile in optlist.fitsfile:
    for findex in range(0, nfiles):
        try:
            hdulist = fits.open(optlist.fitsfile[findex], memmap=optlist.nomemmap)
        except IOError as ioerr:
            logging.error("IOError: %s", ioerr)
            sys.exit(1)
        # info option
        if optlist.info:
            hdulist.info()
            continue

        if optlist.overlay:
            ax = axes[0, 0]
        else:
            logging.debug(
                "ax = np.ravel(axes)[%d + %d]",
                int(findex / npcols) * npcols,
                findex % npcols,
            )
            ax = np.ravel(axes)[int(findex / npcols) * npcols + findex % npcols]

        # construct a list of the HDU's to work on
        hduids = iu.get_requested_image_hduids(
            hdulist, optlist.hduname, optlist.hduindex
        )
        if hduids is None:
            logging.error(
                "No valid HDUs found in %s", optlist.hduname or optlist.hduindex
            )
            sys.exit(1)

        # plot title is truncated filename (w/out path or .fit(s))
        title_str = re.sub(r"^.*/(.*)$", r"\1", optlist.fitsfile[findex])
        title_str = re.sub(r"^(.*)\.fits?(\.fz)*$", r"\1", title_str)
        if npcols < 3:
            title_nchars = 44
            title_fontsize = "x-small"
        else:
            title_nchars = 32
            title_fontsize = "xx-small"
        if len(title_str) > title_nchars:
            title_str = "{}...".format(title_str[:title_nchars])
        else:
            title_str = "{}".format(title_str)
        logging.debug(
            "using title_nchars=%d title_fontsize=%s", title_nchars, title_fontsize
        )
        if optlist.overlay:
            if nfiles > 1:  # trunc'd filename in legend
                ax.plot([], [], " ", label=title_str)
        elif nfiles == 1 and not optlist.title:
            ax.set_title(title_str, fontsize=title_fontsize)
            if optlist.style == "ggplot":
                ax.plot([], [])  # skip the first color
        else:
            ax.set_title(title_str, fontsize=title_fontsize)
            if optlist.style == "ggplot":
                ax.plot([], [])  # skip the first color

        if optlist.xlimits:
            # for ax in np.ravel(axes):
            xbot, xtop = ax.set_xlim(optlist.xlimits[0], optlist.xlimits[1])
            logging.debug("xbot= %.3g xtop= %.3g", xbot, xtop)
        if optlist.ylimits:
            # for ax in np.ravel(axes):
            ybot, ytop = ax.set_ylim(optlist.ylimits[0], optlist.ylimits[1])
            logging.debug("ybot= %.3g ytop= %.3g", ybot, ytop)

        if optlist.logy:
            logging.debug("set_yscale(symlog)")
            ax.set_yscale("symlog")

        # y label depends on offset type
        if not optlist.offset:
            ax.set_ylabel("signal", size="x-small")
        elif optlist.offset == "mean":
            ax.set_ylabel("signal - mean", size="x-small")
        elif optlist.offset == "median":
            ax.set_ylabel("signal - median", size="x-small")
        elif optlist.offset == "delta":
            ax.set_ylabel("signal - mean + 5*j*stdev, j=0,1,..", size="x-small")
        else:
            try:
                offset_value = float(optlist.offset)
            except ValueError as verr:
                logging.error("ValueError: %s", verr)
                logging.error("invalid --offset choice: %s", optlist.offset)
                sys.exit(1)

        # x label
        ax.grid(True)
        if optlist.row is not None:
            ax.set_xlabel("column", size="x-small")
            if optlist.ltype == "timeorder":
                ax.set_xlabel("time ordered series", size="x-small")
        elif optlist.col is not None:
            ax.set_xlabel("row", size="x-small")
            if optlist.ltype == "timeorder":
                ax.set_xlabel("time ordered series", size="x-small")
        else:
            logging.error("must have one of --row or --col")
            sys.exit(1)

        # do the plotting
        pu.plot_hdus(vars(optlist), hduids, hdulist, ax)

        #  done with file, close it
        hdulist.close()

        # end of loop over files

    if optlist.info:  # just print the image info and exit
        sys.exit()

    if optlist.title:  # set the suptitle
        fig.set_tight_layout({"h_pad": 0.50, "w_pad": 1.0, "rect": [0, 0, 1, 0.97]})
    else:
        fig.set_tight_layout({"h_pad": 0.50, "w_pad": 1.0, "rect": [0, 0, 1, 1]})

    # Deal with the legend (ugly)
    # Get list of handles, labels from first plot
    ax = np.ravel(axes)[0]
    handles, labels = ax.get_legend_handles_labels()
    if nfiles == 1 or optlist.overlay:
        ax = np.ravel(axes)[0]
    else:
        ax = np.ravel(axes)[-1]  # put legend in last slot
    pu.mk_legend(optlist.placement, nprows, handles, labels, ax)

    if not optlist.overlay:
        for gidx in range(nfiles, nprows * npcols):
            ax = np.ravel(axes)[int(gidx / npcols) * npcols + gidx % npcols]
            ax.grid(False)
            ax.set_frame_on(False)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

    if optlist.saveplot:
        fig.savefig(f"{optlist.saveplot}", dpi=600)

    plt.show()


if __name__ == "__main__":
    implot()
