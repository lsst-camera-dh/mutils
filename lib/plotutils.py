"""
Common plot functions
"""
import sys
import re
import logging
import math
from astropy import stats
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian1DKernel

# local imports
import imutils as iu
import mutils as mu


def update_rcparams():
    """
    Account for v2.x vs. v3.x differences
    """
    vers = mpl.__version__
    if vers.startswith("3.1"):  # botched hdpi stuff
        # update/override some critical parameters
        logging.debug("found matplotlib version %s", vers)
        plt.rcParams.update({"legend.handlelength": 0.6})
        plt.rcParams.update({"lines.markersize": 2})
        plt.rcParams.update({"figure.dpi": 84})
    elif vers.startswith("2"):
        logging.debug("found matplotlib version %s", vers)
        plt.rcParams.update({"legend.handlelength": 0.6})
        plt.rcParams.update({"lines.markersize": 2})
        plt.rcParams.update({"figure.dpi": 166})
    else:
        logging.debug("found matplotlib version %s", vers)
        plt.rcParams.update({"legend.handlelength": 0.6})
        plt.rcParams.update({"lines.markersize": 2})
        plt.rcParams.update({"figure.dpi": 166})


def get_fig_and_axis(
    nax,
    layout="landscape",
    overlay=False,
    sharex=False,
    sharey=False,
    xdpi=None,
    fsize=None,
):
    """
    Create the figure and subplot axes based on number of axes requested,
    and options.
    Note plt.subplot() is called with squeeze=False which forces axes
    to be a proper array of Axes objects even if there is only one.

    Returns (fig, axes) tuple
    """
    if nax == 1 or overlay:
        npcols = nprows = 1
    else:
        # default is portrait, choose #rows > #cols
        npcols = int(math.ceil(math.sqrt(nax)))
        nprows = int(math.ceil(float(nax) / npcols))
        while float(nprows) / float(npcols) < 1.5:
            nprows += 1
            npcols = int(math.ceil(float(nax) / nprows))

        #
        if layout == "landscape":
            nprows, npcols = npcols, nprows

        if re.match(r"([1-9]+)x([1-9]+)$", layout):
            (rstr, cstr) = re.match(r"([1-9]+)x([1-9]+)$", layout).groups()
            nprows, npcols = int(rstr), int(cstr)

    logging.debug(
        "subplots layout for nax=%s is nprows=%s x npcols=%s", nax, nprows, npcols
    )

    if fsize:
        plt.rcParams["figure.figsize"] = fsize
    else:
        plt.rcParams["figure.figsize"] = (9.0, 5.4)
    fig, axes = plt.subplots(
        nprows,
        npcols,
        sharex=sharex,
        sharey=sharey,
        squeeze=False,
        dpi=xdpi,
    )
    fsize = plt.rcParams["figure.figsize"]
    logging.debug("subplots size is fsize[0]=%s x fsize[1]=%s", fsize[0], fsize[1])

    return (fig, axes)


def set_fig_title(title: str, string_list: list, fig):
    """
    Given the title option string, the file list and the figure,
    generate the suptitle.
    """
    suptitle = None
    if title == "auto":
        logging.debug("using autogenerated title...")
        if string_list:
            suptitle = mu.mkglob(string_list, True)
    else:
        suptitle = title
    logging.debug("using title=%s", suptitle)
    if suptitle:
        fig.suptitle("{}".format(suptitle), size="medium")


def plot_hdus(optdict: dict, hduids: list, hdulist: fits.HDUList, pax: plt.axes):
    """
    For each hdu specified by hduids in hdulist, make a line
    plot on axis pax according the the input parameters

    Parameters
    ----------
    optdict[]: Dictionary with options such as row, col regions, flags
             etc.  Typically the dict version of an argparse namespace.
        key       value
        ---       -----
        col       region spec list
        row       region spec list
        bias      True|False for auto: overrides sbias, pbias
        sbias     mean|median|byrow|byrowsmooth
        pbias     mean|median|bycol|bycolsmooth|lsste2v|lsstitl
        ltype     median|mean|clipped|series
        steps     default|steps-mid
        offset    mean|median|delta
        series    boolean
        wcs       str
        smooth    int
    hduids: List of hdu ids to work on
    hdulist: A fits HDUList object containing the hdu's
    pax: A matplotlib.axes.Axes instance to contain the line plots
    """
    if optdict["row"]:
        map_axis = 0  # first axis (y) converts to scalar
    elif optdict["col"]:
        map_axis = 1  # second axis (x) converts to scalar
    else:
        exit(1)

    if optdict["bias"]:  # auto set [sp]bias, overriding existing
        try:
            optdict["sbias"], optdict["pbias"] = iu.auto_biastype(hdulist)
        except KeyError as kerr:
            logging.error(kerr)
            sys.exit(1)
        except ValueError as verr:
            logging.error(verr)
            sys.exit(1)

    # Process each HDU in the list "hduids"
    for hduid in hduids:
        hdu = hdulist[hduid]
        try:
            name = hdu.name
        except IndexError as ierr:
            logging.debug("IndexError: %s", ierr)
            logging.debug("using name=%s", hduid)
            name = "{}".format(hduid)

        if not optdict["sbias"] and not optdict["pbias"]:
            pass
        else:
            iu.subtract_bias(optdict["sbias"], optdict["pbias"], hdu)

        (datasec, soscan, poscan) = iu.get_data_oscan_slices(hdu)
        wcs = None
        if optdict["wcs"]:
            wcs = WCS(hdu.header, key=optdict["wcs"][0])
            # logging.debug("%s", wcs.printwcs())

        slices = []  # define regions to plot
        for reg in optdict["row"] or optdict["col"]:
            logging.debug("processing %s", reg)
            if re.match(r"data", reg):
                slice_spec = datasec
            elif re.match(r"over", reg):
                slice_spec = soscan
            elif re.match(r"pover", reg):
                slice_spec = poscan
            else:
                slice_spec = iu.parse_region(reg)
            if slice_spec != (None, None):
                slices.append(slice_spec)
            else:
                logging.error("skipping region %s", reg)
        for slice_spec in slices:
            logging.debug(
                "calling line_plot() %s[%s]", name, optdict["row"] or optdict["col"]
            )
            line_plot(
                slice_spec,
                hdu.data,
                optdict["ltype"],
                optdict["steps"],
                optdict["offset"],
                optdict["smooth"],
                wcs,
                map_axis,
                name,
                pax,
            )


def line_plot(
    slice_spec: tuple,
    pix: np.ndarray,
    map_type: str,
    steps: str,
    plot_offset: str,
    smooth: int,
    wcs: WCS,
    map_axis: int,
    hduname: str,
    pax: plt.axes,
):
    """
    make a line plot on the axis given (pax) using the region etc.
    """
    logging.debug("slice_spec: %s", slice_spec)
    logging.debug("shape(pix[slice_spec])=%s", np.shape(pix[slice_spec]))

    #  transform region to a row or column (median, mean, ...)
    logging.debug("map_type=%s", map_type)
    if map_type == "median" or not map_type:
        line1 = np.median(pix[slice_spec], axis=map_axis)
    elif map_type == "mean":
        line1 = np.mean(pix[slice_spec], axis=map_axis)
    elif map_type == "series":
        if map_axis == 0:
            order = "F"
        if map_axis == 1:
            order = "C"
        line1 = pix[slice_spec].flatten(order)
    elif map_type == "clipped":
        l1_avg, line1, l1_std = stats.sigma_clipped_stats(
            pix[slice_spec], axis=map_axis
        )
        logging.debug(
            "shape(l1_avg)=%s, shape(pix[slice_spec])=%s, shape(l1_std)=%s",
            np.shape(l1_avg),
            np.shape(pix[slice_spec]),
            np.shape(l1_std),
        )
    else:
        logging.error("ltype incorrect or programmer error")
        sys.exit(1)
    logging.debug("np.shape(line1)=%s", np.shape(line1))

    #  shift the line (row|column) if needed for plotting
    if plot_offset:
        logging.debug("plot_offset=%s", plot_offset)
        avg, med, std = stats.sigma_clipped_stats(line1)
        if plot_offset == "mean":
            line1 = line1 - avg
        elif plot_offset == "median":
            line1 = line1 - med
        elif plot_offset == "delta":
            # offset by 5 std's per line
            line1 = line1 - med + len(pax.lines) * 5 * std
        else:
            logging.error("invalid offset or programmer error")
            exit(1)

    #  plot the line: np.arrays are 0 indexed, fits is 1 indexed
    xa = np.arange(len(pix[0, :]))
    x = np.arange(xa[slice_spec[1]][0] + 1, xa[slice_spec[1]][-1] + 1 + 1)
    ya = np.arange(len(pix[:, 0]))
    y = np.arange(ya[slice_spec[0]][0] + 1, ya[slice_spec[0]][-1] + 1 + 1)
    logging.debug("xmin=%.2f  xmax=%.2f", x[0], x[-1])
    logging.debug("ymin=%.2f  ymax=%.2f", y[0], y[-1])

    if map_axis == 0:
        s = x
        if wcs:  # can only handle offsets, inversions and rotations mult of 90
            t = np.full(np.shape(s), 0)
            (s, t) = wcs.wcs_pix2world(s, t, 0)
            logging.debug("s[0]=%.2f t[0]=%.2f", s[0], t[0])
            logging.debug("s[-1]=%.2f t[-1]=%.2f", s[-1], t[-1])
            if s[0] == s[-1]:
                s = t
    else:  # map_axis == 1:
        s = y
        if wcs:
            t = np.full(np.shape(s), 0)
            (t, s) = wcs.wcs_pix2world(t, s, 0)
            logging.debug("s[0]=%.2f t[0]=%.2f", s[0], t[0])
            logging.debug("s[-1]=%.2f t[-1]=%.2f", s[-1], t[-1])
            if s[0] == s[-1]:
                s = t

    if map_type == "series":
        s = np.arange(0, len(line1))

    if smooth:  # smooth the line
        kernel = Gaussian1DKernel(smooth)
        line1 = convolve(line1, kernel, boundary="extend")

    slabel = "{}:[{}:{},{}:{}]".format(hduname, y[0], y[-1], x[0], x[-1])
    pax.plot(s, line1, drawstyle="{}".format(steps), label=slabel)
    pax.xaxis.set_tick_params(labelsize="x-small")
    pax.yaxis.set_tick_params(labelsize="x-small")


def mk_legend(placement: str, nrows: int, handles: list, labels: list, pax: plt.axes):
    """
    Generate the legend for the given axis according to options

    Parameters
    ----------
    placement: one of {None, inside, outside, heuristic (default)}
    nrows:     # rows of plots in figure (determine where to truncate)
    handles:   list
    labels:    list
    pax:       axis to use
    """
    logging.debug("mk_legend():")
    trnc = int(39 / nrows)  # max before truncation
    if len(handles) < 5:  # place in the box
        location = "best"
        fsize = "x-small"
        falpha = 0.5
        bb2a = None
    elif len(handles) < trnc:  # place to right, small
        location = "upper left"
        fsize = "x-small"
        falpha = None
        bb2a = (1, 1)
    else:  # to right, truncate legend
        location = "upper left"
        fsize = "xx-small"
        falpha = None
        bb2a = (1, 1)
        if len(handles) >= trnc:
            labels[trnc - 1] = "truncated:[{}--{}]".format(trnc - 1, len(handles))
            labels, handles = labels[:trnc], handles[:trnc]

    if placement.startswith("heur"):
        pass
    elif placement.startswith("ins"):
        location = "best"
        fsize = "xx-small"
        # fsize = 3
        bb2a = None
        falpha = 0.4
    elif placement.startswith("outs"):
        location = "upper left"
        fsize = "xx-small"
        falpha = None
        bb2a = (1, 1)
    else:
        return

    logging.debug("    placement:%s, nrows:%d", placement, nrows)
    logging.debug("    len(handles):%d, len(labels): %d)", len(handles), len(labels))
    logging.debug(
        "    bbox_to_anchor=%s, fontsize=%s, loc=%s, alpha=%s",
        bb2a,
        fsize,
        location,
        falpha,
    )
    leg = pax.legend(
        handles,
        labels,
        loc=location,
        fontsize=fsize,
        bbox_to_anchor=bb2a,
        framealpha=falpha,
    )
    if location == "best":
        try:
            leg.set_draggable(True)
        except AttributeError:
            pass

    return leg
