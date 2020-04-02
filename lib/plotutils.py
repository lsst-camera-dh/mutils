"""
Common plot functions
"""
import sys
import re
import logging
import math
from astropy import stats
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# local imports
import imutils as iu
import mutils as mu


def get_fig_and_axis(nax, layout, overlay, sharex=False, sharey=False):
    """
    Create the figure and subplot axes based on number of axes requested,
    and options.
    Note plt.subplot() is called with squeeze=False which forces axes
    to be a proper array of Axes objects even if there is only one.

    Returns (fig, axes) tuple
    """
    fsize = plt.rcParams["figure.figsize"]
    if nax == 1 or overlay:
        npcols = nprows = 1
    else:
        # default is portrait, choose #rows > #cols
        npcols = int(math.ceil(math.sqrt(nax)))
        nprows = int(math.ceil(float(nax)/npcols))
        while float(nprows)/float(npcols) < 1.5:
            nprows += 1
            npcols = int(math.ceil(float(nax)/nprows))

        # portrait
        if layout == 'landscape':
            nprows, npcols = npcols, nprows
            if fsize[0] < fsize[1]:
                fsize[0], fsize[1] = fsize[1], fsize[0]

        if re.match(r"([1-9]+)x([1-9]+)$", layout):
            (rstr, cstr) = re.match(
                r"([1-9]+)x([1-9]+)$", layout).groups()
            nprows, npcols = int(rstr), int(cstr)

        if nprows > npcols:  # portrait, make y-dimen bigger
            if fsize[0] > fsize[1]:
                fsize[0], fsize[1] = fsize[1], fsize[0]
            if fsize[1] < 1.2 * fsize[0]:
                fsize[1] = 1.2 * fsize[0]
        if nprows < npcols:  # landscape, make x-dimen bigger
            if fsize[1] > fsize[0]:
                fsize[0], fsize[1] = fsize[1], fsize[0]
            if fsize[0] < 1.6 * fsize[1]:
                fsize[0] = 1.6 * fsize[1]

    plt.rcParams.update({'figure.figsize': fsize})
    logging.debug('subplots layout for nax=%s is nprows=%s x npcols=%s',
                  nax, nprows, npcols)

    fig, axes = plt.subplots(nprows, npcols, sharex=sharex, sharey=sharey,
                             squeeze=False)
    return (fig, axes)


def set_fig_title(title, string_list, fig):
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
        fig.suptitle("{}".format(suptitle), size='medium')


def plot_hdus(optdict: dict, hduids: list,
              hdulist: fits.HDUList, pax: plt.axes):
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
        bias      "overscan" or 'x1:x2' column spec
        btype     mean|median|byrow|byrowsmooth|byrowcol|byrowcolsmooth
        ltype     median|mean|clipped
        steps     default|steps-mid
        offset    mean|median|delta
    hduids: List of hdu ids to work on
    hdulist: A fits HDUList object containing the hdu's
    pax: A matplotlib.axes.Axes instance to contain the line plots
    """
    if optdict['row']:
        map_axis = 0  # first axis (y) converts to scalar
    elif optdict['col']:
        map_axis = 1  # second axis (x) converts to scalar
    else:
        exit(1)
    # Process each HDU in the list "hduids"
    for hduid in hduids:
        hdu = hdulist[hduid]
        try:
            name = hdu.name
        except IndexError as ierr:
            logging.debug('IndexError: %s', ierr)
            logging.debug('using name=%s', hduid)
            name = "{}".format(hduid)
        if optdict['bias']:
            iu.subtract_bias(optdict['bias'], optdict['btype'], hdu)
        (datasec, soscan, poscan) = iu.get_data_oscan_slices(hdu)
        slices = []  # define regions to plot
        for reg in optdict['row'] or optdict['col']:
            logging.debug('processing %s', reg)
            if re.match(r"data", reg):
                slice_spec = datasec
            elif re.match(r"over", reg):
                slice_spec = soscan
            elif re.match(r"pover", reg):
                slice_spec = poscan
            else:
                slice_spec = iu.parse_region(reg)
            if slice_spec is not (None, None):
                slices.append(slice_spec)
            else:
                logging.error('skipping region %s', reg)
        for slice_spec in slices:
            logging.debug('calling line_plot() %s[%s]', name, reg)
            line_plot(slice_spec, hdu.data, optdict['ltype'], optdict['steps'],
                      optdict['offset'], map_axis, name, pax)


def line_plot(slice_spec: tuple, pix: np.ndarray,
              map_type: str, steps: str,
              plot_offset: str, map_axis:
              int, hduname: str, pax: plt.axes):
    """
    make a line plot on the axis given (pax) using the region etc.
    """
    logging.debug('slice_spec: %s', slice_spec)
    logging.debug('shape(pix[slice_spec])=%s', np.shape(pix[slice_spec]))

    #  transform region to a row or column (median, mean, ...)
    logging.debug('map_type=%s', map_type)
    if map_type == 'median' or not map_type:
        line1 = np.median(pix[slice_spec], axis=map_axis)
    elif map_type == 'mean':
        line1 = np.mean(pix[slice_spec], axis=map_axis)
    elif map_type == 'clipped':
        l1_avg, line1, l1_std = stats.sigma_clipped_stats(
            pix[slice_spec], axis=map_axis)
        logging.debug(
            'shape(l1_avg)=%s, shape(pix[slice_spec])=%s, shape(l1_std)=%s',
            np.shape(l1_avg), np.shape(pix[slice_spec]), np.shape(l1_std))
    else:
        logging.error('ltype incorrect or programmer error')
        sys.exit(1)

    #  shift the line (row|column) if needed for plotting
    if plot_offset:
        logging.debug('plot_offset=%s', plot_offset)
        avg, med, std = stats.sigma_clipped_stats(line1)
        if plot_offset == 'mean':
            line1 = line1 - avg
        elif plot_offset == 'median':
            line1 = line1 - med
        elif plot_offset == 'delta':
            # offset by 5 std's per line
            line1 = line1 - med + len(pax.lines) * 5 * std
        else:
            logging.error('invalid offset or programmer error')
            exit(1)

    #  plot the line: N.B. arrays are 0 indexed, fits is 1 indexed
    x = np.arange((slice_spec[1].start or 0) + 1,
                  (slice_spec[1].stop or len(pix[0, :])) + 1)
    y = np.arange((slice_spec[0].start or 0) + 1,
                  (slice_spec[0].stop or len(pix[:, 0])) + 1)
    s = x if map_axis == 0 else y
    slabel = "{}:[{}:{},{}:{}]".format(hduname, x[0], x[-1], y[0], y[-1])
    pax.plot(s, line1, drawstyle="{}".format(steps), label=slabel)


def rotateTickLabels(pax, rotation, which, mode='anchor', ha='right'):
    """
    Rotates the ticklabels of a matplotlib Axes
    Parameters
    ----------
    pax : matplotlib Axes
        The Axes object that will be modified.
    rotation : float
        The amount of rotation, in degrees, to be applied to the labels.
    which : string
        The axis whose ticklabels will be rotated. Valid values are 'x',
        'y', or 'both'.
    mode: string, optional
        The rotation point for the ticklabels. Highly recommended to use
        the default value ('anchor').
    ha : string
        The horizontal alignment of the ticks. Again, recommended to use
        the default ('right').
    Returns
    -------
    None
    """

    if which == 'both':
        rotateTickLabels(pax, rotation, 'x', mode=mode, ha=ha)
        rotateTickLabels(pax, rotation, 'y', mode=mode, ha=ha)
    else:
        if which == 'x':
            axis = pax.xaxis

        elif which == 'y':
            axis = pax.yaxis

        for t in axis.get_ticklabels():
            t.set_horizontalalignment(ha)
            t.set_rotation(rotation)
            t.set_rotation_mode(mode)


def mk_legend(placement: str, nrows: int, handles: list,
              labels: list, pax: plt.axes):
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
    logging.debug('mk_legend():')
    trnc = int(50 / nrows)  # max before truncation
    if len(handles) < 5:  # place in the box
        location = 'best'
        fsize = 'xx-small'
        falpha = 0.5
        bb2a = None
    elif len(handles) < trnc:  # place to right, small
        location = 'upper left'
        fsize = 'x-small'
        falpha = None
        bb2a = (1, 1)
    else:                    # to right, truncate legend
        location = 'upper left'
        fsize = 'xx-small'
        falpha = None
        bb2a = (1, 1)
        if len(handles) >= trnc:
            labels[trnc-1] = "truncated:[{}--{}]".format(
                trnc-1, len(handles))
            labels, handles = labels[:trnc], handles[:trnc]

    if placement.startswith('heur'):
        pass
    elif placement.startswith('ins'):
        location = 'best'
        fsize = 'xx-small'
        # fsize = 3
        bb2a = None
        falpha = 0.4
    elif placement.startswith('outs'):
        location = 'upper left'
        fsize = 'xx-small'
        falpha = None
        bb2a = (1, 1)
    else:
        return
    logging.debug('    placement:%s, nrows:%d', placement, nrows)
    logging.debug('    len(handles):%d, len(labels): %d)',
                  len(handles), len(labels))
    logging.debug('    bbox_to_anchor=%s, fontsize=%s, loc=%s, alpha=%s',
                  bb2a, fsize, location, falpha)
    leg = pax.legend(handles, labels, loc=location,
                     fontsize=fsize, bbox_to_anchor=bb2a,
                     framealpha=falpha)
    if location == 'best':
        leg.set_draggable(True)