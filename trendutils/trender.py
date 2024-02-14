#!/usr/bin/env python3
""" trending data app: gets specified channels for requested time period
"""
import re
import os
import argparse
import textwrap
import logging
import sys
import copy
import math
import datetime as dt
from lxml import etree
from dateutil.tz import gettz
import numpy as np
from astropy import stats
import matplotlib.pyplot as plt
import matplotlib.dates as mdate
from astropy.time import Time
import astropy.units as au

# put ../lib directory into sys.path
bp = os.path.dirname(os.path.realpath(__file__)).split(os.sep)
modpath = os.sep.join(bp[:-1] + ["lib"])
sys.path.insert(0, modpath)

#  local imports
try:
    import trendutils as tu
    import mutils as mu
    import plotutils as pu
except ImportError as e:
    logging.error("Import failed: %s", e)
    sys.exit(1)

if sys.version_info[0] < 3 or sys.version_info[1] < 7:
    raise Exception("Must be using Python >=3.7")


def parse_args():
    """
    parse command line
    This is standard argparse stuff which mostly spells out all the
    features of the program but is necessarily long and messy
    """
    style_list = ["default"] + sorted(
        [
            "bmh",
            "classic",
            "dark_background",
            "fast",
            "fivethirtyeight",
            "ggplot",
            "grayscale",
            "seaborn-v0_8",
        ]
    )
    sites_list = tu.get_sites()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """
        Fetch and present trending data from camera database
                                    """
        ),
        epilog=textwrap.dedent(
            """
        This application expects to access the CCS trending database,
        either directly, at lsst-mcm:8080, or via an ssh tunnel which
        must be set up in advance to use localhost:8080.

        Alternately trending data can come from a local file saved
        in an earlier run.  See the "input_file" options.  This allows
        multiple output runs (stats, plots, etc.) using the same trending
        data w/out re-querying the server.

        The "channel_source"(s) are either all files or all patterns (re's).
        Files are best constructed using the sibling application
        "trendingChannels.py" and editing the resuling files to choose
        which channels to activate.  Patterns (regex's) are most useful for
        interactive use.  The "--match" option supports checking which
        channels result from a given set of regex's.

        The interval and start/stop date specifiers can be nearly any
        variant of a full date/time spec as output by the unix
        "date <options>" command.  A suggestd choice is the format
        from "date --iso-8601=seconds".  Formats containing spaces should be
        quoted.  If no timezone is specified, (eg PST, -07:00, etc), the
        local timezone will be used.
                               """
        ),
    )
    # Input args
    parser.add_argument(
        "channel_source", nargs="+", help="filename|regex specifying channels"
    )
    parser.add_argument(
        "--reject", nargs="+", help="filename|regex providing channels to reject"
    )
    tgroup = parser.add_mutually_exclusive_group()  # all to stdout
    tgroup.add_argument(
        "--input_file", nargs="+", help="XML file with trending data, =>no db query"
    )
    tgroup.add_argument(
        "--timebins",
        nargs="?",
        default=None,
        const=0,
        type=int,
        metavar="nBins (blank for autosize)",
        help="retrieve, plot time avg'd bins (esp for long durations)",
    )
    #
    # Output options for text based outputs
    #
    ogroup = parser.add_mutually_exclusive_group()  # all to stdout
    ogroup.add_argument(
        "--xml", action="store_true", help="Print formatted xml from trending to stdout"
    )
    ogroup.add_argument(
        "--text", action="store_true", help="Print (timestamp, value, path) colum text"
    )
    ogroup.add_argument(
        "--stats", action="store_true", help="Print statistics for each channel"
    )
    ogroup.add_argument(
        "--match", action="store_true", help="print list of matching channels and exit"
    )
    parser.add_argument(
        "--rstats",
        action="store_true",
        help="Print additional robust stats per channel",
    )
    #
    # Plotting specifications
    #
    parser.add_argument(
        "--plot", action="store_true", help="produce plots for each channel vs. time"
    )
    parser.add_argument("--logy", action="store_true", help="log y-axis ")
    parser.add_argument("--saveplot", metavar="filename.<pdf|png|..>", help="save as")
    parser.add_argument(
        "--overlay",
        action="store_true",
        help="all channels on a single plot, eg. no subplots",
    )
    parser.add_argument(
        "--overlayunits",
        action="store_true",
        help="channels grouped by units on same (sub)plot",
    )
    parser.add_argument(
        "--overlayregex",
        action="store_true",
        help="channels grouped by regex/units on (sub)plot",
    )
    #
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="normalize chan j as (x-<x>)/std(x) + 5j*std(x)",
    )
    #
    ogroup = parser.add_mutually_exclusive_group()  # all to stdout
    ogroup.add_argument(
        "--overlaytime",
        action="store_true",
        help="time axis is union of all time intervals",
    )
    ogroup.add_argument(
        "--overlaystart",
        action="store_true",
        help="overlay channel data vs time from tstart",
    )
    ogroup.add_argument(
        "--overlaystop",
        action="store_true",
        help="overlay channel data vs time until tstop",
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
    #
    parser.add_argument(
        "--fmt",
        nargs=1,
        metavar="format_str",
        default=None,
        help="Matplotlib format string (eg. 'o-')",
    )
    parser.add_argument(
        "--title",
        nargs="?",
        metavar="Title (or blank)",
        const="auto",
        help="specify Title String or get auto-generated title",
    )
    parser.add_argument(
        "--layout", default="portrait", help='"landscape"|"portrait"|"nxm"'
    )
    parser.add_argument("--dpi", type=int, help="set screen dpi")
    parser.add_argument("--fsize", help="set figsize (x-width,y-height)")
    #
    # Time interval specifications
    #
    igroup = parser.add_mutually_exclusive_group()
    igroup.add_argument(
        "--interval",
        metavar=("tstart", "tstop"),
        nargs=2,
        action="append",
        help="Pair of date/time specifiers",
    )
    igroup.add_argument(
        "--start", metavar="tstart", help="Date/time specifier(s)", nargs="+"
    )
    igroup.add_argument(
        "--stop", metavar="tstop", help="Date/time specifier(s)", nargs="+"
    )
    parser.add_argument(
        "--duration",
        metavar="timespan",
        default=None,
        help="duration *s,*m,*h,*d,*w spec, eg. 1h == 1 hour",
    )
    parser.add_argument(
        "--sharex",
        action="store_true",
        default=True,
        help="use same x-axis limits on all plots",
    )
    #
    # General options
    #
    parser.add_argument(
        "--site",
        required=False,
        choices=sites_list,
        help="Specify trending site",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Print additional debugging info"
    )
    parser.add_argument(
        "--noshow", action="store_true", help="make but don't show the plot"
    )
    parser.add_argument("--itl", action="store_true", help="limit to ITL devices")
    parser.add_argument("--e2v", action="store_true", help="limit to E2V devices")
    parser.add_argument("--science", action="store_true", help="limit to science rafts")
    parser.add_argument("--corner", action="store_true", help="limit to corner rafts")
    parser.add_argument(  # None: use iso date/time
        "--mjd",
        nargs="?",
        # default=argparse.SUPPRESS,
        default=None,
        const=float(math.nan),
        type=float,
        metavar="offset (blank for none)",
        help="use MJD time axis",
    )
    parser.add_argument(
        "--forceupdate",
        action="store_true",
        help="Force update to most recent channel file",
    )
    #
    parser.add_argument("--style", default="ggplot", required=False, choices=style_list)
    #
    return parser.parse_args()


def trender():
    """main logic"""
    # get command args and options
    optlist = parse_args()
    optlistd = vars(optlist)  # dict form of namespace object

    mu.init_logging(optlist.debug)
    mu.init_warnings()

    logging.debug("optlist: %s", optlist)

    # make list of time intervals to process (in milliseconds since the Epoch)
    intervals = tu.get_unique_time_intervals(
        optlist.start, optlist.stop, optlist.interval, optlist.duration
    )
    if intervals:  # interval accounting
        intcnt = len(intervals)
        inttot = int(sum([t[1] - t[0] for t in intervals]) / 1000)  # sec
        tmin = intervals[0][0]
        tmax = intervals[-1][1]
    else:
        logging.error("time interval spec failed")
        sys.exit(1)

    if optlist.timebins is not None:  # using time averaged data
        if optlist.timebins == 0:  # calculate autosize number of bins
            optlist.timebins = tu.autosize_timebins(intervals)

    if optlist.mjd is not None:
        if math.isnan(optlist.mjd):
            mjdoff = math.floor(Time(tmin / 1000, format="unix").mjd)
        else:
            mjdoff = optlist.mjd

    tsite, channel_file = tu.set_trending_source(
        tmin, tmax, optlist.input_file, optlist.site, optlist.forceupdate
    )

    # construct dict oflds of selected channels as {id:path} and store regexes as list
    oflds, regexes = tu.parse_channel_sources(optlist.channel_source, channel_file)
    if not oflds:
        logging.error("no channels found")
        sys.exit(1)

    # remove channels on the reject list (eg bad RTDs etc)
    rflds, rregexes = tu.parse_channel_sources(optlist.reject, channel_file)
    if rflds:
        logging.debug("found %d rejected channels", len(rflds))
        for rid in rflds.keys():
            if rid in oflds.keys():
                removed = oflds.pop(rid)
                logging.debug("removing %s from channels to process", removed)
            else:
                logging.debug("NOT removing %s from channels to process", rflds[rid])
        logging.debug("%d channels remaining", len(oflds))

    # filter on E2V, ITL, science, corner by removing other types
    oflds = tu.filter_raft_types(
        optlist.e2v, optlist.itl, optlist.science, optlist.corner, oflds
    )

    if optlist.debug:
        logging.debug("Found matching channels:")
        for chid in oflds:
            logging.debug("id= %6d  path= %s", int(chid), oflds[chid])

    if optlist.match:  # print matching channels and EXIT
        print("#--- Found matching channels:")
        for chid in oflds:
            print(f"   id: {chid}  path: {oflds[chid]}")
        sys.exit(0)

    # retrieve the requested trending data
    if optlist.input_file:
        responses = tu.get_xml_from_file(optlist.input_file)
    else:
        responses = tu.get_xml_from_trending(tsite, oflds, intervals, optlist.timebins)

    # save data locally and exit
    if optlist.xml:  # dump xml to stdout
        tu.save_to_xml(tsite, responses)
        sys.exit(0)

    # get the data as numpy array
    chanspec, chandata = tu.xml_to_numpy(intervals, oflds, responses)

    # print to stdout a text dump of the data, in time order per channel
    #
    if optlist.text:
        printText(optlist.title, tmin, tmax, inttot, intcnt, tsite, chanspec, chandata)
        sys.exit()

    # print some statistics for each channel
    #
    if optlist.stats:
        printStats(
            optlist.rstats,
            optlist.title,
            tmin,
            tmax,
            inttot,
            intcnt,
            tsite,
            chanspec,
            chandata,
        )

    # Plotting:
    #
    if optlist.plot:
        # make one or more plots of the time series data
        # default is plot per channel per interval
        # option to combine intervals and channels by units
        # and to overlay all

        # update/override some critical parameters
        plt.style.use(optlist.style)
        pu.update_rcparams()

        # figure out how many distinct plots and windows to make
        # subplots layout and shape are determined
        # nax will store the number of actual plots
        # the nxm array of plots may be larger
        #
        nax = len(chanspec)  # default

        if optlist.overlayunits:  # axis set per unit
            unit_map = dict()
            unit_idx = 0  # counts types of units
            # loop over all channels sorted on units then path
            for chid in sorted(
                chanspec, key=lambda x: (chanspec[x]["units"], chanspec[x]["path"])
            ):
                unit = chanspec[chid]["units"]
                if unit not in unit_map:
                    unit_map[unit] = unit_idx
                    unit_idx += 1
            nax = len(unit_map)
            logging.debug("unit_map=%s", unit_map)

        elif optlist.overlayregex:  # axis set per regex and per unit
            # regexes[] is a list of regex's used to select channels
            regex_map = dict()
            axis_idx = 0  # will be the axis index
            # loop over all channels sorted on units then path
            for chid in sorted(
                chanspec, key=lambda x: (chanspec[x]["units"], chanspec[x]["path"])
            ):
                chid_matched = False
                path = chanspec[chid]["path"]
                unit = chanspec[chid]["units"]
                for regex in regexes:
                    if re.search(regex, path):
                        logging.debug("regex_map[%s] matches %s", regex, path)
                        if regex not in regex_map:
                            regex_map[regex] = dict()
                            regex_map[regex][unit] = axis_idx
                            axis_idx += 1
                        elif unit not in regex_map[regex]:
                            regex_map[regex][unit] = axis_idx
                            axis_idx += 1
                        else:  # both regex and unit accounted for
                            pass
                        chid_matched = True
                        break  # found match
                if not chid_matched:
                    logging.error("no regex matches %s", path)
                    sys.exit(1)
            nax = axis_idx  # so now have an axis count, need to re-assign
            # now re-assign axis ids to match command line regexes order
            regex_map_tmp = copy.deepcopy(regex_map)
            aix = 0
            for regex in regexes:
                for unit in sorted(regex_map[regex].keys()):
                    regex_map_tmp[regex][unit] = aix
                    aix += 1
            regex_map = regex_map_tmp
            logging.debug("regex_map=%s", regex_map)
        elif optlist.overlay:
            nax = 1

        if nax == 0:
            logging.error("no data to plot, check inputs?")
            logging.error("try running with --debug")
            sys.exit(1)

        if (
            not optlist.overlaytime
            and not optlist.overlaystart
            and not optlist.overlaystop
        ):
            nax = nax * len(intervals)  # disjoint intervals

        logging.debug("nax=%d", nax)

        if not optlist.overlaytime and len(intervals) > 1 and optlist.sharex:
            sharex = False  # not possible
        else:
            sharex = optlist.sharex

        fig, axes = pu.get_fig_and_axes(
            nax,
            optlist.layout,
            optlist.overlay,
            sharex,
            False,
            optlist.dpi,
            optlist.fsize,
        )

        logging.debug("len(axes)=%d", len(axes))
        logging.debug("axes.shape= %s", axes.shape)
        nrows, ncols = (axes.shape[0], axes.shape[1])

        # loop over data channels and plot them on correct axis
        # chids = list(chanspec.keys())
        chids = sorted(
            chanspec, key=lambda x: (chanspec[x]["units"], chanspec[x]["path"])
        )
        logging.debug("chids=%s", chids)
        unit = None
        for chidx in range(0, len(chids)):  # channels
            chid = chids[chidx]
            unit = chanspec[chid]["units"]
            path = chanspec[chid]["path"]
            tstamp = chandata[chid]["tstamp"]
            mcolor = None
            labeled = False
            for idx in range(0, len(intervals)):
                #
                # choose on which axis to plot
                if optlist.overlayunits:
                    axcnt = unit_map[unit]  # map unit to correct axis
                elif optlist.overlayregex:
                    chid_matched = False
                    for regex in regexes:
                        if re.search(regex, path):
                            logging.debug("regex_map[%s] matches %s", regex, path)
                            axcnt = regex_map[regex][unit]
                            chid_matched = True
                            logging.debug(
                                "using axcnt=%d for regex_map[%s]", axcnt, regex
                            )
                            break
                    if not chid_matched:
                        logging.error("no regex match found for %s", path)
                else:
                    axcnt = chidx
                if not (
                    optlist.overlaytime or optlist.overlaystart or optlist.overlaystop
                ):
                    # stride is number of intervals
                    axcnt = axcnt * len(intervals) + idx

                if optlist.overlay:
                    axcnt = 0
                #
                # now set up this axis
                logging.debug("using axcnt=%d", axcnt)
                ax = np.ravel(axes)[axcnt]
                rowid = int(axcnt / ncols)
                colid = int(axcnt % ncols)
                logging.debug("axcnt= %d  idx= %d", axcnt, idx)
                logging.debug("rowid = %d  colid = %d", rowid, colid)
                ax.grid(True)
                ax.set_frame_on(True)
                ax.get_xaxis().set_visible(True)
                ax.get_yaxis().set_visible(True)
                ax.xaxis.set_tick_params(labelsize="x-small")
                ax.yaxis.set_tick_params(labelsize="x-small")
                if optlist.style == "ggplot":
                    ax.plot([], [])  # consumes the first color (red)
                #
                # mask the tstamps outside of the interval
                mask = (intervals[idx][0] < tstamp) & (tstamp < intervals[idx][1])
                mask_start = intervals[idx][0] / 1000.0
                mask_stop = intervals[idx][1] / 1000.0
                x = chandata[chid]["tstamp"][mask] / 1000.0  # time in seconds
                y = chandata[chid]["value"][mask]
                #
                # deal with point/line format
                if optlist.fmt:
                    fmt = optlist.fmt[0]
                else:
                    if optlist.timebins:
                        fmt = "s-"
                    else:
                        fmt = "-"
                # do the actual plotting
                #
                if not (optlist.overlaystart or optlist.overlaystop):
                    #
                    # convert time axis to matplotlib dates sequence
                    dates = [
                        dt.datetime.fromtimestamp(ts, gettz(tsite["tz"])) for ts in x
                    ]
                    mds = mdate.date2num(dates)
                    if mds.size == 0 and not (
                        optlist.overlay
                        or optlist.overlaytime
                        or optlist.overlayunits
                        or optlist.overlayregex
                    ):
                        #
                        # no data, blank, annotate as empty and skip
                        ax.grid(False)
                        ax.set_frame_on(True)
                        ax.get_xaxis().set_visible(False)
                        ax.get_yaxis().set_visible(False)
                        anno_string = "{}:{} empty".format(chanspec[chid]["path"], idx)
                        ax.annotate(
                            anno_string,
                            xy=(0.03, 0.55),
                            xycoords="axes fraction",
                            horizontalalignment="left",
                            verticalalignment="bottom",
                            fontsize="small",
                        )
                        anno_string = "{} (tstart)".format(
                            dt.datetime.fromtimestamp(
                                intervals[idx][0] / 1000, gettz(tsite["tz"])
                            ).isoformat(timespec="seconds")
                        )
                        ax.annotate(
                            anno_string,
                            xy=(0.03, 0.45),
                            xycoords="axes fraction",
                            horizontalalignment="left",
                            verticalalignment="top",
                            fontsize="small",
                        )
                        continue

                    #   # normalization and shift
                    #   if optlist.normalize:
                    # make label for legend
                    if not labeled:  # label first valid interval, save color
                        mlabel = "{}".format(chanspec[chid]["path"])
                        if optlist.mjd is None:  # use dates
                            line = ax.plot_date(
                                mds, y, fmt, label=mlabel, tz=gettz(tsite["tz"])
                            )
                        else:
                            mjd = Time(mds, format="plot_date").mjd - float(mjdoff)
                            line = ax.plot(mjd, y, fmt, color=mcolor, label=mlabel)
                        mcolor = line[0].get_color()
                        logging.debug("mcolor= %s", mcolor)
                        labeled = True
                    else:  # no label on later intervals, use saved color
                        if optlist.mjd is None:  # use dates
                            line = ax.plot_date(
                                mds,
                                y,
                                fmt,
                                color=mcolor,
                                label=None,
                                tz=gettz(tsite["tz"]),
                            )
                        else:
                            mjd = Time(mds, format="plot_date").mjd - float(mjdoff)
                            line = ax.plot(mjd, y, fmt, color=mcolor, label=None)

                    # set x,y-axis label format
                    if not ax.get_ylabel():
                        ax.set_ylabel("{}".format(unit), size="small")
                        if not optlist.logy:
                            ax.ticklabel_format(
                                axis="y", style="sci", scilimits=(-3, 5)
                            )
                        else:
                            logging.debug("set_yscale(symlog)")
                            ax.set_yscale("symlog")

                    if not ax.get_xlabel():
                        # xlabel and tick labels on bottom plots
                        # only unless multiple intervals
                        if (
                            len(intervals) > 1 and not optlist.overlaytime
                        ) or nax - axcnt - 1 < ncols:
                            xlabel_str = "{} (tstart)".format(
                                dt.datetime.fromtimestamp(
                                    intervals[idx][0] / 1000, gettz(tsite["tz"])
                                ).isoformat(timespec="seconds")
                            )
                            if optlist.mjd is not None:  # nan or a value
                                xlabel_str = "{}  MJD-{}".format(xlabel_str, mjdoff)
                            if optlist.timebins:
                                xlabel_str = "{} [{} bins]".format(
                                    xlabel_str, optlist.timebins
                                )

                            logging.debug("ax.set_xlabel(%s)", xlabel_str)
                            ax.set_xlabel(
                                "{}".format(xlabel_str),
                                position=(0.0, 1e6),
                                size="small",
                                horizontalalignment="left",
                            )

                            ax.tick_params(axis="x", labelbottom=True)
                            # rotate the labels
                            if optlist.mjd is None:
                                for xtick in ax.get_xticklabels():
                                    xtick.set_rotation(30)
                                    xtick.set_horizontalalignment("right")
                        else:
                            ax.tick_params(axis="x", labelbottom=False)
                else:  # overlay start or stop
                    # convert x to duration axis units (s+/-)
                    #
                    if optlist.overlaystart:
                        x = x - mask_start
                    elif optlist.overlaystop:
                        x = x - mask_stop
                    else:
                        logging.error("overlaystart/stop problem")
                        sys.exit(1)

                    mlabel = "{}[{}]".format(chanspec[chid]["path"], idx)
                    line = ax.plot(x, y, fmt, label=mlabel)
                    mcolor = line[0].get_color()
                    logging.debug("mcolor= %s", mcolor)
                    if not ax.get_ylabel():
                        ax.set_ylabel("{}".format(unit), size="small")
                        ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 5))
                    # xlabel for this axis
                    if not ax.get_xlabel():
                        if nax - axcnt - 1 < ncols:
                            if optlist.overlaystart:
                                xstr = "tstart"
                                xid = 0
                            if optlist.overlaystop:
                                xstr = "tstop"
                                xid = 1
                            xlabel_str = "{} ({}[0])".format(
                                dt.datetime.fromtimestamp(
                                    intervals[0][xid] / 1000, gettz(tsite["tz"])
                                ).isoformat(timespec="seconds"),
                                xstr,
                            )
                            if len(intervals) > 1:
                                xlabel_last = "{} ({}[{}])".format(
                                    dt.datetime.fromtimestamp(
                                        intervals[-1][xid] / 1000, gettz(tsite["tz"])
                                    ).isoformat(timespec="seconds"),
                                    xstr,
                                    len(intervals) - 1,
                                )
                                if len(intervals) > 2:
                                    xlabel_last = "...{}".format(xlabel_last)
                                xlabel_str = "{}\n{}".format(xlabel_str, xlabel_last)
                            ax.set_xlabel(
                                "{}".format(xlabel_str),
                                fontsize="small",
                                position=(0.0, 1e6),
                                horizontalalignment="left",
                            )
                            #
                            ax.tick_params(
                                axis="x", labelbottom=True, labelrotation=30.0
                            )
                        else:
                            ax.tick_params(axis="x", labelbottom=False)

        # plot array padded with invisible boxes
        for pcnt in range(nax, nrows * ncols):
            ax = np.ravel(axes)[pcnt]
            ax.grid(False)
            ax.set_frame_on(False)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

        # make the legends for each plot in the array
        # with legend placement/size adjustments
        for pcnt in range(0, nax):
            logging.debug("pcnt= %d, nax= %d", pcnt, nax)
            ax = np.ravel(axes)[pcnt]
            handles, labels = ax.get_legend_handles_labels()
            if not handles or not labels:
                continue
            # sort the labels for easier reading
            # using https://stackoverflow.com/questions/9764298
            labels, handles = (list(t) for t in zip(*sorted(zip(labels, handles))))
            if labels:
                pu.mk_legend(optlist.placement, nrows, handles, labels, ax)

        if optlist.title:  # set the suptitle
            suptitle = optlist.title
            if optlist.title == "auto":
                suptitle = mu.mkglob(
                    [c["path"] for c in list(chanspec.values())], False
                )
            if suptitle:  # mkglob returns None on no solution
                logging.debug("using suptitle=%s", suptitle)
                fig.set_tight_layout(
                    {"h_pad": 0.02, "w_pad": 1.0, "rect": [0, 0, 1, 0.97]}
                )
                fig.suptitle(suptitle, size="medium")
        else:
            fig.set_tight_layout({"h_pad": 0.02, "w_pad": 1.0, "rect": [0, 0, 1, 1]})

        if optlist.saveplot:
            fig.savefig(f"{optlist.saveplot}", dpi=600)

        if optlist.noshow:
            pass
        else:
            plt.show()

    sys.exit(0)
    # end of main()


def printText(
    title: str,
    tmin: int,
    tmax: int,
    inttot: int,
    intcnt: int,
    tsite: dict,
    chanspec: dict,
    chandata: dict,
):
    """
    print out data in text format
    """
    # print a header for the text
    #
    print("#")
    print(f"# {title}")
    print("#")
    tstr = dt.datetime.now(gettz()).isoformat(timespec="seconds")
    print(f"# CCS trending dump at {tstr}")
    tdelta = dt.timedelta(seconds=(tmax / 1000 - tmin / 1000))
    print(
        f"# Data for {inttot} total seconds from {intcnt} intervals",
        f" over {tdelta} (h:m:s) from:",
    )
    tstr = dt.datetime.fromtimestamp(tmin / 1000, gettz(tsite["tz"])).isoformat(
        timespec="seconds"
    )
    print(f'#     tmin={tmin}: "{tstr}"')
    tstr = dt.datetime.fromtimestamp(tmax / 1000, gettz(tsite["tz"])).isoformat(
        timespec="seconds"
    )
    print(f'#     tmax={tmax}: "{tstr}"')
    print(
        f"# {'time (ms)':<13s} {'value':>11s}  {'unit':>5s}  ",
        f"{'channel CCS path':<30s}  {'iso-8601 Date':<30s}",
    )
    # loop over all channels sorted on units then path
    for chid in sorted(
        chanspec.keys(), key=lambda x: (chanspec[x]["units"], chanspec[x]["path"])
    ):
        path = chanspec[chid]["path"]
        unitstr = chanspec[chid]["units"]
        if np.size(chandata[chid]) == 0:
            continue
        for tstamp, value in chandata[chid]:
            try:
                date = dt.datetime.fromtimestamp(
                    tstamp / 1000.0, gettz(tsite["tz"])
                ).isoformat(timespec="milliseconds")
                print(
                    f"{tstamp:<14d} {value:>12.7g} {unitstr:>6s}   ",
                    f"{path:<30s}  {date:<30s}",
                )
            except IOError:
                # 'Broken pipe' IOError when stdout is closed
                pass


def printStats(
    rstats: bool,
    title: str,
    tmin: int,
    tmax: int,
    inttot: int,
    intcnt: int,
    tsite: dict,
    chanspec: dict,
    chandata: dict,
):
    """
    print some statistics for each channel
    """
    #
    # print a header for the stats
    #
    print("#")
    print(f"# {title}")
    print("#")
    tstr = dt.datetime.now(gettz()).isoformat(timespec="seconds")
    print(f"# CCS trending stats at {tstr}")
    tdelta = dt.timedelta(seconds=(tmax / 1000 - tmin / 1000))
    print(
        f"# Data for {inttot} total seconds from {intcnt} intervals",
        f" over {tdelta} (h:m:s) from:",
    )

    tstr = dt.datetime.fromtimestamp(tmin / 1000, gettz(tsite["tz"])).isoformat(
        timespec="seconds"
    )
    print(f'#     tmin="{tstr}"')

    tstr = dt.datetime.fromtimestamp(tmax / 1000, gettz(tsite["tz"])).isoformat(
        timespec="seconds"
    )
    print(f'#     tmax="{tstr}"')

    print(
        f'# {"cnt":>4s} {"mean":>9s} {"median":>8s} {"stddev":>8s} {"min":>8s}',
        f' {"max":>8s} {"d/dt 1/m":>11s}',
        end="",
    )
    if rstats:
        print(f'{"rmean":>8s} {"rmedian":>8s} {"rstddev":>8s}  ', end="")
    print(f'  {"path":<40s} {"units":>6s}')

    # loop over all channels sorted on units then path
    for chid in sorted(
        chanspec.keys(), key=lambda x: (chanspec[x]["units"], chanspec[x]["path"])
    ):
        path = chanspec[chid]["path"]
        unitstr = chanspec[chid]["units"]
        tstamp = chandata[chid]["tstamp"]
        nelem = tstamp.size
        if nelem > 0:
            y = chandata[chid]["value"]
            avg = np.mean(y)
            med = np.median(y)
            std = np.std(y)
            npmin = np.min(y)
            npmax = np.max(y)
            if y.size > 5:
                # silly but better than taking last value
                npgrad = np.gradient(y, tstamp)
                grad = 60 * 1000 * np.sum(npgrad[-4:-1]) / 4.0
            else:
                grad = math.nan
            if rstats:
                rmean, rmedian, rstd = stats.sigma_clipped_stats(y)
        else:
            avg = med = std = npmin = npmax = 0
            grad = rmean = rmedian = rstd = 0
        try:
            print(
                f"{nelem:>6g} {avg:>9.4g} {med:>8.4g} {std:>8.4g} ",
                f"{npmin:>8.3g} {npmax:>8.3g} {grad:>11.3g}",
                end="",
            )
            if rstats:
                print(f"{rmean:>8.4g} {rmedian:>8.4g} {rstd:>8.4g}   ", end="")
            print(f"  {path:<40s} {unitstr:>6s}")
        except IOError:
            # 'Broken pipe' IOError when stdout is closed
            pass


if __name__ == "__main__":
    trender()
