#!/usr/bin/env python
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

# put parent directory into sys.path
bp = os.path.dirname(os.path.realpath(__file__)).split(os.sep)
modpath = os.sep.join(bp[:-1] + ["lib"])
sys.path.insert(0, modpath)

#  local imports
try:
    import trendutils as tu
    import mutils as mu
    import plotutils as pu
    from lsst_camera_data import rafts_of_type
except ImportError as e:
    logging.error("Import failed: %s", e)
    sys.exit(1)

if sys.version_info[0] < 3 or sys.version_info[1] < 7:
    raise Exception("Must be using Python >=3.7")


def parse_args():
    """handle command line"""
    style_list = ["default"] + sorted(
        [
            "fast",
            "ggplot",
            "seaborn-poster",
            "seaborn-notebook",
            "seaborn-darkgrid",
            "fivethirtyeight",
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
        "--input_file", nargs="+", help="XML file with trending data, =>no db query"
    )
    parser.add_argument(
        "--reject", nargs="+", help="filename|regex providing channels to reject"
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
        default="-",
        help="Matplotlib format string (eg. 'o-')",
    )
    parser.add_argument(
        "--title",
        nargs="?",
        metavar="Title (or blank)",
        const="auto",
        help="specify Title String or get auto-generated title",
    )
    parser.add_argument("--style", default="ggplot", required=False, choices=style_list)
    parser.add_argument(
        "--layout", default="portrait", help='"landscape"|"portrait"|"nxm"'
    )
    parser.add_argument("--dpi", type=int, help="set screen dpi")
    parser.add_argument("--fsize", help="set figsize (x-width,y-height)")
    #
    # Time interval specifications
    #
    dgroup = parser.add_mutually_exclusive_group()
    dgroup.add_argument(
        "--interval",
        metavar=("tstart", "tstop"),
        nargs=2,
        action="append",
        help="Pair of date/time specifiers",
    )
    dgroup.add_argument(
        "--start", metavar="tstart", help="Date/time specifier(s)", nargs="+"
    )
    dgroup.add_argument(
        "--stop", metavar="tstop", help="Date/time specifier(s)", nargs="+"
    )
    parser.add_argument(
        "--duration",
        metavar="seconds",
        default=None,
        help="duration [s]|(*s,*m,*h,*d,*w) start/stop spec",
    )
    parser.add_argument(
        "--timebins",
        nargs="?",
        const=0,
        type=int,
        metavar="nBins (blank for autosize)",
        help="retrieve and plot time avg'd bins (esp for long durations)",
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
        "--site", required=False, choices=sites_list, help="Specify trending site",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Print additional debugging info"
    )
    parser.add_argument("--itl", action="store_true", help="limit to ITL devices")
    parser.add_argument("--e2v", action="store_true", help="limit to E2V devices")
    parser.add_argument("--science", action="store_true", help="limit to science rafts")
    parser.add_argument("--corner", action="store_true", help="limit to corner rafts")
    parser.add_argument("--mjd", action="store_true", help="NA, use MJD time axis")
    parser.add_argument(
        "--forceupdate", action="store_true", help="Force update of cached channel file"
    )
    #
    #
    return parser.parse_args()


def main():
    """main logic"""
    # get command args and options
    optlist = parse_args()
    mu.init_logging(optlist.debug)
    mu.init_warnings()

    logging.debug("optlist: %s", optlist)

    # get list of time intervals to process
    intervals = tu.get_unique_time_intervals(
        optlist.start, optlist.stop, optlist.interval, optlist.duration
    )
    if intervals:  # interval accounting
        intcnt = len(intervals)
        inttot = int(sum([t[1] - t[0] for t in intervals]) / 1000)
        tmin = intervals[0][0]
        tmax = intervals[-1][1]
    else:
        logging.error("time interval spec failed")
        sys.exit(1)

    # set up the trending source(chached-on-disk, slac, base, summit etc.)
    if not optlist.input_file:
        tsite = tu.get_trending_server(optlist.site)
        if tsite and tsite["server"]:
            data_url = "http://{}:{}/rest/data/dataserver".format(
                tsite["server"], tsite["port"]
            )
        else:
            logging.error("failed to determine trending server")
            sys.exit(1)

        # access the file with the channel list and update if needed
        if optlist.forceupdate:
            channel_file = tu.update_trending_channels_xml(tsite["name"])
        else:
            channel_file = tu.update_trending_channels_xml(
                tsite["name"], tmin / 1000, tmax / 1000
            )
    else:  # get site and data from input file
        logging.debug("using input file %s", optlist.input_file)
        tsite = tu.init_trending_from_input_xml(optlist.input_file)
        if not tsite:
            logging.error("failed to determine trending server")
            sys.exit(1)
        channel_file = None

    # construct the dict of input channels as {id:path} and store regexes as list
    oflds = dict()  # dict to hold channel information
    regexes = []
    oflds, regexes = tu.parse_channel_sources(optlist.channel_source, channel_file)
    if oflds:
        logging.debug("found %d channels", len(oflds))
    else:
        logging.error("no channels found")
        sys.exit(1)
    if regexes:
        logging.debug("found %d regexes with valid channels", len(regexes))

    # remove channels on the reject list (eg bad RTDs etc)
    rflds, rregexes = tu.parse_channel_sources(optlist.reject, channel_file)
    if rflds:
        logging.debug("found %d channels to reject", len(rflds))
        for rid in rflds.keys():
            if rid in oflds.keys():
                removed = oflds.pop(rid)
                logging.debug("removing %s from channels to process", removed)
            else:
                logging.debug("NOT removing %s from channels to process", rflds[rid])
        logging.debug("%d channels remaining", len(oflds))

    # filter on E2V, ITL, science, corner by removing other types
    rafts_to_reject = []
    if optlist.e2v:
        rafts_to_reject.extend(rafts_of_type["ITL"])
    if optlist.itl:
        rafts_to_reject.extend(rafts_of_type["E2V"])
    if optlist.science:
        rafts_to_reject.extend(rafts_of_type["CORNER"])
    if optlist.corner:
        rafts_to_reject.extend(rafts_of_type["SCIENCE"])
    if rafts_to_reject:
        rids = []
        for chid in oflds:  # loop over paths
            logging.debug("id= %5d  path= %s", int(chid), oflds[chid])
            for raft in set(rafts_to_reject):  # loop over rafts of type
                logging.debug("raft to reject = %s", raft)
                if re.search(f"/{raft}/", oflds[chid]):
                    rids.append(chid)
                    logging.debug("adding %s to channels to reject", oflds[chid])
                    break
            else:
                logging.debug("NOT adding %s to channels to reject", oflds[chid])
        for rid in rids:
            removed = oflds.pop(rid)
            logging.debug("removing %s from channels to process", removed)
    logging.debug("%d channels remaining", len(oflds))

    # now have info needed to query the CCS trending db

    if optlist.match:
        print("#--- Found matching channels:")
        for chid in oflds:
            print("   id: {}  path: {}".format(chid, oflds[chid]))
        sys.exit(0)

    logging.debug("Found matching channels:")
    for chid in oflds:
        logging.debug("id= %5d  path= %s", int(chid), oflds[chid])

    #  Get the trending data either from local saved files or via
    #  trending db queries to the rest service
    if optlist.input_file:
        # get input from files rather than trending service
        # an issue is that the input file need not have the
        # same set of channels or time intervals as requested on command line.
        # The output time intervals will be restricted to the intersection
        # of the intervals present in the input files.
        responses = []
        parser = etree.XMLParser(remove_blank_text=True)
        for ifile in optlist.input_file:
            logging.debug("using %s for input", ifile)
            logging.debug("test for well-formed xml...")
            try:
                tree = etree.parse(ifile, parser)
            except etree.ParseError as e:
                logging.debug("parsing %s failed: %s", ifile, e)
                sys.exit(1)
            except etree.XMLSyntaxError as e:
                logging.debug("parsing %s failed: %s", ifile, e)
                sys.exit(1)
            else:
                logging.debug("successfully parsed %s", ifile)

            logging.debug("appending to responses...")
            responses.append(
                etree.tostring(
                    tree.getroot(),
                    encoding="UTF-8",
                    xml_declaration=True,
                    pretty_print=False,
                )
            )

            logging.debug("deleting the etree")
            del tree

    else:
        # CCS is pre-binned at 5m, 30m, or will rebin on-the-fly
        # default is raw data, timebins triggers stat data
        # query the rest server and place responses into a list
        # join the ids requested as "id0&id=id1&id=id2..." for query
        idstr = "&id=".join(id for id in oflds)
        responses = []
        timebins = 0
        nbins = 0
        if optlist.timebins == 0:  # autosize it
            for ival in intervals:  # only one interval per query allowed
                logging.debug("timebins=0")
                logging.debug("ival[1]= %d, ival[0]= %d", ival[1], ival[0])
                if int((ival[1] - ival[0]) / 1000 / 60) < 5:  # <5m => raw data
                    timebins = None
                elif int((ival[1] - ival[0]) / 1000 / 3600) < 10:  # <10h => 1m bins
                    timebins = int(((ival[1] - ival[0]) / 1000.0) / 60.0)
                elif int((ival[1] - ival[0]) / 1000 / 3600) < 50:  # <50h => 5m bins
                    timebins = int(((ival[1] - ival[0]) / 1000.0) / 300.0)
                else:  # 30m bins
                    timebins = int(((ival[1] - ival[0]) / 1000.0) / 1800.0)
                logging.debug("timebins= %d", timebins)
                if timebins and nbins < timebins:
                    nbins = timebins
        else:
            nbins = optlist.timebins  # is None or an integer

        for ival in intervals:  # only one interval per query allowed
            res = tu.query_rest_server(ival[0], ival[1], data_url, idstr, nbins)
            responses.append(res)
    # Now have the data from trending service

    # Output to stdout a well formed xml tree aggregating the xml received
    # Main use is to save to local file, and re-use for subsequent queries
    # for statistics, plots etc. with subset of channels and time periods
    # Also useful fo debugging and verification of data
    # need to have server info as attribs to get tz correct
    if optlist.xml:
        xml_dec = b'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
        os.write(1, xml_dec)
        datas_str = '<datas {}="{}" {}="{}" {}="{}">\n'.format(
            "trending_server",
            tsite["server"],
            "trending_port",
            tsite["port"],
            "trending_tz",
            tsite["tz"],
        )
        os.write(1, str.encode(datas_str))
        for res in responses:
            root = etree.fromstring(res)
            for data in root.iter("data"):
                os.write(
                    1,
                    etree.tostring(
                        data, encoding="UTF-8", xml_declaration=False, pretty_print=True
                    ),
                )
        try:
            os.write(1, b"</datas>")
        except OSError:
            # 'Broken pipe' OSError when stdout is closed
            pass

        sys.exit(0)

    # Translate the xml responses into internal arrays etc.
    # XML Tree structure looks like this:
    # 1: data [id, path]
    # 2: trendingresult [-]
    #     3: channelmetadata [-]
    #         4: channelmetadatavalue [tstart, tstop, name, value]
    #     3: trendingdata [-]
    #         4: datavalue [name, value]
    #         4: axisvalue [name, value, loweredge, upperedge]
    # where [id, path] could appear multiple times and input time intervals are
    # allowed to overlap
    #
    chanspec = dict()  # where keys are chids, element is also a dict
    chanmd = dict()  # key is chid, elements will be dicts holding arrays
    chandata = dict()  # key is chid, element is list of (time, value) tuples
    datacnt = 0
    for res in responses:
        root = etree.fromstring(res)
        for data in root.iter("data"):
            datacnt += 1
            chid = data.attrib.get("id")
            path = data.attrib.get("path")
            # verify this element's (chid, path) matches the input list
            # logging.debug('id=%s  path=%s', chid, path)
            if chid not in oflds:
                continue
            if path is None or oflds[chid] != path:
                logging.warning(
                    "inputpath(id=%s): %s != %s (xmlpath), using %s",
                    chid,
                    oflds[chid],
                    path,
                    oflds[chid],
                )
                path = oflds[chid]
            # check if chid in
            if chid in chanspec:
                if chanspec[chid]["path"] != path:
                    logging.warning("path mismatch for channel_id= %d", chid)
                    logging.warning(
                        "  %s != %s, skipping....", chanspec[chid]["path"], path
                    )
            else:
                chanspec[chid] = dict()
                chanspec[chid]["path"] = path
                chanspec[chid]["units"] = "none"

            # channelmetadata:
            # each element is a name, value and time interval
            # a name can appear multiple times with distinct time intervals
            # convert to a list, per name, of ordered pairs (value,time)
            # that could be plotted using those points
            #
            if chid not in chanmd:  # check if already exists
                chanmd[chid] = dict()
            # metadata:
            # parse all but only using units for now
            for mdval in data.iter("channelmetadatavalue"):
                if mdval.keys():  # empty sequence is false
                    mdname = mdval.attrib.get("name")  # key
                    mdvalue = mdval.attrib.get("value")  # value
                    mdstart = mdval.attrib.get("tstart")
                    mdstop = mdval.attrib.get("tstop")
                if mdname in chanmd[chid]:
                    chanmd[chid][mdname].append((mdstart, mdvalue))
                    chanmd[chid][mdname].append((mdstop, mdvalue))
                else:  # first assignment
                    chanmd[chid][mdname] = [(mdstart, mdvalue), (mdstop, mdvalue)]
            # trendingdata:
            # extract timestamp, value pairs in axisvalue, datavalue tags
            if chid not in chandata:  # first time
                chandata[chid] = []  # empty list
            for tdval in data.iter("trendingdata"):
                dataval = tdval.find("datavalue")
                if dataval is not None:
                    tvalue = dataval.attrib.get("value")
                else:
                    continue
                axisval = tdval.find("axisvalue")
                if axisval is not None:
                    tstamp = axisval.attrib.get("value")
                else:
                    continue
                # if tstamp is in intervals then append
                for ival in intervals:  # slow, but no other way?
                    if ival[0] < int(tstamp) < ival[1]:
                        chandata[chid].append((tstamp, tvalue))
                        break

    # Done translating the xml responses into internal lists etc.
    # Delete all the raw xml responses
    logging.debug("processed %d xml channel responses", len(responses))
    logging.debug("processed %d uniq channel requests", len(chanspec))
    logging.debug("processed %d total channel queries", datacnt)
    del responses

    # chanspec = dict()  # where keys are chids, values are ccs paths
    # chanmd = dict()  # key is chid, elements will be dicts holding lists
    # chandata = dict() # key is chid, elements are (time, value) pair lists
    # so all responses processed, now have data organized by a set of dicts
    # with the the index on channel id.  Multiple queries for a given channel
    # id are grouped together and there could be duplicate values.
    #
    # To facilitate operating on the data, transform chandat from list[] based
    # (which was easy to append to) to np.array based data.
    chandt = np.dtype({"names": ["tstamp", "value"], "formats": ["int", "float"]})
    trimids = []
    for chid in chanspec:
        path = chanspec[chid]["path"]
        logging.debug("id=%s  path=%s", chid, path)
        for mdname in chanmd[chid]:
            # pick out and process the md's we want
            if mdname == "units" and chanspec[chid]["units"] == "none":
                chanspec[chid]["units"] = chanmd[chid][mdname][-1][1]

        logging.debug("    units=%s", chanspec[chid]["units"])
        # sort and remove duplicates from chandata[chid] where:
        # chandata[chid] = [(t0, v0), (t1, v1), ...., (tn, vn)]
        # and convert to np array
        tmparr = np.array(chandata[chid], dtype=chandt)
        chandata[chid] = np.unique(tmparr)
        logging.debug(
            "    chandata: %d uniq/sorted values from %d entries",
            np.size(chandata[chid]),
            np.size(tmparr),
        )
        del tmparr
        if np.size(chandata[chid]) == 0:  # append chid to trimid list
            logging.debug("%s has no data", chanspec[chid]["path"])
            # arrange to trim empty data
            trimids.append(chid)

    for chid in trimids:
        del chandata[chid]
        del chanmd[chid]
        del chanspec[chid]

    # print to stdout a text dump of the data, in time order per channel
    #
    if optlist.text:
        # print a header for the text
        #
        print("#")
        print("# {}".format(optlist.title))
        print("#")
        print(
            "# CCS trending dump at {}".format(
                dt.datetime.now(gettz()).isoformat(timespec="seconds")
            )
        )
        print(
            "# Data for {} total seconds from {} intervals".format(inttot, intcnt),
            end="",
        )
        print(
            " over {} (h:m:s) from:".format(
                dt.timedelta(seconds=(tmax / 1000 - tmin / 1000))
            )
        )
        print(
            '#     tmin={}: "{}"'.format(
                tmin,
                dt.datetime.fromtimestamp(tmin / 1000, gettz(tsite["tz"])).isoformat(
                    timespec="seconds"
                ),
            )
        )
        print(
            '#     tmax={}: "{}"'.format(
                tmax,
                dt.datetime.fromtimestamp(tmax / 1000, gettz(tsite["tz"])).isoformat(
                    timespec="seconds"
                ),
            )
        )
        print(
            "#{:<{wt}s} {:>{wv}s} {:<{wu}s}  {:<{wp}s}  {:<{wd}s}".format(
                " 'time (ms)'",
                "'value'",
                "'unit'",
                "'channel CCS path'",
                "'iso-8601 Date'",
                wt=13,
                wv=12,
                wu=6,
                wp=30,
                wd=30,
            )
        )
        # loop over all channels sorted on units then path
        for chid in sorted(
            chanspec.keys(), key=lambda x: (chanspec[x]["units"], chanspec[x]["path"])
        ):
            path = chanspec[chid]["path"]
            unitstr = chanspec[chid]["units"]
            if np.size(chandata[chid]) == 0:
                continue
            for (tstamp, value) in chandata[chid]:
                try:
                    date = dt.datetime.fromtimestamp(
                        tstamp / 1000.0, gettz(tsite["tz"])
                    ).isoformat(timespec="milliseconds")
                    print(
                        "{:<{wt}d} {:>{wv}g} {:>{wu}s}   ".format(
                            int(tstamp), float(value), unitstr, wt=14, wv="12.7", wu=6
                        ),
                        end="",
                    )
                    print(
                        "{:<{wp}s}  {:<{wd}s}".format(
                            path, date, wt=14, wv="12.7", wu=6, wp=30, wd=30
                        )
                    )
                except IOError:
                    # 'Broken pipe' IOError when stdout is closed
                    pass

    # print some statistics for each channel
    #
    if optlist.stats:
        # print a header for the stats
        #
        print("#")
        print("# {}".format(optlist.title))
        print("#")
        print(
            "# CCS trending stats at {}".format(
                dt.datetime.now(gettz()).isoformat(timespec="seconds")
            )
        )
        print(
            "# Data for {} total seconds from {} intervals".format(inttot, intcnt),
            end="",
        )
        print(
            " over {} (h:m:s) from:".format(
                dt.timedelta(seconds=(tmax / 1000 - tmin / 1000))
            )
        )
        print(
            '#     tmin="{}"'.format(
                dt.datetime.fromtimestamp(tmin / 1000, gettz(tsite["tz"])).isoformat(
                    timespec="seconds"
                )
            )
        )
        print(
            '#     tmax="{}"'.format(
                dt.datetime.fromtimestamp(tmax / 1000, gettz(tsite["tz"])).isoformat(
                    timespec="seconds"
                )
            )
        )
        print(
            "# {:>4s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>11s}".format(
                "cnt", "mean", "median", "stddev", "min", "max", "d/dt 1/m"
            ),
            end="",
        )
        if optlist.rstats:
            print(
                "{:>8s} {:>8s} {:>8s}  ".format("rmean", "rmedian", "rstddev"), end=""
            )
        print(" {:<{wt}s} {:>{wu}s}".format("path", "units", wt=40, wu=6))

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
                    grad = (
                        60
                        * 1000
                        * (
                            npgrad[-4] * 1.0
                            + npgrad[-3] * 1.0
                            + npgrad[-2] * 1.0
                            + npgrad[-1] * 1.0
                        )
                        / 4.0
                    )
                else:
                    grad = math.nan
                if optlist.rstats:
                    rmean, rmedian, rstd = stats.sigma_clipped_stats(y)
            else:
                avg = med = std = npmin = npmax = 0
                grad = rmean = rmedian = rstd = 0
            try:
                print(
                    "{:>6g} {:>8.4g} {:>8.4g} {:>8.4g} ".format(nelem, avg, med, std,),
                    end="",
                )
                print("{:>8.4g} {:>8.4g} ".format(npmin, npmax), end="")
                print("{:>11.3g} ".format(grad), end="")
                if optlist.rstats:
                    print(
                        "{:>8.4g} {:>8.4g} {:>8.4g}   ".format(rmean, rmedian, rstd),
                        end="",
                    )
                print("{:<{wt}s} {:>{wu}s}".format(path, unitstr, wt=40, wu=6))
            except IOError:
                # 'Broken pipe' IOError when stdout is closed
                pass

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
            nax = nax * len(intervals)  # per interval, per channel

        logging.debug("nax=%d", nax)
        # logging.debug('nrows= %d  ncols=%d', nrows, ncols)

        if not optlist.overlaytime and len(intervals) > 1 and optlist.sharex:
            sharex = False
        else:
            sharex = optlist.sharex

        fig, axes = pu.get_fig_and_axis(
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
                        fmt = "|-"
                    else:
                        fmt = "o-"
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
                        line = ax.plot_date(
                            mds, y, fmt, label=mlabel, tz=gettz(tsite["tz"])
                        )
                        mcolor = line[0].get_color()
                        logging.debug("mcolor= %s", mcolor)
                        labeled = True
                    else:  # no label on later intervals, use saved color
                        line = ax.plot_date(
                            mds, y, fmt, color=mcolor, label=None, tz=gettz(tsite["tz"])
                        )

                    # set x,y-axis label format
                    if not ax.get_ylabel():
                        ax.set_ylabel("{}".format(unit), size="small")
                        ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 5))
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
                            if optlist.timebins:
                                xlabel_str = "{} [{} timebins]".format(
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
            # filename = re.sub(r"(.*).pdf$", r"\1", optlist.saveplot)
            fig.savefig(f"{optlist.saveplot}", dpi=600)
        plt.show()

    # end of main()
    sys.exit(0)


if __name__ == "__main__":
    main()
