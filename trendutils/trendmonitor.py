#!/usr/bin/env python
""" trending snapshot report generator
"""
# code to illustrate parsing of XML files
# importing the required modules
import os
import sys
import re
import argparse
import logging
import time
import requests
from lxml import etree

# put parent directory into sys.path
bp = os.path.dirname(os.path.realpath(__file__)).split(os.sep)
modpath = os.sep.join(bp[:-1] + ["lib"])
sys.path.insert(0, modpath)

#  local imports
try:
    import trendutils as tu
except ImportError as e:
    logging.error("Import failed: %s", e)
    exit(1)

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(description="Display Trending data CCS localdb")
    parser.add_argument(
        "channelFile", nargs="+", help="file specifying channels to obtain"
    )
    parser.add_argument(
        "--debug", action="store_true", help="print additional debugging info"
    )
    parser.add_argument(
        "--date", metavar="datestr", help="where almost any unambiguous format is okay"
    )
    #  note iso-8601 is from unix cmd "date --iso-8601=sconds"
    return parser.parse_args()


def main():
    """main logic"""
    optlist = parse_args()
    if optlist.debug:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
        logging.debug("printing extra debugging...")
    else:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)

    # get time spec for the trending data pull
    if optlist.date:
        logging.debug("datestr=%s", optlist.date)
    t2 = tu.parse_datestr(optlist.date)
    if t2 is None:
        logging.error("Date assignment failed...exiting")
        exit(1)
    t2 = int(t2 * 1000)  # convert to ms
    t1 = t2 - 20 * 1000
    logging.debug("time interval: %d -- %d (%d sec)", t1, t2, (t2 - t1) / 1000)

    trending_server = tu.get_trending_server()
    if trending_server:
        data_url = "http://{}:8080/rest/data/dataserver".format(trending_server)
    else:
        logging.error("trending_server error")
        exit(1)

    # loop over input files containing 1 channel/line
    #
    oflds = []  # list to hold channel information
    for cfile in optlist.channelFile:
        logging.debug("channelFile= %s", cfile)
        try:
            cf = open(cfile, mode="r")
        except OSError as e:
            logging.error("open() failed: %s", e)
            exit(1)

        # populate oflds[] with the channels we want
        for line in cf:
            if re.match(r"^#", line):
                continue
            flds = line.strip().split()
            if len(flds) != 3:
                continue
            if int(flds[0]) == 1:
                oflds.append([flds[1], flds[2]])

    # loop over channels and get xml from restful interface
    # to be parallelized later, sequential for now
    #
    s = requests.Session()
    payload = {"t1": t1, "t2": t2, "flavor": "raw", "n": 1}
    for ch in oflds:
        t_start = time.time()
        uri = "{}/data/?id={}".format(data_url, ch[1])
        try:
            resp = s.get(uri, params=payload)
        except requests.ConnectionError as e:
            logging.error("ConnectionError: %s", e)
            logging.error("channel= %s:%s", ch[0], ch[1])
            logging.error("check status of ssh tunnel to trending server")
        if resp.status_code != 200:
            logging.error(
                "error: invalid response %s from Trending Server", resp.status_code
            )
            exit(1)
        ch.append(resp.content)
        logging.debug("#---->%s: dt= %.3f", ch[0], time.time() - t_start)
    logging.debug("t2=%s", t2)
    s.close()

    print(
        "# CCS trending report at {}".format(
            time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime(t2 / 1000))
        )
    )
    # --- header line
    (chanw, valuew, alarmw, unitw, statew, toffw) = (30, 9, 5, 8, 3, 4)
    hdstr1 = "{:<{w1}s} {:>{w2}s}  {:>{w3}s}".format(
        "# Channel", "Value", "Alo", w1=chanw, w2=valuew, w3=alarmw
    )
    hdstr2 = " {:<{w1}s} {:{w2}s} {:{w3}s}  {:>{w4}s}".format(
        "Ahi", "Units", "St", "Toff", w1=alarmw, w2=unitw, w3=statew, w4=toffw
    )
    print("{} {}".format(hdstr1, hdstr2))

    # loop over channels and get xml from restful interface
    pcols = dict()  # columns to print out for channel
    for ch in oflds:
        logging.debug("processing xml for %s", ch[0])
        # ch[0] channel string
        # ch[1] channel index in TDB
        # ch[2] format(ch[2]) # xml in string
        # loop over channelmetadatavalue elements
        # build columns to print as ordered dictionary
        # columns are: (numbers formated by format attribute)
        # --- channel: ch[1]
        # --- chanval: (last instance of trendingdata(value))
        # --- alarmLow:warnLow
        # --- warnHigh:alarmHigh
        # --- state: NOM|WAR|ALA
        # --- tstamp: convert to seconds UTC, later TAI
        # --- description
        pcols.clear()  # empty out the dictionary and start over
        # channel path goes in first
        pcols["channel"] = ch[0]
        root = etree.fromstring(ch[2])
        # load up meta data name:value pairs
        for mdval in root.iter("channelmetadatavalue"):
            if mdval.keys():
                # check if pcols[name] exists, if so, check tstart and only
                # update if this instance is newer
                mdname = mdval.attrib.get("name")  # key
                mdvalue = mdval.attrib.get("value")  # value
                mdstart = mdval.attrib.get("tstart")
                mdstop = mdval.attrib.get("tstop")
                if mdname not in pcols:  # new entry
                    pcols[mdname] = [mdvalue, mdstart]
                else:  # update only if newer
                    if mdstart > pcols[mdname][1]:
                        pcols[mdname] = [mdvalue, mdstart]
        # iterate over trending data and take the latest value
        for tdval in root.iter("trendingdata"):
            dataval = tdval.find("datavalue")
            if dataval is not None:
                value = dataval.attrib.get("value")
            else:
                continue
            axisval = tdval.find("axisvalue")
            if axisval is not None:
                tstamp = axisval.attrib.get("value")
            else:
                continue
            if "chanval" not in pcols:  # new entry
                pcols["chanval"] = [value, tstamp]
            else:
                if tstamp > pcols["chanval"][1]:
                    pcols["chanval"] = [value, tstamp]

        # print stuff out now
        if "chanval" not in pcols:  # no data, skip
            logging.debug("no data found for %s", pcols["channel"])
            continue
        else:
            chvalstr = "{:>{width}.3g} ".format(
                float(pcols["chanval"][0]), width=valuew
            )
        if len(pcols["channel"]) > chanw:
            chanstr = " ...{:<{width}s}".format(
                pcols["channel"][-(chanw - 4) :], width=(chanw - 4)
            )
        else:
            chanstr = "{:<{width}s}".format(pcols["channel"], width=chanw)
        if "alarmLow" in pcols:
            alarmLowstr = "{:>{width}g}".format(
                float(pcols["alarmLow"][0]), width=alarmw
            )
        else:
            alarmLowstr = "{:>{width}s}".format("", width=alarmw)
        if "alarmHigh" in pcols:
            alarmHighstr = "{:<{width}g}".format(
                float(pcols["alarmHigh"][0]), width=alarmw
            )
        else:
            alarmHighstr = "{:<{width}s}".format("", width=alarmw)
        if "units" in pcols:
            unitstr = "{:<{width}s}".format(pcols["units"][0], width=unitw)
        else:
            unitstr = "{:<{width}s}".format("", width=unitw)
        if "state" in pcols:
            statestr = "{:<{width}.{width}s}".format(pcols["state"][0], width=statew)
        else:
            statestr = "{:<{width}.{width}s}:".format("", width=statew)
        tstampstr = "{:>{width}d}s".format(
            int((int(pcols["chanval"][1]) - t2) / 1000), width=toffw
        )
        print(
            "{0} {1} {2}  {3} {4} {5} {6}".format(
                chanstr,
                chvalstr,
                alarmLowstr,
                alarmHighstr,
                unitstr,
                statestr,
                tstampstr,
            )
        )


if __name__ == "__main__":
    main()
