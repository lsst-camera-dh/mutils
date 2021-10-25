"""
Utility functions for CCS trending access
"""
import sys
import socket
import re
import logging
import collections
import os.path
import time
import stat
from datetime import datetime
import requests
from lxml import etree
import dateutil.parser as dp
from dateutil.tz import gettz
from dateutil.tz import tzutc
from timezone_info import timezone_abbr

# constants? for lack of a better term (default to slac)
# trendnetre = r"134\.79\.[0-9]*\.[0-9]*"
# default_port = "8080"
# default_trending_server = "lsst-mcm.slac.stanford.edu"
# tz_trending = "America/Los_Angeles"
default_req_timeout = 5  # used by requests.*

#  module level dict mapping integer id's to string channel paths
# chid_dict = {"key": "value"}
chid_dict = dict()

# set up known sites for CCS restful servers
# atsccs1.cp.lsst.org has address 139.229.170.33
# ccs-db01.cp.lsst.org has address 139.229.174.2
# lsst-mcm.slac.stanford.edu has address 134.79.209.38
#
sites = dict()
site = None  # module scope variable
#
sites["slac"] = dict()
sites["slac"]["name"] = "slac"
sites["slac"]["netregex"] = r"134\.79\.[0-9]*\.[0-9]*"
sites["slac"]["server"] = "lsst-mcm.slac.stanford.edu"
sites["slac"]["port"] = 8080
sites["slac"]["tz"] = "America/Los_Angeles"
#
sites["lab1"] = dict()
sites["lab1"]["name"] = "lab1"
# sites["lab1"]["netregex"] = r"134.79.100.11"
sites["lab1"]["netregex"] = r"999.999.999.999"  # never matches
sites["lab1"]["server"] = "rddev101.slac.stanford.edu"
sites["lab1"]["port"] = 8080
sites["lab1"]["tz"] = gettz().tzname(datetime.now())
#
sites["lab1tunnel"] = dict()
sites["lab1tunnel"]["name"] = "lab1tunnel"
# sites["lab1tunnel"]["netregex"] = r"134.79.100.11"
sites["lab1tunnel"]["netregex"] = r"999.999.999.999"  # never matches
# sites["lab1tunnel"]["netregex"] = r"[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+"
sites["lab1tunnel"]["server"] = "localhost"
sites["lab1tunnel"]["port"] = 8081
sites["lab1tunnel"]["tz"] = gettz().tzname(datetime.now())
#
sites["ats"] = dict()
sites["ats"]["name"] = "ats"
sites["ats"]["netregex"] = r"139\.229\.170\.[0-9]"
sites["ats"]["server"] = "atsccs1.cp.lsst.org"
sites["ats"]["port"] = 8080
sites["ats"]["tz"] = "UTC"
#
sites["summit"] = dict()
sites["summit"]["name"] = "summit"
sites["summit"]["netregex"] = r"139\.229\.174\.[0-9]"
sites["summit"]["server"] = "ccs-db01.cp.lsst.org"
sites["summit"]["port"] = 8080
sites["summit"]["tz"] = "UTC"
#
sites["comcam"] = dict()
sites["comcam"]["name"] = "comcam"
sites["comcam"]["netregex"] = r"139\.229\.150\.[0-9]"
sites["comcam"]["server"] = "comcam-db01.cp.lsst.org"
sites["comcam"]["port"] = 8080
sites["comcam"]["tz"] = "UTC"
#
# this should be last and matches anything else
# requires a local db or ssh tunnel to the actual server
priv_pat = re.compile(
    r"""
    (^127\.)|(^10\.)|(^172\.1[6-9]\.)|
    (^172\.2[0-9]\.)|(^172\.3[0-1]\.)|(^192\.168\.)
    """,
    re.X,
)
sites["localhost"] = dict()
sites["localhost"]["name"] = "localhost"
# match any address
sites["localhost"]["netregex"] = r"[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+"
# match only private addresses
# sites["localhost"]["netregex"] = priv_pat
sites["localhost"]["server"] = "localhost"
sites["localhost"]["port"] = 8080
sites["localhost"]["tz"] = gettz().tzname(datetime.now())


def set_site(ssite):
    global site
    site = ssite


def print_sites():
    """
    print site keys
    """
    print([*sites])


def get_sites():
    """
    return a list with the site names
    """
    return [*sites]


def get_site(site_name):
    """
    return a dict with the site info
    """
    return sites[site_name]


def deltamtime(pathname):
    """time since last mod"""
    return time.time() - os.stat(pathname)[stat.ST_MTIME]


def get_all_channels(site: str, maxidle: int = -1):
    """get list of channels from the server
    maxidle recovers channels active within maxidle seconds
    """
    trending_server = sites[site]["server"]
    trending_port = sites[site]["port"]

    listpathroot = "/rest/data/dataserver/listchannels?maxIdleSeconds="
    if maxidle:
        listpath = "{}{:d}".format(listpathroot, int(maxidle))
    else:
        listpath = "{}{:s}".format(listpathroot, "-1")
    url = "http://{}:{}{}".format(trending_server, trending_port, listpath)
    logging.debug("url=%s", url)

    # creating HTTP response object from given url
    resp = requests.get(url)
    return resp.content


def parse_datestr(datestr=None):
    """
    convert date string to time in seconds since the epoch

    Returns
    -------
    time in seconds since the epoch (1970) UTC
    """
    if datestr:  # parse the date string almost any format
        logging.debug("parse_datestr(datestr=%s)", datestr)
        try:
            dt = dp.parse(datestr, tzinfos=timezone_abbr)
        except ValueError as e:
            logging.error("ValueError: %s", e)
            logging.error("could not parse the datestring")
            return None
        if dt.tzinfo is None:  # make it tz aware
            dt = dt.astimezone(gettz())  # assume local
    else:
        dt = datetime.now(gettz())  # again local

    return (dt - datetime(1970, 1, 1, tzinfo=tzutc())).total_seconds()


def get_time_interval(startstr, stopstr, duration=None):
    """return the timeinterval boundaries (ms) from datestrings"""
    if startstr:
        t1 = parse_datestr(startstr)
        logging.debug(
            "t1 as UTC: %s",
            datetime.utcfromtimestamp(t1),
        )
        t1 *= 1000

    if stopstr:
        t2 = parse_datestr(stopstr)
        logging.debug(
            "t2 as UTC: %s",
            datetime.utcfromtimestamp(t2),
        )
        t2 *= 1000

    if not duration:
        duration = 600
    else:
        duration = convert_to_seconds(duration)

    # cases for t1, t2
    if startstr and stopstr:
        if t1 > t2:
            logging.error("starttime must be earlier than stop")
            t1 = None
            t2 = None
    elif startstr and not stopstr:  # get durations interval
        t2 = t1 + duration * 1000
    elif not startstr and stopstr:  # get durations interval
        t1 = t2 - duration * 1000
    elif not startstr and not stopstr:
        t2 = int(time.time()) * 1000
        t1 = t2 - duration * 1000
    else:
        t1 = t2 = None

    return (t1, t2)


def get_trending_server(isite: str = None):
    """
    The trending server is one of several in list of sites.
    The default ip address used to access google is checked and
    if it matches one of the known site networks, the server is
    obtained from the list.  If it does not match then the assumption
    is that localhost:8080 is connected via an ssh tunnel and we use that.

    returns dict for trending service
    """
    global site
    if isite:
        site = isite

    if not site:  # try to infer from local network address
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 80))  # connect to google dns server
        res = None
        for ssite in sites:  # break if match
            logging.debug("checking site = %s", ssite)
            res = re.search(sites[ssite]["netregex"], s.getsockname()[0])
            logging.debug(
                "res = re.search(%s, %s)", sites[ssite]["netregex"], s.getsockname()[0]
            )
            logging.debug("res = %s", res)
            if res:
                logging.debug("using site = %s", ssite)
                site = ssite
                break
        s.close()

    if not site:
        logging.error("unable to infer trending server from network, no site given")
        return None

    trending_server = sites[site]["server"]
    trending_port = sites[site]["port"]
    trending_tz = sites[site]["tz"]
    logging.debug(
        "using site: %s: server: %s port: %s tz: %s",
        site,
        trending_server,
        trending_port,
        trending_tz,
    )
    return sites[site]


def get_trending_server_from_xml(xmlfile: str):
    """
    The trending server is one of several in list of sites.
    The default ip address used to access google is checked and
    if it matches one of the known site networks, the server is
    obtained from the list.  If it does not match then the assumption
    is that localhost:8080 is connected via an ssh tunnel and we use that.

    returns dict for trending service
    """
    global site
    trending_server = sites[site]["server"]
    trending_port = sites[site]["port"]
    trending_tz = sites[site]["tz"]
    logging.debug(
        "using site: %s: server: %s port: %s tz: %s",
        site,
        trending_server,
        trending_port,
        trending_tz,
    )
    return sites[site]


def test_trending_server(site: str):
    """
    The trending server is one of several in list of sites.
    The default ip address used to access google is checked and
    if it matches one of the known site networks, the server is
    obtained from the list.  If it does not match then the assumption
    is that localhost:8080 is connected via an ssh tunnel and we use that.

    returns valid service or None on failure
    """
    trending_server = sites[site]["server"]
    trending_port = sites[site]["port"]
    trending_tz = sites[site]["tz"]
    logging.debug(
        "using site: %s: server: %s port: %s tz: %s",
        site,
        trending_server,
        trending_port,
        trending_tz,
    )
    try:
        requests.head(
            "http://{}:{}".format(trending_server, trending_port), timeout=15.0
        )
    except requests.ConnectionError as e:
        logging.error("ConnectionError: %s", e)
        logging.error("check status connection to trending server: %s", trending_server)
        return False
    return True


def geturi(uri):
    """return the xml file from internet or locally"""
    res = re.match(r"http.*", uri)
    if res:
        #  we have a url, use request to return it
        try:
            resp = requests.head(uri, timeout=10.0)
        except requests.ConnectionError as e:
            logging.error("ConnectionError: %s", e)
            return None
        if resp.status_code != 200:
            logging.error("error: invalid response %s from Server", resp.status_code)
            return None
        resp = requests.get(uri)
        return resp.content

    #  it may be a file, try to open as a path
    if not os.path.exists(uri):
        logging.error("%s is not a path", uri)
        return None
    try:
        xmlfile = open(uri, mode="rb")
    except IOError as e:
        logging.error("I/O error(%s): %s", e.errno, e.strerror)
        return None
    xmlstring = xmlfile.read()
    xmlfile.close()
    return xmlstring


def print_channel_content(xmlcontent, ssnames):
    """Walk the tree, find items with subsystem and trending matches
    print out the path and trending-ID
    """
    root = etree.fromstring(xmlcontent)
    for ssname in ssnames:
        for dchan in root.iterfind("datachannel"):
            channel_match = 0
            for md in dchan.iterfind("metadata"):
                if md.attrib.get("name") == "subsystem":
                    if md.attrib.get("value") == ssname:
                        channel_match += 1
                if md.attrib.get("name") == "type":
                    if md.attrib.get("value") == "trending":
                        channel_match += 1

            if channel_match != 2:
                continue
            parr = []
            for pp in dchan.iterfind("path"):
                for pelem in pp.iterfind("pathelement"):
                    if not parr:
                        parr.append("{}".format(pelem.text))
                    else:
                        parr.append("/{}".format(pelem.text))

            for eid in dchan.iterfind("id"):
                parr.append("  {}".format(eid.text))
            print("1 {}".format("".join(parr)))


def print_channel_structure(xmlcontent):
    """print out xml tree structure using algorithm from stackoverflow.com"""
    #
    xml_root = etree.fromstring(xmlcontent)
    raw_tree = etree.ElementTree(xml_root)
    nice_tree = collections.OrderedDict()
    for tag in xml_root.iter():
        path = re.sub(r"\[[0-9]+\]", "", raw_tree.getpath(tag))
        if path not in nice_tree:
            nice_tree[path] = []
            if tag.keys():
                nice_tree[path].extend(
                    attrib for attrib in tag.keys() if attrib not in nice_tree[path]
                )
                for path, attribs in nice_tree.items():
                    indent = int(path.count("/") - 1)
                    print(
                        "{0}{1}: {2} [{3}]".format(
                            "    " * indent,
                            indent,
                            path.split("/")[-1],
                            ", ".join(attribs) if attribs else "-",
                        )
                    )


def update_trending_channels_xml(site, tstart=None, tstop=None):
    """maintain local cache of trending channels in xml file
    arg: tstart -- channels active since tstart (seconds since the epoch)
    """
    logging.debug("update_trending_channels_xml(%s, %s)", tstart, tstop)
    cachedir = "{}/.trender".format(os.environ.get("HOME"))
    channel_file = "{}/.trender/{}_channels.xml".format(os.environ.get("HOME"), site)
    update = True
    # check channel_file exists, get mtime, update if need be
    if not os.path.exists(cachedir):  # make cachdir if not exist
        os.mkdir(cachedir)
    if not os.path.isdir(cachedir):  # is not a directory
        logging.error("%s is not a directory, exiting...", cachedir)
        sys.exit(1)
    #
    # Trigger an update based on whether the interval (tstart,tstop)
    # is covered in the existing cached file based on time of last update
    # and maxIdleSeconds given when the file was fetched.  The attribute
    # 'activeSinceDate' is then maxIdleSeconds prior to last update time
    #
    if os.path.exists(channel_file):  # trigger update based on mtime
        statinfo = os.stat(channel_file)
        mode = statinfo.st_mode
        if not stat.S_IWUSR & mode:  # not writeable
            os.chmod(channel_file, mode | stat.S_IWUSR)
        delta = int(time.time() - statinfo.st_mtime)
        logging.debug("existing cache age: %d (s)", delta)

        chfile = open(channel_file, mode="rb")
        xmlstring = chfile.read()
        chfile.close()

        trending_tz = sites[site]["tz"]
        root = etree.fromstring(xmlstring)
        active_since = root.attrib.get("activeSinceDate")
        if active_since:  # parse, convert and compare to tstart
            xml_start_epoch = parse_datestr(active_since)
            logging.debug(
                "%s channels active_since: %s",
                channel_file,
                datetime.fromtimestamp(xml_start_epoch, gettz(trending_tz)).isoformat(
                    timespec="seconds"
                ),
            )
        # If tstart is inside the interval: [xml_start, last_update]
        # then can guarantee desired channels were being published
        # and hence are already in the cached file.
        if tstart and xml_start_epoch < tstart < statinfo.st_mtime + 43200:
            if not tstop or tstop < (statinfo.st_mtime + 86400):
                update = False

    if update:
        logging.info("updating cached channel_file...")
        if tstart:
            xstart = tstart - 86400  # adjust to 24h earlier
            maxidle = int(time.time() - xstart)
        else:
            maxidle = 3600 * 24 * 7  # give it a week
            xstart = int(time.time() - maxidle)
        xmlstring = get_all_channels(site, maxidle)
        root = etree.fromstring(xmlstring)
        #  attributes in channel_file are from datachannels tag
        #  <datachannels
        #      activeSinceDate="2020-08-21T11:22:35.241-07:00"
        #      host="134.79.209.38"
        #      port="8080">
        #
        active_since = root.attrib.get("activeSinceDate")
        if not active_since:  # needed until service adds this attrib
            active_since_str = datetime.fromtimestamp(
                xstart, gettz(trending_tz)
            ).isoformat(timespec="seconds")
            root.set("activeSinceDate", active_since_str)
            logging.warning(
                "setting activeSinceDate= %s (missing in res)", active_since_str
            )
        else:
            logging.debug("activeSinceDate= %s found", active_since)
        if xmlstring:
            tree = etree.ElementTree(root)
            tree.write(
                channel_file,
                xml_declaration=True,
                encoding="UTF-8",
                pretty_print=False,
                standalone="yes",
            )
    else:
        logging.debug("returning existing channel_file=%s", channel_file)
    return channel_file


def convert_to_seconds(duration_str) -> int:
    """
    for (s)econds, (m)inutes, (h)ours, (d)ays, (w)eeks
    return duration in seconds
    """

    seconds = 0
    if re.match(r"[0-9]+$", duration_str):
        seconds = int(duration_str)
    elif re.match(r"[0-9]+s$", duration_str):
        seconds = int(duration_str[:-1])
    elif re.match(r"[0-9]+m$", duration_str):
        seconds = 60 * int(duration_str[:-1])
    elif re.match(r"[0-9]+h$", duration_str):
        seconds = 3600 * int(duration_str[:-1])
    elif re.match(r"[0-9]+d$", duration_str):
        seconds = 84600 * int(duration_str[:-1])
    elif re.match(r"[0-9]+w$", duration_str):
        seconds = 7 * 84600 * int(duration_str[:-1])

    return seconds


def long_substr(data: list) -> str:
    """
    https://stackoverflow.com/questions/2892931/\
        longest-common-substring-from-more-than-two-strings-python#
    """
    substr = ""
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0]) - i + 1):
                if j > len(substr) and all(data[0][i : i + j] in x for x in data):
                    substr = data[0][i : i + j]
    return substr


def query_rest_server(ts1, ts2, data_url, idstr, nbins):
    """get xml from restful interface for the requested channels
    with a single request
    inputs:
        t1, t2 are start/stop time in ms
        data_url is "server-url:port/default-path"
        idstr is all channel ids as '&id='.join(id for id in oflds)
        nbins: request raw (=None) or binned data from DB
    output: raw xml response from the service is returned
    """
    s = requests.Session()
    if nbins is None:  # raw data
        options = {"t1": int(ts1), "t2": int(ts2), "flavor": "raw", "n": 1}
    else:  # CCS stat data
        options = {"t1": int(ts1), "t2": int(ts2), "flavor": "stat", "n": int(nbins)}
    uri = "{}/data/?id={}".format(data_url, idstr)
    t_start = time.time()
    try:
        resp = s.get(uri, params=options)
    except requests.ConnectionError as e:
        logging.error("ConnectionError: %s", e)
        logging.error("check status of ssh tunnel to trending server")
    if resp.status_code != 200:
        logging.error("invalid response %s from Trending Server", resp.status_code)
        return None
    # logging.debug('URL=%s', resp.url)
    logging.debug("channels: %s", re.sub(r"(id=)?([0-9]+)&*", r"\2 ", idstr))
    logging.debug("dt=%.3f seconds", (time.time() - t_start))
    s.close()
    return resp.content


#   def get_unique_time_intervals(optlist):
def get_unique_time_intervals(starts=None, stops=None, intervalarr=None, duration=None):
    """
    Input: Command line options defining a set of intervals
           using pairs or start/stop with duration,
           empty input will return 1 defaut interval

    The set of intervals as ordered pairs in seconds are processed
    to merge overlapping periods yielding distinct intervals.

    Output: A ordered list of non-overlapping periods are returned
            as [[t00,t01], [t10,t11], ...,[tn0,tn1]]
    """
    intervals = []
    if starts:
        for start in starts:
            (t1, t2) = get_time_interval(start, None, duration)
            intervals.append([t1, t2])
    elif stops:
        for stop in stops:
            (t1, t2) = get_time_interval(None, stop, duration)
            intervals.append([t1, t2])
    elif intervalarr:
        for interval in intervalarr:
            (t1, t2) = get_time_interval(interval[0], interval[1])
            intervals.append([t1, t2])
    else:
        (t1, t2) = get_time_interval(None, None, duration)
        intervals.append([t1, t2])

    for interval in intervals:
        if interval[0] is None or interval[1] is None:
            logging.error("Date assignment failed")
            return None

    i = 0
    for interval in intervals:
        logging.debug(
            "time interval[%d] (before merge): %d -- %d (%d sec)",
            i,
            interval[0],
            interval[1],
            (interval[1] - interval[0]) / 1000,
        )
        i += 1

    # merge overlaps to generate list of distinct intervals
    intervals.sort()  # sorts so that intervals[i][0] <= intervals[i+1][0]
    i = 1
    while i < len(intervals):  # loop over pairs of intervals
        if intervals[i - 1][1] >= intervals[i][0]:
            intervals[i][0] = intervals[i - 1][0]  # move left edge down
            if intervals[i - 1][1] > intervals[i][1]:
                intervals[i][1] = intervals[i - 1][1]  # move right edge up
            del intervals[i - 1]  # delete the 1st of the pair
        else:
            i += 1  # no overlap so move to next pair

    i = 0
    for interval in intervals:
        logging.debug(
            "time interval[%d] (after merge): %d -- %d (%d sec)",
            i,
            interval[0],
            interval[1],
            (interval[1] - interval[0]) / 1000,
        )
        i += 1

    return intervals


def get_channel_dict(channel_file):
    """
    build a channel dictionary for all channels
    struction of channel_file is:
    0: datachannels [-]
      1: datachannel [-]
        2: path [-]
          3: pathelement [-]
      2: id [-]
      2: metadata [name, value]
    """
    logging.debug("building full channel dictionary from %s", channel_file)
    cdict = dict()
    tree = etree.parse(channel_file)
    root = tree.getroot()
    for dchan in root.iterfind("datachannel"):
        chid = dchan.find("id").text
        # build path
        parr = []  # list to hold path elements
        pp = dchan.find("path")
        for pe in pp.iterfind("pathelement"):
            if pe.text:  # work-around for problem xxx
                parr.append(pe.text)
        path = "/".join(parr)
        if path and chid:  # create entry in dict
            cdict[chid] = path
    logging.debug("channel dict contains %d active channels", len(cdict))
    del tree
    return cdict


def parse_channel_sources(sources: list, channel_cache: str) -> tuple:
    """
    Convert list of sources to channels to process
    The sources list
       inputs:
           sources is a list of either filenames or regular expressions
           channel_cache is the filename where the channel id map is cached
       returns:
           fields dict in form {id:path}
           regexes list of regexes from sources that had matching channels
    """
    if not sources:
        return None, None

    oflds_f = get_chanids_from_files(sources)
    if oflds_f:
        logging.debug("found %d valid channels from input files", len(oflds_f))
        return oflds_f, None

    oflds_r, regexes = get_chanids_from_regexes(sources, channel_cache)
    if oflds_r:
        logging.debug("found valid channels from channel patterns")
        return oflds_r, regexes

    return None, None


def get_chanids_from_files(filelist: list) -> dict:
    """
    Convert list of sources to channels to process
    The sources list
       inputs:
           sources is a list of either filenames or regular expressions
       output:
           fields dict in form {id:path}
    """
    # loop over filelist to define the channels for query/output
    # using a file with 4 fields per line (after comment removal) where
    # the line format is '\s+id:\s+<chan_id>\s+path:\s+<path>$'
    # example: " id: 15176  path: focal-plane/R22/Reb0/Temp3"
    #
    # channels listed are used to construct a dict() to be returned
    # oflds[channel_id_num] = trending_full_path
    # eg: oflds[2372] = aliveness-fp/R00/Reb2/RGL
    #
    # from https://stackoverflow.com/questions/16710076
    # regex to split a string preserving quoted fields
    #
    rpat = re.compile(
        r"""
        (?:[^\s"']+)|    # match non-delimiter
        "(?:\\.|[^"]*)"| # match double quoted
        '(?:\\.|[^']*)'  # match single quoted
        """,
        re.X,
    )

    oflds = dict()
    for csource in filelist:
        logging.debug("channel_source= %s", csource)
        #  test to determine type of channel_source
        #
        logging.debug("test for formatted input file...")
        try:
            cf = open(csource, mode="r")
        except OSError as e:
            logging.debug("open(%s) failed: %s", csource, e)
            return None
        else:
            # populate oflds[id] with the corresponding channel path
            for line in cf:
                if re.match(r"^\s*#", line):  # skip block comment
                    continue
                if re.match(r"^\s*$", line):  # skip white space line
                    continue
                # strip inline cmnt
                sline = re.sub(r"""(#[^\'^"]*$)""", "", line)
                # tokenize what remains
                flds = ["".join(t) for t in rpat.findall(sline)]
                if len(flds) != 4:
                    logging.warning("bad input line: %s", line)
                    continue
                # example: " id: 15176  path: focal-plane/R22/Reb0/Temp3"
                if flds[1].isdigit():
                    oflds[flds[1]] = flds[3]  # oflds[id] = path
                else:
                    logging.warning("bad input line format: %s", line)
            cf.close()
    return oflds


def get_chanids_from_regexes(sources: list, channel_cache: str) -> tuple:
    """
    Convert list of sources to channels to process
    The sources list
       inputs:
           sources is a list of either filenames or regular expressions
           channel_cache is the filename where the channel id map is cached
       output:
           fields dict in form {id:path}
           regexes list of regexes from sources that had matching channels
    """
    logging.debug("dir()={}".format(dir()))
    global chid_dict
    oflds = dict()
    regexes = []
    for csource in sources:
        logging.debug("channel_source= %s", csource)
        #  test to determine type of channel_source
        #
        logging.debug("eval pattern for matching channels...")
        if not chid_dict:
            logging.debug("initializing channel id:path dictionary")
            chid_dict = get_channel_dict(channel_cache)

        # add csource as a new parameter in oflds[csource][chid]?
        # search the entire catalog for each pattern, this may be slow
        cpat = re.compile(csource)
        for chid in chid_dict:
            if cpat.search(chid_dict[chid]):
                oflds[chid] = chid_dict[chid]

        regexes.append(csource)

    return oflds, regexes


def init_trending_from_input_xml(input_files: list) -> str:
    """
    Initialize chid_dict and trending site using input data file
    inputs:  xml file(s) with trending data
    returns: site
    outcome: construct the module level dict() chid_dict
             mapping channel ids -> channel paths

    build a channel dictionary for all channels
    structure of input data xml file is:
    0: datas [-] (with attributes trending_server, trending_port, trending_tz
    1: data [id, path]
        2: trendingresult [-]
            3: channelmetadata [-]
                4: channelmetadatavalue [tstart, tstop, name, value]
            3: trendingdata [-]
                4: datavalue [name, value]
                4: axisvalue [name, value, loweredge, upperedge]

    where the top bit has attributes:
        <datas trending_server="localhost" trending_port="8080" trending_tz="PDT>"
    """
    global chid_dict
    global site
    logging.debug("initializing id:path dict from input data (xml) files")
    for ifile in input_files:
        logging.debug("ifile= %s", ifile)
        try:
            iff = open(ifile, mode="rb")
            xmlstring = iff.read()
            iff.close()
        except OSError as e:
            logging.debug("open(%s) failed: %s", ifile, e)
            sys.exit(1)
        else:
            root = etree.fromstring(xmlstring)
            # set the trending site to get the timezone etc
            tserver = root.attrib.get("trending_server")
            logging.debug("tserver = %s", tserver)
            res = None
            if not site:
                for ssite in sites:  # break if match
                    logging.debug("checking site = %s", ssite)
                    res = re.search(sites[ssite]["server"], tserver)
                    logging.debug(
                        "res = re.search(%s, %s)", sites[ssite]["server"], tserver
                    )
                    logging.debug("res = %s", res)
                    if res:
                        logging.debug("using site = %s", ssite)
                        site = ssite
                        break
            for data in root.iter("data"):
                chid = data.attrib.get("id")
                path = data.attrib.get("path")
                if chid and path:
                    chid_dict[chid] = path

    logging.debug("channel dict contains %d active channels", len(chid_dict))
    return sites[site]
