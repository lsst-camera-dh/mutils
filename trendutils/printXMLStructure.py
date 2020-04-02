#!/usr/bin/env python
"""code to illustrate parsing/structure of XML files
"""
import os
import argparse
import logging
import collections
import re
import sys

# put ../lib directory into sys.path
bp = os.path.dirname(os.path.realpath(__file__)).split(os.sep)
modpath = os.sep.join(bp[:-1] + ['lib'])
sys.path.insert(0, modpath)

#  local imports
try:
    import trendutils as tu
except ImportError as e:
    logging.error('Import failed: %s', e)
    exit(1)


if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        description="Parse and display structure of an xml file")
    parser.add_argument("target", metavar="URL|file",
                        help="filepath or URL")
    parser.add_argument("--debug", action='store_true',
                        help="print additional debugging info")
    parser.add_argument("--structure", action='store_true',
                        help="print structure of xml")
    return parser.parse_args()


def main():
    """main logic"""
    optlist = parse_args()
    if optlist.debug:
        logging.basicConfig(format='%(levelname)s: %(message)s',
                            level=logging.DEBUG)
        logging.debug("printing extra debugging...")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s',
                            level=logging.INFO)

    xmlstr = tu.geturi(optlist.target)
    if xmlstr is None:
        logging.error('Failed to open uri: %s', optlist.target)
        exit(1)
    tu.print_channel_structure(xmlstr)


if __name__ ==  '__main__':
    main()
