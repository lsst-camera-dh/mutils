#!/usr/bin/env python3

import os
import sys
import argparse
import textwrap
import logging

# put parent directory into sys.path
bp = os.path.dirname(os.path.realpath(__file__)).split(os.sep)
modpath = os.sep.join(bp[:-1] + ['lib'])
sys.path.insert(0, modpath)

#  local imports
try:
    import trendutils as tu
    import mutils as mu
except ImportError as e:
    logging.error('Import failed: %s', e)
    exit(1)


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        description="list Trending channels from CCS localdb")
    parser.add_argument("ssname", nargs='+',
                        help="CCS subsystem name(s)")
    parser.add_argument("--maxidle", metavar='days', type=int,
                        nargs='?', default=None, const=-1,
                        help="Include only if published within N=maxidle days")
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

    if optlist.maxidle:
        channel_file = tu.update_trending_channels_xml(
            86400 * optlist.maxidle)
    else:
        channel_file = tu.update_trending_channels_xml()
    chfile = open(channel_file, mode='rb')
    logging.debug('got channel_file')
    channels_xml = chfile.read()
    chfile.close()
    if optlist.structure:
        tu.print_channel_structure(channels_xml)
        return
    tu.print_channel_content(channels_xml, optlist.ssname)


if __name__ == '__main__':
    main()
