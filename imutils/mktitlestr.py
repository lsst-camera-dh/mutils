#!/usr/bin/env python
"""
Produce a glob-like pattern from a list of files to indicate
the file set used (nominally to produce a title string)
"""
import os
import sys
import re
import argparse
import logging
import textwrap

# put parent directory into sys.path
bp = os.path.dirname(os.path.realpath(__file__)).split(os.sep)
modpath = os.sep.join(bp[:-1] + ['lib'])
sys.path.insert(0, modpath)

#  local imports
try:
    import trendutils as tu
    import mutils as mu
    import imutils as iu
except ImportError as e:
    logging.error('Import failed: %s', e)
    exit(1)


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        Generate a glob-like pattern from a list of strings.
                                    '''),
        epilog=textwrap.dedent('''\
        The intended use case is a list of filenames (or paths)
        given on the command line via unix globbing.  The returned
        pattern should be a reasonable guess at an equivalent glob
        pattern
                                '''))
    parser.add_argument("file", nargs='+',
                        metavar="file", help="input fits file(s)")
    parser.add_argument("--minsize", default=3, type=int,
                        help="min length substring, default: $(default)s")
    parser.add_argument("--debug", action='store_true',
                        help="print additional debugging messages")
    return parser.parse_args()


def main():
    """
    main logic
    a glob pattern for a given input list of files/paths
    """
    optlist = parse_args()
    mu.init_logging(optlist.debug)
    mu.init_warnings()
    filenames = []
    for fname in optlist.file:
        # fname = re.sub(r"^.*/(.*)$", r"\1", fname)
        fname = re.sub(r"^(.*)\.fits?(\.fz)*$", r"\1", fname)
        filenames.append(fname)
    logging.debug('using %d filenames', len(filenames))
    if len(filenames) < 1:
        logging.error('No filenames to process')
        exit(1)
    if len(filenames) == 1:
        print("{}".format(filenames[0]))
        exit(0)

    ss_arr = []
    iu.get_lcs_array(filenames, ss_arr, 0, '', optlist.minsize)
    if ss_arr:
        logging.debug('%d substrings were found', len(ss_arr))
        suptitle = "{}".format('*'.join(ss_arr))
        if not re.match(ss_arr[0], filenames[0]):
            suptitle = "*{}".format(suptitle)
        if not re.search(r"{}$".format(ss_arr[-0]), filenames[0]):
            suptitle = "{}*".format(suptitle)
        print("{}".format(suptitle))
    else:
        logging.debug("No common substrings found")


if __name__ == '__main__':
    main()
