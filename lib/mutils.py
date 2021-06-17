"""
Common utility functions that are generic (not image or trending)
"""
import re
import logging
import warnings
from astropy.utils.exceptions import AstropyWarning


def init_logging(debug):
    """Set up debug and info level logging"""
    if debug:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
    else:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
    # suppress plotting debug messages
    mpl_logger = logging.getLogger("matplotlib")
    mpl_logger.setLevel(logging.WARNING)


def init_warnings():
    """Block warnings from Astropy"""
    warnings.simplefilter("ignore", category=AstropyWarning)


def long_substr(ilist: list) -> str:
    """
    Find the longest common substring (sstr) from a list of strings.

    This is from:
    https://stackoverflow.com/questions/2892931/\
        longest-common-substring-from-more-than-two-strings-python#
    """
    sstr = ""
    if len(ilist) > 1 and len(ilist[0]) > 0:
        for i in range(len(ilist[0])):
            for j in range(len(ilist[0]) - i + 1):
                if j > len(sstr) and all(ilist[0][i : i + j] in x for x in ilist):
                    sstr = ilist[0][i : i + j]
    return sstr


def cleave(arr: list, substr: str) -> (list, list):
    """
    Split each element in array of strings using a substring.

    Returns
    -------
    pre_arr, post_arr :
        Lists of strings with preceeding/following segments after splitting.
    """

    cleav_pat = re.compile("^(.*){}(.*)$".format(substr))
    pre_arr = []
    post_arr = []
    for istr in arr:
        tl, tr = cleav_pat.match(istr).groups()
        if tl:
            pre_arr.append(tl)
        if tr:
            post_arr.append(tr)

    return (pre_arr, post_arr)


def get_lcs_array(
    s_arr: list, ss_arr: list, index: int, order: str, minsz: int
) -> list:
    """
    Create ordered list of common substrings for list of strings.

    Given a list of filenames or paths or similar, find the set of all
    common substrings longer than 'minsz'.  The result can be used to create
    a 'glob-like' pattern to summarize the original list of strings.

    Parameters
    ----------
    s_arr: list Array of strings to process to find lcs.

    ss_arr: list Substring array -- this is the product.

    index: int Index in ss_arr of parent lcs.

    order: str Insert new lcs before|after index.

    minsz: int Shortest allowed common substring.
    """

    logging.debug("get_lcs_array() on entry: s_arr=%s", s_arr)
    logging.debug("get_lcs_array() on entry: index=%d", index)
    logging.debug("get_lcs_array() on entry: order=%s", order)
    s_arr_length = len(s_arr)
    lcstr = long_substr(s_arr)  # get longest common substring
    logging.debug("lcstr=%s", lcstr)
    if len(lcstr) >= minsz:
        if not ss_arr:  # first pass when ss_arr is None
            ss_arr.append(lcstr)
            logging.debug("first: ss_arr=%s", ss_arr)
        else:
            if order == "before":
                logging.debug("before: on entry ss_arr=%s", ss_arr)
                logging.debug("before: inserting %s at index=%d", lcstr, index)
                ss_arr.insert(index, lcstr)
                logging.debug("before: post insertion ss_arr=%s", ss_arr)
            elif order == "after":
                # treat index as offset from end
                logging.debug("after: on entry ss_arr=%s", ss_arr)
                if index == 0:
                    logging.debug("after: appending %s", lcstr)
                    ss_arr.append(lcstr)
                else:
                    logging.debug("after: inserting %s at index=%d", lcstr, -1 * index)
                    ss_arr.insert(-1 * index, lcstr)
                logging.debug("after: final ss_arr=%s", ss_arr)
            else:
                logging.error("invalid order %s", order)
                exit(1)
    else:
        return

    pre_arr, post_arr = cleave(s_arr, lcstr)
    logging.debug("pre_arr=%s", pre_arr)
    logging.debug("post_arr=%s", post_arr)
    if pre_arr and len(pre_arr) == s_arr_length:
        if order == "after":
            pre_index = len(ss_arr) - index - 1
        else:
            pre_index = index
        logging.debug(
            "calling get_lcs_array(pre_arr, ss_arr, index=%d, before)", pre_index
        )
        get_lcs_array(pre_arr, ss_arr, pre_index, "before", minsz)
    if post_arr and len(post_arr) == s_arr_length:
        if order == "before":
            post_index = len(ss_arr) - index - 1
        else:
            post_index = index
        logging.debug(
            "calling get_lcs_array(post_arr, ss_arr, index=%d, after)", post_index
        )
        get_lcs_array(post_arr, ss_arr, post_index, "after", minsz)


def mkglob(fullpaths: list, trim=False) -> str:
    """
    For an input list of strings or filenames or full paths prepare
    a glob-like pattern that summarizes the list.  The nominal use case
    is to take a list of filenames and produce a glob pattern that
    could be used as a title in a plot or other graphic.  If trim=True
    the leading path and trailing extensions are removed.
    """
    string_list = []
    glob = None
    for fname in fullpaths:
        if trim:
            fname = re.sub(r"^.*/(.*)$", r"\1", fname)
            # fname = re.sub(r"^(.*)\.fits?(\.fz)*$", r"\1", fname)
            fname = re.sub(r"^([^\.]*)\..*$", r"\1", fname)  # trim suffix
        string_list.append(fname)
    logging.debug("string_list[]={}".format(string_list))
    if len(string_list) == 1:
        glob = string_list[0]
    elif len(string_list) > 1:
        # title is longest common substring array
        # joined with *'s to look like a glob pattern
        ss_arr = []
        get_lcs_array(string_list, ss_arr, 0, "", 2)
        if ss_arr:
            glob = "{}".format("*".join(ss_arr))
            if not re.match(ss_arr[0], string_list[0]):
                glob = "*{}".format(glob)
            if not re.search(r"{}$".format(ss_arr[-1]), string_list[0]):
                glob = "{}*".format(glob)
    return glob


def tok_line(line: str) -> list:
    """ """
    #
    # from https://stackoverflow.com/questions/16710076
    # regex to split a string preserving quoted fields
    #
    rpat = re.compile(
        r"""             #
                      (?:[^\s"']+)|    # match non-delimiter
                      "(?:\\.|[^"]*)"| # match double quoted
                      '(?:\\.|[^']*)'  # match single quoted
                      """,
        re.X,
    )
    if re.match(r"^\s*#", line):  # skip block comment
        return []
    if re.match(r"^\s*$", line):  # skip white space line
        return []
    # strip inline cmnt
    sline = re.sub(r"""(#[^\'^"]*$)""", "", line)
    # tokenize what remains
    return ["".join(t) for t in rpat.findall(sline)]


def file_to_tokens(ifiles: list) -> list:
    """
    Input is a list of files (or a single file)
    Each subsequent call tokenizes the next line of content
    skipping comment and blank lines.  Inline comments are filtered
    out as well (within reason).
    """
    #
    # derived from https://stackoverflow.com/questions/16710076
    # regex to split a string preserving quoted fields
    #
    rpat = re.compile(
        r"""             #
                      (?:[^\s"]+)|    # match non-delimiter
                      (?<=\W)"(?:\\.|[^"]*)"(?=\W)| # match double quoted
                      (?<=\W)'(?:\\.|.*?)'(?=\W)  # match single quoted
                      """,
        re.X,
    )
    #
    if not isinstance(ifiles, list):
        ifiles = [ifiles]
    #
    for ffile in ifiles:
        try:
            ff = open(ffile, mode="r")
        except OSError as e:
            logging.debug("open(%s) failed: %s", ffile, e)
        else:
            for line in ff:
                if re.match(r"^\s*#", line):  # skip block comment
                    continue
                if re.match(r"^\s*$", line):  # skip white space line
                    continue
                # strip inline cmnt '<space(s)># to end of line'
                sline = re.sub(r"""(\s*#[^\'"]*$)""", "", line)
                # tokenize what remains
                yield ["".join(t) for t in rpat.findall(sline)]
            ff.close()
