"""
Common image utility functions
"""
import re
import sys
import logging
import warnings
import datetime
import os.path
from astropy.io import fits
from astropy import stats
from astropy import wcs
from astropy.utils.exceptions import AstropyWarning
from astropy.convolution import convolve, Gaussian1DKernel, RickerWavelet1DKernel
import numpy as np
from scipy.ndimage import minimum_filter1d
from scipy.ndimage import percentile_filter
from scipy.interpolate import interp1d


def create_output_hdulist(hdulisti: fits.HDUList, argv: list) -> fits.HDUList:
    """
    Create output HDUList from input HDUList for building new image
    that is the result of processing the inputs (eg. not a blank).

    The Primary header of the input HDUList is used to create the
    Primary header of the output HDUList by appending to the bare
    output HDUList.

    DATE and an HISTORY header cards added to record what was done
    This is generally the first step before subsequent ops to modify
    data arrays and changing additional header keys.
    """
    logging.debug("creating output hdulist")
    # Create the output image, copy and update header comments, history
    hdulisto = fits.HDUList(fits.PrimaryHDU(None, hdulisti[0].header))
    hdu = hdulisto[0]
    hdr = hdu.header
    cstr = hdr.comments["DATE"]  # copy comment
    hdr.rename_keyword("DATE", "DATEORIG", force=True)
    hdr.comments["DATEORIG"] = "Previous file date/time"
    # FITS date format: 'yyyy-mm-ddTHH:MM:SS[.sss]'
    dtstr = datetime.datetime.utcnow().isoformat(timespec="milliseconds")
    hdr.insert("DATEORIG", ("DATE", dtstr, cstr))
    # add HISTORY lines
    hdr.add_history(
        "Header written by {} at: {}".format(os.path.basename(argv[0]), dtstr)
    )
    hdr.add_history("CMD: {} {}".format(os.path.basename(argv[0]), " ".join(argv[1:])))
    return hdulisto


def init_image_hdu(
    hdui: fits.ImageHDU, hdulisto: fits.HDUList, region: tuple = None
) -> fits.ImageHDU:
    """
    Append a new image HDU to output image using input HDU as a template.

    Copy the header and set the size/region specs in preparation for data
    to be added later.

    Returns
    -------
    hduo: fits.ImageHDU That was created during the call.

    """
    # create the output hdu from the master, (primary already exists)
    if not isinstance(hdui, fits.PrimaryHDU):
        hdri = hdui.header.copy()
        hdulisto.append(fits.ImageHDU(None, hdri, hdri["EXTNAME"]))
    hduo = hdulisto[len(hdulisto) - 1]
    hdro = hduo.header
    hdro["NAXIS"] = 2
    hdro.set("NAXIS1", hdri["NAXIS1"], "size of the n'th axis", after="NAXIS")
    hdro.set("NAXIS2", hdri["NAXIS2"], "size of the n'th axis", after="NAXIS1")
    hdro["BITPIX"] = -32
    # make changes to account for region of interest subimage
    if region and region is not (None, None):
        logging.debug("region = {}".format(region))
        naxis2 = (region[0].stop or len(hdui.data[:, 0])) - (region[0].start or 0)
        naxis1 = (region[1].stop or len(hdui.data[0, :])) - (region[1].start or 0)
        hdro.set("NAXIS1", naxis1, "size of the n'th axis", after="NAXIS")
        hdro.set("NAXIS2", naxis2, "size of the n'th axis", after="NAXIS1")
        # update any wcses
        wcses = wcs.find_all_wcs(hdro, fix=False)
        for w in wcses:
            wreg = w.slice(region)
            wreghdro = wreg.to_header()
            for card in wreghdro.cards:
                key = card.keyword
                value = card.value
                comment = card.comment
                hdro.set(key, value, comment)
    # logging.debug('output header:\n%s\n', hdro.tostring())
    return hduo


def parse_region(reg: str) -> tuple:
    """
    Return a pair of slices (slice1, slice2) corresponding
    to the region give as input in numpy slice string format
    If the region can't be parsed sys.exit() is called
    """
    try:
        slices = str_to_slices(reg)
    except ValueError as ve:
        logging.error("ValueError: %s", ve)
        logging.error("Bad region spec: %s", reg)
        sys.exit(1)

    if len(slices) != 2:
        logging.error("Bad region spec: %s", reg)
        sys.exit(1)

    return slices


def parse_iraf_region(reg: str) -> tuple:
    """
    Return a pair of slices (slice1, slice2) corresponding
    to the region give as input in ~IRAF format
    If the region can't be parsed (None, None) is returned
    """
    # peel off any outer brackets
    reg = re.sub(r"^\[([^\]]*)\]$", r"\1", reg)
    #
    # reg = [x1:x2,y1:y2] -- standard rectangle)
    if re.match(r"([0-9]*):([0-9]+),\s*([0-9]+):([0-9]+)$", reg):
        (x1, x2, y1, y2) = re.match(
            r"([0-9]+):([0-9]+),\s*([0-9]+):([0-9]+)$", reg
        ).groups()
        retval = (slice(int(y1) - 1, int(y2)), slice(int(x1) - 1, int(x2)))
    #
    # reg = [x0,y1:y2] -- single column section)
    elif re.match(r"([0-9]+),\s*([0-9]+):([0-9]+)$", reg):
        (x0, y1, y2) = re.match(r"([0-9]+),\s*([0-9]+):([0-9]+)$", reg).groups()
        retval = (slice(int(y1) - 1, int(y2)), slice(int(x0) - 1))
    #
    # reg = [*,y1:y2]) -- row selection
    elif re.match(r"(\*),\s*([0-9]+):([0-9]+)$", reg):
        (x, y1, y2) = re.match(r"(\*),\s*([0-9]+):([0-9]+)$", reg).groups()
        retval = (slice(int(y1) - 1, int(y2)), slice(None, None))
    #
    # reg = [x1:*,y1:y2]) -- row selection w/cols to end
    elif re.match(r"([0-9]+):\s*(\*),\s*([0-9]+):([0-9]+)$", reg):
        (x1, x2, y1, y2) = re.match(
            r"([0-9]+):\s*(\*),\s*([0-9]+):([0-9]+)$", reg
        ).groups()
        retval = (slice(int(y1) - 1, int(y2)), slice(int(x1) - 1, None))
    #
    # reg = [*:x1,y1:y2]) -- row selection w/cols from beginning
    elif re.match(r"(\*):\s*([0-9]+),\s*([0-9]+):([0-9]+)$", reg):
        (x1, x2, y1, y2) = re.match(
            r"(\*):\s*([0-9]+),\s*([0-9]+):([0-9]+)$", reg
        ).groups()
        retval = (slice(int(y1) - 1, int(y2)), slice(None, int(x2) - 1))
    #
    # reg = [x0,y0] -- single pixel
    elif re.match(r"([0-9]+),\s*([0-9]+)$", reg):
        (x0, y0) = re.match(r"([0-9]+),\s*([0-9]+)$", reg).groups()
        retval = (slice(int(y0)), slice(int(x0)))
    #
    # reg = [x1:x2,y0] -- single row section
    elif re.match(r"([0-9]+):([0-9]+),\s*([0-9]+)$", reg):
        (x1, x2, y0) = re.match(r"([0-9]+):([0-9]+),\s*([0-9]+)$", reg).groups()
        retval = (slice(int(y0) - 1), slice(int(x1) - 1, int(x2)))
    #
    # reg = [x1:x2,*] -- column selection
    elif re.match(r"([0-9]+):([0-9]+),\s*(\*)$", reg):
        (x1, x2, y) = re.match(r"([0-9]+):([0-9]+),\s*(\*)$", reg).groups()
        retval = (slice(None, None), slice(int(x1) - 1, int(x2)))
    #
    # reg = [*,*] # redundant, for completeness)
    elif re.match(r"(\*),\s*(\*)$", reg):
        (x, y) = re.match(r"(\*),\s*(\*)$", reg).groups()
        retval = (slice(None, None), slice(None, None))
    #
    # no match found, bad spec
    else:
        logging.error("bad region spec: '%s' no match produced", reg)
        retval = (None, None)
    #
    return retval


def get_requested_hduids(
    hdulist: fits.HDUList, hdunames: list, hduindices: list
) -> list:
    """
    Return a list of image hduids requested in optlist or all by default.

    Check that they exist in hdulist.  Requested hduids that
    don't exist are skipped.  Redundant values are dropped.
    """
    logging.debug("get_requested_hduids() called")
    hduids = []  # list of candidate hduids
    for name in hdunames or []:
        for hdu in hdulist:
            if re.search(name, hdu.name):
                try:
                    hduid = hdulist.index_of(hdu.name)
                    if hduid not in hduids:
                        hduids.append(hduid)
                except KeyError as ke:
                    logging.error("KeyError: %s", ke)
                    logging.error("HDU[%s] not found, skipping", hdu.name)
    for hduid in hduindices or []:
        try:
            hdu = hdulist[hduid]
            if hduid not in hduids:
                hduids.append(hduid)
        except IndexError:
            logging.error("HDU[%d] not found, skipping", hduid)
    if not hduindices and not hdunames:
        for hdu in hdulist:
            hduids.append(hdulist.index(hdu))

    if hduids:
        return hduids
    return None


def get_requested_image_hduids(
    hdulist: fits.HDUList, hdunames: list, hduindices: list
) -> list:
    """
    Return a list of image hduids requested in optlist or all by default.

    Check that they exist in hdulist and have data.  Requested hduids that
    don't exist are skipped.  Redundant values are dropped.
    """
    logging.debug("get_requested_hduids() called")
    chduids = []  # list of candidate hduids
    for name in hdunames or []:
        for hdu in hdulist:
            if re.search(name, hdu.name):
                try:
                    hduid = hdulist.index_of(hdu.name)
                    if hduid not in chduids:
                        chduids.append(hduid)
                except KeyError as ke:
                    logging.error("KeyError: %s", ke)
                    logging.error("HDU[%s] not found, skipping", hdu.name)
    for hduid in hduindices or []:
        try:
            hdu = hdulist[hduid]
            if hduid not in chduids:
                chduids.append(hduid)
        except IndexError:
            logging.error("HDU[%d] not found, skipping", hduid)
    if not hduindices and not hdunames:
        for hdu in hdulist:
            chduids.append(hdulist.index(hdu))

    # Validate the list of candidate HDUs, keep those with pixels
    hduids = []
    for hduid in chduids:
        hdu = hdulist[hduid]
        if isinstance(hdu, fits.PrimaryHDU):  # check for data
            hdr = hdu.header
            if hdr.get("NAXIS") == 2:
                if hdr.get("NAXIS1") and hdr.get("NAXIS2"):
                    naxis1 = hdr.get("NAXIS1")
                    naxis2 = hdr.get("NAXIS2")
                    if naxis1 * naxis2 > 0:
                        logging.debug(
                            "adding %s with index %d to hduid list", hdu.name, hduid
                        )
                        hduids.append(hduid)
        elif isinstance(hdu, (fits.ImageHDU, fits.CompImageHDU)):
            logging.debug("adding %s with index %d to hduid list", hdu.name, hduid)
            hduids.append(hduid)
        else:
            logging.debug(
                "%s with index %d is not type (Comp)ImageHDU", hdu.name, hduid
            )
    if hduids:
        return hduids
    return None


def get_data_oscan_slices(hdu: fits.FitsHDU) -> tuple:
    """
    Get datasec, serial/parallel overscan as slice specifications.

    Given an hdu, uses header keys to infer slice specs.  If a particular
    region cannot be obtained a spec of (None, None) is returned for that
    region.

    Returns a tuple of slice definitions (datasec, soscan, poscan).
    The serial overscan is assumed to be at the end of each row if present.
    """
    # first get serial and parallel overscan region defs
    hdr = hdu.header
    try:
        dstr = hdr["DATASEC"]
    except KeyError as ke:
        logging.debug("KeyError: %s required", ke)
        return (None, None, None)
    logging.debug("DATASEC=%s", dstr)
    try:
        n1 = hdr["NAXIS1"]
    except KeyError as ke:
        logging.error("KeyError: %s required", ke)
        return (None, None, None)
    try:
        n2 = hdr["NAXIS2"]
    except KeyError as ke:
        logging.error("KeyError: %s required", ke)
        return (None, None, None)
    # get DATASEC region
    datasec = parse_iraf_region(dstr)
    if datasec == (None, None):
        return (None, None, None)
    (p1, p2) = (datasec[0].start or 0, datasec[0].stop or len(hdu.data[:, 0]))
    (s1, s2) = (datasec[1].start or 0, datasec[1].stop or len(hdu.data[0, :]))
    if n1 > s2:
        soscan = (slice(0, n2), slice(s2, n1))
    else:  # no serial overscan
        soscan = (slice(None), slice(None))
    if n2 > p2:
        poscan = (slice(p2, n2), slice(0, n1))
    else:
        poscan = (slice(None), slice(None))

    return (datasec, soscan, poscan)


def str_to_slices(sliceStr: str) -> tuple:
    """
    Parse a string containing one or more slice specs separated by commas

    Returns a tuple of slice() objects

    rewrite of:
    https://stackoverflow.com/questions/43089907/
    using-a-string-to-define-numpy-array-slice
    to make it straightforward albeit not nearly as elegant
    """
    # peel off any outer brackets
    sliceStr = re.sub(r"^\[([^\]]*)\]$", r"\1", sliceStr)

    slices = []
    for sspec in sliceStr.split(","):
        if ":" not in sspec:
            slice_args = [int(sspec), int(sspec) + 1]
            slices.append(slice(*tuple(slice_args)))
        else:
            slice_args = []
            for item in sspec.strip().split(":"):
                if item:
                    slice_args.append(int(item))
                else:
                    slice_args.append(None)
            slices.append(slice(*tuple(slice_args)))

    return tuple(slices)


def subtract_bias(stype: str, ptype: str, hdu: fits.ImageHDU):
    """
    Subtract a bias estimate (using overscans) from an hdu.
    Operates in-place on the Image.HDU parameter

    Choices are 'None', 'mean' 'median', 'by(row|col)', 'by(row|col)smooth' and
    applying only to parallel on lsste2v's: 'lsste2v'

    Bias estimates are calculated using DATASEC to infer the overscan regions.

    The fits.ImageHDU is operated on directly
    """
    (datasec, soscan, poscan) = get_data_oscan_slices(hdu)
    logging.debug("bias stype=%s ptype=%s", stype, ptype)
    pcnt = 30.0  # percentile for signal est
    max_rn = 7.0
    rn_est = min(np.std(hdu.data[poscan[0], soscan[1]]), max_rn)

    # First pass uses serial overscan region
    if stype:
        if stype in ("byrow", "byrowsmooth"):
            so_med = np.percentile(hdu.data[soscan], 50, axis=1)
            # clean up any crazy rows (eg overflow in serial from hot column)
            so_med_avg, so_med_med, so_med_std = stats.sigma_clipped_stats(so_med)
            logging.debug(
                "so_med_avg=%.3f, so_med_med=%.3f, so_med_std=%.3f",
                so_med_avg,
                so_med_med,
                so_med_std,
            )
            so_med_bad_ind = np.nonzero(np.abs(so_med - so_med_med) > 100 * rn_est)
            logging.debug("so_med_bad_ind=%s", so_med_bad_ind)
            logging.debug(
                "replacing %d values in so_med[] with %.3f",
                len(so_med_bad_ind),
                so_med_med,
            )
            so_med[so_med_bad_ind] = so_med_med
            if stype in ("byrowsmooth",):
                logging.debug("smoothing serial overscan with Gaussian1DKernel")
                kernel = Gaussian1DKernel(1)
                so_med = convolve(so_med, kernel, boundary="extend")
            # convert shape from (n,) to (n, 1)
            logging.debug("mean serial overscan subtraction: %d", np.median(so_med))
            logging.debug("first 20 rows: \n%s", so_med[0:20])
            so_med = so_med.reshape(np.shape(so_med)[0], 1)
            hdu.data = hdu.data - so_med
        elif stype == "mean":
            hdu.data = hdu.data - np.mean(hdu.data[soscan])
        elif stype == "median":
            hdu.data = hdu.data - np.median(hdu.data[soscan])
        else:
            logging.error("stype: %s not valid", stype)
            sys.exit(1)

    # Second pass uses parallel overscan
    if ptype:
        if ptype in ("bycol", "bycolfilter", "bycolsmooth", "lsste2v"):
            if ptype in ("bycol"):
                bias_row = np.percentile(hdu.data[poscan[0], :], pcnt, axis=0)
            if ptype in ("bycolfilter", "bycolsmooth", "lsste2v"):
                bias_row = get_bias_filtered_est_row(hdu)
                if bias_row is None:
                    logging.warning("could not perform parallel bias subtraction")
                    return
                if ptype in ("bycolsmooth"):
                    logging.debug("smoothing par overscan with Gaussian1DKernel")
                    kernel = Gaussian1DKernel(2)
                    # don't smooth the prescan
                    bias_row[datasec[1].start :] = convolve(
                        bias_row[datasec[1].start :], kernel, boundary="extend"
                    )
            # convert shape from (,n) to (1, n)
            bias_row = bias_row.reshape(1, np.shape(bias_row)[0])
            hdu.data = hdu.data - bias_row.data
            logging.debug("bias_row_median = %.2f", np.median(bias_row.data))
        elif ptype == "mean":
            hdu.data = hdu.data - np.mean(hdu.data[poscan])
        elif ptype == "median":
            hdu.data = hdu.data - np.median(hdu.data[poscan])
        else:
            logging.error("ptype: %s not valid", ptype)
            sys.exit(1)

    # third pass to take out special bias effect on selected e2v CCDs
    if ptype and ptype in ("lsste2v"):
        # subtract an exp decay with amplitude from prescan along each row
        a0 = np.mean(hdu.data[:, 1 : datasec[1].start], axis=1)
        naxis1 = np.shape(hdu.data)[1]
        i0 = (-3.0 / naxis1) * np.arange(naxis1)  # exp decay row vector
        e0 = np.exp(i0)
        # subtract exp decay function along each row
        for i in np.arange(np.size(hdu.data[:, 0])):
            hdu.data[i, :] -= a0[i] * e0
        logging.debug("third_pass_median = %.2f", np.median(a0) * np.median(e0))


def eper_serial(hdu):
    """
    Given datasec and serial overscan as slices, calculate
    eper using the first ecols=3 columns of serial overscan
    """
    datasec, soscan, poscan = get_data_oscan_slices(hdu)
    ecols = 3  # number of columns used for eper signal
    pcnt = 30.0  # percentile for signal est
    ncols = datasec[1].stop - datasec[1].start
    scols = int(0.20 * ncols)

    # signal estimate 1-d array (30% is ~sky)
    sig_est_col = np.percentile(
        hdu.data[datasec[0], (datasec[1].stop - scols) : datasec[1].stop], pcnt, axis=1
    )
    # deferred charge estimate (before bias subtraction)
    dc_sum_col = np.sum(
        hdu.data[datasec[0], soscan[1].start : (soscan[1].start + ecols)], axis=1
    )
    bias_est_col = np.median(hdu.data[datasec[0], (soscan[1].start + ecols) :], axis=1)

    sig_est_col = sig_est_col - bias_est_col
    dc_est_col = dc_sum_col - ecols * bias_est_col

    dc_avg, dc_med, dc_std = stats.sigma_clipped_stats(dc_est_col)
    sig_avg, sig_med, sig_std = stats.sigma_clipped_stats(sig_est_col)

    # if dc_avg > 0 and sig_avg > 0:
    cti_est = dc_avg / sig_avg / ncols

    if cti_est > -0.0001:
        eper = 1 - cti_est
        return eper
    else:
        logging.debug("s-cti est was < 0")
        return None


def get_bias_filtered_est_row(hdu):
    """
    Given hdu, produce a suitable parallel bias estimate for subtraction
    An effort is made to deal with saturation until it gets too high.
    The filtered row attempts to interpolate across regions with bad/hot columns
    """
    datasec, soscan, poscan = get_data_oscan_slices(hdu)
    # parameters
    pcnt = 30.0  # percentile to use for bias estimate
    max_rn = 7.0
    srows = int(0.05 * (datasec[0].stop - datasec[0].start))
    pstart = poscan[0].start
    pstop = poscan[0].stop
    erows = 6
    if pstop - pstart < 2 * erows:
        logging.warning("Not enough parallel overscan rows to estimate pbias row")
        return None
    bias_est_row = np.percentile(hdu.data[poscan[0].start + erows :, :], pcnt, axis=0)
    bias_base = np.percentile(hdu.data[poscan[0], datasec[1]], pcnt)
    bias_floor = np.percentile(hdu.data[poscan[0], soscan[1]], pcnt)
    rn_est = min(np.std(hdu.data[poscan[0], soscan[1]]), max_rn)
    sig_est = np.median(hdu.data[datasec[0].stop - srows : datasec[0].stop, datasec[1]])
    logging.debug("sig_est = %.2f", sig_est)
    logging.debug("bias_floor = %.2f", bias_floor)
    logging.debug("bias_base = %.2f", bias_base)
    logging.debug("rn_est = %.2f", rn_est)
    # bad columns are ones that are above 25% of signal after erows into overscan
    # or if the signal is very small they have improbably high value
    bad_ind = np.array(
        np.nonzero(
            np.median(
                hdu.data[pstart + erows : pstart + 2 * erows, datasec[1]], axis=0,
            )
            > bias_floor + sig_est / 4.0 + 50 * rn_est
        )
    )
    # add +/-3 adjacent columns to the bad list
    bad_ind = np.union1d(np.union1d(bad_ind - 1, bad_ind), bad_ind + 1)
    bad_ind = np.union1d(np.union1d(bad_ind - 1, bad_ind), bad_ind + 1)
    bad_ind = np.union1d(np.union1d(bad_ind - 1, bad_ind), bad_ind + 1)
    bad_ind = np.intersect1d(np.arange(datasec[1].stop - datasec[1].start), bad_ind)
    if np.size(bad_ind) > 0:
        logging.debug(
            "possible hot cols: %s", bad_ind + datasec[1].start,
        )
    if np.size(bad_ind) > 0.5 * (datasec[1].stop - datasec[1].start):
        return None
    # good columns are the complementary set
    good_ind = np.setdiff1d(np.arange(datasec[1].stop - datasec[1].start), bad_ind)
    good_ind = good_ind + datasec[1].start
    # mask off and interpolate across rejected columns, leave pre/overscan alone
    x = np.arange(datasec[1].start, datasec[1].stop)  # x coords
    interp_function = interp1d(
        good_ind, bias_est_row[good_ind], bounds_error=False, fill_value="extrapolate",
    )
    bias_est_row[datasec[1].start : datasec[1].stop] = interp_function(x)
    bias_est_row = bias_est_row - bias_floor  # ensure dbl oscan no change
    logging.debug("bias_est_row = %s", bias_est_row[bad_ind])
    return bias_est_row


def eper_parallel(hdu):
    """
    Given hdu, calculate eper using the first erows (6) rows of parallel overscan
    Note once eper <~ 0.998 it is not accurate as deferred charge then extends
    out beyond the 6 rows used
    """
    datasec, soscan, poscan = get_data_oscan_slices(hdu)
    # need a return None if any of those are missing

    erows = 8  # number of rows used to measure deferred charge
    nrows = datasec[0].stop - datasec[0].start
    srows = int(0.05 * nrows)
    pstart = poscan[0].start
    pstop = poscan[0].stop
    prows = pstop - pstart

    if prows < 2 * erows:
        logging.warning("parallel overscan too small to estimate cte")
        return None

    # bias floor and read noise estimate using double overscan region
    bias_floor = np.percentile(hdu.data[poscan[0], soscan[1]], 30)
    logging.debug("bias_floor = %.2f", bias_floor)
    read_noise_est = min(np.std(hdu.data[poscan[0], soscan[1]]), 7.0)
    logging.debug("read_noise_est = %.2f", read_noise_est)

    # signal estimate 1-d array (use last 5% of rows), not bias subtracted
    sig_est_row = np.median(
        hdu.data[datasec[0].stop - srows : datasec[0].stop, datasec[1]], axis=0
    )
    sig_est0 = np.percentile(sig_est_row, 20) - bias_floor  # estimate
    logging.debug("sig_est0 = %.2f", sig_est0)
    #  get column indices to use in determining p-cti
    if sig_est0 > int(1 << 14) * read_noise_est:  # assuming ~16k dynamic range
        logging.debug("using high signal case")
        dc_est_row = (
            np.sum(hdu.data[poscan[0], datasec[1]], axis=0) - bias_floor * prows
        )
        sig_est_row = sig_est_row - bias_floor
        bad_ind = np.array(np.nonzero(dc_est_row > sig_est0 * (prows - erows)))
    else:
        logging.debug("using NOT high signal case")
        # bias estimate
        bias_est_row = np.percentile(
            hdu.data[pstart + erows :, datasec[1]], 30.0, axis=0
        )
        # deferred charge estimate, use complement of bias_est_rows
        dc_est_row = (
            np.sum(hdu.data[pstart : pstart + erows, datasec[1]], axis=0)
            - bias_est_row * erows
        )
        sig_est_row = sig_est_row - bias_est_row

        # estimate bias level by filtering out high end of distribution
        bias_base = np.percentile(bias_est_row, 20)
        logging.debug("bias_base = %.2f", bias_base)
        # bad columns are ones that are still 25% of signal after erows into overscan
        # or if the signal is very small they have improbably high value
        bad_ind = np.array(
            np.nonzero(
                np.mean(
                    hdu.data[pstart + erows : pstart + 2 * erows, datasec[1],], axis=0,
                )
                > bias_base + sig_est0 / 4.0 + 50 * read_noise_est
            )
        )
        logging.debug(
            "erows * np.sum(bias_est_row - bias_floor) = %.2f",
            erows * (np.sum(bias_est_row - bias_floor)),
        )
    # add +/-3 adjacent columns to the bad list
    bad_ind = np.union1d(np.union1d(bad_ind - 1, bad_ind), bad_ind + 1)
    bad_ind = np.union1d(np.union1d(bad_ind - 1, bad_ind), bad_ind + 1)
    bad_ind = np.union1d(np.union1d(bad_ind - 1, bad_ind), bad_ind + 1)
    bad_ind = np.intersect1d(np.arange(datasec[1].stop - datasec[1].start), bad_ind)
    # good columns are the complementary set
    good_ind = np.setdiff1d(np.arange(datasec[1].stop - datasec[1].start), bad_ind)

    dc_est = np.sum(dc_est_row[good_ind])
    sig_est = np.sum(sig_est_row[good_ind])
    logging.debug("%d cols had usable signal in eper_parallel", np.size(good_ind))
    if np.size(bad_ind) > 0:
        logging.debug(
            "possible hot cols: %s",
            np.intersect1d(np.arange(datasec[1].stop - datasec[1].start), bad_ind)
            + datasec[1].start,
        )
    logging.debug("bias_floor = %.2f", bias_floor)
    logging.debug("read_noise_est = %.2f", read_noise_est)
    logging.debug("dc_est = %.2f  sig_est = %.2f  nrows = %d", dc_est, sig_est, nrows)
    if np.size(good_ind) < 0.5 * (datasec[1].stop - datasec[1].start):
        logging.debug("not enough good columns to determine p-cte")
        return None

    cti_est = dc_est / sig_est / nrows

    if cti_est > -0.0001:
        eper = 1 - cti_est
        return eper
    else:
        logging.warning("p-cti est was < 0")
        return None


def files_to_hdulists(ifiles: list, mmp: bool = True) -> list:
    """
    Given a list of image files return a list of fits.HDUList objects
    that are verified as commensurate for processing as a set (combining etc.)

    The mmp input flag defaults to True to enable memory mapping being used.
    If there are many large files then calling with mmp = False and
    processing by sectioning is another choice.
    """
    # set up the items used to verify file match each other
    #
    list_of_hdulists = []
    for cnt, ffile in enumerate(ifiles):
        try:
            hdulist = fits.open(ffile, memmap=mmp)
        except IOError as ioerr:
            logging.error("IOError: %s", ioerr)
            sys.exit(1)

        # compare selected parameters per hdu per file
        hdu_pars = []  # list of dict()s
        for hdu in hdulist:
            hdr = hdu.header
            hdudict = dict()
            # determine type of HDU
            if isinstance(hdu, fits.PrimaryHDU):  # check for data
                hdudict["type"] = "PrimaryHDU"
            elif isinstance(hdu, (fits.ImageHDU, fits.CompImageHDU)):
                hdudict["type"] = "ImageHDU"
            else:
                hdudict["type"] = "other"
            # get pixel data info for hdu if exists
            hdudict["dimension"] = (None, None)
            if hdudict["type"] in ("ImageHDU", "PrimaryHDU"):
                if hdr.get("NAXIS") == 2:
                    if hdr.get("NAXIS1") and hdr.get("NAXIS2"):
                        naxis1 = hdr.get("NAXIS1")
                        naxis2 = hdr.get("NAXIS2")
                        if naxis1 * naxis2 > 0:
                            hdudict["dimension"] = (naxis1, naxis2)
            hdu_pars.append(hdudict)
            # end of loop overy hdus within file
        if cnt == 0:  # first file defines the valid parameters
            base_pars = hdu_pars
        else:  # compare hdu_pars to first file
            for hpar, bpar in zip(hdu_pars, base_pars):
                for key in bpar.keys():
                    if hpar[key] != bpar[key]:
                        logging.error(
                            "file parameter mismatch: %s: %s != %s",
                            key,
                            hpar[key],
                            bpar[key],
                        )
                        sys.exit(1)
        # end of loop over files
        list_of_hdulists.append(hdulist)
    return list_of_hdulists


def image_combine_hdu(
    iimages: list,
    hduid: int,
    method: list,
    region: tuple,
    bimage: fits.HDUList,
    stype: str,
    ptype: str,
    scaling: tuple,
    hduo: fits.ImageHDU,
):
    """
    From a list of input images (as hdulists) and the id of one extension
    return an ImageHDU.data object containing a pixel-by-pixel "average" of
    the stacked input images.  The processing varies according to the
    additional arguments as to median vs. average, bias subtraction etc.

    Parameters
    ----------
    iimages: list of astropy.io.fits.HDUList objects
    hduid: index specifying a single hdu (present in all iimages) to process
    method: [median], [average], [sigmaclip, sigmaval], [rank, percentile]
    region: (yslice, xslice) specifying ROI to process, full image if None
    bimage: fits.HDUList object with (bias) image to subtract
    stype: param for subtract_bias() function (in this module)
    ptype: param for subtract_bias() function (in this module)
    scaling: (yslice, xslice) specifying ROI to use for scaling
    hduo: a basic ImageHDU object that is modified and is the functions result
    """
    hdudata_list = []
    hdu_scale = []
    for im in iimages:
        hdu = im[hduid].copy()
        if stype | ptype:
            subtract_bias(stype, ptype, hdu)
        if scaling:
            svalue = np.median(hdu.data[scaling[0], scaling[1]])
            hdu_scale.append(svalue)
        if region:
            hdudata_list.append(hdu.data[region[0], region[1]])
            if bimage:
                bdata = bimage[hduid].data[region[0], region[1]]
        else:
            hdudata_list.append(hdu.data)
            if bimage:
                bdata = bimage[hduid].data

    if scaling:  # pass through data and scale it
        hdu_scale_arr = np.asarray(hdu_scale)
        # normalize the scale factors
        hdu_scale_arr = np.mean(hdu_scale_arr) / hdu_scale_arr
        logging.debug(f"applying scale factors: {hdu_scale_arr}")
        for hdudata, hduscale in zip(hdudata_list, hdu_scale_arr):
            hdudata = hdudata * hduscale

    logging.debug(f"using method: {method}")
    if re.match(r"^mea", method[0]):
        hduo.data = np.mean(np.array(hdudata_list), axis=0)
    elif re.match(r"^med", method[0]):
        hduo.data = np.median(np.array(hdudata_list), axis=0)
    elif re.match(r"^sig", method[0]):  # this one is ugly
        hduo.data = np.nanmean(
            stats.sigma_clip(np.array(hdudata_list), method[1], axis=0, masked=False),
            axis=0,
        )
    elif re.match(r"^ran", method[0]):
        hduo.data = np.percentile(np.array(hdudata_list), method[1], axis=0)
    else:
        logging.error("image combine method %s not recognized", method[0])
        sys.exit(1)

    if bimage:
        hduo.data = hduo.data - bdata
