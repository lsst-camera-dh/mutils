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
from astropy.convolution import convolve, Gaussian1DKernel
import numpy as np


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
    to the region give as input in ~IRAF format
    If the region can't be parsed (None, None) is returned
    """
    # peel off any outer brackets
    reg = re.sub(r"^\[*([^\]]*)\]*$", r"\1", reg)
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
    hduids = []  # list of candidate hduids
    for name in hdunames or []:
        try:
            hduid = hdulist.index_of(name)
            if hduid not in hduids:
                hduids.append(hduid)
        except KeyError as ke:
            logging.error("KeyError: %s", ke)
            logging.error("HDU[%s] not found, skipping", name)
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
    chduids = []  # list of candidate hduids
    for name in hdunames or []:
        try:
            hduid = hdulist.index_of(name)
            if hduid not in chduids:
                chduids.append(hduid)
        except KeyError as ke:
            logging.error("KeyError: %s", ke)
            logging.error("HDU[%s] not found, skipping", name)
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
    datasec = parse_region(dstr)
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


def subtract_bias(bstring: str, btype: str, hdu: fits.ImageHDU):
    """
    Subtract a bias estimate (using overscans) from an hdu.

    Choices are mean, median, byrow, byrowsmooth,  byrowcol or byrowcolsmooth
    subtraction of a bias calculated in either a given set of columns or using
    DATASEC to infer the overscan region. The rows-and-columns method uses both
    the serial and parallel oscan to account for the shape of the bias
    (particularly for ITL sensors)
    """
    (datasec, soscan, poscan) = get_data_oscan_slices(hdu)

    if bstring == "overscan":
        logging.debug("use overscan region in hdu:%s", hdu.name)
    elif len(bstring) >= 3:  # user supplies column selection
        # eg. a:b is minimum spec size, peel off outer brackets
        try:
            logging.debug("use %s cols in hdu:%s bias subtraction", bstring, hdu.name)
            reg = re.sub(r"^\[*([0-9]+:[0-9]+)\]*$", r"\1", bstring)
            (x1, x2) = re.match(r"([0-9]+):([0-9]+)$", reg).groups()
            soscan = (soscan[0], slice(int(x1), int(x2)))
        except SyntaxError as se:
            logging.error("SyntaxError: %s", se)
            logging.error("bad bias selection %s", bstring)
            sys.exit(1)
    else:
        logging.error("bad bias selection %s", bstring)
        sys.exit(1)

    logging.debug("biastype=%s", btype)
    # default method is "byrow"
    if btype in ("byrow", "byrowsmooth", "byrowcol", "byrowcolsmooth"):
        # get sigma clipped median per row
        so_avg, so_bias, so_std = stats.sigma_clipped_stats(hdu.data[soscan], axis=1)
        if re.search(r"smooth", btype):  # smooth the row bias vector a bit
            logging.debug("smoothing serial overscan with Gaussian1DKernel")
            gauss_kernel = Gaussian1DKernel(1)
            so_bias = convolve(so_bias, gauss_kernel)
        # convert shape from (n,) to (n, 1)
        so_bias = so_bias.reshape(np.shape(so_bias)[0], 1)
        logging.debug("np.shape(so_bias)=%s", np.shape(so_bias))
        hdu.data = hdu.data - so_bias.data
        if btype in ("byrowcol", "byrowcolsmooth"):
            # get sigma clipped median per column
            logging.debug(
                "poscan=((%d, %d), (%d, %d))",
                poscan[0].start,
                poscan[0].stop,
                poscan[1].start,
                poscan[1].stop,
            )
            po_avg, po_bias, po_std = stats.sigma_clipped_stats(
                hdu.data[poscan], axis=0
            )
            if re.search(r"smooth", btype):  # smooth col bias vector a bit
                logging.debug("smoothing par overscan with Gaussian1DKernel")
                gauss_kernel = Gaussian1DKernel(2)
                po_bias = convolve(po_bias, gauss_kernel)
            # convert shape from (,n) to (1, n)
            po_bias = po_bias.reshape(1, np.shape(po_bias)[0])
            logging.debug("np.shape(po_bias)=%s", np.shape(po_bias))
            hdu.data = hdu.data - po_bias.data
    elif btype == "mean":
        bias = np.mean(hdu.data[soscan])
        hdu.data = hdu.data - bias
    elif btype == "median":
        bias = np.median(hdu.data[soscan])
        hdu.data = hdu.data - bias
    else:
        logging.error("btype: %s not valid", btype)
        sys.exit(1)


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
    bstring: str,
    btype: str,
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
    bstring: param for subtract_bias() function (in this module)
    btype: param for subtract_bias()
    scaling: (yslice, xslice) specifying ROI to use for scaling
    hduo: a basic ImageHDU object that is modified and is the functions result
    """
    hdudata_list = []
    hdu_scale = []
    for im in iimages:
        hdu = im[hduid].copy()
        if bstring:
            subtract_bias(bstring, btype, hdu)
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
    if re.match(r"mea", method[0]):
        hduo.data = np.mean(np.array(hdudata_list), axis=0)
    elif re.match(r"med", method[0]):
        hduo.data = np.median(np.array(hdudata_list), axis=0)
    elif re.match(r"sig", method[0]):  # this one is ugly
        hduo.data = np.nanmean(
            stats.sigma_clip(np.array(hdudata_list), method[1], axis=0, masked=False),
            axis=0,
        )
    elif re.match(r"ran", method[0]):
        hduo.data = np.percentile(np.array(hdudata_list), method[1], axis=0)
    else:
        logging.error("image combine method %s not recognized", method[0])
        sys.exit(1)

    if bimage:
        hduo.data = hduo.data - bdata
