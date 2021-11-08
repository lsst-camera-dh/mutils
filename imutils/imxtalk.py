#!/usr/bin/env python
"""
Calculate xtalk coefs for a single CCD with N segments
in the form of a NxN matrix
"""
import os
import sys
import re
import argparse
import logging
import textwrap
import os.path
from astropy.io import fits
import numpy as np
from sklearn import linear_model
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import binned_statistic


# put parent directory into sys.path
bp = os.path.dirname(os.path.realpath(__file__)).split(os.sep)
modpath = os.sep.join(bp[:-1] + ["lib"])
sys.path.insert(0, modpath)

#  local imports
try:
    import imutils as iu
    import mutils as mu
    import plotutils as pu
except ImportError as e:
    logging.error("Import failed: %s", e)
    exit(1)


def parse_args():
    """handle command line"""
    style_list = ["default"] + sorted(
        ["fast", "ggplot", "seaborn-poster", "seaborn-notebook"]
    )
    #
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """\
            Calculate cross-talk coefs for a single CCD")
                                    """
        ),
        epilog=textwrap.dedent(
            """\
                               """
        ),
    )
    parser.add_argument(
        "fitsfile", nargs="+", metavar="file", help="input fits file(s)"
    )
    # ---------------------------------------------------------------
    # restrict sources to a subset of hdu's
    sgroup = parser.add_mutually_exclusive_group()
    sgroup.add_argument(
        "--srcname", nargs="+", metavar="idn", help="process HDU list by names"
    )
    sgroup.add_argument(
        "--srcindex", nargs="+", type=int, metavar="idx", help="process HDU list by ids"
    )
    # ---------------------------------------------------------------
    # restrict responses to a subset of hdu's
    rgroup = parser.add_mutually_exclusive_group()
    rgroup.add_argument(
        "--rspname", nargs="+", metavar="idn", help="process HDU list by names"
    )
    rgroup.add_argument(
        "--rspindex", nargs="+", type=int, metavar="idx", help="process HDU list by ids"
    )
    # ---------------------------------------------------------------
    parser.add_argument(
        "--intercept", action="store_true", help="fit line with intercept"
    )
    parser.add_argument(
        "--threshold",
        nargs=1,
        metavar="thresh",
        type=float,
        help="threshold for source pixels (defaults to 500 * read-noise-estimate",
    )
    parser.add_argument(
        "--background",
        nargs="?",
        metavar="method_spec",
        const="boxcar",
        help="subtract background: {default=median}",
    )
    parser.add_argument(
        "--plot", action="store_true", help="plot binned (src, rsp) scatter plots"
    )
    parser.add_argument(
        "--ylimits",
        nargs=2,
        type=float,
        required=False,
        help="lower upper",
    )
    parser.add_argument(
        "--predict", action="store_true", help="plot fit prediction lines"
    )
    parser.add_argument(
        "--raw", action="store_true", help="plot raw (src, rsp) scatter plots"
    )
    parser.add_argument(
        "--layout", default="landscape", help='"landscape"|"portrait"|"nxm"'
    )
    parser.add_argument(
        "--style",
        default="ggplot",
        required=False,
        choices=style_list,
        help="default: %(default)s",
    )
    # x-axis matplotlib sharex exclusive
    xgroup = parser.add_mutually_exclusive_group()
    xgroup.add_argument(
        "--sharex",
        dest="sharex",
        action="store_true",
        help="Share X-axis range, this is default",
    )
    xgroup.add_argument(
        "--no-sharex",
        dest="sharex",
        action="store_false",
        help="Don't share Y-axis range",
    )
    xgroup.set_defaults(sharex=False)
    # y-axis matplotlib sharey exclusive
    ygroup = parser.add_mutually_exclusive_group()
    ygroup.add_argument(
        "--sharey",
        dest="sharey",
        action="store_true",
        help="Share Y-axis range, this is default",
    )
    ygroup.add_argument(
        "--no-sharey",
        dest="sharey",
        action="store_false",
        help="Don't share Y-axis range",
    )
    ygroup.set_defaults(sharey=True)
    # title
    parser.add_argument(
        "--title", nargs="?", metavar="Title", const="auto", help="specify Title String"
    )
    parser.add_argument(
        "--info", action="store_true", help="print the info() table summarizing file"
    )
    parser.add_argument(
        "--debug", action="store_true", help="print additional debugging messages"
    )
    return parser.parse_args()


def imxtalk():
    """main logic:"""
    optlist = parse_args()
    mu.init_logging(optlist.debug)
    mu.init_warnings()
    ncalls.counter = 0
    # begin processing -- loop over files
    for ffile in optlist.fitsfile:
        try:
            hdulist = fits.open(ffile)
        except IOError as ioerr:
            logging.error("IOError: %s", ioerr)
            exit(1)
        if optlist.info:  # just print the image info per file
            hdulist.info()
            continue
        # Construct a list of the source HDU's to work on
        srcids = iu.get_requested_image_hduids(
            hdulist, optlist.srcname, optlist.srcindex
        )
        # Construct a list of the response HDU's to work on
        rspids = iu.get_requested_image_hduids(
            hdulist, optlist.rspname, optlist.rspindex
        )
        # do stuff
        logging.debug("calling iu.get_union_of_bad_column_segs(hdulist)")
        bad_segs = iu.get_union_of_bad_column_segs(hdulist)
        logging.debug("bad_segs=%s", bad_segs)
        max_rn = 7.0
        pcnt = 20
        lsst_num = hdulist[0].header.get("LSST_NUM")
        # get_xtalk_coefs(hdulist, srcids, rspids, optlist.threshold)
        for srcid in srcids:
            hdu_s = hdulist[srcid]
            # subtract the bias estimate from the source array
            if lsst_num and re.match(r"^E2V-CCD250", lsst_num):
                sbias = "byrowe2v"
            else:
                sbias = "byrow"
            pbias = "bycolfilter"
            iu.subtract_bias(sbias, pbias, hdu_s, bad_segs)
            logging.info("hdu_s = %s", hdu_s.header["EXTNAME"])
            (datasec_s, soscan_s, poscan_s) = iu.get_data_oscan_slices(hdu_s)
            rn_est = min(np.std(hdu_s.data[poscan_s[0], soscan_s[1]]), max_rn)
            if optlist.threshold:
                thresh = optlist.threshold
            else:
                rn_est = min(np.std(hdu_s.data[poscan_s[0], soscan_s[1]]), max_rn)
                thresh = 500 * rn_est
            # estimate source background level, the threshold will be added to that
            # since we are only interested in source pixels above background by thresh
            thresh += np.percentile(hdu_s.data[datasec_s], pcnt)

            # make the (weights) mask used in response hdu bckgrnd subtraction
            mask_s = np.ones_like(hdu_s.data, dtype=int)
            mask_s[np.nonzero(hdu_s.data > thresh)] = 0
            mask_s[:, bad_segs] = 0  # fold these in

            arr_s = hdu_s.data.flatten("K")
            logging.debug("np.shape(arr_s)= %s", np.shape(arr_s))
            arr_x = arr_s[arr_s > thresh]
            logging.debug("found %d nans in arr_x", np.count_nonzero(np.isnan(arr_x)))
            arr_x = arr_x.reshape(-1, 1)  # infer 1st axis, 2nd axis for 1 "feature"
            logging.debug("np.shape(arr_x)= %s", np.shape(arr_x))
            if np.size(arr_x) < 1000:
                logging.warning(
                    "not enough source points to produce a coef: %d < 100",
                    np.size(arr_x),
                )
                continue

            if optlist.plot:
                plt.style.use(optlist.style)
                pu.update_rcparams()
                fig, axes = pu.get_fig_and_axes(
                    len(rspids),
                    optlist.layout,
                    False,
                    optlist.sharex,
                    optlist.sharey,
                    None,
                )
                nprows, npcols = (axes.shape[0], axes.shape[1])
                pu.set_fig_title(optlist.title, ffile, fig)
                sylim_upper = sylim_lower = 0.0

            for rindex, rspid in enumerate(rspids):
                if rspid == srcid:
                    sindex = rindex
                    continue
                hdu_r = hdulist[rspid]
                logging.info("    hdu_r = %s", hdu_r.header["EXTNAME"])
                if np.shape(hdu_s.data) != np.shape(hdu_r.data):
                    logging.warning(
                        "hdu's %s, %s shapes not commensurate: %s != %s, skipping",
                        hdu_s.header["EXTNAME"],
                        hdu_r.header["EXTNAME"],
                        np.shape(hdu_s.data),
                        np.shape(hdu_r.data),
                    )
                    continue
                (datasec_r, soscan_r, poscan_r) = iu.get_data_oscan_slices(hdu_r)
                iu.subtract_bias(sbias, pbias, hdu_r, bad_segs)
                # need to subtract background level estimate from hdu_s but it may
                # have lots of structure so need somewhat careful estimate
                # ------------
                # redo this with masking and line by line interp across the mask
                logging.debug(
                    "found %d nans in hdu_r", np.count_nonzero(np.isnan(hdu_r.data))
                )
                iu.subtract_background_for_xtalk(hdu_r, mask_s, datasec_r)
                logging.debug(
                    "found %d nans in hdu_r", np.count_nonzero(np.isnan(hdu_r.data))
                )
                arr_r = hdu_r.data.flatten("K")
                logging.debug("np.shape(arr_r)= %s", np.shape(arr_r))
                arr_y = arr_r[arr_s > thresh]

                logging.debug(
                    "found %d nans in arr_y", np.count_nonzero(np.isnan(arr_y))
                )
                arr_xp = arr_x[~np.isnan(arr_y)]
                arr_yp = arr_y[~np.isnan(arr_y)]

                # reject high sources in response channel
                arr_xp = arr_xp[arr_yp < thresh]
                arr_yp = arr_yp[arr_yp < thresh]

                if optlist.intercept:
                    lr = linear_model.LinearRegression()
                    ransac = linear_model.RANSACRegressor()
                else:
                    lr = linear_model.LinearRegression(fit_intercept=False)
                    ransac = linear_model.RANSACRegressor(
                        linear_model.LinearRegression(fit_intercept=False)
                    )

                # lr.fit(arr_xp, arr_yp)
                # ransac.fit(arr_xp, arr_yp)
                # print(f"lr.coef={lr.coef_}")
                # print(f"ransac.estimator.coef={ransac.estimator_.coef_}")
                if np.max(arr_xp) < 0.95 * np.max(arr_x):
                    logging.warning("threshold is too low, raise and re-run")
                nbins = (np.max(arr_xp) - np.min(arr_xp)) / 1000 * rn_est
                logging.debug("np.max(arr_xp) = %.2f", np.max(arr_xp))
                logging.debug("np.min(arr_xp) = %.2f", np.min(arr_xp))
                logging.debug("nbins = %d", nbins)
                s, edges, _ = binned_statistic(arr_xp[:, 0], arr_yp, "median", nbins)
                cnt, cedges, _ = binned_statistic(arr_xp[:, 0], arr_yp, "count", nbins)
                bin_width = edges[1] - edges[0]
                logging.debug("bin_width = %.2f", bin_width)
                binx = edges[1:] - bin_width / 2
                binx = binx[~np.isnan(s)]  # remove the empty bins
                count = cnt[~np.isnan(s)]
                logging.debug(
                    "count: mean: %.2f  median: %.2f  stddev: %.2f min: %.2f  max: %.2f",
                    np.mean(count),
                    np.median(count),
                    np.std(count),
                    np.min(count),
                    np.max(count),
                )
                count = np.sqrt(count) * (np.log10(binx))  # extra weight on high bins
                logging.debug("binx[-10:] = %s", binx[-10:])
                logging.debug("count[-10:] = %s", count[-10:])
                s = s[~np.isnan(s)]
                sylim_upper = np.percentile(arr_yp, 99)
                sylim_lower = np.percentile(arr_yp, 1)
                if optlist.raw:  # expand limits
                    sydel = sylim_upper - sylim_lower
                    sylim_upper += 0.4 * sydel
                    sylim_lower -= 0.3 * sydel

                binx = binx.reshape(-1, 1)  # infer 1st axis, 2nd axis for 1 "feature"
                lr.fit(binx, s, sample_weight=count)
                ransac.fit(binx, s, count)
                inlier_mask = ransac.inlier_mask_
                outlier_mask = np.logical_not(inlier_mask)
                lrcoef = lr.coef_[0]
                lrincpt = lr.intercept_
                rscoef = ransac.estimator_.coef_[0]
                rsincpt = ransac.estimator_.intercept_
                print(f"binned lr.coef={lrcoef:>3g}")
                if optlist.intercept:
                    print(f"binned lr.intercept={lrincpt:>.2f}")
                print(f"binned ransac.estimator.coef={rscoef:>3g}")
                if optlist.intercept:
                    print(f"binned ransac.estimator.intercept={rsincpt:>.2f}")

                # plotting
                if optlist.plot:
                    ax = np.ravel(axes)[int(rindex / npcols) * npcols + rindex % npcols]
                    if optlist.style == "ggplot":
                        ax.scatter([], [])  # skip the first color
                    ax.grid(True)
                    ax.set_xlabel("source signal", size="x-small")
                    ax.set_ylabel("response signal", size="x-small")
                    if optlist.raw:
                        ax.scatter(arr_xp[:, 0], arr_yp, s=1.0, label="raw")
                    si = s[inlier_mask]
                    ci = count[inlier_mask]
                    # ci[-1] = 1.0
                    ax.scatter(
                        binx[inlier_mask, 0],
                        si,
                        # s=np.sqrt(count),
                        s=np.sqrt(ci),
                        color="blue",
                        alpha=0.5,
                        label="inliers",
                    )
                    si = s[outlier_mask]
                    ci = count[outlier_mask]
                    # ci[-1] = 1.0
                    ax.scatter(
                        binx[outlier_mask, 0],
                        si,
                        # s=np.sqrt(count),
                        s=np.sqrt(ci),
                        color="purple",
                        alpha=0.5,
                        label="outliers",
                    )
                    if optlist.predict:  # Predict and plot result of estimated models
                        line_x = np.arange(0.0, binx.max())[:, np.newaxis]
                        line_y = lr.predict(line_x)
                        line_y_ransac = ransac.predict(line_x)
                        lw = 2
                        if optlist.intercept:
                            lbl = f"lr: {lrcoef:>.3g}*x + {lrincpt:>.3g}"
                        else:
                            lbl = f"lr: {lrcoef:>.3g}*x"
                        ax.plot(line_x, line_y, color="navy", linewidth=lw, label=lbl)
                        if optlist.intercept:
                            lbl = f"ransac: {rscoef:>.3g}*x + {rsincpt:>.3g}"
                        else:
                            lbl = f"ransac: {rscoef:>.3g}*x"
                        ax.plot(
                            line_x,
                            line_y_ransac,
                            color="cornflowerblue",
                            linewidth=lw,
                            label=lbl,
                        )

                    if optlist.ylimits:
                        ax.set_ylim(optlist.ylimits[0], optlist.ylimits[1])
                    else:
                        ax.set_ylim(sylim_lower, sylim_upper)
                    ax.xaxis.set_tick_params(labelsize="x-small")
                    ax.xaxis.set_major_formatter(ticker.EngFormatter())
                    ax.yaxis.set_tick_params(labelsize="x-small")
                    ax.set_title(
                        f"SRC:RSP {hdu_s.header['EXTNAME']}:{hdu_r.header['EXTNAME']}",
                        fontsize="xx-small",
                    )
                    handles, labels = ax.get_legend_handles_labels()
                    lgnd = pu.mk_legend("inside", nprows, handles, labels, ax)
                    # big hack
                    lgnd.legendHandles[-2]._sizes = [6]
                    lgnd.legendHandles[-1]._sizes = [6]

            if optlist.plot:
                for gidx in range(rindex + 1, nprows * npcols):
                    # ax = np.ravel(axes)[int(gidx / npcols) * npcols + gidx % npcols]
                    ax = np.ravel(axes)[gidx]
                    ax.grid(False)
                    ax.set_frame_on(False)
                    ax.get_xaxis().set_visible(False)
                    ax.get_yaxis().set_visible(False)

                if srcid in rspids:
                    ax = np.ravel(axes)[sindex]
                    ax.grid(False)
                    ax.set_frame_on(False)
                    ax.get_xaxis().set_visible(False)
                    ax.get_yaxis().set_visible(False)

                fig.set_tight_layout(
                    {"h_pad": 0.50, "w_pad": 1.0, "rect": [0, 0, 1, 0.97]}
                )
                plt.show()

        # end of doing stuff
        ncalls.counter = 0  # reset per file, triggers headers

    ncalls()  # track call count, acts like static variable)


def get_xtalk_coefs(hdulist: fits.HDUList, srcids, rspids, thresh):
    """ """


def ncalls():
    """maintain a counter"""
    ncalls.counter += 1


if __name__ == "__main__":
    imxtalk()
