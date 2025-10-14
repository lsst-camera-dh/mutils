# mutils

This is a collection of command line utilities for working with fits images and CCS trending data.

## Installation

### With uv

1. [Install `uv`](https://docs.astral.sh/uv/getting-started/installation/)
    ```
    curl -LsSf https://astral.sh/uv/install.sh | sh
    ```
2. Create virtual env
    ```
    uv venv
    source .venv/bin/activate
    ```
3. Install `mutils`
    ```
    uv pip install git+https://github.com/lsst-camera-dh/mutils.git@uv-installable
    ```

All executables can be run with
```
uv run <script_name>
```
e.g.
```
uv run trender --fmt 'o' --lay 4x1 --site summit --plot --dur 6h --dpi 300 --overlayreg 'thermal/C.*/A.*' 'chiller1/Chiller/FluidTemperature' 'chiller1/Chiller/FlowRate' 'thermal/Trim.*/CryoTotal_P'
```

To mimic the old behavior of `trender.py`, one could use an alias `alias 'trender.py'='uv run trender'`

### With normal pip

1. Create a virtual env
    ```
    python -m venv .venv
    source .venv/bin/activate
    ```
2. Install `mutils`
    ```
    pip install git+https://github.com/lsst-camera-dh/mutils.git@uv-installable
    ```

## Command Line Utilities

These command line applications are written in python and rely mainly on astropy, numpy and matplotlib to make plots for FITS imagesand time series data from the CCS trending database.  They generally have a usage message obtained by using the '--help' option.

### trendutils

* trender -- This is the main trending application.  It is useful for retrieving and displaying trending data (time series).  It can display as plots, text, summary statistics or raw xml.
* trendmonitor -- Generates a snapshot text report for selected channels with hi/low limits and alarm state.  Originally used for a cron job sending the report via email.
* printXMLStructure -- Just a small script to display the structure of an xml file (sort of like the unix command tree).

### imutils

* imstat -- Display image statistics.  The canonical use is the "--quick" option to summarize stats for images.
* implot -- Generate plots for ROIs in images and display as one or more axes in a matplotlib figure.
* imhead -- Display fits headers as text.
* imfft -- Generate plots of the Image Noise Power Spectrum for fits images.  Useful for diagnosing noise problems.
* imarith -- Perform arithmetic on pairs of images or an image and scalars.
* imcombine -- Combine stacks of images with various methods to make super biases, flats, darks etc.
* hduorder -- Rearrange HDU order within a fits file (originally used for testing).
* imhist -- Not yet available.  Coming soon.  Plot ROI histograms for images.

## Module library

These are several modules with collections of useful functions used by the command line applications.

* mutils.py -- generic functions
* plotutils.py -- plotting related functions
* imutils.py -- image related functions
* trendutils.py -- trending data functions
* timezone_info.py -- timezone data to help allow time spec inputs in many timezones.
