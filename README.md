## mutils

This is a collection of command line utilities for working with fits images and CCS trending data.

### Command Line Utilities

### trendutils

* trender.py -- This is the main trending application.  It is useful for retrieving and displaying trending data.  It can display as plots, text, summary statistics or raw xml.
* trendmonitor.py -- Generates a snapshot text report for selected channels with hi/low limits and alarm state.  Originally used for a cronjob sending the report via email.
* printXMLStructure.py -- Just a small script to display the structure of an xml file (sort of like the unix command tree).

### imutils

* imstat.py -- Display image statistics.  The canonical use is the "--quick" option to summarize stats for images.
* implot.py -- Generate plots for ROIs in images and display as one or more axes in a matplotlib figure.
* imhead.py -- Display fits headers as text.
* imfft.py -- Generate plots of the Image Noise Power Spectrum for fits images.  Useful for diagnosing noise problems.
* imarith.py -- Perform arithmetic on pairs of images or an image and scalars.
* imcombine.py -- Combine stacks of images with various methods to make super biases, flats, darks etc.
* hduorder.py -- Rearrange HDU order within a fits file (originally used for testing).
* imhist.py -- Not yet available.  Coming soon.  Plot ROI histograms for images.
