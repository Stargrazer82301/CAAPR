# CAAPR (the Comprehensive & Adaptable Aperture Photometry Routine)

CAAPR (the Comprehensive & Adaptable Aperture Photometry Routine) is an aperture photometry pipeline, designed with multiwavelength extragalactic astronomy in mind. It's written in Python 2.7. (Yes, I know that Python 2.7 is *soo* 2016; I'll get around to updating it some day. Meanwhile, if you want to avoid compatibility issues, I recommend using a [virtual environment](https://conda.io/docs/user-guide/tasks/manage-environments.html); it's surprisingly easy.)

CAAPR was created to be able to generate aperture-matched photometry that is truly cross-comparable, even in the face of the broad range of data types that are required for multiwavelength astronomy – producing fluxes and **uncertainties** in a consistent manner, despite the enormous variation in the characteristics of observations ranging from ultraviolet to infrared to microwave.

CAAPR is designed to be an end-to-end photometry pipeline; the user sets it running, and it outputs a couple of CSV tables describing the results, along with some thumbnail images illustrating the apertures. Of course, astronomers tend to distrust "black box" code, so CAAPR is designed with a host of options and switches that the user can adjust, if they want. This makes CAAPR more of a "grey box" pipeline.

## Installation

1. [Download](https://github.com/Stargrazer82301/CAAPR/archive/master.zip) the zip archive containing CAAPR from the GitHub repository, and extract it to a temporary location of your choice.
2. Open a terminal and navigate to extracted directory that contains setup.py
3. Install using the terminal command `pip install -e . --user` (Whilst you can try installing using the more common `python setup.py install` route, this sometimes hits compatability issues, whereas the pip local installer seems more robust.)
4. Download and install the ChrisFuncs package, found [here](https://github.com/Stargrazer82301/ChrisFuncs), which CAAPR depends upon. Like CAAPR it can be downloaded, extracted, and installed very straightforwardly.
