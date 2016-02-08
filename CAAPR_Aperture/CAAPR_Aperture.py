# Import smorgasbord
import os
import sys
sys.path.insert(0, '../')
import gc
import pdb
import numpy as np
import astropy.io.fits
import astropy.wcs
import ChrisFuncs



# The aperture-fitting sub-pipeline
def PipelineAperture(source_dict, band_dict, output_dir_path, temp_dir_path):

    # Read in FITS file in question
    in_fitspath = band_dict['band_path']
    in_fitsdata = astropy.io.fits.open(in_fitspath)
    in_image = in_fitsdata[0].data
    in_header = in_fitsdata[0].header
    in_fitsdata.close()
    in_wcs = astropy.wcs.WCS(in_header)
    pix_size = in_wcs.wcs.cdelt

