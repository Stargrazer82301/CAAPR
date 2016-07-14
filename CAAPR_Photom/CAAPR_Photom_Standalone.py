# Import smorgasbord
import os
import sys
sys.path.insert(0, '/home/herdata/spx7cjc/Dropbox/Work/Scripts/')
import pdb
import time
import warnings
import numpy as np
import scipy.ndimage
import astropy.io.fits
import astropy.wcs
import astropy.convolution
import astropy.modeling
import astroquery.irsa_dust
import lmfit
import ChrisFuncs
import CAAPR






# Standalone wrapper around ApCorrect: A function that uses provided beam profile to aperture-correct photometry
def StandaloneApCorrect(psf_path, cutout, pix_arcsec, semimaj_pix, axial_ratio, angle, centre_i, centre_j, annulus_inner, annulus_outer):
    """
    Arguments of CAAPR.CAAPR_Photom_Standalone.ApCorrect:

    psf_path:           Either a string giving the path to FITS file that contains the PSF, or a False boolean (in which case an airy disc PSF will be assumed).
    cutout:             Array upon whcih photometry is being perfomred upon (should alreadt be background-subtracted).
    pix_arcsec:         The width, in arscec, of the pixels in the cutout (this is needed in case there is a pixel size mismatch with PSF).
    semimaj_pix:        Semi-major axis of photometric aperture, in pixels.
    axial_ratio:        Axial ratio of photometryic aperture.
    angle:              Position angle of photometric aperture, in degrees.
    centre_i:           Zero-indexed, 0th-axis coordinate (equivalent to y-axis one-indexed coordinates in FITS terms) of centre position of photometric aperture.
    centre_j:           Zero-indexed, 1st-axis coordinate (equivalent to x-axis one-indexed coordinates in FITS terms) of centre position of photometric aperture.
    annulus_inner:      The semi-major axis of the inner edge of the background annulus, in units of the semi-major axis of the source ellipse.
    annulus_outer:      The semi-major axis of the outer edge of the background annulus, in units of the semi-major axis of the source ellipse.

    Returns:

    ap_correction:      The calculated aperture correction factor
    """


    # Construct staw-man pod
    pod = {'id':'',
           'cutout':cutout,
           'pix_arcsec':pix_arcsec,
           'adj_semimaj_pix':semimaj_pix,
           'adj_axial_ratio':axial_ratio,
           'adj_angle':angle,
           'centre_i':centre_i,
           'centre_j':centre_j,
           'ap_sum':1.0,
           'ap_error':0.1,
           'bg_avg':0.0}

    # Construct straw-man source dictionary
    source_dict = {}

    # Construct straw-man band dictionary
    band_dict = {'beam_correction':psf_path,
                 'annulus_inner':annulus_inner,
                 'annulus_outer':annulus_outer,
                 'subpixel_factor':1.0}

    # Construct straw-man kwargs dictionary
    kwargs_dict = {'verbose':False,
                   'extinction_corr':True}

    # Call actual ApCorrect function
    pod = CAAPR.CAAPR_Photom.ApCorrect(pod, source_dict, band_dict, kwargs_dict)

    # Determine aperture correction factor, and return
    ap_correction = pod['ap_sum']
    return ap_correction





# Standalone wrapper around ExtCorrect: A function that performs extinction correction on photometry, via IRSA dust extinction service (which uses the Schlafly & Finkbeiner 2011 prescription)
def StandaloneExtCorrrct(ra, dec, band_name):
    """
    Arguments of CAAPR.CAAPR_Photom_Standalone.ExtCorrrct:

    ra:                 Right ascenscion of the target coorindate.
    dec:                Declination of the target coordinate.
    band_name:          The name of the band in question.

    Returns:

    excorr:             The determined Galactic extionction correction, in magnititudes, for band in question at the target coordinates
    """



    # Construct staw-man pod
    pod = {'id':'',
           'ap_sum':1.0,
           'ap_error':1.0}

    # Construct straw-man source dictionary
    source_dict = {'ra':ra,
                   'dec':dec}

    # Construct straw-man band dictionary
    band_dict = {'band_name':band_name}

    # Construct straw-man kwargs dictionary
    kwargs_dict = {'verbose':False,
                   'extinction_corr':True}

    # Call actual ExtCorrect function
    pod = CAAPR.CAAPR_Photom.ExtCorrrct(pod, source_dict, band_dict, kwargs_dict)

    # Determine extinction correction, and return
    excorr = np.abs( 2.51 * np.log10( pod['ap_sum'] ) )
    return excorr






#psf_path = '/home/herdata/spx7cjc/Beams/SPIRE_250.fits'
#cutout = astropy.io.fits.getdata('/home/saruman/spx7cjc/DustPedia/SPIRE/Cutouts/DustPedia/NGC4030_SPIRE_250.fits')
#pix_arcsec = 6.0
#semimaj_pix = 41.0
#axial_ratio = 1.1795263352195566
#angle = 115.16660752050387
#centre_i = 300.0
#centre_j = 300.0
#annulus_inner = 1.25
#annulus_outer = 1.601

#ap_corr = CAAPR.CAAPR_Photom.CAAPR_Photom_Standalone.StandaloneApCorrect(psf_path, cutout, pix_arcsec, semimaj_pix, axial_ratio, angle, centre_i, centre_j, annulus_inner, annulus_outer)
