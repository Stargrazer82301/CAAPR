#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.headers Contains functions for extracting information from FITS headers.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import re
import copy
import numpy as np

# Import astronomical modules
from astropy import coordinates
from astropy import units as u

# Import the relevant PTS classes and modules
from ..basics.vector import Extent
from ...core.basics.filter import Filter
from ..basics.coordinatesystem import CoordinateSystem
from ...core.tools.logging import log

# -----------------------------------------------------------------

def flattened(header):

    """
    This function ...
    :param header:
    :return:
    """

    flat_header = copy.deepcopy(header)
    flat_header["NAXIS"] = 2
    if "NAXIS3" in flat_header: del flat_header["NAXIS3"]
    for key in flat_header:
        if "PLANE" in key: del flat_header[key]

    return flat_header

# -----------------------------------------------------------------

def get_pixelscale(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Search for different keywords indicating the pixelscale
    for keyword in ("PIXSCALE", "SECPIX", "PFOV", "PLTSCALE"):

        if keyword in header:

            # Get the scale and try to get the unit
            scale = header[keyword]
            if scale == "N/A": continue

            try: unit = header.comments[keyword].split("[")[1].split("]")[0]
            except IndexError: unit = None

            # Parse the unit with Astropy
            if unit is not None:

                unit = unit.replace("asec", "arcsec")
                if not (unit.endswith("pixel") or unit.endswith("pix")): unit = unit + "/pix"
                try: unit = u.Unit(unit)
                except ValueError: unit = None

            log.debug("pixelscale found in " + str(keyword) + " keyword = " + str(scale))
            log.debug("unit for the pixelscale = " + str(unit))

            # If no unit is found, guess that it's arcseconds / pixel ...
            if unit is None: unit = u.Unit("arcsec/pix")
            scale = scale * unit

            # Return the scale
            return Extent(scale, scale)

    # Search for the 'PXSCAL1' and 'PXSCAL2' keywords
    for keyword_combination in (("PXSCAL1", "PXSCAL2"), ("XPIXSIZE", "YPIXSIZE")):

        if keyword_combination[0] in header and keyword_combination[1] in header:

            scale1 = header[keyword_combination[0]]
            if scale1 == "N/A": continue
            try: unit1 = header.comments[keyword_combination[0]].split("[")[1].split("]")[0]
            except IndexError: unit1 = None

            scale2 = header[keyword_combination[1]]
            if scale2 == "N/A": continue
            try: unit2 = header.comments[keyword_combination[1]].split("[")[1].split("]")[0]
            except IndexError: unit2 = None

            # Parse the unit with Astropy
            if unit1 is not None:

                unit1 = unit1.replace("asec", "arcsec")
                if not (unit1.endswith("pixel") or unit1.endswith("pix")): unit1 = unit1 + "/pix"
                try: unit1 = u.Unit(unit1)
                except ValueError: unit1 = None

            if unit2 is not None:

                unit2 = unit2.replace("asec", "arcsec")
                if not (unit2.endswith("pixel") or unit2.endswith("pix")): unit2 = unit2 + "/pix"
                try: unit2 = u.Unit(unit2)
                except ValueError: unit2 = None

            log.debug("pixelscale found in PXSCAL1 and PXSCAL2 keywords = (" + str(scale1) + "," +str(scale2) + ")")
            log.debug("unit for the pixelscale = (" + str(unit1) + "," + str(unit2) + ")")

            # If no unit is found, guess that it's arcseconds / pixel ...
            if unit1 is None: unit1 = u.Unit("arcsec/pix")
            if unit2 is None: unit2 = u.Unit("arcsec/pix")
            scale1 = scale1 * unit1
            scale2 = scale2 * unit2

            # Return the pixelscale
            return Extent(scale1, scale2)

    # If none of the above keywords were found, return None
    else: return None

# -----------------------------------------------------------------

def get_filter(name, header=None):

    """
    This function ...
    :param name:
    :param header:
    :return:
    """

    filterid = name.lower()
    channel = None
    wavelength = None

    if "kernel" in filterid:
        log.debug("The image represents a kernel, so no filter will be set")
        return None

    # Get information from the header
    if header is not None:

        # Get information regarding the telescope and instrument
        if "TELESCOP" in header: filterid += " " + get_string(header["TELESCOP"]).lower()
        if "INSTRUME" in header: filterid += " " + get_string(header["INSTRUME"]).lower()
        if "ORIGIN" in header: filterid += " " + get_string(header["ORIGIN"]).lower()
        if "OBSERVAT" in header: filterid += " " + get_string(header["OBSERVAT"]).lower()

        # Get a name describing the filter
        if "FILTER" in header:
            try:
                filter = Filter.from_string(header["FILTER"])
                return filter
            except ValueError: pass
            filterid += " " + get_string(header['FILTER']).lower()
        if "FLTRNM" in header:
            try:
                filter = Filter.from_string(header["FLTRNM"])
                return filter
            except ValueError: pass
            filterid += " " + get_string(header['FLTRNM']).lower()

        # Get information about the channel number
        if "CHNLNUM" in header: channel = get_int(header["CHNLNUM"])
        elif "BAND" in header: channel = get_int(header["BAND"])
        else: channel = None

        # Get the wavelength
        if "WAVELEN" in header: wavelength = get_quantity(header["WAVELEN"], "micron")
        elif "WVLNGTH" in header: wavelength = get_quantity(header["WVLNGTH"], "micron")
        else: wavelength = None

    # Debug information
    log.debug("filterid = " + str(filterid))
    log.debug("channel = " + str(channel))
    log.debug("wavelength = " + str(wavelength))

    final_filter_name = None

    # -- UV --

    # GALEX
    if "fuv" in filterid: final_filter_name = "GALEX FUV"
    elif "nuv" in filterid: final_filter_name = "GALEX NUV"

    # SWIFT
    elif "uw2" in filterid: final_filter_name = "SWIFT W2"
    elif "um2" in filterid: final_filter_name = "SWIFT M2"
    elif "uw1" in filterid: final_filter_name = "SWIFT W1"
    elif "swift" in filterid or "uvot" in filterid:

        if "w2" in filterid: final_filter_name = "SWIFT W2"
        elif "m2" in filterid: final_filter_name = "SWIFT M2"
        elif "w1" in filterid: final_filter_name = "SWIFT W1"
        else: log.warning("Could not determine which SWIFT UVOT filter was used for this image")

    # TODO: support other UV instruments

    # -- Optical --

    # SDSS
    elif "sdss" in filterid:

        if "-u" in filterid: final_filter_name = "SDSS u"
        elif "-g" in filterid: final_filter_name = "SDSS g"
        elif "-r" in filterid: final_filter_name = "SDSS r"
        elif "-i" in filterid: final_filter_name = "SDSS i"
        elif "-z" in filterid: final_filter_name = "SDSS z"
        else:

            if "sdss u" in filterid: final_filter_name = "SDSS u"
            elif "sdss g" in filterid: final_filter_name = "SDSS g"
            elif "sdss r" in filterid: final_filter_name = "SDSS r"
            elif "sdss i" in filterid: final_filter_name = "SDSS i"
            elif "sdss z" in filterid: final_filter_name = "SDSS z"
            else:

                if "u" in name: final_filter_name = "SDSS u"
                elif "g" in name: final_filter_name = "SDSS g"
                elif "r" in name: final_filter_name = "SDSS r"
                elif "i" in name: final_filter_name = "SDSS i"
                elif "z" in name: final_filter_name = "SDSS z"
                else: log.warning("Could not determine which SDSS filter was used for this image")

    # R band // not good; H alpha image was also identified as R band ...
    #elif "r" in filterid and "kpno" in filterid: return Filter("KPNO.Mosaic.R")

    # TODO: support other optical instruments

    # -- IR --

    # 2MASS filters
    elif "2mass" in filterid:

        if "h" in filterid: final_filter_name = "2MASS H"
        elif "j" in filterid: final_filter_name = "2MASS J"
        elif "k" in filterid: final_filter_name = "2MASS Ks"
        else: log.warning("Could not determine which 2MASS filter was used for this image")

    # IRAC filters
    elif "irac" in filterid:

        if "3.6" in filterid or "i1" in filterid: final_filter_name = "IRAC I1"
        elif "4.5" in filterid or "i2" in filterid: final_filter_name = "IRAC I2"
        elif "5.8" in filterid or "i3" in filterid: final_filter_name = "IRAC I3"
        elif "8.0" in filterid or "i4" in filterid: final_filter_name = "IRAC I4"
        else:  # Look at the channel number

            if channel is not None:

                if channel == 1: final_filter_name = "IRAC I1"
                elif channel == 2: final_filter_name = "IRAC I2"
                elif channel == 3: final_filter_name = "IRAC I3"
                elif channel == 4: final_filter_name = "IRAC I4"
                else: log.warning("Could not determine which IRAC filter was used for this image")

            elif wavelength is not None:

                if np.isclose(wavelength.to("micron").value, 3.6, rtol=0.05): final_filter_name = "IRAC I1"
                elif np.isclose(wavelength.to("micron").value, 4.5, rtol=0.05): final_filter_name = "IRAC I2"
                elif np.isclose(wavelength.to("micron").value, 5.8, rtol=0.05): final_filter_name = "IRAC I3"
                elif np.isclose(wavelength.to("micron").value, 8.0, rtol=0.05): final_filter_name = "IRAC I4"
                else: log.warning("Could not determine which IRAC filter was used for this image")

            else: log.warning("Could not determine which IRAC filter was used for this image")

    # WISE filters
    elif "wise" in filterid:

        if "w1" in filterid: final_filter_name = "WISE W1"
        elif "w2" in filterid: final_filter_name = "WISE W2"
        elif "w3" in filterid: final_filter_name = "WISE W3"
        elif "w4" in filterid: final_filter_name = "WISE W4"
        else:

            if channel is not None:

                if channel == 1: final_filter_name = "WISE W1"
                elif channel == 2: final_filter_name = "WISE W2"
                elif channel == 3: final_filter_name = "WISE W3"
                elif channel == 4: final_filter_name = "WISE W4"
                else: log.warning("Could not determine which WISE filter was used for this image")

            elif wavelength is not None:

                if np.isclose(wavelength.to("micron").value, 3.4, rtol=0.05): final_filter_name = "WISE W1"
                elif np.isclose(wavelength.to("micron").value, 4.6, rtol=0.05): final_filter_name = "WISE W2"
                elif np.isclose(wavelength.to("micron").value, 12., rtol=0.05): final_filter_name = "WISE W3"
                elif np.isclose(wavelength.to("micron").value, 22., rtol=0.05): final_filter_name = "WISE W4"

            else: log.warning("Could not determine which WISE filter was used for this image")

    # MIPS filters
    elif "mips" in filterid:

        if "24" in filterid: final_filter_name = "MIPS 24"
        elif "70" in filterid: final_filter_name = "MIPS 70"
        elif "160" in filterid: final_filter_name = "MIPS 160"
        else:

            if wavelength is not None:

                if np.isclose(wavelength.to("micron").value, 24., rtol=0.05): final_filter_name = "MIPS 24"
                elif np.isclose(wavelength.to("micron").value, 70., rtol=0.05): final_filter_name = "MIPS 70"
                elif np.isclose(wavelength.to("micron").value, 160., rtol=0.05): final_filter_name = "MIPS 160"
                else: log.warning("Could not determine which MIPS filter was used for this image")

            else: log.warning("Could not determine which MIPS filter was used for this image")

    # Spitzer bands (IRAC and MIPS but "irac" and "mips" are not in filterid)
    elif "spitzer" in filterid:

        if wavelength is not None:

            if np.isclose(wavelength.to("micron").value, 3.6, rtol=0.05): final_filter_name = "IRAC I1"
            elif np.isclose(wavelength.to("micron").value, 4.5, rtol=0.05): final_filter_name = "IRAC I2"
            elif np.isclose(wavelength.to("micron").value, 5.8, rtol=0.05): final_filter_name = "IRAC I3"
            elif np.isclose(wavelength.to("micron").value, 8.0, rtol=0.05): final_filter_name = "IRAC I4"
            elif np.isclose(wavelength.to("micron").value, 24., rtol=0.05): final_filter_name = "MIPS 24"
            elif np.isclose(wavelength.to("micron").value, 70., rtol=0.05): final_filter_name = "MIPS 70"
            elif np.isclose(wavelength.to("micron").value, 160., rtol=0.05): final_filter_name = "MIPS 160"
            else: log.warning("Could not determine which Spitzer filter was used for this image")

        else: log.warning("Could not determine which Spitzer filter was used for this image")

    # PACS filters
    elif "pacs" in filterid:

        if '70' in filterid or 'blue' in filterid: final_filter_name = "Pacs blue"
        elif '100' in filterid or 'green' in filterid: final_filter_name = "Pacs green"
        elif '160' in filterid or 'red' in filterid: final_filter_name = "Pacs red"
        else: log.warning("Could not determine which PACS filter was used for this image")

    # SPIRE filters
    elif "spire" in filterid:

        if "psw" in filterid or "250" in filterid: final_filter_name = "SPIRE PSW"
        elif "pmw" in filterid or "350" in filterid: final_filter_name = "SPIRE PMW"
        elif "plw" in filterid or "500" in filterid: final_filter_name = "SPIRE PLW"
        else:

            if channel is not None:

                if channel == 1: final_filter_name = "SPIRE PSW"
                elif channel == 2: final_filter_name = "SPIRE PMW"
                elif channel == 3: final_filter_name = "SPIRE PLW"
                else: log.warning("Could not determine which SPIRE filter was used for this image")

            elif wavelength is not None:

                if wavelength == 250: final_filter_name = "SPIRE PSW"
                elif wavelength == 350: final_filter_name = "SPIRE PMW"
                elif wavelength == 500: final_filter_name = "SPIRE PLW"
                else: log.warning("Could not determine which SPIRE filter was used for this image")

    # -- H alpha --
    elif "alpha" in filterid or "6561" in filterid or "656_1" in filterid: final_filter_name = "656_1"
    elif "ha" in filterid and "kpno" in filterid: final_filter_name = "Halpha"

    if final_filter_name is None:

        if wavelength is not None:

            value = wavelength.to("micron").value
            five_percent = 0.05 * value
            lower = value - five_percent
            upper = value + five_percent

            if "FILTER" in header: name = header["FILTER"].split("  / ")[0].replace(" ", "")
            else: name = None

            # Create a custom filter around the wavelength
            fltr = Filter((lower, upper), name=name)

        else: fltr = None

    else:

        # Create the filter
        fltr = Filter.from_string(final_filter_name)

        # Inform the user
        log.debug("Filter was identified as " + str(fltr))

    # Create and return a Filter object
    return fltr

# -----------------------------------------------------------------

def get_unit(header):

    """
    This function ...
    :param header:
    :return:
    """

    unit = None

    for keyword in ("BUNIT", "SIGUNIT", "ZUNITS"):

        if keyword in header:

            value = header[keyword].split("   / ")[0].rstrip()

            if value.isupper(): value = value.replace("DN", "count").replace("SEC", "second").lower()
            else: value = value.replace("DN", "count")

            try:
                #print(value)
                unit = u.Unit(value)
                #print(unit)
                break
            except ValueError: continue

    # Return the unit
    return unit

# -----------------------------------------------------------------

def get_fwhm(header):

    """
    This function ...
    :param header:
    :return:
    """

    fwhm = None

    for keyword in ["FWHM"]:

        if keyword in header:

            fwhm = get_quantity(header["FWHM"], default_unit="arcsec")

    # Return the FWHM
    return fwhm

# -----------------------------------------------------------------

def get_zero_point(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Loop over all keys in the header
    for key in header:

        if "MAGZP" in key: return header[key]

    # If no keyword is found that states the zero-point, return None
    return None

# -----------------------------------------------------------------

def is_sky_subtracted(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Initially, set the boolean to False
    subtracted = False

    if 'BACK_SUB' in header: subtracted = header['BACK_SUB']

    # Return the boolean value
    return subtracted

# -----------------------------------------------------------------

def get_number_of_frames(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Initially, set the boolean to False
    nframes = 1

    if 'NAXIS' in header:

        # If there are 3 axes, get the size of the third
        if header['NAXIS'] == 3: nframes = header['NAXIS3']

    # Return the boolean value
    return nframes

# -----------------------------------------------------------------

def get_frame_name_and_description(header, i, always_call_first_primary=True):

    """
    This function ...
    :param header:
    :param i:
    :param always_call_first_primary:
    :return:
    """

    planeX = "PLANE" + str(i)

    # Return the description
    if planeX in header: description = header[planeX]
    else: description = None

    plane_type = "frame"

    # FITS file created by AstroMagic
    if description is not None and "[" in description and "]" in description:

        name = description.split(" [")[0]
        plane_type = description.split("[")[1].split("]")[0]

    elif i == 0 and always_call_first_primary:

        # Get the name of this frame, but the first frame always gets the name 'primary' unless the
        # 'always_call_first_primary' flag is disabled

        description = "the primary signal map"
        name = "primary"

    elif description is not None:

        # Convert spaces to underscores and ignore things between parentheses
        name = description.split("(")[0].rstrip(" ").replace(" ", "_")

        # If the frame name contains 'error', use the standard name "errors" for this frame
        if 'error' in name: name = "errors"

    else: ## description is None

        description = ""
        name = "frame"+str(i)

    # Return the name and description
    return name, description, plane_type

# -----------------------------------------------------------------

def get_frame_description(header, i):

    """
    This function ...
    :param header:
    :param i:
    :return:
    """

    planeX = "PLANE" + str(i)

    # Return the description
    if planeX in header: return header[planeX]
    else: return None

# -----------------------------------------------------------------

def get_frame_index(header, name):

    """
    This function ...
    """

    for key in header:

        # Skip keys not ...
        if not "PLANE" in key: continue

        if header[key] == name: return int(key.split("PLANE")[1])

    return None

# -----------------------------------------------------------------

def get_frame_name(description):

    """
    This function ...
    :param description:
    :return:
    """

    # Convert spaces to underscores and ignore things between parentheses
    name = description.split("(")[0].rstrip(" ").replace(" ", "_")

    # If the frame name contains 'error', use the standard name "errors" for this frame
    if 'error' in name: name = "errors"

    # Return the frame name
    return name

# -----------------------------------------------------------------

def check_header_matches_image(image, header):

    """
    This function ...
    :param image:
    :param header:
    :return:
    """

    wcs = CoordinateSystem(header)

    # wcs.naxis attributes are deprecated, so we perform this check conditionally
    if ((hasattr(wcs,'naxis1') and hasattr(wcs,'naxis2')) and not
            (wcs.naxis1 == image.shape[1] and wcs.naxis2 == image.shape[0])):
        raise Exception("Image shape must match header shape.")

# -----------------------------------------------------------------

def get_pixel_mapping(header1, header2):

    """
    This function determines the mapping from pixel coordinates in header1 to pixel coordinates in header2 (the
    reference header). It takes the following arguments:
    :param header1:
    :param header2:
    :return: a NumPy array describing a grid of y,x pixel locations in the input header's pixel units but the output
    header's world units. It raises a TypeError if neither header is not a Header or WCS instance, and a
    NotImplementedError if the CTYPE in the header is not recognized.
    """

    # Get the WCS from the two headers
    wcs1 = CoordinateSystem(header1)
    wcs2 = CoordinateSystem(header2)

    # Convert the coordinates
    if not all([w1 == w2 for w1,w2 in zip(wcs1.wcs.ctype,wcs2.wcs.ctype)]):
        allowed_coords = ('GLON', 'GLAT', 'RA', 'DEC')
        if all([(any(word in w1 for word in allowed_coords) and
                 any(word in w2 for word in allowed_coords))
                for w1,w2 in zip(wcs1.wcs.ctype,wcs2.wcs.ctype)]):
            csys1 = ctype_to_csys(wcs1.wcs)
            csys2 = ctype_to_csys(wcs2.wcs)
            convert_coordinates = True
        else:
            # do unit conversions
            raise NotImplementedError("Unit conversions between {0} and {1} have not yet been implemented.".format(wcs1.wcs.ctype,wcs2.wcs.ctype))

    else: convert_coordinates = False

    # sigh... why does numpy use matrix convention?  Makes everything so
    # much harder...
    # WCS has naxis attributes because it is loaded with
    # _load_wcs_from_header
    outshape = [wcs2.naxis2,wcs2.naxis1]

    yy2,xx2 = np.indices(outshape)

    # get the world coordinates of the output image
    lon2,lat2 = wcs2.wcs_pix2world(xx2, yy2, 0)

    # Alternative
    #x = np.arange(wcs2.naxis1)
    #y = np.arange(wcs2.naxis2)
    #X, Y = np.meshgrid(x, y)
    #lon2, lat2 = wcs2.wcs_pix2world(X, Y, 0)

    if convert_coordinates:

        # Transform the world coordinates from the output image into the coordinate
        # system of the input image
        C2 = coordinates.SkyCoord(lon2,lat2,unit=(u.deg,u.deg),frame=csys2)
        C1 = C2.transform_to(csys1)
        lon2,lat2 = C1.spherical.lon.deg,C1.spherical.lat.deg

    xx1,yy1 = wcs1.wcs_world2pix(lon2, lat2, 0)
    grid = np.array([yy1.reshape(outshape),xx1.reshape(outshape)])

    # Return the grid
    return grid

# -----------------------------------------------------------------

def ctype_to_csys(wcs):

    """
    This function ...
    :param wcs:
    :return:
    """

    ctype = wcs.ctype[0]
    if 'RA' in ctype or 'DEC' in ctype:
        if wcs.equinox == 2000 or wcs.equinox == 2000.:
            return 'fk5'
        elif wcs.equinox == 1950 or wcs.equinox == 1950.:
            return 'fk4'
        else:
            raise NotImplementedError("Non-fk4/fk5 equinoxes are not allowed")
    elif 'GLON' in ctype or 'GLAT' in ctype:
        return 'galactic'

# -----------------------------------------------------------------

def get_float(entry):

    """
    This function ...
    :param entry:
    :return:
    """

    try: return float(entry)
    except ValueError:
        value = entry.split("   / ")[0].rstrip()
        return float(value)

# -----------------------------------------------------------------

def get_quantity(entry, default_unit=None):

    """
    This function ...
    :param entry:
    :param default_unit:
    :return:
    """

    if isinstance(entry, basestring): value = entry.split("   / ")[0].rstrip()
    else: value = entry

    try:

        num_value = float(value)
        if default_unit is None: raise RuntimeError("Default unit is not provided")
        unit = u.Unit(default_unit)

    except ValueError:

        #floats = re.findall("[-+]?\d*\.\d+|\d+", value)
        #assert len(floats) == 1
        #num_value = floats[0]

        #unit_description = value.split(str(floats[0])[-1:])[1].rstrip()

        #print("unit:::", unit_description)

        #unit = u.Unit(unit_description)

        composite_unit = u.Unit(value)

        num_value = composite_unit.to("micron")
        unit = u.Unit("micron")

    return num_value * unit

# -----------------------------------------------------------------

def get_string(entry):

    """
    This function ...
    :param entry:
    :return:
    """

    value = entry.split("   / ")[0].rstrip()
    return value

# -----------------------------------------------------------------

def get_int(entry):

    """
    This function ...
    :param entry:
    :return:
    """

    try: return int(entry)
    except ValueError:
        value = entry.split("   / ")[0].rstrip()
        return int(value)

# -----------------------------------------------------------------
