#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.io Tools for IO...

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import OrderedDict

# Import astronomical modules
from astropy.io import fits

# Import the relevants PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ..basics.coordinatesystem import CoordinateSystem
from ..tools import headers
from ..basics.mask import Mask
from .frame import Frame

# -----------------------------------------------------------------

def get_frame_names(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Load the header
    header = fits.getheader(path)

    # Get the number of planes
    nplanes = headers.get_number_of_frames(header)

    # Initialize a dictionary to contain the frame names and corresponding descriptions
    frames = dict()

    # Look at the properties of each plane
    for i in range(nplanes):

        # Get name and description of plane
        name, description, plane_type = headers.get_frame_name_and_description(header, i, always_call_first_primary=False)
        if plane_type == "frame": frames[name] = description

    # Return the frames with their name and description
    return frames

# -----------------------------------------------------------------

def load_frames(path, index=None, name=None, description=None, always_call_first_primary=True, rebin_to_wcs=False, hdulist_index=0, no_filter=False):

    """
    This function ...
    :param path:
    :param index:
    :param name:
    :param description:
    :param always_call_first_primary:
    :param rebin_to_wcs:
    :param hdulist_index:
    :return:
    """

    frames = OrderedDict()
    masks = OrderedDict()
    metadata = dict()

    filename = fs.strip_extension(fs.name(path))

    # Check if the file exists
    if not fs.is_file(path): raise IOError("File '" + path + "' does not exist")

    # Show which image we are importing
    log.debug("Reading in file '" + path + "' ...")

    # Open the HDU list for the FITS file
    hdulist = fits.open(path)

    # Get the primary HDU
    hdu = hdulist[hdulist_index]

    # Get the image header
    original_header = hdu.header

    # Get flattened form of the header
    flattened_header = headers.flattened(original_header)

    # Obtain the world coordinate system
    wcs = CoordinateSystem(flattened_header)

    # Set the filter
    if no_filter: fltr = None
    else:

        # Obtain the filter for this image
        fltr = headers.get_filter(filename, original_header)

        # Inform the user on the filter
        if fltr is not None: log.debug("The filter for the '" + filename + "' image is " + str(fltr))
        else: log.warning("Could not determine the filter for the image '" + filename + "'")

    # Obtain the units of this image
    unit = headers.get_unit(original_header)

    # Obtain the FWHM of this image
    fwhm = headers.get_fwhm(original_header)

    # Get the magnitude zero-point
    zero_point = headers.get_zero_point(original_header)

    # Check whether the image is sky-subtracted
    sky_subtracted = headers.is_sky_subtracted(original_header)

    # Check whether multiple planes are present in the FITS image
    nframes = headers.get_number_of_frames(original_header)
    if nframes > 1:

        # For each frame
        for i in range(nframes):

            # If only a frame with specific index needs to be imported, skip this frame if it does not correspond
            if index is not None and i != index: continue

            # Get name and description of frame
            name, description, plane_type = headers.get_frame_name_and_description(original_header, i, always_call_first_primary)

            # The sky-subtracted flag should only be set for the primary frame
            subtracted = sky_subtracted if i == 0 else False

            # Add this frame to the frames dictionary
            if plane_type == "frame":

                # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
                frame = Frame(hdu.data[i],
                              wcs=wcs,
                              name=name,
                              description=description,
                              unit=unit,
                              zero_point=zero_point,
                              filter=fltr,
                              sky_subtracted=subtracted,
                              fwhm=fwhm)
                frames[name] = frame

            elif plane_type == "mask":

                #data, name=None, description=None
                mask = Mask(hdu.data[i], name=name, description=description)
                masks[name] = mask

            else: raise ValueError("Unrecognized type (must be frame or mask)")

    else:

        # Sometimes, the 2D frame is embedded in a 3D array with shape (1, xsize, ysize)
        if len(hdu.data.shape) == 3: hdu.data = hdu.data[0]

        if name is None: name = "primary"
        if description is None: description = "the primary signal map"

        dummy_name, dummy_description, plane_type = headers.get_frame_name_and_description(original_header, 0)

        if plane_type == "frame":

            # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
            frame = Frame(hdu.data,
                          wcs=wcs,
                          name=name,
                          description=description,
                          unit=unit,
                          zero_point=zero_point,
                          filter=fltr,
                          sky_subtracted=sky_subtracted,
                          fwhm=fwhm)
            # Add the primary image frame
            frames[name] = frame

        elif plane_type == "mask":

            mask = Mask(hdu.data, name=name, description=description)
            # Add the mask
            masks[name] = mask

        else: raise ValueError("Unrecognized type (must be frame or mask)")

    # Add meta information
    for key in original_header: metadata[key.lower()] = original_header[key]

    # Close the FITS file
    hdulist.close()

    # Frames, masks and meta data
    return frames, masks, metadata

# -----------------------------------------------------------------

def load_frame(cls, path, index=None, name=None, description=None, plane=None, hdulist_index=None, no_filter=False, fwhm=None, add_meta=True):

    """
    This function ...
    :param path:
    :param index:
    :param name:
    :param description:
    :param plane:
    :param hdulist_index:
    :param no_filter:
    :param fwhm:
    :param add_meta:
    :return:
    """

    metadata = dict()

    # Open the HDU list for the FITS file
    hdulist = fits.open(path)

    # Look for the first HDU with data
    if hdulist_index is None:

        index = 0
        while True:

            if hdulist[index].data is not None:
                hdulist_index = index
                break
            index += 1

        if hdulist_index is None: raise ValueError("The FITS file does not contain any data")

    # Get the primary HDU
    hdu = hdulist[hdulist_index]

    # Get the image header
    header = hdu.header

    # Add meta information
    if add_meta:
        for key in header: metadata[key.lower()] = header[key]

    # Check whether multiple planes are present in the FITS image
    nframes = headers.get_number_of_frames(header)

    # Remove references to a potential third axis
    flat_header = headers.flattened(header)

    # Obtain the world coordinate system from the 'flattened' header
    wcs = CoordinateSystem(flat_header)

    # Load the frames
    header_pixelscale = headers.get_pixelscale(header)  # NOTE: SOMETIMES PLAIN WRONG IN THE HEADER !!
    pixelscale = wcs.pixelscale

    # Check whether pixelscale defined in the header is correct
    if header_pixelscale is not None:

        x_isclose = np.isclose(header_pixelscale.x.to("arcsec/pix").value, pixelscale.x.to("arcsec/pix").value)
        y_isclose = np.isclose(header_pixelscale.y.to("arcsec/pix").value, pixelscale.y.to("arcsec/pix").value)

        if not (x_isclose or y_isclose):
            print("WARNING: the pixel scale defined in the header is WRONG:")
            print("           - header pixelscale: (", header_pixelscale.x.to("arcsec/pix"),
                  header_pixelscale.y.to("arcsec/pix"), ")")
            print("           - actual pixelscale: (", pixelscale.x.to("arcsec/pix"), pixelscale.y.to("arcsec/pix"),
                  ")")

    if no_filter:
        fltr = None
    else:

        # Obtain the filter for this image
        fltr = headers.get_filter(fs.name(path[:-5]), header)

    # Obtain the units of this image
    unit = headers.get_unit(header)

    # Obtain the FWHM of this image
    if fwhm is None: fwhm = headers.get_fwhm(header)

    # Get the magnitude zero-point
    zero_point = headers.get_zero_point(header)

    # Check whether the image is sky-subtracted
    sky_subtracted = headers.is_sky_subtracted(header)

    if nframes > 1:

        if plane is not None:

            for i in range(nframes):

                # Get name and description of frame
                name, description, plane_type = headers.get_frame_name_and_description(header, i,
                                                                                       always_call_first_primary=False)

                if plane == name and plane_type == "frame":
                    index = i
                    break

            # If a break is not encountered, a matching plane name is not found
            else:
                raise ValueError("Plane with name '" + plane + "' not found")

        elif index is not None:

            name, description, plane_type = headers.get_frame_name_and_description(header, index,
                                                                                   always_call_first_primary=False)

        else:  # index and plane is None

            for i in range(nframes):
                # Get name and description of frame
                name, description, plane_type = headers.get_frame_name_and_description(header, i,
                                                                                       always_call_first_primary=False)
                if name == "primary": index = i
                break

            if index is None: index = 0  # if index is still None, set it to zero (take the first plane)

        # Get the name from the file path
        if name is None: name = fs.name(path[:-5])

        # Return the frame
        return cls(hdu.data[index],
                   wcs=wcs,
                   name=name,
                   description=description,
                   unit=unit,
                   zero_point=zero_point,
                   filter=fltr,
                   sky_subtracted=sky_subtracted,
                   fwhm=fwhm, meta=metadata)

    else:

        # Sometimes, the 2D frame is embedded in a 3D array with shape (1, xsize, ysize)
        if len(hdu.data.shape) == 3: hdu.data = hdu.data[0]

        # Get the name from the file path
        if name is None: name = fs.name(path[:-5])

        # Return the frame
        return cls(hdu.data,
                   wcs=wcs,
                   name=name,
                   description=description,
                   unit=unit,
                   zero_point=zero_point,
                   filter=fltr,
                   sky_subtracted=sky_subtracted,
                   fwhm=fwhm, meta=metadata)

# -----------------------------------------------------------------

def write_frame(data, header, path):

    """
    This function ...
    :param data:
    :param header:
    :return:
    """

    # Create the HDU
    hdu = fits.PrimaryHDU(data, header)

    # Write the HDU to a FITS file
    hdu.writeto(path, clobber=True)

# -----------------------------------------------------------------

def write_datacube(datacube, header, path):

    """
    This function ...
    :param datacube:
    :param header:
    :param path:
    :return:
    """

    # Create the HDU from the data array and the header
    hdu = fits.PrimaryHDU(np.array(datacube), header)

    # Write the HDU to a FITS file
    hdu.writeto(path, clobber=True)

    # Inform the user that the file has been created
    log.debug("File " + path + " created")

# -----------------------------------------------------------------
