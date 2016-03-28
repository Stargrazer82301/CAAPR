#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.preparation.plotting Contains the PreparationPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools.logging import log

# -----------------------------------------------------------------

class PreparationPlotter(ModelingComponent):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(PreparationPlotter, self).__init__(config)

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new PreparationPlotter instance
        plotter = cls()

        # ...

        # Return the plotter
        return plotter

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

    # -----------------------------------------------------------------

    def plot(self):

        """
        An example script where aplpy is used to plot FITS images with matplotlib :

        #################
        # PREPARE IMAGE #
        #################

        # Open WISE 4 file. Units are DN (datanumbers), needs conversion to Jy
        hdulist = pyfits.open('NGC3628_W4.fits')
        im  = hdulist[0].data[0:,0:] - 233.65
        hdr = hdulist[0].header

        # For WISE 4, the DN-to-Jy conversion factor is 5.2269e-5
        im_Jy_pix = im * 5.2269e-5
        # The pixel size is 1.375 arcsec
        im_MJy_sr = im_Jy_pix * 1.e-6 / (1.375/206264.806247)**2
        hdr['BUNIT'] = 'MJySr'

        # Write out
        hdu = pyfits.PrimaryHDU(im_MJy_sr,hdr)
        hdu.writeto('NGC3628_W4_MJy_sr.fits',clobber=True)

        ##############
        # PLOT IMAGE #
        ##############

        fig = pyplot.figure(figsize=(5,3.5))

        # Set standard setup (in case multiple subfigures need to be plotted)
        def standard_setup(sp):
            sp.set_frame_color('white')
            sp.set_tick_labels_font(size='8')
            sp.set_axis_labels_font(size='10')
            sp.set_xaxis_coord_type('longitude')
            sp.set_yaxis_coord_type('longitude')
            sp.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
            sp.set_tick_color('white')
            sp.recenter(x=170.07069,y=13.589137,width=0.22,height=0.13)
            #sp.set_tick_xspacing(0.25)
            #sp.set_tick_yspacing(0.25)
            sp.set_system_latex(True)
            sp.tick_labels.hide()
            sp.axis_labels.hide()
            sp.tick_labels.set_xposition('top')
            sp.axis_labels.set_xposition('top')

        # Location of the subfigure(s)
        plotloc = [[0.16,0.10,0.78,0.80]]

        # SUBFIGURE 1
        # Load fits file
        f1 = aplpy.FITSFigure('NGC3628_W4_MJy_sr.fits',figure=fig, subplot=plotloc[0])
        # Apply standard setup for NGC 3628
        standard_setup(f1)
        # Apply other settings
        f1.show_colorscale(vmax=100, vmin=-1, cmap='hot', stretch='arcsinh')
        #f1.show_beam(major=18./3600., minor=18./3600., angle=0, fill=True, color='green')
        f1.tick_labels.show()
        f1.axis_labels.show()

        # Add colorbar
        f1.add_colorbar()
        f1.colorbar.set_location('bottom')
        f1.colorbar.set_pad(0.1)  # arbitrary units, default is 0.05
        f1.colorbar.set_width(0.2)  # arbitrary units, default is 0.2
        f1.colorbar.set_axis_label_text("WISE $S_{22}$ [MJy/sr]")
        f1.colorbar.set_font(size=10)
        f1.colorbar.set_ticks([0,10,20,40, 80])

        # plot and save
        fig.patch.set_facecolor('#3f3f3f')
        fig.canvas.draw()
        fig.savefig("NGC3628_W4.pdf", dpi=150)
        """

        pass

# -----------------------------------------------------------------