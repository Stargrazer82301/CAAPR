#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.data Contains the DataPlotter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import PlottingComponent
from ..data.component import DataComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.core.frame import Frame
from ...magic.plot.imagegrid import StandardImageGridPlotter
from ...core.plot.sed import SEDPlotter
from ..core.sed import ObservedSED

# -----------------------------------------------------------------

class DataPlotter(PlottingComponent, DataComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        #super(DataPlotter, self).__init__(config)  # not sure this works
        PlottingComponent.__init__(self, config)
        DataComponent.__init__(self)

        # -- Attributes --

        # The observed SED
        self.sed = None

        # The dictionary of image frames
        self.images = dict()

    # -----------------------------------------------------------------

    def run(self, features=None):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the observed SED
        self.load_sed()

        # 3. Load the images
        self.load_images()

        # 4. Plot
        self.plot()

    # -----------------------------------------------------------------

    def load_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed SED ...")

        # Load the sed
        self.sed = ObservedSED.from_file(self.observed_sed_path)

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the images ...")

        # Loop over all subdirectories of the preparation directory
        for directory_path, directory_name in fs.directories_in_path(self.prep_path, returns=["path", "name"]):

            # Debugging
            log.debug("Opening " + directory_name + " image ...")

            # Look if an initialized image file is present
            image_path = fs.join(directory_path, "initialized.fits")
            if not fs.is_file(image_path):
                log.warning("Initialized image could not be found for " + directory_name)
                continue

            # Open the prepared image frame
            frame = Frame.from_file(image_path)

            # Set the image name
            frame.name = directory_name

            # Add the image to the dictionary
            self.images[directory_name] = frame

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the observed SED
        self.plot_sed()

        # Plot the images
        self.plot_images()

    # -----------------------------------------------------------------

    def plot_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the observed SED ...")

        # Create the SED plotter
        plotter = SEDPlotter()

        # Set properties
        plotter.transparent = True

        # Add the observed SED
        plotter.add_observed_sed(self.sed, "DustPedia")

        # Determine the path to the plot file
        path = fs.join(self.plot_data_path, "sed.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the images ...")

        # Create the image plotter
        plotter = StandardImageGridPlotter()

        # Sort the image labels based on wavelength
        sorted_labels = sorted(self.images.keys(), key=lambda key: self.images[key].filter.pivotwavelength())

        # Add the images
        for label in sorted_labels: plotter.add_image(self.images[label], label)

        # Determine the path to the plot file
        path = fs.join(self.plot_data_path, "images.pdf")

        # Set the plot title
        plotter.set_title("Images")

        # Make the plot
        plotter.run(path)

# -----------------------------------------------------------------
