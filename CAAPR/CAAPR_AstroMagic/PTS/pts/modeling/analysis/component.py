#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.component Contains the AnalysisComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools import filesystem as fs
from ...core.launch.timing import TimingTable
from ...core.launch.memory import MemoryTable

# -----------------------------------------------------------------

class AnalysisComponent(ModelingComponent):
    
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
        super(AnalysisComponent, self).__init__(config)

        # -- Attributes --

        # The path to the analysis/in directory
        self.analysis_in_path = None

        # The path to the analysis/out directory
        self.analysis_out_path = None

        # The path to the analysis/extr directory
        self.analysis_extr_path = None

        # The path to the analysis/plot directory
        self.analysis_plot_path = None

        # The path to the analysis/misc directory
        self.analysis_misc_path = None

        # The path to the analysis/attenuation directory
        self.analysis_attenuation_path = None

        # The path to the analysis/colours directory
        self.analysis_colours_path = None

        # The path to the analysis/residuals directory
        self.analysis_residuals_path = None

        # The path to the analysis/heating directory
        self.analysis_heating_path = None

        # The path to the ski file and the wavelength grid file
        self.analysis_ski_path = None
        self.analysis_wavelengths_path = None

        # The path to the timing table
        self.timing_table_path = None

        # The path to the memory table
        self.memory_table_path = None

        # The path to the reference image
        self.reference_path = None

        # The path to the analysis/scripts directory
        self.analysis_scripts_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisComponent, self).setup()

        # Set the path to the analysis/in path
        self.analysis_in_path = fs.create_directory_in(self.analysis_path, "in")

        # Set the path to the analysis/out path
        self.analysis_out_path = fs.create_directory_in(self.analysis_path, "out")

        # Set the path to the analysis/extr path
        self.analysis_extr_path = fs.create_directory_in(self.analysis_path, "extr")

        # Set the path to the analysis/plot path
        self.analysis_plot_path = fs.create_directory_in(self.analysis_path, "plot")

        # Set the path to the analysis/misc path
        self.analysis_misc_path = fs.create_directory_in(self.analysis_path, "misc")

        # Set the path to the analysis/attenuation path
        self.analysis_attenuation_path = fs.create_directory_in(self.analysis_path, "attenuation")

        # Set the path to the analysis/colours path
        self.analysis_colours_path = fs.create_directory_in(self.analysis_path, "colours")

        # Set the path to the analysis/residuals path
        self.analysis_residuals_path = fs.create_directory_in(self.analysis_path, "residuals")

        # Set the path to the analysis/heating path
        self.analysis_heating_path = fs.create_directory_in(self.analysis_path, "heating")

        # Set the path to the analysis/scripts directory
        self.analysis_scripts_path = fs.create_directory_in(self.analysis_path, "scripts")

        # Set the path to the ski file and wavelength grid file
        self.analysis_ski_path = fs.join(self.analysis_path, self.galaxy_name + ".ski")
        self.analysis_wavelengths_path = fs.join(self.analysis_in_path, "wavelengths.txt")

        # Set the path to the timing table
        self.timing_table_path = fs.join(self.analysis_path, "timing.dat")

        # Initialize the timing table
        timing_table = TimingTable(self.timing_table_path)

        # Set the path to the memory table
        self.memory_table_path = fs.join(self.analysis_path, "memory.dat")

        # Initialize the memory table
        memory_table = MemoryTable(self.memory_table_path)

        # Set the path to the reference image
        self.reference_path = fs.join(self.truncation_path, self.reference_image + ".fits")

# -----------------------------------------------------------------
