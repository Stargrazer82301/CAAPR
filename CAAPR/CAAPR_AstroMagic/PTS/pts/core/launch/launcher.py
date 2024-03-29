#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.launcher Contains the SkirtLauncher class, which can be used to launch SKIRT simulations
#  in a convenient way.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .analyser import SimulationAnalyser
from ..simulation.execute import SkirtExec
from ..simulation.arguments import SkirtArguments
from ..basics.configurable import OldConfigurable
from ..test.resources import ResourceEstimator
from ..tools import monitoring
from ..tools import filesystem as fs
from ..tools.logging import log

# -----------------------------------------------------------------

class SkirtLauncher(OldConfigurable):

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
        super(SkirtLauncher, self).__init__(config, "core")

        # -- Attributes --

        # Create the SKIRT execution context
        self.skirt = SkirtExec()

        # Create a SimulationAnalyser instance
        self.analyser = SimulationAnalyser()

        # Set the simulation instance to None initially
        self.simulation = None

        # Set the paths to None initialy
        self.base_path = None
        self.input_path = None
        self.output_path = None
        self.extr_path = None
        self.plot_path = None
        self.misc_path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new SkirtLauncher instance
        launcher = cls()

        ## Adjust the configuration settings according to the command-line arguments

        # Ski file
        launcher.config.arguments.ski_pattern = arguments.filepath

        # Simulation logging
        launcher.config.arguments.logging.brief = arguments.brief
        launcher.config.arguments.logging.verbose = arguments.verbose
        launcher.config.arguments.logging.memory = arguments.memory
        launcher.config.arguments.logging.allocation = arguments.allocation

        # Parallelization
        if arguments.parallel is not None:
            launcher.config.arguments.parallel.processes = arguments.parallel[0]
            launcher.config.arguments.parallel.threads = arguments.parallel[1]
        else:
            launcher.config.arguments.parallel.processes = None
            launcher.config.arguments.parallel.threads = None

        # Other simulation arguments
        launcher.config.arguments.emulate = arguments.emulate
        launcher.config.arguments.single = True  # For now, we only allow single simulations

        # Extraction
        launcher.config.analysis.extraction.progress = arguments.extractprogress
        launcher.config.analysis.extraction.timeline = arguments.extracttimeline
        launcher.config.analysis.extraction.memory = arguments.extractmemory

        # Plotting
        launcher.config.analysis.plotting.seds = arguments.plotseds
        launcher.config.analysis.plotting.grids = arguments.plotgrids
        launcher.config.analysis.plotting.progress = arguments.plotprogress
        launcher.config.analysis.plotting.timeline = arguments.plottimeline
        launcher.config.analysis.plotting.memory = arguments.plotmemory
        launcher.config.analysis.plotting.reference_sed = arguments.refsed

        # Miscellaneous
        launcher.config.analysis.misc.rgb = arguments.makergb
        launcher.config.analysis.misc.wave = arguments.makewave
        launcher.config.analysis.misc.fluxes = arguments.fluxes
        launcher.config.analysis.misc.images = arguments.images
        launcher.config.analysis.misc.observation_filters = arguments.filters
        launcher.config.analysis.misc.observation_instruments = arguments.instruments
        launcher.config.analysis.misc.images_wcs = arguments.wcs
        launcher.config.analysis.misc.images_unit = arguments.unit

        # Return the new launcher
        return launcher

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Set the parallelization scheme
        if not self.has_parallelization: self.set_parallelization()

        # 3. Run the simulation
        self.simulate()

        # 4. Analyse the simulation output
        self.analyse()

    # -----------------------------------------------------------------

    @property
    def has_parallelization(self):

        """
        This function ...
        :return:
        """

        # Check whether the number of processes and the number of threads are both defined
        return self.config.parameters.parallel.processes is not None and self.config.parameters.parallel.threads is not None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SkirtLauncher, self).setup()

        # Set the paths
        self.base_path = fs.directory_of(self.config.parameters.ski_pattern) if "/" in self.config.parameters.ski_pattern else fs.cwd()
        self.input_path = fs.join(self.base_path, "in")
        self.output_path = fs.join(self.base_path, "out")
        self.extr_path = fs.join(self.base_path, "extr")
        self.plot_path = fs.join(self.base_path, "plot")
        self.misc_path = fs.join(self.base_path, "misc")

        # Check if an input directory exists
        if not fs.is_directory(self.input_path): self.input_path = None

        # Set the paths for the simulation
        self.config.parameters.input_path = self.input_path
        self.config.parameters.output_path = self.output_path

        # Create the output directory if necessary
        if not fs.is_directory(self.output_path): fs.create_directory(self.output_path, recursive=True)

        # Create the extraction directory if necessary
        if self.config.extraction.progress or self.config.extraction.timeline or self.config.extraction.memory:
            if not fs.is_directory(self.extr_path): fs.create_directory(self.extr_path, recursive=True)

        # Create the plotting directory if necessary
        if self.config.plotting.seds or self.config.plotting.grids or self.config.plotting.progress \
            or self.config.plotting.timeline or self.config.plotting.memory:
            if not fs.is_directory(self.plot_path): fs.create_directory(self.plot_path, recursive=True)

        # Create the misc directory if necessary
        if self.config.misc.fluxes or self.config.misc.images:
            if not fs.is_directory(self.misc_path): fs.create_directory(self.misc_path, recursive=True)

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining the parallelization scheme by estimating the memory requirements...")

        # Create and run a ResourceEstimator instance
        estimator = ResourceEstimator()
        estimator.run(self.config.arguments.ski_pattern)

        # Calculate the maximum number of processes based on the memory requirements
        processes = int(monitoring.free_memory() / estimator.memory)

        # If there is too little free memory for the simulation, the number of processes will be smaller than one
        if processes < 1:

            # Exit with an error
            log.error("Not enough memory available to run this simulation")
            exit()

        # Calculate the maximum number of threads per process based on the current cpu load of the system
        threads = int(monitoring.free_cpus() / processes)

        # If there are too little free cpus for the amount of processes, the number of threads will be smaller than one
        if threads < 1:

            processes = int(monitoring.free_cpus())
            threads = 1

        # Set the parallelization options
        self.config.arguments.parallel.processes = processes
        self.config.arguments.parallel.threads = threads

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing the simulation...")

        # Run the simulation
        arguments = SkirtArguments.from_config(self.config.arguments)
        self.simulation = self.skirt.run(arguments, silent=True)

    # -----------------------------------------------------------------

    def analyse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the simulation output...")

        # Set simulation analysis flags
        self.simulation.set_analysis_options(self.config.analysis)

        # Set simulation analysis paths
        self.simulation.analysis.extraction.path = self.extr_path  # or self.config.analysis.extraction.path ?
        self.simulation.analysis.plotting.path = self.plot_path    # or self.config.analysis.plotting.path ?
        self.simulation.analysis.misc.path = self.misc_path        # or self.config.analysis.misc.path ?

        # Run the analyser on the simulation
        self.analyser.run(self.simulation)

# -----------------------------------------------------------------
