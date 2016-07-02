#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.scaling Test the scaling of SKIRT on a particular system.
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.test.scaling import ScalingTest
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration()

# Required arguments
config.add_required("ski_path", "absolute_path", "the name of the ski file to be used for the scaling test")
config.add_required("remote", str, "the name of the remote host")
config.add_required("mode", str, "the parallelization mode for the scaling test", choices=["mpi", "hybrid", "threads"])

# Optional arguments
config.add_positional_optional("maxnodes", float, "the maximum number of nodes", 1)
config.add_positional_optional("minnodes", float, "the minimum number of nodes. In hybrid mode, this also defines the number of threads per process", 0)
config.add_optional("cluster", str, "the name of the cluster", None)
config.add_optional("wavelengths", float, "boost the number of wavelengths by a certain factor", None)
config.add_optional("packages", float, "boost the number of photon packages by a certain factor", None)

# Flags
config.add_flag("manual", "launch and inspect job scripts manually")
config.add_flag("keep", "keep the output generated by the different SKIRT simulations")

# Read the configuration settings from the provided command-line arguments
config.read()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), time.unique_name("scaling") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting scaling ...")

# -----------------------------------------------------------------

# Create a ScalingTest instance
test = ScalingTest(config.get_settings())

# Run the scaling test
test.run()

# -----------------------------------------------------------------
