#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.initializeinput Initialize the input directory for the fitting procedure.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.modeling.fitting import InputInitializer
from pts.core.tools import logging, time, filesystem, parsing

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("path", type=str, nargs='?', help="the modeling path")
parser.add_argument("--nlambda", type=int, help="the total number of wavelengths")
parser.add_argument("--lambda_minmax", type=parsing.float_tuple, help="the minimum and maximum wavelength of the total grid")
parser.add_argument("--lambda_minmax_zoom", type=parsing.float_tuple, help="the minimum and maximum wavelength of the zoomed-in grid")
parser.add_argument("--packages", type=float, help="the number of photon packages per wavelength")
parser.add_argument("--selfabsorption", action="store_true", help="enable dust self-absorption")
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument("--report", action='store_true', help='write a report file')
parser.add_argument("--config", type=str, help="the name of a configuration file")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Set the modeling path
if arguments.path is None: arguments.path = filesystem.cwd()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = filesystem.join(arguments.path, time.unique_name("initialization") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
logging.setup_log(level=level, path=logfile_path)
logging.log.info("Starting input initialization ...")

# -----------------------------------------------------------------

# Create a InputInitializer object
initializer = InputInitializer.from_arguments(arguments)

# Run the input initialization
initializer.run()

# -----------------------------------------------------------------