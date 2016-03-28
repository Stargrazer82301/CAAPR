#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.photometry Calculate the photometry of the different input images

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import the relevant PTS classes and modules
from pts.modeling import PhotoMeter
from pts.core.tools import logging, time

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("image", type=str, nargs='?', help="the name of the image for which to calculate the photometry")
parser.add_argument("path", type=str, nargs='?', help="the modeling path")
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument("--report", action='store_true', help="write a report file")
parser.add_argument("--config", type=str, help="the name of a configuration file")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Set the modeling path
if arguments.path is None: arguments.path = os.getcwd()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = os.path.join(arguments.path, time.unique_name("photometry") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
logging.setup_log(level=level, path=logfile_path)
logging.log.info("Starting photometry ...")

# -----------------------------------------------------------------

# Create a PhotoMeter object
photometer = PhotoMeter.from_arguments(arguments)

# Run the photometry
photometer.run()

# -----------------------------------------------------------------