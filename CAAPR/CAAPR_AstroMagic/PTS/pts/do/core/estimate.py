#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.estimate Estimate the resource requirements for a certain ski file

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import the relevant PTS classes and modules
from pts.core.test.resources import ResourceEstimator
from pts.core.tools.logging import log

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("skifile", type=str, help="the name of the ski file")
parser.add_argument("threads", type=int, help="the number of parallel threads")
parser.add_argument("processes", type=int, help="the number of parallel processes")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Initialize the logger
log.info("Starting estimate script ...")

# -----------------------------------------------------------------

# Determine the full path to the parameter file
ski_path = os.path.abspath(arguments.file)

# Create and run a ResourceEstimator oject
estimator = ResourceEstimator()
estimator.run(ski_path, arguments.processes, arguments.threads)

# Inform the user about the resource requirements
log.info("This simulation requires " + estimator.memory + " GB of virtual memory")
#log.info("This simulation requires a walltime of " + estimator.walltime)

# -----------------------------------------------------------------
