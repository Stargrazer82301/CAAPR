#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.relaunch Relaunch previous simulations

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.launch.synchronizer import RemoteSynchronizer
from pts.core.tools import logging, time, parsing

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("relaunch", type=parsing.simulation_ids, help="the ID's of the simulations to relaunch")
parser.add_argument("--debug", action="store_true", help="add this option to enable debug output")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Create a RemoteSynchronizer instance
synchronizer = RemoteSynchronizer.from_arguments(arguments)

# Run the synchronizer
synchronizer.run()

# -----------------------------------------------------------------
