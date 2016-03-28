#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.fitsed Fit the observed SED to the best RT model

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.modeling.core import SEDFitter

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Create a SEDFitter object
fitter = SEDFitter.from_arguments(arguments)

# Run the fitting procedure
fitter.run()

# -----------------------------------------------------------------