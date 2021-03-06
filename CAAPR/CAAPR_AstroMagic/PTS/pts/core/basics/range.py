#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.range Contains the IntegerRange, RealRange and QuantityRange classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.units import Unit

# -----------------------------------------------------------------

# TODO: what to do with the combination of inclusive=False and invert=True ??
# Define inclusive for minimum value and maximum value seperately??

# -----------------------------------------------------------------

class Range(object):
    
    """
    This class ...
    """
    
    def __init__(self, min_value, max_value, inclusive=True, invert=False):
        
        """
        The constructor ...
        :param min_value:
        :param max_value:
        :param inclusive:
        """

        self.min = min_value
        self.max = max_value
        self.inclusive = inclusive
        self.invert = invert

    # -----------------------------------------------------------------

    def linear(self, npoints):

        """
        This function ...
        """

        values = np.linspace(self.min, self.max, npoints, endpoint=self.inclusive)
        if self.invert: values = np.flipud(values)
        return values

    # -----------------------------------------------------------------

    def log(self, npoints):

        """
        This function ...
        :param npoints:
        :return:
        """

        values = np.logspace(self.min, self.max, npoints, endpoint=self.inclusive)
        if self.invert: values = np.flipud(values)
        return values

    # -----------------------------------------------------------------

    def sqrt(self, npoints):

        """
        This function ...
        :param npoints:
        :return:
        """

        width = self.max - self.min
        normalized = np.linspace(0.0, 1.0, npoints, endpoint=self.inclusive)
        values = self.min + normalized * width
        if self.invert: values = np.flipud(values)
        return values

# -----------------------------------------------------------------

class IntegerRange(Range):

    """
    This class ...
    """

    def __init__(self, min_value, max_value, inclusive=True, invert=False):

        """
        This function ...
        :param min_value:
        :param max_value:
        :param inclusive:
        :param invert:
        """

        assert isinstance(min_value, int)
        assert isinstance(max_value, int)

        # Call the constructor of the base class
        super(IntegerRange, self).__init__(min_value, max_value, inclusive=inclusive, invert=invert)

    # -----------------------------------------------------------------

    def linear(self, npoints):

        """
        This function ...
        :return:
        """

        real = super(IntegerRange, self).linear(npoints)
        integers = list(set(map(int, real)))
        return np.array(integers)

    # -----------------------------------------------------------------

    def log(self, npoints):

        """
        This function ...
        :param npoints:
        :return:
        """

        real = super(IntegerRange, self).log(npoints)
        integers = list(set(map(int, real)))
        return np.array(integers)

    # -----------------------------------------------------------------

    def sqrt(self, npoints):

        """
        This function ...
        :param npoints:
        :return:
        """

        real = super(IntegerRange, self).sqrt(npoints)
        integers = list(set(map(int, real)))
        return np.array(integers)

# -----------------------------------------------------------------

class RealRange(Range):

    """
    This class ...
    """

    def __init__(self, min_value, max_value, inclusive=True, invert=False):

        """
        The constructor ...
        :param min_value:
        :param max_value:
        :param inclusive:
        """

        assert isinstance(min_value, float)
        assert isinstance(max_value, float)

        # Call the constructor of the base class
        super(RealRange, self).__init__(min_value, max_value, inclusive=inclusive, invert=invert)

# -----------------------------------------------------------------

class QuantityRange(Range):

    """
    This function ...
    """

    def __init__(self, min_value, max_value, unit=None, inclusive=True, invert=False):

        """
        This function ...
        :param min_value:
        :param max_value:
        :param unit:
        :param inclusive:
        :param invert:
        """

        # Convert everything so that min_value and max_value are floats in the same unit, and so that 'unit' is the corresponding Unit
        min_is_quantity = hasattr(min_value, "unit")
        max_is_quantity = hasattr(max_value, "unit")

        if min_is_quantity and max_is_quantity:

            unit = min_value.unit
            min_value = min_value.value
            max_value = max_value.to(unit).value

        elif (not min_is_quantity) and (not max_is_quantity):

            if unit is None: raise ValueError("Unit must be specified if min_value and max_value are not quantities")
            elif isinstance(unit, basestring): unit = Unit(unit)

        else: raise ValueError("min_value and max_value must be either both quantities or both floats (with unit specified seperately)")

        # Call the constructor of the base class
        super(QuantityRange, self).__init__(min_value, max_value, inclusive=inclusive, invert=invert)

        # Set the unit
        self.unit = unit

    # -----------------------------------------------------------------

    def linear(self, npoints, as_list=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :return:
        """

        real = super(QuantityRange, self).linear(npoints)

        if as_list:
            result = []
            for num in real: result.append(num * self.unit)
        else: result = real * self.unit

        # Return the result (list or quantity)
        return result

    # -----------------------------------------------------------------

    def log(self, npoints, as_list=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :return:
        """

        real = super(QuantityRange, self).log(npoints)

        if as_list:
            result = []
            for num in real: result.append(num * self.unit)
        else: result = real * self.unit

        # Return the result (list or quantity)
        return result

    # -----------------------------------------------------------------

    def sqrt(self, npoints, as_list=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :return:
        """

        real = super(QuantityRange, self).sqrt(npoints)

        if as_list:
            result = []
            for num in real: result.append(num * self.unit)
        else: result = real * self.unit

        # Return the result (list or quantity)
        return result

# -----------------------------------------------------------------

def zip_linear(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    npoints = kwargs.pop("npoints")

    temp = []
    for arg in args:
        if isinstance(arg, QuantityRange): temp.append(arg.linear(npoints, as_list=True))
        else: temp.append(arg.linear(npoints))

    # Zip
    result = zip(*temp)
    return result

# -----------------------------------------------------------------

def zip_log(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    npoints = kwargs.pop("npoints")

    temp = []
    for arg in args:
        if isinstance(arg, QuantityRange): temp.append(arg.log(npoints, as_list=True))
        else: temp.append(arg.log(npoints))

    # Zip
    result = zip(*temp)
    return result

# -----------------------------------------------------------------

def zip_sqrt(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    npoints = kwargs.pop("npoints")

    temp = []
    for arg in args:
        if isinstance(arg, QuantityRange): temp.append(arg.sqrt(npoints, as_list=True))
        else: temp.append(arg.sqrt(npoints))

    # Zip
    result = zip(*temp)
    return result

# -----------------------------------------------------------------
