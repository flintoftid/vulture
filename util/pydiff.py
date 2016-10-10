#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of Vulture.
#
# Vulture finite-difference time-domain electromagnetic solver.
# Copyright (C) 2011-2016 Ian David Flintoft
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
#
# Author: Ian Flintoft <ian.flintoft@googlemail.com>
#

import numpy 
import string
import sys
import os

if __name__ == '__main__':

    """                                                         """
    """ Compare two ASCII data files for equality within        """
    """ a relative tolerance and absolute tolerance.            """
    """                                                         """
    """ Usage: pydiff <valFileName> <testFileName>              """
    """                                                         """
    """ <valFileName>  - path to validation data file.          """
    """ <testFileName> - path to test data file.                """
    """                                                         """
    """ Returns success if the data table (a and b ) in the     """
    """ files have the same extents and                         """
    """                                                         """
    """      absolute(a - b) <= (atol + rtol * absolute(b))     """
    """                                                         """
    """ where rtol = 1e-5 and atol = 1e-2.                      """
    """                                                         """ 
    """ Differences in data that are close to the numerical     """
    """ noise floor should  not cause failure.                  """
    """                                                         """

    # Path to validated data file.
    validFileName = sys.argv[1]

    # Path to test data file.
    testFileName = sys.argv[2]

    # Load data.
    validData = numpy.loadtxt( validFileName )
    validShape = numpy.shape( validData )

    testData = numpy.loadtxt( testFileName )
    testShape = numpy.shape( testData )

    if( not numpy.allclose( validData , testData , 1e-5 , 1e-2 ) ):
        if numpy.any( validShape != testShape ):
            print "Test data and validation data are different shapes!"
        else:
            if( not numpy.allclose( validData , testData , 1e-4 , 1e-1 ) ):
                print "Mismatch in values of validation and test data!"
	        sys.exit( 1 )
            else:
		print "Marginal mismatch in values of validation and test data!"
	        sys.exit( 2 )

