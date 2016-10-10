/* 
 * This file is part of Vulture.
 *
 * Vulture finite-difference time-domain electromagnetic solver.
 * Copyright (C) 2011-2016 Ian David Flintoft
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 *
 * Author: Ian Flintoft <ian.flintoft@googlemail.com>
 *
 */

#include <math.h>

#include "util.h"
#include "physical.h"


/* Return true if two real numbers are equal within relative tolerance. */
bool isEqualRel( real x , real y , real rtol )
{

  return ( fabs( x - y ) <= ( rtol * fmax( fabs( x ) , fabs( y ) ) + REAL_EPSILON ) );

}

/* Convert degrees to radians. */
real degrees2radians( real angle )
{
  return pi / 180.0 * angle;
}

/* Convert radians to degrees. */
real radians2degrees( real angle )
{
  return 180.0 / pi * angle;
}
