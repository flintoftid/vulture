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

#ifndef _GNUPLOT_H_
#define _GNUPLOT_H_

#include <stdbool.h>
#include <stdio.h>

#include "fdtd_types.h"

void gnuplotMesh( bool isPhysUnits , bool isExternalSurfaces );
void gnuplotProblemSize( FILE *outputFile , int ibbox[6] );
void gnuplotBoundingBox( FILE *fp , int bbox[6] );
void gnuplotBoundingBoxNodes( FILE *outputFile , int bbox[6] , int step[3] );
void gnuplotBoundingBoxArrow( FILE *outputFile , int bbox[6] , FieldComponent field );
void gnuplotArrow( FILE *outputFile , real start[3] , real end[3] , real norm[3] , int headStyle );

#endif
