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

#ifndef _PLANEWAVE_H_
#define _PLANEWAVE_H_

#include "fdtd_types.h"

/* Index type used for counting and iterating over plane waves. */
typedef unsigned int PlaneWaveIndex;
#define MAX_PLANE_WAVE UINT_MAX

/*
 * Public method interfaces.
 */

bool parsePW( char *line );
void initPlaneWaves( void );
void updatePlaneWavesEfield( real timeE );
void updatePlaneWavesHfield( real timeH );
void reportPlaneWaves( void );
void deallocPlaneWaves( void );
void gnuplotPlaneWaves( void );
void gmshPlaneWaves( void );
bool thereArePlaneWaves( void );

#endif
