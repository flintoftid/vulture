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

#ifndef _SOURCE_H_
#define _SOURCE_H_

#include "fdtd_types.h"

/* Index type used for counting and iterating over sources. */
typedef unsigned int SourceIndex;
#define MAX_SOURCE     UINT_MAX

/* 
 * Source types.
 * 
 * Source types must begin at zero, be contigous and end with ST_UNDEFINED, which
 * *is not* included in the number NUM_SOURCE_TYPES.
 */
 
#define NUM_SOURCE_TYPES 13
 
typedef enum {

  ST_EFIELD = 0,
  ST_HFIELD,
  ST_ELEC_CURR_DENSITY, 
  ST_MAGN_CURR_DENSITY, 
  ST_ELEC_SURF_CURR_DENSITY, 
  ST_MAGN_SURF_CURR_DENSITY, 
  ST_ELEC_CURRENT, 
  ST_MAGN_CURRENT, 
  ST_ELEC_CURRENT_MOMENT, 
  ST_MAGN_CURRENT_MOMENT, 
  ST_VOLTAGE,
  ST_THEVENIN_VOLTAGE,
  ST_NORTON_CURRENT,
  ST_UNDEFINED

} SourceType;

/*
 * Public method interfaces.
 */

bool parseEX( char *line );
void initSources( void );
void updateSourcesEfield( real timeE );
void updateSourcesHfield( real timeH );
void reportSources( void );
void deallocSources( void );
void gnuplotSources( void );
void gmshSources( void );
bool thereAreSources( SourceType );

#endif
