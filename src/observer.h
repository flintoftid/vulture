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

#ifndef _OBSERVER_H_
#define _OBSERVER_H_

#include "fdtd_types.h"

/* Index type used for counting and iterating over observers. */
typedef unsigned int ObserverIndex;
#define MAX_OBSERVER UINT_MAX

/* 
 * Observer formats.
 * 
 * Observer formats must begin at zero, be contigous and end with OF_UNDEFINED, which
 * *is not* included in the number NUM_OBSERVER_FORMATS.
 */
#define NUM_OBSERVER_FORMATS 3

/* Observer formats. */
typedef enum {

  OF_ASCII,
  OF_BINARY,
  OF_HDF5,
  OF_UNDEFINED

} ObserverFormat;

/* 
 * Observer domains.
 * 
 * Observer domains must begin at zero, be contigous and end with OD_UNDEFINED, which
 * *is not* included in the number NUM_OBSERVER_DOMAINS.
 */
#define NUM_OBSERVER_DOMAINS 2

typedef enum {

  OD_TIME,
  OD_FREQ,
  OD_UNDEFINED

} ObserverDomain;

/* 
 * Observer quantities.
 * 
 * Observer quantities must begin at zero, be contigous and end with OQ_UNDEFINED, which
 * *is not* included in the number NUM_OBSERVER_QUANTITIES.
 */
#define NUM_OBSERVER_QUANTITIES 9

typedef enum {

  OQ_WF,            // Waveform.
  OQ_E,             // electricField.
  OQ_H,             // magneticField.
  OQ_EH,            // EM field.
  OQ_S,             // Poynting vector.
  OQ_P,             // (volume) power density
  OQ_V,             // Voltage.
  OQ_I,             // Current.
  OQ_Z,             // Impedance.
  OQ_UNDEFINED

} ObserverQuantity;

/*
 * Public method interfaces.
 */

bool parseOP( char *line );
bool parseFF( char *line );
bool parseOT( char *line );
bool parseOF( char *line );
void initObservers( void );
void deallocObservers( void );
void updateObservers( unsigned long tstepNum , real t );
void reportObservers( void );
void gnuplotObservers( void );
void gmshObservers( void );
void gmshObservers( void );
bool thereAreObservers( void );
bool thereAreObserversDomain( ObserverDomain domain );
bool thereAreObserversFormat( ObserverFormat format );

#endif
