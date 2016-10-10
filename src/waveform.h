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

#ifndef _WAVEFORM_H_
#define _WAVEFORM_H_

#include <stdbool.h>

#include "fdtd_types.h"

/* Index type used for counting and iterating over waveforms. */
typedef unsigned int WaveformIndex;
#define MAX_WAVEFORM UINT_MAX

/* 
 * Waveform types.
 * 
 * Waveform types must begin at zero, be contigous and end with WT_UNDEFINED, which
 * *is not* included in the number NUM_WAVEFORM_TYPES.
 */

#define NUM_WAVEFORM_TYPES 10

/* Waveform types. */
typedef enum {

  WT_GAUSSIAN_PULSE = 0,
  WT_NARROW_GAUSSIAN_PULSE,
  WT_DIFFERENTIATED_GAUSSIAN_PULSE,
  WT_RICKER_WAVELET,
  WT_MODULATED_GAUSSIAN_PULSE,
  WT_COMPACT_PULSE,
  WT_DIFFERENTIATED_COMPACT_PULSE,
  WT_MODULATED_COMPACT_PULSE,
  WT_RAMPED_SINUSOID,
  WT_EXTERNAL,
  WT_UNDEFINED

} WaveformType;

/*
 * Public method interfaces.
 */

bool parseWF( char *line );
void initWaveforms( void );
void reportWaveforms( void );
void updateWaveforms( unsigned long tstepNum , real t );
void deallocWaveforms( void );
real getWaveformValue( real t , WaveformIndex waveformNumber , real delay );
WaveformIndex getNumberOfWaveforms( void );
bool isWaveform( char *name , WaveformIndex *number );
bool thereAreWaveforms( WaveformType );
char *getWaveformName( WaveformIndex waveformNumber );

#endif
