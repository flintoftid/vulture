/* 
 * This file is part of Vulture.
 *
 * Vulture finite-difference time-domain electromagnetic solver.
 * Copyright (C) 2011-2013 Ian David Flintoft
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
 * Author: Ian Flintoft <ian.flintoft@york.ac.uk>
 *
 */

#include <stdio.h>
#include <math.h>

#include "memory.h"
#include "message.h"

/* Global data. */
Memory memory = { 0UL , 0UL , 0UL , 0UL , 0UL , 0UL , 0UL , 0UL , 0UL , 0UL , 0UL , 0UL };

/* Private functions. */
double normMemory( unsigned long numBytes , int n );


/* Report memory usage. */
void reportMemory( void  )
{

  message( MSG_LOG , 0 , "\nMemory usage:\n\n" );

  message( MSG_LOG , 0 , "  Grid auxiliary arrays:          %g kiB\n" , normMemory( memory.grid , 1 ) );
  message( MSG_LOG , 0 , "  Grid E/H field arrays:          %g MiB\n" , normMemory( memory.ehFields , 2 ) );
  message( MSG_LOG , 0 , "  Grid update coefficient arrays: %g MiB\n" , normMemory( memory.ehCoeffs , 2 ) );
  message( MSG_LOG , 0 , "  PML field arrays:               %g MiB\n" , normMemory( memory.pmlFields , 2 ) );
  message( MSG_LOG , 0 , "  PML coefficients arrays:        %g kiB\n" , normMemory( memory.pmlCoeffs , 1 ) );
  message( MSG_LOG , 0 , "  Waveforms:                      %g kiB\n" , normMemory( memory.waveforms , 1 ) );
  message( MSG_LOG , 0 , "  Sources:                        %g kiB\n" , normMemory( memory.sources , 1 ) );
  message( MSG_LOG , 0 , "  Observers:                      %g kiB\n" , normMemory( memory.observers , 1 ) );
  message( MSG_LOG , 0 , "  Boundaries:                     %g kiB\n" , normMemory( memory.boundaries , 1 ) );
  message( MSG_LOG , 0 , "  Surfaces:                       %g kiB\n" , normMemory( memory.surfaces , 1 ) );
  message( MSG_LOG , 0 , "  Media:                          %g kiB\n" , normMemory( memory.media , 1 ) );
  message( MSG_LOG , 0 , "  Blocks:                         %g kiB\n" , normMemory( memory.blocks , 1 ) );
  message( MSG_LOG , 0 , "  Wires:                          %g kiB\n" , normMemory( memory.wires , 1 ) );
  message( MSG_LOG , 0 , "  Lines:                          %g kiB\n" , normMemory( memory.lines , 1 ) );

  return;

}

/* Normalise memory size. */
double normMemory( unsigned long numBytes , int n )
{

  return (double) numBytes / pow( 1024 , n );

}
