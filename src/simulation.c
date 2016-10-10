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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "simulation.h"
#include "message.h"
#include "waveform.h"
#include "source.h"
#include "planewave.h"
#include "observer.h"
#include "block.h"
#include "line.h"
#include "surface.h"
#include "mur.h"  
#include "grid.h"
#include "timer.h"

/* 
 * Private data.
 */

/* Number of time steps. */
static unsigned long numTimeSteps;    

/* Courant stability factor. */
static real courantNumber = -1.0;

/* 
 * Private method interfaces. 
 */

void setCourantNumber( real counrantNum );

/*
 * Method Implementations.
 */

/* Parse Courant number. */
bool parseCN( char *line )
{

  double courantNum;

  if( sscanf( line , "%lf" , &courantNum ) != 1 )
    return false;  

  if( courantNum >= 1.0 )
  {
    message( MSG_LOG , 0 , "  Courant number must less than 1!\n" );
    return false;
  }
  else
  {
    setCourantNumber( (real) courantNum );
  }

  return true;

}

/* Parse number of time-steps. */
bool parseNT( char *line )
{

  unsigned long numSteps = 0UL;

  if( sscanf( line , "%lu" , &numSteps ) != 1 )
    return false;  

  if( numSteps < 0UL )
  {
    message( MSG_LOG , 0 , "  Number of time steps must be >= 0!\n" );
    return false;
  }
  else
  {
    setNumTimeSteps( numSteps );
  }

  return true;

}

/* Initialise simulation. */
/* Depends: */
void initSimulation( void )
{

  message( MSG_LOG , 0 , "\nInitialising simulation...\n\n" );

  /* Set default Courant stability factor. */
  if( courantNumber < 0.0 )
    courantNumber = sqrt( 3.0 ) / 2.0;

  return;

}

/* Report simulation. */
void reportSimulation()
{

  message( MSG_LOG , 0 , "  Number of time steps: %lu\n" , numTimeSteps );
  message( MSG_LOG , 0 , "  Courant number: %g\n" , courantNumber );

  return;

}

/* Propgate fields. */
void propagate( void )
{

  /* Time step counter. */
  unsigned long timeStepNumber = 0UL;

  /* Time step. */
  real dt = 0.0;
    
  /* Current physical time for electric and magnetic fields. */
  real timeE = 0.0;
  real timeH = 0.0;

  dt = getGridTimeStep();
  
  /* Time loop. */
  message( MSG_LOG , 0 , "\nStarting time stepping loop...\n" );
	
  startTimer( 0 , numTimeSteps );

  for ( timeStepNumber = 0 ; timeStepNumber <= numTimeSteps - 1 ; timeStepNumber++ )  {

    /* Electric field time. */
    timeE = timeStepNumber * dt;

    /* Magnetic field time. */
    timeH = ( timeStepNumber + 0.5 ) * dt;
    
    updateTimer( timeStepNumber , numTimeSteps );

    /* Update observers. */
    updateObservers( timeStepNumber , timeE );

    /* Update waveforms - currently a no-op. */
    updateWaveforms( timeStepNumber , timeE );

    /* Update external surface E fields. */
    updateExternalSurfacesEfield();

    /* Update main grid E fields. */
    updateGridEfield();
    
    /* Update block E fields. */
    updateBlocksEfield(); 

    /* Update internal surface H field. */
    updateInternalSurfacesEfield();

    /* Thin wire E field update goes here. */
    updateLinesEfield();

#ifndef CHECK_LIMITS
    /* Update electric field sources. */
    updateSourcesEfield( timeE );

    /* TFSF E field update. */
    updatePlaneWavesEfield( timeE );    
#endif

    updateGhostEfield();
    
    /* Update external surface H fields. */
    updateExternalSurfacesHfield();
    
    /* Update main grid H fields. */
    updateGridHfield();

    /* Update block H fields. */
    updateBlocksHfield();
    
    /* Update internal surface H field. */
    updateInternalSurfacesHfield();

    /* Thin wire H field update goes here. */
    updateLinesHfield();

#ifndef CHECK_LIMITS    
    /* Update magnetic field sources. */
    updateSourcesHfield( timeH );

    /* TFSF H field update. */
    updatePlaneWavesHfield( timeH );    
#endif
    
    updateGhostHfield();
    
  } /* for */
  
  stopTimer( numTimeSteps , numTimeSteps );

  message( MSG_LOG , 0 , "\nCompleted time stepping loop.\n\n" );

  return;

}

/* Deallocate simulation. */
void deallocSimulation( void )
{

  message( MSG_DEBUG1 , 0 , "Dealloctaing simulation...\n" );

  return;

}

/* Get the Courant number. */
real getCourantNumber( void )
{

  return courantNumber;

}

/* Set the Courant number. */
void setCourantNumber( real counrantNum )
{

  courantNumber = counrantNum;

  return;

}

/* Get number of time steps. */
unsigned long getNumTimeSteps( void )
{

  return numTimeSteps;

}

/* Set number of time steps. */
void setNumTimeSteps( unsigned long numSteps )
{

  numTimeSteps = numSteps;

  return;

}
