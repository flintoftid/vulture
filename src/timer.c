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

/* 
 * Determine time of completion for iterative process.
 *
 * Refactored from eta.c in Falcon source code by
 * Stuart Porter and John Dawson.
 *
 */

#include <stdio.h>
#include <time.h>
#include <math.h>

#include "timer.h"
#include "message.h"
#include "grid.h"

#define MIN_WAIT_TIME      10       // Minimum time in seconds to wait before calculating time per iteration.
#define MIN_WAIT_STEPS     10       // Minimum number of steps to wait before calculating time per iteration.
#define NUM_FAST_SAMPLES   10       // Number of samples taken using fast filter.
#define TPI_STABLE_TOL   0.01       // Proportion change of FTI allowed between samples for 'stable' FTI.
#define END_TIME_TOL     0.01       // Proportion change in end time before new time printed.
#define FAST_COEFF        0.5       // Fast filter coefficient.
#define SLOW_COEFF        0.9       // Slow filter coefficient.

/* Private data. */
static time_t startTime;            // Start time.
static time_t lastTime;             // Last time-step time. 
static time_t nowTime;              // Current time-step time.
static time_t estEndTime;           // Current estimated end time.
static time_t lastEstEndTime;       // Last estimated end time.
static double timePerIteration;     // Current time-per-iteration.
static double lastTimePerIteration; // Last time-per-iteration.
static unsigned long lastTimeStep;  // Last time-step number.
static unsigned long numSamples;    // Number of samples taken.

/* Private functions.*/
double timeFilter( double coeff , double timePerIteration , double timeDiff , unsigned long numStep );


/* Initialise time. */
void startTimer( unsigned long timeStep , unsigned long numTimeSteps )
{

  timePerIteration = lastTimePerIteration = 0.0;

  startTime = nowTime = lastTime = lastEstEndTime = time( NULL );

  numSamples = lastTimeStep = 0;

  message( MSG_LOG , 0 , "\n  %lu Iterations - Start time: %s" , numTimeSteps , ctime ( &startTime ) );

  return;

}

/* Stop timer. */
void stopTimer( unsigned long timeStep , unsigned long numTimeSteps )
{

  double averageTimePerIteration;
  unsigned long numCells;
  int innerBox[6];
  int outerBox[6];

  nowTime = time( NULL );

  averageTimePerIteration = difftime( nowTime , startTime ) / (double) numTimeSteps;

  getGridBoundingBox( innerBox , outerBox );

  numCells = ( outerBox[XHI] - outerBox[XLO] ) * ( outerBox[YHI] - outerBox[YLO] ) * ( outerBox[ZHI] - outerBox[ZLO] );

  message( MSG_LOG , 0 , "\n  %lu/%lu Iterations - Completed: %s" , timeStep , numTimeSteps , ctime( &nowTime ) );

  message( MSG_LOG , 0 , "\n  Average spi %g, average spi/cell %g ns\n" , averageTimePerIteration , averageTimePerIteration / (double) numCells / 1e-9 );

  return;

}

/* Update timer. */
void updateTimer( unsigned long timeStep , unsigned long numTimeSteps )
{

  double timeDiff;  

  /* Time now. */
  nowTime = time( NULL );

  /* Time difference between last and current sample. */
  timeDiff = difftime( nowTime , lastTime );

  /* Wait at least MIN_WAIT_TIME before estimating iteration period. */
  if( timeDiff >= MIN_WAIT_TIME && timeStep > MIN_WAIT_STEPS )
  {

    /* Count samples. */
    numSamples++; 

    /* Calculate time per iteration. */
    if( numSamples == 1 )
    { 
      /* Prime filter with first sample. */
      timePerIteration = timeDiff / (double) ( timeStep - lastTimeStep ); 
    } 
    else
    { 
      /* Update filter and calculate completion time. */
      if( numSamples <= NUM_FAST_SAMPLES )
      {	
        /* Fast filter first NUM_FAST_SAMPLES samples. */
        timePerIteration = timeFilter( FAST_COEFF , timePerIteration , timeDiff , timeStep - lastTimeStep );
      } 
      else
      {	
        /* Slow filter after NUM_FAST_SAMPLES samples. */
       timePerIteration = timeFilter( SLOW_COEFF , timePerIteration , timeDiff , timeStep - lastTimeStep );
      }
	
      /* Estimate completion time if time-per-iteration stable. */
      if( TPI_STABLE_TOL > fabs( ( timePerIteration - lastTimePerIteration ) / timePerIteration ) )
      {	
        /* Stable. */
        estEndTime = nowTime + (time_t) ( timePerIteration * (float) ( numTimeSteps - timeStep ) );

        /* Print completion time if it differs from last completion time by more than END_TIME_TOL seconds. */
        if( END_TIME_TOL < fabs( (double) ( estEndTime - lastEstEndTime ) / (double) ( estEndTime - startTime ) ) )
        {
          message( MSG_LOG , 0 , "\n  %g spi, %lu/%lu iterations, Est. end: %s" , timePerIteration , timeStep , numTimeSteps , ctime( &estEndTime ) );
          lastEstEndTime = estEndTime;
        }
      }
    }

    /* Save current values. */
    lastTime = nowTime;
    lastTimePerIteration = timePerIteration;
    lastTimeStep = timeStep;

  }

  return;

}

/* Filter. */
double timeFilter( double coeff , double timePerIteration , double timeDiff , unsigned long numStep )
{

 return ( timePerIteration * coeff ) + ( 1.0 - coeff ) * ( timeDiff / (double) numStep );

}

