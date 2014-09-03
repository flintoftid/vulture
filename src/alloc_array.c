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

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "message.h"

#define EXTENTSIZE 12

/* Total memory allocated. */
static unsigned long totalMemory = 0UL;

/* 
 * Allocate an array of objects of given dimension and extents.
 * Extents are zero based.
 * 
 * Example
 *
 * double ***a;
 * a = (double ***) allocArray( &bytes , sizeof( double ) , 3 , 10 , 20 , 30 )
 * 
 * is equivalent to the static array declaration
 *
 * double a[10][20][30]:
 *
 * The array can be accessed using a[i][j][j] and free'd using free( a ).
 * The array is contiguous, hence column based offsets work. 
 * The number of bytes allocated is returned in variable bytes.
*/
void *allocArray( unsigned long *bytes , size_t size , int dimension , ... )
{

  void ***retval, ***currentPointer, **destinationPointer;
  unsigned long i, j, prevLength, stepsize , inNodes , blockSize;
  va_list marker;
  unsigned long *lengthOfLayer, *extent;
  char *extentString;
  char thisExtent[EXTENTSIZE] = "";

  extentString = (char *) malloc( sizeof( char ) * EXTENTSIZE * dimension + 3 );
  extentString[0] = '\0';

  lengthOfLayer = (unsigned long *) malloc( 2 * dimension * sizeof(unsigned long) );
  if( lengthOfLayer == 0 )
    message( MSG_ERROR , 0 , "allocArray: Failed to alloc %lu -D array temporary array!\n" , dimension );

  extent = lengthOfLayer + dimension;

  va_start( marker , dimension );

  inNodes = 0;

  strncat( extentString , "(" , 1 );
  prevLength = lengthOfLayer[0] = extent[0] = va_arg( marker , int );
  snprintf( thisExtent , EXTENTSIZE , "%lu" , extent[0] );
  strncat( extentString , thisExtent , 32 );

  for( i = 1; i < dimension ; i++ )
  {
    inNodes += prevLength;
    prevLength = lengthOfLayer[i] = prevLength * ( extent[i] = va_arg( marker , int ));
    snprintf( thisExtent , EXTENTSIZE , "x%lu" , extent[i] );
    strncat( extentString , thisExtent , EXTENTSIZE );
  }

  strncat( extentString , ")" , 1 );

  va_end( marker );

  blockSize = inNodes * sizeof( void * ) + lengthOfLayer[dimension-1] * size;
  currentPointer = retval = (void***)malloc( blockSize );
  if( retval == 0 )
    message( MSG_ERROR , 0 , "  allocArray: Failed to allocate %.3lf MiB %s %lu-D array!\n" , 
             blockSize / 1024.0 / 1024.0 , extentString , dimension );
  else
#ifdef USE_MINGW_FIX
  message( MSG_DEBUG3 , 0 , "  allocArray: Allocated %.3lf MiB () %d-D array!\n" , 
            blockSize / 1024.0 / 1024.0 , dimension );
#else
    message( MSG_DEBUG3 , 0 , "  allocArray: Allocated %.3lf MiB %s %d-D array!\n" , 
             blockSize / 1024.0 / 1024.0 , extentString , dimension );
#endif

  totalMemory += blockSize;
 
  *bytes = blockSize;

  destinationPointer = (void **)retval + extent[0];

  if( dimension > 1 )
  {

    for( i = 0 , dimension -=2 ; i < dimension ; i++ )
    {
      for( j = 0 ; j < lengthOfLayer[i] ; j++ ) 
      {
	*currentPointer++ = destinationPointer;
	destinationPointer += extent[i+1];  
      }
    }

    stepsize = extent[i+1] * size;

    for( j = 0 ; j < lengthOfLayer[i] ; j++ )
    {
      *currentPointer++ = destinationPointer;
      destinationPointer = (void **)((unsigned long)destinationPointer + stepsize );
    }

  }

  free( (void *)lengthOfLayer );
  free( extentString );

  return (void *)retval;

}

/* Deallocate array. */
void deallocArray( void *array , int dimension , ... )
{

  free( array );

  return;

}

/* Report total memory allocation. */
void allocArrayReport()
{

  if( totalMemory < 1024 )
    message( MSG_LOG , 0 , "\n  Total array allocation %lu bytes\n\n" , totalMemory );
  else if( totalMemory < 1024 * 1024 )
    message( MSG_LOG , 0 , "\n  Total array allocation %.1lf kiB\n\n" , (double)totalMemory / 1024.0 );
  else if( totalMemory < 1024 * 1024 * 1024 )
    message( MSG_LOG , 0 , "\n  Total array allocation %.1lf MiB\n\n" , (double)totalMemory / 1024.0 / 1024.0 );
  else
    message( MSG_LOG , 0 , "\n  Total array allocation %.1lf GiB\n\n" , (double)totalMemory / 1024.0 / 1024.0 / 1024.0 );

  return;

}

