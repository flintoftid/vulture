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

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "waveform.h"
#include "utlist.h"
#include "uthash.h"
#include "alloc_array.h"
#include "physical.h"
#include "message.h"
#include "simulation.h"
#include "grid.h"

/* 
 * Waveform class. 
 */
                                                         
typedef struct WaveformItem_t {

  WaveformIndex number;          // Waveform number, assigned in order found.

  /* Parameters from mesh. */
  char name[TAG_SIZE];           // Waveform name from mesh.
  WaveformType type;             // Waveform type.
  real size;                     // Waveform ampltude [-].
  real delay;                    // Waveform delay [s].
  real width;                    // Waveform width [s].
  real frequency;                // Frequench for sinusoids [Hz].

  /* Derived parameters. */
  unsigned long tableSize;     // Number of samples in interpolation table.
  unsigned long lastIdx;       // Index of last table entry used.
  real **table;                  // Interpolation table for external waveform.
  char fileName[PATH_SIZE];     // Name of input file for external waveforms.

  /* UT list/hash. */

  struct WaveformItem_t *prev;
  struct WaveformItem_t *next;

  UT_hash_handle hh;             // Waveform name hash.
  UT_hash_handle hhint;          // Waveform number hash.

} WaveformItem;

/* 
 * Private data. 
 */

/* Waveform type strings. */
static char WAVEFORM_TYPE_STR[NUM_WAVEFORM_TYPES][22] = { "GAUSSIAN_PULSE", 
                                                          "NARROW_GAUSSIAN_PULSE",
                                                          "DIFF_GAUSSIAN_PULSE",
                                                          "RICKER_WAVELET",
                                                          "MOD_GAUSSIAN_PULSE",
                                                          "COMPACT_PULSE" , 
                                                          "DIFF_COMPACT_PULSE" , 
                                                          "MOD_COMPACT_PULSE" , 
                                                          "RAMPED_SINUSOID" ,
                                                          "EXTERNAL" };

/* Number of waveforms. */
static WaveformIndex numWaveform = 0;

/* Existance Flags for waveforms of each type, including undefined. */
static bool isWaveformType[NUM_WAVEFORM_TYPES+1] = { false };

/* List of waveforms. */
static WaveformItem *waveformList = NULL;

/* Hash of waveforms using name. */
static WaveformItem *waveformHash = NULL;

/* Hash of waveforms using number. */
static WaveformItem *waveformNumberHash = NULL;

/* 
 * Private method interfaces. 
 */

void addWaveform( char *name , WaveformType type , real size , real width , real delay , real frequency , char *fileName );
real gaussianPulse( real time , real width );
real differentiatedGaussianPulse( real time , real width );
real rickerWavelet( real time , real width );
real modulatedGaussianPulse( real time , real width , real frequency );
real compactPulse( real time , real width );
real differentiatedCompactPulse( real time , real width );
real modulatedCompactPulse( real time , real width , real frequency );
real rampedSinusoid( real time , real width , real frequency );
real externalWaveform( real time , unsigned long tableSize , real **table , unsigned long *pLastIdx );
unsigned long loadExternalWaveform( char *fileName , real ***table );
double evalSpline( real **table , unsigned long tableSize , real time , unsigned long *lastIdx );
void createSplines( char *fileName , real **table , unsigned long tableSize , real deriv1 , real deriv2 , real del_t );

/*
 * Method implementations.
 */

/* Add waveform to lists. */
void addWaveform( char *name , WaveformType type , real size , real width , real delay , real frequency , char *fileName )
{

  WaveformItem *item = NULL;

  if( numWaveform == MAX_WAVEFORM )
    message( MSG_ERROR , 0 , "Maximum number of waveforms exceeded!\n" );

  item = (WaveformItem *) malloc( sizeof( WaveformItem ) );
  if( !item )
    message( MSG_ERROR , 0 , "Failed to allocate waveform %s\n" , name );

  strncpy( item->name , name , TAG_SIZE );
  item->number = numWaveform;
  item->type = type;
  item->size = size;
  item->width = width;
  item->delay = delay;
  item->frequency = (real) frequency;
  strncpy( item->fileName , fileName , PATH_SIZE );

  /* Add to list and hash. */
  DL_APPEND( waveformList , item );
  HASH_ADD_STR( waveformHash , name , item );
  HASH_ADD( hhint , waveformNumberHash , number , sizeof( numWaveform ) , item );
  numWaveform++;
  isWaveformType[type] = true;
  isWaveformType[WT_UNDEFINED] = true;
  
  return;

}

/* Parse waveform. */
bool parseWF( char *line )
{

  int numScanned = 0;
  char name[TAG_SIZE] = "";
  char fileName[PATH_SIZE] = "";
  char typeStr[TAG_SIZE] = "";
  WaveformType type = WT_GAUSSIAN_PULSE;
  WaveformIndex waveformNumber;
  bool foundType = false;
  double size = -1.0;
  double delay = -1.0;
  double width = -1.0;
  double frequency = -1.0;

  numScanned = sscanf( line , "%31s %31s" , name , typeStr );
  if( numScanned < 2 )
    return false;  

  /* Check waveform is not already defined. */
  if( isWaveform( name , &waveformNumber ) )
  {
    message( MSG_LOG , 0 , "  Waveform %s already defined\n" , name );
    return false;
  }

  /* Find type. */
  for( int waveform = 0 ; waveform < NUM_WAVEFORM_TYPES ; waveform++ )
    if( strncmp( typeStr , WAVEFORM_TYPE_STR[waveform] , TAG_SIZE ) == 0 )
    {
      type = (WaveformType)waveform;      
      foundType = true;
    }

  if( !foundType )
  {
    message( MSG_LOG , 0 , "  Invalid waveform: %s\n" , type );
    return false;
  }

  /* Re-parse and validate  depending on type. */
  if( type == WT_EXTERNAL )
  {

    numScanned = sscanf( line , "%31s %31s \"%1023[^\"]\" %lf %lf" , name , typeStr , fileName ,  &size , &delay );

    if( numScanned < 3 )
    {
      message( MSG_LOG , 0 , "  Invalid waveform card:\n" );
      return false;
    } 

    if( numScanned >= 4 && size < 0.0 )
    {
      message( MSG_LOG , 0 , "  Waveform size must be positive:\n" );
      return false;
    }

    if(  numScanned >= 5 && delay < 0.0 )
      message( MSG_WARN , 0 , "  Waveform delay negative:\n" );

  }
  else
  {

    numScanned = sscanf( line , "%31s %31s %lf %lf %lf %lf" , name , typeStr , &size , &delay , &width, &frequency );

    if( numScanned >= 3 && size < 0.0 )
    {
      message( MSG_LOG , 0 , "  Waveform size must be positive:\n" );
      return false;
    }

    if(  numScanned >= 4 && delay < 0.0 )
      message( MSG_WARN , 0 , "  Waveform delay negative:\n" );

    if( numScanned >= 5 && width <= 0.0 )
    {
      message( MSG_LOG , 0 , "  Waveform width must be positive:\n" );
      return false;
    }

    if( numScanned == 6 && frequency <= 0.0 )
    {
      message( MSG_LOG , 0 , "  Waveform frequency must be positive:\n" );
      return false;
    }

  }

  addWaveform( name , type , (real)size , (real)width , (real)delay , (real)frequency , fileName );

  return true;

}

/* Initialise waveforms. */
/* Depends: initGrid */
void initWaveforms( void )
{

  real del_t;
  WaveformItem *item;

  message( MSG_LOG , 0 , "\nInitialising waveforms...\n\n" );
   
  /* Nominal time step of largest cell in mesh - used for default parameters. */
  del_t = getGridTimeStep();

 /* Set waveforms parameters. */
  DL_FOREACH( waveformList , item ) 
  {
    
    /* Set default parameters. */
    if( item->size < 0 )
       item->size = 1.0;

    item->table = NULL;
    
    switch( item->type )
    {
    case WT_GAUSSIAN_PULSE:
    case WT_DIFFERENTIATED_GAUSSIAN_PULSE:
    case WT_RICKER_WAVELET:
      if( item->width < 0 )
         item->width = 5.0 * sqrt( 2.0 ) * del_t;
      if( item->delay < 0 )
         item->delay = 40.0 * del_t;
      item->frequency = 0.0;
      break;
    case WT_NARROW_GAUSSIAN_PULSE:
      if( item->width < 0 )
         item->width = 8.0 * del_t;
      if( item->delay < 0 )
         item->delay = 12.0 * del_t;
      item->frequency = 0.0;
      break;
    case WT_MODULATED_GAUSSIAN_PULSE:
      if( item->width < 0 )
         item->width = 20.0 * sqrt( 2.0 ) * del_t;
      if( item->delay < 0 )
         item->delay = 120.0 * del_t;
      if( item->frequency < 0 )
         item->frequency = 0.05 / del_t;
      break;
    case WT_COMPACT_PULSE:
    case WT_DIFFERENTIATED_COMPACT_PULSE:
      if( item->width < 0 )
         item->width = 20.0 * del_t; 
      if( item->delay < 0 )
         item->delay = 0.0; 
      item->frequency = 0.0;
      break;
    case WT_MODULATED_COMPACT_PULSE:
      if( item->width < 0 )
         item->width = 80.0 * del_t;   
      if( item->delay < 0 )
         item->delay = 0.0;
      if( item->frequency < 0 )
         item->frequency = 0.05 / del_t;
      break;
    case WT_RAMPED_SINUSOID:
      if( item->width < 0 )
         item->width = 20.0 * del_t;   
      if( item->delay < 0 )
         item->delay = 0.0;
      if( item->frequency < 0 )
         item->frequency = 0.05 / del_t;
      break;
    case WT_EXTERNAL:
      if( item->delay < 0 )
         item->delay = 0.0;
      item->tableSize = loadExternalWaveform( item->fileName , &(item->table) );
      item->lastIdx = 0UL;
      message( MSG_LOG , 0 , "  Read %ld entries from external waveform table in file %s\n" , item->tableSize , item->fileName );
      createSplines( item->fileName , item->table , item->tableSize , 0.0 , 0.0 , del_t );
     break;
    default:
        assert( 0 );
      break;
    }

    message( MSG_DEBUG3 , 0 , "  Setting %s waveform: size=%g, delay=%g, width=%g, freq=%g\n" , WAVEFORM_TYPE_STR[item->type] ,
             item->size , item->delay , item->width , item->frequency );
	       
  }

  return;

}

/* Report waveform parameters. */
void reportWaveforms( void )
{

  WaveformItem *item;

  message( MSG_LOG , 0 , "  Number of waveforms: %lu\n" , (unsigned long) numWaveform );

  DL_FOREACH( waveformList , item ) 
  {
    message( MSG_DEBUG3 , 0 , "    Waveform #%lu: Name=%s Type=%s size=%e delay=%e width=%e frequency=%e\n" , 
             (unsigned long) item->number , item->name , WAVEFORM_TYPE_STR[item->type] , 
             item->size , item->delay , item->width , item->frequency );
  }

  return;

}

/* Return true if there are waveforms are given type. */
bool thereAreWaveforms( WaveformType type )
{
  
  return isWaveformType[type];

}

/* Update all waveforms to current grid time. */
void updateWaveforms( unsigned long tstepNum , real t )
{

  return;

}

/* Deallocate waveform data. */
void deallocWaveforms( void )
{

  WaveformItem *item , *tmp;

  message( MSG_DEBUG1 , 0 , "Deallocating waveforms...\n" );

  /* Free waveform number hash. */
  HASH_ITER( hhint , waveformNumberHash , item , tmp )
  {
    if( item->type == WT_EXTERNAL )
      deallocArray( item->table , 2 , item->tableSize , 3 ); 

    HASH_DELETE( hhint , waveformNumberHash , item );
  }
 
  /* Free waveform name hash and the waveforms. */
  HASH_ITER( hh , waveformHash , item , tmp )
  {
    HASH_DELETE( hh , waveformHash , item );
    free( item );
  }
 
  return;

}

/* Gaussian pulse. */
real gaussianPulse( real time , real width )
{

  return exp( -0.5 * ( pow( time / width , 2.0 ) ) );

}

/* Differentiated Gaussian pulse. */
real differentiatedGaussianPulse( real time , real width )
{

  return -time / width * exp( -0.5 * ( pow( time / width , 2.0 ) ) );

}

/* Double-differentiated Gaussian pulse - Ricker wavelet. */
real rickerWavelet( real time , real width )
{

  return ( 1.0 - pow( time / width , 2 ) ) * exp( -0.5 * ( pow( time / width , 2.0 ) ) );

}

/* Modulated Gaussian pulse. */
real modulatedGaussianPulse( real time , real width , real frequency )
{

  return exp( -0.5 * ( pow( time / width , 2.0 ) ) ) * sin( 2.0 * pi * frequency * time );

}

/* Compact pulse. */
real compactPulse( real time , real width )
{

  if( time <= 0.0 )
    return 0.0;
  else if( time < 2.0 * width )
    return 1.0 / 32.0 * ( 10.0 - 15.0 * cos( pi / width * time ) + 
           6.0 * cos( 2.0 * pi / width * time ) - 
           cos( 3.0 * pi / width * time ) );
  else
    return 0.0;

}

/* Differentiated compact pulse. */
real differentiatedCompactPulse( real time , real width )
{

  if( time <= 0.0 )
    return 0.0;
  else if( time < 2.0 * width )
    return 1.0 / 32.0 * ( 15.0 * sin( pi / width * time ) - 
           12.0 * sin( 2.0 * pi / width * time ) + 
           3.0 * sin( 3.0 * pi / width * time ) );
  else
    return 0.0;

}

/* Modulated compact pulse. */
real modulatedCompactPulse( real time , real width , real frequency )
{

  real envelope;

  if( time <= 0.0 )
    envelope = 0.0;
  else if( time < 2.0 * width )
    envelope = 1.0 / 32.0 * ( 10.0 - 15.0 * cos( pi / width * time ) + 
               6.0 * cos( 2.0 * pi / width * time ) - 
               cos( 3.0 * pi / width * time ) );
  else
    envelope = 0.0;

  return envelope * sin( 2.0 * pi * frequency * time );
}

/* Sinusoid ramped by rising edge of compact pulse. */
real rampedSinusoid( real time , real width , real frequency )
{
  real ramp;

  if( time <= 0.0 )
    ramp = 0.0;
  else if( time < 1.0 * width )
    ramp = 1.0 / 32.0 * ( 10.0 - 15.0 * cos( pi / width * time ) + 
           6.0 * cos( 2.0 * pi / width * time ) - 
           cos( 3.0 * pi / width * time ) );
  else
    ramp = 1.0;

  return ramp * sin( 2.0 * pi * frequency * time );

}

/* Read external waveform from file. */
unsigned long loadExternalWaveform( char *fileName , real ***table )
{

  FILE *fp;
  double tmp1 = 0.0;
  double tmp2 = 0.0;
  unsigned long tableSize = 0;
  unsigned long bytes;

  /* Open file. */
  fp = fopen( fileName , "r" );
  if( fp == NULL ) 
     message( MSG_ERROR , 0 , "  ***Error: Cannot open external waveform file %s\n" , fileName );

  /* Preview file to get number of elements. */
  while( !feof( fp ) ) 
  {
    if( fscanf( fp , "%lf %lf" , &tmp2 , &tmp2 ) == 2 )
      tableSize++;
    else
      if( feof( fp ) )
        break;
      else
        message( MSG_ERROR , 0 , "  *** Error reading from %s\n" , fileName );
  }

  rewind( fp );

  *table = allocArray( &bytes , sizeof( real )  , 2 , tableSize , 3 );

  for( int i = 0 ; i < tableSize ; i++ )
  {
    if( fscanf( fp , "%lf %lf" , &tmp1, &tmp2 ) != 2 )
      message( MSG_ERROR , 0 , "  *** Error reading from %s\n" , fileName );
    (*table)[i][0] = (real)tmp1;
    (*table)[i][1] = (real)tmp2;   
    (*table)[i][2] = 0.0;     
  }

  fclose( fp );  
 
  return tableSize;

}

/* Interpolate waveform from table. */
real externalWaveform( real time , unsigned long tableSize , real **table , unsigned long *pLastIdx )
{

  if( time < table[0][0] || time > table[tableSize-1][0] )
    return 0.0;
  else
    return evalSpline( table , tableSize , time , pLastIdx );

  return 0.0;

}

/* Get waveform number from name. */
bool isWaveform( char *name , WaveformIndex *number )
{

  WaveformItem *item;

  HASH_FIND_STR( waveformHash , name , item );
  if( item )
  {
    *number = item->number;
    return true;
  }
  else
  {
    *number = 0;
    return false;
  }
  
}

/* Get value of waveform by number. */
real getWaveformValue( real t , WaveformIndex waveformNumber , real delay )
{
  
  real value = 0.0;
  WaveformItem *item;

  /* Get waveform by number. */
  HASH_FIND( hhint , waveformNumberHash , &waveformNumber , sizeof( waveformNumber ) , item );
  if( !item )
    assert( 0 ); /* Parser has failed if this happens. */

  switch( item->type )
  {
  case WT_GAUSSIAN_PULSE:
  case WT_NARROW_GAUSSIAN_PULSE:
    value = item->size * gaussianPulse( t - delay - item->delay , item->width );
    break;
  case WT_DIFFERENTIATED_GAUSSIAN_PULSE:
    value = item->size * differentiatedGaussianPulse( t - delay - item->delay , item->width );
    break;
  case WT_RICKER_WAVELET:
    value = item->size * rickerWavelet( t - delay - item->delay , item->width );
    break;
  case WT_MODULATED_GAUSSIAN_PULSE:
    value = item->size * modulatedGaussianPulse( t - delay - item->delay , item->width , item->frequency );
    break;
  case WT_COMPACT_PULSE:
    value = item->size * compactPulse( t - delay - item->delay , item->width );
    break;
  case WT_DIFFERENTIATED_COMPACT_PULSE:
    value = item->size * differentiatedCompactPulse( t - delay - item->delay , item->width );
    break;
  case WT_MODULATED_COMPACT_PULSE:
    value = item->size * modulatedCompactPulse( t - delay - item->delay , item->width , item->frequency );
    break;
  case WT_RAMPED_SINUSOID:
    value = item->size * rampedSinusoid( t - delay - item->delay , item->width , item->frequency );
    break;
  case WT_EXTERNAL:
    value = item->size * externalWaveform( t - delay - item->delay , item->tableSize , item->table , &(item->lastIdx) );
    break;
  default:
    assert( 0 );
    break;
  }

  return value;
  
}


/* Get pointer to name of waveform by number. */
char *getWaveformName( WaveformIndex waveformNumber )
{
  
  WaveformItem *item;

  /* Get waveform by number. */
  HASH_FIND( hhint , waveformNumberHash , &waveformNumber , sizeof( waveformNumber ) , item );
  if( !item )
    assert( 0 ); /* Parser has failed if this happens. */

  return item->name;
  
}

/* Get number of waveforms. */
WaveformIndex getNumberOfWaveforms( void )
{

  return numWaveform;

}

/* Create second derivative entries in spline table. */
void createSplines( char *fileName , real **table , unsigned long tableSize , real deriv1 , real deriv2 , real del_t )
{

  unsigned long bytes;
  unsigned long i;
  unsigned long k;
  real p;
  real qnm1;
  real sign;
  real unm1;
  real *u;
  real maxDiff; 
  real thisDiff;
  real samplingRatio;

  /* Check we have enough points. */  
  if( tableSize < 2UL )
    message( MSG_ERROR , 0 , "  Insufficient points for spline evaluation in file %s.\n" , fileName );

  /* Check sampling rate is adequate. */
  maxDiff = 0.0; 
  for( i = 1UL ; i < tableSize ; i++ )
  { 
    thisDiff = table[i][0] - table[i-1][0];  
    if( thisDiff <= 0.0 )
      message( MSG_ERROR , 0 , "  Time data in %s is not monotonically increasing.\n" , fileName );
    else if( thisDiff > maxDiff )
      maxDiff = thisDiff;
  }
 
  samplingRatio = maxDiff / del_t;
 
  if( samplingRatio < 1.5 )
    message( MSG_LOG , 0 , "  External waveform in %s is well sampled (ratio %e).\n" , fileName , samplingRatio );
  else if( ( samplingRatio > 1.5 ) && ( samplingRatio < 3.0 ) )
    message( MSG_WARN , 0 , "  External waveform in %s may be undersampled (ratio %e)!\n" , fileName , samplingRatio );
  else
    message( MSG_WARN , 0 , "  *** External waveform in %s greatly undersampled (ratio %e) ***.\n" , fileName , samplingRatio );

  /* Temporary array for ? */
  u = allocArray( &bytes , sizeof( real ) , 1 , tableSize - 1UL );

  /* Low side boundary point. */
  if( deriv1 > 0.99e30 )
  {
    table[0][2] = 0.0;
    u[0] = 0.0;
  }
  else 
  {
    table[0][2] = -0.5;
    u[0] = ( 3.0 / ( table[1][0] - table[0][0] ) ) * ( ( table[1][1] - table[0][1] ) / ( table[1][0] - table[0][0] ) - deriv1 );
  }

  /* Internal points. */
  for( i = 1UL ; i <= tableSize - 2UL ; i++ ) 
  {
    sign = ( table[i][0] - table[i-1][0] ) / ( table[i+1][0] - table[i-1][0] );
    p = sign * table[i-1][2] + 2.0;
    table[i][2] = ( sign - 1.0 ) / p;
    u[i] = ( table[i+1][1] - table[i][1] ) / ( table[i+1][0] - table[i][0] ) - 
           ( table[i][1] - table[i-1][1] ) / ( table[i][0] - table[i-1][0] );
    u[i] = ( 6.0 * u[i] / ( table[i+1][0] - table[i-1][0] ) - sign * u[i-1] ) / p;
  }

  /* High side boundary point. */
  if( deriv2 > 0.99e30 )
  {
    qnm1 = 0.0;
    unm1 = 0.0;
  } 
  else 
  {
    qnm1 = 0.5;
    unm1 = ( 3.0 / ( table[tableSize-1][0] - table[tableSize-2][0] ) ) * 
           ( deriv2 - ( table[tableSize-1][1] - table[tableSize-2][1] ) / ( table[tableSize-1][0] - table[tableSize-2][0] ) );
  }

  /* Back recursion of second derivatives. */
  table[tableSize-1][2] = ( unm1 - qnm1 * u[tableSize-2] ) / ( qnm1 * table[tableSize-2][2] + 1.0 );

  /* Iterate down from (tableSize - 1) to 0 - beware unsigned type */
  /* so for( k = tableSize - 2UL ; k>= 0 ; k-- ) fails to stop! */
  for( k = tableSize - 1UL ; k-- > 0 ; )
    table[k][2] = table[k][2] * table[k+1][2] + u[k];

  deallocArray( u , 1 , tableSize - 1UL );

  return;

}

/* Evaluate splines. Return index of which interval is used. */
double evalSpline( real **table , unsigned long tableSize , real time , unsigned long *lastIdx )
{

  unsigned long k;
  unsigned long klo;
  unsigned long khi;
  real h;
  real b;
  real a;

  /* Set inteval to last one used. */
  klo = *lastIdx;
  khi = klo + 1UL;
  if( klo >= tableSize - 1UL ) klo = 0;
  if( klo < 0 )  klo = 0;

  /* If not still in last interval perform binary search for required interval. */
  if( ( table[klo][0] > time ) || ( table[klo+1][0] < time ) )
  {
    klo = 0;
    khi = tableSize - 1UL;

    while( khi - klo > 1UL ) 
    {
      k = ( khi + klo ) >> 1;
      if( table[k][0] > time ) 
        khi = k;
      else 
        klo = k;
    }
  }

  /* Update last interval used. */
  *lastIdx = klo;

  /* Evaluate spline. */
  h = table[khi][0] - table[klo][0];
	
  if( h <= 0.0 ) 
    message( MSG_ERROR , 0 , "\n  Bisection failure in spline evaluation.\n\n" );

  a = ( table[khi][0] - time ) / h;
  b = ( time - table[klo][0] ) / h;

  return a * table[klo][1] + b * table[khi][1] +( ( a * a * a - a ) * 
          table[klo][2] + ( b * b * b - b ) * table[khi][2] ) * ( h * h ) / 6.0;

}

