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
#include <string.h>
#include <math.h>

#include "observer.h"
#include "utlist.h"
#include "alloc_array.h"
#include "message.h"
#include "bounding_box.h"
#include "simulation.h"
#include "grid.h"
#include "physical.h"
#include "gnuplot.h"
#include "gmsh.h"
#include "waveform.h"
#include "mesh.h"
#include "memory.h"
  
/* 
 * Observer class. 
 */

#define CACHE_SIZE 5

/* Maximum number of components in an observer. */
#define MAX_COMP 6     

typedef struct ObserverItem_t {

  ObserverIndex number;             // Observer number, assigned in order found.

  /* Parameters from mesh. */
  char name[TAG_SIZE];              // Observer name from mesh.
  ObserverFormat format;            // Observer format.
  ObserverDomain domain;            // Observer domain.
  ObserverQuantity quantity;        // Observer quantity.        
  int mbbox[6];                     // Bounding box on mesh.
  int step[3];                      // Step along each axis.
  WaveformIndex waveformNumber;     // Number of reference waveform.
  int cacheSize;                    // Length of cache. 
  bool isInternal;                  // Boolean indicating if observer is private.
  
  /* Derived parameters. */
  int numComp;                            // Number of components.
  int gbbox[6];                           // Bounding box on grid.
  FILE *outputFile;                       // Output file for ACSII types.
  real **dft_real;                        // Real part of running DFT - 1D.
  real **dft_imag;                        // Imaginary part of running DFT - 1D.
  struct ObserverItem_t *waveformObserver;// Pointer to reference waveform for DFT types.
  //real *****var_real;                   // Cache/DFT array. var_real[ii][jj][kk][comp][t/f][ii][jj][kk]
  //real *****var_imag;                   // Cache/DFT array. var_imag[ii][jj][kk][comp][1/f][ii][jj][kk]
  //struct ObserverItem_t *subObs;        // Array of sub-observers.
  //char dataPath[PATH_SIZE];             // Path name to data in HDF file.
  
  /* UT list. */

  struct ObserverItem_t *prev;
  struct ObserverItem_t *next;
  
  UT_hash_handle hh;                // Observer name hash.

} ObserverItem;

/* 
 * Private data. 
 */

/* Observer format strings. */
static char OBSERVER_FORMAT_STR[NUM_OBSERVER_FORMATS][7] = { "ASCII" , "BINARY" , "HDF5" };

/* Observer domain strings. */
static char OBSERVER_DOMAIN_STR[NUM_OBSERVER_DOMAINS][5] = { "TIME" , "FREQ" };

/* Observer quantity strings. */
static char OBSERVER_QUANTITY_STR[NUM_OBSERVER_QUANTITIES][10] = { "WAVEFORM" , "EFIELD" , "HFIELD" ,
                                                                   "EHFIELD" , "POYNTING" , "POWDEN" ,
                                                                   "VOLTAGE" , "CURRENT" , "IMPEDANCE" };

/* Observer quantity number of components map. */                                                                   
static int observerCompMap[NUM_OBSERVER_QUANTITIES] = { 1 , 3 , 3 ,
                                                        6 , 3 , 1 ,
                                                        1 , 1 , 1 };
                                                        
/* Number of observers. */
static ObserverIndex numObserver = 0;

/* Number of observers of each format. */
static ObserverIndex numObserverFormat[NUM_OBSERVER_FORMATS] = { 0 };

/* Existance flag for observers of each format, including undefined. */
static bool isObserverFormat[NUM_OBSERVER_FORMATS+1] = { false };

/* Number of observers of each domain. */
static ObserverIndex numObserverDomain[NUM_OBSERVER_DOMAINS] = { 0 };

/* Existance flag for observers of each domain, including undefined. */
static bool isObserverDomain[NUM_OBSERVER_DOMAINS+1] = { false };

/* List of observers. */
static ObserverItem *observerList = NULL;

/* Hash of observers using name. */
static ObserverItem *observerHash = NULL;

/* Existance, start and stop time-steps and times and number of output time-steps. */
static bool isOT = false;
static unsigned long startTimeStep = 0;
static unsigned long stopTimeStep = 0;
static real startTime = -1.0;
static real stopTime = -1.0;
static unsigned long numOutTimeSteps = 0;

/* Existance, start, stop, step and number of frequencies for DFTs. */
static bool isOF = false;
static real startFreq = -1.0;
static real stopFreq = -1.0;
static real stepFreq = -1.0;
static unsigned long numFreq = 0; 

/* Array of angular frequenceis for DFTs. */
static real *omega = NULL;

/* 
 * Private method interfaces. 
 */

ObserverItem *addObserver( int mbbox[6] , int step[3] , char name[TAG_SIZE] , ObserverFormat format , ObserverDomain domain , 
                           ObserverQuantity quantity , unsigned long cacheSize , bool isInternal , WaveformIndex waveformNumber );
void getNumberOfOutputNodes( int *nx , int *ny , int *nz , int bbox[6] , int step[3] );
bool isObserver( char *name , ObserverIndex *number );
void initObserverAsciiTime( ObserverItem *item );
void initObserverAsciiFreq( ObserverItem *item );
void updateObserverAsciiTime( ObserverItem *item , unsigned long tstepNum , real t );
void updateObserverAsciiFreq( ObserverItem *item , unsigned long tstepNum , real t );
void deallocObserverAsciiTime( ObserverItem *item );
void deallocObserverAsciiFreq( ObserverItem *item );
void initObserverDft( ObserverItem *item );
void flushObserverDft( ObserverItem *item );
void deallocObserverDft( ObserverItem *item );
void initBinaryObservers( real dt );
void initExciteDat( void );
void writeProcessDat( void );
void initImpulseDat( real dt );
void updateExciteDat( real value );
void updateImpulseDat(  ObserverItem *item );
void deallocBinaryObservers( void );

/*
 * Method Implementations.
 */

/* Parse observers. */
bool parseOP( char *line )
{
  char OBSERVER_TYPE_STR[3][24]   = { "TDOM_ASCII" , "FDOM_ASCII" , "TDOM_BINARY" };
  ObserverFormat obsFormat[3]     = { OF_ASCII     , OF_ASCII     , OF_BINARY     };
  ObserverDomain obsDomain[3]     = { OD_TIME      , OD_FREQ      , OD_TIME       };
  ObserverQuantity obsQuantity[3] = { OQ_EH        , OQ_EH        , OQ_EH         };    
  int numScanned = 0;
  char typeStr[TAG_SIZE] = "";
  char waveformName[TAG_SIZE] = "";
  char name[TAG_SIZE] = "";
  int mbbox[6] = { 0 , 0 , 0 , 0 , 0 , 0 };
  int step[3] = { 1 , 1 , 1 };
  ObserverFormat format = OF_ASCII;
  ObserverDomain domain = OD_TIME;
  ObserverQuantity quantity = OQ_EH;
  bool foundType = false;
  WaveformIndex waveformNumber = 0;
  ObserverIndex observerNumber = 0;
  
  /* Get the generic part of the card. */
  numScanned = sscanf( line , "%d %d %d %d %d %d %31s %31s" , 
                       &mbbox[XLO] , &mbbox[XHI] , &mbbox[YLO] , &mbbox[YHI] , &mbbox[ZLO] , &mbbox[ZHI] , 
                       name , typeStr );

  if( numScanned < 8 )
    return false;  

  /* Check observer name is not already defined. */
  if( isObserver( name , &observerNumber ) )
  {
    message( MSG_LOG , 0 , "  Observer %s already defined\n" , name );
    return false;
  }
  
  /* Validate bounding box. */ 
  if( !bboxIsNormal( mbbox ) )
  {
    message( MSG_LOG , 0 , "  Bounding box is abnormal:\n" );
    return false;
  }
  else if( !bboxIsWithin( mbbox , mbox ) )
  {
    message( MSG_LOG , 0 , "  Bounding box is outside mesh:\n" );
    return false;
  }

  /* Find observer type. */
  for( int observer = 0 ; observer < 3 ; observer++ )
    if( strncmp( typeStr , OBSERVER_TYPE_STR[observer] , TAG_SIZE ) == 0 )
    {
      format = obsFormat[observer];     
      domain = obsDomain[observer];
      quantity = obsQuantity[observer];  
      foundType = true;
    }

  if( !foundType )
  {
    message( MSG_LOG , 0 , "  Invalid observer type: %s\n" , typeStr );
    return false;
  }

  /* Get remaining parameters depending on type. */
  if( format == OF_ASCII )
  {
    if( bboxType( mbbox) != BB_POINT )
    {
      message( MSG_LOG , 0 , "  ASCII observers only valid for single node bounding boxes!\n" );
      return false;
    }
    numScanned = sscanf( line , "%d %d %d %d %d %d %31s %31s %31s" , 
                         &mbbox[XLO] , &mbbox[XHI] , &mbbox[YLO] , &mbbox[YHI] , &mbbox[ZLO] , &mbbox[ZHI] , 
                         name , typeStr , waveformName ); 
  }    
  else if( format == OF_BINARY )
  {
    numScanned = sscanf( line , "%d %d %d %d %d %d %31s %31s %d %d %d" , 
                         &mbbox[XLO] , &mbbox[XHI] , &mbbox[YLO] , &mbbox[YHI] , &mbbox[ZLO] , &mbbox[ZHI] , 
                         name , typeStr , &step[XDIR] , &step[YDIR] , &step[ZDIR] );
    if( step[XDIR] < 0 ||  step[YDIR] < 0 ||  step[ZDIR] < 0 )
    {
      message( MSG_LOG , 0 , "  Steps must be positive or zero:\n" );
      return false;
    }
  }
  else
  {
    message( MSG_ERROR , 0 , "  unsupported observer format %s!\n" , OBSERVER_FORMAT_STR[format] );     
  }

  /* If waveform given check it is defined and its get number. */
  if( strncmp( waveformName , "" , TAG_SIZE ) != 0 )
  {
    if( !isWaveform( waveformName , &waveformNumber ) )
    {
      message( MSG_LOG , 0 , "  Waveform %s not defined in field excitation card\n" , waveformName );
      return false;
    }
  }

  addObserver( mbbox , step , name , format , domain , quantity , CACHE_SIZE , false , waveformNumber );

  return true; 

}

/* Add observer to lists. */
ObserverItem *addObserver( int mbbox[6] , int step[3] , char name[TAG_SIZE] , ObserverFormat format , ObserverDomain domain , 
                           ObserverQuantity quantity , unsigned long cacheSize , bool isInternal , WaveformIndex waveformNumber )
{

  ObserverItem *item = NULL;

  if( numObserver == MAX_OBSERVER )
    message( MSG_ERROR , 0 , "*** Error: Maximum number of observers exceeded!\n" );

  item = (ObserverItem *) malloc( sizeof( ObserverItem ) );
  if( !item )
    message( MSG_ERROR , 0 , "*** Error: Failed to allocate observer\n" );

  strncpy( item->name , name , TAG_SIZE );
  item->domain = domain;
  item->format = format;
  item->quantity = quantity;
  item->number = numObserver;
  item->cacheSize = cacheSize;
  item->isInternal = isInternal;
  item->waveformNumber = waveformNumber;
  for( int boundary = XLO ; boundary <= ZHI ; boundary++ ) item->mbbox[boundary] = mbbox[boundary];
  for( int direction = XDIR ; direction <= ZDIR ; direction++ ) item->step[direction] = step[direction];

  /* Add to list. */
  DL_APPEND( observerList , item );
  HASH_ADD_STR( observerHash , name , item );
  numObserver++;
  numObserverDomain[domain]++;
  isObserverDomain[domain] = true;
  numObserverFormat[format]++;
  isObserverFormat[format] = true;
 
  return item;

}

/* Parse far-field requests. */
bool parseFF( char *line )
{

  message( MSG_WARN , 0 , "*** Warning: Parsing of boundary far-fields not implemeneted yet - ignoring FF card.\n" );

  return true;

}

/* Output times. */
bool parseOT( char *line )
{

  if( sscanf( line , "%lu %lu" , &startTimeStep , &stopTimeStep ) != 2 )
    return false;

  if( stopTimeStep < startTimeStep )
  {
    message( MSG_LOG , 0 , "  Stop time-step must be >= start time-step!\n" );
    return false; 
  }

  isOT = true;

  return true;

}

/* Output frequencies. */
bool parseOF( char *line )
{

  double fstart = -1.0;
  double fstop = -1.0;
  unsigned long num = 0;

  if( sscanf( line , "%lf %lf %lu" , &fstart , &fstop , &num ) != 3 )
    return false;

  if( fstart < 0 )
  {
    message( MSG_LOG , 0 , "  Start frequency must be >=0 !\n" );
    return false; 
  }
  else
  {
    startFreq = (real) fstart;
  }

  if( fstop < fstart )
  {
    message( MSG_LOG , 0 , "  Stop frequency must be >= start frequency!\n" );
    return false; 
  }
  else
  {
    stopFreq = (real) fstop;
  }

  if( num < 1UL )
  {
    message( MSG_LOG , 0 , "  Number of frequencies must be >= 1!\n" );
    return false; 
  }
  else
  {
    numFreq = num;
  }

  isOF = true;
    
  return true;

}

/* Initialise observers. */
/* Depends: initGrid, initSimulation, initWaveforms */
void initObservers( void )
{

  ObserverItem *item = NULL;
  ObserverItem **waveformObserverFreqList = NULL;
  real dt = 0.0;
  unsigned long numTimeSteps = 0;
  unsigned long bytes = 0;
  WaveformIndex numWaveformObserver = 0; 
  int bbox[6] = { 0 };
  int step[3] = { 0 };
  
  message( MSG_LOG , 0 , "\nInitialising observers...\n\n" );

  dt = (float) getGridTimeStep();
  numTimeSteps = getNumTimeSteps();
  
  /* Output times. */
  if( !isOT )
  {
    /* Default to all time-steps. */
    startTimeStep = 0UL;
    stopTimeStep = numTimeSteps - 1UL;
  }
  
  startTime = startTimeStep * dt;
  stopTime = stopTimeStep * dt;
  numOutTimeSteps = stopTimeStep - startTimeStep + 1UL; 

  message( MSG_LOG , 0 , "  Observer times: tstart=%ul (%g ns), tstop=%ul (%g ns), numsteps=%lu\n" , 
                            startTimeStep , startTime / 1e-9 , stopTimeStep , stopTime / 1e-9 , numOutTimeSteps );


  /* Output frequencies. */
  if( !isOF )
  {
    /* Set default frequencies if no OF card was given. */
    startFreq = 0.0;
    numFreq = fmax( getNumTimeSteps() / 10 , 1) ;
    stepFreq = 1.0 / ( getNumTimeSteps() * getGridTimeStep() );
    stopFreq = startFreq + ( numFreq - 1 ) * stepFreq;
  }
  else
  {
    stepFreq = ( stopFreq - startFreq ) / ( numFreq - 1 );
  }

  message( MSG_LOG , 0 , "  Observer freqs: fstart=%g MHz, fstop=%g MHz, fstep=%g MHz, fnumber=%lu\n" , 
                            startFreq / 1e6 , stopFreq / 1e6 , stepFreq , numFreq );
 
  /* Add time and frequency domain waveform observers for every waveform. */
  /* Keep mapping from waveform number to its frequency domain observer. */
  numWaveformObserver = getNumberOfWaveforms();
  waveformObserverFreqList =  allocArray( &bytes , sizeof( ObserverItem* ) , 1 , numWaveformObserver );
  for( WaveformIndex wf = 0 ; wf < numWaveformObserver ; wf++ )
  {  
    addObserver( bbox , step , getWaveformName( wf ) , OF_ASCII , OD_TIME , OQ_WF , CACHE_SIZE , true , wf );
    waveformObserverFreqList[wf] = addObserver( bbox , step , getWaveformName( wf ) , OF_ASCII , OD_FREQ , OQ_WF , CACHE_SIZE , true , wf );
  }

  /* Set up observers. */
  DL_FOREACH( observerList , item ) 
  {

    offsetBoundingBox( item->gbbox , item->mbbox , gibox );

    message( MSG_DEBUG3 , 0 , "  Setting %s %s %s observer \"%s\" on [%d,%d,%d,%d,%d,%d]/[%d,%d,%d,%d,%d,%d]\n" , 
             OBSERVER_FORMAT_STR[item->format] , OBSERVER_DOMAIN_STR[item->domain] , 
             OBSERVER_QUANTITY_STR[item->quantity] , item->name ,
             item->mbbox[XLO] , item->mbbox[XHI] , item->mbbox[YLO] , item->mbbox[YHI] , 
	     item->mbbox[ZLO] , item->mbbox[ZHI] , item->gbbox[XLO] , item->gbbox[XHI] , 
	     item->gbbox[YLO] , item->gbbox[YHI] , item->gbbox[ZLO] , item->gbbox[ZHI] );

    /* Number of components. */
    item->numComp = observerCompMap[item->quantity];

    if( item->format == OF_ASCII )   
    {
      if( item->domain == OD_TIME )
      {
        initObserverAsciiTime( item );
      }
      else if( item->domain == OD_FREQ )
      {
        initObserverAsciiFreq( item );
        initObserverDft( item );
      }      
      else
      {
        message( MSG_ERROR , 0 , "*** Error: Unsupported observer format/domain for observer number %lu!\n" , (unsigned long) item->number );
      }
    }
    else if( item->format == OF_BINARY )
    {
      /* No-op. */;    
    }
    else
    {
      message( MSG_ERROR , 0 , "*** Error: Unsupported observer format/domain for observer number %lu!\n" , (unsigned long) item->number );
    }
    
    /* Link frequency domain observers to their reference waveform observers. */
    if( item->domain == OD_FREQ && item->quantity != OQ_WF )
      item->waveformObserver = waveformObserverFreqList[item->waveformNumber];

  } // DL_FOREACH

  /* Deallocate waveform number map. */
  deallocArray( waveformObserverFreqList , 1 , numWaveformObserver );

  /* Set up DFT frequencies. */
  omega = allocArray( &bytes , sizeof( real ) , 1 , numFreq );
  memory.observers += bytes;
  for( int f = 0; f < numFreq ; f++ )
    omega[f] = 2.0 * pi * ( startFreq + f * stepFreq ); 

  if( thereAreObserversFormat( OF_BINARY ) )
    initBinaryObservers( dt );

  return;

}

/* Get number of outputs in bounding box. */
void getNumberOfOutputNodes( int *nx , int *ny , int *nz , int bbox[6] , int step[3] )
{

  *nx = ( 1 + ( bbox[XHI] - bbox[XLO] ) / step[XDIR] );
  *ny = ( 1 + ( bbox[YHI] - bbox[YLO] ) / step[YDIR] );
  *nz = ( 1 + ( bbox[ZHI] - bbox[ZLO] ) / step[ZDIR] );

  return;

}

/* Return true if there are observers. */
bool thereAreObservers( void )
{
  
  bool isTrue = false;
  
  if( numObserver > 0 )
    isTrue = true;
  else
    isTrue = false;
  
  return isTrue;

}

/* Return true if there are observers of given format. */
bool thereAreObserversFormat( ObserverFormat format )
{
  
  bool isTrue = false;
  
  if( format == OF_UNDEFINED )
  {
    if( numObserver > 0 )
      isTrue = true;
    else
      isTrue = false;
  }
  else
  {
    if( numObserverFormat[format] > 0 )
      isTrue = true;
    else
      isTrue = false;
  }
  
  return isTrue;

}

/* Return true if there are observers of given domain. */
bool thereAreObserversDomain( ObserverDomain domain )
{
  
  bool isTrue = false;
  
  if( domain == OD_UNDEFINED )
  {
    if( numObserver > 0 )
      isTrue = true;
    else
      isTrue = false;
  }
  else
  {
    if( numObserverDomain[domain] > 0 )
      isTrue = true;
    else
      isTrue = false;
  }
  
  return isTrue;

}

/* Report observers. */
void reportObservers( void )
{

  ObserverItem *observerItem;

  message( MSG_LOG , 0 , "  Number of observers: %lu\n" , (unsigned long) numObserver );

  DL_FOREACH( observerList , observerItem ) 
  {
    message( MSG_DEBUG3 , 0 , "    Observer \"%s\" (#%lu): Waveform#=%d Format=%s Domain=%s Quantity=%s BBOX=[%d,%d,%d,%d,%d,%d] step=[%d,%d,%d]\n" , 
             observerItem->name , (unsigned long) observerItem->number , observerItem->waveformNumber , 
             OBSERVER_FORMAT_STR[observerItem->format] , OBSERVER_DOMAIN_STR[observerItem->domain] , OBSERVER_QUANTITY_STR[observerItem->quantity] ,
             observerItem->mbbox[XLO] , observerItem->mbbox[XHI] , observerItem->mbbox[YLO] , 
             observerItem->mbbox[YHI] , observerItem->mbbox[ZLO] , observerItem->mbbox[ZHI] ,
             observerItem->step[XDIR] , observerItem->step[YDIR] , observerItem->step[ZDIR] );
  }

  return;

}

/* Update observers. */
void updateObservers( unsigned long tstepNum , real t )
{

  ObserverItem *item;
  bool isOTValid = tstepNum >= startTimeStep && tstepNum <= stopTimeStep;
  
  /* Processing compatible output format for first waveform. */
  if( thereAreObserversFormat( OF_BINARY ) )
    updateExciteDat( getWaveformValue( t , 0 , 0.0 ) );

  DL_FOREACH( observerList , item ) 
  {
    if( item->quantity == OQ_WF && item->domain == OD_TIME )
      updateObserverAsciiTime( item , tstepNum , t );
    else if( item->quantity == OQ_WF && item->domain == OD_FREQ )
      updateObserverAsciiFreq( item , tstepNum , t );     
    else if( item->format == OF_BINARY && isOTValid )
      updateImpulseDat( item );
    else if( item->domain == OD_TIME && item->format == OF_ASCII && isOTValid )
      updateObserverAsciiTime( item , tstepNum , t );
    else if( item->domain == OD_FREQ && item->format == OF_ASCII && isOTValid )
      updateObserverAsciiFreq( item , tstepNum , t );
    else
      /* No-op. Maybe outside OT limits. */;
  }

  return;

}

/* Deallocate observers. */
void deallocObservers( void )
{

  ObserverItem *item , *tmp;

  message( MSG_DEBUG1 , 0 , "Deallocating observers...\n" );


  /* First flush DFT observers. Need to do this before risk of any required */
  /* reference waveform DFT observer being deallocated. */
  DL_FOREACH( observerList , item ) 
  {
    if( item->format == OF_ASCII && item->domain == OD_FREQ )
      flushObserverDft( item );
  }
    
  /* Now deallocate observer hash ansd all observers. */
  HASH_ITER( hh , observerHash , item , tmp )
  {

    if( item->format == OF_ASCII )
    {
      if( item->domain == OD_FREQ )
      {
        deallocObserverAsciiFreq( item );  
        deallocObserverDft( item );
      }
      else
      {
        deallocObserverAsciiTime( item );  
      }
    }
    
    HASH_DELETE( hh , observerHash , item );
    free( item );
  }

  /* Deallocate DFT angular frequency array. */ 
  deallocArray( omega , 1 , numFreq );

  /* Deallocate binary observers. */
  if( thereAreObserversFormat( OF_BINARY ) )
    deallocBinaryObservers();
  
  return;

}

/* Get observer number from name. */
bool isObserver( char *name , ObserverIndex *number )
{

  ObserverItem *item;

  HASH_FIND_STR( observerHash , name , item );
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

/*
 * ASCII observer methods.
 */

/* Initialise ASCII time observer. */
void initObserverAsciiTime( ObserverItem *item )
{
  
  char fileName[PATH_SIZE];
  real physbbox[6];

  switch( item->quantity )
  {
    case OQ_WF:
      sprintf( fileName , "wf_%s_td.asc", item->name );
      item->outputFile = fopen( fileName , "w" );
      if( !item->outputFile )
        message( MSG_ERROR , 0 , "*** Error: Failed to open time domain output file for waveform number %lu\n" , (unsigned long)item->waveformNumber );  
      fprintf( item->outputFile , "# Waveform# %lu\n" , (unsigned long)item->waveformNumber );
      fprintf( item->outputFile , "# %6s %16s %16s\n" , "ts (-)" , "t (s)" , "wf (-)" );
      break;
    case OQ_EH:
      bboxInPhysicalUnits( physbbox , item->mbbox );
      sprintf( fileName , "eh_%s_td.asc", item->name );
      item->outputFile = fopen( fileName , "w" );
      if( !item->outputFile )
        message( MSG_ERROR , 0 , "*** Error: Failed to open output file for observer number %lu\n" , (unsigned long) item->number );
      fprintf( item->outputFile , "# (%d,%d,%d)->(%g,%g,%g)\n" , item->mbbox[XLO] , item->mbbox[YLO] , item->mbbox[ZLO] , physbbox[XLO] , physbbox[YLO] , physbbox[ZLO] );
      fprintf( item->outputFile , "# %6s %16s %16s %16s %16s %16s %16s %16s\n" , "ts (-)" , "t (s)" , "Ex (V/m)" , "Ey (V/m)" ,"Ez (V/m)" ,"Hx (A/m)" ,"Hy (A/m)" ,"Hz (A/m)" );
      break;
    default:
      assert( 0 );
      break;
  }

  return;
}

/* Update ASCII time observer. */
void updateObserverAsciiTime( ObserverItem *item , unsigned long tstepNum , real t )
{
  
  int i;
  int j;
  int k;
  real value[MAX_COMP];

  /* Get observable value. */
  switch( item->quantity )
  {
    case OQ_WF:       
      value[0] = getWaveformValue( t , item->waveformNumber , 0.0 );
      break;
    case OQ_EH:
      i = item->gbbox[XLO];
      j = item->gbbox[YLO];
      k = item->gbbox[ZLO];
      value[EX] = UNSCALE_Ex( Ex[i][j][k] , i );
      value[EY] = UNSCALE_Ey( Ey[i][j][k] , j );
      value[EZ] = UNSCALE_Ez( Ez[i][j][k] , k );
      value[HX] = UNSCALE_Hx( Hx[i][j][k] , i );
      value[HY] = UNSCALE_Hy( Hy[i][j][k] , j );
      value[HZ] = UNSCALE_Hz( Hz[i][j][k] , k );
      break;
    default:
      assert( 0 );
      break;
  }
  
  /* Write out value. */
  fprintf( item->outputFile , "%8lu %16.8e ", tstepNum , t );
  for( int comp = 0; comp < item->numComp ; comp++ )
    fprintf( item->outputFile , "%16.8e ", value[comp] );
  fprintf( item->outputFile , "\n" );

  return;

}

/* Deallocate ASCII time observer. */
void deallocObserverAsciiTime( ObserverItem *item )
{
  
  switch( item->quantity )
  {
    case OQ_WF:
    case OQ_EH:
      fclose( item->outputFile );
      break;
    default:
      assert( 0 );
      break;
  }

  return;
}

/* Initialise ASCII frequency observer. */
void initObserverAsciiFreq( ObserverItem *item )
{
  
  char fileName[PATH_SIZE];
  real physbbox[6];
  
  switch( item->quantity )
  {
    case OQ_WF:
      sprintf( fileName , "wf_%s_fd.asc", item->name );
      item->outputFile = fopen( fileName , "w" );
      if( !item->outputFile )
        message( MSG_ERROR , 0 , "*** Error: Failed to open frequency domain output file for waveform number %lu\n" , (unsigned long)item->waveformNumber );  
      fprintf( item->outputFile , "# Waveform# %lu\n" , (unsigned long)item->waveformNumber );
      fprintf( item->outputFile , "# %14s %16s %16s\n" , "f (Hz)" , "Re(wf) (-)" , "Im(wf) (-)" );
      break;
    case OQ_EH:
      bboxInPhysicalUnits( physbbox , item->mbbox );
      sprintf( fileName , "eh_%s_fd.asc", item->name );
      item->outputFile = fopen( fileName , "w" );
      if( !item->outputFile )
        message( MSG_ERROR , 0 , "*** Error: Failed to open output file for observer number %lu\n" , (unsigned long) item->number );
      fprintf( item->outputFile , "# (%d,%d,%d)->(%g,%g,%g)\n" , item->mbbox[XLO] , item->mbbox[YLO] , item->mbbox[ZLO] , physbbox[XLO] , physbbox[YLO] , physbbox[ZLO] );
      fprintf( item->outputFile , "# %14s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s\n" , 
                              "f (Hz)" , "Re(Ex) (V/m)" , "Im(Ex) (V/m)" , "Re(Ey) (V/m)" , "Im(Ey) (V/m)", 
                              "Re(Ez) (V/m)" , "Im(Ez) (V/m)", "Re(Hx) (V/m)" , "Im(Hx) (V/m)", "Re(Hy) (V/m)" , 
                              "Im(Hy) (V/m)", "Re(Hz) (V/m)" , "Im(Hz) (V/m)");
      break;
    default:
      assert( 0 );
      break;
  }
  
  return;
  
}

/* Update ASCII frequency observer. */
void updateObserverAsciiFreq( ObserverItem *item , unsigned long tstepNum , real t )
{
  int i;
  int j;
  int k;
  real value[MAX_COMP];
  real omegat;
  
  /* Get observable value. */
  switch( item->quantity )
  {
    case OQ_WF:       
      value[0] = getWaveformValue( t , item->waveformNumber , 0.0 );
      break;
    case OQ_EH:
      i = item->gbbox[XLO];
      j = item->gbbox[YLO];
      k = item->gbbox[ZLO];
      value[EX] = UNSCALE_Ex( Ex[i][j][k] , i );
      value[EY] = UNSCALE_Ey( Ey[i][j][k] , j );
      value[EZ] = UNSCALE_Ez( Ez[i][j][k] , k );
      value[HX] = UNSCALE_Hx( Hx[i][j][k] , i );
      value[HY] = UNSCALE_Hy( Hy[i][j][k] , j );
      value[HZ] = UNSCALE_Hz( Hz[i][j][k] , k );
      break;
    default:
      assert( 0 );
      break;
  }

  /* Add value to DFT. */
  for( int comp = 0; comp < item->numComp ; comp++ )
  {
    for( int f = 0; f < numFreq ; f++ )
    {
      omegat = omega[f] * t;
      item->dft_real[comp][f] += value[comp] * cos( omegat );
      item->dft_imag[comp][f] -= value[comp] * sin( omegat );
    }
  }
  
  return;
}

/* Deallocate ASCII frequency observer. */
void deallocObserverAsciiFreq( ObserverItem *item )
{
  
  switch( item->quantity )
  {
    case OQ_WF:
    case OQ_EH:
      fclose( item->outputFile );
      break;
    default:
      assert( 0 );
      break;
  }

  return;
}

/* Initialise observer DFT arrays. */
void initObserverDft( ObserverItem *item )
{

  unsigned long bytes;

  switch( item->quantity )
  {
    case OQ_WF:
    case OQ_EH:
      item->dft_real = allocArray( &bytes , sizeof( real ) , 2 , item->numComp , numFreq );
      item->dft_imag = allocArray( &bytes , sizeof( real ) , 2 , item->numComp , numFreq );
      memory.observers += bytes;
      for( int comp = 0; comp < item->numComp ; comp++ )
        for( int f = 0; f < numFreq ; f++ )
        {
          item->dft_real[comp][f] = 0.0;
          item->dft_imag[comp][f] = 0.0;
        }
      break;
    default:
      assert( 0 );
      break;
  }
  
  return;

}

/* Write out observer DFT normalised by waveform DFT. */
void flushObserverDft( ObserverItem *item )
{

  real denom;
  real realPart;
  real imagPart;
  real wf_r;
  real wf_i;
  real comp_r;
  real comp_i;

  switch( item->quantity )
  {
    case OQ_WF:
      for( int f = 0; f < numFreq ; f++ )
      {
        fprintf( item->outputFile , "%16.8e " , startFreq + f * stepFreq );
        fprintf( item->outputFile , "%16.8e %16.8e " , item->dft_real[0][f] , item->dft_imag[0][f] );
        fprintf( item->outputFile , "\n" );
      }
      break;
    case OQ_EH:
      for( int f = 0; f < numFreq ; f++ )
      {
        wf_r = item->waveformObserver->dft_real[0][f];
        wf_i = item->waveformObserver->dft_imag[0][f];
        denom = wf_r * wf_r + wf_i * wf_i;
        fprintf( item->outputFile , "%16.8e " , startFreq + f * stepFreq );
        for( int comp = 0; comp < item->numComp ; comp++ )
        {
          comp_r = item->dft_real[comp][f];
          comp_i = item->dft_imag[comp][f];
          realPart = ( comp_r * wf_r + comp_i * wf_i ) / denom;
          imagPart = ( comp_i * wf_r - comp_r * wf_i ) / denom;
          fprintf( item->outputFile , "%16.8e %16.8e " , realPart , imagPart );
        }
        fprintf( item->outputFile , "\n" );
      }
      break;
    default:
      assert( 0 );
      break;
  }
  
  return;

}

/* Deallocate observer DFT arrays. */
void deallocObserverDft( ObserverItem *item )
{

  switch( item->quantity )
  {
    case OQ_WF:
    case OQ_EH:
      deallocArray( item->dft_real , 2 , item->numComp , numFreq );
      deallocArray( item->dft_imag , 2 , item->numComp , numFreq );
      break;
    default:
      assert( 0 );
      break;
  }

  return;

}

/*
 * Binary format methods.
 */

/* Maximum size of comment line in process.dat file. */
#define COMMENT_BUFFER_SIZE 70

/* impulse.dat file pointer. */
static FILE *impulseDatFile = NULL;

/* excite.dat file pointer. */
static FILE *exciteFile = NULL;

/* Initialise binary observers. */
void initBinaryObservers( real dt )
{

  /* Initialise impulse.dat file. */
  initImpulseDat( dt );
    
  /* Write post processing file. */
  initExciteDat();

  /* Write processing data file. */
  writeProcessDat();
    
  return;

}

/* Write processing tools comptabile excitation data file. */
void initExciteDat( void )
{

  exciteFile = fopen ( "excite.dat" , "w" );
  if( NULL == exciteFile ) 
    message( MSG_ERROR , 0 , "*** Error: Failed to open %s file\n" , "excite.dat" );

  fprintf ( exciteFile, " %d\n", 8 );

  return;

}

/* Write processing data file. */
void writeProcessDat( void )
{

  char commentBuffer[COMMENT_BUFFER_SIZE] = "";
  FILE *processFile;
  real fstart;
  real fstop;
  real fstep;
  real meshSize;
  ObserverItem *item;
  ObserverItem *firstTimeDomBinaryItem = NULL;

  processFile = fopen ( "process.dat" , "w" );
  if( NULL == processFile ) 
    message( MSG_ERROR , 0 , "*** Error: Failed to open %s file\n" , "process.dat" );

  meshSize = getGridTimeStep() * ( 2.0 * c0 ); 
  fstart = 0.0;
  fstep = 1.0 / ( getNumTimeSteps() * getGridTimeStep() );
  fstop = getNumTimeSteps() * fstep;

  /* Processing tools bjork if comment too long! */
  snprintf ( commentBuffer , COMMENT_BUFFER_SIZE - 1 , "%s" , getCommentReference() );
  fprintf ( processFile , "CE %s\n" , commentBuffer );

  fprintf ( processFile , "%lu\n" , (unsigned long) numObserverFormat[OF_BINARY] );

  DL_FOREACH( observerList , item ) 
  {
    if( item->domain == OD_TIME && item->format == OF_BINARY )
    {
      if( firstTimeDomBinaryItem == NULL )
        firstTimeDomBinaryItem = item;

      fprintf( processFile , " %d", item->gbbox[XLO] - gibox[XLO] );
      fprintf( processFile , " %d", item->gbbox[XHI] - gibox[XLO] );
      fprintf( processFile , " %d", item->step[XDIR] );
      fprintf( processFile , " %d", item->gbbox[YLO] - gibox[YLO] );
      fprintf( processFile , " %d", item->gbbox[YHI] - gibox[YLO] );
      fprintf( processFile , " %d", item->step[YDIR] );
      fprintf( processFile , " %d", item->gbbox[ZLO] - gibox[ZLO] );
      fprintf( processFile , " %d", item->gbbox[ZHI] - gibox[ZLO] );
      fprintf( processFile , " %d", item->step[ZDIR] );
      fprintf (processFile, "\n");

    }
  }

  /* Note add one to time step number for processing compatibility! */
  fprintf ( processFile , "%lu %lu\n" , startTimeStep + 1 , stopTimeStep + 1  );
  fprintf ( processFile , "%g %g %g\n" , fstart , fstop , fstep );
  fprintf ( processFile , "%g\n" , meshSize );

  if( firstTimeDomBinaryItem )
  {
    fprintf( processFile , " %d", firstTimeDomBinaryItem->gbbox[XLO] - gibox[XLO] );
    fprintf( processFile , " %d", firstTimeDomBinaryItem->gbbox[XHI] - gibox[XLO] );
    fprintf( processFile , " %d", firstTimeDomBinaryItem->gbbox[YLO] - gibox[YLO] );
    fprintf( processFile , " %d", firstTimeDomBinaryItem->gbbox[YHI] - gibox[YLO] );
    fprintf( processFile , " %d", firstTimeDomBinaryItem->gbbox[ZLO] - gibox[ZLO] );
    fprintf( processFile , " %d", firstTimeDomBinaryItem->gbbox[ZHI] - gibox[ZLO] );
    fprintf ( processFile , " %d\n" , EX + 1 );
  }
  else
  {
    message( MSG_WARN , 0 , "*** Warning: process.dat contains no valid observers!\n" );  
  }

  fclose( processFile );

  return;

}

/* Open and initialise impulse.dat file. */
void initImpulseDat( real dt )
{

  /* Open binary impulse dat file. */
  impulseDatFile = fopen( "impulse.dat" , "wb" );
  if( !impulseDatFile )
    message( MSG_ERROR , 0 , "*** Error: Failed to open binary data file %s\n" , "impulse.dat" );

  /* Write header. */
  fwrite( &numOutTimeSteps , sizeof( int ) , (size_t) 1 , impulseDatFile );
  fwrite( &dt , sizeof( float ) , (size_t) 1 , impulseDatFile );    

  return;

}   

/* Update excite.dat file. */
void updateExciteDat( real value )
{

  fprintf( exciteFile , "%16.8e\n", value ); 

  return;

}

/* Update impulse.dat file. */
void updateImpulseDat( ObserverItem *item )
{

  /* Types for impulse.dat - must match processing tools on given platform. */
  int i, j, k;
  float outEx , outEy , outEz , outHx , outHy , outHz;
  int iout , jout , kout;

  for( k = item->gbbox[ZLO]; k <= item->gbbox[ZHI] ; k += item->step[ZDIR] )
    for( j = item->gbbox[YLO] ; j <= item->gbbox[YHI] ; j += item->step[YDIR] )
      for( i = item->gbbox[XLO] ; i <= item->gbbox[XHI] ; i += item->step[XDIR] )
        {
          outEx = UNSCALE_Ex( Ex[i][j][k] , i );
          outEy = UNSCALE_Ey( Ey[i][j][k] , j );
          outEz = UNSCALE_Ez( Ez[i][j][k] , k );
          outHx = UNSCALE_Hx( Hx[i][j][k] , i );
          outHy = UNSCALE_Hy( Hy[i][j][k] , j );
          outHz = UNSCALE_Hz( Hz[i][j][k] , k );
          iout = i - gibox[XLO];
          jout = j - gibox[YLO];
          kout = k - gibox[ZLO];
          fwrite( &iout , sizeof( iout ) , (size_t) 1 , impulseDatFile );
          fwrite( &jout , sizeof( jout ) , (size_t) 1 , impulseDatFile );
          fwrite( &kout , sizeof( kout ) , (size_t) 1 , impulseDatFile );
          fwrite( &outEx , sizeof( outEx ) , (size_t) 1 , impulseDatFile );
          fwrite( &outEy , sizeof( outEy ) , (size_t) 1 , impulseDatFile );
          fwrite( &outEz , sizeof( outEz ) , (size_t) 1 , impulseDatFile );
          fwrite( &outHx , sizeof( outHx ) , (size_t) 1 , impulseDatFile );
          fwrite( &outHy , sizeof( outHy ) , (size_t) 1 , impulseDatFile );
          fwrite( &outHz , sizeof( outHz ) , (size_t) 1 , impulseDatFile );
        }
    
  return;

}

/* Dealloc binary observers. */
void deallocBinaryObservers( void )
{

  /* Close processing tools excitation file. */
  fclose( exciteFile );
 
  /* Close impulse data file. */ 
  fclose( impulseDatFile );
  
  return;

}

/*
 * Graphics methods.
 */

/* Output gnuplot compatible data file for observers. */
void gnuplotObservers( void )
{

  char observerFileName[] = "gnuplot-observer.dat";
  FILE *outputFile; 
  ObserverItem *item;

  outputFile = fopen( observerFileName , "w" );
  if( !outputFile )
    message( MSG_ERROR , 0 , "*** Error: Failed to open observer output file %s\n" , observerFileName );

  gnuplotProblemSize( outputFile , mbox );

  DL_FOREACH( observerList , item ) 
  {
    gnuplotBoundingBoxNodes( outputFile , item->mbbox , item->step );
  }
			
  fclose( outputFile );

  return;

}

/* Output gmsh compatible data file for observers. */
void gmshObservers( void )
{

  ObserverItem *item;
  unsigned long entityNumber;
  char name[GMSH_NAME_LENGTH];
  
  DL_FOREACH( observerList , item ) 
  {
    /* gmsh doesn't seem to like multiple identical nodes elements with different entity 
     * and different physical numbers. OPs likely to have overlapping nodes to give them
     * entity number zero (no tag).
     */
    entityNumber = 0UL;
    snprintf( name , GMSH_NAME_LENGTH - 1 , "OP_%s" , item->name );
    gmshAddEntity( entityNumber , BB_POINT , name , item->mbbox , item->step );
  }
                        
  return;

}
