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
#include <assert.h>
#include <string.h>

#include "line.h"
#include "utlist.h"
#include "alloc_array.h"
#include "grid.h"
#include "message.h"
#include "wire.h"
#include "medium.h"
#include "bounding_box.h"
#include "gnuplot.h"
#include "gmsh.h"
#include "memory.h"

/* 
 * Private data.
 */

/* Number of lines. */
static LineIndex numLine = 0;

/* Existance flag for lines of each wire type, including undefined. */
static bool isLineType[NUM_WIRE_TYPES+1] = { false };

/* List of lines. */
static LineItem *lineList = NULL;

/* 
 * Private method interfaces. 
 */

void addLine( int mbbox[6] , char *wireName , WireIndex wireNumber , WireEndType lowEndType , WireEndType highEndType );

/*
 * Method Implementations.
 */

/* Add a line to the lists. */
void addLine( int mbbox[6] , char *wireName , WireIndex wireNumber , WireEndType lowEndType , WireEndType highEndType )
{

  LineItem *item = NULL;

  if( numLine == MAX_LINE )
    message( MSG_ERROR , 0 , "*** Error: Maximum number of lines exceeded!\n" );

  item = (LineItem *) malloc( sizeof( LineItem ) );
  if( !item )
    message( MSG_ERROR , 0 , "*** Error: Failed to allocate line.\n" );

  for( int i = XLO ; i <= ZHI ; i++ ) item->mbbox[i] = mbbox[i];
  item->wireNumber = wireNumber;
  strncpy( item->wireName , wireName , TAG_SIZE );
  item->lowEndType = lowEndType;
  item->highEndType = highEndType;

  /* Add to list. */
  DL_APPEND( lineList , item );
  numLine++;

  return;

}

/* Parse internal surfaces. */
bool parseTW( char *line )
{

  int numScanned = 0;
  char wireName[TAG_SIZE] = "";
  int mbbox[6];
  WireIndex wireNumber = 0;
  char lowEndTag[TAG_SIZE] = "";
  char highEndTag[TAG_SIZE] = "";
  WireEndType lowEndType = WE_UNDEFINED;
  WireEndType highEndType = WE_UNDEFINED;

  numScanned = sscanf( line , "%d %d %d %d %d %d %31s %31s %31s" , 
                       &mbbox[XLO] , &mbbox[XHI] , &mbbox[YLO] , &mbbox[YHI] , &mbbox[ZLO] , &mbbox[ZHI] , 
                       wireName , lowEndTag , highEndTag );

  if( numScanned < 7 )
    return false;  

  /* Validate bounding box. */ 
  if( !bboxIsNormal( mbbox ) )
  {
    message( MSG_LOG , 0 , "  Bounding box is abnormal:\n" );
    return false;
  }
  else if( !bboxIsWithin( mbbox , mbox ) )
  {
    message( MSG_LOG , 0 , "Bounding box is outside mesh:\n" );
    return false;
  }
  else if( bboxType( mbbox ) != BB_LINE )
  {
    message( MSG_LOG , 0 , "  Bounding box is not a line!\n" );
    return false;
  }

  /* Check boundary exists. */
  if( !isWire( wireName , &wireNumber ) )
  {
    message( MSG_LOG , 0 , "  Wire %s not defined in TW card\n" , wireName );
    return false;
  }

  /* Process end types. */
  if( numScanned >= 8 )
  {
    for( WireEndType endType = 0 ; endType < NUM_WIRE_END_TYPES ; endType++ )
      if( strncmp( WIRE_END_TYPE_STR[endType] , lowEndTag , TAG_SIZE ) == 0 )
        lowEndType = endType;
    if( lowEndType == WE_UNDEFINED )
      message( MSG_ERROR , 0 , "  Invalid wire low end type %s in TW card\n" , lowEndTag );
  }
 
  if( numScanned >= 9 )
  {
    for( WireEndType endType = 0 ; endType < NUM_WIRE_END_TYPES ; endType++ )
      if( strncmp( WIRE_END_TYPE_STR[endType] , highEndTag , TAG_SIZE ) == 0 )
        highEndType = endType;
    if( highEndType == WE_UNDEFINED )
      message( MSG_ERROR , 0 , "  Invalid wire high end type %s in TW card\n" , highEndTag );
  }

  addLine( mbbox , wireName , wireNumber , lowEndType , highEndType );

  /* This is needed by mesh renderers in gvulture.*/
  isLineType[TW_UNDEFINED] = true;

  return true;

}

/* Initialise lines. */
/* Depends:  */
void initLines( void )
{

  LineItem *item;
  WireType type;

  message( MSG_LOG , 0 , "\nInitialising lines...\n\n" );

  DL_FOREACH( lineList , item ) 
  {

    type = getWireType( item->wireNumber );
    isLineType[type] = true;

    switch( type )
    {
    case TW_PEC:
      offsetBoundingBox( item->gbbox , item->mbbox , gibox );
      message( MSG_DEBUG3 , 0 , "  Setting PEC line medium#%lu on [%d,%d,%d,%d,%d,%d]/[%d,%d,%d,%d,%d,%d]\n" , 
               MT_PEC , item->mbbox[XLO] , item->mbbox[XHI] , item->mbbox[YLO] , item->mbbox[YHI] , 
	       item->mbbox[ZLO] , item->mbbox[ZHI] , item->gbbox[XLO] , item->gbbox[XHI] , 
	       item->gbbox[YLO] , item->gbbox[YHI] , item->gbbox[ZLO] , item->gbbox[ZHI] );
      setMediumOnGrid( item->gbbox , MT_PEC , FACE_MASK_ALL );
      break;
    case TW_FREE_SPACE:
      offsetBoundingBox( item->gbbox , item->mbbox , gibox );
      message( MSG_DEBUG3 , 0 , "  Setting FREE_SPACE line medium#%lu on [%d,%d,%d,%d,%d,%d]/[%d,%d,%d,%d,%d,%d]\n" , 
               MT_PEC , item->mbbox[XLO] , item->mbbox[XHI] , item->mbbox[YLO] , item->mbbox[YHI] , 
               item->mbbox[ZLO] , item->mbbox[ZHI] , item->gbbox[XLO] , item->gbbox[XHI] , 
               item->gbbox[YLO] , item->gbbox[YHI] , item->gbbox[ZLO] , item->gbbox[ZHI] );
      setMediumOnGrid( item->gbbox , MT_FREE_SPACE , FACE_MASK_ALL );
      break;
    default:
      break;
    }
  }
  
  return;

}

/* Report lines. */
void reportLines( void )
{

  LineItem *lineItem;
  LineIndex counter = 0;

  message( MSG_LOG , 0 , "  Number of lines: %lu\n" , numLine );
  DL_FOREACH( lineList , lineItem ) 
  {
    message( MSG_DEBUG3 , 0 , "    Line #%lu: Wire=%s Wire#=%lu BBOX=[%d,%d,%d,%d,%d,%d]\n" , 
             (unsigned long) counter , lineItem->wireName , (unsigned long) lineItem->wireNumber ,
             lineItem->mbbox[XLO] , lineItem->mbbox[XHI] , lineItem->mbbox[YLO] , 
             lineItem->mbbox[YHI] , lineItem->mbbox[ZLO] , lineItem->mbbox[ZHI] );
    counter++;
  }

  return;

}

/* Return true if there are internal surfaces of given boundary type. */
bool thereAreLines( WireType type )
{
  
  return isLineType[type];

}

/* Update lines E fields. */
void updateLinesEfield( void )
{

  return;

}

/* Update lines H fields. */
void updateLinesHfield( void )
{

  return;

}

/* Deallocate lines. */
void deallocLines( void )
{

  LineItem *item , *tmp;

  message( MSG_DEBUG1 , 0 , "Deallocating lines...\n" );

  /* Free line list. */
  DL_FOREACH_SAFE( lineList , item , tmp ) 
  {
    DL_DELETE( lineList , item );
    free( item );
  }

  return;

}

/* Draw lines. */
void gnuplotLines( void )
{

  char lineFileName[] = "gnuplot-wires.dat";
  FILE *outputFile;
  LineItem *item;

  outputFile = fopen( lineFileName , "w" );
  if( !outputFile )
    message( MSG_ERROR , 0 , "*** Error: Failed to open line output file %s\n" , lineFileName );

  gnuplotProblemSize( outputFile , mbox );

  DL_FOREACH( lineList , item ) 
  {
    gnuplotBoundingBox( outputFile , item->mbbox );
  }

  fclose( outputFile );

  return;

}

/* Draw lines. */
void gmshLines( void )
{

  int step[3] = { 1 , 1 , 1 };
  LineItem *item;
  unsigned long entityNumber;
  char name[GMSH_NAME_LENGTH];
    
  DL_FOREACH( lineList , item ) 
  {
    entityNumber = gmshGetEntityNumber();
    snprintf( name , GMSH_NAME_LENGTH - 1 , "WT_%s" , getWireName( item->wireNumber ) );
    gmshAddEntity( entityNumber , BB_LINE , name , item->mbbox , step );
  }
  
  return;

}

