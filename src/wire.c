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
#include <assert.h>
#include <stdio.h>

#include "wire.h"
#include "alloc_array.h"
#include "message.h"
#include "grid.h"
#include "pml.h"
#include "physical.h"

/*
 * Global data.
 */

/* Wire type strings. */
char WIRE_TYPE_STR[NUM_WIRE_TYPES][11] = { "PEC" , "FREE_SPACE" };

/* Wire end type strings. */
char WIRE_END_TYPE_STR[NUM_WIRE_END_TYPES][7] = { "THRU" , "END" , "CORNER" };

/* Mapping of wire type to medium type. */
MediumType WIRE_MEDIUM_TYPE[NUM_WIRE_TYPES] = { MT_PEC , MT_FREE_SPACE };

/* 
 * Private data. 
 */

/* Number of wire types. */
static WireIndex numWire = 0;

/* Existance flag for wires of each type, including undefined. */
static bool isWireType[NUM_WIRE_TYPES+1] = { false };

/* List of wire types. */
static WireItem *wireList = NULL;

/* Hash of wire types by name. */
static WireItem *wireHash = NULL;

/* Hash of wire types by number. */
static WireItem *wireNumberHash = NULL;

/* 
 * Private method interfaces. 
 */


/*
 * Method Implementations.
 */

/* Parse boundaries. */
bool parseWT( char *line )
{

  int numScanned = 0;
  char name[TAG_SIZE] = "";
  char typeStr[TAG_SIZE] = "";
  WireType type = TW_PEC;
  bool foundType = false;
  int radius = 0.0;
  WireIndex number;

  numScanned = sscanf( line , "%31s %31s" , name , typeStr );
  if( numScanned < 2 )
    return false;  
 
  /* Check not already defined. */
  if( isWire( name , &number ) )
  {
    message( MSG_LOG , 0 , "  Wire %s already defined\n" , name );
    return false;
  }

  /* Find type. */
  for( int wire = 0 ; wire < NUM_WIRE_TYPES ; wire++ )
    if( strncmp( typeStr , WIRE_TYPE_STR[wire] , TAG_SIZE ) == 0 )
    {
      type = (WireType)wire;      
      foundType = true;
    }

  if( !foundType )
  {
    message( MSG_LOG , 0 , "  Invalid wire type: %s\n" , typeStr );
    return false;
  }

  /* Validate parameters. */
  switch( type )
  {
  case TW_PEC:
    break;
  default:
    assert( 0 );
    break;
  }

  addWire( name , type , radius );

  return true;

}

/* Add wire to lists. */
void addWire( char *name , WireType type , real radius )
{

  WireItem *item = NULL;

  if( numWire == MAX_WIRE )
    message( MSG_ERROR , 0 , "*** Error: Maximum number of wires exceeded!\n" );

  item = (WireItem *) malloc( sizeof( WireItem ) );
  if( !item )
    message( MSG_ERROR , 0 , "*** Error: Failed to allocate wire!\n" );

  strncpy( item->name , name , TAG_SIZE );
  item->type = type;
  item->number = numWire;
  item->radius = radius;

  /* Add to list and hash. */
  DL_APPEND( wireList , item );
  HASH_ADD_STR( wireHash , name , item );
  HASH_ADD( hhint , wireNumberHash , number , sizeof( numWire ) , item );
  numWire++;
  isWireType[type] = true;
  isWireType[TW_UNDEFINED] = true;
  
  return;

}

/* Get wire number from name. */
bool isWire( char *name , WireIndex *number )
{

  WireItem *item;

  HASH_FIND_STR( wireHash , name , item );
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

/* Initialise wires. */
/* Depends: initGrid */
void initWires( void )
{

  WireItem *item;
  
  message( MSG_LOG , 0 , "  Initialising wires...\n" );
  
  /* Determine update coefficient for each material. */
  DL_FOREACH( wireList , item ) 
  {
    switch( item->type )
    {
    default:
      break;
    }
  }

  return;

}

/* Deallocate wires. */
void deallocWires( void )
{

  WireItem *item , *tmp;

  message( MSG_DEBUG1 , 0 , "Deallocating wires...\n" );

  /* Free the number hash. */
  HASH_ITER( hhint , wireNumberHash , item , tmp )
  {
    HASH_DELETE( hhint , wireNumberHash , item );
  }

  /* Free the name hash and wires. */
  HASH_ITER( hh , wireHash , item , tmp )
  {
    HASH_DELETE( hh , wireHash , item );
    free( item );
  }
 
  return;

}

/* Report wires. */
void reportWires( void )
{

  WireItem *item;

  message( MSG_LOG , 0 , "  Number of wires: %lu\n" , numWire );

  DL_FOREACH( wireList , item ) 
  {
    message( MSG_DEBUG3 , 0 , "    Wire #%lu: Name=%s Type=%s Radius=%e\n" , 
             (unsigned long) item->number , item->name , WIRE_TYPE_STR[item->type] , item->radius );
  }

  return;

}

/* Return true if there are wires of given type. */
bool thereAreWires( WireType type )
{
  
  return isWireType[type];

}

/* Get wire pointer. */
WireItem *getWire( WireIndex number )
{

  WireItem *item;

  HASH_FIND( hhint , wireNumberHash , &number , sizeof( number ) , item );
  if( !item)
    assert( 0 );

  return item;

}

/* Get wire types by number. */
WireType getWireType( WireIndex number )
{

  WireItem *item;

  HASH_FIND( hhint , wireNumberHash , &number , sizeof( number ) , item );
  if( !item)
    assert( 0 );

  return item->type;

}

/* Get pointer to wire name by number. */
char *getWireName( WireIndex number )
{

  WireItem *item;

  HASH_FIND( hhint , wireNumberHash , &number , sizeof( number ) , item );
  if( !item)
    assert( 0 );

  return item->name;

}

