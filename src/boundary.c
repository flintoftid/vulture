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

#include "boundary.h"
#include "alloc_array.h"
#include "message.h"
#include "grid.h"
#include "pml.h"
#include "physical.h"

#ifdef WITH_SIBC
#include "sibc.h"
#endif

/*
 * Global data.
 */

/* Boundary type strings. */
char BOUNDARY_TYPE_STR[NUM_BOUNDARY_TYPES][11] = { "PMC" , "PML" , "PEC" , "FREE_SPACE" , "PERIODIC" , "MUR" , "SIBC" };

/* Mapping of boundary type to medium type. */
MediumType BOUNDARY_MEDIUM_TYPE[NUM_BOUNDARY_TYPES] = { MT_UNDEFINED , MT_PEC , MT_PEC , MT_FREE_SPACE , MT_UNDEFINED , MT_UNDEFINED , MT_PEC };

/* 
 * Private data. 
 */

/* Number of boundary types. */
static BoundaryIndex numBoundary = 0;

/* Existance flag for boundaries of each type, including undefined. */
static bool isBoundaryType[NUM_BOUNDARY_TYPES+1] = { false };

/* List of boundary types. */
static BoundaryItem *boundaryList = NULL;

/* Hash of boundary types by name. */
static BoundaryItem *boundaryHash = NULL;

/* Hash of boun dary types by number. */
static BoundaryItem *boundaryNumberHash = NULL;

/* 
 * Private method interfaces. 
 */


/*
 * Method Implementations.
 */

/* Parse extra boundary parameters - PML depths. */
bool parseBE( char *line )
{

  BoundaryIndex number;
  int depth[6];

  if( sscanf( line , "%d %d %d %d %d %d" , &depth[XLO] , &depth[XHI] , &depth[YLO] , &depth[YHI] , &depth[ZLO] , &depth[ZHI] ) != 6 )
    return false;  

  for( int surface = XLO ; surface <= ZHI ; surface++ )
  {
    if( depth[surface] < 0 )
    {
      message( MSG_LOG , 0 , "  Invalid depth, %d, for PML on %s boundary\n" , depth[surface] , BOUNDARY_TYPE_STR[surface] );
      return false;
    }
    else
    {
      if( !isBoundary( FACE[surface] , &number ) )
        message( MSG_LOG , 0 , "  Failed to set depth for PML on %s boundary\n" , BOUNDARY_TYPE_STR[surface] );
      setBoundaryNumLayers( number , depth[surface] );
    }

  }

  return true;

}

/* Parse boundaries. */
bool parseBT( char *line )
{

  int numScanned = 0;
  char name[TAG_SIZE] = "";
  char typeStr[TAG_SIZE] = "";
  BoundaryType type = BT_PEC;
  bool foundType = false;
  int numLayers = 0;
  int order = 0;
  real n_eff = 0.0;
  real refCoeff = 0.0;
  real kmax = 0.0;
  char fileName[PATH_SIZE] = "";
  BoundaryIndex number;
  real S_TM[2][2] = { { -1.0 , 0.0 } , { 0.0 , -1.0 } };
  real S_TE[2][2] = { { -1.0 , 0.0 } , { 0.0 , -1.0 } };

  numScanned = sscanf( line , "%31s %31s" , name , typeStr );
  if( numScanned < 2 )
    return false;  
 
  /* Check not already defined. */
  if( isBoundary( name , &number ) )
  {
    message( MSG_LOG , 0 , "  Boundary %s already defined\n" , name );
    return false;
  }

  /* Find type. */
  for( int boundary = 0 ; boundary < NUM_BOUNDARY_TYPES ; boundary++ )
    if( strncmp( typeStr , BOUNDARY_TYPE_STR[boundary] , TAG_SIZE ) == 0 )
    {
      type = (BoundaryType)boundary;      
      foundType = true;
    }

  if( !foundType )
  {
    message( MSG_LOG , 0 , "  Invalid boundary type: %s\n" , typeStr );
    return false;
  }

  /* Validate parameters. */
  switch( type )
  {
  case BT_PEC:
    numLayers = 0;
    refCoeff = -1.0;
    break;
  case BT_PMC:
    numLayers = 0;
    refCoeff = +1.0;
    break;
  case BT_FREE_SPACE:
    numLayers = 0;
    refCoeff = 0.0;   
    break;
  case BT_PERIODIC:
    numLayers = 0;
    refCoeff = +1.0;
    break;
  case BT_PML:
    setPmlDefaults( &numLayers , &order , &n_eff , &refCoeff , &kmax );
    numScanned = sscanf( line , "%31s %31s %d %d %e %e %e" , name , typeStr , &numLayers , 
                         &order , &n_eff , &refCoeff , &kmax  );
    if( numScanned < 2 )
      assert( 0 );   
    break;
  case BT_MUR:
    numLayers = 0;
    refCoeff = 0.0;   
    break;
  case BT_SIBC:
    /* See if scattering parameters given. */
    numScanned = sscanf( line , "%31s %31s %e %e %e %e %e %e %e %e" , name , typeStr , 
                         &S_TM[0][0] , &S_TM[0][1] , &S_TM[1][0] , &S_TM[1][1] , 
                         &S_TE[0][0] , &S_TE[0][1] , &S_TE[1][0] , &S_TE[1][1] );
    if( numScanned == 10 )
    {
      /* Anisotropic SIBC. */;
    }
    else if( numScanned == 6 )    
    {
      /* Isotropic SIBC. */
      S_TE[0][0] = S_TM[0][0];
      S_TE[0][1] = S_TM[0][1];
      S_TE[1][0] = S_TM[1][0];
      S_TE[1][1] = S_TM[1][1];
    }
    else
    {
      /* See if external file. */
      numScanned = sscanf( line , "%31s %31s \"%1023[^\"]\"" , name , typeStr , fileName );
      if( numScanned != 3 )
      {
        message( MSG_LOG , 0 , "  Invalid/missing file name in boundary card:\n" );
        return false;
      }
    }
    break;
  default:
    assert( 0 );
    break;
  }

  addBoundary( name , type , numLayers , order , n_eff , refCoeff , kmax , fileName , S_TM , S_TE );

  return true;

}

/* Add boundary to lists. */
void addBoundary( char *name , BoundaryType type , int numLayers , int order , 
                  real n_eff , real refCoeff , real kmax , char *fileName , real S_TM[2][2] , real S_TE[2][2] )
{

  BoundaryItem *item = NULL;

  if( numBoundary == MAX_BOUNDARY )
    message( MSG_ERROR , 0 , "*** Error: Maximum number of boundaries exceeded!\n" );

  item = (BoundaryItem *) malloc( sizeof( BoundaryItem ) );
  if( !item )
    message( MSG_ERROR , 0 , "*** Error: Failed to allocate boundary!\n" );

  strncpy( item->name , name , TAG_SIZE );
  item->type = type;
  item->number = numBoundary;
  item->numLayers = numLayers;
  item->order = order;
  item->n_eff = n_eff;
  item->refCoeff = refCoeff;
  item->kmax = kmax;
  strncpy( item->fileName , fileName , PATH_SIZE );
  if( S_TM != NULL )
    for( int i = 0 ; i <= 1 ; i++ )
      for( int j = 0 ; j <= 1 ; j++ )
        item->S_TM[i][j] = S_TM[i][j];
  if( S_TE != NULL )
    for( int i = 0 ; i <= 1 ; i++ )
      for( int j = 0 ; j <= 1 ; j++ )
        item->S_TE[i][j] = S_TE[i][j];

  /* Add to list and hash. */
  DL_APPEND( boundaryList , item );
  HASH_ADD_STR( boundaryHash , name , item );
  HASH_ADD( hhint , boundaryNumberHash , number , sizeof( numBoundary ) , item );
  numBoundary++;
  isBoundaryType[type] = true;
  isBoundaryType[BT_UNDEFINED] = true;
  
  return;

}

/* Get boundary number from name. */
bool isBoundary( char *name , BoundaryIndex *number )
{

  BoundaryItem *item;

  HASH_FIND_STR( boundaryHash , name , item );
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

/* Initialise boundaries. */
/* Depends: initGrid */
void initBoundaries( void )
{

  BoundaryItem *item;
  
  message( MSG_LOG , 0 , "  Initialising boundaries...\n" );
  
  /* Determine update coefficient for each material. */
  DL_FOREACH( boundaryList , item ) 
  {
    switch( item->type )
    {
    case BT_SIBC:
#ifdef WITH_SIBC
      initSibcBoundary( item );
#else
      message( MSG_ERROR , 0 , "SIBCs not included in this executable!\n" );
#endif
      break;
    default:
      break;
    }
  }

  return;

}

/* Deallocate boundaries. */
void deallocBoundaries( void )
{

  BoundaryItem *item , *tmp;

  message( MSG_DEBUG1 , 0 , "Deallocating boundaries...\n" );
  
  /* Free the number hash. */
  HASH_ITER( hhint , boundaryNumberHash , item , tmp )
  {
    HASH_DELETE( hhint , boundaryNumberHash , item );
  }

  /* Free the name hash and boundaries. */
  HASH_ITER( hh , boundaryHash , item , tmp )
  {
#ifdef WITH_SIBC
    if( item->type == BT_SIBC )
      deallocSibcBoundary( item );
#endif
    HASH_DELETE( hh , boundaryHash , item );
    free( item );
  }
 
  return;

}

/* Report boundaries. */
void reportBoundaries( void )
{

  BoundaryItem *item;

  message( MSG_LOG , 0 , "  Number of boundaries: %lu\n" , numBoundary );

  DL_FOREACH( boundaryList , item ) 
  {
    message( MSG_DEBUG3 , 0 , "    Boundary #%lu: Name=%s Type=%s Layers=%d Order=%d n_eff=%e rho=%e kmax=%e file=%s\n" , 
             (unsigned long) item->number , item->name , BOUNDARY_TYPE_STR[item->type] , 
             item->numLayers , item->order , item->n_eff , item->refCoeff , item->kmax , item->fileName );
  }

  return;

}

/* Return true if there are boundaries of given type. */
bool thereAreBoundaries( BoundaryType type )
{
  
  return isBoundaryType[type];

}

/* Get boundary pointer. */
BoundaryItem *getBoundary( BoundaryIndex number )
{

  BoundaryItem *item;

  HASH_FIND( hhint , boundaryNumberHash , &number , sizeof( number ) , item );
  if( !item)
    assert( 0 );

  return item;

}

/* Get boundary types by number. */
BoundaryType getBoundaryType( BoundaryIndex number )
{

  BoundaryItem *item;

  HASH_FIND( hhint , boundaryNumberHash , &number , sizeof( number ) , item );
  if( !item)
    assert( 0 );

  return item->type;

}

/* Get pointer to boundary name by number. */
char *getBoundaryName( BoundaryIndex number )
{

  BoundaryItem *item;

  HASH_FIND( hhint , boundaryNumberHash , &number , sizeof( number ) , item );
  if( !item)
    assert( 0 );

  return item->name;

}

/* Get number of layers by number. */
int getBoundaryNumLayers( BoundaryIndex number )
{

  BoundaryItem *item;

  HASH_FIND( hhint , boundaryNumberHash , &number , sizeof( number ) , item );
  if( !item)
    assert( 0 );

  return item->numLayers;

}

/* Get reflection coefficient by number. */
real getBoundaryRefCoeff( BoundaryIndex number )
{

  BoundaryItem *item;

  HASH_FIND( hhint , boundaryNumberHash , &number , sizeof( number ) , item );
  if( !item)
    assert( 0 );

  return item->refCoeff;

}

/* Get external boundary parameters by boundary number. */
void getExternalBoundaryParams( BoundaryIndex number , int *order , real *n_eff , real *refCoeff , real *kmax )
{
  BoundaryItem *item;

  HASH_FIND( hhint , boundaryNumberHash , &number , sizeof( number ) , item );
  if( !item)
    assert( 0 );
  
  *order = item->order;
  *n_eff = item->n_eff;
  *refCoeff = item->refCoeff;
  *kmax = item->kmax;

  return;
}

/* Set the number of layers by number. */
void setBoundaryNumLayers( BoundaryIndex number , int numLayers )
{

  BoundaryItem *item;

  HASH_FIND( hhint , boundaryNumberHash , &number , sizeof( number ) , item );
  if( !item)
    assert( 0 );

  item->numLayers = numLayers;

  return;

}
