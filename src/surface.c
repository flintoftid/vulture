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

#include "surface.h"
#include "utlist.h"
#include "alloc_array.h"
#include "grid.h"
#include "message.h"
#include "boundary.h"
#include "medium.h"
#include "pml.h"
#include "bounding_box.h"
#include "gnuplot.h"
#include "gmsh.h"
#include "memory.h"

#include "mur.c"

#ifdef WITH_SIBC
  #include "sibc.h"
#endif
    
/* 
 * Private data.
 */

/* Number of interal surfaces. */
static SurfaceIndex numInternalSurface = 0;

/* Existance flag for internal surfaces of each boundary type, including undefined. */
static bool isInternalSurfaceType[NUM_BOUNDARY_TYPES+1] = { false };

/* List of internal surfaces. */
static SurfaceItem *surfaceList = NULL;

/* Number of exteral surfaces. */
static SurfaceIndex numExternalSurface = 0;

/* Existance flag for external surfaces of each boundary type, including undefined. */
static bool isExternalSurfaceType[NUM_BOUNDARY_TYPES+1] = { false };

/* Array of external surfaces. */
static SurfaceItem externalSurfaceList[6];

/* 
 * Private method interfaces. 
 */

void addExternalSurface( SurfaceIndex surfaceNumber , char *boundaryName , BoundaryIndex boundaryNumber , int orient , real angle );
void addInternalSurface( int mbbox[6] , char *boundaryName , BoundaryIndex boundaryNumber , int orient , real angle );

/*
 * Method Implementations.
 */

/* Add an external surface to the array. */
void addExternalSurface( SurfaceIndex surfaceNumber , char *boundaryName , BoundaryIndex boundaryNumber , int orient , real angle )
{

  externalSurfaceList[surfaceNumber].boundaryNumber = boundaryNumber;
  strncpy( externalSurfaceList[surfaceNumber].boundaryName , boundaryName , TAG_SIZE );
  getFaceOfBoundingBox( externalSurfaceList[surfaceNumber].mbbox , mbox , surfaceNumber );
  externalSurfaceList[surfaceNumber].orientation = orient;
  externalSurfaceList[surfaceNumber].angle = angle;
  numExternalSurface++;
  
  return;

}

/* Add an internal surface to the lists. */
void addInternalSurface( int mbbox[6] , char *boundaryName , BoundaryIndex boundaryNumber , int orient , real angle )
{

  SurfaceItem *item = NULL;

  if( numInternalSurface == MAX_SURFACE )
    message( MSG_ERROR , 0 , "*** Error: Maximum number of surfaces exceeded!\n" );

  item = (SurfaceItem *) malloc( sizeof( SurfaceItem ) );
  if( !item )
    message( MSG_ERROR , 0 , "*** Error: Failed to allocate surface.\n" );

  for( int i = XLO ; i <= ZHI ; i++ ) item->mbbox[i] = mbbox[i];
  item->boundaryNumber = boundaryNumber;
  strncpy( item->boundaryName , boundaryName , TAG_SIZE );
  item->orientation = orient;
  item->angle = angle;

  /* Add to list. */
  DL_APPEND( surfaceList , item );
  numInternalSurface++;

  return;

}

/* Parse external surfaces. */
bool parseBR( char *line )
{
  
  int intType[6];
  int numLayers = 0;
  int order = 0;
  real n_eff = 0.0;
  real refCoeff = 0.0;
  real kmax = 0.0;
  char fileName[PATH_SIZE] = "";
  BoundaryType type;
  
  if( sscanf( line , "%d %d %d %d %d %d" , &intType[XLO] , &intType[XHI] , &intType[YLO] , &intType[YHI] , &intType[ZLO] , &intType[ZHI] ) != 6 )
    return false;  

  message( MSG_WARN , 0 , "  BR is obsolete - please specify external boundaries using BT\n" );
  
  for( int boundary = XLO ; boundary <= ZHI ; boundary++ )
  {
    switch( intType[boundary] )
    {
    case -1:
      type = BT_PEC;
      numLayers = 0;
      refCoeff = -1.0;
      break;
    case 0:
      type = BT_PML;
      setPmlDefaults( &numLayers , &order , &n_eff , &refCoeff , &kmax );
      break;
    case +1:
      type = BT_PMC;
      numLayers = 0;
      refCoeff = 1.0;
      break;
    default:
      message( MSG_LOG , 0 , "  Invalid type, %d, for %s boundary\n" , intType[boundary] , FACE[boundary] );
      return false;
      break;
    }  

    addBoundary( FACE[boundary] , type , numLayers , order , n_eff , refCoeff , kmax , fileName , NULL , NULL );

  }

  /* This is needed by mesh renderers in gvulture.*/
  isExternalSurfaceType[BT_UNDEFINED] = true;

  return true;

}

/* Parse internal surfaces. */
bool parseTB( char *line )
{

  int numScanned = 0;
  char boundaryName[TAG_SIZE] = "";
  int mbbox[6];
  BoundaryIndex boundaryNumber = 0;
  int orient = 1;
  double angle = 0.0;

  numScanned = sscanf( line , "%d %d %d %d %d %d %31s %d %lf" , 
                       &mbbox[XLO] , &mbbox[XHI] , &mbbox[YLO] , &mbbox[YHI] , &mbbox[ZLO] , &mbbox[ZHI] , boundaryName , &orient , &angle );

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
  else if( bboxType( mbbox ) != BB_SURFACE )
  {
    message( MSG_LOG , 0 , "  Bounding box is not a surface!\n" );
    return false;
  }

  /* Check boundary exists. */
  if( !isBoundary( boundaryName , &boundaryNumber ) )
  {
    message( MSG_LOG , 0 , "  Boundary %s not defined in TB card\n" , boundaryName );
    return false;
  }

  if( orient != 1 && orient != -1 )
  {
    message( MSG_LOG , 0 , "  Invalid orientation %d in TB card\n" , orient );   
    return false;
  }

  if( angle < -180.0 || angle > +180.0 )
  {
    message( MSG_LOG , 0 , "  Invalid angle %g in TB card\n" , angle );  
    return false;
  }

  addInternalSurface( mbbox , boundaryName , boundaryNumber , orient , (real) angle );

  /* This is needed by mesh renderers in gvulture.*/
  isInternalSurfaceType[BT_UNDEFINED] = true;

  return true;

}

/* Initialise external surface parameters. */
/* Depends: initBoundary */
void initExternalSurfaceParameters( void )
{

  int numLayers;
  BoundaryIndex number;
  BoundaryType type;
  BoundaryType oppositeSurfaceType;
  int oppositeSurface;
  
  message( MSG_LOG , 0 , "  Initialising the external surface arrays...\n" );
  
  /* Check external boundaries all defined - add default if not. */
  for( int boundary = XLO ; boundary <= ZHI ; boundary++ )
  {
    if( !isBoundary( FACE[boundary] , &number ) )
    {
      int numLayers;
      int order;
      real n_eff;
      real refCoeff;
      real kmax;
      setPmlDefaults( &numLayers , &order , &n_eff , &refCoeff , &kmax );
      addBoundary( FACE[boundary] , BT_PML , numLayers , order , n_eff , refCoeff , kmax , "" , NULL , NULL );
      if( !isBoundary( FACE[boundary] , &number ) ) assert( 0 );
    }
    
    addExternalSurface( boundary , FACE[boundary] , number , 1 , 0.0 );
  }
 
 
 
  /* Validate external surface type consistency. */
  for( int surface = XLO ; surface <= ZHI ; surface++ )
  {

    numLayers = getBoundaryNumLayers( externalSurfaceList[surface].boundaryNumber );

    type = getBoundaryType( externalSurfaceList[surface].boundaryNumber );
    isExternalSurfaceType[type] = true;
    
    switch( type )
    {
    case BT_PEC:
      if( numLayers != 0 ) 
      {
        message( MSG_WARN , 0 , "*** Warning: PEC on %s surface has %d PML layers - reset to zero\n" , FACE[surface] , numLayers );
        setBoundaryNumLayers( externalSurfaceList[surface].boundaryNumber , 0 );
      }
      break;
    case BT_PMC:
      if( numLayers != 0 ) 
      {
        message( MSG_WARN , 0 , "*** Warning: PMC on %s surface has %d PML layers - reset to zero\n" , FACE[surface] , numLayers );
        setBoundaryNumLayers( externalSurfaceList[surface].boundaryNumber , 0 );
      }
      break;
    case BT_PERIODIC:
      if( numLayers != 0 ) 
      {
        message( MSG_WARN , 0 , "*** Warning: PERIODIC on %s surface has %d PML layers - reset to zero\n" , FACE[surface] , numLayers );
        setBoundaryNumLayers( externalSurfaceList[surface].boundaryNumber , 0 );
      }
      oppositeSurface = surface - 2 * ( surface % 2 ) + 1;
      oppositeSurfaceType = getBoundaryType( externalSurfaceList[oppositeSurface].boundaryNumber );
      if( oppositeSurfaceType != BT_PERIODIC ) 
        message( MSG_ERROR , 0 , "*** PERIODIC boundary on %s surface doesn't match that on surface %s\n" , FACE[surface] , FACE[oppositeSurface] );
      break;
    case BT_PML:
      if( numLayers < 1 )
        message( MSG_ERROR , 0 , "*** Warning: PML on %s surface has less than one (%d) layers\n" , FACE[surface] , numLayers );
      break;
    case BT_MUR: 
      if( numLayers != 0 ) 
      {
        message( MSG_WARN , 0 , "*** Warning: MUR on %s surface has %d PML layers - reset to zero\n" , FACE[surface] , numLayers );
        setBoundaryNumLayers( externalSurfaceList[surface].boundaryNumber , 0 );
      }       
      break;
    default:
      assert( 0 );
      break;
    }
  }

  return;

}

/* Initialise external surfaces. */
/* Depends: initInternalSurfaces, initMedia */
void initExternalSurfaces( void )
{

  message( MSG_LOG , 0 , "\nInitialising non-PEC/PMC external surfaces...\n\n" );

  /* Initialise external PEC/PMC surfaces. */
  initExternalPecPmcSurfaces();

  /* Initialise PML boundaries. */
  initPmlBoundaries();

  /* Initialise Mur boundaries. */
  initMurBoundaries();
  
  return;

}

/* Initialise internal surfaces. */
/* Depends: initExternalPecPmcSurfaces */
void initInternalSurfaces( void )
{

  SurfaceItem *item;
  SurfaceIndex numSibc = 0;
  BoundaryType type;
  
  message( MSG_LOG , 0 , "\nInitialising internal surfaces...\n\n" );

  DL_FOREACH( surfaceList , item ) 
  {
    
    type = getBoundaryType( item->boundaryNumber );
    isInternalSurfaceType[type] = true;

    switch( type )
    {
    case BT_PEC:
      offsetBoundingBox( item->gbbox , item->mbbox , gibox );
      message( MSG_DEBUG3 , 0 , "  Setting PEC surface medium#%lu on [%d,%d,%d,%d,%d,%d]/[%d,%d,%d,%d,%d,%d]\n" , 
               MT_PEC , item->mbbox[XLO] , item->mbbox[XHI] , item->mbbox[YLO] , item->mbbox[YHI] , 
	       item->mbbox[ZLO] , item->mbbox[ZHI] , item->gbbox[XLO] , item->gbbox[XHI] , 
	       item->gbbox[YLO] , item->gbbox[YHI] , item->gbbox[ZLO] , item->gbbox[ZHI] );
      /* Apply PEC surfaces using medium coefficients. */
      setMediumOnGrid( item->gbbox , MT_PEC , FACE_MASK_ALL );
      break;
    case BT_FREE_SPACE:
      offsetBoundingBox( item->gbbox , item->mbbox , gibox );
      message( MSG_DEBUG3 , 0 , "  Setting FREE_SPACE surface medium#%lu on [%d,%d,%d,%d,%d,%d]/[%d,%d,%d,%d,%d,%d]\n" , 
               MT_FREE_SPACE , item->mbbox[XLO] , item->mbbox[XHI] , item->mbbox[YLO] , item->mbbox[YHI] , 
               item->mbbox[ZLO] , item->mbbox[ZHI] , item->gbbox[XLO] , item->gbbox[XHI] , 
               item->gbbox[YLO] , item->gbbox[YHI] , item->gbbox[ZLO] , item->gbbox[ZHI] );
      setMediumOnGrid( item->gbbox , MT_FREE_SPACE , FACE_MASK_ALL );
      break;
    case BT_SIBC:
      offsetBoundingBox( item->gbbox , item->mbbox , gibox );
      message( MSG_DEBUG3 , 0 , "  Setting SIBC surface medium#%lu on [%d,%d,%d,%d,%d,%d]/[%d,%d,%d,%d,%d,%d]\n" , 
               MT_PEC , item->mbbox[XLO] , item->mbbox[XHI] , item->mbbox[YLO] , item->mbbox[YHI] , 
	       item->mbbox[ZLO] , item->mbbox[ZHI] , item->gbbox[XLO] , item->gbbox[XHI] , 
	       item->gbbox[YLO] , item->gbbox[YHI] , item->gbbox[ZLO] , item->gbbox[ZHI] );
      setMediumOnGrid( item->gbbox , MT_PEC , FACE_MASK_ALL );
      numSibc++;
      break;
    default:
      break;
    }
  }

  /* Initialise SIBCs. */ 
#ifdef WITH_SIBC
    initSibcSurfaces( numSibc , surfaceList );
#endif
  
  return;

}

/* Initialise PEC and PMC boundaries, including PEC backing of PML. */
/* Depends: initGrid, initMedia */
void initExternalPecPmcSurfaces( void )
{

  BoundaryType boundaryType;
  MediumType mediumType;
  int bbox[6];

  message( MSG_LOG , 0 , "\nInitialising PEC/PMC surfaces...\n\n" );

  /* Mesh bounding box. */
  for( int face = XLO ; face <= ZHI ; face++ )
    offsetBoundingBox( externalSurfaceList[face].gbbox , externalSurfaceList[face].mbbox , gibox ); 
  
  /* Set material type for each external grid surface - CAREFUL HERE! 
   * For non-PML boundaries the external grid and external mesh surfaces coincide.
   * This code sets the PEC backing of PML boundaries.
   * For PEC external mesh surfaces we must overwrite any material parameters from MB/TB/TW. 
   */
  for( int boundary = XLO ; boundary <= ZHI ; boundary++ )
  {
  
    /* Type of boundary. */
    boundaryType = getBoundaryType( externalSurfaceList[boundary].boundaryNumber );

    /* Get medium type associated with boundary type. */
    mediumType = BOUNDARY_MEDIUM_TYPE[boundaryType];
    
    /* Get surface - for PML this is the *back* PEC face. For others gibox and gobox coincide. */
    getFaceOfBoundingBox( bbox , gobox , boundary );

    message( MSG_LOG , 0 , "  %s (#=%lu) [%d,%d,%d,%d,%d,%d]: %s (#=%lu) -> %s\n" , FACE[boundary] , 
             (unsigned long)externalSurfaceList[boundary].boundaryNumber ,
             bbox[XLO] , bbox[XHI] , bbox[YLO] , bbox[YHI] , bbox[ZLO] , bbox[ZHI] ,           
             BOUNDARY_TYPE_STR[boundaryType] , 
             (unsigned long)boundaryType , MEDIUM_TYPE_STR[mediumType] );

    switch( mediumType )
    {
      case MT_PEC:
        setMediumOnGrid( bbox , mediumType , FACE_MASK_ALL );
        break;
      default:
        break;        
    }

  }

  return;

}

/* Report surfaces. */
void reportSurfaces( void )
{

  SurfaceItem *surfaceItem;
  SurfaceIndex counter = 0;

  message( MSG_LOG , 0 , "  Number of external surfaces: 6\n" );  

  for( int surface = XLO ; surface <= ZHI ; surface++ )
  {
    message( MSG_DEBUG3 , 0 , "    Surface #%d: Boundary=%s Boundary#=%lu BBOX=[%d,%d,%d,%d,%d,%d] orient=%d angle=%e\n" , 
             surface , externalSurfaceList[surface].boundaryName , (unsigned long) externalSurfaceList[surface].boundaryNumber ,
             externalSurfaceList[surface].mbbox[XLO] , externalSurfaceList[surface].mbbox[XHI] , externalSurfaceList[surface].mbbox[YLO] , 
             externalSurfaceList[surface].mbbox[YHI] , externalSurfaceList[surface].mbbox[ZLO] , externalSurfaceList[surface].mbbox[ZHI] ,
             externalSurfaceList[surface].orientation , externalSurfaceList[surface].angle );

  }

  message( MSG_LOG , 0 , "  Number of internal surfaces: %lu\n" , numInternalSurface );
  DL_FOREACH( surfaceList , surfaceItem ) 
  {
    message( MSG_DEBUG3 , 0 , "    Surface #%lu: Boundary=%s Boundary#=%lu BBOX=[%d,%d,%d,%d,%d,%d] orient=%d angle=%e\n" , 
             (unsigned long) counter , surfaceItem->boundaryName , (unsigned long) surfaceItem->boundaryNumber ,
             surfaceItem->mbbox[XLO] , surfaceItem->mbbox[XHI] , surfaceItem->mbbox[YLO] , 
             surfaceItem->mbbox[YHI] , surfaceItem->mbbox[ZLO] , surfaceItem->mbbox[ZHI] ,
             surfaceItem->orientation , surfaceItem->angle );
    counter++;
  }

  return;

}

/* Return true if there are internal surfaces of given boundary type. */
bool thereAreInternalSurfaces( BoundaryType type )
{
  
  return isInternalSurfaceType[type];

}

/* Return true if there are external surfaces of given boundary type. */
bool thereAreExternalSurfaces( BoundaryType type )
{
  
  return isExternalSurfaceType[type];

}

/* Update external surface E fields. */
void updateExternalSurfacesEfield( void )
{

  /* Update E field in PML regions. */
  updatePmlEfield();
  
  /* Update E field in ghost regions. */
  //updateGhostEfield();
  
  /* Update E field on Mur boundaries. */
  updateMurEfield();
  
  return;
}

/* Update external surface H fields. */
void updateExternalSurfacesHfield( void )
{

  /* Update H field in PML regions. */
  updatePmlHfield();
  
  /* Update H field in ghost regions. */
  //updateGhostHfield();

  /* Update H field on Mur boundaries. */
  updateMurHfield();
  
  return;

}

/* Update internal surface E fields. */
void updateInternalSurfacesEfield( void )
{

#ifdef WITH_SIBC
  updateSibcSurfacesEfield();
#endif
  
  return;

}

/* Update internal surface H fields. */
void updateInternalSurfacesHfield( void )
{

#ifdef WITH_SIBC
  updateSibcSurfacesHfield();
#endif
  
  return;

}

/* Update ghost electric fields. */
/* Used by observers when averaging on external surfaces. */
void updateGhostEfield( void )
{

  int i , j , k;

  /* XLO boundary. */
  if( outerSurfaceType( XLO ) == BT_PEC )
    for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
        Ex[gobox[XLO]-1][j][k] = Ex[gobox[XLO]][j][k];
  else if( outerSurfaceType( XLO ) == BT_PMC )
    for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
        Ex[gobox[XLO]-1][j][k] = -Ex[gobox[XLO]][j][k];
  else if( outerSurfaceType( XLO ) == BT_PERIODIC )
    for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
        Ex[gobox[XLO]-1][j][k] = Ex[gobox[XHI]-1][j][k];
 
  /* XHI boundary. */
  if( outerSurfaceType( XHI ) == BT_PEC )
    for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
        Ex[gobox[XHI]][j][k] = Ex[gobox[XHI]-1][j][k];
  else if( outerSurfaceType( XHI ) == BT_PMC )
    for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
       Ex[gobox[XHI]][j][k] = -Ex[gobox[XHI]-1][j][k];
  else if( outerSurfaceType( XHI ) == BT_PERIODIC )
    for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
       Ex[gobox[XHI]][j][k] = Ex[gobox[XLO]][j][k];
        
  /* YLO boundary. */
  if( outerSurfaceType( YLO ) == BT_PEC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
        Ey[i][gobox[YLO]-1][k] = Ey[i][gobox[YLO]][k];
  else if( outerSurfaceType( YLO ) == BT_PMC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
        Ey[i][gobox[YLO]-1][k] = -Ey[i][gobox[YLO]][k];
  else if( outerSurfaceType( YLO ) == BT_PERIODIC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
        Ey[i][gobox[YLO]-1][k] = Ey[i][gobox[YHI]-1][k];

  /* YHI boundary. */
  if( outerSurfaceType( YHI ) == BT_PEC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
        Ey[i][gobox[YHI]][k] = Ey[i][gobox[YHI]-1][k];
  else if( outerSurfaceType( YHI ) == BT_PMC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
        Ey[i][gobox[YHI]][k] = -Ey[i][gobox[YHI]-1][k];
  else if( outerSurfaceType( YHI ) == BT_PERIODIC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
        Ey[i][gobox[YHI]][k] = Ey[i][gobox[YLO]][k];
      
  /* ZLO boundary. */
  if( outerSurfaceType( ZLO ) == BT_PEC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
        Ez[i][j][gobox[ZLO]-1] = Ez[i][j][gobox[ZLO]];
  else if( outerSurfaceType( ZLO ) == BT_PMC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
        Ez[i][j][gobox[ZLO]-1] = -Ez[i][j][gobox[ZLO]];
  else if( outerSurfaceType( ZLO ) == BT_PERIODIC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
        Ez[i][j][gobox[ZLO]-1] = Ez[i][j][gobox[ZHI]-1];
    
  /* ZHI boundary. */
  if( outerSurfaceType( ZHI ) == BT_PEC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
        Ez[i][j][gobox[ZHI]] = Ez[i][j][gobox[ZHI]-1];
  else if( outerSurfaceType( ZHI ) == BT_PMC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
        Ez[i][j][gobox[ZHI]] = -Ez[i][j][gobox[ZHI]-1];
  else if( outerSurfaceType( ZHI ) == BT_PERIODIC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
        Ez[i][j][gobox[ZHI]] = Ez[i][j][gobox[ZLO]];

  return;

}

/* Update external PMC boundary ghost fields. */
/* Tangential magnetic fields in ghost cells are set to minus */
/* the tangential magnetic field just inside the boundary. The */
/* fields on the boundary are updated by either the inner mesh time-stepping */
/* or PML time-stepping funcions. */
void updateGhostHfield( void )
{

  int i , j , k;

  /* XLO boundary. */
  if( outerSurfaceType( XLO ) == BT_PMC )
    for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
      {
        Hy[gobox[XLO]-1][j][k] = -Hy[gobox[XLO]][j][k];
        Hz[gobox[XLO]-1][j][k] = -Hz[gobox[XLO]][j][k];
      }
  else if( outerSurfaceType( XLO ) == BT_PERIODIC )
    for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
      {
        Hy[gobox[XLO]-1][j][k] = Hy[gobox[XHI]-1][j][k];
        Hz[gobox[XLO]-1][j][k] = Hz[gobox[XHI]-1][j][k];
      }
  else if( outerSurfaceType( XLO ) == BT_PEC )
    for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
      {
        Hy[gobox[XLO]-1][j][k] = Hy[gobox[XLO]][j][k];
        Hz[gobox[XLO]-1][j][k] = Hz[gobox[XLO]][j][k];
      }
    
  /* XHI boundary. */
  if( outerSurfaceType( XHI ) == BT_PMC )
    for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
      {
        Hy[gobox[XHI]][j][k] = -Hy[gobox[XHI]-1][j][k];
        Hz[gobox[XHI]][j][k] = -Hz[gobox[XHI]-1][j][k];
      }
  else if( outerSurfaceType( XHI ) == BT_PERIODIC )
    for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
      {
        Hy[gobox[XHI]][j][k] = Hy[gobox[XLO]][j][k];
        Hz[gobox[XHI]][j][k] = Hz[gobox[XLO]][j][k];
      }
  else if( outerSurfaceType( XHI ) == BT_PEC )
    for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
      {
        Hy[gobox[XHI]][j][k] = Hy[gobox[XHI]-1][j][k];
        Hz[gobox[XHI]][j][k] = Hz[gobox[XHI]-1][j][k];
      }

  /* YLO boundary. */
  if( outerSurfaceType( YLO ) == BT_PMC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
      {
        Hz[i][gobox[YLO]-1][k] = -Hz[i][gobox[YLO]][k];
        Hx[i][gobox[YLO]-1][k] = -Hx[i][gobox[YLO]][k];
      }
  else if( outerSurfaceType( YLO ) == BT_PERIODIC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
      {
        Hz[i][gobox[YLO]-1][k] = Hz[i][gobox[YHI]-1][k];
        Hx[i][gobox[YLO]-1][k] = Hx[i][gobox[YHI]-1][k];
      }
  else if( outerSurfaceType( YLO ) == BT_PEC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
      {
        Hz[i][gobox[YLO]-1][k] = Hz[i][gobox[YLO]][k];
        Hx[i][gobox[YLO]-1][k] = Hx[i][gobox[YLO]][k];
      }

  /* YHI boundary. */
  if( outerSurfaceType( YHI ) == BT_PMC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
      {
        Hz[i][gobox[YHI]][k] = -Hz[i][gobox[YHI]-1][k];
        Hx[i][gobox[YHI]][k] = -Hx[i][gobox[YHI]-1][k];
      }
  else if( outerSurfaceType( YHI ) == BT_PERIODIC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
      {
        Hz[i][gobox[YHI]][k] = Hz[i][gobox[YLO]][k];
        Hx[i][gobox[YHI]][k] = Hx[i][gobox[YLO]][k];
      }
  else if( outerSurfaceType( YHI ) == BT_PEC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
      {
        Hz[i][gobox[YHI]][k] = Hz[i][gobox[YHI]-1][k];
        Hx[i][gobox[YHI]][k] = Hx[i][gobox[YHI]-1][k];
      }

  /* ZLO boundary. */
  if( outerSurfaceType( ZLO ) == BT_PMC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      {
        Hx[i][j][gobox[ZLO]-1] = -Hx[i][j][gobox[ZLO]];
        Hy[i][j][gobox[ZLO]-1] = -Hy[i][j][gobox[ZLO]];
      }
  else if( outerSurfaceType( ZLO ) == BT_PERIODIC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      {
        Hx[i][j][gobox[ZLO]-1] = Hx[i][j][gobox[ZHI]-1];
        Hy[i][j][gobox[ZLO]-1] = Hy[i][j][gobox[ZHI]-1];
      }
  else if( outerSurfaceType( ZLO ) == BT_PEC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      {
        Hx[i][j][gobox[ZLO]-1] = Hx[i][j][gobox[ZLO]];
        Hy[i][j][gobox[ZLO]-1] = Hy[i][j][gobox[ZLO]];
      }

  /* ZHI boundary. */
  if( outerSurfaceType( ZHI ) == BT_PMC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      {
        Hx[i][j][gobox[ZHI]] = -Hx[i][j][gobox[ZHI]-1];
        Hy[i][j][gobox[ZHI]] = -Hy[i][j][gobox[ZHI]-1];
      }
  else if( outerSurfaceType( ZHI ) == BT_PERIODIC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      {
        Hx[i][j][gobox[ZHI]] = Hx[i][j][gobox[ZLO]];
        Hy[i][j][gobox[ZHI]] = Hy[i][j][gobox[ZLO]];
      }
  else if( outerSurfaceType( ZHI ) == BT_PEC )
    for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
      {
        Hx[i][j][gobox[ZHI]] = Hx[i][j][gobox[ZHI]-1];
        Hy[i][j][gobox[ZHI]] = Hy[i][j][gobox[ZHI]-1];
      }

  return;

}

/* Deallocate external boundary arrays. */
void deallocExternalSurfaces( void )
{

  message( MSG_DEBUG1 , 0 , "Deallocating external surfaces...\n" );

  deallocMurArrays();
  deallocPmlArrays();

  return;

}

/* Deallocate grid internal surfaces. */
void deallocInternalSurfaces( void )
{

  SurfaceItem *item , *tmp;

  message( MSG_DEBUG1 , 0 , "Deallocating internal surfaces...\n" );

  /* Free internal surface list. */
  DL_FOREACH_SAFE( surfaceList , item , tmp ) 
  {
    DL_DELETE( surfaceList , item );
    free( item );
  }

#ifdef WITH_SIBC
  deallocSibcSurfaces();
#endif
  
  return;

}

/* Draw external surfaces. */
void gnuplotExternalSurfaces( void )
{

  char externalFileName[] = "gnuplot-external.dat";
  FILE *outputFile;
  int bbox[6];

  outputFile = fopen( externalFileName , "w" );
  if( !outputFile )
    message( MSG_ERROR , 0 , "*** Error: Failed to external surface output file %s\n" , externalFileName );

  gnuplotProblemSize( outputFile , mbox );

  getFaceOfBoundingBox( bbox , mbox , ZLO );
  gnuplotBoundingBox( outputFile , bbox );
  getFaceOfBoundingBox( bbox , mbox , XLO );
  gnuplotBoundingBox( outputFile , bbox );
  getFaceOfBoundingBox( bbox , mbox , YLO );
  gnuplotBoundingBox( outputFile , bbox );

  getFaceOfBoundingBox( bbox , mbox , ZHI );
  gnuplotBoundingBox( outputFile , bbox );
  getFaceOfBoundingBox( bbox , mbox , XHI );
  gnuplotBoundingBox( outputFile , bbox );
  getFaceOfBoundingBox( bbox , mbox , YHI );
  gnuplotBoundingBox( outputFile , bbox );
 
  fclose( outputFile );

  return;

}

/* Draw internal surfaces. */
void gnuplotInternalSurfaces( void )
{

  char surfaceFileName[] = "gnuplot-surface.dat";
  FILE *outputFile;
  SurfaceItem *item;

  outputFile = fopen( surfaceFileName , "w" );
  if( !outputFile )
    message( MSG_ERROR , 0 , "*** Error: Failed to open surface output file %s\n" , surfaceFileName );

  gnuplotProblemSize( outputFile , mbox );

  DL_FOREACH( surfaceList , item ) 
  {
    gnuplotBoundingBox( outputFile , item->mbbox );
  }

  fclose( outputFile );

  return;

}

/* Draw external surfaces. */
void gmshExternalSurfaces( void )
{

  int bbox[6];
  int step[3] = { 1 , 1 , 1 };
  unsigned long entityNumber;

  entityNumber = gmshGetEntityNumber();
  getFaceOfBoundingBox( bbox , mbox , ZLO );
  gmshAddEntity( entityNumber , BB_SURFACE , "BT_ZLO" , bbox , step );
  entityNumber = gmshGetEntityNumber();
  getFaceOfBoundingBox( bbox , mbox , XLO );
  gmshAddEntity( entityNumber , BB_SURFACE , "BT_XLO" , bbox , step );
  entityNumber = gmshGetEntityNumber();
  getFaceOfBoundingBox( bbox , mbox , YLO );
  gmshAddEntity( entityNumber , BB_SURFACE , "BT_YLO" , bbox , step );

  entityNumber = gmshGetEntityNumber();
  getFaceOfBoundingBox( bbox , mbox , ZHI );
  gmshAddEntity( entityNumber , BB_SURFACE , "BT_ZHI" , bbox , step );
  entityNumber = gmshGetEntityNumber();
  getFaceOfBoundingBox( bbox , mbox , XHI );
  gmshAddEntity( entityNumber , BB_SURFACE , "BT_XHI" , bbox , step );
  entityNumber = gmshGetEntityNumber();
  getFaceOfBoundingBox( bbox , mbox , YHI );
  gmshAddEntity( entityNumber , BB_SURFACE , "BT_YHI" , bbox , step );

  return;

}

/* Draw internal surfaces. */
void gmshInternalSurfaces( void )
{

  int step[3] = { 1 , 1 , 1 };
  SurfaceItem *item;
  unsigned long entityNumber;
  char name[GMSH_NAME_LENGTH];
    
  DL_FOREACH( surfaceList , item ) 
  {
    entityNumber = gmshGetEntityNumber();
    snprintf( name , GMSH_NAME_LENGTH - 1 , "BT_%s" , getBoundaryName( item->boundaryNumber ) );
    gmshAddEntity( entityNumber , BB_SURFACE , name , item->mbbox , step );    
  }
  
  return;

}

/* Check if material arrays on boundaries still have their designated */
/* values after all initialisation code is complete. */
void checkExternalSurfaces( void )
{

  int obbox[6];
  int ibbox[6];
  
  message( MSG_LOG , 0 , "\nChecking boundary material parameters...\n\n" );

  for( int boundary = XLO ; boundary <= ZHI ; boundary++ )
  {

    getFaceOfBoundingBox( obbox , gobox , boundary );
    getFaceOfBoundingBox( ibbox , gibox , boundary );
    
    switch( outerSurfaceType( boundary ) )
    {
    case BT_PEC:
    case BT_PML:
      /* PEC and PML backing should still be PEC. */
      checkMediumOnGrid( obbox , MT_PEC );
      break;
    case BT_MUR:
      /* Mur only valid with free-soace at the moment. */ 
      checkMediumOnGrid( ibbox , MT_FREE_SPACE );
      break;
    case BT_PMC:
      /* Should be OK with any material. */
      break;
    default:
      assert( 0 );
      break;
    }
  }

  return;

}

/* Get external boundary type for given external surface. */
BoundaryType outerSurfaceType( MeshFace face )
{  
  return getBoundaryType( externalSurfaceList[face].boundaryNumber );
}

/* Get external boundary number of layer for given external surface. */
int outerSurfaceNumLayers( MeshFace face )
{
  return getBoundaryNumLayers( externalSurfaceList[face].boundaryNumber );
}

/* Get external boundary reflection coefficient for given external surface. */
real outerSurfaceReflectCoeff( MeshFace face )
{
  return getBoundaryRefCoeff( externalSurfaceList[face].boundaryNumber );
}

void getOuterSurfaceParams( MeshFace face , int *order , real *n_eff , real *refCoeff , real *kmax )
{
 
  getExternalBoundaryParams( externalSurfaceList[face].boundaryNumber , order , n_eff , refCoeff , kmax );

  return;
}

/* Determine if edge is on a PMC. */
bool isPmcEdge( CoordAxis direction , int index )
{

  switch( direction )
  {
    case XDIR:
      return ( ( ( outerSurfaceType( XLO ) == BT_PMC || outerSurfaceType( XLO ) == BT_PERIODIC ) && index == gibox[XLO] ) || 
                ( ( outerSurfaceType( XHI ) == BT_PMC || outerSurfaceType( XHI ) == BT_PERIODIC ) && index == gibox[XHI] ) );
      break;
    case YDIR:
      return ( ( ( outerSurfaceType( YLO ) == BT_PMC || outerSurfaceType( YLO ) == BT_PERIODIC ) && index == gibox[YLO] ) || 
                ( ( outerSurfaceType( YHI ) == BT_PMC || outerSurfaceType( YHI ) == BT_PERIODIC ) && index == gibox[YHI] ) );
      break;
    case ZDIR:
      return ( ( ( outerSurfaceType( ZLO ) == BT_PMC || outerSurfaceType( ZLO ) == BT_PERIODIC ) && index == gibox[ZLO] ) || 
                ( ( outerSurfaceType( ZHI ) == BT_PMC || outerSurfaceType( ZHI ) == BT_PERIODIC ) && index == gibox[ZHI] ) );
      break;
    default:
      assert( false );
      return false;
      break;
  }

}

