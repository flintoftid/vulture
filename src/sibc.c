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
#include <string.h>
#include <math.h>

#include "surface.h"
#include "utlist.h"
#include "alloc_array.h"
#include "grid.h"
#include "message.h"
#include "boundary.h"
#include "bounding_box.h"
#include "memory.h"
#include "physical.h"
#include "filter.h"

/* 
 * SIBC class. 
 */

typedef struct Sibc_t {

  int gbbox[6];                      // Bounding box on grid.
  int orientation;                   // Orientation of boundary on surface.
  CoordAxis normal;                  // Normal direction. 
  real cosa;                         // Cosine of angle of boundary on surface.
  real sina;                         // Sine of angle of boundary on surface.
  BoundaryItem *boundary;            // Pointer to boundary model.
  yfRecConvStateM ***rcm_s;          // Array of SIBC states.
  real ****Etan;                     // Tangential fields on each face.
  bool ****isAdjA;                   // Face adjacency array, side a.
  bool ****isAdjB;                   // Face adjacency array, side b.

} SibcItem;

/* 
 * Private data.
 */

/* Number of internal SIBC surfaces. */
static SurfaceIndex numSibcSurface = 0;

/* Internal SIBC surface array. */
static SibcItem *sibcArray = NULL;         

/* Orientation matrices. */
static real Acp[4][4] = { { 0.0 , 0.0 , 1.0 , 0.0 } ,
                          { 0.0 , 0.0 , 0.0 ,-1.0 } ,
                          {-1.0 , 0.0 , 0.0 , 0.0 } ,
                          { 0.0 , 1.0 , 0.0 , 0.0 } };

static real Acm[4][4] = { { 0.0 , 0.0 , 0.0 ,-1.0 } ,
                          { 0.0 , 0.0 , 1.0 , 0.0 } ,
                          { 0.0 , 1.0 , 0.0 , 0.0 } ,
                          {-1.0 , 0.0 , 0.0 , 0.0 } };

static real Asp[4][4] = { {-1.0 , 0.0 , 0.0 , 0.0 } ,
                          { 0.0 , 1.0 , 0.0 , 0.0 } ,
                          { 0.0 , 0.0 ,-1.0 , 0.0 } ,
                          { 0.0 , 0.0 , 0.0 , 1.0 } };

static real Asm[4][4] = { { 0.0 , 1.0 , 0.0 , 0.0 } ,
                          {-1.0 , 0.0 , 0.0 , 0.0 } ,
                          { 0.0 , 0.0 , 0.0 , 1.0 } ,
                          { 0.0 , 0.0 ,-1.0 , 0.0 } };

static real Bcp[4][4] = { { 1.0 , 0.0 , 0.0 , 0.0 } ,
                          { 0.0 , 1.0 , 0.0 , 0.0 } ,
                          { 0.0 , 0.0 , 1.0 , 0.0 } ,
                          { 0.0 , 0.0 , 0.0 , 1.0 } };

static real Bcm[4][4] = { { 0.0 , 1.0 , 0.0 , 0.0 } ,
                          { 1.0 , 0.0 , 0.0 , 0.0 } ,
                          { 0.0 , 0.0 , 0.0 , 1.0 } ,
                          { 0.0 , 0.0 , 1.0 , 0.0 } };

static real Bsp[4][4] = { { 0.0 , 0.0 ,-1.0 , 0.0 } ,
                          { 0.0 , 0.0 , 0.0 ,-1.0 } ,
                          { 1.0 , 0.0 , 0.0 , 0.0 } ,
                          { 0.0 , 1.0 , 0.0 , 0.0 } };

static real Bsm[4][4] = { { 0.0 , 0.0 , 0.0 ,-1.0 } ,
                          { 0.0 , 0.0 ,-1.0 , 0.0 } ,
                          { 0.0 , 1.0 , 0.0 , 0.0 } ,
                          { 1.0 , 0.0 , 0.0 , 0.0 } };

/* 
 * Private method interfaces. 
 */

void matMulVector( real y[4] , real A[4][4] , real x[4] );
void matLinearComb( real C[4][4] , real a , real A[4][4] , real b , real B[4][4] );
void tportS2Z( real Z[2][2] , real S[2][2] );
bool isPassiveS( real S[2][2] );
void setSibcFace( bool ****isSibcFace , int gbbox[6] , CoordAxis dir , bool value );

/*
 * Method Implementations.
 */

/* Intialise a single SIBC boundary type. */ 
void initSibcBoundary( BoundaryItem *item )
{

  yfPoleResidueM prm;
  real dt;
  real Z_TM[2][2];
  real Z_TE[2][2];

  dt = getGridTimeStep();  
  
  if( strlen( item->fileName ) == 0 )
  {
    if( !isPassiveS( item->S_TM ) || !isPassiveS( item->S_TE ) )
      message( MSG_ERROR , 0 , "  Passivity violation for SIBC boundary %s\n", item->name );
    tportS2Z( Z_TM , item->S_TM );
    tportS2Z( Z_TE , item->S_TE );
    item->prm = yfAllocPoleResidueM( 4 , 4 );
    for( int i = 0 ; i < 2 ; i++ )
    {
      for( int j = 0 ; j < 2 ; j++ )
      {  
        item->prm.pr[i][j]     = yfAllocPoleResidue( 0 , Z_TM[i][j] , NULL , NULL );
        item->prm.pr[i+2][j]   = yfAllocPoleResidue( 0 , 0.0 , NULL , NULL );
        item->prm.pr[i][j+2]   = yfAllocPoleResidue( 0 , 0.0 , NULL , NULL );
        item->prm.pr[i+2][j+2] = yfAllocPoleResidue( 0 , Z_TE[i][j] , NULL , NULL );
      }
    }       
  }
  else
  {
    /* Read in pole residue matrix. */
    prm = yfReadPoleResidueM( item->fileName );
    if( prm.m == 4 && prm.n == 4 )
    {    
      /* Anisotropic SIBC. */
      message( MSG_DEBUG3 , 0 , "  Setting anisotropic SIBC boundary %s RC model from %s\n", item->name, item->fileName );
      item->prm = prm;
    }
    else if( prm.m == 2 && prm.n == 2 )
    {
      /* Isotropic SIBC. */
      message( MSG_DEBUG3 , 0 , "  Setting isotropic SIBC boundary %s RC model from %s\n", item->name, item->fileName );
      item->prm = yfAllocPoleResidueM( 4 , 4 );
      for( int i = 0 ; i < 2 ; i++ )
      {
        for( int j = 0 ; j < 2 ; j++ )
        {  
          item->prm.pr[i][j]     = yfAllocPoleResidue( prm.pr[i][j].numPoles , prm.pr[i][j].asymp , prm.pr[i][j].poles , prm.pr[i][j].residues );
          item->prm.pr[i+2][j]   = yfAllocPoleResidue( 0 , 0.0 , NULL , NULL );
          item->prm.pr[i][j+2]   = yfAllocPoleResidue( 0 , 0.0 , NULL , NULL );
          item->prm.pr[i+2][j+2] = yfAllocPoleResidue( prm.pr[i][j].numPoles , prm.pr[i][j].asymp , prm.pr[i][j].poles , prm.pr[i][j].residues );
        }
      }  
      yfDeallocPoleResidueM( prm );
    }
    else
    {
      message( MSG_ERROR , 0 , "  Pole-residue model must be 2x2 or 4x4!\n" );
    }
  }

  /* Determine recursive convolution coefficients. */
  item->rcm = yfPoleResidueM2RecConvM( item->prm , (double)dt );

  #ifdef DEBUG
    yfPrintPoleResidueM( item->prm );
    yfPrintRecConvM( item->rcm );
  #endif
        
  return;

}

/* Deallocate SIBC boundary data. */
void deallocSibcBoundary( BoundaryItem *item )
{
  
  yfDeallocRecConvM( item->rcm );
  yfDeallocPoleResidueM( item->prm );

  return;
  
}

/* Initialise internal SIBC surfaces. */
void initSibcSurfaces( SurfaceIndex number , SurfaceItem *surfaceList )
{

  int i , j , k;
  int ii , jj , kk;
  SurfaceItem *item;
  SurfaceIndex surface;
  unsigned long bytes;
  
  message( MSG_LOG , 0 , "\nInitialising SIBC surfaces...\n\n" );

  numSibcSurface = number;

  /* Allocate SIBC surfaces. */
  message( MSG_DEBUG1 , 0 , "  Allocating SIBC surface array\n" );
  sibcArray = allocArray( &bytes , sizeof( SibcItem ) , 1 , numSibcSurface );
  memory.surfaces += bytes;

  /* Temporary SIBC face utilisation array. */
  bool ****isSibcFace;
  message( MSG_DEBUG1 , 0 , "  Allocating SIBC utilisation array\n" );
  isSibcFace = allocArray( &bytes , sizeof( bool ) , 4 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] , 3 );
  setSibcFace( isSibcFace , gobox , XDIR , false );
  setSibcFace( isSibcFace , gobox , YDIR , false );
  setSibcFace( isSibcFace , gobox , ZDIR , false );

  /* Set up SIBC surfaces. */
  surface = 0;
  DL_FOREACH( surfaceList , item ) 
  {
    switch( getBoundaryType( item->boundaryNumber ) )
    {
    case BT_SIBC:
      sibcArray[surface].boundary = getBoundary( item->boundaryNumber );
      sibcArray[surface].orientation = item->orientation;
      sibcArray[surface].normal = bboxDirection( item->mbbox );
      sibcArray[surface].cosa = cos( pi * item->angle / 180.0 );
      sibcArray[surface].sina = sin( pi * item->angle / 180.0 );
      offsetBoundingBox( sibcArray[surface].gbbox , item->mbbox , gibox );
      /* CARE! Add one in normal direction for convenient management of filter arrays. */
      switch( sibcArray[surface].normal )
      {
      case XDIR:
        sibcArray[surface].gbbox[XHI] = sibcArray[surface].gbbox[XHI] + 1;
        break;
      case YDIR:
        sibcArray[surface].gbbox[YHI] = sibcArray[surface].gbbox[YHI] + 1;
        break;
      case ZDIR:
        sibcArray[surface].gbbox[ZHI] = sibcArray[surface].gbbox[ZHI] + 1;
        break;
      default:
        assert( 0 );
        break;
      }
      
      message( MSG_DEBUG3 , 0 , "  Setting SIBC type #%lu on [%d,%d,%d,%d,%d,%d]/[%d,%d,%d,%d,%d,%d]: orient=%2d, norm=%s, angle=%e [deg.]\n" ,
        item->boundaryNumber ,
        item->mbbox[XLO] , item->mbbox[XHI] , item->mbbox[YLO] , item->mbbox[YHI] , 
        item->mbbox[ZLO] , item->mbbox[ZHI] , sibcArray[surface].gbbox[XLO] , sibcArray[surface].gbbox[XHI] , 
        sibcArray[surface].gbbox[YLO] , sibcArray[surface].gbbox[YHI] , sibcArray[surface].gbbox[ZLO] , 
        sibcArray[surface].gbbox[ZHI] , sibcArray[surface].orientation , AXIS[sibcArray[surface].normal] , item->angle );
      
      /* Allocate arrays for filter states and tangential electric fields. This would be more complicated if we hadn't added one in normal direction above. */
      sibcArray[surface].rcm_s = allocArray( &bytes , sizeof( yfRecConvStateM ) , 3 , sibcArray[surface].gbbox[XHI] - sibcArray[surface].gbbox[XLO] , 
                                                                                      sibcArray[surface].gbbox[YHI] - sibcArray[surface].gbbox[YLO] , 
                                                                                      sibcArray[surface].gbbox[ZHI] - sibcArray[surface].gbbox[ZLO] );
      memory.surfaces += bytes;
      sibcArray[surface].Etan = allocArray( &bytes , sizeof( real ) , 4 , sibcArray[surface].gbbox[XHI] - sibcArray[surface].gbbox[XLO] , 
                                                                          sibcArray[surface].gbbox[YHI] - sibcArray[surface].gbbox[YLO] , 
                                                                          sibcArray[surface].gbbox[ZHI] - sibcArray[surface].gbbox[ZLO] , 4 );
      memory.surfaces += bytes;

      sibcArray[surface].isAdjA = allocArray( &bytes , sizeof( bool ) , 4 , sibcArray[surface].gbbox[XHI] - sibcArray[surface].gbbox[XLO] , 
                                                                            sibcArray[surface].gbbox[YHI] - sibcArray[surface].gbbox[YLO] , 
                                                                            sibcArray[surface].gbbox[ZHI] - sibcArray[surface].gbbox[ZLO] , 4 );
      memory.surfaces += bytes;
        
      sibcArray[surface].isAdjB = allocArray( &bytes , sizeof( bool ) , 4 , sibcArray[surface].gbbox[XHI] - sibcArray[surface].gbbox[XLO] , 
                                                                            sibcArray[surface].gbbox[YHI] - sibcArray[surface].gbbox[YLO] , 
                                                                            sibcArray[surface].gbbox[ZHI] - sibcArray[surface].gbbox[ZLO] , 4 );
      memory.surfaces += bytes;
      
      /* Loop over *faces* of surface - note limits. */
      for( i = sibcArray[surface].gbbox[XLO] , ii = 0 ; i < sibcArray[surface].gbbox[XHI] ; i++ , ii++ )
        for( j = sibcArray[surface].gbbox[YLO] , jj = 0 ; j < sibcArray[surface].gbbox[YHI] ; j++ , jj++ )
          for( k = sibcArray[surface].gbbox[ZLO] , kk = 0 ; k < sibcArray[surface].gbbox[ZHI] ; k++ , kk++ )
          {
            sibcArray[surface].rcm_s[ii][jj][kk] = yfAllocRecConvStateM( sibcArray[surface].boundary->rcm );
            sibcArray[surface].Etan[ii][jj][kk][0] = 0.0;
            sibcArray[surface].Etan[ii][jj][kk][1] = 0.0;
            sibcArray[surface].Etan[ii][jj][kk][2] = 0.0;
            sibcArray[surface].Etan[ii][jj][kk][3] = 0.0;
            sibcArray[surface].isAdjA[ii][jj][kk][0] = false;
            sibcArray[surface].isAdjA[ii][jj][kk][1] = false;
            sibcArray[surface].isAdjA[ii][jj][kk][2] = false;
            sibcArray[surface].isAdjA[ii][jj][kk][3] = false;
            sibcArray[surface].isAdjB[ii][jj][kk][0] = false;
            sibcArray[surface].isAdjB[ii][jj][kk][1] = false;
            sibcArray[surface].isAdjB[ii][jj][kk][2] = false;
            sibcArray[surface].isAdjB[ii][jj][kk][3] = false;            
          }

      setSibcFace( isSibcFace , sibcArray[surface].gbbox , sibcArray[surface].normal , true );
        
      surface++;
          
      break;
    default:
      break;
    }
  }

  /* Set adjacency flags for all SIBC faces. */
  surface = 0;
  DL_FOREACH( surfaceList , item ) 
  {
    switch( getBoundaryType( item->boundaryNumber ) )
    {
    case BT_SIBC:   
      /* Loop over *faces* of surface - note limits. */
      switch( sibcArray[surface].normal )
      {
      case XDIR:
        for( i = sibcArray[surface].gbbox[XLO] , ii = 0 ; i < sibcArray[surface].gbbox[XHI] ; i++ , ii++ )
          for( j = sibcArray[surface].gbbox[YLO] , jj = 0 ; j < sibcArray[surface].gbbox[YHI] ; j++ , jj++ )
            for( k = sibcArray[surface].gbbox[ZLO] , kk = 0 ; k < sibcArray[surface].gbbox[ZHI] ; k++ , kk++ )
            {
              sibcArray[surface].isAdjA[ii][jj][kk][0] = isSibcFace[i-1][j][k][YDIR];
              sibcArray[surface].isAdjA[ii][jj][kk][1] = isSibcFace[i-1][j+1][k][YDIR];
              sibcArray[surface].isAdjA[ii][jj][kk][2] = isSibcFace[i-1][j][k][ZDIR];
              sibcArray[surface].isAdjA[ii][jj][kk][3] = isSibcFace[i-1][j][k+1][ZDIR];
              sibcArray[surface].isAdjB[ii][jj][kk][0] = isSibcFace[i][j][k][YDIR];
              sibcArray[surface].isAdjB[ii][jj][kk][1] = isSibcFace[i][j+1][k][YDIR];
              sibcArray[surface].isAdjB[ii][jj][kk][2] = isSibcFace[i][j][k][ZDIR];
              sibcArray[surface].isAdjB[ii][jj][kk][3] = isSibcFace[i][j+1][k+1][ZDIR];       
            }
        break;
      case YDIR:
        for( i = sibcArray[surface].gbbox[XLO] , ii = 0 ; i < sibcArray[surface].gbbox[XHI] ; i++ , ii++ )
          for( j = sibcArray[surface].gbbox[YLO] , jj = 0 ; j < sibcArray[surface].gbbox[YHI] ; j++ , jj++ )
            for( k = sibcArray[surface].gbbox[ZLO] , kk = 0 ; k < sibcArray[surface].gbbox[ZHI] ; k++ , kk++ )
            {
              sibcArray[surface].isAdjA[ii][jj][kk][0] = isSibcFace[i][j-1][k][ZDIR];
              sibcArray[surface].isAdjA[ii][jj][kk][1] = isSibcFace[i][j-1][k+1][ZDIR];
              sibcArray[surface].isAdjA[ii][jj][kk][2] = isSibcFace[i][j-1][k][XDIR];
              sibcArray[surface].isAdjA[ii][jj][kk][3] = isSibcFace[i+1][j-1][k][XDIR];
              sibcArray[surface].isAdjB[ii][jj][kk][0] = isSibcFace[i][j][k][ZDIR];
              sibcArray[surface].isAdjB[ii][jj][kk][1] = isSibcFace[i][j][k+1][ZDIR];
              sibcArray[surface].isAdjB[ii][jj][kk][2] = isSibcFace[i][j][k][XDIR];
              sibcArray[surface].isAdjB[ii][jj][kk][3] = isSibcFace[i+1][j][k+1][XDIR];       
            }
        break;
      case ZDIR:
        for( i = sibcArray[surface].gbbox[XLO] , ii = 0 ; i < sibcArray[surface].gbbox[XHI] ; i++ , ii++ )
          for( j = sibcArray[surface].gbbox[YLO] , jj = 0 ; j < sibcArray[surface].gbbox[YHI] ; j++ , jj++ )
            for( k = sibcArray[surface].gbbox[ZLO] , kk = 0 ; k < sibcArray[surface].gbbox[ZHI] ; k++ , kk++ )
            {
              sibcArray[surface].isAdjA[ii][jj][kk][0] = isSibcFace[i][j][k-1][XDIR];
              sibcArray[surface].isAdjA[ii][jj][kk][1] = isSibcFace[i+1][j][k-1][XDIR];
              sibcArray[surface].isAdjA[ii][jj][kk][2] = isSibcFace[i][j][k-1][YDIR];
              sibcArray[surface].isAdjA[ii][jj][kk][3] = isSibcFace[i][j+1][k-1][YDIR];
              sibcArray[surface].isAdjB[ii][jj][kk][0] = isSibcFace[i][j][k][XDIR];
              sibcArray[surface].isAdjB[ii][jj][kk][1] = isSibcFace[i+1][j][k][XDIR];
              sibcArray[surface].isAdjB[ii][jj][kk][2] = isSibcFace[i][j][k][YDIR];
              sibcArray[surface].isAdjB[ii][jj][kk][3] = isSibcFace[i+1][j+1][k][YDIR];
            }
        break;
      default:
        assert( 0 );
        break;
      }
         
      surface++;
          
      break;
    default:
      break;
    }
  }
  
  /* Deallocate temporary adjacency array. */
  deallocArray( isSibcFace , 4 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] , 3 );
      
  return;

}

/* Deallocate SIBC surfaces. */
void deallocSibcSurfaces( void )
{

  int i , j , k;
  int ii , jj , kk;

  message( MSG_DEBUG1 , 0 , "Deallocating SIBC surfaces...\n" );

  for( SurfaceIndex surface = 0 ; surface < numSibcSurface ; surface++ )
  {

    /* Deallocate RC filter states. */ 
    for( i = sibcArray[surface].gbbox[XLO] , ii = 0 ; i < sibcArray[surface].gbbox[XHI] ; i++ , ii++ )
      for( j = sibcArray[surface].gbbox[YLO] , jj = 0 ; j < sibcArray[surface].gbbox[YHI] ; j++ , jj++ )
        for( k = sibcArray[surface].gbbox[ZLO] , kk = 0 ; k < sibcArray[surface].gbbox[ZHI] ; k++ , kk++ )
          yfDeallocRecConvStateM( sibcArray[surface].rcm_s[ii][jj][kk] );

    /* Deallocate face adjacency arrays. */
    deallocArray( sibcArray[surface].isAdjA , 4 , sibcArray[surface].gbbox[XHI] - sibcArray[surface].gbbox[XLO] , 
                                                  sibcArray[surface].gbbox[YHI] - sibcArray[surface].gbbox[YLO] , 
                                                  sibcArray[surface].gbbox[ZHI] - sibcArray[surface].gbbox[ZLO] , 4 );

    deallocArray( sibcArray[surface].isAdjB , 4 , sibcArray[surface].gbbox[XHI] - sibcArray[surface].gbbox[XLO] , 
                                                  sibcArray[surface].gbbox[YHI] - sibcArray[surface].gbbox[YLO] , 
                                                  sibcArray[surface].gbbox[ZHI] - sibcArray[surface].gbbox[ZLO] , 4 );
    
    /* Deallocate face state variables. */
    deallocArray( sibcArray[surface].Etan , 4 , sibcArray[surface].gbbox[XHI] - sibcArray[surface].gbbox[XLO] , 
                                                sibcArray[surface].gbbox[YHI] - sibcArray[surface].gbbox[YLO] , 
                                                sibcArray[surface].gbbox[ZHI] - sibcArray[surface].gbbox[ZLO] , 4 );

    /* Deallocate RC filters. */
    deallocArray( sibcArray[surface].rcm_s , 3 , sibcArray[surface].gbbox[XHI] - sibcArray[surface].gbbox[XLO] , 
                                                 sibcArray[surface].gbbox[YHI] - sibcArray[surface].gbbox[YLO] , 
                                                 sibcArray[surface].gbbox[ZHI] - sibcArray[surface].gbbox[ZLO] );
  }

  deallocArray( sibcArray , 1 , numSibcSurface );

  return;

}

/* SIBC E field update. */
void updateSibcSurfacesEfield( void )
{

  int i , j , k , p , q;
  int ii , jj , kk;
  real Eout[4];
  real Htan[4];
  real Hin[4];
  real A[4][4] = { { 0.0 } };
  real B[4][4] = { { 0.0 } };

  /* Find orientation and polarisation matrices for surface. 
   * This could be done at init phase and saved but could be a huge memory
     hog if there are many one face TBs in the mesh.
   */ 
  for( SurfaceIndex surface = 0 ; surface < numSibcSurface ; surface++ )
  { 

    /* Calculate pre and post filtering axes transformation matrices. */
    switch( sibcArray[surface].orientation )
    {
    case 1:
      matLinearComb( A , sibcArray[surface].cosa , Acp , sibcArray[surface].sina , Asp );
      matLinearComb( B , sibcArray[surface].cosa , Bcp , sibcArray[surface].sina , Bsp );
      break;
    case -1:
      matLinearComb( A , sibcArray[surface].cosa , Acm , sibcArray[surface].sina , Asm );
      matLinearComb( B , sibcArray[surface].cosa , Bcm , sibcArray[surface].sina , Bsm );
      break;
    default:
      assert( 0 );
      break;
    }

    switch( sibcArray[surface].normal )
    {
    case XDIR:
      for( i = sibcArray[surface].gbbox[XLO] , ii = 0 ; i < sibcArray[surface].gbbox[XHI] ; i++ , ii++ )
        for( j = sibcArray[surface].gbbox[YLO] , jj = 0 ; j < sibcArray[surface].gbbox[YHI] ; j++ , jj++ )
          for( k = sibcArray[surface].gbbox[ZLO] , kk = 0 ; k < sibcArray[surface].gbbox[ZHI] ; k++ , kk++ )
          { 
            /* Zero tangential electric field on mesh. These should be unnecessary. */
            Ey[i][j][k] = 0.0;
            Ey[i][j][k+1] = 0.0;
            Ez[i][j][k] = 0.0;
            Ez[i][j+1][k] = 0.0;
            /* Find magnetc field at face centre. */
            Htan[0] = 0.5 * ( 1 + sibcArray[surface].isAdjA[ii][jj][kk][0] + sibcArray[surface].isAdjA[ii][jj][kk][1] ) 
                          * ( UNSCALE_Hy( Hy[i-1][j][k] , j ) + UNSCALE_Hy( Hy[i-1][j+1][k] , j + 1 ) );
            Htan[1] = 0.5 * ( 1 + sibcArray[surface].isAdjB[ii][jj][kk][0] + sibcArray[surface].isAdjB[ii][jj][kk][1] ) 
                          * ( UNSCALE_Hy( Hy[i][j][k]   , j ) + UNSCALE_Hy( Hy[i][j+1][k]   , j + 1 ) );
            Htan[2] = 0.5 * ( 1 + sibcArray[surface].isAdjA[ii][jj][kk][2] + sibcArray[surface].isAdjA[ii][jj][kk][3] ) 
                          * ( UNSCALE_Hz( Hz[i-1][j][k] , k ) + UNSCALE_Hz( Hz[i-1][j][k+1] , k + 1 ) );
            Htan[3] = 0.5 * ( 1 + sibcArray[surface].isAdjB[ii][jj][kk][2] + sibcArray[surface].isAdjB[ii][jj][kk][3] )
                          * ( UNSCALE_Hz( Hz[i][j][k]   , k ) + UNSCALE_Hz( Hz[i][j][k+1]   , k + 1 ) );
            /* Transform magnetic field vectors from mesh to principal axes. */
            matMulVector( Hin , A , Htan );

            /* Apply SIBC. */
            for( p = 0 ; p < 4 ; p++ )
            {
              Eout[p] = 0.0;
              for( q = 0 ; q < 4 ; q++ )
                Eout[p] = Eout[p] + yfRecConvStep( sibcArray[surface].boundary->rcm.rc[p][q] , &sibcArray[surface].rcm_s[ii][jj][kk].rc_s[p][q] , Hin[q] );
            } 
            /* Transform output electric field from principal to mesh axes and store in surface. */ 
            matMulVector( sibcArray[surface].Etan[ii][jj][kk] , B , Eout );            
          }
      break;
    case YDIR:
      for( i = sibcArray[surface].gbbox[XLO] , ii = 0 ; i < sibcArray[surface].gbbox[XHI] ; i++ , ii++ )
        for( j = sibcArray[surface].gbbox[YLO] , jj = 0 ; j < sibcArray[surface].gbbox[YHI] ; j++ , jj++ )
          for( k = sibcArray[surface].gbbox[ZLO] , kk = 0 ; k < sibcArray[surface].gbbox[ZHI] ; k++ , kk++ )
          {
            /* Zero tangential electric field on mesh. These should be unnecessary. */
            Ez[i][j][k] = 0.0;
            Ez[i+1][j][k] = 0.0;
            Ex[i][j][k] = 0.0;
            Ex[i][j][k+1] = 0.0;
            /* Find magnetc field at face centre. */
            Htan[0] = 0.5 * ( 1 + sibcArray[surface].isAdjA[ii][jj][kk][0] + sibcArray[surface].isAdjA[ii][jj][kk][1] ) 
                          * ( UNSCALE_Hz( Hz[i][j-1][k] , k ) + UNSCALE_Hz( Hz[i][j-1][k+1] , k + 1 ) );
            Htan[1] = 0.5 * ( 1 + sibcArray[surface].isAdjB[ii][jj][kk][0] + sibcArray[surface].isAdjB[ii][jj][kk][1] )
                          * ( UNSCALE_Hz( Hz[i][j][k]   , k ) + UNSCALE_Hz( Hz[i][j][k+1]   , k + 1 ) );
            Htan[2] = 0.5 * ( 1 + sibcArray[surface].isAdjA[ii][jj][kk][2] + sibcArray[surface].isAdjA[ii][jj][kk][3] ) 
                          * ( UNSCALE_Hx( Hx[i][j-1][k] , i ) + UNSCALE_Hx( Hx[i+1][j-1][k] , i + 1 ) );
            Htan[3] = 0.5 * ( 1 + sibcArray[surface].isAdjB[ii][jj][kk][2] + sibcArray[surface].isAdjB[ii][jj][kk][3] )
                          * ( UNSCALE_Hx( Hx[i][j][k]   , i ) + UNSCALE_Hx( Hx[i+1][j][k]   , i + 1 ) );
            /* Transform magnetic field vectors from mesh to principal axes. */
            matMulVector( Hin , A , Htan );   
            /* Apply SIBC. */       
            for( p = 0 ; p < 4 ; p++ )
            {
              Eout[p] = 0.0;
              for( q = 0 ; q < 4 ; q++ )
                Eout[p] = Eout[p] + yfRecConvStep( sibcArray[surface].boundary->rcm.rc[p][q] , &sibcArray[surface].rcm_s[ii][jj][kk].rc_s[p][q] , Hin[q] );
            }
            /* Transform output electric field from principal to mesh axes and store in surface. */ 
            matMulVector( sibcArray[surface].Etan[ii][jj][kk] , B , Eout );
          }
      break;
    case ZDIR:
      for( i = sibcArray[surface].gbbox[XLO] , ii = 0 ; i < sibcArray[surface].gbbox[XHI] ; i++ , ii++ )
        for( j = sibcArray[surface].gbbox[YLO] , jj = 0 ; j < sibcArray[surface].gbbox[YHI] ; j++ , jj++ )
          for( k = sibcArray[surface].gbbox[ZLO] , kk = 0 ; k < sibcArray[surface].gbbox[ZHI] ; k++ , kk++ )
          {
            /* Zero tangential electric field on mesh. These should be unnecessary. */
            Ex[i][j][k] = 0.0;
            Ex[i][j+1][k] = 0.0;
            Ey[i][j][k] = 0.0;
            Ey[i+1][j][k] = 0.0;
            /* Find magnetc field at face centre. */
            Htan[0] = 0.5 * ( 1 + sibcArray[surface].isAdjA[ii][jj][kk][0] + sibcArray[surface].isAdjA[ii][jj][kk][1] )
                          * ( UNSCALE_Hx( Hx[i][j][k-1] , i ) + UNSCALE_Hx( Hx[i+1][j][k-1] , i + 1 ) );
            Htan[1] = 0.5 * ( 1 + sibcArray[surface].isAdjB[ii][jj][kk][0] + sibcArray[surface].isAdjB[ii][jj][kk][1] )
                          * ( UNSCALE_Hx( Hx[i][j][k]   , i ) + UNSCALE_Hx( Hx[i+1][j][k]   , i + 1 ) );
            Htan[2] = 0.5 * ( 1 + sibcArray[surface].isAdjA[ii][jj][kk][2] + sibcArray[surface].isAdjA[ii][jj][kk][3] )
                          * ( UNSCALE_Hy( Hy[i][j][k-1] , j ) + UNSCALE_Hy( Hy[i][j+1][k-1] , j + 1 ) );
            Htan[3] = 0.5 * ( 1 + sibcArray[surface].isAdjB[ii][jj][kk][2] + sibcArray[surface].isAdjB[ii][jj][kk][3] )
                          * ( UNSCALE_Hy( Hy[i][j][k]   , j ) + UNSCALE_Hy( Hy[i][j+1][k]   , j + 1 ) );
            /* Transform magnetic field vectors from mesh to principal axes. */
            matMulVector( Hin , A , Htan );
            /* Apply SIBC. */
            for( p = 0 ; p < 4 ; p++ )
            {
              Eout[p] = 0.0;
              for( q = 0 ; q < 4 ; q++ )
                Eout[p] = Eout[p] + yfRecConvStep( sibcArray[surface].boundary->rcm.rc[p][q] , &sibcArray[surface].rcm_s[ii][jj][kk].rc_s[p][q] , Hin[q] );
            }
            /* Transform output electric field from principal to mesh axes and store in surface. */ 
            matMulVector( sibcArray[surface].Etan[ii][jj][kk] , B , Eout );
          }         
      break;
    default:
      assert( 0 );
      break;
    } /* switch */

  } /* for */

  return;

}

/* Pre-multiply 4-vector with 4x4 matrix. */
void matMulVector( real y[4] , real A[4][4] , real x[4] )
{

  for( int p = 0 ; p < 4 ; p++ )
  {
    y[p] = 0.0;
    for( int q = 0 ; q < 4 ; q++ )
      y[p] = y[p] + A[p][q] * x[q]; 
  }

  return;

}

/* Linear combination of two 4x4 matrices. */
void matLinearComb( real C[4][4] , real a , real A[4][4] , real b , real B[4][4] )
{

  for( int p = 0 ; p < 4 ; p++ )
    for( int q = 0 ; q < 4 ; q++ )
      C[p][q] = a * A[p][q] + b * B[p][q];

  return;

}

/* SIBC H field correction. */
void updateSibcSurfacesHfield( void )
{

  int i , j , k;
  int ii , jj , kk;
  real Exa , Exb , Eya , Eyb , Eza , Ezb;
  real edgeWeightxl , edgeWeightyl , edgeWeightzl;
  real edgeWeightxh , edgeWeightyh , edgeWeightzh;
  
  for( SurfaceIndex surface = 0 ; surface < numSibcSurface ; surface++ )
  { 
    switch( sibcArray[surface].normal )
    {
    case XDIR:
      for( i = sibcArray[surface].gbbox[XLO] , ii = 0 ; i < sibcArray[surface].gbbox[XHI] ; i++ , ii++ )
        for( j = sibcArray[surface].gbbox[YLO] , jj = 0 ; j < sibcArray[surface].gbbox[YHI] ; j++ , jj++ )
          for( k = sibcArray[surface].gbbox[ZLO] , kk = 0 ; k < sibcArray[surface].gbbox[ZHI] ; k++ , kk++ )
          {
            /* Zero normal magnetic field on mesh. */
            Hx[i][j][k] = 0.0;
            /* Tangential electric fields. Half due to edge being shared by two faces. */   
            edgeWeightyl = 0.5 * ( 1 + isPmcEdge( ZDIR , k ) );
            edgeWeightyh = 0.5 * ( 1 + isPmcEdge( ZDIR , k + 1 ) );
            edgeWeightzl = 0.5 * ( 1 + isPmcEdge( YDIR , j ) );
            edgeWeightzh = 0.5 * ( 1 + isPmcEdge( YDIR , j + 1 ) );
            Eya = SCALE_Ey( sibcArray[surface].Etan[ii][jj][kk][0] , j );
            Eyb = SCALE_Ey( sibcArray[surface].Etan[ii][jj][kk][1] , j );
            Eza = SCALE_Ez( sibcArray[surface].Etan[ii][jj][kk][2] , k );
            Ezb = SCALE_Ez( sibcArray[surface].Etan[ii][jj][kk][3] , k );            
            /* Apply H field correction at each edge. */
            Hy[i-1][j][k]   = Hy[i-1][j][k]   + GAMMA_HY(i-1,j,k)   * edgeWeightzl * dEz_dx( Eza , i - 1 );
            Hy[i-1][j+1][k] = Hy[i-1][j+1][k] + GAMMA_HY(i-1,j+1,k) * edgeWeightzh * dEz_dx( Eza , i - 1 );
            Hy[i][j][k]     = Hy[i][j][k]     - GAMMA_HY(i,j,k)     * edgeWeightzl * dEz_dx( Ezb , i );
            Hy[i][j+1][k]   = Hy[i][j+1][k]   - GAMMA_HY(i,j+1,k)   * edgeWeightzh * dEz_dx( Ezb , i );
            Hz[i-1][j][k]   = Hz[i-1][j][k]   - GAMMA_HZ(i-1,j,k)   * edgeWeightyl * dEy_dx( Eya , i - 1 );
            Hz[i-1][j][k+1] = Hz[i-1][j][k+1] - GAMMA_HZ(i-1,j,k+1) * edgeWeightyh * dEy_dx( Eya , i - 1 );
            Hz[i][j][k]     = Hz[i][j][k]     + GAMMA_HZ(i,j,k)     * edgeWeightyl * dEy_dx( Eyb , i );
            Hz[i][j][k+1]   = Hz[i][j][k+1]   + GAMMA_HZ(i,j,k+1)   * edgeWeightyh * dEy_dx( Eyb , i );

            /* Apply parallel adjacent magnetic field correction */
            Hx[i][j+1][k] = Hx[i][j+1][k] - GAMMA_HX(i,j+1,k) * 0.5 * edgeWeightzh * (!sibcArray[surface].isAdjA[ii][jj][kk][0] * dEz_dx( Eza, i-1) + !sibcArray[surface].isAdjB[ii][jj][kk][0] * dEz_dx( Ezb, i));
            Hx[i][j-1][k] = Hx[i][j-1][k] + GAMMA_HX(i,j-1,k) * 0.5 * edgeWeightzl * (!sibcArray[surface].isAdjA[ii][jj][kk][1] * dEz_dx( Eza, i-1) + !sibcArray[surface].isAdjB[ii][jj][kk][1] * dEz_dx( Ezb, i));
            Hx[i][j][k+1] = Hx[i][j][k+1] - GAMMA_HX(i,j,k+1) * 0.5 * edgeWeightyh * (!sibcArray[surface].isAdjA[ii][jj][kk][2] * dEy_dx( Eya, i-1) + !sibcArray[surface].isAdjB[ii][jj][kk][2] * dEy_dx( Eyb, i));
            Hx[i][j][k-1] = Hx[i][j][k-1] + GAMMA_HX(i,j,k-1) * 0.5 * edgeWeightyl * (!sibcArray[surface].isAdjA[ii][jj][kk][3] * dEy_dx( Eya, i-1) + !sibcArray[surface].isAdjB[ii][jj][kk][3] * dEy_dx( Eyb, i));
          }
      break;
    case YDIR:
      for( i = sibcArray[surface].gbbox[XLO] , ii = 0 ; i < sibcArray[surface].gbbox[XHI] ; i++ , ii++ )
        for( j = sibcArray[surface].gbbox[YLO] , jj = 0 ; j < sibcArray[surface].gbbox[YHI] ; j++ , jj++ )
          for( k = sibcArray[surface].gbbox[ZLO] , kk = 0 ; k < sibcArray[surface].gbbox[ZHI] ; k++ , kk++ )
          {
            /* Zero normal magnetic field on mesh. */
            Hy[i][j][k] = 0.0;
            /* Tangential electric fields. Half due to edge being shared by two faces. */
            edgeWeightzl = 0.5 * ( 1 + isPmcEdge( XDIR , i ) );
            edgeWeightzh = 0.5 * ( 1 + isPmcEdge( XDIR , i + 1 ) ); 
            edgeWeightxl = 0.5 * ( 1 + isPmcEdge( ZDIR , k ) );
            edgeWeightxh = 0.5 * ( 1 + isPmcEdge( ZDIR , k + 1 ) );     
            Eza = SCALE_Ez( sibcArray[surface].Etan[ii][jj][kk][0] , k );
            Ezb = SCALE_Ez( sibcArray[surface].Etan[ii][jj][kk][1] , k );
            Exa = SCALE_Ex( sibcArray[surface].Etan[ii][jj][kk][2] , i );
            Exb = SCALE_Ex( sibcArray[surface].Etan[ii][jj][kk][3] , i );
            /* Apply H field correction at each edge. */
            Hz[i][j-1][k]   = Hz[i][j-1][k]   + GAMMA_HZ(i,j-1,k)   * edgeWeightxl * dEx_dy( Exa , j - 1 );
            Hz[i][j-1][k+1] = Hz[i][j-1][k+1] + GAMMA_HZ(i,j-1,k+1) * edgeWeightxh * dEx_dy( Exa , j - 1 );
            Hz[i][j][k]     = Hz[i][j][k]     - GAMMA_HZ(i,j,k)     * edgeWeightxl * dEx_dy( Exb , j );
            Hz[i][j][k+1]   = Hz[i][j][k+1]   - GAMMA_HZ(i,j,k+1)   * edgeWeightxh * dEx_dy( Exb , j );
            Hx[i][j-1][k]   = Hx[i][j-1][k]   - GAMMA_HX(i,j-1,k)   * edgeWeightzl * dEz_dy( Eza , j - 1 );
            Hx[i+1][j-1][k] = Hx[i+1][j-1][k] - GAMMA_HX(i+1,j-1,k) * edgeWeightzh * dEz_dy( Eza , j - 1 );
            Hx[i][j][k]     = Hx[i][j][k]     + GAMMA_HX(i,j,k)     * edgeWeightzl * dEz_dy( Ezb , j );
            Hx[i+1][j][k]   = Hx[i+1][j][k]   + GAMMA_HX(i+1,j,k)   * edgeWeightzh * dEz_dy( Ezb , j );

            /* Apply parallel adjacent magnetic field correction */
            Hy[i][j][k+1] = Hy[i][j][k+1] + GAMMA_HY(i,j,k+1) * 0.5 * edgeWeightzh * (!sibcArray[surface].isAdjA[ii][jj][kk][0] * dEz_dy( Eza, j-1) + !sibcArray[surface].isAdjB[ii][jj][kk][0] * dEz_dy( Exb, j));
            Hy[i][j][k-1] = Hy[i][j][k-1] - GAMMA_HY(i,j,k-1) * 0.5 * edgeWeightzl * (!sibcArray[surface].isAdjA[ii][jj][kk][1] * dEz_dy( Eza, j-1) + !sibcArray[surface].isAdjB[ii][jj][kk][1] * dEz_dy( Exb, j));
            Hy[i+1][j][k] = Hy[i+1][j][k] + GAMMA_HY(i+1,j,k) * 0.5 * edgeWeightxh * (!sibcArray[surface].isAdjA[ii][jj][kk][2] * dEx_dy( Exa, j-1) + !sibcArray[surface].isAdjB[ii][jj][kk][2] * dEx_dy( Exb, j));
            Hy[i-1][j][k] = Hy[i-1][j][k] - GAMMA_HY(i-1,j,k) * 0.5 * edgeWeightxl * (!sibcArray[surface].isAdjA[ii][jj][kk][3] * dEx_dy( Exa, j-1) + !sibcArray[surface].isAdjB[ii][jj][kk][3] * dEx_dy( Exb, j));
          }
      break;
    case ZDIR:
      for( i = sibcArray[surface].gbbox[XLO] , ii = 0 ; i < sibcArray[surface].gbbox[XHI] ; i++ , ii++ )
        for( j = sibcArray[surface].gbbox[YLO] , jj = 0 ; j < sibcArray[surface].gbbox[YHI] ; j++ , jj++ )
          for( k = sibcArray[surface].gbbox[ZLO] , kk = 0 ; k < sibcArray[surface].gbbox[ZHI] ; k++ , kk++ )
          {
            /* Zero normal magnetic field on mesh. */
            Hz[i][j][k] = 0.0;
            /* Tangential electric fields. Half due to edge being shared by two faces. */
            edgeWeightxl = 0.5 * ( 1 + isPmcEdge( YDIR , j ) );
            edgeWeightxh = 0.5 * ( 1 + isPmcEdge( YDIR , j + 1 ) );   
            edgeWeightyl = 0.5 * ( 1 + isPmcEdge( XDIR , i ) );
            edgeWeightyh = 0.5 * ( 1 + isPmcEdge( XDIR , i + 1 ) );  
            Exa = SCALE_Ex( sibcArray[surface].Etan[ii][jj][kk][0] , i );
            Exb = SCALE_Ex( sibcArray[surface].Etan[ii][jj][kk][1] , i );
            Eya = SCALE_Ey( sibcArray[surface].Etan[ii][jj][kk][2] , j );
            Eyb = SCALE_Ey( sibcArray[surface].Etan[ii][jj][kk][3] , j );
            /* Apply H field correction at each edge. */
            Hx[i][j][k-1]   = Hx[i][j][k-1]   + GAMMA_HX(i,j,k-1)   * edgeWeightyl * dEy_dz( Eya , k - 1 );
            Hx[i+1][j][k-1] = Hx[i+1][j][k-1] + GAMMA_HX(i+1,j,k-1) * edgeWeightyh * dEy_dz( Eya , k - 1 );
            Hx[i][j][k]     = Hx[i][j][k]     - GAMMA_HX(i,j,k)     * edgeWeightyl * dEy_dz( Eyb , k );
            Hx[i+1][j][k]   = Hx[i+1][j][k]   - GAMMA_HX(i+1,j,k)   * edgeWeightyh * dEy_dz( Eyb , k );
            Hy[i][j][k-1]   = Hy[i][j][k-1]   - GAMMA_HY(i,j,k-1)   * edgeWeightxl * dEx_dz( Exa , k - 1 );
            Hy[i][j+1][k-1] = Hy[i][j+1][k-1] - GAMMA_HY(i,j+1,k-1) * edgeWeightxh * dEx_dz( Exa , k - 1 );
            Hy[i][j][k]     = Hy[i][j][k]     + GAMMA_HY(i,j,k)     * edgeWeightxl * dEx_dz( Exb , k );
            Hy[i][j+1][k]   = Hy[i][j+1][k]   + GAMMA_HY(i,j+1,k)   * edgeWeightxh * dEx_dz( Exb , k );

            /* Apply parallel adjacent magnetic field correction */
            Hz[i+1][j][k] = Hz[i+1][j][k] + GAMMA_HZ(i+1,j,k) * 0.5 * edgeWeightyh * (!sibcArray[surface].isAdjA[ii][jj][kk][0] * dEy_dz( Eya, k-1) + !sibcArray[surface].isAdjB[ii][jj][kk][0] * dEy_dz( Eyb, k));
            Hz[i-1][j][k] = Hz[i-1][j][k] - GAMMA_HZ(i-1,j,k) * 0.5 * edgeWeightyl * (!sibcArray[surface].isAdjA[ii][jj][kk][1] * dEy_dz( Eya, k-1) + !sibcArray[surface].isAdjB[ii][jj][kk][1] * dEy_dz( Eyb, k));
            Hz[i][j+1][k] = Hz[i][j+1][k] + GAMMA_HZ(i,j+1,k) * 0.5 * edgeWeightxh * (!sibcArray[surface].isAdjA[ii][jj][kk][2] * dEx_dz( Exa, k-1) + !sibcArray[surface].isAdjB[ii][jj][kk][2] * dEx_dz( Exb, k));
            Hz[i][j-1][k] = Hz[i][j-1][k] - GAMMA_HZ(i,j-1,k) * 0.5 * edgeWeightxl * (!sibcArray[surface].isAdjA[ii][jj][kk][3] * dEx_dz( Exa, k-1) + !sibcArray[surface].isAdjB[ii][jj][kk][3] * dEx_dz( Exb, k));
          }
      break;
    default:
      assert( 0 );
      break;
    } /* switch */

  } /* for */

  
  for( SurfaceIndex surface = 0 ; surface < numSibcSurface ; surface++ )
  { 
    switch( sibcArray[surface].normal )
    {
    case XDIR:
      for( i = sibcArray[surface].gbbox[XLO] , ii = 0 ; i < sibcArray[surface].gbbox[XHI] ; i++ , ii++ )
        for( j = sibcArray[surface].gbbox[YLO] , jj = 0 ; j < sibcArray[surface].gbbox[YHI] ; j++ , jj++ )
          for( k = sibcArray[surface].gbbox[ZLO] , kk = 0 ; k < sibcArray[surface].gbbox[ZHI] ; k++ , kk++ )
          {
            /* Zero normal magnetic field on mesh. */
            Hx[i][j][k] = 0.0;
          }
      break;
    case YDIR:
      for( i = sibcArray[surface].gbbox[XLO] , ii = 0 ; i < sibcArray[surface].gbbox[XHI] ; i++ , ii++ )
        for( j = sibcArray[surface].gbbox[YLO] , jj = 0 ; j < sibcArray[surface].gbbox[YHI] ; j++ , jj++ )
          for( k = sibcArray[surface].gbbox[ZLO] , kk = 0 ; k < sibcArray[surface].gbbox[ZHI] ; k++ , kk++ )
          {
            /* Zero normal magnetic field on mesh. */
            Hy[i][j][k] = 0.0;
          }
      break;
    case ZDIR:
      for( i = sibcArray[surface].gbbox[XLO] , ii = 0 ; i < sibcArray[surface].gbbox[XHI] ; i++ , ii++ )
        for( j = sibcArray[surface].gbbox[YLO] , jj = 0 ; j < sibcArray[surface].gbbox[YHI] ; j++ , jj++ )
          for( k = sibcArray[surface].gbbox[ZLO] , kk = 0 ; k < sibcArray[surface].gbbox[ZHI] ; k++ , kk++ )
          {
            /* Zero normal magnetic field on mesh. */
            Hz[i][j][k] = 0.0;
           }
      break;
    default:
      assert( 0 );
      break;
    } /* switch */

  } /* for */
  
  return;

}

/* Return true if there are SIBC internal surfaces. */
bool thereAreSibcSurfaces( void )
{
  
  if( numSibcSurface > 0 )
    return true;
  else
    return false;

}

/* Convert 2x2 scattering matrix with free-space port impedances to impedance matrix. */
void tportS2Z( real Z[2][2] , real S[2][2] )
{

  real DeltaS;

  DeltaS = ( 1.0 - S[0][0] ) * ( 1.0 - S[1][1] ) - S[0][1] * S[1][0];
  Z[0][0] = ( ( 1.0 + S[0][0] ) * ( 1.0 - S[1][1] ) + S[0][1] * S[1][0] ) / DeltaS * eta0;
  Z[0][1] = 2.0 * S[0][1] * eta0 / DeltaS;
  Z[1][0] = 2.0 * S[1][0] * eta0 / DeltaS;
  Z[1][1] = ( ( 1.0 - S[0][0] ) * ( 1.0 + S[1][1] ) + S[0][1] * S[1][0] ) / DeltaS * eta0;
  
  return;

}

/* Verify scattering matrix is passive. */
bool isPassiveS( real S[2][2] )
{

  bool elem00;
  bool elem01;
  bool elem10;
  bool elem11;

  /* Compare elements of I - transpose(conj(S)) * S to zero matrix. */
  elem00 = ( 1 - ( S[0][0] * S[0][0] + S[1][0] * S[1][0] ) ) >= 0.0;
  elem01 = ( S[0][0] * S[0][1] - S[1][0] * S[1][1] ) >= 0.0;
  elem10 = ( S[0][1] * S[0][0] - S[1][1] * S[1][0] ) >= 0.0;
  elem11 = ( 1 - ( S[0][1] * S[0][1] + S[1][1] * S[1][1] ) ) >= 0.0;

  /* Passive if positive definite. */
  return ( elem00 && elem01 && elem10 && elem11 );

}

/* Set SIBC adjacency on array. */
void setSibcFace( bool ****isSibcFace , int gbbox[6] , CoordAxis dir , bool value )
{

  for( int i = gbbox[XLO]  ; i < gbbox[XHI] ; i++ )
    for( int j = gbbox[YLO]  ; j < gbbox[YHI] ; j++ )
      for( int k = gbbox[ZLO]  ; k < gbbox[ZHI] ; k++ )
        isSibcFace[i][j][k][dir] = value;

  return;

}