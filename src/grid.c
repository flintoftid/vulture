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
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>

#include "grid.h"
#include "boundary.h"
#include "surface.h"
#include "alloc_array.h"
#include "physical.h"
#include "message.h"
#include "bounding_box.h"
#include "simulation.h"
#include "medium.h"
#include "mesh.h"
#include "gnuplot.h"
#include "pml.h"
#include "memory.h"
#include "util.h"

/* Tolerance on grid type test */
#define GRID_TYPE_TOL 1e-5    

/* Tolerance on medium checking test. */
#define CHECK_LIMITS_RTOL    1e-8

/* 
 * Global variables. 
 */

/* Number of cells in each direction, including ghost cells. */
int numCells[3];

/* Mesh bounding box. */
int mbox[6];

/* Inner grid bounding box - holds the mesh. */
int gibox[6];

/* Outer grid boundaring box - include the PML cells. */
int gobox[6];

/* Ghost bounding box - include the PML and ghost cells. */
int ggbox[6];                 

/* Field array limits for inner grid. */
int gfilim[6][6];

/* Field array limits for outer grid. */
int gfolim[6][6];

/* Field array limits for ghost grid. */
int gfglim[6][6];

/* EM field arrays. */
real ***Ex;                   
real ***Ey;
real ***Ez;
real ***Hx;
real ***Hy;
real ***Hz;

/* Update coefficient arrays. */
#ifdef USE_INDEXED_MEDIA
  MediumIndex ***mediumEx;
  MediumIndex ***mediumEy;
  MediumIndex ***mediumEz;
  MediumIndex ***mediumHx;
  MediumIndex ***mediumHy;
  MediumIndex ***mediumHz;
#else
  real ***alphaEx;             
  real ***alphaEy;
  real ***alphaEz;
  real ***betaEx;
  real ***betaEy;
  real ***betaEz;
  real ***gammaHx;
  real ***gammaHy;
  real ***gammaHz;
#endif

/* Primary grid edge lengths. */
real *dex;                    
real *dey;
real *dez;

/* Secondary grid edge lengths. */
real *dhx;                    
real *dhy;
real *dhz;

/* Inverse primary grid edge lengths. */
real *idex;
real *idey;
real *idez;

/* Inverse secondary grid edge lengths. */
real *idhx;
real *idhy;
real *idhz;

/* Time step interval. */
real dt;

/* 
 * Private data.
 */

/* Grid type strings. */
char GRID_TYPE[4][11]  = { "CUBIC" , "UNIFORM" , "NONUNIFORM" , "UNDEFINED" };

/* Type of grid. */
static GridType gridType;

/* Number of mesh lines in each direction. */
static int numLines[3];

/* Line coordinates. */
static real *xlines = NULL;
static real *ylines = NULL;
static real *zlines = NULL;

/* Minimum and maximum edge lengths in each direction. */
static real dmin[3];
static real dmax[3];

/* Edge lengths for uniform grid type. */
static real duni[3];

/* 
 * Private method interfaces. 
 */

void setGridExtents( void );
void initCellEdges( void );
void initFieldArrayLimits( void );
void allocGridArrays( void );
void clearGrid( void );
void writeLines( char *fileName , int n , real *v );
void setUniformMesh( real v[] , int numLines , real del );
char decodeAlphaBeta( real alpha , real beta );
char decodeGamma( real gamma );
void setGridType( void );
real numPhaseVelocityFunc( real k , real A[3] , real B );

/*
 * Method Implementations.
 */

/* Parse mesh extents. */
bool parseDM( char *line )
{

  if( sscanf( line , "%d %d %d" , &numCells[XDIR] , &numCells[YDIR] , &numCells[ZDIR] ) != 3 )
    return false;  

  for( int direction = XDIR ; direction <= ZDIR ; direction++ )
  {
    if( numCells[direction] < 1 )
    {
      message( MSG_LOG , 0 , "  Number of mesh lines in %s direction must be greater than 1\n" , AXIS[direction] );
      return false;
    }
    else
    {
      numLines[direction] = numCells[direction] + 1;
    }
  }

  mbox[XLO] = 0;
  mbox[XHI] = numCells[XDIR];
  mbox[YLO] = 0;
  mbox[YHI] = numCells[YDIR];
  mbox[ZLO] = 0;
  mbox[ZHI] = numCells[ZDIR];

  return true;

}

/* Parse uniform mesh card. */
bool parseMS( char *line )
{

  int numScanned = 0;
  double delta[3];
  unsigned long bytes;

  numScanned = sscanf( line , "%lf %lf %lf" , &delta[XDIR] , &delta[YDIR] , &delta[ZDIR] );

  switch( numScanned )
  {
  case 1:
    delta[YDIR] = delta[XDIR];
    delta[ZDIR] = delta[XDIR];
    break;
  case 3:
    break;
  default:
    message( MSG_LOG , 0 , "  Invalid MS card - must have one or three parameters\n" );  
    return false;
    break;
  }

  /* Validate parameters. */
  for( int direction = XDIR ; direction <= ZDIR ; direction++ )
  {
    if( delta[direction] <= 0.0 )
    {
      message( MSG_LOG , 0 , "  Mesh increment in %s direction must be greater than 0\n" , AXIS[direction] );
      return false;
    }
  }

  /* Allocate and set mesh lines. */
  xlines = allocArray( &bytes , sizeof( real ) , 1 , numLines[XDIR] );
  memory.grid += bytes;
  ylines = allocArray( &bytes , sizeof( real ) , 1 , numLines[YDIR] );
  memory.grid += bytes;
  zlines = allocArray( &bytes , sizeof( real ) , 1 , numLines[ZDIR] );
  memory.grid += bytes;
  setUniformMesh( xlines , numLines[XDIR] , (real) (delta[XDIR]) );
  setUniformMesh( ylines , numLines[YDIR] , (real) (delta[YDIR]) );
  setUniformMesh( zlines , numLines[ZDIR] , (real) (delta[ZDIR]) );

  return true;

}

/* Parse non-uniform mesh card. */
bool parseXL( char *line )
{

  unsigned long bytes;

  xlines = allocArray( &bytes , sizeof( real ) , 1 , numLines[XDIR] );
  memory.grid += bytes;

  if( !meshReadRealArray( numLines[XDIR] , xlines ) )
  {
    message( MSG_LOG , 0 , "  Could not read x-lines from mesh file\n" );  
    return false; 
  }

  return true;

}

/* Parse non-uniform mesh card. */
bool parseYL( char *line )
{

  unsigned long bytes;

  ylines = allocArray( &bytes , sizeof( real ) , 1 , numLines[YDIR] );
  memory.grid += bytes;

  if( !meshReadRealArray( numLines[YDIR] , ylines ) )
  {
    message( MSG_LOG , 0 , "  Could not read y-lines from mesh file\n" );  
    return false; 
  }

  return true;

}

/* Parse non-uniform mesh card. */
bool parseZL( char *line )
{

  unsigned long bytes;

  zlines = allocArray( &bytes , sizeof( real ) , 1 , numLines[ZDIR] );
  memory.grid += bytes;

  if( !meshReadRealArray( numLines[ZDIR] , zlines ) )
  {
    message( MSG_LOG , 0 , "  Could not read z-lines from mesh file\n" );  
    return false; 
  }

  return true;

}

/* Create a uniform mesh in one dimension. */
void setUniformMesh( real v[] , int numLines , real del )
{

  for( int i = 0 ; i < numLines ; i++ )
    v[i] = i * del;

  return;

}

/* Initialise the grid fields. */
/* Depnds: initSimulation(). */
void initGrid( void )
{

  message( MSG_LOG , 0  , "\nInitialising the grid...\n\n" );

  /* Need to get the external boundary and surface information to set the grid extents. */
  initExternalSurfaceParameters();

  /* Set the grid extents. */
  setGridExtents();
 
  /* Find field array limits. */
  initFieldArrayLimits();

  /* Allocate grid arrays. */
  allocGridArrays();

  /* Determine cell edge lengths and time step. */
  initCellEdges();

  /* Determine the grid type. */
  setGridType();
  
  /* Clear the E and H field arrays. */
  clearGrid();

  return;

}

/* Set the extents of the inner and outer grid. */
void setGridExtents( void )
{

  /* Grid bounding boxes. */
  ggbox[XLO] = 0;
  ggbox[YLO] = 0;
  ggbox[ZLO] = 0;
  gobox[XLO] = NUM_GHOST_CELLS;
  gobox[YLO] = NUM_GHOST_CELLS;
  gobox[ZLO] = NUM_GHOST_CELLS;
  gibox[XLO] = gobox[XLO] + outerSurfaceNumLayers( XLO );
  gibox[YLO] = gobox[YLO] + outerSurfaceNumLayers( YLO );
  gibox[ZLO] = gobox[ZLO] + outerSurfaceNumLayers( ZLO );
  gibox[XHI] = gibox[XLO] + numLines[XDIR] - 1;
  gibox[YHI] = gibox[YLO] + numLines[YDIR] - 1;
  gibox[ZHI] = gibox[ZLO] + numLines[ZDIR] - 1;
  gobox[XHI] = gibox[XHI] + outerSurfaceNumLayers( XHI );
  gobox[YHI] = gibox[YHI] + outerSurfaceNumLayers( YHI );
  gobox[ZHI] = gibox[ZHI] + outerSurfaceNumLayers( ZHI );
  ggbox[XHI] = gobox[XHI] + 1;
  ggbox[YHI] = gobox[YHI] + 1;
  ggbox[ZHI] = gobox[ZHI] + 1;

  /* Cell allocation extents. */
  numCells[XDIR] = gobox[XHI] - gobox[XLO] + 2 * NUM_GHOST_CELLS;
  numCells[YDIR] = gobox[YHI] - gobox[YLO] + 2 * NUM_GHOST_CELLS;
  numCells[ZDIR] = gobox[ZHI] - gobox[ZLO] + 2 * NUM_GHOST_CELLS;

  return;

}

/* Initialise the field array limits for the grid. */
void initFieldArrayLimits( void )
{
                                  /* XLO     XHI     YLO     YHI     ZLO     ZHI */
  bool includeInnerBoundary[6] = {  true ,  true ,  true ,  true ,  true ,  true };
  bool includeOuterBoundary[6] = {  true ,  true ,  true ,  true ,  true ,  true };
  bool includeGhostBoundary[6] = { false , false , false , false , false , false };

  for( int boundary = XLO ; boundary <= ZHI ; boundary++ )
    if( outerSurfaceType( boundary ) == BT_MUR ) includeInnerBoundary[boundary] = false;
 
  message( MSG_LOG , 0 , "  Initialising array limits ...\n" );

  /* Field limits for inner, outer and ghost grid boundaries. */
  setFieldLimits( gibox , gfilim , includeInnerBoundary );
  setFieldLimits( gobox , gfolim , includeOuterBoundary );
  setFieldLimits( ggbox , gfglim , includeGhostBoundary );

  return;

}

/* Set the field array limits for a volume defined by cell limits, */
/* with option to include or exclude the tangnetial electric and */
/* normal magnetic fields on each face of the volume. */
void setFieldLimits( int bbox[6] , int fieldLimits[6][6] , bool includeBoundary[6] )
{

  for( FieldComponent field = EX ; field <= HZ ; field++ )
  {
    for( MeshFace face = XLO ; face <= ZHI ; face++ )
    {
      switch( face )
      {
      case XLO:
      case YLO: 
      case ZLO:
        if( fieldIsInBoundary( field , face ) )
          if( includeBoundary[face] )
            fieldLimits[field][face] = bbox[face];
          else
            fieldLimits[field][face] = bbox[face] + 1;
        else
          fieldLimits[field][face] = bbox[face];
        break;
      case XHI:
      case YHI:
      case ZHI:
        if( fieldIsInBoundary( field , face ) )
          if( includeBoundary[face] )
            fieldLimits[field][face] = bbox[face];
          else
            fieldLimits[field][face] = bbox[face] - 1;
        else
          fieldLimits[field][face] = bbox[face] - 1;   
        break;
      default:
        assert( 0 );
        break;
      }
    }
  }

  return;

}

/* Allocate the grid arrays. */
void allocGridArrays( void )
{

  unsigned long bytes;

  message( MSG_LOG , 0 , "  Allocating grid arrays...\n" );

  /* Allocate primary and secondary grid edge length arrays. */
  message( MSG_DEBUG1 , 0 , "  Allocating grid dex array\n" );
  dex = allocArray( &bytes , sizeof( real ) , 1 , numCells[XDIR] );
  memory.grid += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid dey array\n" );
  dey = allocArray( &bytes , sizeof( real ) , 1 , numCells[YDIR] );
  memory.grid += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid dez array\n" );
  dez = allocArray( &bytes , sizeof( real ) , 1 , numCells[ZDIR] );
  memory.grid += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid idex array\n" );
  idex = allocArray( &bytes , sizeof( real ) , 1 , numCells[XDIR] );
  memory.grid += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid idey array\n" );
  idey = allocArray( &bytes , sizeof( real ) , 1 , numCells[YDIR] );
  memory.grid += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid idez array\n" );
  idez = allocArray( &bytes , sizeof( real ) , 1 , numCells[ZDIR]);
  memory.grid += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid dhx array\n" );
  dhx = allocArray( &bytes , sizeof( real ) , 1 , numCells[XDIR] );
  memory.grid += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid dhy array\n" );
  dhy = allocArray( &bytes , sizeof( real ) , 1 , numCells[YDIR] );
  memory.grid += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid dhz array\n" );
  dhz = allocArray( &bytes , sizeof( real ) , 1 , numCells[ZDIR] );
  memory.grid += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid idhx array\n" );
  idhx = allocArray( &bytes , sizeof( real ) , 1 , numCells[XDIR] );
  memory.grid += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid idhy array\n" );
  idhy = allocArray( &bytes , sizeof( real ) , 1 , numCells[YDIR] );
  memory.grid += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid idhz array\n" );
  idhz = allocArray( &bytes , sizeof( real ) , 1 , numCells[ZDIR] );
  memory.grid += bytes;

  /* Allocate main#include "fdtd_types.h" field arrays. */
  message( MSG_DEBUG1 , 0 , "  Allocating grid Ex array\n" );
  Ex = allocArray( &bytes , sizeof( real ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
  memory.ehFields += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid Ey array\n" );
  Ey = allocArray( &bytes , sizeof( real ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
  memory.ehFields += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid Ez array\n" );
  Ez = allocArray( &bytes , sizeof( real ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
  memory.ehFields += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid Hx array\n" );
  Hx = allocArray( &bytes , sizeof( real ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
  memory.ehFields += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid Hy array\n" );
  Hy = allocArray( &bytes , sizeof( real ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
  memory.ehFields += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid Hz array\n" );
  Hz = allocArray( &bytes , sizeof( real ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
  memory.ehFields += bytes;

  /* Allocate media arrays. */
  #ifdef USE_INDEXED_MEDIA

    message( MSG_DEBUG1 , 0 , "  Allocating grid mediumEx array\n" );
    mediumEx = allocArray( &bytes , sizeof( MediumIndex ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    memory.ehCoeffs += bytes;
    message( MSG_DEBUG1 , 0 , "  Allocating grid mediumEy array\n" );
    mediumEy = allocArray( &bytes , sizeof( MediumIndex ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    memory.ehCoeffs += bytes;
    message( MSG_DEBUG1 , 0 , "  Allocating grid mediumEz array\n" );
    mediumEz = allocArray( &bytes , sizeof( MediumIndex ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    memory.ehCoeffs += bytes;
    message( MSG_DEBUG1 , 0 , "  Allocating grid mediumHx array\n" );
    mediumHx = allocArray( &bytes , sizeof( MediumIndex ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    memory.ehCoeffs += bytes;
    message( MSG_DEBUG1 , 0 , "  Allocating grid mediumHy array\n" );
    mediumHy = allocArray( &bytes , sizeof( MediumIndex ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    memory.ehCoeffs += bytes;
    message( MSG_DEBUG1 , 0 , "  Allocating grid mediumHz array\n" );
    mediumHz = allocArray( &bytes , sizeof( MediumIndex ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    memory.ehCoeffs += bytes;

  #else
  
    message( MSG_DEBUG1 , 0 , "  Allocating grid alphaEx array\n" );
    alphaEx = allocArray( &bytes , sizeof( real ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    memory.ehCoeffs += bytes;
    message( MSG_DEBUG1 , 0 , "  Allocating grid alphaEy array\n" );
    alphaEy = allocArray( &bytes , sizeof( real ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    memory.ehCoeffs += bytes;
    message( MSG_DEBUG1 , 0 , "  Allocating grid alphaEz array\n" );
    alphaEz = allocArray( &bytes , sizeof( real ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    memory.ehCoeffs += bytes;
    message( MSG_DEBUG1 , 0 , "  Allocating grid betaEx array\n" );
    betaEx = allocArray( &bytes , sizeof( real ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    memory.ehCoeffs += bytes;
    message( MSG_DEBUG1 , 0 , "  Allocating grid betaEy array\n" );
    betaEy = allocArray( &bytes , sizeof( real ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    memory.ehCoeffs += bytes;
    message( MSG_DEBUG1 , 0 , "  Allocating grid betaEz array\n" );
    betaEz = allocArray( &bytes , sizeof( real ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    memory.ehCoeffs += bytes;
    message( MSG_DEBUG1 , 0 , "  Allocating grid gammaHx array\n" );
    gammaHx = allocArray( &bytes , sizeof( real ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    memory.ehCoeffs += bytes;
    message( MSG_DEBUG1 , 0 , "  Allocating grid gammaHy array\n" );
    gammaHy = allocArray( &bytes , sizeof( real ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    memory.ehCoeffs += bytes;
    message( MSG_DEBUG1 , 0 , "  Allocating grid gammaHz array\n" );
    gammaHz = allocArray( &bytes , sizeof( real ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    memory.ehCoeffs += bytes;

  #endif
  
  message( MSG_DEBUG1 , 0 , "\n" );

  return;

}

/* Initialise the cell edge length arrays and time-step. */
void initCellEdges( void )
{

  /* Counters. */
  int i , j , k;

  message( MSG_LOG , 0 , "  Calculating cell edge lengths and time step...\n" );

  /* Primary grid edge lengths */
  dmin[XDIR] = REAL_MAX;
  dmin[YDIR] = REAL_MAX;
  dmin[ZDIR] = REAL_MAX;
  dmax[XDIR] = 0.0;
  dmax[YDIR] = 0.0;
  dmax[ZDIR] = 0.0;

  for( i = gobox[XLO] ; i <= gibox[XLO] - 1 ; i++ ) dex[i] = xlines[1] - xlines[0];
  for( i = gibox[XLO] ; i <= gibox[XHI] - 1 ; i++ ) dex[i] = xlines[i-gibox[XLO]+1] - xlines[i-gibox[XLO]];
  for( i = gibox[XHI] ; i <= gobox[XHI] - 1 ; i++ ) dex[i] = xlines[gibox[XHI]-gibox[XLO]] - xlines[gibox[XHI]-gibox[XLO]-1];
  dex[gobox[XLO]-1] = dex[gobox[XLO]];
  dex[gobox[XHI]] = dex[gobox[XHI]-1];

  for( j = gobox[YLO] ; j <= gibox[YLO] - 1 ; j++ ) dey[j] = ylines[1] - ylines[0];
  for( j = gibox[YLO] ; j <= gibox[YHI] - 1 ; j++ ) dey[j] = ylines[j-gibox[YLO]+1] - ylines[j-gibox[YLO]];
  for( j = gibox[YHI] ; j <= gobox[YHI] - 1 ; j++ ) dey[j] = ylines[gibox[YHI]-gibox[YLO]] - ylines[gibox[YHI]-gibox[YLO]-1];
  dey[gobox[YLO]-1] = dey[gobox[YLO]];
  dey[gobox[YHI]] = dey[gobox[YHI]-1];

  for( k = gobox[ZLO] ; k <= gibox[ZLO] - 1 ; k++ ) dez[k] = zlines[1] - zlines[0];
  for( k = gibox[ZLO] ; k <= gibox[ZHI] - 1 ; k++ ) dez[k] = zlines[k-gibox[ZLO]+1] - zlines[k-gibox[ZLO]];
  for( k = gibox[ZHI] ; k <= gobox[ZHI] - 1 ; k++ ) dez[k] = zlines[gibox[ZHI]-gibox[ZLO]] - zlines[gibox[ZHI]-gibox[ZLO]-1];
  dez[gobox[ZLO]-1] = dez[gobox[ZLO]];
  dez[gobox[ZHI]] = dez[gobox[ZHI]-1];

  /* Find minimum and maximum primary grid edge lengths. */
  message( MSG_DEBUG3 , 0 , "  Edge lengths along x-axis:\n" );
  for( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ ) {
    idex[i] = 1.0 / dex[i];
    if( dex[i] < dmin[XDIR] ) dmin[XDIR] = dex[i];
    if( dex[i] > dmax[XDIR] ) dmax[XDIR] = dex[i];
    message( MSG_DEBUG3 , 0 , "  dex[%d] = %e\n" , i , dex[i] );
  }
  message( MSG_DEBUG3 , 0 , "  Edge lengths along y-axis:\n" );
  for( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ ) {
    idey[j] = 1.0 / dey[j];
    if( dey[j] < dmin[YDIR] ) dmin[YDIR] = dey[j];
    if( dey[j] > dmax[YDIR] ) dmax[YDIR] = dey[j];
    message( MSG_DEBUG3 , 0 , "  dey[%d] = %e\n" , j , dey[j] );
  }
  message( MSG_DEBUG3 , 0 , "  Edge lengths along z-axis:\n" );
  for( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ ) {
    idez[k] = 1.0 / dez[k];
    if( dez[k] < dmin[ZDIR] ) dmin[ZDIR] = dez[k];
    if( dez[k] > dmax[ZDIR] ) dmax[ZDIR] = dez[k];
    message( MSG_DEBUG3 , 0 , "  dez[%d] = %e\n" , k , dez[k] );
  }
  message( MSG_DEBUG3 , 0 , "\n" );

  /* Secondary grid edge lengths */
  for( i = gobox[XLO] ; i <= gobox[XHI] ; i++ ) dhx[i] = 0.5 * ( dex[i] + dex[i-1] );
  for( j = gobox[YLO] ; j <= gobox[YHI] ; j++ ) dhy[j] = 0.5 * ( dey[j] + dey[j-1] );
  for( k = gobox[ZLO] ; k <= gobox[ZHI] ; k++ ) dhz[k] = 0.5 * ( dez[k] + dez[k-1] );
  dhx[gobox[XLO]-1] = dhx[gobox[XLO]];
  dhy[gobox[YLO]-1] = dhy[gobox[YLO]];
  dhz[gobox[ZLO]-1] = dhz[gobox[ZLO]];

  /* Find minimum secondary grid edge lengths. */
  for( i = gobox[XLO] ; i <= gobox[XHI] ; i++ ) {
    idhx[i] = 1.0 / dhx[i];
    if( dhx[i] < dmin[XDIR] ) dmin[XDIR] = dhx[i];
    if( dhx[i] > dmax[XDIR] ) dmax[XDIR] = dhx[i];
    message( MSG_DEBUG3 , 0 , "  dhx[%d] = %e\n" , i , dhx[i] );
  }
  for( j = gobox[YLO] ; j <= gobox[YHI] ; j++ ) {
    idhy[j] = 1.0 / dhy[j];
    if( dhy[j] < dmin[YDIR] ) dmin[YDIR] = dhy[j];
    if( dhy[j] > dmax[YDIR] ) dmax[YDIR] = dhy[j];
    message( MSG_DEBUG3 , 0 , "  dhy[%d] = %e\n" , j , dhy[j] );
  }
  for( k = gobox[ZLO] ; k <= gobox[ZHI] ; k++ ) {
    idhz[k] = 1.0 / dhz[k];
    if( dhz[k] < dmin[ZDIR] ) dmin[ZDIR] = dhz[k];
    if( dhz[k] > dmax[ZDIR] ) dmax[ZDIR] = dhz[k];
    message( MSG_DEBUG3 , 0 , "  dhz[%d] = %e\n" , k , dhz[k] );
  }
  idhx[gobox[XLO]-1] = idhx[gobox[XLO]];
  idhy[gobox[YLO]-1] = idhy[gobox[YLO]];
  idhz[gobox[ZLO]-1] = idhz[gobox[ZLO]];

  /* Time step. */
  dt = getCourantNumber() / c0 / sqrt( ( 1.0 / dmin[XDIR] ) * ( 1.0 / dmin[XDIR] )
                                     + ( 1.0 / dmin[YDIR] ) * ( 1.0 / dmin[YDIR] )
                                     + ( 1.0 / dmin[ZDIR] ) * ( 1.0 / dmin[ZDIR] ) );

  return;

}

/* Set initial field value througout the grid, including PML. */
void clearGrid( void )
{

  int i , j , k;

  message( MSG_LOG , 0 , "  Clearing the grid...\n" );

  /* Clear E/H over entire grid, including PML and mirror cells. */
  /* INITIAL_FIELD_VALUE is zero unless limit checking mode is enabled. */
  for ( i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ ) {
    for ( j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ ) {
      for ( k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ ) {
        Ex[i][j][k] = INITIAL_FIELD_VALUE;
        Ey[i][j][k] = INITIAL_FIELD_VALUE;
        Ez[i][j][k] = INITIAL_FIELD_VALUE;
        Hx[i][j][k] = INITIAL_FIELD_VALUE;
        Hy[i][j][k] = INITIAL_FIELD_VALUE;
        Hz[i][j][k] = INITIAL_FIELD_VALUE;
      }
    }
  }

  return;

}

/* Step electric fields in inner grid. */
void updateGridEfield( void )
{

  int i , j , k;
  real **Ex_i, **Ey_i, **Ez_i, **Hx_i, **Hy_i, **Hz_i;
  real *Ex_ij, *Ey_ij, *Ez_ij, *Hx_ij, *Hy_ij, *Hz_ij;
  real **Hz_i1 , **Hy_i1 , *Hz_ij1 , *Hz_i1j , *Hy_i1j , *Hx_ij1;

  /* Update Ex. */
  #ifdef WITH_OPENMP
    #pragma omp parallel for private( i , j , k , Ex_i , Ex_ij , Hy_i , Hy_ij , Hz_i , Hz_ij , Hz_ij1 )
  #endif
  for ( i = gfilim[EX][XLO] ; i <= gfilim[EX][XHI] ; i++ ) 
  {
    Ex_i = Ex[i];
    Hy_i = Hy[i];
    Hz_i = Hz[i];
    for ( j = gfilim[EX][YLO] ; j <= gfilim[EX][YHI] ; j++ ) 
    {
      Ex_ij = Ex_i[j];
      Hy_ij = Hy_i[j];
      Hz_ij = Hz_i[j];
      Hz_ij1 = Hz_i[j-1];
      for ( k = gfilim[EX][ZLO] ; k <= gfilim[EX][ZHI] ; k++ ) 
      {
	CHECK_NOT_VISITED( Ex_ij[k] );
        Ex_ij[k] = ALPHA_EX(i,j,k) * Ex_ij[k] + BETA_EX(i,j,k)
	  * curl_Hx( Hz_ij[k] , Hz_ij1[k] , Hy_ij[k-1] , Hy_ij[k] , i , j , k );
	MARK_AS_VISITED( Ex_ij[k] );  
      }
    }
  }

  /* Update Ey. */
  #ifdef WITH_OPENMP
    #pragma omp parallel for private( i , j , k , Ey_i , Ey_ij , Hx_i , Hx_ij , Hz_i , Hz_ij , Hz_i1 , Hz_i1j )
  #endif
  for ( i = gfilim[EY][XLO] ; i <= gfilim[EY][XHI] ; i++ ) 
  {
    Ey_i = Ey[i];
    Hx_i = Hx[i];
    Hz_i = Hz[i];
    Hz_i1 = Hz[i-1];   
    for ( j = gfilim[EY][YLO] ; j <= gfilim[EY][YHI] ; j++ ) 
    {
      Ey_ij = Ey_i[j];
      Hx_ij = Hx_i[j];
      Hz_ij = Hz_i[j];  
      Hz_i1j = Hz_i1[j];
      for ( k = gfilim[EY][ZLO] ; k <= gfilim[EY][ZHI] ; k++ ) 
      {
	CHECK_NOT_VISITED( Ey_ij[k] );
        Ey_ij[k] = ALPHA_EY(i,j,k) * Ey_ij[k] + BETA_EY(i,j,k)
	  * curl_Hy( Hx_ij[k] , Hx_ij[k-1] , Hz_i1j[k] , Hz_ij[k] , i , j , k ); 
	MARK_AS_VISITED( Ey_ij[k] );
      }
    }
  }

  /* Update Ez. */
  #ifdef WITH_OPENMP
    #pragma omp parallel for private( i , j , k , Ez_i , Ez_ij , Hx_i , Hx_ij , Hy_i , Hy_ij , Hy_i1 , Hy_i1j , Hx_ij1 )
  #endif
  for ( i = gfilim[EZ][XLO] ; i <= gfilim[EZ][XHI] ; i++ ) 
  {
    Ez_i = Ez[i];
    Hx_i = Hx[i];
    Hy_i = Hy[i];
    Hy_i1 = Hy[i-1];
    for ( j = gfilim[EZ][YLO] ; j <= gfilim[EZ][YHI] ; j++ ) 
    {
      Ez_ij = Ez_i[j];
      Hx_ij = Hx_i[j];
      Hy_ij = Hy_i[j];       
      Hy_i1j = Hy_i1[j];
      Hx_ij1 = Hx_i[j-1];
      for ( k = gfilim[EZ][ZLO] ; k <= gfilim[EZ][ZHI] ; k++ ) 
      {
	CHECK_NOT_VISITED( Ez_ij[k] );
        Ez_ij[k] = ALPHA_EZ(i,j,k) * Ez_ij[k] + BETA_EZ(i,j,k)
	  * curl_Hz( Hy_ij[k] , Hy_i1j[k] , Hx_ij1[k] , Hx_ij[k] , i , j , k );
	MARK_AS_VISITED( Ez_ij[k] );
      }
    }
  }

  return;

}

/* Step magnetic fields in inner grid. */
void updateGridHfield( void )
{

  int i , j , k;
  real **Ex_i, **Ey_i, **Ez_i, **Hx_i, **Hy_i, **Hz_i;
  real *Ex_ij, *Ey_ij, *Ez_ij, *Hx_ij, *Hy_ij, *Hz_ij;
  real **Ez_i1 , **Ey_i1 , *Ez_ij1 , *Ez_i1j , *Ex_ij1 , *Ey_i1j;

  /* Update Hx. */
  #ifdef WITH_OPENMP
    #pragma omp parallel for private( i , j , k , Ey_i , Ey_ij , Ez_i , Ez_ij , Hx_i , Hx_ij , Ez_ij1 )
  #endif
  for ( i = gfilim[HX][XLO] ; i <= gfilim[HX][XHI] ; i++ ) 
  {
    Hx_i = Hx[i];
    Ey_i = Ey[i];
    Ez_i = Ez[i];   
    for ( j = gfilim[HX][YLO] ; j <= gfilim[HX][YHI] ; j++ ) 
    {
      Hx_ij = Hx_i[j];
      Ey_ij = Ey_i[j];
      Ez_ij = Ez_i[j];   
      Ez_ij1 = Ez_i[j+1];
      for ( k = gfilim[HX][ZLO] ; k <= gfilim[HX][ZHI] ; k++ ) 
      {
        CHECK_NOT_VISITED( Hx_ij[k] );
        Hx_ij[k] = Hx_ij[k] + GAMMA_HX(i,j,k)
	  * curl_Ex( Ey_ij[k+1] , Ey_ij[k] , Ez_ij[k] , Ez_ij1[k] , i , j , k ); 
	MARK_AS_VISITED( Hx_ij[k] );
      }
    }
  }

  /* Update Hy. */
  #ifdef WITH_OPENMP
    #pragma omp parallel for private( i , j , k , Ex_i , Ex_ij , Ez_i , Ez_ij , Hy_i , Hy_ij , Ez_i1 , Ez_i1j )
  #endif
  for ( i = gfilim[HY][XLO] ; i <= gfilim[HY][XHI] ; i++ ) 
  {
    Hy_i = Hy[i];
    Ex_i = Ex[i];
    Ez_i = Ez[i];
    Ez_i1 = Ez[i+1];
    for ( j = gfilim[HY][YLO] ; j <= gfilim[HY][YHI] ; j++ ) 
    {
      Hy_ij = Hy_i[j];
      Ex_ij = Ex_i[j];
      Ez_ij = Ez_i[j];   
      Ez_i1j = Ez_i1[j];
      for ( k = gfilim[HY][ZLO] ; k <= gfilim[HY][ZHI] ; k++ ) 
      {
	CHECK_NOT_VISITED( Hy_ij[k] );
        Hy_ij[k] = Hy_ij[k] + GAMMA_HY(i,j,k)
	  * curl_Ey( Ez_i1j[k] , Ez_ij[k] , Ex_ij[k] , Ex_ij[k+1] , i , j , k );
	MARK_AS_VISITED( Hy_ij[k] );
      }
    }
  }

  /* Update Hz. */
  #ifdef WITH_OPENMP
    #pragma omp parallel for private( i , j , k , Ex_i , Ex_ij , Ey_i , Ey_ij , Hz_i , Hz_ij , Ex_ij1 , Ey_i1 , Ey_i1j )
  #endif
  for ( i = gfilim[HZ][XLO] ; i <= gfilim[HZ][XHI] ; i++ ) 
  {
    Hz_i = Hz[i];
    Ex_i = Ex[i];
    Ey_i = Ey[i];
    Ey_i1 = Ey[i+1];
    for ( j = gfilim[HZ][YLO] ; j <= gfilim[HZ][YHI] ; j++ ) 
    {
      Hz_ij = Hz_i[j];
      Ex_ij = Ex_i[j];
      Ey_ij = Ey_i[j];   
      Ex_ij1 = Ex_i[j+1];
      Ey_i1j = Ey_i1[j];
      for ( k = gfilim[HZ][ZLO] ; k <= gfilim[HZ][ZHI] ; k++ ) 
      {
        CHECK_NOT_VISITED( Hz_ij[k] );
        Hz_ij[k] = Hz_ij[k] + GAMMA_HZ(i,j,k)
	  * curl_Ez( Ex_ij1[k] , Ex_ij[k] , Ey_ij[k] , Ey_i1j[k] , i , j , k );
	MARK_AS_VISITED( Hz_ij[k] );
      }
    }
  }

  return;

}

/* Report grid. */
void reportGrid( void )
{

  message( MSG_LOG , 0 , "\nGrid characteristics:\n\n" );

  message( MSG_LOG , 0 , "  Grid is %s\n" , GRID_TYPE[gridType] );
    
  message( MSG_LOG , 0 , "  Number of lines x: %d y: %d z: %d\n" , numLines[XDIR] , numLines[YDIR] , numLines[ZDIR] );

  message( MSG_LOG , 0 , "  Mesh BBOX=[%d,%d,%d,%d,%d,%d]\n" ,
           mbox[XLO] , mbox[XHI] , mbox[YLO] , mbox[YHI] , mbox[ZLO] , mbox[ZHI] );

  message( MSG_LOG , 0 , "  Grid dimensions [cells]: %d x %d x %d\n" , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
  message( MSG_LOG , 0 , "  Grid size: %lu cells\n" , ( gobox[XHI] - gobox[XLO] ) * ( gobox[YHI] - gobox[YLO] ) * ( gobox[ZHI] - gobox[ZLO] ) );

  message( MSG_LOG , 0 , "  Inner grid: BBOX=[%d,%d,%d,%d,%d,%d]\n",
           gibox[XLO] , gibox[XHI] , gibox[YLO] ,
           gibox[YHI] , gibox[ZLO] , gibox[ZHI] );

  message( MSG_LOG , 0 , "  Outer grid: BBOX=[%d,%d,%d,%d,%d,%d]\n",
           gobox[XLO] , gobox[XHI] , gobox[YLO] ,
           gobox[YHI] , gobox[ZLO] , gobox[ZHI] );

    message( MSG_LOG , 0 , "  Ghost grid: BBOX=[%d,%d,%d,%d,%d,%d]\n",
           ggbox[XLO] , ggbox[XHI] , ggbox[YLO] ,
           ggbox[YHI] , ggbox[ZLO] , ggbox[ZHI] );
    
  message( MSG_LOG , 0 , "  Minimum edge lengths: DXMIN = %e, DYMIN = %e, DZMIN = %e\n" , dmin[XDIR] , dmin[YDIR] , dmin[ZDIR] );

  message( MSG_LOG , 0 , "  Maximum edge lengths: DXMAX = %e, DYMAX = %e, DZMAX = %e\n" , dmax[XDIR] , dmax[YDIR] , dmax[ZDIR] );

  message( MSG_LOG , 0 , "  Time step [s]: %e\n" , dt );
  message( MSG_LOG , 0 , "  CFLN [-]: %e\n" ,  getCourantNumber() );

  for( FieldComponent field = EX ; field <= HZ ; field++ )
  {
    message( MSG_DEBUG1 , 0 , "  Inner Grid %s field limits: [%d,%d,%d,%d,%d,%3d]\n", FIELD[field] ,
             gfilim[field][XLO] , gfilim[field][XHI] , gfilim[field][YLO] ,
             gfilim[field][YHI] , gfilim[field][ZLO] , gfilim[field][ZHI] );

    message( MSG_DEBUG1 , 0 , "  Outer Grid %s field limits: [%d,%d,%d,%d,%d,%3d]\n", FIELD[field] ,
             gfolim[field][XLO] , gfolim[field][XHI] , gfolim[field][YLO] , 
             gfolim[field][YHI] , gfolim[field][ZLO] , gfolim[field][ZHI] );
  }

  reportPml();

  /* Log total memory usage. */
  allocArrayReport();

  /* Write mesh lines to files. */
  writeLines( "xlines.dat" , numLines[XDIR] , xlines );
  writeLines( "ylines.dat" , numLines[YDIR] , ylines );
  writeLines( "zlines.dat" , numLines[ZDIR] , zlines );

  return;

}

/* Verify all fields are updated when in limit checking mode. */
void checkGrid( void )
{
  
  int i , j , k;
  int field;
  bool numError = 0;
  
  message( MSG_LOG , 0 , "  Checking all are fields updated...\n" ); 
    
  field = EX;
  for ( i = gfolim[field][XLO] ; i <= gfolim[field][XHI] ; i++ )
    for ( j = gfolim[field][YLO] ; j <= gfolim[field][YHI] ; j++ )
      for ( k = gfolim[field][ZLO] ; k <= gfolim[field][ZHI] ; k++ )
        if( Ex[i][j][k] != VISITED_FIELD_VALUE )
	{  
	  message( MSG_WARN , 0 , "*** Warning: %s[%d][%d][%d] = %e != %g\n" , FIELD[field] , i , j , k , Ex[i][j][k] , VISITED_FIELD_VALUE );
	  numError++;
	}
	
  field = EY;
  for ( i = gfolim[field][XLO] ; i <= gfolim[field][XHI] ; i++ )
    for ( j = gfolim[field][YLO] ; j <= gfolim[field][YHI] ; j++ )
      for ( k = gfolim[field][ZLO] ; k <= gfolim[field][ZHI] ; k++ )
        if( Ey[i][j][k] != VISITED_FIELD_VALUE )
	{
	  message( MSG_WARN , 0 , "*** Warning: %s[%d][%d][%d] = %e != %g\n" , FIELD[field] , i , j , k , Ey[i][j][k] , VISITED_FIELD_VALUE );
	  numError++;
	}
	
  field = EZ;
  for ( i = gfolim[field][XLO] ; i <= gfolim[field][XHI] ; i++ )
    for ( j = gfolim[field][YLO] ; j <= gfolim[field][YHI] ; j++ )
      for ( k = gfolim[field][ZLO] ; k <= gfolim[field][ZHI] ; k++ )
        if( Ez[i][j][k] != VISITED_FIELD_VALUE )
	{
	  message( MSG_WARN , 0 , "*** Warning: %s[%d][%d][%d] = %e != %g\n" , FIELD[field] , i , j , k , Ez[i][j][k] , VISITED_FIELD_VALUE );
	  numError++;
	}
      
  field = HX;
  for ( i = gfolim[field][XLO] ; i <= gfolim[field][XHI] ; i++ )
    for ( j = gfolim[field][YLO] ; j <= gfolim[field][YHI] ; j++ )
      for ( k = gfolim[field][ZLO] ; k <= gfolim[field][ZHI] ; k++ )
        if( Hx[i][j][k] != VISITED_FIELD_VALUE )
	{
	  message( MSG_WARN , 0 , "*** Warning: %s[%d][%d][%d] = %e != %g\n" , FIELD[field] , i , j , k , Hx[i][j][k] , VISITED_FIELD_VALUE );
	  numError++;
	}
	
  field = HY;
  for ( i = gfolim[field][XLO] ; i <= gfolim[field][XHI] ; i++ )
    for ( j = gfolim[field][YLO] ; j <= gfolim[field][YHI] ; j++ )
      for ( k = gfolim[field][ZLO] ; k <= gfolim[field][ZHI] ; k++ )
        if( Hy[i][j][k] != VISITED_FIELD_VALUE )
	{
	  message( MSG_WARN , 0 , "*** Warning: %s[%d][%d][%d] = %e != %g\n" , FIELD[field] , i , j , k , Hy[i][j][k] , VISITED_FIELD_VALUE );
	  numError++;
	}
	
  field = HZ;
  for ( i = gfolim[field][XLO] ; i <= gfolim[field][XHI] ; i++ )
    for ( j = gfolim[field][YLO] ; j <= gfolim[field][YHI] ; j++ )
      for ( k = gfolim[field][ZLO] ; k <= gfolim[field][ZHI] ; k++ )
        if( Hz[i][j][k] != VISITED_FIELD_VALUE )
	{
	  message( MSG_WARN , 0 , "*** Warning: %s[%d][%d][%d] = %e != %g\n" , FIELD[field] , i , j , k , Hz[i][j][k] , VISITED_FIELD_VALUE );
	  numError++;
	}
	
  if( numError > 0 )
    message( MSG_ERROR , 0 , "  ** Found %d field elements that have not been updated! ** \n" , numError );    

  return;
  
}

/* Deallocate the grid arrays. */
void deallocGridArrays( void )
{

  message( MSG_DEBUG1 , 0 , "Deallocating the grid...\n" );

  message( MSG_DEBUG1 , 0 , "  Deallocating grid dex array\n" );
  deallocArray( dex , 1 , numCells[XDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid dey array\n" );
  deallocArray( dey , 1 , numCells[YDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid dez array\n" );
  deallocArray( dez , 1 , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid dhx array\n" );
  deallocArray( dhx , 1 , numCells[XDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid dhy array\n" );
  deallocArray( dhy , 1 , numCells[YDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid dhz array\n" );
  deallocArray( dhz , 1 , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid idex array\n" );
  deallocArray( idex , 1 , numCells[XDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid idey array\n" );
  deallocArray( idey , 1 , numCells[YDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid idez array\n" );
  deallocArray( idez , 1 , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid idhx array\n" );
  deallocArray( idhx , 1 , numCells[XDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid idhy array\n" );
  deallocArray( idhy , 1 , numCells[YDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid idhz array\n" );
  deallocArray( idhz , 1 , numCells[ZDIR] );

  message( MSG_DEBUG1 , 0 , "  Deallocating grid Ex array\n" );
  deallocArray( Ex , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid Ey array\n" );
  deallocArray( Ey , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid Ez array\n" );
  deallocArray( Ez , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid Hx array\n" );
  deallocArray( Hx , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid Hy array\n" );
  deallocArray( Hy , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid Hz array\n" );
  deallocArray( Hz , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );

  #ifdef USE_INDEXED_MEDIA
    message( MSG_DEBUG1 , 0 , "  Deallocating grid mediumEx array\n" );
    deallocArray( mediumEx , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    message( MSG_DEBUG1 , 0 , "  Deallocating grid mediumEy array\n" );
    deallocArray( mediumEy , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    message( MSG_DEBUG1 , 0 , "  Deallocating grid mediumEz array\n" );
    deallocArray( mediumEz , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    message( MSG_DEBUG1 , 0 , "  Deallocating grid mediumHx array\n" );
    deallocArray( mediumHx , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    message( MSG_DEBUG1 , 0 , "  Deallocating grid mediumHy array\n" );
    deallocArray( mediumHy , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    message( MSG_DEBUG1 , 0 , "  Deallocating grid mediumHz array\n" );
    deallocArray( mediumHz , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );  
  #else
    message( MSG_DEBUG1 , 0 , "  Deallocating grid alphaEx array\n" );
    deallocArray( alphaEx , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    message( MSG_DEBUG1 , 0 , "  Deallocating grid alphaEy array\n" );
    deallocArray( alphaEy , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    message( MSG_DEBUG1 , 0 , "  Deallocating grid alphaEz array\n" );
    deallocArray( alphaEz , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    message( MSG_DEBUG1 , 0 , "  Deallocating grid betaEx array\n" );
    deallocArray( betaEx , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    message( MSG_DEBUG1 , 0 , "  Deallocating grid betaEy array\n" );
    deallocArray( betaEy , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    message( MSG_DEBUG1 , 0 , "  Deallocating grid betaEz array\n" );
    deallocArray( betaEz , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    message( MSG_DEBUG1 , 0 , "  Deallocating grid gammaHx array\n" );
    deallocArray( gammaHx , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    message( MSG_DEBUG1 , 0 , "  Deallocating grid gammaHy array\n" );
    deallocArray( gammaHy , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
    message( MSG_DEBUG1 , 0 , "  Deallocating grid gammaHz array\n" );
    deallocArray( gammaHz , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
  #endif
    
  message( MSG_DEBUG1 , 0 , "  Deallocating xlines array\n" );
  deallocArray( xlines , 1 , numLines[XDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating ylines array\n" );
  deallocArray( ylines , 1 , numLines[YDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating zlines array\n" );
  deallocArray( zlines , 1 , numLines[ZDIR] );

  return;

}

/* Get maximum edge length in requested direction. */
real getGridMaxEdgeLength( CoordAxis direction )
{
  
  return dmax[direction];

}

/* Get time step. */
real getGridTimeStep( void )
{

  return dt;

}

/* Get number of cells in mesh. */
void getGridNumCells( int numMeshCells[3] )
{

  numMeshCells[XDIR] = numCells[XDIR];
  numMeshCells[YDIR] = numCells[YDIR];
  numMeshCells[ZDIR] = numCells[ZDIR];

  return;

}

/* Initialise media arrays to free space. */
/* Depends: initMedia */
void initMediaArrays( void )
{

  FaceMask mask;
  
  mask = FACE_MASKS[XLO] | FACE_MASKS[YLO] | FACE_MASKS[ZLO];
  assert( mask == 0x2A );
  //bool includeBoundary[6] = {  true ,  false ,  true ,  false ,  true ,  false };

  message( MSG_LOG , 0 , "\nInitialising grid media arrays ...\n\n" );

  setMediumOnGrid( ggbox , MT_FREE_SPACE , mask );
  
  return;

}

/* Write mesh lines to file for processing tools. */
void writeLines( char *fileName , int n , real *v )
{

  FILE *fp = NULL;

  fp = fopen( fileName , "w" );
  if( NULL == fp ) 
    message( MSG_ERROR , 0 , "*** Error: Failed to open %s file.\n" , fileName );
  for( int i = 0 ; i < n ; i++ )
    fprintf( fp , "%d %e\n" , i , v[i] );
  fclose( fp );

  return;

}

/* Convert bounding box to physical units. */
void bboxInPhysicalUnits( real physbbox[6] , int bbox[6] )
{

  physbbox[XLO] = xlines[bbox[XLO]];
  physbbox[XHI] = xlines[bbox[XHI]];
  physbbox[YLO] = ylines[bbox[YLO]];
  physbbox[YHI] = ylines[bbox[YHI]];
  physbbox[ZLO] = zlines[bbox[ZLO]];
  physbbox[ZHI] = zlines[bbox[ZHI]];

  return;

}

/* Get physical coordinate of x mesh line. */
real getMeshLineCoordX( int lineNumber )
{
  return xlines[lineNumber];
}

/* Get physical coordinate of y mesh line. */
real getMeshLineCoordY( int lineNumber )
{
  return ylines[lineNumber];
}

/* Get physical coordinate of z mesh line. */
real getMeshLineCoordZ( int lineNumber )
{
  return zlines[lineNumber];
}

/* Get physical coordinates of node. */
void getMeshNodeCoords( real nodeCoords[3] , int nodeIndices[3] )
{

  nodeCoords[XDIR] = xlines[nodeIndices[XDIR]];
  nodeCoords[YDIR] = ylines[nodeIndices[YDIR]];
  nodeCoords[ZDIR] = zlines[nodeIndices[ZDIR]];
    
  return;
  
}

/* Convert coordinate line index to physical units. */
real indexInPhysicalUnits( int index , CoordAxis dir )
{

  real value = 0.0;
  
  switch( dir )
  {
  case XDIR:
    value = xlines[index];
    break;
  case YDIR:
    value = ylines[index];
    break;
  case ZDIR:
    value = zlines[index];
    break;
  default:
    assert( 0 );
    break;
  }

  return value;
  
}

/* Convert coordinate line index as real to physical units. */
void nodeInPhysicalUnits( real r[3] , real rijk[3] )
{

  int ijk[3];
  real frac[3];
  
  for( CoordAxis dir = XDIR ; dir <= ZDIR ; dir++ )
  {
    ijk[dir] = floor( rijk[dir] );
    frac[dir] = rijk[dir] - ijk[dir];
  }
  
  if( ijk[XDIR] < 0 )
    r[XDIR] = xlines[0] + rijk[XDIR] * ( xlines[1] - xlines[0] );
  else if( ijk[XDIR] > numLines[XDIR] - 2 )
    r[XDIR] = xlines[numLines[XDIR]-1] + ( rijk[XDIR] - numLines[XDIR] + 1 ) * ( xlines[numLines[XDIR]-1] - xlines[numLines[XDIR]-2] );
  else
     r[XDIR] = ( 1.0 - frac[XDIR] ) * xlines[ijk[XDIR]] + frac[XDIR] * xlines[ijk[XDIR]+1];

  if( ijk[YDIR] < 0 )
    r[YDIR] = ylines[0] + rijk[YDIR] * ( ylines[1] - ylines[0] );
  else if( ijk[YDIR] > numLines[YDIR] - 2 )
    r[YDIR] = ylines[numLines[YDIR]-1] + ( rijk[YDIR] - numLines[YDIR] + 1 ) * ( ylines[numLines[YDIR]-1] - ylines[numLines[YDIR]-2] );
  else
    r[YDIR] = ( 1.0 - frac[YDIR] ) * ylines[ijk[YDIR]] + frac[YDIR] * ylines[ijk[YDIR]+1];
  
  if( ijk[ZDIR] < 0 )
    r[ZDIR] = zlines[0] + rijk[ZDIR] * ( zlines[1] - zlines[0] );
  else if( ijk[ZDIR] > numLines[ZDIR] - 2 )
    r[ZDIR] = zlines[numLines[ZDIR]-1] + ( rijk[ZDIR] - numLines[ZDIR] + 1 ) * ( zlines[numLines[ZDIR]-1] - zlines[numLines[ZDIR]-2] );
  else
     r[ZDIR] = ( 1.0 - frac[ZDIR] ) * zlines[ijk[ZDIR]] + frac[ZDIR] * zlines[ijk[ZDIR]+1];
  
  //r[XDIR] = ( 1.0 - frac[XDIR] ) * xlines[ijk[XDIR]] + frac[XDIR] * xlines[ijk[XDIR]+1];
  //r[YDIR] = ( 1.0 - frac[YDIR] ) * ylines[ijk[YDIR]] + frac[YDIR] * ylines[ijk[YDIR]+1];
  //r[ZDIR] = ( 1.0 - frac[ZDIR] ) * zlines[ijk[ZDIR]] + frac[ZDIR] * zlines[ijk[ZDIR]+1];
   
  return;
  
}

/* Output gnuplot compatible data for mesh lines. */
void gnuplotGridLines( void )
{

  char linesFileName[] = "gnuplot-lines.dat";
  FILE *outputFile;
  int ibbox[6];

  outputFile = fopen( linesFileName , "w" );
  if( !outputFile )
    message( MSG_ERROR , 0 , "*** Error: Failed to open lines output file %s\n" , linesFileName );

  gnuplotProblemSize( outputFile , mbox );

  /* xlines */
  ibbox[YLO] = mbox[YLO];
  ibbox[YHI] = mbox[YHI];
  ibbox[ZLO] = mbox[ZLO];
  ibbox[ZHI] = mbox[ZLO];
  for( int i = 0 ; i < numLines[XDIR] ; i++ )
  {
    ibbox[XLO] = i;
    ibbox[XHI] = i;
    gnuplotBoundingBox( outputFile , ibbox );
  }

  ibbox[YLO] = mbox[YHI];
  ibbox[YHI] = mbox[YHI];
  ibbox[ZLO] = mbox[ZLO];
  ibbox[ZHI] = mbox[ZHI];
  for( int i = 0 ; i < numLines[XDIR] ; i++ )
  {
    ibbox[XLO] = i;
    ibbox[XHI] = i;
    gnuplotBoundingBox( outputFile , ibbox );
  }

  /* ylines */
  ibbox[XLO] = mbox[XLO];
  ibbox[XHI] = mbox[XHI];
  ibbox[ZLO] = mbox[ZLO];
  ibbox[ZHI] = mbox[ZLO];
  for( int j = 0 ; j < numLines[YDIR] ; j++ )
  {
    ibbox[YLO] = j;
    ibbox[YHI] = j;
    gnuplotBoundingBox( outputFile , ibbox );
  }

  ibbox[XLO] = mbox[XLO];
  ibbox[XHI] = mbox[XLO];
  ibbox[ZLO] = mbox[ZLO];
  ibbox[ZHI] = mbox[ZHI];
  for( int j = 0 ; j < numLines[YDIR] ; j++ )
  {
    ibbox[YLO] = j;
    ibbox[YHI] = j;
    gnuplotBoundingBox( outputFile , ibbox );
  }

  /* zlines */
  ibbox[XLO] = mbox[XLO];
  ibbox[XHI] = mbox[XHI];
  ibbox[YLO] = mbox[YHI];
  ibbox[YHI] = mbox[YHI];
  for( int k = 0 ; k < numLines[ZDIR] ; k++ )
  {
    ibbox[ZLO] = k;
    ibbox[ZHI] = k;
    gnuplotBoundingBox( outputFile , ibbox );
  }

  ibbox[XLO] = mbox[XLO];
  ibbox[XHI] = mbox[XLO];
  ibbox[YLO] = mbox[YLO];
  ibbox[YHI] = mbox[YHI];
  for( int k = 0 ; k < numLines[ZDIR] ; k++ )
  {
    ibbox[ZLO] = k;
    ibbox[ZHI] = k;
    gnuplotBoundingBox( outputFile , ibbox );
  }

  fclose( outputFile );

  return;

}

/* Get grid inner and outer bounding boxes. */
void getGridBoundingBox( int innerBox[6] , int outerBox[6] )
{

  for( MeshFace face = XLO ; face <= ZHI ; face++ )
  {
    innerBox[face] = gibox[face];
    outerBox[face] = gobox[face];
  }

  return;

}

/* Get location of field point (in grid indices.) in physical units. */
void getFieldPhysicalLocation( real r[3] , FieldComponent field , int ig , int jg , int kg )
{

  int im, jm, km;   // Mesh indices.

  im = ig - gibox[XLO];
  jm = jg - gibox[YLO];
  km = kg - gibox[ZLO];

  switch( field )
  {
  case EX:
    r[XDIR] = xlines[im] + 0.5 * dex[ig];
    r[YDIR] = ylines[jm];
    r[ZDIR] = zlines[km];
    break;
  case EY:
    r[XDIR] = xlines[im];
    r[YDIR] = ylines[jm] + 0.5 * dey[jg];
    r[ZDIR] = zlines[km];
    break;
  case EZ:
    r[XDIR] = xlines[im];
    r[YDIR] = ylines[jm];
    r[ZDIR] = zlines[km] + 0.5 * dez[kg];
    break;
  case HX:
    r[XDIR] = xlines[im];
    r[YDIR] = ylines[jm] + 0.5 * dey[jg];
    r[ZDIR] = zlines[km] + 0.5 * dez[kg];
    break;
  case HY:
    r[XDIR] = xlines[im] + 0.5 * dex[ig];
    r[YDIR] = ylines[jm];
    r[ZDIR] = zlines[km] + 0.5 * dez[kg];
    break;
  case HZ:
    r[XDIR] = xlines[im] + 0.5 * dex[ig];
    r[YDIR] = ylines[jm] + 0.5 * dey[jg];
    r[ZDIR] = zlines[km];
    break;
  default:
    assert( 0 );
    break;
  }

  return;

}

/* Get location of field point in grid indices. */
void getFieldIndexLocation( real ijk[3] , FieldComponent field , int ig , int jg , int kg )
{

  switch( field )
  {
  case EX:
    ijk[XDIR] = ig + 0.5;
    ijk[YDIR] = jg;
    ijk[ZDIR] = kg;
    break;
  case EY:
    ijk[XDIR] = ig;
    ijk[YDIR] = jg + 0.5;
    ijk[ZDIR] = kg;
    break;
  case EZ:
    ijk[XDIR] = ig;
    ijk[YDIR] = jg;
    ijk[ZDIR] = kg + 0.5;
    break;
  case HX:
    ijk[XDIR] = ig;
    ijk[YDIR] = jg + 0.5;
    ijk[ZDIR] = kg + 0.5;
    break;
  case HY:
    ijk[XDIR] = ig + 0.5;
    ijk[YDIR] = jg;
    ijk[ZDIR] = kg + 0.5;
    break;
  case HZ:
    ijk[XDIR] = ig + 0.5;
    ijk[YDIR] = jg + 0.5;
    ijk[ZDIR] = kg;
    break;
  default:
    assert( 0 );
    break;
  }

  return;

}

/* Get location of node (in grid indices) in physical units. */
void getNodeLocation( real r[3] , int ig , int jg , int kg )
{

  r[XDIR] = xlines[ig - gibox[XLO]];
  r[YDIR] = ylines[jg - gibox[YLO]];
  r[ZDIR] = zlines[kg - gibox[ZLO]];

  return;

}

/* Determine and set grid type. */
void setGridType( void )
{

  bool xIsUniform = true;
  bool yIsUniform = true;
  bool zIsUniform = true;

  for( int i = gibox[XLO] + 1 ; i <= gibox[XHI] - 1 ; i++ ) 
    if( dex[i] - dex[i-1] >  GRID_TYPE_TOL )
      xIsUniform = false;

  for( int j = gibox[YLO] + 1 ; j <= gibox[YHI] - 1 ; j++ ) 
    if( dey[j] - dey[j-1] >  GRID_TYPE_TOL )
      yIsUniform = false;

  for( int k = gibox[ZLO] + 1 ; k <= gibox[ZHI] - 1 ; k++ ) 
    if( dez[k] - dez[k-1] >  GRID_TYPE_TOL )
      zIsUniform = false;
    
  if( xIsUniform && yIsUniform && zIsUniform )
  {
    if( ( fabs( dex[gibox[XLO]] - dey[gibox[YLO]] ) < GRID_TYPE_TOL ) &&
        ( fabs( dey[gibox[YLO]] - dez[gibox[ZLO]] ) < GRID_TYPE_TOL ) )
    {
      gridType = GT_CUBIC;
    }
    else
    {
      gridType = GT_UNIFORM;
    }
    duni[XDIR] = dex[gibox[XLO]];
    duni[YDIR] = dey[gibox[YLO]];
    duni[ZDIR] = dez[gibox[ZLO]];
  }
  else
  {
    gridType = GT_NONUNIFORM;
    duni[XDIR] = -1.0;
    duni[YDIR] = -1.0;
    duni[ZDIR] = -1.0;
  }

  return;

}

/* Get grid type */
GridType getGridType( void )
{
 
  return gridType;
  
}

/* Get grid size for uniform grid */
void getUniformGridSize( real d[3] )
{
 
  d[XDIR] = duni[XDIR];
  d[YDIR] = duni[YDIR];
  d[ZDIR] = duni[ZDIR];

  return ;
  
}

/* Equation for numerical phase velocity determination. */
real numPhaseVelocityFunc( real k , real A[3] , real B )
{
  real tmp[3];
  real func;
  real derivative;
  
  tmp[XDIR] = sin( A[XDIR] * k ) / duni[XDIR];
  tmp[YDIR] = sin( A[YDIR] * k ) / duni[YDIR];
  tmp[ZDIR] = sin( A[ZDIR] * k ) / duni[ZDIR];
     
  func = tmp[XDIR] * tmp[XDIR] + tmp[YDIR] * tmp[YDIR] + tmp[ZDIR] * tmp[ZDIR] - B * B;
  
  derivative = A[XDIR] * sin( 2 * A[XDIR] * k ) / ( duni[XDIR] * duni[XDIR] ) + 
               A[YDIR] * sin( 2 * A[YDIR] * k ) / ( duni[YDIR] * duni[YDIR] ) +
               A[ZDIR] * sin( 2 * A[ZDIR] * k ) / ( duni[ZDIR] * duni[ZDIR] );
  
  return func / derivative;
  
}

/* Determine numerical phase velocity using Newton-Rhapson iteration. */
real numericalPhaseVelocity( real theta , real phi )
{
  real A[3];
  real B;
  real w;
  real k;

  /* Shouldn't be doing this on a nonuniform grid. */
  assert( gridType != GT_NONUNIFORM );
  
  /* Choose frequency based on mesh size. */
  w = 2.0 * pi * ( 1.0 / dt / 23 );

  A[XDIR] = 0.5 * duni[XDIR] * sin( theta ) * cos( phi );
  A[YDIR] = 0.5 * duni[YDIR] * sin( theta ) * sin( phi );
  A[ZDIR] = 0.5 * duni[ZDIR] * cos( theta );
  B = ( sin( 0.5 * w * dt ) / ( c0 * dt ) ); 

  k = w / c0;
  for( int i = 0 ; i < 10 ; i++ )
    k = k - numPhaseVelocityFunc( k , A , B );

  return w / k; 
  
}

/* Print ASCII dump of the material arrays. */
void dumpMediaOnGrid( FieldComponent field )
{
 
  FILE *fp;
  char fileName[PATH_SIZE];

  message( MSG_LOG , 0 , "\nPrinting the grid %s field media array...\n\n" , FIELD[field] ); 
  
  /* Determine file name for each field. */
  switch( field )
  {
  case EX:
    strncpy( fileName , "mediaEx.dat" , PATH_SIZE );
    break;
  case EY:
    strncpy( fileName , "mediaEy.dat" , PATH_SIZE );
    break;
  case EZ:
    strncpy( fileName , "mediaEz.dat" , PATH_SIZE );
    break;
  case HX:
    strncpy( fileName , "mediaHx.dat" , PATH_SIZE );
    break;
  case HY:
    strncpy( fileName , "mediaHy.dat" , PATH_SIZE );
    break;
  case HZ:
    strncpy( fileName , "mediaHz.dat" , PATH_SIZE );
    break;
  default:
   assert( 0 );
  }

  fp = fopen( fileName , "w" );
 
  for ( int k = gobox[ZHI] ; k >= gobox[ZLO] - 1; k-- )
  {
    for ( int j = gobox[YHI] ; j >= gobox[YLO] - 1; j-- )
    {
      for( int space = gobox[YLO] - 1; space < j ; space++ ) fprintf( fp , " " );
      for ( int i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      {
        switch( field )
        {
        case EX:
          fprintf( fp , "%c" , decodeAlphaBeta( ALPHA_EX(i,j,k) , UNSCALE_betaEx( BETA_EX(i,j,k) , i , j , k ) ) );
          break;
        case EY:
          fprintf( fp , "%c" , decodeAlphaBeta( ALPHA_EY(i,j,k) , UNSCALE_betaEy( BETA_EY(i,j,k) , i , j , k ) ) );
          break;
        case EZ:
          fprintf( fp , "%c" , decodeAlphaBeta( ALPHA_EZ(i,j,k) , UNSCALE_betaEz( BETA_EZ(i,j,k) , i , j , k ) ) );
          break;
        case HX:
          fprintf( fp , "%c" , decodeGamma( UNSCALE_gammaHx( GAMMA_HX(i,j,k) , i , j , k ) ) );
          break;
        case HY:
          fprintf( fp , "%c" , decodeGamma( UNSCALE_gammaHy( GAMMA_HY(i,j,k) , i , j , k ) ) );
          break;
        case HZ:
          fprintf( fp , "%c" , decodeGamma( UNSCALE_gammaHz( GAMMA_HZ(i,j,k) , i , j , k ) ) );
          break;
        default:
          assert( 0 );
        }
      }
      if( j == gobox[YHI] - 1 ) fprintf( fp , " y=%d" , j );
      if( j == ( gobox[YHI] + gobox[YLO] - 1 ) / 2 ) fprintf( fp , "       z=%d" , k );
      if( j == gobox[YLO] ) fprintf( fp , " y=%d" , j );
      fprintf( fp , "\n" );
    }
  }

  fclose( fp );

  return;
  
}

/* Decode electric field update coefficient. */
/* '.' for free space, '*' for PEC and 'o' for other. */
char decodeAlphaBeta( real alpha , real beta )
{

  if( fabs( alpha - 1.0 ) < 1e-6 && fabs( beta - dt / eps0 ) < 1e-6 )
    return '.'; /* Free space. */
  else if( fabs( alpha + 1.0 ) < 1e-6 && fabs( beta ) < 1e-6 ) 
    return '*'; /* PEC. */
  else
    return 'o'; /* Other. */

}

/* Decode magnetic field update coefficient. */
/* '.' for free space and 'o' for other. */
char decodeGamma( real gamma )
{

  if( fabs( gamma - dt / mu0 ) < 1e-6 )
    return '.'; /* Free space. */
  else
    return 'o'; /* Other. */

}

/* Set update coefficients in bbox on grid to those of medium. */
void setMediumOnGrid( int gbbox[6] , MediumIndex medium , FaceMask mask )
{

  int flim[6][6];
  FieldComponent field;
  bool includeBoundary[6];
  
  faceMask2boolArray( includeBoundary , mask );
    
  /* Exclude outer surfaces which lie in/on mesh boundaries unless PML/PMC? */
  //for( int surface = XLO ; surface <= ZHI ; surface++ )
  //  if( ( gbbox[surface] == gibox[surface] ) && ( outerSurfaceType( surface ) == BT_PML ) ) includeBoundary[surface] = false;

  switch( bboxType( gbbox ) )
  {
  case BB_VOLUME:
    setFieldLimits( gbbox , flim , includeBoundary );
    break;
  case BB_SURFACE:
    setFieldLimits( gbbox , flim , includeBoundary );
    break;
  case BB_LINE:
    setFieldLimits( gbbox , flim , includeBoundary );
    break;
  case BB_POINT:
    return;
    break;
  default:
    assert( 0 );
    break;
  }

  message( MSG_DEBUG3 , 0 , "    EX FLIM=[%d,%d,%d,%d,%d,%d]\n" , flim[EX][XLO] , flim[EX][XHI] , flim[EX][YLO] , flim[EX][YHI] , flim[EX][ZLO] , flim[EX][ZHI] );
  message( MSG_DEBUG3 , 0 , "    EY FLIM=[%d,%d,%d,%d,%d,%d]\n" , flim[EY][XLO] , flim[EY][XHI] , flim[EY][YLO] , flim[EY][YHI] , flim[EY][ZLO] , flim[EY][ZHI] );
  message( MSG_DEBUG3 , 0 , "    EZ FLIM=[%d,%d,%d,%d,%d,%d]\n" , flim[EZ][XLO] , flim[EZ][XHI] , flim[EZ][YLO] , flim[EZ][YHI] , flim[EZ][ZLO] , flim[EZ][ZHI] );
  message( MSG_DEBUG3 , 0 , "    HX FLIM=[%d,%d,%d,%d,%d,%d]\n" , flim[HX][XLO] , flim[HX][XHI] , flim[HX][YLO] , flim[HX][YHI] , flim[HX][ZLO] , flim[HX][ZHI] );
  message( MSG_DEBUG3 , 0 , "    HY FLIM=[%d,%d,%d,%d,%d,%d]\n" , flim[HY][XLO] , flim[HY][XHI] , flim[HY][YLO] , flim[HY][YHI] , flim[HY][ZLO] , flim[HY][ZHI] );
  message( MSG_DEBUG3 , 0 , "    HZ FLIM=[%d,%d,%d,%d,%d,%d]\n" , flim[HZ][XLO] , flim[HZ][XHI] , flim[HZ][YLO] , flim[HZ][YHI] , flim[HZ][ZLO] , flim[HZ][ZHI] );

  #ifdef USE_INDEXED_MEDIA

    message( MSG_DEBUG3 , 0 , "      Medium number#=%lu\n" , medium );
    
    /* Set update coefficients. */
    field = EX;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
          mediumEx[i][j][k] = medium;

    field = EY;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
          mediumEy[i][j][k] = medium;

    field = EZ;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
          mediumEz[i][j][k] = medium;

    field = HX;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
          mediumHx[i][j][k] = medium;

    field = HY;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
          mediumHy[i][j][k] = medium;

    field = HZ;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
          mediumHz[i][j][k] = medium;
        
  #else
  
    real alpha;
    real beta;
    real gamma;
  
    /* Get medium update coefficients. */
    getSimpleMediumCoefficients( &alpha , &beta , &gamma , medium );
  
    message( MSG_DEBUG3 , 0 , "      Medium#=%lu: alpha=%e beta=%e gamma=%e\n" , medium , alpha , beta , gamma );
        
    /* Scale and set update coefficients. */
    field = EX;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
        {
          alphaEx[i][j][k] = alpha;
          betaEx[i][j][k] = SCALE_betaEx( beta , i , j , k );
        }

    field = EY;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
        {
          alphaEy[i][j][k] = alpha;
          betaEy[i][j][k] = SCALE_betaEy( beta , i , j , k );
        }

    field = EZ;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
        {
          alphaEz[i][j][k] = alpha;
          betaEz[i][j][k] = SCALE_betaEz( beta , i , j , k );
        }

    field = HX;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
        {
          gammaHx[i][j][k] = SCALE_gammaHx( gamma , i , j , k );
        }

    field = HY;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
        {
          gammaHy[i][j][k] = SCALE_gammaHy( gamma , i , j , k );
        }

    field = HZ;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
        {
          gammaHz[i][j][k] = SCALE_gammaHz( gamma , i , j , k );
        }

  #endif

  return;

}

/* Verify update coefficients in bbox on grid are those of medium. */
void checkMediumOnGrid( int gbbox[6] , MediumIndex medium )
{

  bool includeBoundary[6] = {  true ,  true ,  true ,  true ,  true ,  true };
  int flim[6][6];
  FieldComponent field;

  switch( bboxType( gbbox ) )
  {
  case BB_VOLUME:
    setFieldLimits( gbbox , flim , includeBoundary );
    break;
  case BB_SURFACE:
    setFieldLimits( gbbox , flim , includeBoundary );
    break;
  case BB_LINE:
    setFieldLimits( gbbox , flim , includeBoundary );
    break;
  case BB_POINT:
    return;
    break;
  default:
    assert( 0 );
    break;
  }

  #ifdef USE_INDEXED_MEDIA

    /* Check update coefficients. */
    field = EX;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
          assert( mediumEx[i][j][k] == medium );

    field = EY;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
          assert( mediumEy[i][j][k] == medium );

    field = EZ;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
          assert( mediumEz[i][j][k] == medium );

    field = HX;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
          assert( mediumHx[i][j][k] == medium );

    field = HY;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
          assert( mediumHy[i][j][k] == medium );

    field = HZ;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
          assert( mediumHz[i][j][k] == medium );

  #else

    real alpha;
    real beta;
    real gamma;
  
    /* Get medium update coefficients. */
    getSimpleMediumCoefficients( &alpha , &beta , &gamma , medium );
  
    /* Scale and set update coefficients. */
    field = EX;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
        {
          assert( isEqualRel( alphaEx[i][j][k] , alpha , CHECK_LIMITS_RTOL ) );
          assert( isEqualRel( betaEx[i][j][k] , SCALE_betaEx( beta , i , j , k ) , CHECK_LIMITS_RTOL ) );
        }

    field = EY;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
        {
          assert( isEqualRel( alphaEy[i][j][k] , alpha , CHECK_LIMITS_RTOL ) );
          assert( isEqualRel( betaEy[i][j][k] , SCALE_betaEy( beta , i , j , k ) , CHECK_LIMITS_RTOL ) );
        }

    field = EZ;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
        {
          assert( isEqualRel( alphaEz[i][j][k] , alpha , CHECK_LIMITS_RTOL ) );
          assert( isEqualRel( betaEz[i][j][k] , SCALE_betaEz( beta , i , j , k ) , CHECK_LIMITS_RTOL ) );
        }

    field = HX;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
        {
          assert( isEqualRel( gammaHx[i][j][k] , SCALE_gammaHx( gamma , i , j , k ) , CHECK_LIMITS_RTOL ) );
        }

    field = HY;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
        {
          assert( isEqualRel( gammaHy[i][j][k] , SCALE_gammaHy( gamma , i , j , k ) , CHECK_LIMITS_RTOL ) );
        }

    field = HZ;
    for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
      for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
        for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
        {
          assert( isEqualRel( gammaHz[i][j][k] , SCALE_gammaHz( gamma , i , j , k ) , CHECK_LIMITS_RTOL ) );
        }

  #endif

  return;

}

#ifdef USE_AVERAGED_MEDIA

/* Apply voxels to the grid using media averaging. */
void applyVoxelsToGrid( MediumIndex ***blockArray )
{
  
  int field , i , j , k;
  int bbox[6];
  int flim[6][6];
  //bool includeBoundary[6] = { true  ,  true  ,  true  ,  true  ,  true  ,  true  };
  bool includeBoundary[6] = { false  ,  false ,  false  ,  false  ,  false ,  false  };
  real eps_r1 , sigma1 , mu_r1 , eps_r2 , sigma2 , mu_r2;
  real eps_r3 , sigma3 , mu_r3 , eps_r4 , sigma4 , mu_r4;
  real eps_r , sigma , mu_r , alpha , beta , gamma;
  
  /* Copy voxels on grid inner faces into cell within boundary */
  /* so that media on outer surface of mesh are unaveraged by */
  /* the general algorithm below. Note that at least one cell */
  /* is guaranteed to exist outside the inner grid - ghost cell. */
  getFaceOfBoundingBox( bbox , gibox , XLO );
  i = bbox[XLO];
  for( j = bbox[YLO]  ; j <= bbox[YHI] ; j++ )
    for( k = bbox[ZLO]  ; k <= bbox[ZHI] ; k++ )
      blockArray[i-1][j][k] = blockArray[i][j][k];
    
  getFaceOfBoundingBox( bbox , gibox , XHI );
  i = bbox[XHI];
  for( j = bbox[YLO]  ; j <= bbox[YHI] ; j++ )
    for( k = bbox[ZLO]  ; k <= bbox[ZHI] ; k++ )
      blockArray[i][j][k] = blockArray[i-1][j][k];
    
  getFaceOfBoundingBox( bbox , gibox , YLO );
  j = bbox[YLO];
  for( i = bbox[XLO]  ; i <= bbox[XHI] ; i++ )
    for( k = bbox[ZLO]  ; k <= bbox[ZHI] ; k++ )
      blockArray[i][j-1][k] = blockArray[i][j][k];
    
  getFaceOfBoundingBox( bbox , gibox , YHI );
  j = bbox[YHI];
  for( i = bbox[XLO]  ; i <= bbox[XHI] ; i++ )
    for( k = bbox[ZLO]  ; k <= bbox[ZHI] ; k++ )
      blockArray[i][j][k] = blockArray[i][j-1][k];
 
  getFaceOfBoundingBox( bbox , gibox , ZLO );
  k = bbox[ZLO];
  for( i = bbox[XLO]  ; i <= bbox[XHI] ; i++ )
    for( j = bbox[YLO]  ; j <= bbox[YHI] ; j++ )
      blockArray[i][j][k-1] = blockArray[i][j][k];
  
  getFaceOfBoundingBox( bbox , gibox , ZHI );
  k = bbox[ZHI];
  for( i = bbox[XLO]  ; i <= bbox[XHI] ; i++ )
    for( j = bbox[YLO]  ; j <= bbox[YHI] ; j++ )
      blockArray[i][j][k] = blockArray[i][j][k-1];

  /* Field limits.*/
  setFieldLimits( gibox , flim , includeBoundary );
  
  /* Scale and set update coefficients. */
  field = EX;
  for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
    for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
      for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
      {    
        getMediumParameters( &eps_r1 , &sigma1 , &mu_r1 , blockArray[i][j-1][k]   );
        getMediumParameters( &eps_r2 , &sigma2 , &mu_r2 , blockArray[i][j][k]     );
        getMediumParameters( &eps_r3 , &sigma3 , &mu_r3 , blockArray[i][j][k-1]   );
        getMediumParameters( &eps_r4 , &sigma4 , &mu_r4 , blockArray[i][j-1][k-1] );
        eps_r = 0.25 * ( eps_r1 * dey[j-1] * dez[k] + eps_r2 * dey[j] * dez[k] + 
                         eps_r3 * dey[j] * dez[k-1] + eps_r4 * dey[j-1] * dez[k-1] ) / ( dhy[j] * dhz[k] ); 
        sigma = 0.25 * ( sigma1 * dey[j-1] * dez[k] + sigma2 * dey[j] * dez[k] + 
                         sigma3 * dey[j] * dez[k-1] + sigma4 * dey[j-1] * dez[k-1] ) / ( dhy[j] * dhz[k] );
        calcCoeffFromParam( dt , eps_r , sigma , 1.0 , &alpha , &beta , &gamma );
        alphaEx[i][j][k] = alpha;
        betaEx[i][j][k] = SCALE_betaEx( beta , i , j , k );
      }

  field = EY;
  for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
    for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
      for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
      {
        getMediumParameters( &eps_r1 , &sigma1 , &mu_r1 , blockArray[i][j][k-1]   );
        getMediumParameters( &eps_r2 , &sigma2 , &mu_r2 , blockArray[i][j][k]     );
        getMediumParameters( &eps_r3 , &sigma3 , &mu_r3 , blockArray[i-1][j][k]   );
        getMediumParameters( &eps_r4 , &sigma4 , &mu_r4 , blockArray[i-1][j][k-1] );
        eps_r = 0.25 * ( eps_r1 * dez[k-1] * dex[i] + eps_r2 * dez[k] * dex[i] + 
                         eps_r3 * dez[k] * dex[i-1] + eps_r4 * dez[k-1] * dex[i-1] ) / ( dhz[k] * dhx[i] ); 
        sigma = 0.25 * ( sigma1 * dez[k-1] * dex[i] + sigma2 * dez[k] * dex[i] + 
                         sigma3 * dez[k] * dex[i-1] + sigma4 * dez[k-1] * dex[i-1] ) / ( dhz[k] * dhx[i] ); 
        calcCoeffFromParam( dt , eps_r , sigma , 1.0 , &alpha , &beta , &gamma );
        alphaEy[i][j][k] = alpha;
        betaEy[i][j][k] = SCALE_betaEy( beta , i , j , k );
      }

  field = EZ;
  for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
    for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
      for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
      {
        getMediumParameters( &eps_r1 , &sigma1 , &mu_r1 , blockArray[i-1][j][k]   );
        getMediumParameters( &eps_r2 , &sigma2 , &mu_r2 , blockArray[i][j][k]     );
        getMediumParameters( &eps_r3 , &sigma3 , &mu_r3 , blockArray[i][j-1][k]   );
        getMediumParameters( &eps_r4 , &sigma4 , &mu_r4 , blockArray[i-1][j-1][k] );
        eps_r = 0.25 * ( eps_r1 * dex[i-1] * dey[j] + eps_r2 * dex[i] * dey[j] + 
                         eps_r3 * dex[i] * dey[j-1] + eps_r4 * dex[i-1] * dey[j-1] ) / ( dhx[i] * dhy[j] ); 
        sigma = 0.25 * ( sigma1 * dex[i-1] * dey[j] + sigma2 * dex[i] * dey[j] + 
                         sigma3 * dex[i] * dey[j-1] + sigma4 * dex[i-1] * dey[j-1] ) / ( dhx[i] * dhy[j] );
        calcCoeffFromParam( dt , eps_r , sigma , 1.0 , &alpha , &beta , &gamma );
        alphaEz[i][j][k] = alpha;
        betaEz[i][j][k] = SCALE_betaEz( beta , i , j , k );
      }

  field = HX;
  for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
    for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
      for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
      {
        getMediumParameters( &eps_r1 , &sigma1 , &mu_r1 , blockArray[i-1][j][k] );
        getMediumParameters( &eps_r2 , &sigma2 , &mu_r2 , blockArray[i][j][k] );
        mu_r = 2.0 * mu_r1 * mu_r2 / ( mu_r1 + mu_r2 ); 
        calcCoeffFromParam( dt , 1.0 , 0.0 , mu_r , &alpha , &beta , &gamma );
        gammaHx[i][j][k] = SCALE_gammaHx( gamma , i , j , k );
      }

  field = HY;
  for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
    for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
      for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
      {
        getMediumParameters( &eps_r1 , &sigma1 , &mu_r1 , blockArray[i][j-1][k] );
        getMediumParameters( &eps_r2 , &sigma2 , &mu_r2 , blockArray[i][j][k] );
        mu_r = 2.0 * mu_r1 * mu_r2 / ( mu_r1 + mu_r2 ); 
        calcCoeffFromParam( dt , 1.0 , 0.0 , mu_r , &alpha , &beta , &gamma );
        gammaHy[i][j][k] = SCALE_gammaHy( gamma , i , j , k );
      }

  field = HZ;
  for( int i = flim[field][XLO]  ; i <= flim[field][XHI] ; i++ )
    for( int j = flim[field][YLO]  ; j <= flim[field][YHI] ; j++ )
      for( int k = flim[field][ZLO]  ; k <= flim[field][ZHI] ; k++ )
      {
        getMediumParameters( &eps_r1 , &sigma1 , &mu_r1 , blockArray[i][j][k-1] );
        getMediumParameters( &eps_r2 , &sigma2 , &mu_r2 , blockArray[i][j][k] );
        mu_r = 2.0 * mu_r1 * mu_r2 / ( mu_r1 + mu_r2 ); 
        calcCoeffFromParam( dt , 1.0 , 0.0 , mu_r , &alpha , &beta , &gamma );
        gammaHz[i][j][k] = SCALE_gammaHz( gamma , i , j , k );
      }

  return;

}

#endif // USE_AVERAGED_MEDIA
