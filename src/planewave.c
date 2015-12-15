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

#include <assert.h>
#include <math.h>
#include <string.h>

#include "planewave.h"
#include "utlist.h"
#include "waveform.h"
#include "alloc_array.h"
#include "message.h"
#include "bounding_box.h"
#include "gnuplot.h"
#include "gmsh.h"
#include "grid.h"
#include "medium.h"
#include "surface.h"
#include "physical.h"
#include "util.h"

/* 
 * Plane wave class. 
 */

typedef struct PlaneWaveItem_t {

  PlaneWaveIndex number;           // Source number, assigned in order found.

  /* Parameters from mesh. */
  char name[TAG_SIZE];              // Observer name from mesh.
  int mbbox[6];                   // Bounding box on mesh.
  WaveformIndex waveformNumber;   // Waveform number.
  bool isActive[6];               // Mask for which box faces are active.
  real theta;                     // Indicent altitude angles [degrees].
  real phi;                       // Incident azimuthal angle [degress].
  real eta;                       // Polarisation angle [rad].
  real size;                      // Amplitude.
  real delay;                     // Delay [s].
 
  /* Derived parameters. */
  int gbbox[6];                   // Bounding box on grid.
  int flim[6][6][6];              // Field limits for each face of TF/SF box.
  real kinc[3];                   // Incident unit wave vector.
  real Finc[6];                   // Incident electric and magnetic field vectors.
  real ijk0[3];                   // Origin of planewave in grid indices.
  real r0[3];                     // Origin of planewave in physical units.
  real phaseVelocity;             // Free-space numerical phase velocity on grid.

  /* Auxiliary grid parameters. */
  int nx;                         // Length of grid.
  real *Eyi;                      // Field arrays.
  real *Hzi;
  real betaEyi;                   // Material property arrays.
  real gammaHzi;
  real *Pyi;                      // PML EM Field Arrays.
  real *PPyi;
  real *Bzi; 
  real *adx;                      // PML Loss Profile Arrays.
  real *bdx;
  real *ahx;
  real *bhx;
  int xb;                         // PML Boundary position.

  /* UT list. */

  struct PlaneWaveItem_t *prev;
  struct PlaneWaveItem_t *next;

  UT_hash_handle hh;              // Plane-wave name hash.
  
} PlaneWaveItem;

/* 
 * Private data. 
 */

/* Number of plane waves. */
static PlaneWaveIndex numPlaneWave = 0;

/* List of plane waves. */
static PlaneWaveItem *planeWaveList;    

/* Hash of plane-waves using name. */
static PlaneWaveItem *planeWaveHash = NULL;

/* Use auxiliary grid if true. */
bool useAuxGrid = false;

/* Null phase point of auxiliary grids. */
static int m0 = 2;

/* PML Depth. */
static int npml = 10;

/* Function pointer for indicent field function. */
static real (*incidentField)( FieldComponent field , int i , int j , int k , real time , PlaneWaveItem *item ) = NULL;

/* 
 * Private method interfaces. 
 */

void addPlaneWave( int mbbox[6] , char name[TAG_SIZE] , bool isActive[6] , real theta , real phi , real eta , 
                   real size , real delay , WaveformIndex waveformNumber );
bool decodeFaceMask( bool isActive[6] , char maskStr[] );
real incidentFieldAnalytic( FieldComponent field , int i , int j , int k , real time , PlaneWaveItem *item );
real incidentFieldAuxGrid( FieldComponent field , int i , int j , int k , real time , PlaneWaveItem *item );
void updateAuxGridHfield( PlaneWaveItem *item , real time );
void updateAuxGridEfield( PlaneWaveItem *item , real time);
void deallocAuxGrid( PlaneWaveItem *item );
bool isPlaneWave( char *name , PlaneWaveIndex *number );
void initAuxGrid( PlaneWaveItem *item );
void calcIncidentFieldVectors( real kinc[3] , real Finc[6] , real ijk0[3] , int gbbox[6] ,
                               real size , real theta , real phi , real eta );

/*
 * Method Implementations.
 */

/* Add plane wave to lists. */
void addPlaneWave( int mbbox[6] , char name[TAG_SIZE] , bool isActive[6] , real theta , real phi , real eta , 
                   real size , real delay , WaveformIndex waveformNumber )
{

  PlaneWaveItem *item = NULL;

  if( numPlaneWave == MAX_PLANE_WAVE )
    message( MSG_ERROR , 0 , "*** Error: Maximum number of plane waves exceeded!\n" );

  item = (PlaneWaveItem *) malloc( sizeof( PlaneWaveItem ) );
  if( !item )
    message( MSG_ERROR , 0 , "*** Error: Failed to allocate plane wave item" );

  strncpy( item->name , name , TAG_SIZE );
  item->number = numPlaneWave;
  item->waveformNumber = waveformNumber;
  for( int face = XLO ; face <= ZHI ; face++ ) item->mbbox[face] = mbbox[face];
  for( int face = XLO ; face <= ZHI ; face++ ) item->isActive[face] = isActive[face];
  item->theta = theta;
  item->phi = phi;
  item->eta = eta;
  item->size = size;
  item->delay = delay;

  /* Add to list. */
  DL_APPEND( planeWaveList , item );
  HASH_ADD_STR( planeWaveHash , name , item );
  numPlaneWave++;

  return;

}

/* Parse plane waves. */
bool parsePW( char *line )
{
  
  int numScanned = 0;
  char waveformName[TAG_SIZE] = "";
  int mbbox[6];
  char name[TAG_SIZE];
  char maskStr[TAG_SIZE] = "111111";
  bool isActive[6] = { true , true , true , true , true , true };
  double size = 1.0;
  double delay = 0.0;
  double theta = 0.0;
  double phi = 0.0;
  double eta = 0.0;
  WaveformIndex waveformNumber;
  PlaneWaveIndex planeWaveNumber;

  numScanned = sscanf( line , "%d %d %d %d %d %d %31s %31s %lf %lf %lf %31s %lf %lf" , 
                       &mbbox[XLO] , &mbbox[XHI] , &mbbox[YLO] , &mbbox[YHI] , &mbbox[ZLO] , &mbbox[ZHI] , 
                       name , waveformName , &theta , &phi , &eta , maskStr , &size , &delay );

  if( numScanned < 11 )
    return false;  

  /* Check plane-wave is not already defined. */
  if( isPlaneWave( name , &planeWaveNumber ) )
  {
    message( MSG_LOG , 0 , "  Plane-wave %s already defined\n" , name );
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

  /* Check waveform exists. */
  if( !isWaveform( waveformName , &waveformNumber ) )
  {
    message( MSG_LOG , 0 , "  Waveform %s not defined in plane wave card\n" , waveformName );
    return false;
  }

  /* Validate angles. */
  if( theta < 0.0 || theta > 180.0 )
  {
    message( MSG_LOG , 0 , "\n  Incident angle theta must be in range [0,180] degrees.\n" );
    return false;
  }
  if( phi < 0.0 || phi > 360.0 )
  {
    message( MSG_LOG , 0 , "\n  Incident angle phi must be in range [0,360) degrees.\n" );
    return false;
  }
  if( eta < 0 || eta >= 360.0 )
  {
    message( MSG_LOG , 0 , "\n  Polarisation angle eta must be in range [0,360) degrees.\n" );
    return false;
  }

  /* Decode mask. */
  if( numScanned >= 12 )
  {
    if( strlen( maskStr ) != 6 )
    {
      message( MSG_LOG , 0 , "  Face mask %s in plane wave card must have exactly six bits!\n" , maskStr );
      return false;    
    }
    else if( !decodeFaceMask( isActive , maskStr ) )
    {
      message( MSG_LOG , 0 , "  Face mask %s invalid in plane wave card\n" , maskStr );
      return false;
    }
  }

  /* Validate size. */
  if( numScanned >= 13 && size < 0.0 )
  {
    message( MSG_LOG , 0 , "  Waveform size must be positive:\n" );
    return false;
  }

  /* Validate delay. */
  if(  numScanned >= 14 && delay < 0.0 )
    message( MSG_WARN , 0 , "Waveform delay negative:\n" );

  addPlaneWave( mbbox , name , isActive, theta , phi , eta , size , delay , waveformNumber );

  return true;

}

/* Initialise sources. */
/* Depends: initGrid, initWaveforms */
void initPlaneWaves( void )
{

  PlaneWaveItem *item;
  GridType gridType;       // Type of grid.
  int bbox[6];             // Region bounding box.
  bool includeBoundary[6]; // Region boundary inclusion flags. 
  bool edgeIsActive[6];    // Edge boundary inclusion flags. 

  message( MSG_LOG , 0 , "\nInitialising plane waves...\n\n" );

  message( MSG_DEBUG1 , 0 , "  Allocating plane wave array\n" );

  /* Set phase velocity of plane wave. */
  gridType = getGridType();
  switch( gridType )
  {
  case GT_CUBIC:
    useAuxGrid = true;
    incidentField = &incidentFieldAuxGrid;
    message( MSG_LOG , 0 , "  Setting plane wave auxiliary grid incident field calculation\n" );
    break;
  case GT_UNIFORM:
    incidentField = &incidentFieldAnalytic;
    message( MSG_LOG , 0 , "  Setting plane wave analytic incident field calculation\n" );
    break;
  case GT_NONUNIFORM:   
    incidentField = &incidentFieldAnalytic;
    message( MSG_LOG , 0 , "  Setting plane wave analytic incident field calculation\n" );
    break;
  default:
    assert( 0 );
    break;
  }

  /* Iterate over plane waves. */
  DL_FOREACH( planeWaveList , item ) 
  {
            
    offsetBoundingBox( item->gbbox , item->mbbox , gibox );

    message( MSG_DEBUG3 , 0 , "  Setting plane wave \"%s\" on [%d,%d,%d,%d,%d,%d]/[%d,%d,%d,%d,%d,%d]: mask=[%d,%d,%d,%d,%d,%d] dir=(%.0f,%.0f) pol=%.0f size=%g, delay=%g\n" ,
             item->name ,
             item->mbbox[XLO] , item->mbbox[XHI] , item->mbbox[YLO] , item->mbbox[YHI] , 
             item->mbbox[ZLO] , item->mbbox[ZHI] , item->gbbox[XLO] , item->gbbox[XHI] , 
             item->gbbox[YLO] , item->gbbox[YHI] , item->gbbox[ZLO] , item->gbbox[ZHI] , 
             item->isActive[XLO] , item->isActive[XHI] , item->isActive[YLO] , 
             item->isActive[YHI] , item->isActive[ZLO] , item->isActive[ZHI] , 
             item->theta , item->phi , item->eta , item->size, item->delay );

    /* Determine incident field origin and vectors. */
    calcIncidentFieldVectors( item->kinc , item->Finc , item->ijk0 , item->gbbox , 
                              item->size , item->theta , item->phi , item->eta );

    /* Origin in physical units. */
    getNodeLocation( item->r0 , item->ijk0[XDIR] , item->ijk0[YDIR] , item->ijk0[ZDIR] );
    
    message( MSG_DEBUG3 , 0 , "    uinc=(%e,%e,%e) [-], r0=(%e,%e,%e) [m]\n" ,
             item->kinc[XDIR] , item->kinc[YDIR] , item->kinc[ZDIR] ,  
             item->r0[XDIR]   , item->r0[YDIR]   , item->r0[ZDIR]   );
    message( MSG_DEBUG3 , 0 , "    Einc =(%e,%e,%e) [V/m], Hinc=(%e,%e,%e) [A/m]\n" ,
             item->Finc[EX] , item->Finc[EY] , item->Finc[EZ] , 
             item->Finc[HX] , item->Finc[HY] , item->Finc[HZ] );
    
    /* Determine the method used to calculate incident field. */
    switch( gridType )
    {
    case GT_CUBIC:
      initAuxGrid( item );    
      break;
    case GT_UNIFORM:
      item->phaseVelocity = numericalPhaseVelocity( degrees2radians( item->theta ) , degrees2radians( item->phi ) );
      message( MSG_DEBUG3 , 0 , "    Numerical phase velocity=%g*c0\n" , item->phaseVelocity / c0 );
      break;
    case GT_NONUNIFORM:   
      item->phaseVelocity = c0;   
      message( MSG_DEBUG3 , 0 , "    Numerical phase velocity=%g*c0\n" , item->phaseVelocity / c0 );
      break;
    default:
      assert( 0 );
      break;
    }
         
    /* Edge activity flags. If edge fields are on PMC external surface we want to keep them on. */
    for( MeshFace face = XLO; face <= ZHI ; face++ )
    {
      edgeIsActive[face] = item->isActive[face];
      if( ( item->gbbox[face] == gibox[face] ) && ( ( outerSurfaceType( face ) == BT_PMC ) || ( outerSurfaceType( face ) == BT_PERIODIC ) ) )
        edgeIsActive[face] = true;
    }
    
    /* Field limits for field corrections. */
    /* XLO. */
    setBoundingBoxFromNodes( bbox , item->gbbox[XLO] - 1 , item->gbbox[XLO] , 
                                    item->gbbox[YLO]     , item->gbbox[YHI] , 
                                    item->gbbox[ZLO]     , item->gbbox[ZHI] );
    setBoundingBoxBoundaryFlags( includeBoundary , false             , true              ,
                                                   edgeIsActive[YLO] , edgeIsActive[YHI] , 
                                                   edgeIsActive[ZLO] , edgeIsActive[ZHI] );
    setFieldLimits( bbox , item->flim[XLO] , includeBoundary );

    /* XHI. */
    setBoundingBoxFromNodes( bbox , item->gbbox[XHI] , item->gbbox[XHI] + 1 , 
                                    item->gbbox[YLO] , item->gbbox[YHI]     , 
                                    item->gbbox[ZLO] , item->gbbox[ZHI]     );
    setBoundingBoxBoundaryFlags( includeBoundary , true              , false             ,
                                                   edgeIsActive[YLO] , edgeIsActive[YHI] , 
                                                   edgeIsActive[ZLO] , edgeIsActive[ZHI] );
    setFieldLimits( bbox , item->flim[XHI] , includeBoundary );
    
     /* YLO. */
    setBoundingBoxFromNodes( bbox , item->gbbox[XLO]     , item->gbbox[XHI] , 
                                    item->gbbox[YLO] - 1 , item->gbbox[YLO] , 
                                    item->gbbox[ZLO]     , item->gbbox[ZHI] );
    setBoundingBoxBoundaryFlags( includeBoundary , edgeIsActive[XLO] , edgeIsActive[XHI] , 
                                                   false             , true              , 
                                                   edgeIsActive[ZLO] , edgeIsActive[ZHI] );
    setFieldLimits( bbox , item->flim[YLO] , includeBoundary );
 
    /* YHI. */
    setBoundingBoxFromNodes( bbox , item->gbbox[XLO] , item->gbbox[XHI]     , 
                                    item->gbbox[YHI] , item->gbbox[YHI] + 1 , 
                                    item->gbbox[ZLO] , item->gbbox[ZHI]     );
    setBoundingBoxBoundaryFlags( includeBoundary , edgeIsActive[XLO] , edgeIsActive[XHI] , 
                                                   true              , false             , 
                                                   edgeIsActive[ZLO] , edgeIsActive[ZHI] );
    setFieldLimits( bbox , item->flim[YHI] , includeBoundary );
    
     /* ZLO. */
    setBoundingBoxFromNodes( bbox , item->gbbox[XLO]     , item->gbbox[XHI] , 
                                    item->gbbox[YLO]     , item->gbbox[YHI] , 
                                    item->gbbox[ZLO] - 1 , item->gbbox[ZLO] );
    setBoundingBoxBoundaryFlags( includeBoundary , edgeIsActive[XLO] , edgeIsActive[XHI] , 
                                                   edgeIsActive[YLO] , edgeIsActive[YHI] ,
                                                   false             , true              );
    setFieldLimits( bbox , item->flim[ZLO] , includeBoundary );
    
     /* ZHI. */
    setBoundingBoxFromNodes( bbox , item->gbbox[XLO] , item->gbbox[XHI]     , 
                                    item->gbbox[YLO] , item->gbbox[YHI]     , 
                                    item->gbbox[ZHI] , item->gbbox[ZHI] + 1 );
    setBoundingBoxBoundaryFlags( includeBoundary , edgeIsActive[XLO] , edgeIsActive[XHI] , 
                                                   edgeIsActive[YLO] , edgeIsActive[YHI] ,
                                                   true              , false             );
    setFieldLimits( bbox , item->flim[ZHI] , includeBoundary );    

    /* Log field limits. */
    for( MeshFace face = XLO ; face <= ZHI ; face++ )
    {
      for( FieldComponent field = EX ; field <= HZ ; field++ )
      {
        if( fieldIsParallelToBoundary( field , face ) ) 
        {
          message( MSG_DEBUG1 , 0 , "    Face %s Field %s Limits: [%d,%d,%d,%d,%d,%d]\n",
                   FACE[face] , FIELD[field] ,
                   item->flim[face][field][XLO] , item->flim[face][field][XHI] ,
                   item->flim[face][field][YLO] , item->flim[face][field][YHI] ,
                   item->flim[face][field][ZLO] , item->flim[face][field][ZHI] );
        }
        else
        {
          continue; 
        }
      }
    }
              
  } //   DL_FOREACH( planeWaveList , item ) 
    
  return;

}

/* Calculate plane wave origin and vectors. */
void calcIncidentFieldVectors( real kinc[3] , real Finc[6] , real ijk0[3] , int gbbox[6] , 
                               real size , real theta , real phi , real eta )
{

  real theta_rad;    // Incident angle in radians.
  real phi_rad;      // Incident angle in radians.
  real eta_rad;      // Incident polarisation angle in radians.
  
  /* Incident angles and polarisation angle in radians. */
  theta_rad = degrees2radians( theta );
  phi_rad = degrees2radians(  phi );
  eta_rad = degrees2radians(  eta );

  /* Incident direction - eqn. (5.60). */
  kinc[XDIR] = sin( theta_rad ) * cos( phi_rad );
  kinc[YDIR] = sin( theta_rad ) * sin( phi_rad );
  kinc[ZDIR] = cos( theta_rad );
  
  /* Incident fields vector components - eqn. (5.63) . */
  Finc[EX] = size * (  cos( eta_rad ) * sin( phi_rad ) - sin( eta_rad ) * cos( theta_rad ) * cos( phi_rad ) );
  Finc[EY] = size * ( -cos( eta_rad ) * cos( phi_rad ) - sin( eta_rad ) * cos( theta_rad ) * sin( phi_rad ) );
  Finc[EZ] = size * (  sin( eta_rad ) * sin( theta_rad ) );
  Finc[HX] = size / eta0 * (  sin( eta_rad ) * sin( phi_rad ) + cos( eta_rad ) * cos( theta_rad ) * cos( phi_rad ) );
  Finc[HY] = size / eta0 * ( -sin( eta_rad ) * cos( phi_rad ) + cos( eta_rad ) * cos( theta_rad ) * sin( phi_rad ) );
  Finc[HZ] = size / eta0 * ( -cos( eta_rad ) * sin( theta_rad ) );

  /* Determine origin of plane wave. */
  if( theta >= 0 && theta <= 90 )
  {
    if( phi >= 0 && phi <= 90 )
    {
      /* Eqn. (5.61a). */
      ijk0[XDIR] = gbbox[XLO];
      ijk0[YDIR] = gbbox[YLO];
      ijk0[ZDIR] = gbbox[ZLO];
    }
    else if( phi > 90 && phi <= 180 )
    {
      /* Eqn. (5.61b). */
      ijk0[XDIR] = gbbox[XHI];
      ijk0[YDIR] = gbbox[YLO];
      ijk0[ZDIR] = gbbox[ZLO];
    }
    else if( phi > 180 && phi <= 270 )
    {
      /* Eqn. (5.61c). */
      ijk0[XDIR] = gbbox[XHI];
      ijk0[YDIR] = gbbox[YHI];
      ijk0[ZDIR] = gbbox[ZLO];
    }
    else if( phi > 270 && phi < 360 )
    {
      /* Eqn. (5.61d). */
      ijk0[XDIR] = gbbox[XLO];
      ijk0[YDIR] = gbbox[YHI];
      ijk0[ZDIR] = gbbox[ZLO];
    }
    else
    {
      assert( 0 );
    }
  }
  else if( theta > 90 && theta <= 180 )
  {
    if( phi >= 0 && phi <= 90 )
    {
      /* Eqn. (5.62a). */
      ijk0[XDIR] = gbbox[XLO];
      ijk0[YDIR] = gbbox[YLO];
      ijk0[ZDIR] = gbbox[ZHI];
    }
    else if( phi > 90 && phi <= 180 )
    {
      /* Eqn. (5.62b). */
      ijk0[XDIR] = gbbox[XHI];
      ijk0[YDIR] = gbbox[YLO];
      ijk0[ZDIR] = gbbox[ZHI];
    }
    else if( phi > 180 && phi <= 270 )
    {
      /* Eqn. (5.62c). */
      ijk0[XDIR] = gbbox[XHI];
      ijk0[YDIR] = gbbox[YHI];
      ijk0[ZDIR] = gbbox[ZHI];
    }
    else if( phi > 270 && phi < 360 )
    {
      /* Eqn. (5.62d). */
      ijk0[XDIR] = gbbox[XLO];
      ijk0[YDIR] = gbbox[YHI];
      ijk0[ZDIR] = gbbox[ZHI];
    }
    else
    {
      assert( 0 );
    }
  }
  else
  {
    assert( 0 );
  }
    
  return;

}

/* Determine incient field component in cell (i,j,k) at time t. */
real incidentFieldAnalytic( FieldComponent field , int i , int j , int k , real time , PlaneWaveItem *item )
{

  real rcomp[3];
  real d;

  /* Physical location of field point. */
  getFieldPhysicalLocation( rcomp , field , i , j , k );

  /* Eqn. (5.41). */
  d = item->kinc[XDIR] * ( rcomp[XDIR] - item->r0[XDIR] ) + 
      item->kinc[YDIR] * ( rcomp[YDIR] - item->r0[YDIR] ) + 
      item->kinc[ZDIR] * ( rcomp[ZDIR] - item->r0[ZDIR] );

  return item->Finc[field] * getWaveformValue( time - d / item->phaseVelocity , item->waveformNumber , item->delay );

}

/* Apply electric field plane wave correction. */
void updatePlaneWavesEfield( real timeE )
{

  PlaneWaveItem *item;
  int i , j , k;
  real incField;

  DL_FOREACH( planeWaveList , item ) 
  {

    if( useAuxGrid )
      updateAuxGridEfield( item , timeE );
      
    if( item->isActive[YLO] )
    {
      /* YLO face, EX - eqn. (5.48a). */
      j = item->flim[YLO][EX][YLO];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , k , incField )
      #endif
      for ( i = item->flim[YLO][EX][XLO] ; i <= item->flim[YLO][EX][XHI] ; i++ ) 
      {
        for ( k = item->flim[YLO][EX][ZLO] ; k <= item->flim[YLO][EX][ZHI] ; k++ )
        {
          incField = SCALE_Hz( incidentField( HZ , i , j - 1 , k , timeE , item ) , k );
          Ex[i][j][k] = ALPHA_EX(i,j,k) * Ex[i][j][k] - BETA_EX(i,j,k) * dHz_dy( incField , j );  
        }
      }

      /* YLO face, EZ - eqn. (5.48b). */
      j = item->flim[YLO][EZ][YLO];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , k , incField )
      #endif
      for ( i = item->flim[YLO][EZ][XLO] ; i <= item->flim[YLO][EZ][XHI] ; i++ ) 
      {
        for ( k = item->flim[YLO][EZ][ZLO] ; k <= item->flim[YLO][EZ][ZHI] ; k++ )
        {
          incField = SCALE_Hx( incidentField( HX , i , j - 1 , k , timeE , item ) , i );
          Ez[i][j][k] = ALPHA_EZ(i,j,k) * Ez[i][j][k] + BETA_EZ(i,j,k) * dHx_dy( incField , j );
        }
      }

    } // if YLO.

    if( item->isActive[YHI] )
    {
      /* YHI face, EX - eqn. (5.49a). */
      j = item->flim[YHI][EX][YHI];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , k , incField )
      #endif
      for ( i = item->flim[YHI][EX][XLO] ; i <= item->flim[YHI][EX][XHI] ; i++ ) 
      {
        for ( k = item->flim[YHI][EX][ZLO] ; k <= item->flim[YHI][EX][ZHI] ; k++ )
        {
          incField = SCALE_Hz( incidentField( HZ , i , j , k , timeE , item ) , k );
          Ex[i][j][k] = ALPHA_EX(i,j,k) * Ex[i][j][k] + BETA_EX(i,j,k) * dHz_dy( incField , j );  
        }
      }

      /* YHI face, EZ - eqn. (5.49b). */
      j = item->flim[YHI][EZ][YHI];      
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , k , incField )
      #endif
      for ( i = item->flim[YHI][EZ][XLO] ; i <= item->flim[YHI][EZ][XHI] ; i++ ) 
      {
        for ( k = item->flim[YHI][EZ][ZLO] ; k <= item->flim[YHI][EZ][ZHI] ; k++ )
        {
          incField = SCALE_Hx( incidentField( HX , i , j , k , timeE , item ) , i );
          Ez[i][j][k] = ALPHA_EZ(i,j,k) * Ez[i][j][k] - BETA_EZ(i,j,k) * dHx_dy( incField , j );  
        }
      }

    } // if YHI.

    if( item->isActive[ZLO] )
    {

      /* ZLO face, EX - eqn. (5.50a). */
      k = item->flim[ZLO][EX][ZLO];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , j , incField )
      #endif
      for ( i = item->flim[ZLO][EX][XLO] ; i <= item->flim[ZLO][EX][XHI] ; i++ ) 
      {
        for ( j = item->flim[ZLO][EX][YLO] ; j <= item->flim[ZLO][EX][YHI] ; j++ )
        {
          incField = SCALE_Hy( incidentField( HY , i , j , k - 1  , timeE , item ) , j );
          Ex[i][j][k] = ALPHA_EX(i,j,k) * Ex[i][j][k] + BETA_EX(i,j,k) * dHy_dz( incField , k );
        }
      }

      /* ZLO face, EY - eqn. (5.50b). */
      k = item->flim[ZLO][EY][ZLO];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , j , incField )
      #endif
      for ( i = item->flim[ZLO][EY][XLO] ; i <= item->flim[ZLO][EY][XHI] ; i++ ) 
      {
        for ( j = item->flim[ZLO][EY][YLO] ; j <= item->flim[ZLO][EY][YHI] ; j++ )
        {
          incField = SCALE_Hx( incidentField( HX , i , j , k - 1  , timeE , item ) , i );
          Ey[i][j][k] = ALPHA_EY(i,j,k) * Ey[i][j][k] - BETA_EY(i,j,k) * dHx_dz( incField , k ); 
        }
      }

    } // if ZLO.

    if( item->isActive[ZHI] )
    {

      /* ZHI face, EX - eqn. (5.51a). */
      k = item->flim[ZHI][EX][ZHI];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , j , incField )
      #endif
      for ( i = item->flim[ZHI][EX][XLO] ; i <= item->flim[ZHI][EX][XHI] ; i++ ) 
      {
        for ( j = item->flim[ZHI][EX][YLO] ; j <= item->flim[ZHI][EX][YHI] ; j++ )
        {
          incField = SCALE_Hy( incidentField( HY , i , j , k , timeE , item ) , j );
          Ex[i][j][k] = ALPHA_EX(i,j,k) * Ex[i][j][k] - BETA_EX(i,j,k) * dHy_dz( incField , k );
        }
      }

      /* ZHI face, EY - eqn. (5.51b). */
      k = item->flim[ZHI][EY][ZHI];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , j , incField )
      #endif
      for ( i = item->flim[ZHI][EY][XLO] ; i <= item->flim[ZHI][EY][XHI] ; i++ ) 
      {
        for ( j = item->flim[ZHI][EY][YLO] ; j <= item->flim[ZHI][EY][YHI] ; j++ )
        {
          incField = SCALE_Hx( incidentField( HX , i , j , k , timeE , item ) , i );
          Ey[i][j][k] = ALPHA_EY(i,j,k) * Ey[i][j][k] + BETA_EY(i,j,k) * dHx_dz( incField , k ); 
        }
      }
  
    } // if ZHI.

    if( item->isActive[XLO] )
    {     

      /* XLO face, EY - eqn. (5.52a). */
      i = item->flim[XLO][EY][XLO];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( j , k , incField )
      #endif
      for ( j = item->flim[XLO][EY][YLO] ; j <= item->flim[XLO][EY][YHI] ; j++ )
      {
        for ( k = item->flim[XLO][EY][ZLO] ; k <= item->flim[XLO][EY][ZHI] ; k++ )
        {
          incField = SCALE_Hz( incidentField( HZ , i - 1  , j , k , timeE , item ) , k );
          Ey[i][j][k] = ALPHA_EY(i,j,k) * Ey[i][j][k] + BETA_EY(i,j,k) * dHz_dx( incField , i );
        }
      }

      /* XLO face, EZ - eqn. (5.52b). */
      i = item->flim[XLO][EZ][XLO];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( j , k , incField )
      #endif
      for ( j = item->flim[XLO][EZ][YLO] ; j <= item->flim[XLO][EZ][YHI] ; j++ )
      {
        for ( k = item->flim[XLO][EZ][ZLO] ; k <= item->flim[XLO][EZ][ZHI] ; k++ )
        {
          incField = SCALE_Hy( incidentField( HY , i - 1  , j , k , timeE , item ) , j );
          Ez[i][j][k] = ALPHA_EZ(i,j,k) * Ez[i][j][k] - BETA_EZ(i,j,k) * dHy_dx( incField , i );
        }
      }

    } // if XLO.

    if( item->isActive[XHI] )
    {  

      /* XHI face, EY - eqn. (5.53a). */
      i = item->flim[XHI][EY][XHI];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( j , k , incField )
      #endif
      for ( j = item->flim[XHI][EY][YLO] ; j <= item->flim[XHI][EY][YHI] ; j++ )
      {
        for ( k = item->flim[XHI][EY][ZLO] ; k <= item->flim[XHI][EY][ZHI] ; k++ )
        {
          incField = SCALE_Hz( incidentField( HZ , i , j , k , timeE , item ) , k );
          Ey[i][j][k] = ALPHA_EY(i,j,k) * Ey[i][j][k] - BETA_EY(i,j,k) * dHz_dx( incField , i );
        }
      }

      /* XHI face, EZ - eqn. (5.53b). */
      i = item->flim[XHI][EZ][XHI];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( j , k , incField )
      #endif
      for ( j = item->flim[XHI][EZ][YLO] ; j <= item->flim[XHI][EZ][YHI] ; j++ )
      {
        for ( k = item->flim[XHI][EZ][ZLO] ; k <= item->flim[XHI][EZ][ZHI] ; k++ )
        {
          incField = SCALE_Hy( incidentField( HY , i , j , k , timeE , item ) , j );
          Ez[i][j][k] = ALPHA_EZ(i,j,k) * Ez[i][j][k] + BETA_EZ(i,j,k) * dHy_dx( incField , i );
        }
      }

    } // if XHI.

  } // DL_FOREACH( planeWaveList , item ) 

  return;

}

/* Apply magnetic field plane wave correction. */
void updatePlaneWavesHfield( real timeH )
{

  PlaneWaveItem *item;
  int i , j , k;
  real incField;

  DL_FOREACH( planeWaveList , item ) 
  {

    if( useAuxGrid )
      updateAuxGridHfield( item , timeH );
    
    if( item->isActive[YLO] )
    {
      
      /* YLO face, HZ - eqn. (5.54a). */
      j = item->flim[YLO][HZ][YLO];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , k , incField )
      #endif
      for ( i = item->flim[YLO][HZ][XLO] ; i <= item->flim[YLO][HZ][XHI] ; i++ ) 
      {
        for ( k = item->flim[YLO][HZ][ZLO] ; k <= item->flim[YLO][HZ][ZHI] ; k++ )
        {
          incField = SCALE_Ex( incidentField( EX , i , j + 1 , k , timeH , item ) , i );       
          Hz[i][j][k] = Hz[i][j][k] - GAMMA_HZ(i,j,k) * dEx_dy( incField , j );
        }
      }

      /* YLO face, HX - eqn. (5.54b). */
      j = item->flim[YLO][HX][YLO];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , k , incField )
      #endif
      for ( i = item->flim[YLO][HX][XLO] ; i <= item->flim[YLO][HX][XHI] ; i++ ) 
      {
        for ( k = item->flim[YLO][HX][ZLO] ; k <= item->flim[YLO][HX][ZHI] ; k++ )
        {
          incField = SCALE_Ez( incidentField( EZ , i , j + 1 , k , timeH , item ) , k );
          Hx[i][j][k] = Hx[i][j][k] + GAMMA_HX(i,j,k) * dEz_dy( incField , j ); 
        }
      }

    } // if YLO.
 
    if( item->isActive[YHI] )
    {
      
      /* YHI face, HZ - eqn. (5.55a). */
      j = item->flim[YHI][HZ][YHI];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , k , incField )
      #endif
      for ( i = item->flim[YHI][HZ][XLO] ; i <= item->flim[YHI][HZ][XHI] ; i++ ) 
      {
        for ( k = item->flim[YHI][HZ][ZLO] ; k <= item->flim[YHI][HZ][ZHI] ; k++ )
        {
          incField = SCALE_Ex( incidentField( EX , i , j , k , timeH , item ) , i );
          Hz[i][j][k] = Hz[i][j][k] + GAMMA_HZ(i,j,k) * dEx_dy( incField , j );
        }
      }

      /* YHI face, HX - eqn. (5.55b). */
      j = item->flim[YHI][HX][YHI];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , k , incField )
      #endif
      for ( i = item->flim[YHI][HX][XLO] ; i <= item->flim[YHI][HX][XHI] ; i++ ) 
      {
        for ( k = item->flim[YHI][HX][ZLO] ; k <= item->flim[YHI][HX][ZHI] ; k++ )
        {
          incField = SCALE_Ez( incidentField( EZ , i , j , k , timeH , item ) , k );
          Hx[i][j][k] = Hx[i][j][k] - GAMMA_HX(i,j,k) * dEz_dy( incField , j ); 
        }
      }

    } // if YHI.
    
    if( item->isActive[ZLO] )
    {
      
      /* ZLO face, HY - eqn. (5.56a). */
      k = item->flim[ZLO][HY][ZLO];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , j , incField )
      #endif
      for ( i = item->flim[ZLO][HY][XLO] ; i <= item->flim[ZLO][HY][XHI] ; i++ ) 
      {
        for ( j = item->flim[ZLO][HY][YLO] ; j <= item->flim[ZLO][HY][YHI] ; j++ )
        {
          incField = SCALE_Ex( incidentField( EX , i , j , k + 1 , timeH , item ) , i );
          Hy[i][j][k] = Hy[i][j][k] + GAMMA_HY(i,j,k) * dEx_dz( incField , k );
        }
      }

      /* ZLO face, HX - eqn. (5.56b). */
      k = item->flim[ZLO][HX][ZLO];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , j , incField )
      #endif
      for ( i = item->flim[ZLO][HX][XLO] ; i <= item->flim[ZLO][HX][XHI] ; i++ ) 
      {
        for ( j = item->flim[ZLO][HX][YLO] ; j <= item->flim[ZLO][HX][YHI] ; j++ )
        {
          incField = SCALE_Ey( incidentField( EY , i , j , k + 1 , timeH , item ) , j );
          Hx[i][j][k] = Hx[i][j][k] - GAMMA_HX(i,j,k) * dEy_dz( incField , k );
        }
      }

    } // if ZLO.
    
    if( item->isActive[ZHI] )
    {
      
      /* ZHI face, HY - eqn. (5.57a). */
      k = item->flim[ZHI][HY][ZHI];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , j , incField )
      #endif
      for ( i = item->flim[ZHI][HY][XLO] ; i <= item->flim[ZHI][HY][XHI] ; i++ ) 
      {
        for ( j = item->flim[ZHI][HY][YLO] ; j <= item->flim[ZHI][HY][YHI] ; j++ )
        {
          incField = SCALE_Ex( incidentField( EX , i , j , k , timeH , item ) , i );
          Hy[i][j][k] = Hy[i][j][k] - GAMMA_HY(i,j,k) * dEx_dz( incField , k );
        }
      }

      /* ZHI face, HX - eqn. (5.57b). */
      k = item->flim[ZHI][HX][ZHI];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , j , incField )
      #endif
      for ( i = item->flim[ZHI][HX][XLO] ; i <= item->flim[ZHI][HX][XHI] ; i++ ) 
      {
        for ( j = item->flim[ZHI][HX][YLO] ; j <= item->flim[ZHI][HX][YHI] ; j++ )
        {
          incField = SCALE_Ey( incidentField( EY , i , j , k , timeH , item ) , j );
          Hx[i][j][k] = Hx[i][j][k] + GAMMA_HX(i,j,k) * dEy_dz( incField , k );
        }
      }

    } // if ZHI.
    
    if( item->isActive[XLO] )
    {
      
      /* XLO face, HZ - eqn. (5.58a). */
      i = item->flim[XLO][HZ][XLO];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( j , k , incField )
      #endif
      for ( j = item->flim[XLO][HZ][YLO] ; j <= item->flim[XLO][HZ][YHI] ; j++ )
      {
        for ( k = item->flim[XLO][HZ][ZLO] ; k <= item->flim[XLO][HZ][ZHI] ; k++ )
        {
          incField = SCALE_Ey( incidentField( EY , i + 1 , j , k , timeH , item ) , j );
          Hz[i][j][k] = Hz[i][j][k] + GAMMA_HZ(i,j,k) * dEy_dx( incField , i );
        }
      }

      /* XLO face, HY - eqn. (5.58b). */
      i = item->flim[XLO][HY][XLO];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( j , k , incField )
      #endif
      for ( j = item->flim[XLO][HY][YLO] ; j <= item->flim[XLO][HY][YHI] ; j++ )
      {
        for ( k = item->flim[XLO][HY][ZLO] ; k <= item->flim[XLO][HY][ZHI] ; k++ )
        {
          incField = SCALE_Ez( incidentField( EZ , i + 1 , j , k , timeH , item ) , k );
          Hy[i][j][k] = Hy[i][j][k] - GAMMA_HY(i,j,k) * dEz_dx( incField , i );
        }
      }

    } // if XLO.
    
    if( item->isActive[XHI] )
    {
      
      /* XHI face, HZ - eqn. (5.59a). */
      i = item->flim[XHI][HZ][XHI];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( j , k , incField )
      #endif
      for ( j = item->flim[XHI][HZ][YLO] ; j <= item->flim[XHI][HZ][YHI] ; j++ )
      {
        for ( k = item->flim[XHI][HZ][ZLO] ; k <= item->flim[XHI][HZ][ZHI] ; k++ )
        {
          incField = SCALE_Ey( incidentField( EY , i , j , k , timeH , item ) , j );
          Hz[i][j][k] = Hz[i][j][k] - GAMMA_HZ(i,j,k) * dEy_dx( incField , i );
        }
      }

      /* XHI face, HY - eqn. (5.59b). */
      i = item->flim[XHI][HY][XHI];
      #ifdef WITH_OPENMP
      #pragma omp parallel for private( j , k , incField )
      #endif
      for ( j = item->flim[XHI][HY][YLO] ; j <= item->flim[XHI][HY][YHI] ; j++ )
      {
        for ( k = item->flim[XHI][HY][ZLO] ; k <= item->flim[XHI][HY][ZHI] ; k++ )
        {
          incField = SCALE_Ez( incidentField( EZ , i , j , k , timeH , item ) , k );
          Hy[i][j][k] = Hy[i][j][k] + GAMMA_HY(i,j,k) * dEz_dx( incField , i );
        }
      }

    } // if XHI.
    
  }

  return;

}

/* Return true if there are plane waves. */
bool thereArePlaneWaves( void )
{
  
  if( numPlaneWave > 0 )
    return true;
  else
    return false;

}

/* Report plane waves. */
void reportPlaneWaves( void )
{
  
  PlaneWaveItem *item;

  message( MSG_LOG , 0 , "  Number of plane waves: %lu\n" , numPlaneWave );

  DL_FOREACH( planeWaveList , item ) 
  {
    message( MSG_DEBUG3 , 0 , "    Plane wave \"%s\" (#%lu): Waveform#=%d Direction=(%.0f,%.0f) Pol. =%.0f BBOX=[%d,%d,%d,%d,%d,%d] mask=[%d,%d,%d,%d,%d,%d] size=%e delay=%e\n" , 
             item->name , (unsigned long) item->number , item->waveformNumber , 
             item->theta , item->phi , item->eta ,
             item->mbbox[XLO] , item->mbbox[XHI] , item->mbbox[YLO] , 
             item->mbbox[YHI] , item->mbbox[ZLO] , item->mbbox[ZHI] ,
             item->isActive[XLO] , item->isActive[XHI] , item->isActive[YLO] , item->isActive[YHI] , item->isActive[ZLO] , item->isActive[ZHI] , item->size , item->delay );
  }

  return;

}

/* Deallocate plane waves. */
void deallocPlaneWaves( void )
{

  PlaneWaveItem *item , *tmp;

  message( MSG_DEBUG1 , 0 , "Deallocating plane waves...\n" );

  /* Free plane-wave name hash and the plane-waves. */
  HASH_ITER( hh , planeWaveHash , item , tmp )
  {
     if( useAuxGrid )
       deallocAuxGrid( item );
      
    HASH_DELETE( hh , planeWaveHash , item );
    free( item );
  }

  return;

}

/* Get plane-wave number from name. */
bool isPlaneWave( char *name , PlaneWaveIndex *number )
{

  PlaneWaveItem *item;

  HASH_FIND_STR( planeWaveHash , name , item );
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

/* Output gnuplot compatiable plot data for plane waves. */
void gnuplotPlaneWaves( void )
{

  char planeWaveFileName[] = "gnuplot-planewave.dat";
  FILE *outputFile; 
  real kinc[3];
  real Finc[6];
  real Einc[3];
  real Hinc[3];
  real ijk0[3];
  real start[3];
  real end[3];
  real magnitude;
  real minSide;
  real scale[3];
      
  PlaneWaveItem *item;

  outputFile = fopen( planeWaveFileName , "w" );
  if( !outputFile )
    message( MSG_ERROR , 0 , "*** Error: Failed to open plane wave output file %s\n" , planeWaveFileName );

  gnuplotProblemSize( outputFile , mbox );

  DL_FOREACH( planeWaveList , item ) 
  {
    /* Bounding box of TF/Sf region. */
    gnuplotBoundingBox( outputFile , item->mbbox );

    /* Length of arrows is 10 % of bounding box diagonal. */
    //scale = 0.1 * sqrt( ( item->mbbox[XHI] - item->mbbox[XLO] ) * ( item->mbbox[XHI] - item->mbbox[XLO] ) +
    //                    ( item->mbbox[YHI] - item->mbbox[YLO] ) * ( item->mbbox[YHI] - item->mbbox[YLO] ) +
    //                    ( item->mbbox[ZHI] - item->mbbox[ZLO] ) * ( item->mbbox[ZHI] - item->mbbox[ZLO] ) );
    
    /* Find smallest side length of plane-wave bounding box. */
    minSide = fmin( fmin ( item->mbbox[XHI] - item->mbbox[XLO] , item->mbbox[YHI] - item->mbbox[YLO] ) , item->mbbox[ZHI] - item->mbbox[ZLO] );
    
    /* Scale arrows according to smallest side. */
    scale[XDIR] = 0.3 * minSide; //( item->mbbox[XHI] - item->mbbox[XLO] );
    scale[YDIR] = 0.3 * minSide; //( item->mbbox[YHI] - item->mbbox[YLO] );
    scale[ZDIR] = 0.3 * minSide; //( item->mbbox[ZHI] - item->mbbox[ZLO] );
    
    /* Unit vectors for kinc and Einc. */
    calcIncidentFieldVectors( kinc , Finc , ijk0 , item->mbbox , 
                              1.0 , item->theta , item->phi , item->eta );
    magnitude = sqrt( Finc[EX] * Finc[EX] + Finc[EY] * Finc[EY] + Finc[EZ] * Finc[EZ] );
    Einc[XDIR] = Finc[EX] / magnitude;
    Einc[YDIR] = Finc[EY] / magnitude;
    Einc[ZDIR] = Finc[EZ] / magnitude;
    magnitude = sqrt( Finc[HX] * Finc[HX] + Finc[HY] * Finc[HY] + Finc[HZ] * Finc[HZ] );
    Hinc[XDIR] = Finc[HX] / magnitude;
    Hinc[YDIR] = Finc[HY] / magnitude;
    Hinc[ZDIR] = Finc[HZ] / magnitude;
    
    /* kinc vector - open arrow. */
    start[XDIR] = ijk0[XDIR];
    start[YDIR] = ijk0[YDIR];
    start[ZDIR] = ijk0[ZDIR];
    end[XDIR] = ijk0[XDIR] + scale[XDIR] * kinc[XDIR];
    end[YDIR] = ijk0[YDIR] + scale[YDIR] * kinc[YDIR];
    end[ZDIR] = ijk0[ZDIR] + scale[ZDIR] * kinc[ZDIR];
    //printf( "*** kinc [%.2f,%.2f,%.2f]->[%.2f,%.2f,%.2f]\n" , start[XDIR],start[YDIR],start[ZDIR],end[XDIR],end[YDIR],end[ZDIR]);
    gnuplotArrow( outputFile , start , end , Einc , 2 );

    /* Einc vector - closed arrow. */
    end[XDIR] = ijk0[XDIR] + scale[XDIR] * Einc[XDIR];
    end[YDIR] = ijk0[YDIR] + scale[YDIR] * Einc[YDIR];
    end[ZDIR] = ijk0[ZDIR] + scale[ZDIR] * Einc[ZDIR];
    //printf( "*** Einc [%.2f,%.2f,%.2f]->[%.2f,%.2f,%.2f]\n" , start[XDIR],start[YDIR],start[ZDIR],end[XDIR],end[YDIR],end[ZDIR]);
    gnuplotArrow( outputFile , start , end , kinc , 1 );

    /* Hinc vector - closed arrow. */
    end[XDIR] = ijk0[XDIR] + scale[XDIR] * Hinc[XDIR];
    end[YDIR] = ijk0[YDIR] + scale[YDIR] * Hinc[YDIR];
    end[ZDIR] = ijk0[ZDIR] + scale[ZDIR] * Hinc[ZDIR];
    //printf( "*** Hinc [%.2f,%.2f,%.2f]->[%.2f,%.2f,%.2f]\n" , start[XDIR],start[YDIR],start[ZDIR],end[XDIR],end[YDIR],end[ZDIR]);
    gnuplotArrow( outputFile , start , end , Einc , 3 );

  }

  fclose (outputFile);

  return;

}

/* Output gmsh compatiable plot data for plane waves. */
void gmshPlaneWaves( void )
{

  int step[3] = { 1 , 1 , 1 };
  int bbox[6];
  PlaneWaveItem *item;
  unsigned long entityNumber;
  char name[GMSH_NAME_LENGTH];
  
  DL_FOREACH( planeWaveList , item ) 
  {
    
    snprintf( name , GMSH_NAME_LENGTH - 1 , "PW_%s" , item->name );

    for( MeshFace face = XLO ; face <= ZHI ; face++ )
      if( item->isActive[face] )
      {
        entityNumber = gmshGetEntityNumber();
        getFaceOfBoundingBox( bbox , item->mbbox , face );   
        gmshAddEntity( entityNumber , BB_SURFACE , name , bbox , step );
      }
  }

  return;

}

/* Decode face mask. */
bool decodeFaceMask( bool isActive[6] , char maskStr[] )
{
  
  bool ok = true;
  
  for( int face = XLO ; face <= ZHI ; face++ )
  {
    switch( maskStr[face] )
    {
      case '0':
        isActive[face] = false;
        break;
      case '1':
        isActive[face] = true;
        break;
      default:
        ok = false;      
        break;
    }
  }
  
  return ok;

}

/* Initialise auxiliary grid. */
void initAuxGrid( PlaneWaveItem *item )
{

  real relPhaseVelocity;
  int i;
  real dt;
  real d[3];
  real depth, sprof;
  
  dt = getGridTimeStep();
  getUniformGridSize( d );
  
  /* vp(0,0) / vp(item->theta,item->phi) */
  relPhaseVelocity = numericalPhaseVelocity( 0.0 , 0.0 ) / numericalPhaseVelocity( pi * item->theta / 180.0 , pi * item->phi / 180.0 );

  message( MSG_DEBUG3 , 0 , "    Relative numerical phase velocity=%g\n" , relPhaseVelocity );

  /* Free space parameters for incident field buffer. */
  item->betaEyi = dt / ( eps0 * d[0] ) / relPhaseVelocity;
  item->gammaHzi = dt / ( mu0 * d[0] ) / relPhaseVelocity;

  /* Grid size - diagonal of bbox plus 10. */
  item->nx = 6 + npml + sqrt( ( item->gbbox[XHI] - item->gbbox[XLO] ) * ( item->gbbox[XHI] - item->gbbox[XLO] ) +
                              ( item->gbbox[YHI] - item->gbbox[YLO] ) * ( item->gbbox[YHI] - item->gbbox[YLO] ) +
                              ( item->gbbox[ZHI] - item->gbbox[ZLO] ) * ( item->gbbox[ZHI] - item->gbbox[ZLO] ) );
  message( MSG_DEBUG3 , 0 , "    Aux. grid length=%d\n" , item->nx );

  /* Determine position of PML boundary. */

  item->xb = item->nx-npml;
    
  /* Allocate arrays. */
  item->Eyi  = (real *) malloc( sizeof( real ) * ( item->nx + 1 ) );
  item->Hzi  = (real *) malloc( sizeof( real ) * ( item->nx + 1 ) );
  item->Pyi  = (real *) malloc( sizeof( real ) * ( npml ) );
  item->PPyi  = (real *) malloc( sizeof( real ) * ( npml ) );
  item->Bzi  = (real *) malloc( sizeof( real ) * ( npml ) );
  item->adx  = (real *) malloc( sizeof( real ) * ( npml ) );
  item->bdx  = (real *) malloc( sizeof( real ) * ( npml ) );
  item->ahx  = (real *) malloc( sizeof( real ) * ( npml ) );
  item->bhx  = (real *) malloc( sizeof( real ) * ( npml ) );

  /* Clear the mesh and initialise to free space. */
  for ( i = 0 ; i <= item->nx ; i++ ) 
  {
      item->Eyi[i] = 0.0;
      item->Hzi[i] = 0.0;
  }

  /* Initialise PML arrays and create loss profiles. */
  for ( i = 0 ; i < npml ; i++)
  {
      item->Pyi[i] = 0.0;
      item->PPyi[i] = 0.0;
      item->Bzi[i] = 0.0;
      depth = abs( i ) / (real)npml;
      sprof = pow( depth , 4.0 ) * 0.8 * (5) / d[0] / eta0;
      item->bdx[i] = 1.0 / ( 1.0 + sprof );
      item->adx[i] = ( 1.0 - sprof ) / ( 1.0 + sprof );

      depth = (abs( i ) + 0.5) / (real)(npml);
      sprof = pow( depth , 4.0 ) * 0.8 * (5) / d[0] / eta0;
      item->bhx[i] = 1.0 / ( 1.0 + sprof );
      item->ahx[i] = ( 1.0 - sprof ) / ( 1.0 + sprof );

  }


  return;

}

/* Deallocate auxiliary grid. */
void deallocAuxGrid( PlaneWaveItem *item )
{

  free( item->Eyi );
  free( item->Hzi );
  free( item->PPyi );
  free( item->Pyi );
  free( item->Bzi );
  free( item->adx );
  free( item->bdx );
  
  return;
  
}

/* Update electric field in auxiliary grid. */
void updateAuxGridEfield( PlaneWaveItem *item , real time )
{

  int i,Lp;
  real waveform;
  real oldPPyi, oldPyi;

  /* Update incident field buffer - eqn. (5.45a). */
  for ( i = 1 ; i < item->xb ; i++ ) 
    item->Eyi[i] = item->Eyi[i] + item->betaEyi * ( item->Hzi[i-1] - item->Hzi[i] );


  /* Update incident field Eyi from Hzi, PPyi and Pyi in PML region. */
  /* PEC at the edges implemented by not updating those E fields. */
  /* Must come after standard update!*/
  for ( i = item->xb ; i < item->nx; i++ ) 
  {

    Lp = i-item->xb;
    oldPPyi = item->PPyi[Lp];
    item->PPyi[Lp] = item->PPyi[Lp] + item->betaEyi * ( item->Hzi[i-1] - item->Hzi[i] );
    oldPyi = item->Pyi[Lp];
    item->Pyi[Lp] = item->Pyi[Lp] + ( item->PPyi[Lp] - oldPPyi );
    item->Eyi[i] = item->adx[Lp] * item->Eyi[i] + item->bdx[Lp] * ( item->Pyi[Lp] - oldPyi );

  }

  /* Get waveform. */
  waveform = getWaveformValue( time , item->waveformNumber , item->delay );

  /* Add Electric field excitation to incident field buffer - eqn. (5.44). */
  /* can put source strength here too. */
  item->Eyi[m0-2] = waveform;
  
}

/* Update magnetic field in auxiliary grid. */
void updateAuxGridHfield( PlaneWaveItem *item , real time )
{

  int i, Lp;
  real oldBzi;

  /* Update incident field buffer - eqn. (5.45b). */
  for ( i = 0 ; i < item->xb ; i++ ) 
    item->Hzi[i] = item->Hzi[i] + item->gammaHzi * ( item->Eyi[i] - item->Eyi[i+1]);

  /* Update incident field Hzi from Ei and Bzi in PML region. */
  /* Must come after standard update!*/
  for ( i = item->xb  ; i < item->nx ; i++ ) 
  {
    Lp = i-item->xb;
    oldBzi = item->Bzi[Lp];
    item->Bzi[Lp] = item->ahx[Lp] * item->Bzi[Lp] + item->gammaHzi * item->bhx[Lp] * ( -item->Eyi[i+1] + item->Eyi[i] );
    item->Hzi[i] = item->Hzi[i] + ( item->Bzi[Lp] - oldBzi );
  }

}

/* Determine incient field component in cell (i,j,k) at time t using auxiliary grid. */
real incidentFieldAuxGrid( FieldComponent field , int i , int j , int k , real time , PlaneWaveItem *item )
{

  real rcomp[3];
  real d;
  int id;
  real dp;
  real value = 0.0;
  
  /* Location of field point in grid units. */
  getFieldIndexLocation( rcomp , field , i , j , k );

  /* Eqn. (5.41). */
  d = item->kinc[XDIR] * ( rcomp[XDIR] - item->ijk0[XDIR] ) + 
      item->kinc[YDIR] * ( rcomp[YDIR] - item->ijk0[YDIR] ) + 
      item->kinc[ZDIR] * ( rcomp[ZDIR] - item->ijk0[ZDIR] );

  switch( field )
  {
  case EX:
  case EY:
  case EZ:
    /* Eqn. (5.46a). */
    id = floor( d );
    dp = d - id;
    value = ( 1 - dp ) * item->Eyi[m0+id] + dp * item->Eyi[m0+id+1]; 
    break;
  case HX:
  case HY:
  case HZ:
    /* Eqn. (5.46b). */
    id = floor( d + 0.5 );
    dp = d + 0.5 - id;
    /* Insert eta here as we have already included 1/eta0 in Finc! */
    value = eta0 * ( ( 1 - dp ) * item->Hzi[m0-1+id] + dp * item->Hzi[m0+id] );  
    break;
  default:
    assert( 0 );
    break;
  }
  
  /* Eqn. (5.63). */
  return item->Finc[field] * value;

}
