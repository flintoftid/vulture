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

#ifndef _GRID_H_
#define _GRID_H_

#include "vulture.h"
#include "fdtd_types.h"
#include "medium.h"

/* Cells used for PMC on XLO/YLO/ZLO and PMC and tangential fields on XHI/YHI/ZHI. */
#define NUM_GHOST_CELLS 1

/* 
 * Grid types.
 * 
 * Grid types must begin at zero, be contigous and end with GT_UNDEFINED, which
 * *is not* included in the number NUM_GRID_TYPES.
 */

#define NUM_GRID_TYPES 3

/* Grid types. */
typedef enum {

  /* The order and values are significant! */
  GT_CUBIC,
  GT_UNIFORM,
  GT_NONUNIFORM,
  GT_UNDEFINED

} GridType;

/*
 * Global variables.
 */

/* Number of cells in each direction, including ghost cells. */
extern int numCells[3];

/* Mesh bounding box. */
extern int mbox[6];

/* Inner grid (mesh) bounding box. */
extern int gibox[6];
 
/* Outer grid boundaring box. */
extern int gobox[6];

/* Field array limits for inner grid. */
extern int gfilim[6][6];

/* Field array limits for outer grid. */
extern int gfolim[6][6];

/* EM field arrays. */
extern real ***Ex;
extern real ***Ey;
extern real ***Ez;
extern real ***Hx;
extern real ***Hy;
extern real ***Hz;
/* Primary grid edge lengths. */
extern real *dex;
extern real *dey;
extern real *dez;

/* Secondary grid edge lengths. */
extern real *dhx;
extern real *dhy;
extern real *dhz;

/* Inverse primary grid edge lengths. */
extern real *idex;
extern real *idey;
extern real *idez;

/* Inverse secondary grid edge lengths. */
extern real *idhx;
extern real *idhy;
extern real *idhz;

/* 
 * Update coefficient array macros to support indexed and unindexed media. 
 */

#ifdef USE_INDEXED_MEDIA

  extern MediumIndex ***mediumEx;
  extern MediumIndex ***mediumEy;
  extern MediumIndex ***mediumEz;
  extern MediumIndex ***mediumHx;
  extern MediumIndex ***mediumHy;
  extern MediumIndex ***mediumHz;
  
  #define ALPHA_EX(i,j,k) mediumArray[mediumEx[i][j][k]]->alpha
  #define ALPHA_EY(i,j,k) mediumArray[mediumEy[i][j][k]]->alpha  
  #define ALPHA_EZ(i,j,k) mediumArray[mediumEz[i][j][k]]->alpha
  #define BETA_EX(i,j,k)  mediumArray[mediumEx[i][j][k]]->beta
  #define BETA_EY(i,j,k)  mediumArray[mediumEy[i][j][k]]->beta
  #define BETA_EZ(i,j,k)  mediumArray[mediumEz[i][j][k]]->beta
  #define GAMMA_HX(i,j,k) mediumArray[mediumHx[i][j][k]]->gamma
  #define GAMMA_HY(i,j,k) mediumArray[mediumHy[i][j][k]]->gamma
  #define GAMMA_HZ(i,j,k) mediumArray[mediumHz[i][j][k]]->gamma
  
  #define COPY_ALPHA_EX(i1,j1,k1,i0,j0,k0) mediumEx[i1][j1][k1] = mediumEx[i0][j0][k0]
  #define COPY_ALPHA_EY(i1,j1,k1,i0,j0,k0) mediumEy[i1][j1][k1] = mediumEy[i0][j0][k0]
  #define COPY_ALPHA_EZ(i1,j1,k1,i0,j0,k0) mediumEz[i1][j1][k1] = mediumEz[i0][j0][k0]
  #define COPY_BETA_EX(i1,j1,k1,i0,j0,k0)  mediumEx[i1][j1][k1] = mediumEx[i0][j0][k0]
  #define COPY_BETA_EY(i1,j1,k1,i0,j0,k0)  mediumEy[i1][j1][k1] = mediumEy[i0][j0][k0]
  #define COPY_BETA_EZ(i1,j1,k1,i0,j0,k0)  mediumEz[i1][j1][k1] = mediumEz[i0][j0][k0]
  #define COPY_GAMMA_HX(i1,j1,k1,i0,j0,k0) mediumHx[i1][j1][k1] = mediumHx[i0][j0][k0]
  #define COPY_GAMMA_HY(i1,j1,k1,i0,j0,k0) mediumHy[i1][j1][k1] = mediumHy[i0][j0][k0]
  #define COPY_GAMMA_HZ(i1,j1,k1,i0,j0,k0) mediumHz[i1][j1][k1] = mediumHz[i0][j0][k0]
  
#else // USE_INDEXED_MEDIA

  extern real ***alphaEx;
  extern real ***alphaEy;
  extern real ***alphaEz;
  extern real ***betaEx;
  extern real ***betaEy;
  extern real ***betaEz;
  extern real ***gammaHx;
  extern real ***gammaHy;
  extern real ***gammaHz;
  
  #define ALPHA_EX(i,j,k) alphaEx[i][j][k]  
  #define ALPHA_EY(i,j,k) alphaEy[i][j][k] 
  #define ALPHA_EZ(i,j,k) alphaEz[i][j][k] 
  #define BETA_EX(i,j,k)  betaEx[i][j][k] 
  #define BETA_EY(i,j,k)  betaEy[i][j][k] 
  #define BETA_EZ(i,j,k)  betaEz[i][j][k] 
  #define GAMMA_HX(i,j,k) gammaHx[i][j][k] 
  #define GAMMA_HY(i,j,k) gammaHy[i][j][k] 
  #define GAMMA_HZ(i,j,k) gammaHz[i][j][k] 

  #define COPY_ALPHA_EX(i1,j1,k1,i0,j0,k0) alphaEx[i1][j1][k1] = alphaEx[i0][j0][k0]
  #define COPY_ALPHA_EY(i1,j1,k1,i0,j0,k0) alphaEy[i1][j1][k1] = alphaEy[i0][j0][k0]
  #define COPY_ALPHA_EZ(i1,j1,k1,i0,j0,k0) alphaEz[i1][j1][k1] = alphaEz[i0][j0][k0]
  #define COPY_BETA_EX(i1,j1,k1,i0,j0,k0)  betaEx[i1][j1][k1]  = betaEx[i0][j0][k0]
  #define COPY_BETA_EY(i1,j1,k1,i0,j0,k0)  betaEy[i1][j1][k1]  = betaEy[i0][j0][k0]
  #define COPY_BETA_EZ(i1,j1,k1,i0,j0,k0)  betaEz[i1][j1][k1]  = betaEz[i0][j0][k0]
  #define COPY_GAMMA_HX(i1,j1,k1,i0,j0,k0) gammaHx[i1][j1][k1] = gammaHx[i0][j0][k0]
  #define COPY_GAMMA_HY(i1,j1,k1,i0,j0,k0) gammaHy[i1][j1][k1] = gammaHy[i0][j0][k0]
  #define COPY_GAMMA_HZ(i1,j1,k1,i0,j0,k0) gammaHz[i1][j1][k1] = gammaHz[i0][j0][k0]

#endif // USE_INDEXED_MEDIA

/* 
 * Limit checking mode macros.
 *
 * Fields are initialised to special value and propagated for one
 * time step. Before updating each field a check is made that it hasn't
 * already been updated. Each field is marked as visited after update using
 * another special value. After the time iteration a check is performed to 
 * make sure all fields have been updated.
 */

#ifdef CHECK_LIMITS
  #define INITIAL_FIELD_VALUE -1.0
  #define VISITED_FIELD_VALUE 0.0
  #define CHECK_NOT_VISITED( arg ) assert( (arg ) == INITIAL_FIELD_VALUE )
  #define MARK_AS_VISITED( arg ) ( (arg) = VISITED_FIELD_VALUE )
#else
  #define INITIAL_FIELD_VALUE 0.0
  #define VISITED_FIELD_VALUE 0.0
  #define CHECK_NOT_VISITED( arg ) 
  #define MARK_AS_VISITED( arg )  
#endif

/*
 * Macros for scaling and unscaling fields, curl operators and update coefficients.
 * Implemented as macros so there is no performance penalty.
 * Scaled fields are usually preferred as the basic update has one less multiplication.
 */

#ifdef USE_SCALE_FIELDS

  /* Curl operators - without inverse edge lengths. */
  #define curl_Hx( Hz_ijk , Hz_ij1k , Hy_ijk1 , Hy_ijk , i , j , k ) ( Hz_ijk  - Hz_ij1k + Hy_ijk1 - Hy_ijk  )
  #define curl_Hy( Hx_ijk , Hx_ijk1 , Hz_i1jk , Hz_ijk , i , j , k ) ( Hx_ijk  - Hx_ijk1 + Hz_i1jk - Hz_ijk  )
  #define curl_Hz( Hy_ijk , Hy_i1jk , Hx_ij1k , Hx_ijk , i , j , k ) ( Hy_ijk  - Hy_i1jk + Hx_ij1k - Hx_ijk  )
  #define curl_Ex( Ey_ijk1 , Ey_ijk , Ez_ijk , Ez_ij1k , i , j , k ) ( Ey_ijk1 - Ey_ijk  + Ez_ijk  - Ez_ij1k )
  #define curl_Ey( Ez_i1jk , Ez_ijk , Ex_ijk , Ex_ijk1 , i , j , k ) ( Ez_i1jk - Ez_ijk  + Ex_ijk  - Ex_ijk1 )
  #define curl_Ez( Ex_ij1k , Ex_ijk , Ey_ijk , Ey_i1jk , i , j , k ) ( Ex_ij1k - Ex_ijk  + Ey_ijk  - Ey_i1jk )
  
  /* Derivative operators - without inverse edge lengths. */
  #define dHz_dy( Hz_ijk , j ) ( Hz_ijk )
  #define dHy_dz( Hy_ijk , k ) ( Hy_ijk )
  #define dHx_dz( Hx_ijk , k ) ( Hx_ijk )
  #define dHz_dx( Hz_ijk , i ) ( Hz_ijk )
  #define dHy_dx( Hy_ijk , i ) ( Hy_ijk )
  #define dHx_dy( Hx_ijk , j ) ( Hx_ijk )
  #define dEy_dz( Ey_ijk , k ) ( Ey_ijk )
  #define dEz_dy( Ez_ijk , j ) ( Ez_ijk )
  #define dEz_dx( Ez_ijk , i ) ( Ez_ijk )
  #define dEx_dz( Ex_ijk , k ) ( Ex_ijk )
  #define dEx_dy( Ex_ijk , j ) ( Ex_ijk )
  #define dEy_dx( Ey_ijk , i ) ( Ey_ijk )

  /* Scale fields. */
  #define SCALE_Ex( ex , i ) ( dex[i] * ex )
  #define SCALE_Ey( ey , j ) ( dey[j] * ey ) 
  #define SCALE_Ez( ez , k ) ( dez[k] * ez )
  #define SCALE_Hx( hx , i ) ( dhx[i] * hx )
  #define SCALE_Hy( hy , j ) ( dhy[j] * hy )
  #define SCALE_Hz( hz , k ) ( dhz[k] * hz )
  #define SCALE_Jx( jx , i ) ( dhy[j] * dhz[k] * jx )
  #define SCALE_Jy( jy , j ) ( dhx[i] * dhz[k] * jy ) 
  #define SCALE_Jz( jz , k ) ( dhx[i] * dhy[j] * jz )
  #define SCALE_JMx( jmx , i ) ( dey[j] * dez[k] * jmx )
  #define SCALE_JMy( jmy , j ) ( dex[i] * dez[k] * jmy )
  #define SCALE_JMz( jmz , k ) ( dex[i] * dey[j] * jmz )

  /* Unscale fields. */
  #define UNSCALE_Ex( ex , i ) ( idex[i] * ex )
  #define UNSCALE_Ey( ey , j ) ( idey[j] * ey )
  #define UNSCALE_Ez( ez , k ) ( idez[k] * ez )
  #define UNSCALE_Hx( hx , i ) ( idhx[i] * hx )
  #define UNSCALE_Hy( hy , j ) ( idhy[j] * hy )
  #define UNSCALE_Hz( hz , k ) ( idhz[k] * hz )

  /* Scale update ceofficients. */
  #define SCALE_betaEx( betaEx , i , j , k ) (  dex[i] * idhy[j] * idhz[k] * betaEx )
  #define SCALE_betaEy( betaEy , i , j , k ) ( idhx[i] *  dey[j] * idhz[k] * betaEy )
  #define SCALE_betaEz( betaEz , i , j , k ) ( idhx[i] * idhy[j] *  dez[k] * betaEz )
  #define SCALE_gammaHx( gammaHx , i , j , k ) (  dhx[i] * idey[j] * idez[k] * gammaHx )
  #define SCALE_gammaHy( gammaHy , i , j , k ) ( idex[i] *  dhy[j] * idez[k] * gammaHy )
  #define SCALE_gammaHz( gammaHz , i , j , k ) ( idex[i] * idey[j] *  dhz[k] * gammaHz )

  /* Unscale update ceofficients. */
  #define UNSCALE_betaEx( betaEx , i , j , k ) ( idex[i] *  dhy[j] *  dhz[k] * betaEx )
  #define UNSCALE_betaEy( betaEy , i , j , k ) (  dhx[i] * idey[j] *  dhz[k] * betaEy )
  #define UNSCALE_betaEz( betaEz , i , j , k ) (  dhx[i] *  dhy[j] * idez[k] * betaEz )
  #define UNSCALE_gammaHx( gammaHx , i , j , k ) ( idhx[i] *  dey[j] *  dez[k] * gammaHx )
  #define UNSCALE_gammaHy( gammaHy , i , j , k ) (  dex[i] * idhy[j] *  dez[k] * gammaHy )
  #define UNSCALE_gammaHz( gammaHz , i , j , k ) (  dex[i] *  dey[j] *  idhz[k] * gammaHz )

#else // USE_SCALE_FIELDS

  /* Curl operators - includes inverse edge lengths. */  
  #define curl_Hx( Hz_ijk , Hz_ij1k , Hy_ijk1 , Hy_ijk , i , j , k ) ( idhy[j] * ( Hz_ijk - Hz_ij1k ) + idhz[k] * ( Hy_ijk1 - Hy_ijk ) )
  #define curl_Hy( Hx_ijk , Hx_ijk1 , Hz_i1jk , Hz_ijk , i , j , k ) ( idhz[k] * ( Hx_ijk - Hx_ijk1 ) + idhx[i] * ( Hz_i1jk - Hz_ijk ) )
  #define curl_Hz( Hy_ijk , Hy_i1jk , Hx_ij1k , Hx_ijk , i , j , k ) ( idhx[i] * ( Hy_ijk - Hy_i1jk ) + idhy[j] * ( Hx_ij1k - Hx_ijk ) )
  #define curl_Ex( Ey_ijk1 , Ey_ijk , Ez_ijk , Ez_ij1k , i , j , k ) ( idez[k] * ( Ey_ijk1 - Ey_ijk ) + idey[j] * ( Ez_ijk - Ez_ij1k ) )
  #define curl_Ey( Ez_i1jk , Ez_ijk , Ex_ijk , Ex_ijk1 , i , j , k ) ( idex[i] * ( Ez_i1jk - Ez_ijk ) + idez[k] * ( Ex_ijk - Ex_ijk1 ) )
  #define curl_Ez( Ex_ij1k , Ex_ijk , Ey_ijk , Ey_i1jk , i , j , k ) ( idey[j] * ( Ex_ij1k - Ex_ijk ) + idex[i] * ( Ey_ijk - Ey_i1jk ) )

  /* Derivative operators - includes inverse edge lengths. */  
  #define dHz_dy( Hz_ijk , j ) ( idhy[j] * Hz_ijk )
  #define dHy_dz( Hy_ijk , k ) ( idhz[k] * Hy_ijk )
  #define dHx_dz( Hx_ijk , k ) ( idhz[k] * Hx_ijk )
  #define dHz_dx( Hz_ijk , i ) ( idhx[i] * Hz_ijk )
  #define dHy_dx( Hy_ijk , i ) ( idhx[i] * Hy_ijk )
  #define dHx_dy( Hx_ijk , j ) ( idhy[j] * Hx_ijk )
  #define dEy_dz( Ey_ijk , k ) ( idez[k] * Ey_ijk )
  #define dEz_dy( Ez_ijk , j ) ( idey[j] * Ez_ijk )
  #define dEz_dx( Ez_ijk , i ) ( idex[i] * Ez_ijk )
  #define dEx_dz( Ex_ijk , k ) ( idez[k] * Ex_ijk )
  #define dEx_dy( Ex_ijk , j ) ( idey[j] * Ex_ijk )
  #define dEy_dx( Ey_ijk , i ) ( idex[i] * Ey_ijk )

  /* Scale fields - no-op. */
  #define SCALE_Ex( ex , i ) ( ex )
  #define SCALE_Ey( ey , j ) ( ey )
  #define SCALE_Ez( ez , k ) ( ez )
  #define SCALE_Hx( hx , i ) ( hx )
  #define SCALE_Hy( hy , j ) ( hy )
  #define SCALE_Hz( hz , k ) ( hz )
  #define SCALE_Jx( jx , i ) ( jx )
  #define SCALE_Jy( jy , j ) ( jy )
  #define SCALE_Jz( jz , k ) ( jz )
  #define SCALE_JMx( jmx , i ) ( jmx )
  #define SCALE_JMy( jmy , j ) ( jmy )
  #define SCALE_JMz( jmz , k ) ( jmz )

  /* Unscale fields - no-op. */
  #define UNSCALE_Ex( ex , i ) ( ex )
  #define UNSCALE_Ey( ey , j ) ( ey )
  #define UNSCALE_Ez( ez , k ) ( ez )
  #define UNSCALE_Hx( hx , i ) ( hx )
  #define UNSCALE_Hy( hy , j ) ( hy )
  #define UNSCALE_Hz( hz , k ) ( hz )

  /* Scale update ceofficients - no-op. */
  #define SCALE_betaEx( betaEx , i , j , k ) ( betaEx )
  #define SCALE_betaEy( betaEy , i , j , k ) ( betaEy )
  #define SCALE_betaEz( betaEz , i , j , k ) ( betaEz )
  #define SCALE_gammaHx( gammaHx , i , j , k ) ( gammaHx )
  #define SCALE_gammaHy( gammaHy , i , j , k ) ( gammaHy )
  #define SCALE_gammaHz( gammaHz , i , j , k ) ( gammaHz )

  /* Unscale update ceofficients - no-op. */
  #define UNSCALE_betaEx( betaEx , i , j , k ) ( betaEx )
  #define UNSCALE_betaEy( betaEy , i , j , k ) ( betaEy )
  #define UNSCALE_betaEz( betaEz , i , j , k ) ( betaEz )
  #define UNSCALE_gammaHx( gammaHx , i , j , k ) ( gammaHx )
  #define UNSCALE_gammaHy( gammaHy , i , j , k ) ( gammaHy )
  #define UNSCALE_gammaHz( gammaHz , i , j , k ) ( gammaHz )

#endif // USE_SCALE_FIELDS

/*
 * Public method interfaces.
 */
 
bool parseDM( char *line );
bool parseMS( char *line );
bool parseXL( char *line );
bool parseYL( char *line );
bool parseZL( char *line );
void initGrid( void );
void reportGrid( void );
void updateGridEfield( void );
void updateGridHfield( void  );
void deallocGridArrays( void );
void gnuplotGridLines( void );
void checkGrid( void );
void setFieldLimits( int cellLimits[6] , int fieldLimits[6][6] , bool includeBoundary[6] );
real getGridMaxEdgeLength( CoordAxis direction );
real getGridTimeStep( void );
void getGridNumCells( int numMeshCells[3] );
void setMediumOnGrid( int bbox[6] , MediumIndex medium , FaceMask mask ); 
void initMediaArrays( void  );
void dumpMediaOnGrid( FieldComponent field );
void bboxInPhysicalUnits( real physbbox[6] , int bbox[6] );
real indexInPhysicalUnits( int index , CoordAxis dir );
real getMeshLineCoordX( int lineNumber );
real getMeshLineCoordY( int lineNumber );
real getMeshLineCoordZ( int lineNumber );
void getMeshNodeCoords( real nodeCoords[3] , int nodeIndices[3] );
void bboxInPhysicalUnits( real physbbox[6] , int bbox[6] );
void getGridBoundingBox( int innerBox[6] , int outerBox[6] );
void getFieldPhysicalLocation( real r[3] , FieldComponent field , int i , int j , int k );
void getFieldIndexLocation( real ijk[3] , FieldComponent field , int i , int j , int k );
void getNodeLocation( real r[3] , int i , int j , int k );
GridType getGridType();
void getUniformGridSize( real d[3] );
real numericalPhaseVelocity( real theta , real phi );
void nodeInPhysicalUnits( real r[3] , real ijk[3] );
void checkMediumOnGrid( int gbbox[6] , MediumIndex medium );
void applyVoxelsToGrid( MediumIndex ***blockArray );

#endif
