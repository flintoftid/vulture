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
#include <stdbool.h>

#include "pml.h"
#include "grid.h"
#include "alloc_array.h"
#include "message.h"
#include "surface.h"
#include "medium.h"
#include "memory.h"
#include "physical.h"

/* 
 * Private data.
 */

static int pbox[6][6];       // Bounding boxes for PML regions.
static int fplim[6][6][6];   // Field array limits for PML regions.
static real ***PPx[6];       // Auxiliary PML field arrays.
static real ***PPy[6];
static real ***PPz[6];
static real ***Px[6];
static real ***Py[6];
static real ***Pz[6];
static real ***Bx[6];
static real ***By[6];
static real ***Bz[6];
static real *adx;            // PML loss profile arrays.
static real *bdx;
static real *gdx;
static real *kdx;
static real *ady;
static real *bdy;
static real *gdy;
static real *kdy;
static real *adz;
static real *bdz;
static real *gdz;
static real *kdz;
static real *ahx;
static real *bhx;
static real *ghx;
static real *khx;
static real *ahy;
static real *bhy;
static real *ghy;
static real *khy;
static real *ahz;
static real *bhz;
static real *ghz;
static real *khz;
static real *ibdx;
static real *ibdy;
static real *ibdz;
static real *ibhx;
static real *ibhy;
static real *ibhz;

/* 
 * Private method interfaces. 
 */

void setPmlLimits( void );
void setPmlParameters( void );
void allocPmlArrays( void );
void initPmlMaterialArrays( void );
void initPmlMaterialArrays2( void );
void clearPml( void );
void setProfile( real *de , real *dh , int low , int high , int dir , 
                 int order , real n_eff , real refCoeff , real kmax ,
                 real *gd , real *kd , real *bd , real *ad ,
                 real *gh , real *kh , real *bh , real *ah ,
                 real *ibd , real *ibh );
real pmlSigmaProfile( real x , real totalDepth , real meshSize , int order , real n_eff , real refCoeff );
real pmlKappaProfile( real x , real totalDepth , real meshSize , int order , real kmax );
/*
 * Method Implementations.
 */

/* Initialise PML boundaries. */
void initPmlBoundaries( void )
{

  message( MSG_LOG , 0 , "\nInitialising PML...\n\n" );

  /* Set PML region array limits. */
  setPmlLimits();
 
  /* Allocate PML arrays. */
  allocPmlArrays();

  /* Initialise PML parameters to free space. */
  setPmlParameters();

  /* Carry materials on inner boundaries into PML.*/
  initPmlMaterialArrays();

  /* Clear the PML arrays. */
  clearPml();

  return;

}

/* Set field and array limits for PML regions. */
void setPmlLimits( void )
{
                                  // XLO      XHI     YLO     YHI     ZLO     ZHI
  bool includeBoundary[6][6] = { {  true ,  false ,  true ,  true ,  true ,  true } ,  // XLO region.
                                 { false ,   true ,  true ,  true ,  true ,  true } ,  // XHI region.
                                 {  true ,   true ,  true , false ,  true ,  true } ,  // YLO region.
                                 {  true ,   true , false ,  true ,  true ,  true } ,  // YHI region.
                                 {  true ,   true ,  true ,  true ,  true , false } ,  // ZLO region.
                                 {  true ,   true ,  true ,  true , false ,  true } }; // ZHI region.

  message( MSG_LOG , 0 , "  Setting PML limits ...\n" );

  /* Set cell limits for PML regions. */
  /* XLO region. */
  pbox[XLO][XLO] = gobox[XLO];
  pbox[XLO][XHI] = gibox[XLO];
  pbox[XLO][YLO] = gobox[YLO];
  pbox[XLO][YHI] = gobox[YHI];
  pbox[XLO][ZLO] = gobox[ZLO];
  pbox[XLO][ZHI] = gobox[ZHI];
  /* XHI region. */
  pbox[XHI][XLO] = gibox[XHI];
  pbox[XHI][XHI] = gobox[XHI];
  pbox[XHI][YLO] = gobox[YLO];
  pbox[XHI][YHI] = gobox[YHI];
  pbox[XHI][ZLO] = gobox[ZLO];
  pbox[XHI][ZHI] = gobox[ZHI];
  /* YLO region. */
  pbox[YLO][XLO] = gibox[XLO];
  pbox[YLO][XHI] = gibox[XHI];
  pbox[YLO][YLO] = gobox[YLO];
  pbox[YLO][YHI] = gibox[YLO];
  pbox[YLO][ZLO] = gobox[ZLO];
  pbox[YLO][ZHI] = gobox[ZHI];
  /* YHI region. */
  pbox[YHI][XLO] = gibox[XLO];
  pbox[YHI][XHI] = gibox[XHI];
  pbox[YHI][YLO] = gibox[YHI];
  pbox[YHI][YHI] = gobox[YHI];
  pbox[YHI][ZLO] = gobox[ZLO];
  pbox[YHI][ZHI] = gobox[ZHI];
  /* ZLO region. */
  pbox[ZLO][XLO] = gibox[XLO];
  pbox[ZLO][XHI] = gibox[XHI];
  pbox[ZLO][YLO] = gibox[YLO];
  pbox[ZLO][YHI] = gibox[YHI];
  pbox[ZLO][ZLO] = gobox[ZLO];
  pbox[ZLO][ZHI] = gibox[ZLO];
  /* ZHI region. */
  pbox[ZHI][XLO] = gibox[XLO];
  pbox[ZHI][XHI] = gibox[XHI];
  pbox[ZHI][YLO] = gibox[YLO];
  pbox[ZHI][YHI] = gibox[YHI];
  pbox[ZHI][ZLO] = gibox[ZHI];
  pbox[ZHI][ZHI] = gobox[ZHI];

  /* Field limits for PML regions. */
  for( MeshFace region = XLO ; region <= ZHI ; region++ )
    setFieldLimits( pbox[region] , fplim[region] , includeBoundary[region] );

  return;

}

/* Set PML parameter profiles in PML regions. */
void setProfile( real *de , real *dh , int low , int high , int dir , 
                 int order , real n_eff , real refCoeff , real kmax ,
                 real *gd , real *kd , real *bd , real *ad ,
                 real *gh , real *kh , real *bh , real *ah ,
                 real *ibd , real *ibh )
{

  real total_depth = 0.0;
  real depthE , depthH;
  int i;
  real sprof, kprof;
  real dt;
  
  dt = getGridTimeStep();
  
  for ( i = low ; i <= high ; i++ )
    total_depth =  total_depth + de[i];

  if( dir == -1 )
  {
    depthE = total_depth;
    depthH = total_depth - 0.5 * de[low];
  }
  else
  {
    depthE = 0.0;
    depthH = 0.5 * de[low];
  }

  for ( i = low ; i <= high ; i++ )
  {
    sprof = 0.5 * dt / eps0 * pmlSigmaProfile( depthE / total_depth , total_depth , de[low] , order , n_eff , refCoeff );
    kprof = pmlKappaProfile( depthE / total_depth , total_depth , de[low] , order , kmax );
    gh[i] = sprof;
    kh[i] = kprof;
    bd[i] = 1.0 / ( kprof + sprof );
    ibd[i] = kprof + sprof;
    ad[i] = ( kprof - sprof ) / ( kprof + sprof );
    depthE = depthE + dir * de[i];
  }

  if( dir == + 1 )
  {
    sprof = 0.5 * dt / eps0 * pmlSigmaProfile( depthE / total_depth , total_depth , de[low] , order , n_eff , refCoeff );
    kprof = pmlKappaProfile( depthE / total_depth , total_depth , de[low] , order , kmax );
    gh[high+1] = sprof;
    kh[high+1] = kprof;
    bd[high+1] = 1.0 / ( kprof + sprof );
    ibd[high+1] = kprof + sprof;
    ad[high+1] = ( kprof - sprof ) / ( kprof + sprof );  
  }

  for ( i = low ; i <= high ; i++ )
  {
    sprof = 0.5 * dt / eps0 * pmlSigmaProfile( depthH / total_depth , total_depth , de[low] , order , n_eff , refCoeff );
    kprof = pmlKappaProfile( depthH / total_depth , total_depth , de[low] , order , kmax );
    gd[i] = sprof;
    kd[i] = kprof;
    bh[i] = 1.0 / ( kprof + sprof );
    ibh[i] = kprof + sprof;
    ah[i] = ( kprof - sprof) / ( kprof + sprof );
    depthH = depthH + dir * dh[i];
  }

  return;

}

/* Loss profile for PML. x=0 is PML to internal space boundary and */
/* x=1 is the PML to PEC boundary. */
real pmlSigmaProfile( real x , real totalDepth , real meshSize , int order , real n_eff , real refCoeff )
{

  real sigmaMax;

  assert( x >= 0.0 && x <= 1.0 );
  
  if( refCoeff > 0 )
  {
    /* Theoretical profile for given reflection coefficient. Taflove eqn. (7.55a) & (7.57). */
    sigmaMax = -( order + 1.0 ) / ( 2.0  * eta0 * n_eff * totalDepth ) * log( refCoeff );
  }
  else
  {
    /* Optimum profile. */
    /* Taflove eqn. (7.61). |R| = exp(-16) for N = 10 and |R| = exp(-8) for N = 5. */
    sigmaMax = 4.0 * ( order + 1.0 ) / ( 5.0  * eta0 * n_eff * meshSize );
  }

  return sigmaMax * pow( x , order );

}

/* Kappa profile for PML. x=0 is PML to internal space boundary and */
/* x=1 is the PML to PEC boundary. */
real pmlKappaProfile( real x , real totalDepth , real meshSize , int order , real kmax )
{
  
  /* Polynmial profile, Taflove eqn. (7.55b). */
  return 1.0 + ( kmax - 1.0 ) * pow( x , order );

}

/* Set PML grading parameters. */
void setPmlParameters( void )
{

  int order;
  real n_eff;
  real refCoeff;
  real kmax;
  
  for ( int i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ ) 
  {
    kdx[i] = 1.0;
    khx[i] = 1.0;
    gdx[i] = 0.0;
    ghx[i] = 0.0;
    bdx[i] = 1.0;
    bhx[i] = 1.0;
    ibdx[i] = 1.0;
    ibhx[i] = 1.0;
    adx[i] = 1.0;
    ahx[i] = 1.0;
  }

  for ( int j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
  {
    kdy[j] = 1.0;
    khy[j] = 1.0;
    gdy[j] = 0.0;
    ghy[j] = 0.0;
    bdy[j] = 1.0;
    bhy[j] = 1.0;
    ibdy[j] = 1.0;
    ibhy[j] = 1.0;
    ady[j] = 1.0;
    ahy[j] = 1.0;
  }

  for ( int k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
  {
    kdz[k] = 1.0;
    khz[k] = 1.0;
    gdz[k] = 0.0;
    ghz[k] = 0.0;
    bdz[k] = 1.0;
    bhz[k] = 1.0;
    ibdz[k] = 1.0;
    ibhz[k] = 1.0;
    adz[k] = 1.0;
    ahz[k] = 1.0;
  }

  
  /* Set up PML loss profiles for each face. */
  if( outerSurfaceType( XLO ) == BT_PML )
  {
    getOuterSurfaceParams( XLO , &order , &n_eff , &refCoeff , &kmax );
    message( MSG_DEBUG3 , 0 , "  Setting PML profile, XLO: order=%d n_eff=%e refCeoff=%e kmax=%e\n" , 
             order , n_eff , refCoeff , kmax );
    setProfile( dex , dhx  , gobox[XLO] , gibox[XLO] - 1 , -1 , 
                order , n_eff , refCoeff , kmax ,
                gdx , kdx , bdx , adx , ghx , khx , bhx , ahx , ibdx , ibhx );
  }
  
  if( outerSurfaceType( XHI ) == BT_PML )
  {
    getOuterSurfaceParams( XHI , &order , &n_eff , &refCoeff , &kmax );
    message( MSG_DEBUG3 , 0 , "  Setting PML profile, XHI: order=%d n_eff=%e refCeoff=%e kmax=%e\n" , 
             order , n_eff , refCoeff , kmax );
    setProfile( dex , dhx , gibox[XHI] , gobox[XHI] - 1 , +1 , 
                order , n_eff , refCoeff , kmax ,
                gdx , kdx , bdx , adx , ghx , khx , bhx , ahx , ibdx , ibhx );
  }
  
  if( outerSurfaceType( YLO ) == BT_PML )
  {
    getOuterSurfaceParams( YLO , &order , &n_eff , &refCoeff , &kmax );
    message( MSG_DEBUG3 , 0 , "  Setting PML profile, YLO: order=%d n_eff=%e refCeoff=%e kmax=%e\n" , 
             order , n_eff , refCoeff , kmax );
    setProfile( dey , dhy  , gobox[YLO] , gibox[YLO] - 1 , -1 , 
                order , n_eff , refCoeff , kmax ,
                gdy , kdy , bdy , ady , ghy , khy , bhy , ahy , ibdy , ibhy );
  }
  
  if( outerSurfaceType( YHI ) == BT_PML )
  {
    getOuterSurfaceParams( YHI , &order , &n_eff , &refCoeff , &kmax );
    message( MSG_DEBUG3 , 0 , "  Setting PML profile, YHI: order=%d n_eff=%e refCeoff=%e kmax=%e\n" , 
             order , n_eff , refCoeff , kmax );
    setProfile( dey , dhy , gibox[YHI] , gobox[YHI] - 1 , +1 , 
                order , n_eff , refCoeff , kmax ,
                gdy , kdy , bdy , ady , ghy , khy , bhy , ahy , ibdy , ibhy );
  }
  
  if( outerSurfaceType( ZLO ) == BT_PML )
  {
    getOuterSurfaceParams( ZLO , &order , &n_eff , &refCoeff , &kmax );
    message( MSG_DEBUG3 , 0 , "  Setting PML profile, ZLO: order=%d n_eff=%e refCeoff=%e kmax=%e\n" , 
             order , n_eff , refCoeff , kmax );
    setProfile( dez , dhz  , gobox[ZLO] , gibox[ZLO] - 1 , -1 , 
                order , n_eff , refCoeff , kmax ,
                gdz , kdz , bdz , adz , ghz , khz , bhz , ahz , ibdz , ibhz );
  }
  
  if( outerSurfaceType( ZHI ) == BT_PML )
  {
    getOuterSurfaceParams( ZHI , &order , &n_eff , &refCoeff , &kmax );
    message( MSG_DEBUG3 , 0 , "  Setting PML profile, ZHI: order=%d n_eff=%e refCeoff=%e kmax=%e\n" , 
             order , n_eff , refCoeff , kmax );
    setProfile( dez , dhz , gibox[ZHI] , gobox[ZHI] - 1 , +1 , 
                order , n_eff , refCoeff , kmax ,
                gdz , kdz , bdz , adz , ghz , khz , bhz , ahz , ibdz , ibhz );
  }
  
  message( MSG_DEBUG3 , 0 , "  PML profile, XDIR:\n" );
  for ( int i = gobox[XLO]; i <= gobox[XHI] ; i++ ) 
  {
    message( MSG_DEBUG3 , 0 , "    %5.1f %6.4f %6.4f %6.4f %6.4f\n", 
             (real)i, khx[i], adx[i], bdx[i], ghx[i] );
    message( MSG_DEBUG3 , 0  ,"    %5.1f %6.4f %6.4f %6.4f %6.4f\n", 
             (real)i + 0.5, kdx[i], ahx[i], bhx[i], gdx[i] );
  }

  message( MSG_DEBUG3 , 0 , "  PML profile, YDIR:\n" );
  for ( int j = gobox[YLO]; j <= gobox[YHI] ; j++ )
  {
    message( MSG_DEBUG3 , 0 , "    %5.1f %6.4f %6.4f %6.4f %6.4f\n", 
             (real)j, khy[j], ady[j], bdy[j], ghy[j] );
    message( MSG_DEBUG3 , 0 , "    %5.1f %6.4f %6.4f %6.4f %6.4f\n", 
             (real)j + 0.5, kdy[j], ahy[j], bhy[j], gdy[j] );
  }

  message( MSG_DEBUG3 , 0 , "  PML profile, ZDIR:\n");
  for ( int k = gobox[ZLO] ; k <= gobox[ZHI] ; k++ ) 
  {
    message( MSG_DEBUG3 , 0 , "    %5.1f %6.4f %6.4f %6.4f %6.4f\n", 
             (real)k, khz[k], adz[k], bdz[k], ghz[k] );
    message( MSG_DEBUG3 , 0 , "    %5.1f %6.4f %6.4f %6.4f %6.4f\n", 
             (real)k+0.5, kdz[k], ahz[k], bhz[k], gdz[k] );
  }

  return;

}

/* Copy material arrays into PML.
 * Operates on the scaled update coefficients.
 *
 * The order of the faces - from Z to Y to X is significant!
 * The ZLO/ZHI faces are done first to carry the necessary material
 * parameters from the inner mesh faces onto their side faces - these
 * are subsequently used by the YLO/YHI/XLO/XHI faces to carry into
 * their own regions. Similarly for the YLO/YHI faces with regard to
 * XLO/XHI. If a PML is present then it is responsible for updating the
 * tangential E and normal H on its outer faces (i.e. on the PEC.
 * Otherwise there are updated by the main grid updating functions.
 */ 
void initPmlMaterialArrays( void )
{
                                 // XLO      XHI     YLO     YHI     ZLO     ZHI
  //bool includeBoundary[6][6] = { { false ,  false , false , false , false , false } ,  // XLO region.
  //                               { false ,  false , false , false , false , false } ,  // XHI region.
  //                               {  true ,   true , false , false , false , false } ,  // YLO region.
  //                              {  true ,   true , false , false , false , false } ,  // YHI region.
  //                               {  true ,   true ,  true ,  true , false , false } ,  // ZLO region.
  //                               {  true ,   true ,  true ,  true , false , false } }; // ZHI region.
  bool includeBoundary[6][6] = { { false ,  false ,  true ,  true ,  true ,  true } ,  // XLO region.
                                 { false ,  false ,  true ,  true ,  true ,  true } ,  // XHI region.
                                 {  true ,   true , false , false ,  true ,  true } ,  // YLO region.
                                 {  true ,   true , false , false ,  true ,  true } ,  // YHI region.
                                 {  true ,   true ,  true ,  true , false , false } ,  // ZLO region.
                                 {  true ,   true ,  true ,  true , false , false } }; // ZHI region.
                                 
  int fpmlim[6][6][6];

  /* Field limits for PML regions. */
  for( MeshFace region = XLO ; region <= ZHI ; region++ )
    setFieldLimits( pbox[region] , fpmlim[region] , includeBoundary[region] );

  /* Counters. */
  int i , j , k;
  int region , field;
  int offset;

  message( MSG_LOG , 0 , "  Initialising PML materials...\n" );

  for( region = ZLO ; region <= ZHI ; region++ )
  {
    if( region == ZLO )
      offset = +1;
    else
      offset = -1;
    field = EX;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
        {
          COPY_ALPHA_EX( i , j , k , i , j , gibox[region] );
          COPY_BETA_EX( i , j , k , i , j , gibox[region] );
        }
    field = EY;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
        {
          COPY_ALPHA_EY( i , j , k , i , j , gibox[region] );
          COPY_BETA_EY( i , j , k , i , j , gibox[region] );
        }
    field = EZ;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
        {
          COPY_ALPHA_EZ( i , j , k , i , j , gibox[region] + offset );
          COPY_BETA_EZ( i , j , k , i , j , gibox[region] + offset );
        }
    field = HX;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
          COPY_GAMMA_HX( i , j , k , i , j , gibox[region] + offset );
    field = HY;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
          COPY_GAMMA_HY( i , j , k , i , j , gibox[region] + offset );
    field = HZ;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
          COPY_GAMMA_HZ( i , j , k , i , j , gibox[region] );
  }

  for( region = YLO ; region <= YHI ; region++ )
  {
    if( region == YLO )
      offset = +1;
    else
      offset = -1;
    field = EX;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
        {
          COPY_ALPHA_EX( i , j , k , i , gibox[region] , k );
          COPY_BETA_EX( i , j , k , i , gibox[region] , k );
        }
    field = EY;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
        {
          COPY_ALPHA_EY( i , j , k , i , gibox[region] + offset , k );
          COPY_BETA_EY( i , j , k , i , gibox[region] + offset , k );
        }
    field = EZ;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
        {
          COPY_ALPHA_EZ( i , j , k , i , gibox[region] , k );
          COPY_BETA_EZ( i , j , k , i , gibox[region] , k );
        }
    field = HX;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
          COPY_GAMMA_HX( i , j , k , i , gibox[region] + offset , k );
    field = HY;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
          COPY_GAMMA_HY( i , j , k , i , gibox[region] , k );
    field = HZ;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
          COPY_GAMMA_HZ( i , j , k , i , gibox[region] + offset , k );
  }

  for( region = XLO ; region <= XHI ; region++ )
  {
    if( region == XLO )
      offset = +1;
    else
      offset = -1;
    field = EX;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
        {
          COPY_ALPHA_EX( i , j , k , gibox[region] + offset , j , k );
          COPY_BETA_EX( i , j , k , gibox[region] + offset , j , k );
        }
    field = EY;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
        {
          COPY_ALPHA_EY( i , j , k , gibox[region] , j , k );
          COPY_BETA_EY( i , j , k , gibox[region] , j , k );
        }
    field = EZ;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
          {
            COPY_ALPHA_EZ( i , j , k , gibox[region] , j , k );
            COPY_BETA_EZ( i , j , k , gibox[region] , j , k );
          }
    field = HX;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
          COPY_GAMMA_HX( i , j , k , gibox[region] , j , k );
    field = HY;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
          COPY_GAMMA_HY( i , j , k , gibox[region] + offset , j , k );
    field = HZ;
    for ( i = fpmlim[region][field][XLO] ; i <= fpmlim[region][field][XHI] ; i++ )
      for ( j = fpmlim[region][field][YLO] ; j <= fpmlim[region][field][YHI] ; j++ )
        for ( k = fpmlim[region][field][ZLO] ; k <= fpmlim[region][field][ZHI] ; k++ )
          COPY_GAMMA_HZ( i , j , k , gibox[region] + offset , j , k );
  }

  return;

}

/* Allocate PML arrays. */
void allocPmlArrays( void )
{

  int region;
  unsigned long bytes;

  message( MSG_LOG , 0 , "  Allocating PML arrays...\n" );

  for( region = XLO ; region <= ZHI ; region++ )
  {
    if( outerSurfaceType( region ) == BT_PML )
    {
      message( MSG_DEBUG1 , 0 , "  Allocating grid PML Px[%s] array\n" , FACE[region] );
      Px[region] = allocArray( &bytes , sizeof( real ) , 3 , fplim[region][EX][XHI] - fplim[region][EX][XLO] + 1 , 
                                                        fplim[region][EX][YHI] - fplim[region][EX][YLO] + 1 , 
                                                        fplim[region][EX][ZHI] - fplim[region][EX][ZLO] + 1 );
      memory.pmlFields += bytes;
      message( MSG_DEBUG1 , 0 , "  Allocating grid PML Py[%s] array\n" , FACE[region] );
      Py[region] = allocArray( &bytes , sizeof( real ) , 3 , fplim[region][EY][XHI] - fplim[region][EY][XLO] + 1 , 
                                                        fplim[region][EY][YHI] - fplim[region][EY][YLO] + 1 , 
                                                        fplim[region][EY][ZHI] - fplim[region][EY][ZLO] + 1 );
      memory.pmlFields += bytes;
      message( MSG_DEBUG1 , 0 , "  Allocating grid PML Pz[%s] array\n" , FACE[region] );
      Pz[region] = allocArray( &bytes , sizeof( real ) , 3 , fplim[region][EZ][XHI] - fplim[region][EZ][XLO] + 1 , 
                                                        fplim[region][EZ][YHI] - fplim[region][EZ][YLO] + 1 , 
                                                        fplim[region][EZ][ZHI] - fplim[region][EZ][ZLO] + 1 );
      memory.pmlFields += bytes;
      message( MSG_DEBUG1 , 0 , "  Allocating grid PML PPx[%s] array\n" , FACE[region] );
      PPx[region] = allocArray( &bytes , sizeof( real ) , 3 , fplim[region][EX][XHI] - fplim[region][EX][XLO] + 1 , 
                                                         fplim[region][EX][YHI] - fplim[region][EX][YLO] + 1 , 
                                                         fplim[region][EX][ZHI] - fplim[region][EX][ZLO] + 1 );
      memory.pmlFields += bytes;
      message( MSG_DEBUG1 , 0 , "  Allocating grid PML PPy[%s] array\n" , FACE[region] );
      PPy[region] = allocArray( &bytes , sizeof( real ) , 3 , fplim[region][EY][XHI] - fplim[region][EY][XLO] + 1 , 
                                                         fplim[region][EY][YHI] - fplim[region][EY][YLO] + 1 , 
                                                         fplim[region][EY][ZHI] - fplim[region][EY][ZLO] + 1 );
      memory.pmlFields += bytes;
      message( MSG_DEBUG1 , 0 , "  Allocating grid PML PPz[%s] array\n" , FACE[region] );
      PPz[region] = allocArray( &bytes , sizeof( real ) , 3 , fplim[region][EZ][XHI] - fplim[region][EZ][XLO] + 1 , 
                                                         fplim[region][EZ][YHI] - fplim[region][EZ][YLO] + 1 , 
                                                         fplim[region][EZ][ZHI] - fplim[region][EZ][ZLO] + 1 );
      memory.pmlFields += bytes;
      message( MSG_DEBUG1 , 0 , "  Allocating grid PML Bx[%s] array\n" , FACE[region] );
      Bx[region] = allocArray( &bytes , sizeof( real ) , 3 , fplim[region][HX][XHI] - fplim[region][HX][XLO] + 1 , 
                                                        fplim[region][HX][YHI] - fplim[region][HX][YLO] + 1 , 
                                                        fplim[region][HX][ZHI] - fplim[region][HX][ZLO] + 1 );
      memory.pmlFields += bytes;
      message( MSG_DEBUG1 , 0 , "  Allocating grid PML By[%s] array\n" , FACE[region] );
      By[region] = allocArray( &bytes , sizeof( real ) , 3 , fplim[region][HY][XHI] - fplim[region][HY][XLO] + 1 , 
                                                        fplim[region][HY][YHI] - fplim[region][HY][YLO] + 1 , 
                                                        fplim[region][HY][ZHI] - fplim[region][HY][ZLO] + 1 );
      memory.pmlFields += bytes;
      message( MSG_DEBUG1 , 0 , "  Allocating grid PML Bz[%s] array\n" , FACE[region] );
      Bz[region] = allocArray( &bytes , sizeof( real ) , 3 , fplim[region][HZ][XHI] - fplim[region][HZ][XLO] + 1 , 
                                                        fplim[region][HZ][YHI] - fplim[region][HZ][YLO] + 1 , 
                                                        fplim[region][HZ][ZHI] - fplim[region][HZ][ZLO] + 1 );
      memory.pmlFields += bytes;
    }
  }

  message( MSG_DEBUG1 , 0 , "  Allocating grid PML adx array\n" );
  adx = allocArray( &bytes , sizeof( real ) , 1 , numCells[XDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML bdx array\n" );
  bdx = allocArray( &bytes , sizeof( real ) , 1 , numCells[XDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML gdx array\n" );
  gdx = allocArray( &bytes , sizeof( real ) , 1 , numCells[XDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML kdx array\n" );
  kdx = allocArray( &bytes , sizeof( real ) , 1 , numCells[XDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML ahx array\n" );
  ahx = allocArray( &bytes , sizeof( real ) , 1 , numCells[XDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML bhx array\n" );
  bhx = allocArray( &bytes , sizeof( real ) , 1 , numCells[XDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML ghx array\n" );
  ghx = allocArray( &bytes , sizeof( real ) , 1 , numCells[XDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML khx array\n" );
  khx = allocArray( &bytes , sizeof( real ) , 1 , numCells[XDIR] );
  memory.pmlCoeffs += bytes;

  message( MSG_DEBUG1 , 0 , "  Allocating grid PML ady array\n" );
  ady = allocArray( &bytes , sizeof( real ) , 1 , numCells[YDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML bdy array\n" );
  bdy = allocArray( &bytes , sizeof( real ) , 1 , numCells[YDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML gdy array\n" );
  gdy = allocArray( &bytes , sizeof( real ) , 1 , numCells[YDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML kdy array\n" );
  kdy = allocArray( &bytes , sizeof( real ) , 1 , numCells[YDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML ahy array\n" );
  ahy = allocArray( &bytes , sizeof( real ) , 1 , numCells[YDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML bhy array\n" );
  bhy = allocArray( &bytes , sizeof( real ) , 1 , numCells[YDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML ghy array\n" );
  ghy = allocArray( &bytes , sizeof( real ) , 1 , numCells[YDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML khy array\n" );
  khy = allocArray( &bytes , sizeof( real ) , 1 , numCells[YDIR] );
  memory.pmlCoeffs += bytes;

  message( MSG_DEBUG1 , 0 , "  Allocating grid PML adz array\n" );
  adz = allocArray( &bytes , sizeof( real ) , 1 , numCells[ZDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML bdz array\n" );
  bdz = allocArray( &bytes , sizeof( real ) , 1 , numCells[ZDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML gdz array\n" );
  gdz = allocArray( &bytes , sizeof( real ) , 1 , numCells[ZDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML kdz array\n" );
  kdz = allocArray( &bytes , sizeof( real ) , 1 , numCells[ZDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML ahz array\n" );
  ahz = allocArray( &bytes , sizeof( real ) , 1 , numCells[ZDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML bhz array\n" );
  bhz = allocArray( &bytes , sizeof( real ) , 1 , numCells[ZDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML ghz array\n" );
  ghz = allocArray( &bytes , sizeof( real ) , 1 , numCells[ZDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML khz array\n" );
  khz = allocArray( &bytes , sizeof( real ) , 1 , numCells[ZDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML ibdx array\n" );
  ibdx = allocArray( &bytes , sizeof( real ) , 1 , numCells[XDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML ibhx array\n" );
  ibhx = allocArray( &bytes , sizeof( real ) , 1 , numCells[XDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML ibdy array\n" );
  ibdy = allocArray( &bytes , sizeof( real ) , 1 , numCells[YDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML ibhy array\n" );
  ibhy = allocArray( &bytes , sizeof( real ) , 1 , numCells[YDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML ibdz array\n" );
  ibdz = allocArray( &bytes , sizeof( real ) , 1 , numCells[ZDIR] );
  memory.pmlCoeffs += bytes;
  message( MSG_DEBUG1 , 0 , "  Allocating grid PML ibhz array\n" );
  ibhz = allocArray( &bytes , sizeof( real ) , 1 , numCells[ZDIR] );
  memory.pmlCoeffs += bytes;

  return;

}

/* Set initial PML field values. */
void clearPml( void )
{

  int i , j , k;

  message( MSG_LOG , 0 , "  Clearing the PML...\n" );

  /* Clear PML arrays. */
  for( int region = XLO ; region <= ZHI ; region++ )
  {

    for ( i = 0 ; i <= fplim[region][EX][XHI] - fplim[region][EX][XLO] ; i++ ) {
      for ( j = 0 ; j <= fplim[region][EX][YHI] - fplim[region][EX][YLO] ; j++ ) {
        for ( k = 0 ; k <= fplim[region][EX][ZHI] - fplim[region][EX][ZLO] ; k++ ) {
          PPx[region][i][j][k] = INITIAL_FIELD_VALUE;
          Px[region][i][j][k] = INITIAL_FIELD_VALUE;
        }
      }
    }

    for ( i = 0 ; i <= fplim[region][EY][XHI] - fplim[region][EY][XLO] ; i++ ) {
      for ( j = 0 ; j <= fplim[region][EY][YHI] - fplim[region][EY][YLO] ; j++ ) {
        for ( k = 0 ; k <= fplim[region][EY][ZHI] - fplim[region][EY][ZLO] ; k++ ) {
          PPy[region][i][j][k] = INITIAL_FIELD_VALUE;
          Py[region][i][j][k] = INITIAL_FIELD_VALUE;
        }
      }
    }

    for ( i = 0 ; i <= fplim[region][EZ][XHI] - fplim[region][EZ][XLO] ; i++ ) {
      for ( j = 0 ; j <= fplim[region][EZ][YHI] - fplim[region][EZ][YLO] ; j++ ) {
        for ( k = 0 ; k <= fplim[region][EZ][ZHI] - fplim[region][EZ][ZLO] ; k++ ) {
          PPz[region][i][j][k] = INITIAL_FIELD_VALUE;
          Pz[region][i][j][k] = INITIAL_FIELD_VALUE;
        }
      }
    }

    for ( i = 0 ; i <= fplim[region][HX][XHI] - fplim[region][HX][XLO] ; i++ ) {
      for ( j = 0 ; j <= fplim[region][HX][YHI] - fplim[region][HX][YLO] ; j++ ) {
        for ( k = 0 ; k <= fplim[region][HX][ZHI] - fplim[region][HX][ZLO] ; k++ ) {
          Bx[region][i][j][k] = INITIAL_FIELD_VALUE;
        }
      }
    }

    for ( i = 0 ; i <= fplim[region][HY][XHI] - fplim[region][HY][XLO] ; i++ ) {
      for ( j = 0 ; j <= fplim[region][HY][YHI] - fplim[region][HY][YLO] ; j++ ) {
        for ( k = 0 ; k <= fplim[region][HY][ZHI] - fplim[region][HY][ZLO] ; k++ ) {
          By[region][i][j][k] = INITIAL_FIELD_VALUE;
        }
      }
    }

    for ( i = 0 ; i <= fplim[region][HZ][XHI] - fplim[region][HZ][XLO] ; i++ ) {
      for ( j = 0 ; j <= fplim[region][HZ][YHI] - fplim[region][HZ][YLO] ; j++ ) {
        for ( k = 0 ; k <= fplim[region][HZ][ZHI] - fplim[region][HZ][ZLO] ; k++ ) {
          Bz[region][i][j][k] = INITIAL_FIELD_VALUE;
        }
      }
    }

  }

  return;

}

/* Update electric field in PML regions. */
void updatePmlEfield( void )
{

  int i , j , k , region;
  int ir , jr , kr;
  int *fplim_rf;
  real **Ex_i, **Ey_i, **Ez_i, **Hx_i, **Hy_i, **Hz_i;
  real *Ex_ij, *Ey_ij, *Ez_ij, *Hx_ij, *Hy_ij, *Hz_ij;
  real **Hz_i1 , **Hy_i1 , *Hz_ij1 , *Hz_i1j , *Hy_i1j , *Hx_ij1;
  real **Px_i, **Py_i, **Pz_i;
  real *Px_ij, *Py_ij, *Pz_ij;
  real **PPx_i, **PPy_i, **PPz_i;
  real *PPx_ij, *PPy_ij, *PPz_ij;
  real ahx_i, ibhx_i , ady_j, bdy_j;
  real adx_i , bdx_i , ahy_j , ibhy_j;

  /* Temporary storage for field arrays. */
  real oldPx, oldPy, oldPz, oldPPx, oldPPy, oldPPz;

  /* Since the 6 PML regions do not coincide with the boundaries of each face we cannot */
  /* easily check which loops to run from the border flags. The update of the PML is */
  /* therefore controlled by the array limits only. */
  for( region = XLO ; region <= ZHI ; region++ )
  {
 
    if( Px[region] != NULL )
    {

      /* Update Ex from H, PPx and Px. */
      fplim_rf = fplim[region][EX];
      #ifdef WITH_OPENMP
        #pragma omp parallel for private( i , j , k , ir , jr , kr , Ex_i , Ex_ij , Hy_i , Hy_ij , Hz_i , Hz_ij , Hz_ij1 , Px_i , Px_ij , PPx_i , PPx_ij , oldPx , oldPPx , ahx_i, ibhx_i , ady_j, bdy_j )
      #endif
      for ( i = fplim_rf[XLO] ; i <= fplim_rf[XHI] ; i++ ) 
      {
        ir = i - fplim_rf[XLO];
        Ex_i = Ex[i];
        Hy_i = Hy[i];
        Hz_i = Hz[i];
        Px_i = Px[region][ir];
        PPx_i = PPx[region][ir];
        ahx_i = ahx[i];  
        ibhx_i = ibhx[i];
        for ( j = fplim_rf[YLO] ; j <= fplim_rf[YHI] ; j++ ) 
        {
          jr = j - fplim_rf[YLO];
          Ex_ij = Ex_i[j];
          Hy_ij = Hy_i[j];
          Hz_ij = Hz_i[j];
          Hz_ij1 = Hz_i[j-1];
          Px_ij = Px_i[jr];
          PPx_ij = PPx_i[jr];
          ady_j = ady[j];
          bdy_j = bdy[j];
          for ( k = fplim_rf[ZLO] ; k <= fplim_rf[ZHI] ; k++ ) 
          {
            kr = k - fplim_rf[ZLO];
            oldPPx = PPx_ij[kr];
            PPx_ij[kr] = ALPHA_EX(i,j,k) * PPx_ij[kr] + BETA_EX(i,j,k) 
              * curl_Hx( Hz_ij[k] , Hz_ij1[k] , Hy_ij[k-1] , Hy_ij[k] , i , j , k );
            oldPx = Px_ij[kr];
            Px_ij[kr] = ady_j * Px_ij[kr] + bdy_j * ( PPx_ij[kr] - oldPPx );
	    CHECK_NOT_VISITED( Ex_ij[k] );
	    Ex_ij[k] = adz[k] * Ex_ij[k] + bdz[k] * ibhx_i * ( Px_ij[kr] - ahx_i * oldPx );
	    MARK_AS_VISITED( Ex_ij[k] ); 
	  } // for k
        } // for j
      } // for i
    } // if
    
    /* Update Ey from H, PPy and Py. */
    if( Py[region] != NULL )
    {
      fplim_rf = fplim[region][EY]; 
      #ifdef WITH_OPENMP
        #pragma omp parallel for private( i , j , k , ir , jr , kr , Ey_i , Ey_ij , Hx_i , Hx_ij , Hz_i , Hz_ij , Hz_i1 , Hz_i1j , Py_i , Py_ij , PPy_i , PPy_ij , oldPy , oldPPy , adx_i , bdx_i , ahy_j , ibhy_j )
      #endif
      for ( i = fplim_rf[XLO] ; i <= fplim_rf[XHI] ; i++ ) 
      {
        ir = i - fplim_rf[XLO];
        Ey_i = Ey[i];
        Hx_i = Hx[i];
        Hz_i = Hz[i];
        Hz_i1 = Hz[i-1];
        Py_i = Py[region][ir];
        PPy_i = PPy[region][ir];
        adx_i = adx[i];
        bdx_i = bdx[i];
        for ( j = fplim_rf[YLO] ; j <= fplim_rf[YHI] ; j++ ) 
        {
          jr = j - fplim_rf[YLO];
          Ey_ij = Ey_i[j];
          Hx_ij = Hx_i[j];
          Hz_ij = Hz_i[j];  
          Hz_i1j = Hz_i1[j];
          Py_ij = Py_i[jr];
          PPy_ij = PPy_i[jr];
          ahy_j = ahy[j];
          ibhy_j = ibhy[j];
          for ( k = fplim_rf[ZLO] ; k <= fplim_rf[ZHI] ; k++ ) 
	  {
            kr = k - fplim_rf[ZLO];
            oldPPy = PPy_ij[kr];
            PPy_ij[kr] = ALPHA_EY(i,j,k) * PPy_ij[kr] + BETA_EY(i,j,k)
	      * curl_Hy( Hx_ij[k] , Hx_ij[k-1] , Hz_i1j[k] , Hz_ij[k] , i , j , k );
            oldPy = Py_ij[kr];
            Py_ij[kr] = adz[k] * Py_ij[kr] + bdz[k] * ( PPy_ij[kr] - oldPPy );
	    CHECK_NOT_VISITED( Ey_ij[k] );
	    Ey_ij[k] = adx_i * Ey_ij[k] + bdx_i * ibhy_j * ( Py_ij[kr] - ahy_j * oldPy );
	    MARK_AS_VISITED( Ey_ij[k] ); 
	  } // for k
        } // for j
      } // for i
    } // if
    
    /* Update Ez from H, PPz and Pz. */
    if( Pz[region] != NULL )
      {
      fplim_rf = fplim[region][EZ]; 
      #ifdef WITH_OPENMP
        #pragma omp parallel for private( i , j , k , ir , jr , kr , Ez_i , Ez_ij , Hx_i , Hx_ij , Hy_i , Hy_ij , Hy_i1 , Hy_i1j , Hx_ij1 , Pz_i , Pz_ij , PPz_i , PPz_ij , oldPz , oldPPz , adx_i , bdx_i , ady_j , bdy_j )
      #endif
      for ( i = fplim_rf[XLO] ; i <= fplim_rf[XHI] ; i++ ) 
      {
        ir = i - fplim_rf[XLO];
        Ez_i = Ez[i];
        Hx_i = Hx[i];
        Hy_i = Hy[i];
        Hy_i1 = Hy[i-1];
        Pz_i = Pz[region][ir];
        PPz_i = PPz[region][ir];
        adx_i = adx[i];
        bdx_i = bdx[i];
        for ( j = fplim_rf[YLO] ; j <= fplim_rf[YHI] ; j++ ) 
        {
          jr = j - fplim_rf[YLO];
          Ez_ij = Ez_i[j];
          Hx_ij = Hx_i[j];
          Hy_ij = Hy_i[j];       
          Hy_i1j = Hy_i1[j];
          Hx_ij1 = Hx_i[j-1];
          Pz_ij = Pz_i[jr];
          PPz_ij = PPz_i[jr];
          ady_j = ady[j];
          bdy_j = bdy[j];
          for ( k = fplim_rf[ZLO] ; k <= fplim_rf[ZHI] ; k++ ) 
	  {
            kr = k - fplim_rf[ZLO];
            oldPPz = PPz_ij[kr];
            PPz_ij[kr] = ALPHA_EZ(i,j,k) * PPz_ij[kr] + BETA_EZ(i,j,k)
     	      * curl_Hz( Hy_ij[k] , Hy_i1j[k] , Hx_ij1[k] , Hx_ij[k] , i , j , k );
            oldPz = Pz_ij[kr];
            Pz_ij[kr] = adx_i * Pz_ij[kr] + bdx_i * ( PPz_ij[kr] - oldPPz );
	    CHECK_NOT_VISITED( Ez_ij[k] );
	    Ez_ij[k] = ady_j * Ez_ij[k] + bdy_j * ibhz[k] * ( Pz_ij[kr] - ahz[k] * oldPz );
	    MARK_AS_VISITED( Ez_ij[k] ); 
	  } // for i
        } // for j
      } // for i
    } // if
  }

  return;

}

/* Update magnetic field in PML regions. */
void updatePmlHfield( void )
{

  int i , j , k , region;
  int ir , jr , kr;
  int *fplim_rf;
  real **Ex_i, **Ey_i, **Ez_i, **Hx_i, **Hy_i, **Hz_i;
  real *Ex_ij, *Ey_ij, *Ez_ij, *Hx_ij, *Hy_ij, *Hz_ij;
  real **Ez_i1 , **Ey_i1 , *Ez_ij1 , *Ez_i1j , *Ex_ij1 , *Ey_i1j;
  real **Bx_i, **By_i, **Bz_i;
  real *Bx_ij, *By_ij, *Bz_ij;
  real ahy_j , bhy_j , adx_i , ibdx_i;
  real ahx_i , bhx_i , ady_j , ibdy_j;
  
  /* Temporary storage for field arrays. */
  real oldBx, oldBy, oldBz; 

  /* Since the 6 PML regions do not coincide with the boundaries of each face we cannot */
  /* easily check which loops to run from the border flags. The update of the PML is */
  /* therefore controlled by the array limits only. */
  for( region = XLO ; region <= ZHI ; region++ )
  {

    /* Update Hx from E and Bx. */
    if( Bx[region] != NULL )
    {
      fplim_rf = fplim[region][HX];
      #ifdef WITH_OPENMP
        #pragma omp parallel for private( i , j , k , ir , jr , kr , Ey_i , Ey_ij , Ez_i , Ez_ij , Hx_i , Hx_ij , Ez_ij1 , Bx_i , Bx_ij , oldBx , ahy_j , bhy_j , adx_i , ibdx_i )
      #endif
      for ( i = fplim_rf[XLO] ; i <= fplim_rf[XHI] ; i++ ) 
      {
        ir = i - fplim_rf[XLO];
        Hx_i = Hx[i];
        Ey_i = Ey[i];
        Ez_i = Ez[i];
        Bx_i = Bx[region][ir];
        adx_i = adx[i];
        ibdx_i = ibdx[i];
        for ( j = fplim_rf[YLO] ; j <= fplim_rf[YHI] ; j++ ) 
        {
          jr = j - fplim_rf[YLO];
          Hx_ij = Hx_i[j];
          Ey_ij = Ey_i[j];
          Ez_ij = Ez_i[j];   
          Ez_ij1 = Ez_i[j+1];
          Bx_ij = Bx_i[jr];
          ahy_j = ahy[j];
          bhy_j = bhy[j];
          for ( k = fplim_rf[ZLO] ; k <= fplim_rf[ZHI] ; k++ ) 
	  {
            kr = k - fplim_rf[ZLO];
            oldBx = Bx_ij[kr];
            Bx_ij[kr] = ahy_j * Bx_ij[kr] + GAMMA_HX(i,j,k) * bhy_j
      	      * curl_Ex( Ey_ij[k+1] , Ey_ij[k] , Ez_ij[k] , Ez_ij1[k] , i , j , k ); 
	    CHECK_NOT_VISITED( Hx_ij[k] );
            Hx_ij[k] = ahz[k] * Hx_ij[k] + bhz[k] * ibdx_i * ( Bx_ij[kr] - adx_i * oldBx );
	    MARK_AS_VISITED( Hx_ij[k] ); 
	  } // for k
        } // for j
      } // for i
    } // if
    
    /* Update Hy from E and By. */
    if( By[region] != NULL )
    {
      fplim_rf = fplim[region][HY];
      #ifdef WITH_OPENMP
        #pragma omp parallel for private( i , j , k , ir , jr , kr , Ex_i , Ex_ij , Ez_i , Ez_ij , Hy_i , Hy_ij , Ez_i1 , Ez_i1j , By_i , By_ij , oldBy , ahx_i , bhx_i , ady_j , ibdy_j )
      #endif
      for ( i = fplim_rf[XLO] ; i <= fplim_rf[XHI] ; i++ ) 
      {
        ir = i - fplim_rf[XLO];
        Hy_i = Hy[i];
        Ex_i = Ex[i];
        Ez_i = Ez[i];
        Ez_i1 = Ez[i+1];
        By_i = By[region][ir];   
        ahx_i = ahx[i];
        bhx_i = bhx[i];
        for ( j = fplim_rf[YLO] ; j <= fplim_rf[YHI] ; j++ ) 
        {
          jr = j - fplim_rf[YLO];
          Hy_ij = Hy_i[j];
          Ex_ij = Ex_i[j];
          Ez_ij = Ez_i[j];   
          Ez_i1j = Ez_i1[j];
          By_ij = By_i[jr];
          ady_j = ady[j];
          ibdy_j = ibdy[j];
          for ( k = fplim_rf[ZLO] ; k <= fplim_rf[ZHI] ; k++ ) 
	  {
            kr = k - fplim_rf[ZLO];
            oldBy = By_ij[kr];
            By_ij[kr] = ahz[k] * By_ij[kr] + GAMMA_HY(i,j,k) * bhz[k]
    	      * curl_Ey( Ez_i1j[k] , Ez_ij[k] , Ex_ij[k] , Ex_ij[k+1] , i , j , k );
	    CHECK_NOT_VISITED( Hy_ij[k] );
	    Hy_ij[k] = ahx_i * Hy_ij[k] + bhx_i * ibdy_j * ( By_ij[kr] - ady_j * oldBy );
	    MARK_AS_VISITED( Hy_ij[k] ); 
  	  } // for k
        } // for j
      } // for i
    } // if
    
    /* Update Hz from E and Bz. */
    if( Bz[region] != NULL )
    {
      fplim_rf = fplim[region][HZ];
      #ifdef WITH_OPENMP
        #pragma omp parallel for private( i , j , k , ir , jr , kr , Ex_i , Ex_ij , Ey_i , Ey_ij , Hz_i , Hz_ij , Ex_ij1 , Ey_i1 , Ey_i1j , Bz_i , Bz_ij , oldBz , ahx_i , bhx_i , ahy_j , bhy_j )
      #endif
      for ( i = fplim_rf[XLO] ; i <= fplim_rf[XHI] ; i++ ) 
      {
        ir = i - fplim_rf[XLO];
        Hz_i = Hz[i];
        Ex_i = Ex[i];
        Ey_i = Ey[i];
        Ey_i1 = Ey[i+1];
        Bz_i = Bz[region][ir];
        ahx_i = ahx[i];
        bhx_i = bhx[i];
        for ( j = fplim_rf[YLO] ; j <= fplim_rf[YHI] ; j++ ) 
        {
          jr = j - fplim_rf[YLO];
          Hz_ij = Hz_i[j];
          Ex_ij = Ex_i[j];
          Ey_ij = Ey_i[j];   
          Ex_ij1 = Ex_i[j+1];
          Ey_i1j = Ey_i1[j];
          Bz_ij = Bz_i[jr];
          ahy_j = ahy[j];
          bhy_j = bhy[j];
          for ( k = fplim_rf[ZLO] ; k <= fplim_rf[ZHI] ; k++ ) 
  	  {
            kr = k - fplim_rf[ZLO];
            oldBz = Bz_ij[kr];
            Bz_ij[kr] = ahx_i * Bz_ij[kr] + GAMMA_HZ(i,j,k) * bhx_i
 	      * curl_Ez( Ex_ij1[k] , Ex_ij[k] , Ey_ij[k] , Ey_i1j[k] , i , j , k );   
	    CHECK_NOT_VISITED( Hz_ij[k] );
	    Hz_ij[k] = ahy_j * Hz_ij[k] + bhy_j * ibdz[k] * ( Bz_ij[kr] - adz[k] * oldBz );
	    MARK_AS_VISITED( Hz_ij[k] ); 
	  } // for k
        } // for j
      } // for i 
    } // if

  }

  return;

}

/* Deallocate PML arrays. */
void deallocPmlArrays( void )
{

  int region;

  message( MSG_DEBUG1 , 0 , "Deallocating the PML...\n" );

  for( region = XLO ; region <= ZHI ; region++ )
  {
    if( outerSurfaceType( region ) == BT_PML )
    {
      message( MSG_DEBUG1 , 0 , "  Deallocating grid PML Px array\n" );
      deallocArray( Px[region] , 3 , fplim[region][EX][XHI] - fplim[region][EX][XLO] + 1 , 
                                         fplim[region][EX][YHI] - fplim[region][EX][YLO] + 1 , 
                                         fplim[region][EX][ZHI] - fplim[region][EX][ZLO] + 1 );
      message( MSG_DEBUG1 , 0 , "  Deallocating grid PML Py array\n" );
      deallocArray( Py[region] , 3 , fplim[region][EY][XHI] - fplim[region][EY][XLO] + 1 , 
                                         fplim[region][EY][YHI] - fplim[region][EY][YLO] + 1 , 
                                         fplim[region][EY][ZHI] - fplim[region][EY][ZLO] + 1 );
      message( MSG_DEBUG1 , 0 , "  Deallocating grid PML Pz array\n" );
      deallocArray( Pz[region] , 3 , fplim[region][EZ][XHI] - fplim[region][EZ][XLO] + 1 , 
                                         fplim[region][EZ][YHI] - fplim[region][EZ][YLO] + 1 , 
                                         fplim[region][EZ][ZHI] - fplim[region][EZ][ZLO] + 1 );
      message( MSG_DEBUG1 , 0 , "  Deallocating grid PML PPx array\n" );
      deallocArray( PPx[region] , 3 , fplim[region][EX][XHI] - fplim[region][EX][XLO] + 1 , 
                                          fplim[region][EX][YHI] - fplim[region][EX][YLO] + 1 , 
                                          fplim[region][EX][ZHI] - fplim[region][EX][ZLO] + 1 );
      message( MSG_DEBUG1 , 0 , "  Deallocating grid PML PPy array\n" );
      deallocArray( PPy[region] , 3 , fplim[region][EY][XHI] - fplim[region][EY][XLO] + 1 , 
                                          fplim[region][EY][YHI] - fplim[region][EY][YLO] + 1 , 
                                          fplim[region][EY][ZHI] - fplim[region][EY][ZLO] + 1 );
      message( MSG_DEBUG1 , 0 , "  Deallocating grid PML PPz array\n" );
      deallocArray( PPz[region] , 3 , fplim[region][EZ][XHI] - fplim[region][EZ][XLO] + 1 , 
                                          fplim[region][EZ][YHI] - fplim[region][EZ][YLO] + 1 , 
                                          fplim[region][EZ][ZHI] - fplim[region][EZ][ZLO] + 1 );
      message( MSG_DEBUG1 , 0 , "  Deallocating grid PML Bx array\n" );
      deallocArray( Bx[region] , 3 , fplim[region][HX][XHI] - fplim[region][HX][XLO] + 1 , 
                                         fplim[region][HX][YHI] - fplim[region][HX][YLO] + 1 , 
                                         fplim[region][HX][ZHI] - fplim[region][HX][ZLO] + 1 );
      message( MSG_DEBUG1 , 0 , "  Deallocating grid PML By array\n" );
      deallocArray( By[region] , 3 , fplim[region][HY][XHI] - fplim[region][HY][XLO] + 1 , 
                                         fplim[region][HY][YHI] - fplim[region][HY][YLO] + 1 , 
                                         fplim[region][HY][ZHI] - fplim[region][HY][ZLO] + 1 );
      message( MSG_DEBUG1 , 0 , "  Deallocating grid PML Bz array\n" );
      deallocArray( Bz[region] , 3 , fplim[region][HZ][XHI] - fplim[region][HZ][XLO] + 1 , 
                                         fplim[region][HZ][YHI] - fplim[region][HZ][YLO] + 1 , 
                                         fplim[region][HZ][ZHI] - fplim[region][HZ][ZLO] + 1 );
    }
  }

  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML adx array\n" );
  deallocArray( adx , numCells[XDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML bdx array\n" );
  deallocArray( bdx , numCells[XDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML gdx array\n" );
  deallocArray( gdx , numCells[XDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML kdx array\n" );
  deallocArray( kdx , numCells[XDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML ahx array\n" );
  deallocArray( ahx , numCells[XDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML bhx array\n" );
  deallocArray( bhx , numCells[XDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML ghx array\n" );
  deallocArray( ghx , numCells[XDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML khx array\n" );
  deallocArray( khx , numCells[XDIR] );

  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML ady array\n" );
  deallocArray( ady , numCells[YDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML bdy array\n" );
  deallocArray( bdy , numCells[YDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML gdy array\n" );
  deallocArray( gdy , numCells[YDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML kdy array\n" );
  deallocArray( kdy , numCells[YDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML ahy array\n" );
  deallocArray( ahy , numCells[YDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML bhy array\n" );
  deallocArray( bhy , numCells[YDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML ghy array\n" );
  deallocArray( ghy , numCells[YDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML khy array\n" );
  deallocArray( khy , numCells[YDIR] );

  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML adz array\n" );
  deallocArray( adz , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML bdz array\n" );
  deallocArray( bdz , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML gdz array\n" );
  deallocArray( gdz , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML kdz array\n" );
  deallocArray( kdz , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML ahz array\n" );
  deallocArray( ahz , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML bhz array\n" );
  deallocArray( bhz , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML ghz array\n" );
  deallocArray( ghz , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML khz array\n" );
  deallocArray( khz , numCells[ZDIR] );
  
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML ibdx array\n" );
  deallocArray( ibdx , numCells[XDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML ibhx array\n" );
  deallocArray( ibhx , numCells[XDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML ibdy array\n" );
  deallocArray( ibdy , numCells[YDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML ibhy array\n" );
  deallocArray( ibhy , numCells[YDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML ibdz array\n" );
  deallocArray( ibdz , numCells[ZDIR] );
  message( MSG_DEBUG1 , 0 , "  Deallocating grid PML ibhz array\n" );
  deallocArray( ibhz , numCells[ZDIR] );

  return;

}

/* Report PML. */
void reportPml( void )
{

  for( MeshFace region = XLO ; region <= ZHI ; region++ )
  {
    message( MSG_LOG , 0 , "  PML %s Region: BBOX=[%d,%d,%d,%d,%d,%d]\n",
             FACE[region] , pbox[region][XLO] , pbox[region][XHI] - 1,
                            pbox[region][YLO] , pbox[region][YHI] - 1,
                            pbox[region][ZLO] , pbox[region][ZHI] - 1 );
  }

  for( MeshFace region = XLO ; region <= ZHI ; region++ )
  {
    for( FieldComponent field = EX ; field <= HZ ; field++ )
    {
      message( MSG_DEBUG1 , 0 , "  PML %s Region, %s Field limits: [%d,%d,%d,%d,%d,%d]\n",
               FACE[region] , FIELD[field] ,
               fplim[region][field][XLO] , fplim[region][field][XHI] ,
               fplim[region][field][YLO] , fplim[region][field][YHI] ,
               fplim[region][field][ZLO] , fplim[region][field][ZHI] );
    }
  }

  return;

}

/* Default PML parameters. */
void setPmlDefaults( int *numLayers , int *order , real *n_eff , real *refCoeff , real *kmax )
{

  *numLayers = 6;
  *order = 4;
  *n_eff = 1.0;
  *refCoeff = -1.0;  // Use optimal profile.
  *kmax = 1.0;

  return;

}
