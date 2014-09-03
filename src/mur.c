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

#include "mur.h"
#include "grid.h"
#include "alloc_array.h"
#include "bounding_box.h"
#include "message.h"
#include "surface.h"
#include "medium.h"
#include "memory.h"
#include "physical.h"

/* 
 * Private data.
 */

/* Bounding boxes for MUR faces. */
static int murbox[6][6];      

/* Field array limits for MUR faces. */
static int fmlim[6][6][6];   

/* MUR phase velocity factor arrays. */
static real zeta[6];

/* 
 * Private method interfaces. 
 */

void deselectAdjacentEdgesByType( bool includeBoundary[6] , MeshFace face , BoundaryType type );

/*
 * Method Implementations.
 */

/* Deselect adjacent edges in include flags if boundary is of given type. */
void deselectAdjacentEdgesByType( bool includeBoundary[6] , MeshFace face , BoundaryType type )
{

  switch( face )
  {
    case XLO:
    case XHI:
      if( outerSurfaceType( YLO ) == type ) includeBoundary[YLO] = false;
      if( outerSurfaceType( YHI ) == type ) includeBoundary[YHI] = false;
      if( outerSurfaceType( ZLO ) == type ) includeBoundary[ZLO] = false;
      if( outerSurfaceType( ZHI ) == type ) includeBoundary[ZHI] = false;            
      break;
    case YLO:
    case YHI:
      if( outerSurfaceType( XLO ) == type ) includeBoundary[XLO] = false;
      if( outerSurfaceType( XHI ) == type ) includeBoundary[XHI] = false;
      if( outerSurfaceType( ZLO ) == type ) includeBoundary[ZLO] = false;
      if( outerSurfaceType( ZHI ) == type ) includeBoundary[ZHI] = false;            
      break;
    case ZLO:
    case ZHI:
      if( outerSurfaceType( XLO ) == type ) includeBoundary[XLO] = false;
      if( outerSurfaceType( XHI ) == type ) includeBoundary[XHI] = false;
      if( outerSurfaceType( YLO ) == type ) includeBoundary[YLO] = false;
      if( outerSurfaceType( YHI ) == type ) includeBoundary[YHI] = false;            
      break;
    default:
      assert( 0 );
      break;
  }
  
  return;
  
}
  
/* Initialise Mur boundaries. */
void initMurBoundaries( void )
{

  real dt;

  message( MSG_LOG , 0 , "\nInitialising Mur boundaries...\n\n" );

  for( MeshFace face = XLO ; face <= ZHI ; face++ )
  {
    bool includeBoundary[6] = { true  , true  , true  , true  , true , true };

    if( outerSurfaceType( face ) == BT_MUR )
    {
      getFaceOfBoundingBox( murbox[face] , gibox , face ); 
      deselectAdjacentEdgesByType( includeBoundary , face , BT_MUR );
      setFieldLimits( murbox[face] , fmlim[face] , includeBoundary );
    }
  }

  /* Set phase velocity parameters on boundaries to free-space.*/
  dt = getGridTimeStep();
  zeta[XLO] = ( c0 * dt - dex[murbox[XLO][XLO]] )   / (  c0 * dt + dex[murbox[XLO][XLO]] );
  zeta[XHI] = ( c0 * dt - dex[murbox[XHI][XHI-1]] ) / (  c0 * dt + dex[murbox[XHI][XHI-1]] );
  zeta[YLO] = ( c0 * dt - dey[murbox[YLO][YLO]] )   / (  c0 * dt + dey[murbox[YLO][YLO]] );
  zeta[YHI] = ( c0 * dt - dey[murbox[YHI][YHI-1]] ) / (  c0 * dt + dey[murbox[YHI][YHI-1]] );
  zeta[ZLO] = ( c0 * dt - dez[murbox[ZLO][ZLO]] )   / (  c0 * dt + dez[murbox[ZLO][ZLO]] );
  zeta[ZHI] = ( c0 * dt - dez[murbox[ZHI][ZHI-1]] ) / (  c0 * dt + dez[murbox[ZHI][ZHI-1]] );
  
  return;

}

/* Update electric field on Mur boundaries. */
/* Must be done before E field updates. */
void updateMurEfield( void )
{

  int i , j , k;
  real Extemp , Eytemp , Eztemp;
  
  /* Mur ABC at XLO. */  
  if( outerSurfaceType( XLO ) == BT_MUR )
  {
  
    i = fmlim[XLO][EY][XLO];
    #ifdef WITH_OPENMP
      #pragma omp parallel for private( j , k , Eytemp )
    #endif
    for ( j = fmlim[XLO][EY][YLO] ; j <= fmlim[XLO][EY][YHI] ; j++ )
    {
      for ( k = fmlim[XLO][EY][ZLO] ; k <= fmlim[XLO][EY][ZHI] ; k++ )
      {
        //zeta = ( sqrt( BETA_EY(i+1,j,k) * GAMMA_EY(i+1,j,k) ) - 1.0 ) /
        //       ( sqrt( BETA_EY(i+1,j,k) * GAMMA_EY(i+1,j,k) ) + 1.0 );
        Eytemp = ALPHA_EY(i+1,j,k) * Ey[i+1][j][k] + BETA_EY(i+1,j,k)
          * curl_Hy( Hx[i+1][j][k] , Hx[i+1][j][k-1] , Hz[i][j][k] , Hz[i+1][j][k] , i + 1 , j , k ); 
        CHECK_NOT_VISITED( Ey[i][j][k] );
        Ey[i][j][k] = Ey[i+1][j][k] + zeta[XLO] * ( Eytemp - Ey[i][j][k] );
        MARK_AS_VISITED( Ey[i][j][k] );
      }
    }

    i = fmlim[XLO][EZ][XLO];
    #ifdef WITH_OPENMP
      #pragma omp parallel for private( j , k , Eztemp )
    #endif
    for ( j = fmlim[XLO][EZ][YLO] ; j <= fmlim[XLO][EZ][YHI] ; j++ )
    {
      for ( k = fmlim[XLO][EZ][ZLO] ; k <= fmlim[XLO][EZ][ZHI] ; k++ )
      {
        //zeta = ( sqrt( BETA_EZ(i+1,j,k) * GAMMA_EZ(i+1,j,k) ) - 1.0 ) /
        //       ( sqrt( BETA_EZ(i+1,j,k) * GAMMA_EZ(i+1,j,k) ) + 1.0 );
        Eztemp = ALPHA_EZ(i+1,j,k) * Ez[i+1][j][k] + BETA_EZ(i+1,j,k)
          * curl_Hz( Hy[i+1][j][k] , Hy[i][j][k] , Hx[i+1][j-1][k] , Hx[i+1][j][k] , i + 1 , j , k );
        CHECK_NOT_VISITED( Ez[i][j][k] );
        Ez[i][j][k] = Ez[i+1][j][k] + zeta[XLO] * ( Eztemp - Ez[i][j][k] );
        MARK_AS_VISITED( Ez[i][j][k] );        
      }
    }

  } // if

  /* Mur ABC at XHI. */
  if( outerSurfaceType( XHI ) == BT_MUR )
  {
  
    i = fmlim[XHI][EY][XHI];
    #ifdef WITH_OPENMP
      #pragma omp parallel for private( j , k , Eytemp )
    #endif
    for ( j = fmlim[XHI][EY][YLO] ; j <= fmlim[XHI][EY][YHI] ; j++ )
    {
      for ( k = fmlim[XHI][EY][ZLO] ; k <= fmlim[XHI][EY][ZHI] ; k++ )
      {
        //zeta = ( sqrt( BETA_EY(i-1,j,k) * GAMMA_EY(i-1,j,k) ) - 1.0 ) /
        //       ( sqrt( BETA_EY(i-1,j,k) * GAMMA_EY(i-1,j,k) ) + 1.0 );
        Eytemp = ALPHA_EY(i-1,j,k) * Ey[i-1][j][k] + BETA_EY(i-1,j,k)
          * curl_Hy( Hx[i-1][j][k] , Hx[i-1][j][k-1] , Hz[i-2][j][k] , Hz[i-1][j][k] , i - 1 , j , k ); 
        CHECK_NOT_VISITED( Ey[i][j][k] );
        Ey[i][j][k] = Ey[i-1][j][k] + zeta[XHI] * ( Eytemp - Ey[i][j][k] );
        MARK_AS_VISITED( Ey[i][j][k] );         
      }
    }

    i = fmlim[XHI][EZ][XHI];
    #ifdef WITH_OPENMP
      #pragma omp parallel for private( j , k , Eztemp )
    #endif
    for ( j = fmlim[XHI][EZ][YLO] ; j <= fmlim[XHI][EZ][YHI] ; j++ )
    {
      for ( k = fmlim[XHI][EZ][ZLO] ; k <= fmlim[XHI][EZ][ZHI] ; k++ )
      {
        //zeta = ( sqrt( BETA_EZ(i-1,j,k) * GAMMA_EZ(i-1,j,k) ) - 1.0 ) /
        //       ( sqrt( BETA_EZ(i-1,j,k) * GAMMA_EZ(i-1,j,k) ) + 1.0 );
        Eztemp = ALPHA_EZ(i-1,j,k) * Ez[i-1][j][k] + BETA_EZ(i-1,j,k)
          * curl_Hz( Hy[i-1][j][k] , Hy[i-2][j][k] , Hx[i-1][j-1][k] , Hx[i-1][j][k] , i - 1 , j , k );
        CHECK_NOT_VISITED( Ez[i][j][k] );
        Ez[i][j][k] = Ez[i-1][j][k] + zeta[XHI] * ( Eztemp - Ez[i][j][k] );
        MARK_AS_VISITED( Ez[i][j][k] );          
      }
    }

  } // if

  /* Mur ABC at YLO. */ 
  if( outerSurfaceType( YLO ) == BT_MUR )
  {
    
    j = fmlim[YLO][EX][YLO];
    #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , k , Extemp )
    #endif
    for ( i = fmlim[YLO][EX][XLO] ; i <= fmlim[YLO][EX][XHI] ; i++ )
    {
      for ( k = fmlim[YLO][EX][ZLO] ; k <= fmlim[YLO][EX][ZHI] ; k++ )
      {
        //zeta = ( sqrt( BETA_EX(i,j+1,k) * GAMMA_EX(i,j+1,k) ) - 1.0 ) /
        //       ( sqrt( BETA_EX(i,j+1,k) * GAMMA_EX(i,j+1,k) ) + 1.0 );
        Extemp = ALPHA_EX(i,j+1,k) * Ex[i][j+1][k] + BETA_EX(i,j+1,k)
          * curl_Hx( Hz[i][j+1][k] , Hz[i][j][k] , Hy[i][j+1][k-1] , Hy[i][j+1][k] , i , j + 1 , k );
        CHECK_NOT_VISITED( Ex[i][j][k] );          
        Ex[i][j][k] = Ex[i][j+1][k] + zeta[YLO] * ( Extemp - Ex[i][j][k] );
        MARK_AS_VISITED( Ex[i][j][k] );          
      }
    }

    j = fmlim[YLO][EZ][YLO];
    #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , k , Eztemp )
    #endif
    for ( i = fmlim[YLO][EZ][XLO] ; i <= fmlim[YLO][EZ][XHI] ; i++ )
    {
      for ( k = fmlim[YLO][EZ][ZLO] ; k <= fmlim[YLO][EZ][ZHI] ; k++ )
      {
        //zeta = ( sqrt( BETA_EZ(i,j+1,k) * GAMMA_EZ(i,j+1,k) ) - 1.0 ) /
        //       ( sqrt( BETA_EZ(i,j+1,k) * GAMMA_EZ(i,j+1,k) ) + 1.0 );
        Eztemp = ALPHA_EZ(i,j+1,k) * Ez[i][j+1][k] + BETA_EZ(i,j+1,k)
          * curl_Hz( Hy[i][j+1][k] , Hy[i-1][j+1][k] , Hx[i][j][k] , Hx[i][j+1][k] , i , j + 1 , k );
        CHECK_NOT_VISITED( Ez[i][j][k] );
        Ez[i][j][k] = Ez[i][j+1][k] + zeta[YLO] * ( Eztemp - Ez[i][j][k] );
        MARK_AS_VISITED( Ez[i][j][k] );          
      }
    }
 
  } // if
  
  /* Mur ABC at YHI. */
  if( outerSurfaceType( YHI ) == BT_MUR )
  {
    
    j = fmlim[YHI][EX][YHI];
    #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , k , Extemp )
    #endif
    for ( i = fmlim[YHI][EX][XLO] ; i <= fmlim[YHI][EX][XHI] ; i++ )
    {
      for (k = fmlim[YHI][EX][ZLO] ; k <= fmlim[YHI][EX][ZHI] ; k++ )
      {
        //zeta = ( sqrt( BETA_EX(i,j-1,k) * GAMMA_EX(i,j-1,k) ) - 1.0 ) /
        //       ( sqrt( BETA_EX(i,j-1,k) * GAMMA_EX(i,j-1,k) ) + 1.0 );
        Extemp = ALPHA_EX(i,j-1,k) * Ex[i][j-1][k] + BETA_EX(i,j-1,k)
          * curl_Hx( Hz[i][j-1][k] , Hz[i][j-2][k] , Hy[i][j-1][k-1] , Hy[i][j-1][k] , i , j - 1 , k );
        CHECK_NOT_VISITED( Ex[i][j][k] ); 
        Ex[i][j][k] = Ex[i][j-1][k] + zeta[YHI] * ( Extemp - Ex[i][j][k] );
        MARK_AS_VISITED( Ex[i][j][k] );              
      }
    }

    j = fmlim[YHI][EZ][YHI]; 
    #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , k , Eztemp )
    #endif
    for ( i = fmlim[YHI][EZ][XLO] ; i <= fmlim[YHI][EZ][XHI] ; i++ )
    {
      for (k = fmlim[YHI][EZ][ZLO] ; k <= fmlim[YHI][EZ][ZHI] ; k++ )
      {
        //zeta = ( sqrt( BETA_EZ(i,j-1,k) * GAMMA_EZ(i,j-1,k) ) - 1.0 ) /
        //       ( sqrt( BETA_EZ(i,j-1,k) * GAMMA_EZ(i,j-1,k) ) + 1.0 );
        Eztemp = ALPHA_EZ(i,j-1,k) * Ez[i][j-1][k] + BETA_EZ(i,j-1,k)
          * curl_Hz( Hy[i][j-1][k] , Hy[i-1][j-1][k] , Hx[i][j-2][k] , Hx[i][j-1][k] , i , j - 1 , k );
        CHECK_NOT_VISITED( Ez[i][j][k] );
        Ez[i][j][k] = Ez[i][j-1][k] + zeta[YHI] * ( Eztemp - Ez[i][j][k] );
        MARK_AS_VISITED( Ez[i][j][k] );
      }
    }

  }
  
  /* Mur ABC at ZLO. */ 
  if( outerSurfaceType( ZLO ) == BT_MUR )
  {
    
    k = fmlim[ZLO][EY][ZLO];
    #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , j , Eytemp )
    #endif
    for ( i = fmlim[ZLO][EY][XLO] ; i <= fmlim[ZLO][EY][XHI] ; i++ )
    {
      for ( j = fmlim[ZLO][EY][YLO] ; j <= fmlim[ZLO][EY][YHI] ; j++ )
      {
        //zeta = ( sqrt( BETA_EY(i,j,k+1) * GAMMA_EY(i,j,k+1) ) - 1.0 ) /
        //       ( sqrt( BETA_EY(i,j,k+1) * GAMMA_EY(i,j,k+1) ) + 1.0 );
        Eytemp = ALPHA_EY(i,j,k+1) * Ey[i][j][k+1] + BETA_EY(i,j,k+1)
          * curl_Hy( Hx[i][j][k+1] , Hx[i][j][k] , Hz[i-1][j][k+1] , Hz[i][j][k+1] , i , j , k + 1 );
        CHECK_NOT_VISITED( Ey[i][j][k] );         
        Ey[i][j][k] = Ey[i][j][k+1] + zeta[ZLO] * ( Eytemp - Ey[i][j][k] );
        MARK_AS_VISITED( Ey[i][j][k] );        
      }
    }

    k = fmlim[ZLO][EX][ZLO];
    #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , j , Extemp )
    #endif
    for ( i = fmlim[ZLO][EX][XLO] ; i <= fmlim[ZLO][EX][XHI] ; i++ )
    {
      for ( j = fmlim[ZLO][EX][YLO] ; j <= fmlim[ZLO][EX][YHI] ; j++ )
      {
        //zeta = ( sqrt( BETA_EX(i,j,k+1) * GAMMA_EX(i,j,k+1) ) - 1.0 ) /
        //       ( sqrt( BETA_EX(i,j,k+1) * GAMMA_EX(i,j,k+1) ) + 1.0 );
        Extemp = ALPHA_EX(i,j,k+1) * Ex[i][j][k+1] + BETA_EX(i,j,k+1)
          * curl_Hx( Hz[i][j][k+1] , Hz[i][j-1][k+1] , Hy[i][j][k] , Hy[i][j][k+1] , i , j , k + 1 );
        CHECK_NOT_VISITED( Ex[i][j][k] );  
        Ex[i][j][k] = Ex[i][j][k+1] + zeta[ZLO] * ( Extemp - Ex[i][j][k] );
        MARK_AS_VISITED( Ex[i][j][k] );          
      }
    }

  }

  /* Mur ABC at ZHI. */  
  if( outerSurfaceType( ZHI ) == BT_MUR )
  {
    
    k = fmlim[ZHI][EY][ZHI];
    #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , j , Eytemp )
    #endif
    for ( i = fmlim[ZHI][EY][XLO] ; i <= fmlim[ZHI][EY][XHI]; i++ )
    {
      for ( j = fmlim[ZHI][EY][YLO] ; j <= fmlim[ZHI][EY][YHI] ; j++ )
      {
        //zeta = ( sqrt( BETA_EY(i,j,k-1) * GAMMA_EY(i,j,k-1) ) - 1.0 ) /
        //       ( sqrt( BETA_EY(i,j,k-1) * GAMMA_EY(i,j,k-1) ) + 1.0 );
        Eytemp = ALPHA_EY(i,j,k-1) * Ey[i][j][k-1] + BETA_EY(i,j,k-1)
          * curl_Hy( Hx[i][j][k-1] , Hx[i][j][k-2] , Hz[i-1][j][k-1] , Hz[i][j][k-1] , i , j , k - 1 ); 
        CHECK_NOT_VISITED( Ey[i][j][k] );          
        Ey[i][j][k] = Ey[i][j][k-1] + zeta[ZHI] * ( Eytemp - Ey[i][j][k] );
        MARK_AS_VISITED( Ey[i][j][k] );         
      }
    }

    k = fmlim[ZHI][EX][ZHI];
    #ifdef WITH_OPENMP
      #pragma omp parallel for private( i , j , Extemp )
    #endif
    for ( i = fmlim[ZHI][EX][XLO] ; i <= fmlim[ZHI][EX][XHI]; i++ )
    {
      for ( j = fmlim[ZHI][EX][YLO] ; j <= fmlim[ZHI][EX][YHI] ; j++ )
      {
        //zeta = ( sqrt ( BETA_EX(i,j,k-1) * GAMMA_EX(i,j,k-1) ) - 1.0 ) /
        //       ( sqrt ( BETA_EX(i,j,k-1) * GAMMA_EX(i,j,k-1) ) + 1.0 );
        Extemp = ALPHA_EX(i,j,k-1) * Ex[i][j][k-1] + BETA_EX(i,j,k-1)
          * curl_Hx( Hz[i][j][k-1] , Hz[i][j-1][k-1] , Hy[i][j][k-2] , Hy[i][j][k-1] , i , j , k - 1 );
        CHECK_NOT_VISITED( Ex[i][j][k] ); 
        Ex[i][j][k] = Ex[i][j][k-1] + zeta[ZHI] * ( Extemp - Ex[i][j][k] );
        MARK_AS_VISITED( Ex[i][j][k] );
      }
    }
   
  }
  
  return;

}

/* Update magnetic field on Mur boundaries. */
/* These are only necessary to make the normal fields have correct values. */
/* They can be omitted without affecting the behaviour of the boundary as */
/* they are not used to update any other fields that can enter the grid. */
void updateMurHfield( void )
{
  
  int i , j , k;
  real **Ex_i, **Ey_i, **Ez_i, **Hx_i, **Hy_i, **Hz_i;
  real *Ex_ij, *Ey_ij, *Ez_ij, *Hx_ij, *Hy_ij, *Hz_ij;
  real **Ez_i1 , **Ey_i1 , *Ez_ij1 , *Ez_i1j , *Ex_ij1 , *Ey_i1j;

  /* Mur ABC at XLO. */  
  if( outerSurfaceType( XLO ) == BT_MUR )
  {
    i = fmlim[XLO][HX][XLO];
    Hx_i = Hx[i];
    Ey_i = Ey[i];
    Ez_i = Ez[i];   
    for ( j = fmlim[XLO][HX][YLO] ; j <= fmlim[XLO][HX][YHI] ; j++ ) 
    {
      Hx_ij = Hx_i[j];
      Ey_ij = Ey_i[j];
      Ez_ij = Ez_i[j];   
      Ez_ij1 = Ez_i[j+1];
      for ( k = fmlim[XLO][HX][ZLO] ; k <= fmlim[XLO][HX][ZHI] ; k++ ) 
      {
        CHECK_NOT_VISITED( Hx_ij[k] );
        Hx_ij[k] = Hx_ij[k] + GAMMA_HX(i,j,k)
          * curl_Ex( Ey_ij[k+1] , Ey_ij[k] , Ez_ij[k] , Ez_ij1[k] , i , j , k ); 
        MARK_AS_VISITED( Hx_ij[k] );
      }
    }
  } // if

  /* Mur ABC at XHI. */  
  if( outerSurfaceType( XHI ) == BT_MUR )
  {
    i = fmlim[XHI][HX][XLO];
    Hx_i = Hx[i];
    Ey_i = Ey[i];
    Ez_i = Ez[i];   
    for ( j = fmlim[XHI][HX][YLO] ; j <= fmlim[XHI][HX][YHI] ; j++ ) 
    {
      Hx_ij = Hx_i[j];
      Ey_ij = Ey_i[j];
      Ez_ij = Ez_i[j];   
      Ez_ij1 = Ez_i[j+1];
      for ( k = fmlim[XHI][HX][ZLO] ; k <= fmlim[XHI][HX][ZHI] ; k++ ) 
      {
        CHECK_NOT_VISITED( Hx_ij[k] );
        Hx_ij[k] = Hx_ij[k] + GAMMA_HX(i,j,k)
          * curl_Ex( Ey_ij[k+1] , Ey_ij[k] , Ez_ij[k] , Ez_ij1[k] , i , j , k ); 
        MARK_AS_VISITED( Hx_ij[k] );
      }
    }    
  } // if

  /* Mur ABC at YLO. */  
  if( outerSurfaceType( YLO ) == BT_MUR )
  {
    for ( i = fmlim[YLO][HY][XLO] ; i <= fmlim[YLO][HY][XHI] ; i++ ) 
    {
      Hy_i = Hy[i];
      Ex_i = Ex[i];
      Ez_i = Ez[i];
      Ez_i1 = Ez[i+1];
      j = fmlim[YLO][HY][YLO];
      Hy_ij = Hy_i[j];
      Ex_ij = Ex_i[j];
      Ez_ij = Ez_i[j];   
      Ez_i1j = Ez_i1[j];
      for ( k = fmlim[YLO][HY][ZLO] ; k <= fmlim[YLO][HY][ZHI] ; k++ ) 
      {
        CHECK_NOT_VISITED( Hy_ij[k] );
        Hy_ij[k] = Hy_ij[k] + GAMMA_HY(i,j,k)
          * curl_Ey( Ez_i1j[k] , Ez_ij[k] , Ex_ij[k] , Ex_ij[k+1] , i , j , k );
        MARK_AS_VISITED( Hy_ij[k] );
      }
    }    
  } // if
  
  /* Mur ABC at YHI. */  
  if( outerSurfaceType( YHI ) == BT_MUR )
  {
    for ( i = fmlim[YHI][HY][XLO] ; i <= fmlim[YHI][HY][XHI] ; i++ ) 
    {
      Hy_i = Hy[i];
      Ex_i = Ex[i];
      Ez_i = Ez[i];
      Ez_i1 = Ez[i+1];
      j = fmlim[YHI][HY][YHI];
      Hy_ij = Hy_i[j];
      Ex_ij = Ex_i[j];
      Ez_ij = Ez_i[j];   
      Ez_i1j = Ez_i1[j];
      for ( k = fmlim[YHI][HY][ZLO] ; k <= fmlim[YHI][HY][ZHI] ; k++ ) 
      {
        CHECK_NOT_VISITED( Hy_ij[k] );
        Hy_ij[k] = Hy_ij[k] + GAMMA_HY(i,j,k)
          * curl_Ey( Ez_i1j[k] , Ez_ij[k] , Ex_ij[k] , Ex_ij[k+1] , i , j , k );
        MARK_AS_VISITED( Hy_ij[k] );
      }
    }    
  } // if

  /* Mur ABC at ZLO. */  
  if( outerSurfaceType( ZLO ) == BT_MUR )
  {
    for ( i = fmlim[ZLO][HZ][XLO] ; i <= fmlim[ZLO][HZ][XHI] ; i++ ) 
    {
      Hz_i = Hz[i];
      Ex_i = Ex[i];
      Ey_i = Ey[i];
      Ey_i1 = Ey[i+1];
      for ( j = fmlim[ZLO][HZ][YLO] ; j <= fmlim[ZLO][HZ][YHI] ; j++ ) 
      {
        Hz_ij = Hz_i[j];
        Ex_ij = Ex_i[j];
        Ey_ij = Ey_i[j];   
        Ex_ij1 = Ex_i[j+1];
        Ey_i1j = Ey_i1[j];
        k = fmlim[ZLO][HZ][ZLO];
        CHECK_NOT_VISITED( Hz_ij[k] );
        Hz_ij[k] = Hz_ij[k] + GAMMA_HZ(i,j,k)
          * curl_Ez( Ex_ij1[k] , Ex_ij[k] , Ey_ij[k] , Ey_i1j[k] , i , j , k );
        MARK_AS_VISITED( Hz_ij[k] );
      }
    }  
  } // if

  /* Mur ABC at ZHI. */  
  if( outerSurfaceType( ZHI ) == BT_MUR )
  {
    for ( i = fmlim[ZHI][HZ][XLO] ; i <= fmlim[ZHI][HZ][XHI] ; i++ ) 
    {
      Hz_i = Hz[i];
      Ex_i = Ex[i];
      Ey_i = Ey[i];
      Ey_i1 = Ey[i+1];
      for ( j = fmlim[ZHI][HZ][YLO] ; j <= fmlim[ZHI][HZ][YHI] ; j++ ) 
      {
        Hz_ij = Hz_i[j];
        Ex_ij = Ex_i[j];
        Ey_ij = Ey_i[j];   
        Ex_ij1 = Ex_i[j+1];
        Ey_i1j = Ey_i1[j];
        k = fmlim[ZHI][HZ][ZHI];
        CHECK_NOT_VISITED( Hz_ij[k] );
        Hz_ij[k] = Hz_ij[k] + GAMMA_HZ(i,j,k)
          * curl_Ez( Ex_ij1[k] , Ex_ij[k] , Ey_ij[k] , Ey_i1j[k] , i , j , k );
        MARK_AS_VISITED( Hz_ij[k] );
      }
    }  
  } // if

  return;
  
}

/* Deallocate Mur arrays. */
void deallocMurArrays( void )
{

  int region;

  message( MSG_DEBUG1 , 0 , "Deallocating Mur boundaries...\n" );

  for( region = XLO ; region <= ZHI ; region++ )
  {
    if( outerSurfaceType( region ) == BT_MUR )
    {
      ;  
    }
  }

  return;

}
