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

#include "debye.h"
#include "utlist.h"
#include "alloc_array.h"
#include "grid.h"
#include "message.h"
#include "bounding_box.h"
#include "memory.h"
#include "physical.h"


/* 
 * Debye class. 
 */

typedef struct Debye_t {

  int flim[3][6];             // Electric field limits on grid.
  MediumItem *medium;         // Pointer to medium model.
  double complex ****Jpolx;   // Array of Jx[i][j][k][p] polarisation currents.
  double complex ****Jpoly;   // Array of Jy[i][j][k][p] polarisation currents.
  double complex ****Jpolz;   // Array of Jz[i][j][k][p] polarisation currents.
  real ***lastEx;             // Cached Ex[i][j][k] from last time-step.
  real ***lastEy;             // Cached Ey[i][j][k] from last time-step.
  real ***lastEz;             // Cached Ez[i][j][k] from last time-step.

} DebyeItem;

/* 
 * Private data.
 */

/* Number of Debye blocks. */
static BlockIndex numDebyeBlock = 0;

/* Debye block array. */
static DebyeItem *debyeArray = NULL;         

/* 
 * Private method interfaces. 
 */

/*
 * Method Implementations.
 */

/* Initialise Debye blocks. */
void initDebyeBlocks( BlockIndex number , BlockItem *blockList )
{

  int i , j , k;
  int ii , jj , kk;
  FieldComponent field;
  BlockItem *item = NULL;
  BlockIndex block = 0;
  unsigned long bytes = 0;
  int gbbox[6];
  int flim[6][6];
  bool includeBoundary[6] = { true , true , true , true , true , true };
  int poleIdx = 0;

  message( MSG_LOG , 0 , "\nInitialising Debye blocks...\n\n" );

  numDebyeBlock = number;

  /* Allocate Debye blocks. */
  message( MSG_DEBUG1 , 0 , "  Allocating Debye block array\n" );
  debyeArray = allocArray( &bytes , sizeof( DebyeItem ) , 1 , numDebyeBlock );
  memory.blocks += bytes;

  /* Set up Debye blocks. */
  block = 0;
  DL_FOREACH( blockList , item ) 
  {
    switch( getMediumType( item->mediumNumber ) )
    {
    case MT_DEBYE:
      debyeArray[block].medium = getMedium( item->mediumNumber );
      offsetBoundingBox( gbbox , item->mbbox , gibox );
      faceMask2boolArray( includeBoundary , item->mask );
      setFieldLimits( gbbox , flim , includeBoundary );
      for( field = EX ; field <= EZ ; field++ )
        for( MeshFace face = XLO ; face <= ZHI ; face++ )
          debyeArray[block].flim[field][face] = flim[field][face];

      message( MSG_DEBUG3 , 0 , "  Setting Debye block on [%d,%d,%d,%d,%d,%d]/[%d,%d,%d,%d,%d,%d]\n" ,
        item->mbbox[XLO] , item->mbbox[XHI] , item->mbbox[YLO] , item->mbbox[YHI] , 
        item->mbbox[ZLO] , item->mbbox[ZHI] , gbbox[XLO] , gbbox[XHI] , 
        gbbox[YLO] , gbbox[YHI] , gbbox[ZLO] , gbbox[ZHI] );
      message( MSG_DEBUG3 , 0 , "    EX FLIM=[%d,%d,%d,%d,%d,%d]\n" , flim[EX][XLO] , flim[EX][XHI] , flim[EX][YLO] , flim[EX][YHI] , flim[EX][ZLO] , flim[EX][ZHI] );
      message( MSG_DEBUG3 , 0 , "    EY FLIM=[%d,%d,%d,%d,%d,%d]\n" , flim[EY][XLO] , flim[EY][XHI] , flim[EY][YLO] , flim[EY][YHI] , flim[EY][ZLO] , flim[EY][ZHI] );
      message( MSG_DEBUG3 , 0 , "    EZ FLIM=[%d,%d,%d,%d,%d,%d]\n" , flim[EZ][XLO] , flim[EZ][XHI] , flim[EZ][YLO] , flim[EZ][YHI] , flim[EZ][ZLO] , flim[EZ][ZHI] );
  
      /* Allocate and initialise arrays for polarisation currents and last electric field value caches. */
      field = EX;
      debyeArray[block].Jpolx = allocArray( &bytes , sizeof( double complex ) , 4 , 
                                            debyeArray[block].flim[field][XHI] - debyeArray[block].flim[field][XLO] + 1 ,
                                            debyeArray[block].flim[field][YHI] - debyeArray[block].flim[field][YLO] + 1 ,
                                            debyeArray[block].flim[field][ZHI] - debyeArray[block].flim[field][ZLO] + 1 ,
                                            debyeArray[block].medium->numPoles );
      memory.blocks += bytes;
      debyeArray[block].lastEx = allocArray( &bytes , sizeof( real ) , 3 , 
                                            debyeArray[block].flim[field][XHI] - debyeArray[block].flim[field][XLO] + 1 ,
                                            debyeArray[block].flim[field][YHI] - debyeArray[block].flim[field][YLO] + 1 ,
                                            debyeArray[block].flim[field][ZHI] - debyeArray[block].flim[field][ZLO] + 1 );
      memory.blocks += bytes;  
      for( i = debyeArray[block].flim[field][XLO] , ii = 0 ; i <= debyeArray[block].flim[field][XHI] ; i++ , ii++ )
        for( j = debyeArray[block].flim[field][YLO] , jj = 0 ; j <= debyeArray[block].flim[field][YHI] ; j++ , jj++ )
          for( k = debyeArray[block].flim[field][ZLO] , kk = 0 ; k <= debyeArray[block].flim[field][ZHI] ; k++ , kk++ )
          {
            for( poleIdx = 0 ; poleIdx < debyeArray[block].medium->numPoles; poleIdx++ )
              debyeArray[block].Jpolx[ii][jj][kk][poleIdx] = 0.0 + I * 0.0;
            debyeArray[block].lastEx[ii][jj][kk] = 0.0;
          }
          
      field = EY;
      debyeArray[block].Jpoly = allocArray( &bytes , sizeof( double complex ) , 4 , 
                                            debyeArray[block].flim[field][XHI] - debyeArray[block].flim[field][XLO] + 1 ,
                                            debyeArray[block].flim[field][YHI] - debyeArray[block].flim[field][YLO] + 1 ,
                                            debyeArray[block].flim[field][ZHI] - debyeArray[block].flim[field][ZLO] + 1 ,
                                            debyeArray[block].medium->numPoles );
      memory.blocks += bytes;
      debyeArray[block].lastEy = allocArray( &bytes , sizeof( real ) , 3 , 
                                            debyeArray[block].flim[field][XHI] - debyeArray[block].flim[field][XLO] + 1 ,
                                            debyeArray[block].flim[field][YHI] - debyeArray[block].flim[field][YLO] + 1 ,
                                            debyeArray[block].flim[field][ZHI] - debyeArray[block].flim[field][ZLO] + 1 );
      memory.blocks += bytes; 
      for( i = debyeArray[block].flim[field][XLO] , ii = 0 ; i <= debyeArray[block].flim[field][XHI] ; i++ , ii++ )
        for( j = debyeArray[block].flim[field][YLO] , jj = 0 ; j <= debyeArray[block].flim[field][YHI] ; j++ , jj++ )
          for( k = debyeArray[block].flim[field][ZLO] , kk = 0 ; k <= debyeArray[block].flim[field][ZHI] ; k++ , kk++ )
          {
            for( poleIdx = 0 ; poleIdx < debyeArray[block].medium->numPoles; poleIdx++ )
              debyeArray[block].Jpoly[ii][jj][kk][poleIdx] = 0.0 + I * 0.0;
            debyeArray[block].lastEy[ii][jj][kk] = 0.0;
          }
          
      field = EZ;
      debyeArray[block].Jpolz = allocArray( &bytes , sizeof( double complex ) , 4 , 
                                            debyeArray[block].flim[field][XHI] - debyeArray[block].flim[field][XLO] + 1 ,
                                            debyeArray[block].flim[field][YHI] - debyeArray[block].flim[field][YLO] + 1 ,
                                            debyeArray[block].flim[field][ZHI] - debyeArray[block].flim[field][ZLO] + 1 ,
                                            debyeArray[block].medium->numPoles );
      memory.blocks += bytes;        
      debyeArray[block].lastEz = allocArray( &bytes , sizeof( real ) , 3 , 
                                            debyeArray[block].flim[field][XHI] - debyeArray[block].flim[field][XLO] + 1 ,
                                            debyeArray[block].flim[field][YHI] - debyeArray[block].flim[field][YLO] + 1 ,
                                            debyeArray[block].flim[field][ZHI] - debyeArray[block].flim[field][ZLO] + 1 );
      memory.blocks += bytes; 
      for( i = debyeArray[block].flim[field][XLO] , ii = 0 ; i <= debyeArray[block].flim[field][XHI] ; i++ , ii++ )
        for( j = debyeArray[block].flim[field][YLO] , jj = 0 ; j <= debyeArray[block].flim[field][YHI] ; j++ , jj++ )
          for( k = debyeArray[block].flim[field][ZLO] , kk = 0 ; k <= debyeArray[block].flim[field][ZHI] ; k++ , kk++ )
          {
            for( poleIdx = 0 ; poleIdx < debyeArray[block].medium->numPoles; poleIdx++ )
              debyeArray[block].Jpolz[ii][jj][kk][poleIdx] = 0.0 + I * 0.0;
            debyeArray[block].lastEz[ii][jj][kk] = 0.0;
          }
    
      block++;
          
      break;
    default:
      break;
    }
  }

  return;

}

/* Deallocate Debye blocks. */
void deallocDebyeBlocks( void )
{

  FieldComponent field;

  message( MSG_DEBUG1 , 0 , "Deallocating Debye blocks...\n" );

  for( BlockIndex block = 0 ; block < numDebyeBlock ; block++ )
  {
    /* Deallocate arrays for polarisation currents and last electric field value caches. */
    field = EX;
    deallocArray( debyeArray[block].Jpolx , 4 , debyeArray[block].flim[field][XHI] - debyeArray[block].flim[field][XLO] + 1 ,
                                                debyeArray[block].flim[field][YHI] - debyeArray[block].flim[field][YLO] + 1 ,
                                                debyeArray[block].flim[field][ZHI] - debyeArray[block].flim[field][ZLO] + 1 ,
                                                debyeArray[block].medium->numPoles );
    deallocArray( debyeArray[block].lastEx , 3 , debyeArray[block].flim[field][XHI] - debyeArray[block].flim[field][XLO] + 1 ,
                                                 debyeArray[block].flim[field][YHI] - debyeArray[block].flim[field][YLO] + 1 ,
                                                 debyeArray[block].flim[field][ZHI] - debyeArray[block].flim[field][ZLO] + 1 );
    field = EY;
    deallocArray( debyeArray[block].Jpoly , 4 , debyeArray[block].flim[field][XHI] - debyeArray[block].flim[field][XLO] + 1 ,
                                                debyeArray[block].flim[field][YHI] - debyeArray[block].flim[field][YLO] + 1 ,
                                                debyeArray[block].flim[field][ZHI] - debyeArray[block].flim[field][ZLO] + 1 ,
                                                debyeArray[block].medium->numPoles );
    deallocArray( debyeArray[block].lastEy , 3 , debyeArray[block].flim[field][XHI] - debyeArray[block].flim[field][XLO] + 1 ,
                                                 debyeArray[block].flim[field][YHI] - debyeArray[block].flim[field][YLO] + 1 ,
                                                 debyeArray[block].flim[field][ZHI] - debyeArray[block].flim[field][ZLO] + 1 );
    field = EZ;
    deallocArray( debyeArray[block].Jpolz , 4 , debyeArray[block].flim[field][XHI] - debyeArray[block].flim[field][XLO] + 1 ,
                                                debyeArray[block].flim[field][YHI] - debyeArray[block].flim[field][YLO] + 1 ,
                                                debyeArray[block].flim[field][ZHI] - debyeArray[block].flim[field][ZLO] + 1 ,
                                                debyeArray[block].medium->numPoles );
    deallocArray( debyeArray[block].lastEz , 3 , debyeArray[block].flim[field][XHI] - debyeArray[block].flim[field][XLO] + 1 ,
                                                 debyeArray[block].flim[field][YHI] - debyeArray[block].flim[field][YLO] + 1 ,
                                                 debyeArray[block].flim[field][ZHI] - debyeArray[block].flim[field][ZLO] + 1 );

  }

  deallocArray( debyeArray , 1 , numDebyeBlock );

  return;

}

/* Debye E field update. Must come before standatd E field update. */
void updateDebyeBlocksEfield( void )
{

  int i , j , k;
  int ii , jj , kk;
  int poleIdx;
  double complex Jsum = 0.0 + I * 0.0;
  FieldComponent field;

  for( BlockIndex block = 0 ; block < numDebyeBlock ; block++ )
  { 

    /* Update Jpolx. */
    field = EX;
    for( i = debyeArray[block].flim[field][XLO] , ii = 0 ; i <= debyeArray[block].flim[field][XHI] ; i++ , ii++ )
      for( j = debyeArray[block].flim[field][YLO] , jj = 0 ; j <= debyeArray[block].flim[field][YHI] ; j++ , jj++ )
        for( k = debyeArray[block].flim[field][ZLO] , kk = 0 ; k <= debyeArray[block].flim[field][ZHI] ; k++ , kk++ )
        {
          Jsum = 0.0 + I * 0.0;
          for( poleIdx = 0 ; poleIdx < debyeArray[block].medium->numPoles; poleIdx++ )
            Jsum = Jsum + ( 1 + debyeArray[block].medium->dalpha[poleIdx] ) * debyeArray[block].Jpolx[ii][jj][kk][poleIdx];    
          Ex[i][j][k] = Ex[i][j][k] - BETA_EX(i,j,k) * SCALE_Jx( creal( Jsum ) , i );
          for( poleIdx = 0 ; poleIdx < debyeArray[block].medium->numPoles; poleIdx++ )
            debyeArray[block].Jpolx[ii][jj][kk][poleIdx] = debyeArray[block].medium->dalpha[poleIdx] * debyeArray[block].Jpolx[ii][jj][kk][poleIdx] + 
              debyeArray[block].medium->dbeta[poleIdx] * UNSCALE_Ex( Ex[i][j][k] - debyeArray[block].lastEx[ii][jj][kk] , i );
          debyeArray[block].lastEx[ii][jj][kk] = Ex[i][j][k];
        }
      
    /* Update Jpoly. */         
    field = EY;
    for( i = debyeArray[block].flim[field][XLO] , ii = 0 ; i <= debyeArray[block].flim[field][XHI] ; i++ , ii++ )
      for( j = debyeArray[block].flim[field][YLO] , jj = 0 ; j <= debyeArray[block].flim[field][YHI] ; j++ , jj++ )
        for( k = debyeArray[block].flim[field][ZLO] , kk = 0 ; k <= debyeArray[block].flim[field][ZHI] ; k++ , kk++ )
        {
          Jsum = 0.0 + I * 0.0;
          for( poleIdx = 0 ; poleIdx < debyeArray[block].medium->numPoles; poleIdx++ )
            Jsum = Jsum + ( 1 + debyeArray[block].medium->dalpha[poleIdx] ) * debyeArray[block].Jpoly[ii][jj][kk][poleIdx];    
          Ey[i][j][k] = Ey[i][j][k] - BETA_EY(i,j,k) * SCALE_Jy( creal( Jsum ) , j );  
          for( poleIdx = 0 ; poleIdx < debyeArray[block].medium->numPoles; poleIdx++ )
            debyeArray[block].Jpoly[ii][jj][kk][poleIdx] = debyeArray[block].medium->dalpha[poleIdx] * debyeArray[block].Jpoly[ii][jj][kk][poleIdx] + 
              debyeArray[block].medium->dbeta[poleIdx] * UNSCALE_Ey( Ey[i][j][k] - debyeArray[block].lastEy[ii][jj][kk] , j );
          debyeArray[block].lastEy[ii][jj][kk] = Ey[i][j][k]; 
        }
      
    /* Update Jpolz. */
    field = EZ; 
      for( i = debyeArray[block].flim[field][XLO] , ii = 0 ; i <= debyeArray[block].flim[field][XHI] ; i++ , ii++ )
        for( j = debyeArray[block].flim[field][YLO] , jj = 0 ; j <= debyeArray[block].flim[field][YHI] ; j++ , jj++ )
          for( k = debyeArray[block].flim[field][ZLO] , kk = 0 ; k <= debyeArray[block].flim[field][ZHI] ; k++ , kk++ )
          {
          Jsum = 0.0 + I * 0.0;
          for( poleIdx = 0 ; poleIdx < debyeArray[block].medium->numPoles; poleIdx++ )
            Jsum = Jsum + ( 1 + debyeArray[block].medium->dalpha[poleIdx] ) * debyeArray[block].Jpolz[ii][jj][kk][poleIdx];    
          Ez[i][j][k] = Ez[i][j][k] - BETA_EZ(i,j,k) * SCALE_Jz( creal( Jsum ) , k );  
          for( poleIdx = 0 ; poleIdx < debyeArray[block].medium->numPoles; poleIdx++ )
            debyeArray[block].Jpolz[ii][jj][kk][poleIdx] = debyeArray[block].medium->dalpha[poleIdx] * debyeArray[block].Jpolz[ii][jj][kk][poleIdx] + 
              debyeArray[block].medium->dbeta[poleIdx] * UNSCALE_Ez( Ez[i][j][k] - debyeArray[block].lastEz[ii][jj][kk] , k );
          debyeArray[block].lastEz[ii][jj][kk] = Ez[i][j][k]; 
          }

  }

  return;

}
