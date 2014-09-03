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
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <assert.h>

#include "filter.h"
#include "message.h"
#include "alloc_array.h"


/* Allocate pole residue filter. */
yfPoleResidue yfAllocPoleResidue( int numPoles , double asymp , double complex poles[] , double complex residues[] )
{

   yfPoleResidue pr;
  
   pr.numPoles = numPoles;
   pr.asymp = asymp;
   pr.residues = (double complex *) malloc( numPoles * sizeof( double complex ) );
   pr.poles = (double complex *) malloc( numPoles * sizeof( double complex ) );
 
   if( residues != NULL )
     for( int i = 0 ; i < numPoles ; i++ )
       pr.residues[i] = residues[i];
 
   if( poles != NULL )
     for( int i = 0 ; i < numPoles ; i++ )
       pr.poles[i] = poles[i];
     
   return pr;

}

/* Deallocate pole residue filter. */
void yfDeallocPoleResidue( yfPoleResidue pr )
{

  free( pr.residues );
  free( pr.poles );

  return;

}

/* Read a pole residue filter from ASCII pr file. */
yfPoleResidue yfReadPoleResidue( char fileName[] )
{

  int k;
  FILE *fp;
  int numPoles;
  int numScanned;
  double asymp;
  double pole_r;
  double pole_i;
  double residue_r;
  double residue_i;
  yfPoleResidue pr;

  fp = fopen( fileName , "r" );

  numScanned = fscanf( fp , "%d" , &numPoles );
  if( numScanned != 1 )
    message( MSG_ERROR , 0 , "error reading number of poles from %s.\n" , fileName );

  numScanned = fscanf( fp , "%lf" , &asymp );
  if( numScanned != 1 )
    message( MSG_ERROR , 0 , "error reading asymptote from %s.\n" , fileName );

  pr = yfAllocPoleResidue( numPoles , asymp , NULL , NULL );

  for( k = 0 ; k < pr.numPoles ; k++ )
  {
    numScanned = fscanf( fp , "%lf %lf %lf %lf" , &pole_r , &pole_i , &residue_r , &residue_i );
    pr.poles[k] = pole_r + I * pole_i;
    pr.residues[k] = residue_r + I * residue_i;
  }

  fclose( fp );

  return pr;

}

/* Allocate pole-residue filter matrix. */
yfPoleResidueM yfAllocPoleResidueM( int m , int n )
{

  yfPoleResidueM prm;
  unsigned long bytes;
  
  prm.m = m;
  prm.n = n;

  prm.pr = (yfPoleResidue **) allocArray( &bytes, sizeof( yfPoleResidue ) , 2 , m , n );

  return prm;

}

/* Deallocate pole-residue filter matrix. */
void yfDeallocPoleResidueM( yfPoleResidueM prm )
{

  for( int row = 0 ; row < prm.m ; row++ )
    for( int col = 0 ; col < prm.n ; col++ )
      yfDeallocPoleResidue( prm.pr[row][col] );

  deallocArray( prm.pr , 2 , prm.m , prm.n );

  return;

}

/* Read pole-residue filter matrix from disk. */
yfPoleResidueM yfReadPoleResidueM( char prmFileName[] )
{

  int m;
  int n;
  FILE *fp;
  int numPoles;
  int numScanned;
  double asymp;
  double pole_r;
  double pole_i;
  double residue_r;
  double residue_i;
  yfPoleResidueM prm;

  fp = fopen( prmFileName , "r" );
  if( !fp )
    message( MSG_ERROR , 0 , "cannot open pole-residue matrix file %s.\n" , prmFileName );

  numScanned = fscanf( fp , "%d %d" , &m , &n );
  if( numScanned != 2 )
    message( MSG_ERROR , 0 , "error reading extents from %s.\n" , prmFileName );

  prm = yfAllocPoleResidueM( m , n );
 
  for( int row = 0 ; row < m ; row++ )
    for( int col = 0 ; col < n ; col++ )
    {

      numScanned = fscanf( fp , "%d" , &numPoles );
      numScanned = fscanf( fp , "%lf" , &asymp );
      prm.pr[row][col] = yfAllocPoleResidue( numPoles , asymp , NULL , NULL );

      for( int k = 0 ; k < numPoles ; k++ )
      {
        numScanned = fscanf( fp , "%lf %lf %lf %lf" , &pole_r , &pole_i , &residue_r , &residue_i );
        prm.pr[row][col].poles[k] = pole_r + I * pole_i;
        prm.pr[row][col].residues[k] = residue_r + I * residue_i;
      }

    }

  fclose( fp );

  return prm;

}

/* Allocate a RC filter. */
yfRecConv yfAllocRecConv( int numPoles , double asymp , double complex **q )
{

   yfRecConv rc;
   unsigned long bytes;
   
   rc.numPoles = numPoles;
   rc.asymp = asymp;

   rc.q = (double complex **) allocArray( &bytes , sizeof( double complex ) , 2 , 3 , numPoles );

   if( q != NULL )
     for( int i = 0 ; i < 3 ; i++ )
        for( int j = 0 ; j < numPoles ; j++ )
           rc.q[i][j] = q[i][j];

   return rc;
  
}

/* Deallocate a RC filter. */
void yfDeallocRecConv( yfRecConv rc )
{

  deallocArray( rc.q , 2 , 3 , rc.numPoles );

  return;
  
}

/* Allocate recursive convolution filter matrix. */
yfRecConvM yfAllocRecConvM( int m , int n )
{

  yfRecConvM rcm;
  unsigned long bytes;
   
  rcm.m = m;
  rcm.n = n;

  rcm.rc = (yfRecConv **) allocArray( &bytes , sizeof( yfRecConv ) , 2 , m , n ); 

  return rcm;

}

/* Deallocate recursive convolution filter matrix. */
void yfDeallocRecConvM( yfRecConvM rcm )
{

  for( int row = 0 ; row < rcm.m ; row++ )
    for( int col = 0 ; col < rcm.n ; col++ )
      yfDeallocRecConv( rcm.rc[row][col] );

  deallocArray( rcm.rc , 2 , rcm.m , rcm.n );

  return;

}

/* Allocate RC filter state. */
yfRecConvState yfAllocRecConvState( yfRecConv rc )
{
   
  yfRecConvState rc_s;

  rc_s.zeta = (double complex *) calloc( rc.numPoles , sizeof( double complex ) );
  rc_s.old = 0.0;
 
  return rc_s;
  
}

/* Deallocate RC filter state. */
void yfDeallocRecConvState( yfRecConvState rc_s )
{

  free( rc_s.zeta );

  return;
  
}

/* Allocate recursive convolution state matrix. */
yfRecConvStateM yfAllocRecConvStateM( yfRecConvM rcm )
{

  yfRecConvStateM rcm_s;
  unsigned long bytes;
  
  rcm_s.m =rcm.m;
  rcm_s.n =rcm.n;

  rcm_s.rc_s = (yfRecConvState **) allocArray( &bytes , sizeof( yfRecConvState ) , 2 , rcm.m , rcm.n );

  for( int row = 0 ; row < rcm.m ; row++ )
    for( int col = 0 ; col < rcm.n ; col++ )
      rcm_s.rc_s[row][col] = yfAllocRecConvState( rcm.rc[row][col] );

  return rcm_s;

}
 
/* Deallocate recursive convolution state matrix. */
void yfDeallocRecConvStateM( yfRecConvStateM rcm_s )
{

  for( int row = 0 ; row < rcm_s.m ; row++ )
    for( int col = 0 ; col < rcm_s.n ; col++ )
      yfDeallocRecConvState( rcm_s.rc_s[row][col] );

  deallocArray( rcm_s.rc_s , 2 , rcm_s.m , rcm_s.n );

  return;

}

/* Step RC filter one time-step. */
double yfRecConvStep( yfRecConv rc , yfRecConvState *rc_s , double x )
{

  int k;
  double y;

  y = rc.asymp * x;
  
  for( k = 0 ; k < rc.numPoles ; k++ )
  {
    rc_s->zeta[k] = rc.q[0][k] * rc_s->zeta[k] + rc.q[1][k] * rc_s->old + rc.q[2][k] * x;
    y = y + rc_s->zeta[k];
  } 

  rc_s->old = x;

  return y;

}

/* Filter a time-series using a RC filter. */
void yfRecConvFiltSeq( yfRecConv rc , int numPoints , double *x , double *y )
{

  yfRecConvState rc_s;

  rc_s = yfAllocRecConvState( rc );

  for( int k = 0 ; k < numPoints ; k++ )
    y[k] = yfRecConvStep( rc , &rc_s , x[k] );

  yfDeallocRecConvState( rc_s );

  return;

}

/* Determine RC coefficients from pole residue expansion. */
yfRecConv yfPoleResidue2RecConv( yfPoleResidue pr , double T )
{

  yfRecConv rc;
  int k;
  double complex alpha;
  double complex beta;

  rc = yfAllocRecConv( pr.numPoles , pr.asymp , NULL );

  for( k = 0 ; k < rc.numPoles ; k++ ) 
  {
      alpha = pr.residues[k] / pr.poles[k];
      beta  = pr.poles[k] * T;
      rc.q[0][k] = cexp( beta );
      rc.q[1][k] = alpha / beta * ( 1.0 + ( beta - 1.0 ) * cexp( beta ) );
      rc.q[2][k] = alpha / beta * ( cexp( beta ) - beta - 1.0 );
  }

  return rc;

}

/* Convert pole-residue coefficent matrix to recursive convolution coefficient matrix. */   
yfRecConvM yfPoleResidueM2RecConvM( yfPoleResidueM prm , double dt )
{

  yfRecConvM rcm;

  rcm = yfAllocRecConvM( prm.m , prm.n );

  rcm.m = prm.m;
  rcm.n = prm.n; 

  for( int row = 0 ; row < prm.m ; row++ )
    for( int col = 0 ; col < prm.n ; col++ )
       rcm.rc[row][col] = yfPoleResidue2RecConv( prm.pr[row][col] , dt );

  return rcm;

}

/* Print pole-residue matrix. */
void yfPrintPoleResidueM( yfPoleResidueM prm )
{

  printf( "\nPRM: %dx%d\n" , prm.m , prm.n );
  for( int row = 0 ; row < prm.m ; row++ )
    for( int col = 0 ; col < prm.n ; col++ )
    {
      printf( "pr(%d,%d):\n" , row , col );
      yfPrintPoleResidue( prm.pr[row][col] ); 
    }
  printf( "\n" );

  return;

}

/* Print pole-residue. */
void yfPrintPoleResidue( yfPoleResidue pr )
{

  printf( "PR: asymp = %e, numPoles = %d\n" , pr.asymp , pr.numPoles );
  for( int p = 0 ; p < pr.numPoles ; p++ )
  {
    printf( "[%d]: %16.8e %16.8e %16.8e %16.8e\n" , p , creal( pr.poles[p] ) , cimag( pr.poles[p] ) , creal( pr.residues[p] ) , cimag( pr.residues[p] ) ); 
  }

  return;

}

/* Print recursive convolution matrix. */
void yfPrintRecConvM( yfRecConvM rcm )
{

  printf( "\nRCM: %dx%d\n" , rcm.m , rcm.n );
  for( int row = 0 ; row < rcm.m ; row++ )
    for( int col = 0 ; col < rcm.n ; col++ )
    {
      printf( "rc(%d,%d):\n" , row , col );
      yfPrintRecConv( rcm.rc[row][col] ); 
    }
  printf( "\n" );

  return;

}

/* Print recursive convolution. */
void yfPrintRecConv( yfRecConv rc )
{

  printf( "RC: asymp = %e, numPoles = %d\n" , rc.asymp , rc.numPoles );
  for( int p = 0 ; p < rc.numPoles ; p++ )
  {
    printf( "[%d]: %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n" , p , creal( rc.q[0][p] ) , cimag( rc.q[0][p] ) , 
                                                                      creal( rc.q[1][p] ) , cimag( rc.q[1][p] ) , 
                                                                      creal( rc.q[2][p] ) , cimag( rc.q[2][p] ) ); 
  }

  return;

}

