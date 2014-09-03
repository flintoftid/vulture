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

#ifndef FILTER_H
#define FILTER_H

#include <complex.h>
#include <stdbool.h>

/* Single filter pole-residue form. */
typedef struct yfPoleResidueType {

  int numPoles;   
  double asymp;
  double complex *residues;
  double complex *poles;
  
} yfPoleResidue;

/* Single RC filter. */
typedef struct yfRecConvFilterType {

  int numPoles;
  double asymp;
  double complex **q;

} yfRecConv;

/* Single RC filter state. */
typedef struct yfRecConvStateType {

  double complex *zeta;
  double old;

} yfRecConvState;

/* Pole-residue coefficient matrix. */
typedef struct yfPoleResidueTypeM {

  int m;
  int n;

  yfPoleResidue **pr;

} yfPoleResidueM;

/* Recursive convolution coefficient matrix. */
typedef struct yfRecConvTypeM {

  int m;
  int n;

  yfRecConv **rc;

} yfRecConvM;

/* Recursive convolution state matrix. */
typedef struct yfRecConvStateTypeM {

  int m;
  int n;

  yfRecConvState **rc_s;

} yfRecConvStateM;

/* Pole residue filters. */
yfPoleResidue yfAllocPoleResidue( int numPoles , double asymp , double complex poles[] , double complex residues[] );
void yfDeallocPoleResidue( yfPoleResidue pr );
yfPoleResidue yfReadPoleResidue( char fileName[] );
yfPoleResidueM yfAllocPoleResidueM( int m , int n );
void yfDeallocPoleResidueM( yfPoleResidueM prm );
yfPoleResidueM yfReadPoleResidueM( char prmFileName[] );
void yfPrintPoleResidue( yfPoleResidue pr );
void yfPrintPoleResidueM( yfPoleResidueM prm );

/* RC filters. */
yfRecConv yfAllocRecConv( int numPoles , double asymp , double complex **q );
void yfDeallocRecConv( yfRecConv rc );
yfRecConvState yfAllocRecConvState( yfRecConv rc );
void yfDeallocRecConvState( yfRecConvState rc_s );
double yfRecConvStep( yfRecConv rc , yfRecConvState *rc_s , double x );
void yfRecConvFiltSeq( yfRecConv rc , int length , double *x , double *y );
yfRecConvM yfAllocRecConvM( int m , int n );
void yfDeallocRecConvM( yfRecConvM rcm );
yfRecConvStateM yfAllocRecConvStateM( yfRecConvM rcm );
void yfDeallocRecConvStateM( yfRecConvStateM rcm_s );
void yfPrintRecConvM( yfRecConvM rcm );
void yfPrintRecConv( yfRecConv rc );

/* Filter conversions. */
yfRecConv yfPoleResidue2RecConv( yfPoleResidue pr , double T );
yfRecConvM yfPoleResidueM2RecConvM( yfPoleResidueM prm , double dt );

#endif
