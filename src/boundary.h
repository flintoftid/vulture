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

#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "fdtd_types.h"
#include "medium.h"
#include "utlist.h"
#include "uthash.h"

#ifdef WITH_SIBC
#include "filter.h"
#endif

/* Index type used for counting and iterating over boundary types. */
typedef unsigned long BoundaryIndex;
#define MAX_BOUNDARY ULONG_MAX

/* 
 * Boundary types.
 * 
 * Boundary types must begin at zero, be contigous and end with BT_UNDEFINED, which
 * *is not* included in the number NUM_BOUNDARY_TYPES.
 */

#define NUM_BOUNDARY_TYPES 7

/* Boundary types. */
typedef enum {

  /* The order and values are significant for the explicitly assigned types. */

  BT_PMC = 0,
  BT_PML = 1,
  BT_PEC = 2,
  BT_FREE_SPACE = 3,
  BT_PERIODIC,
  BT_MUR,
  BT_SIBC,
  BT_UNDEFINED
  
} BoundaryType;

/* 
 * Boundary class. 
 */

typedef struct BoundaryItem_t {

  BoundaryIndex number;           // Boundary number, assigned in order.

  /* Parameters from mesh. */
  char name[TAG_SIZE];            // Boundary name.
  BoundaryType type;              // Boundary type.
  int numLayers;                  // Number of layers (PML).
  int order;                      // Order of boundary (PML).
  real n_eff;                     // Effective refractive index of dominant modes (PML).
  real refCoeff;                  // Theoretical reflection coefficient (PML).
  real kmax;                      // Maximum permittivity (kappa) of PML (PML).
  char fileName[PATH_SIZE];       // Name of model file (SIBC).
  real S_TM[2][2];                 // Scattering matrix of TM mode. 
  real S_TE[2][2];                 // Scattering matrix of TE mode. 

  /* Derived parameters.*/
#ifdef WITH_SIBC
  yfPoleResidueM prm;             // Pole-residue model.
  yfRecConvM rcm;                 // Recursive convolution coefficients.
#endif

  /* UT list/hash. */

  struct BoundaryItem_t *prev;
  struct BoundaryItem_t *next;

  UT_hash_handle hh;             // Boundary name hash.
  UT_hash_handle hhint;          // Boundary number hash.

} BoundaryItem;

/*
 * Global data.
 */

/* Mapping of boundary types to medium types. */
extern MediumType BOUNDARY_MEDIUM_TYPE[NUM_BOUNDARY_TYPES];
extern char BOUNDARY_TYPE_STR[NUM_BOUNDARY_TYPES][11];

/*
 * Public method interfaces.
 */

bool parseBE( char *line );
bool parseBT( char *line );
void initBoundaries( void );
void deallocBoundaries( void );
void reportBoundaries( void );
bool isBoundary( char *name , BoundaryIndex *number );
void addBoundary( char *name , BoundaryType type , int numLayers , int order , 
                  real n_eff , real refCoeff , real kmax , char *fileName , real S_TM[2][2] , real S_TE[2][2] );
BoundaryItem * getBoundary( BoundaryIndex number );
BoundaryType getBoundaryType( BoundaryIndex number );
char *getBoundaryName( BoundaryIndex number );
int getBoundaryNumLayers( BoundaryIndex number );
real getBoundaryRefCoeff( BoundaryIndex number );
void setBoundaryNumLayers( BoundaryIndex number , int numLayers );
bool thereAreBoundaries( BoundaryType );
void getExternalBoundaryParams( BoundaryIndex number , int *order , real *n_eff , real *refCoeff , real *kmax );

#endif
