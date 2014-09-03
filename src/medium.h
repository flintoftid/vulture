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

#ifndef _MEDIA_H_
#define _MEDIA_H_

#include <complex.h>

#include "fdtd_types.h"
#include "utlist.h"
#include "uthash.h"

/* Index type used for counting and iterating over media. */
typedef unsigned int MediumIndex;
#define MAX_MEDIA UINT_MAX
//typedef unsigned char MediumIndex;
//#define MAX_MEDIA UCHAR_MAX


/* 
 * Medium types.
 * 
 * Medium types must begin at zero, be contigous and end with MT_UNDEFINED, which
 * *is not* included in the number NUM_MEDIUM_TYPES.
 */

#define NUM_MEDIUM_TYPES 4

typedef enum {

  /* The order and values of the explicitly assigned elements are significant! */
  MT_FREE_SPACE = 0,
  MT_PEC = 1,
  MT_SIMPLE = 2,
  MT_DEBYE,
  MT_UNDEFINED

} MediumType;

/* Medium type strings. */
extern char MEDIUM_TYPE_STR[NUM_MEDIUM_TYPES][TAG_SIZE];

/* 
 * Medium class. 
 */

typedef struct MediumItem_t {

  MediumIndex number;               // Medium number, assigned in order found.

  /* Parameters from the mesh. */
  char name[TAG_SIZE];              // Medium name.
  MediumType type;                  // Medium type.
  real eps_r;                       // Relative permittivity [-].
  real sigma;                       // Conductivity [S/m].
  real mu_r;                        // Relative permeability [-].
  int numPoles;                     // Number of Debye poles.
  double complex *residues;         // Debye residues.
  double complex *poles;            // Debye poles.
  char fileName[PATH_SIZE];         // Debye material parameter input file.

  /* Dervied parameters. */
  real alpha;                       // Simple media update coefficients.
  real beta;
  real gamma;
  double complex *dalpha;           // Debye update coefficients.
  double complex *dbeta;

  /* UT list/hash. */

  struct MediumItem_t *prev;
  struct MediumItem_t *next;

  UT_hash_handle hh;
  UT_hash_handle hhint;
  
} MediumItem;


/* 
 * Global data. 
 */

/* This has to be exposed in order to allow indexed media to pick up */
/* the medium coeffficients from the media list. */
extern MediumItem **mediumArray;

/*
 * Public method interfaces.
 */

bool parseMT( char *line );
void initMedia( void );
void deallocMedia( void );
void addMedium( char *name , MediumType type , real eps_r , real sigma , real mu_r , 
                int numPoles , double residues[3] , double poles[3] , char fileName[PATH_SIZE] );
void reportMedia( void );
bool isMedium( char *name , MediumIndex *number );
bool mediumTypeByName( char *name , MediumType *type );
MediumType getMediumType( MediumIndex number );
MediumItem * getMedium( MediumIndex number );
char *getMediumName( MediumIndex number );
void getSimpleMediumCoefficients( real *alpha , real *beta , real *gamma , MediumIndex medium );
void getSimpleMediumParameters( real *eps_r , real *sigma , real *mu_r , MediumIndex medium );
void calcCoeffFromParam( real *alpha , real *beta , real *gamma , double complex *dalpha , double complex *dbeta ,
                         real dt , real eps_r , real sigma , real mu_r , 
                         int numPoles , double complex residues[] , double complex poles[] );
void updateSimpleMedium( MediumIndex index , real eps_r , real sigma , real mu_r );
bool thereAreMedia( MediumType );

#endif
