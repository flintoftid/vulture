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

#ifndef _SURFACE_H_
#define _SURFACE_H_

#include "fdtd_types.h"
#include "boundary.h"

/* Index type used for counting and iterating over surfaces. */
typedef unsigned long SurfaceIndex;
#define MAX_SURFACE ULONG_MAX

/* 
 * Surface class. 
 */

typedef struct SurfaceItem_t {

  /* Parameters from mesh. */
  int mbbox[6];                      // Bounding box on mesh.
  char boundaryName[TAG_SIZE];       // Boundary name.
  int orientation;                   // Orientation of boundary on surface.                  
  real angle;                        // Angle of boundary on surface.

  /* Derived parameters. */
  BoundaryIndex boundaryNumber;      // Boundary number.
  int gbbox[6];                      // Bounding box on grid.

  /* UT list. */

  struct SurfaceItem_t *prev;
  struct SurfaceItem_t *next;

} SurfaceItem;

/*
 * Public method interfaces.
 */

bool parseBR( char *line );
bool parseTB( char *line );
void initExternalSurfaceParameters( void );
void initExternalSurfaces( void );
void initInternalSurfaces( void );
void reportSurfaces( void );
void updateExternalSurfacesEfield( void );
void updateExternalSurfacesHfield( void );
void updateInternalSurfacesEfield( void );
void updateInternalSurfacesHfield( void );
void updateGhostEfield( void );
void updateGhostHfield( void );
void deallocExternalSurfaces( void );
void deallocInternalSurfaces( void );
void gnuplotExternalSurfaces( void );
void gnuplotInternalSurfaces( void );
void gmshExternalSurfaces( void );
void gmshInternalSurfaces( void );
void checkExternalSurfaces( void );
BoundaryType outerSurfaceType( MeshFace face );
int outerSurfaceNumLayers( MeshFace face );
real outerSurfaceReflectCoeff( MeshFace face );
void getOuterSurfaceParams( MeshFace face , int *order , real *n_eff , real *refCoeff , real *kmax );
void initExternalPecPmcSurfaces( void );
bool thereAreInternalSurfaces( BoundaryType type );
bool thereAreExternaSurfaces( BoundaryType type );
bool isPmcEdge( CoordAxis direction , int index );

#endif
