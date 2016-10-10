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

#ifndef _GMSH_H_
#define _GMSH_H_

#include <stdbool.h>
#include <stdio.h>

#include "fdtd_types.h"
 
/* Length of physical name buffers. */
#define GMSH_NAME_LENGTH 48

/* Element types. */
#define NUM_ELEMENT_TYPES 4

typedef enum {

  ET_NODE1,
  ET_BAR2,
  ET_QUAD4,
  ET_HEX8,
  ET_UNDEFINED

} ElementType;

void gmshMesh( bool isPhysUnits , bool isExternalSurfaces );
unsigned long gmshGetEntityNumber( void );
void gmshAddEntity( unsigned long entityNumber , BoundingBoxType typeToAdd , char *name , int mbox[6] , int step[3] );

#endif
