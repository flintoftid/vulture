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

#ifndef _BLOCK_H_
#define _BLOCK_H_

#include "fdtd_types.h"
#include "medium.h"

/* Index type used for counting and iterating over block. */
typedef unsigned long BlockIndex;
#define MAX_BLOCK ULONG_MAX

/* 
 * Block class. 
 */

typedef struct BlockItem_t {

  /* Parameters from mesh. */

  int mbbox[6];                   // Bounding box on mesh.
  char mediumName[TAG_SIZE];      // Medium name.
  FaceMask mask;                  // Face inclusion mask.

  /* Dervied parameters. */
  MediumIndex mediumNumber;       // Medium number.
 
  /* UT list. */

  struct BlockItem_t *prev;
  struct BlockItem_t *next;

} BlockItem;

/*
 * Public method interfaces.
 */

bool parseMB( char *line );
void initBlocks( void );
void deallocBlocks( void );
void addBlock( int mbbox[6] , char *mediumName , MediumIndex mediumNumber , FaceMask mask );
void gnuplotBlocks( void );
void gmshBlocks( void );
void reportBlocks( void );
bool thereAreBlocks( MediumType type );
void updateBlocksEfield( void );
void updateBlocksHfield( void );

#endif
