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

#ifndef _LINE_H_
#define _LINE_H_

#include "fdtd_types.h"
#include "wire.h"

/* Index type used for counting and iterating over lines. */
typedef unsigned long LineIndex;
#define MAX_LINE ULONG_MAX

/* 
 * Line class. 
 */

typedef struct LineItem_t {

  /* Parameters from mesh. */
  int mbbox[6];                      // Bounding box on mesh.
  char wireName[TAG_SIZE];           // Wire name.
  WireEndType lowEndType;            // End type at low end.
  WireEndType highEndType;           // End type at high end.
  
  /* Derived parameters. */
  WireIndex wireNumber;              // Wire number.
  int gbbox[6];                      // Bounding box on grid.

  /* UT list. */

  struct LineItem_t *prev;
  struct LineItem_t *next;

} LineItem;

/*
 * Public method interfaces.
 */

bool parseTW( char *line );
void initLines( void );
void reportLines( void );
void updateLinesEfield( void );
void updateLinesHfield( void );
void deallocLines( void );
bool thereAreLines( WireType type );
void gmshLines( void );
void gnuplotLines( void );

#endif

