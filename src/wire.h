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

#ifndef _WIRE_H_
#define _WIRE_H_

#include "fdtd_types.h"
#include "medium.h"
#include "utlist.h"
#include "uthash.h"

/* Index type used for counting and iterating over wire types. */
typedef unsigned long WireIndex;
#define MAX_WIRE ULONG_MAX

/* 
 * Wire types.
 * 
 * Wire types must begin at zero, be contigous and end with TW_UNDEFINED, which
 * *is not* included in the number NUM_WIRE_TYPES.
 */

#define NUM_WIRE_TYPES 2

/* Wire types. */
typedef enum {
  
  TW_PEC,
  TW_FREE_SPACE,
  TW_UNDEFINED

} WireType;

/* Wire end types. */

#define NUM_WIRE_END_TYPES 3

typedef enum {

  WE_THRU,
  WE_END,
  WE_CORNER,
  WE_UNDEFINED

} WireEndType;

extern char WIRE_END_TYPE_STR[NUM_WIRE_END_TYPES][7];

/* 
 * Wire class. 
 */

typedef struct WireItem_t {

  WireIndex number;           // Wire number, assigned in order.

  /* Parameters from mesh. */
  char name[TAG_SIZE];            // Wire name.
  WireType type;                  // Wire type.
  real radius;                    // Wire radius.

  /* Derived parameters.*/

  /* UT list/hash. */

  struct WireItem_t *prev;
  struct WireItem_t *next;

  UT_hash_handle hh;             // Wire name hash.
  UT_hash_handle hhint;          // Wire number hash.

} WireItem;

/*
 * Global data.
 */

/* Mapping of boundary types to medium types. */
extern MediumType WIRE_MEDIUM_TYPE[NUM_WIRE_TYPES];
extern char WIRE_TYPE_STR[NUM_WIRE_TYPES][11];

/*
 * Public method interfaces.
 */

bool parseWT( char *line );
void initWires( void );
void deallocWires( void );
void reportWires( void );
bool isWire( char *name , WireIndex *number );
void addWire( char *name , WireType type , real radius );
WireItem * getWire( WireIndex number );
WireType getWireType( WireIndex number );
char *getWireName( WireIndex number );
bool thereAreWires( WireType );

#endif
