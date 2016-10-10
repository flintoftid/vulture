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

#ifndef _FDTD_TYPES_H_
#define _FDTD_TYPES_H_

#include <stdbool.h>
#include <float.h>
#include <limits.h>

/* Maximum length of filename. */
#define FMAX 1024

/* Size of real type.*/
typedef float real;
#define REAL_MAX     FLT_MAX
#define REAL_EPSILON FLT_EPSILON

/* Maxium size of tag strings. */
#define TAG_SIZE 32

/* Maxium size of path name. */
#define PATH_SIZE 1024

/* Maximum size of comment string. */
#define COMMENT_SIZE 1024

/*
 * Mesh entities.
 */

typedef enum coord_axis {

  XDIR = 0,
  YDIR = 1,
  ZDIR = 2,
  CA_UNDEFINED

} CoordAxis;

extern char AXIS[3][2];

typedef enum coord_direction {

  XPOS = 0,
  XNEG = 1,
  YPOS = 2,
  YNEG = 3,
  ZPOS = 4,
  ZNEG = 5

} CoordDirection;

extern char DIR[6][5];

typedef struct coord {

  int x;
  int y;
  int z;

} Coord;

typedef enum mesh_faces {

  /* The order and values are significant! */
  XLO = 0,
  XHI = 1,
  YLO = 2,
  YHI = 3,
  ZLO = 4,
  ZHI = 5

} MeshFace;

extern char FACE[6][4];

typedef enum field_component {

  /* The order and values are significant! */
  EX = 0,
  EY = 1,
  EZ = 2,
  HX = 3,
  HY = 4,
  HZ = 5

} FieldComponent;

extern char FIELD[6][3];

typedef enum {

  BB_INVALID = 0,
  BB_POINT,
  BB_LINE,
  BB_SURFACE,
  BB_VOLUME

} BoundingBoxType;

extern char BBOX_STR[5][8];

extern char BOOL[2][6];

typedef unsigned char FaceMask;
#define FACE_MASK_ERROR 0xFF
#define FACE_MASK_ALL   0x3F
extern FaceMask FACE_MASKS[6];

#endif

