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

#include "fdtd_types.h"

char AXIS[3][2]  = { "X" , "Y" , "Z" };
char DIR[6][5]   = { "XPOS" , "XNEG" , "YPOS" , "YNEG" , "ZPOS" , "ZNEG" };
char FACE[6][4]  = { "XLO" , "XHI" , "YLO" , "YHI" , "ZLO" , "ZHI" };
char FIELD[6][3] = { "EX" , "EY" , "EZ" , "HX" , "HY" , "HZ" };
char BOOL[2][6]  = { "FALSE" , "TRUE" };
char BBOX_STR[5][8]= { "INVALID" , "POINT" , "LINE" , "SURFACE" , "VOLUME" };
FaceMask FACE_MASKS[6] = { 0x20 , 0x10 , 0x8 , 0x4 , 0x2 , 0x1 };