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

#ifndef _VULTURE_H_
#define _VULTURE_H_

// These are now passed as compiler options by cmake.

/* If defined, all the fields are scaled by their respectivew edge lengths. */
//#define USE_SCALE_FIELDS

/* If defined, media at interfaces are averaged. */
//#define USE_AVERAGED_MEDIA

/* If defined indexed media are used. Incompatible with scaled fields and averaged media! */
//#define USE_INDEXED_MEDIA

#endif

