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

#ifndef _BOUNDING_BOX_H_
#define _BOUNDING_BOX_H_

#include "fdtd_types.h"

/* True if XLO <= XHI etc. */
inline bool bboxIsNormal( int bbox[6] )
{
  return ( bbox[XLO] <= bbox[XHI] && bbox[YLO] <= bbox[YHI] && bbox[ZLO] <= bbox[ZLO] );
}

/* True if (normal) bbox1 is within (normal) bbox2. */
inline bool bboxIsWithin( int bbox1[6] , int bbox2[6] )
{
  return ( bbox1[XLO] >= bbox2[XLO] && bbox1[XHI] <= bbox2[XHI] 
        && bbox1[YLO] >= bbox2[YLO] && bbox1[YHI] <= bbox2[YHI] 
        && bbox1[ZLO] >= bbox2[ZLO] && bbox1[ZHI] <= bbox2[ZHI] );
}

/* Offset the bounding box bbox by the lower limits of bounding box offset. */
inline void offsetBoundingBox( int offsetBbox[6] , int bbox[6] , int offset[6] )
{

  offsetBbox[XLO] = bbox[XLO] + offset[XLO];
  offsetBbox[YLO] = bbox[YLO] + offset[YLO];
  offsetBbox[ZLO] = bbox[ZLO] + offset[ZLO];
  offsetBbox[XHI] = bbox[XHI] + offset[XLO];
  offsetBbox[YHI] = bbox[YHI] + offset[YLO];
  offsetBbox[ZHI] = bbox[ZHI] + offset[ZLO];

  return;

}

/* Prototypes. */
BoundingBoxType bboxType( int bbox[6] );
void getFaceOfBoundingBox( int faceBox[6] , int bbox[6] , MeshFace face );
bool bboxIsElemental( int bbox[6] );
CoordAxis bboxDirection( int bbox[6] );
bool fieldIsParallelToBoundary( FieldComponent field , MeshFace boundary );
bool fieldIsInBoundary( FieldComponent field , MeshFace boundary );
bool fieldIsParallelToAxis( FieldComponent field , CoordAxis axis );
void setBoundingBoxFromNodes( int bbox[6] , int ilo , int ihi , int jlo , int jhi , int klo , int khi );
void setBoundingBoxBoundaryFlags( bool includeBoundary[6] , bool isXlo ,  bool isXhi , bool isYlo , bool isYhi , bool isZlo , bool isZhi );
FaceMask setFaceMaskFromString( char *maskStr );
bool isFaceActive( FaceMask mask , MeshFace face );
void faceMask2boolArray( bool includeBoundary[6] , FaceMask mask );

#endif

