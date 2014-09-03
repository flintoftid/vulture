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

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "bounding_box.h"

/* Ensure inlines are mentioned once. */
extern inline bool bboxIsNormal( int bbox[6] );
extern inline bool bboxIsWithin( int bbox1[6] , int bbox2[6] );
extern inline void offsetBoundingBox( int offsetBbox[6] , int bbox[6] , int offset[6] );

/* Determine type of normal bounding box. */
BoundingBoxType bboxType( int bbox[6] )
{
  int ix , iy , iz , mask;
  BoundingBoxType type[8] = { BB_VOLUME  , BB_SURFACE , BB_SURFACE , BB_LINE , BB_SURFACE , BB_LINE , BB_LINE , BB_POINT };

  ix = ( bbox[XLO] == bbox[XHI] );
  iy = ( bbox[YLO] == bbox[YHI] );
  iz = ( bbox[ZLO] == bbox[ZHI] );
  
  mask = ix + 2 * iy + 4 * iz;
    
  return type[mask];

}

/* Determine bounding box direction. */
CoordAxis bboxDirection( int bbox[6] )
{
  int ix , iy , iz , mask;
  CoordAxis type[8] = { CA_UNDEFINED , XDIR , YDIR , ZDIR , ZDIR , YDIR , XDIR , CA_UNDEFINED };

  ix = ( bbox[XLO] == bbox[XHI] );
  iy = ( bbox[YLO] == bbox[YHI] );
  iz = ( bbox[ZLO] == bbox[ZHI] );
  
  mask = ix + 2 * iy + 4 * iz;
    
  return type[mask];

}

/* Determine type of normal bounding box. */
void getFaceOfBoundingBox( int faceBox[6] , int bbox[6] , MeshFace face )
{

  switch( face )
  {
  case XLO:
    faceBox[XLO] = bbox[XLO];
    faceBox[XHI] = bbox[XLO];
    faceBox[YLO] = bbox[YLO];
    faceBox[YHI] = bbox[YHI];
    faceBox[ZLO] = bbox[ZLO];
    faceBox[ZHI] = bbox[ZHI];
    break;
  case XHI:
    faceBox[XLO] = bbox[XHI];
    faceBox[XHI] = bbox[XHI];
    faceBox[YLO] = bbox[YLO];
    faceBox[YHI] = bbox[YHI];
    faceBox[ZLO] = bbox[ZLO];
    faceBox[ZHI] = bbox[ZHI];
    break;
  case YLO:
    faceBox[XLO] = bbox[XLO];
    faceBox[XHI] = bbox[XHI];
    faceBox[YLO] = bbox[YLO];
    faceBox[YHI] = bbox[YLO];
    faceBox[ZLO] = bbox[ZLO];
    faceBox[ZHI] = bbox[ZHI];
    break;
  case YHI:
    faceBox[XLO] = bbox[XLO];
    faceBox[XHI] = bbox[XHI];
    faceBox[YLO] = bbox[YHI];
    faceBox[YHI] = bbox[YHI];
    faceBox[ZLO] = bbox[ZLO];
    faceBox[ZHI] = bbox[ZHI];
    break;
  case ZLO:
    faceBox[XLO] = bbox[XLO];
    faceBox[XHI] = bbox[XHI];
    faceBox[YLO] = bbox[YLO];
    faceBox[YHI] = bbox[YHI];
    faceBox[ZLO] = bbox[ZLO];
    faceBox[ZHI] = bbox[ZLO];
    break;
  case ZHI:
    faceBox[XLO] = bbox[XLO];
    faceBox[XHI] = bbox[XHI];
    faceBox[YLO] = bbox[YLO];
    faceBox[YHI] = bbox[YHI];
    faceBox[ZLO] = bbox[ZHI];
    faceBox[ZHI] = bbox[ZHI];
    break;
  default:
    assert( 0 );
    break;
  }

  return;

}

/* True if (normal) bbox is elemental. */
bool bboxIsElemental( int bbox[6] )
{

  int maxSize;

  maxSize = fmax( bbox[XHI] - bbox[XLO] , bbox[YHI] - bbox[YLO] );
  maxSize = fmax( maxSize , bbox[ZHI] - bbox[ZLO] );
  
  switch( bboxType( bbox ) )
  {
  case BB_POINT:
    if( maxSize == 0 )
      return true;
    break;
  case BB_LINE:
  case BB_SURFACE:
  case BB_VOLUME:
    if( maxSize == 1 )
      return true;
    return true;
    break;  
  default:
     break;
  }

  return false;
}

/* True if field component is parallel to mesh face. */
bool fieldIsInBoundary( FieldComponent field , MeshFace boundary )
{
  
  bool ans = false;

  switch( field )
  {
  case EX:
    if( boundary == XLO || boundary == XHI )
      ans = false;
    else 
      ans = true;
    break;
  case EY:
    if( boundary == YLO || boundary == YHI )
      ans = false;
    else 
      ans = true; 
    break;
  case EZ:
    if( boundary == ZLO || boundary == ZHI )
      ans = false;
    else 
      ans = true; 
    break;
  case HX:
    if( boundary == XLO || boundary == XHI )
      ans = true;
    else 
      ans = false; 
    break;
  case HY:
    if( boundary == YLO || boundary == YHI )
      ans = true;
    else 
      ans = false; 
    break;
  case HZ:
    if( boundary == ZLO || boundary == ZHI )
      ans = true;
    else 
      ans = false; 
    break;
  default:
    assert( 0 );
    break;
  }

  return ans;
  
}


/* True if field component is parallel to mesh face. */
bool fieldIsParallelToBoundary( FieldComponent field , MeshFace boundary )
{
  
  bool ans = false;

  switch( field )
  {
  case EX:
  case HX:
    if( boundary == XLO || boundary == XHI )
      ans = false;
    else 
      ans = true;
    break;
  case EY:
  case HY:
    if( boundary == YLO || boundary == YHI )
      ans = false;
    else 
      ans = true; 
    break;
  case EZ:
  case HZ:
    if( boundary == ZLO || boundary == ZHI )
      ans = false;
    else 
      ans = true; 
    break;
  default:
    assert( 0 );
    break;
  }

  return ans;
  
}

/* True if field component is parallel to direciton. */
bool fieldIsParallelToAxis( FieldComponent field , CoordAxis axis )
{

  bool ans = false;
  
  switch( field )
  {
  case EX:
  case HX:
    if( axis == XDIR )
      ans = true;
    else 
      ans = false;
    break;
  case EY:
  case HY:
    if( axis == YDIR )
      ans = true;
    else 
      ans = false;
    break;
  case EZ:
  case HZ:
    if( axis == ZDIR )
      ans = true;
    else 
      ans = false;
    break;
  default:
    assert( 0 );
    break;
  }

  return ans;

}

/* Set a bounding box from its nodes. */
void setBoundingBoxFromNodes( int bbox[6] , int ilo , int ihi , int jlo , int jhi , int klo , int khi )    
{

  bbox[XLO] = ilo;
  bbox[XHI] = ihi;
  bbox[YLO] = jlo;
  bbox[YHI] = jhi;
  bbox[ZLO] = klo;
  bbox[ZHI] = khi;

  return;
    
}
    
/* Set bounding box boundary inclusion flags. */
void setBoundingBoxBoundaryFlags( bool includeBoundary[6] , bool isXlo ,  bool isXhi , bool isYlo , bool isYhi , bool isZlo , bool isZhi )    
{
  
  includeBoundary[XLO] = isXlo;
  includeBoundary[XHI] = isXhi;
  includeBoundary[YLO] = isYlo;
  includeBoundary[YHI] = isYhi;
  includeBoundary[ZLO] = isZlo;
  includeBoundary[ZHI] = isZhi;

  return;
  
}

/* Set a face mask from a string. */ 
FaceMask setFaceMaskFromString( char *maskStr )
{
  
  FaceMask mask = 0;
  long lmask = 0;
  char *junk = NULL;

  lmask = strtol( maskStr , &junk , 2 );

  if( lmask < 0 || lmask > 63L || *junk )
    mask = FACE_MASK_ERROR;  
  else
    mask = (FaceMask)lmask;
  
  return mask;

}

/* Test if mask is active on given face. */
bool isFaceActive( FaceMask mask , MeshFace face )
{

  return (bool)( mask & FACE_MASKS[face] );

}

/* Convert face bit mask to bool activity array. */
void faceMask2boolArray( bool includeBoundary[6] , FaceMask mask )
{

  for( MeshFace face = XLO ; face <= ZHI ; face++ )
    includeBoundary[face] = (bool)( mask & FACE_MASKS[face] );
  
  return;
  
}
