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

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "gnuplot.h"
#include "bounding_box.h"
#include "message.h"
#include "source.h"
#include "planewave.h"
#include "observer.h"
#include "block.h"
#include "line.h"
#include "grid.h"
#include "mesh.h"
#include "surface.h"

static bool isPhysicalUnits = true;

/* Private functions. */
void gnuplotScript( bool isExternalSurfaces );
void gnuplotPoint( FILE *outputFile , int bbox[6] );
void gnuplotLine( FILE *outputFile , int bbox[6] , CoordAxis direction );
void gnuplotSurface( FILE *outputFile , int bbox[6] , CoordAxis direction );
void gnuplotVolume( FILE *outputFile , int bbox[6] );
void bboxInRealUnits( real physbbox[6] , int bbox[6] );
real indexInRealUnits( int index , CoordAxis dir );
bool equalBoundaries();

/* Render mesh to gnuplot compatible files. */
void gnuplotMesh( bool isPhysUnits , bool isExternalSurfaces )
{

  isPhysicalUnits = isPhysUnits;

  gnuplotExternalSurfaces();
  gnuplotInternalSurfaces();
  gnuplotBlocks();
  gnuplotLines();
  gnuplotSources();
  gnuplotPlaneWaves();
  gnuplotObservers(); 
  gnuplotGridLines();
  gnuplotScript( isExternalSurfaces );
  
  return;

}

/* Write gnuplot script. */
void gnuplotScript( bool isExternalSurfaces )
{

  char scriptFileName[] = "mesh.gnp";
  FILE *scriptFile;
  char unitStr[2][2] = { "m" , "-" };
  int units = 0;

  scriptFile = fopen( scriptFileName , "w" );
  if( !scriptFile )
    message( MSG_ERROR , 0 , "*** Error: Failed to open script output file %s\n" , scriptFileName );
 
  if( isPhysicalUnits )  
    units = 0;
  else
    units = 1;

  printf("%i %i %i %i %i %i surface\n", mbox[XLO], mbox[XHI], mbox[YLO], mbox[YHI], mbox[ZLO], mbox[ZHI]);

  fprintf( scriptFile , "set term push\n" );
  fprintf( scriptFile , "set term post eps enhanced color \"Helvetica\" 16\n" );
  fprintf( scriptFile , "set output 'mesh.eps'\n" );
  fprintf( scriptFile , "# Grid lines.\n" );
  fprintf( scriptFile , "set style line  1 lt  0 lc rgb \"#BBBBBB\" lw 1\n" );
  fprintf( scriptFile , "# External surfaces.\n" );
  fprintf( scriptFile , "set style line  2 lt  3 lc rgb \"#FF8C00\" lw 2\n" );
  fprintf( scriptFile , "# Internal surface.\n" );
  fprintf( scriptFile , "set style line  3 lt  1 lc rgb \"#FF0000\" lw 2\n" );
  fprintf( scriptFile , "# Blocks.\n" );
  fprintf( scriptFile , "set style line  4 lt  1 lc rgb \"#0000FF\" lw 2\n" );
  fprintf( scriptFile , "# Sources.\n" );
  fprintf( scriptFile , "set style line  5 lt  5 lc rgb \"#FF00FF\" lw 2\n" );
  fprintf( scriptFile , "# Obervers.\n" );
  fprintf( scriptFile , "set style line  6 pt  6 lc rgb \"#000000\" ps 0.2\n" );
  if ( equalBoundaries() )
    fprintf(scriptFile , "set view equal xyz\n");

  fprintf( scriptFile , "\n" );
  fprintf( scriptFile , "set ticslevel 0\n" );
  fprintf( scriptFile , "set xlabel 'x [%s]'\n" , unitStr[units] );
  fprintf( scriptFile , "set ylabel 'y [%s]'\n" , unitStr[units] );
  fprintf( scriptFile , "set zlabel 'z [%s]'\n" , unitStr[units] );
  fprintf( scriptFile , "set title '%s'\n" , getCommentReference() );
  fprintf( scriptFile , "splot 'gnuplot-lines.dat'    ti ''                  w l   ls 1 , \\\n" );
  if( isExternalSurfaces ) 
    fprintf( scriptFile , "      'gnuplot-external.dat' ti 'External surfaces' w l   ls 2 , \\\n" );
  else
    fprintf( scriptFile , "#      'gnuplot-external.dat' ti 'External surfaces' w l   ls 2 , \\\n" );    
  if( thereAreInternalSurfaces( BT_UNDEFINED ) )
    fprintf( scriptFile , "      'gnuplot-surface.dat'  ti 'Internal surfaces' w l   ls 3 , \\\n" );
  if( thereAreBlocks( MT_UNDEFINED ) )
    fprintf( scriptFile , "      'gnuplot-block.dat'    ti 'Blocks'            w l   ls 4 , \\\n" );
  if( thereAreLines( TW_UNDEFINED ) )
    fprintf( scriptFile , "      'gnuplot-wires.dat'    ti 'Wires'            w l   ls 7 , \\\n" );
  if( thereAreSources( ST_UNDEFINED ) )
    fprintf( scriptFile , "      'gnuplot-source.dat'      ti 'Sources'           w l   ls 5 , \\\n" );
  if( thereArePlaneWaves() )
    fprintf( scriptFile , "      'gnuplot-planewave.dat'   ti 'Plane waves'       w l   ls 5 , \\\n" );
  if( thereAreObservers() )
    fprintf( scriptFile , "      'gnuplot-observer.dat' ti 'Observers'         w poi ls 6\n" );
  fprintf( scriptFile , "set output\n" );
  fprintf( scriptFile , "set term pop\n" );
  fprintf( scriptFile , "replot\n");

  fclose( scriptFile );

  return;

}

/* Define problem size and irregular spacing to prevent mesh generation. */
void gnuplotProblemSize( FILE *outputFile , int ibbox[6] )
{

  real bbox[6];

  bboxInRealUnits( bbox , ibbox );

  fprintf( outputFile , "%g %g %g\n" , bbox[XLO] , bbox[YLO] , bbox[ZLO] );
  fprintf( outputFile , "%g %g %g\n" , bbox[XLO] , bbox[YLO] , bbox[ZLO] );
  fprintf( outputFile , "%g %g %g\n\n", 0.99995 * bbox[XLO] , 0.99995 * bbox[YLO] , 0.99995 * bbox[ZLO] );

  fprintf( outputFile , "%g %g %g\n" , bbox[XHI] , bbox[YHI] , bbox[ZHI] );
  fprintf( outputFile , "%g %g %g\n\n", 1.00005 * bbox[XHI] , 1.00005 * bbox[YHI] , 1.00005 * bbox[ZHI] );

  return;

}

/* Render a bounding box. */
void gnuplotBoundingBox( FILE *outputFile , int bbox[6] )
{

  switch( bboxType( bbox ) )
  {
  case BB_POINT:
    gnuplotPoint( outputFile , bbox );
    break;
  case BB_LINE:
    gnuplotLine( outputFile , bbox , bboxDirection( bbox ) );
    break;
  case BB_SURFACE:
    gnuplotSurface( outputFile , bbox , bboxDirection( bbox ) );
    break;
  case BB_VOLUME:
    gnuplotVolume( outputFile , bbox );
    break;
  default:
    break;
  }

  return;

}

/* Render bounding box nodes as points. */
void gnuplotBoundingBoxNodes( FILE *outputFile , int ibbox[6] , int step[3] )
{

  real ri , rj , rk;

  for( int k = ibbox[ZLO] ; k <= ibbox[ZHI] ; k += step[ZDIR] )
    for( int j = ibbox[YLO] ; j <= ibbox[YHI] ; j += step[YDIR] ) 
      for( int i = ibbox[XLO] ; i <= ibbox[XHI] ; i += step[XDIR] )
        {
           ri = indexInRealUnits( i , XDIR );
           rj = indexInRealUnits( j , YDIR );
           rk = indexInRealUnits( k , ZDIR );
           fprintf( outputFile , "%g %g %g\n\n" , ri , rj , rk );
        }

  return;

}

/* Render a point. */
void gnuplotPoint( FILE *outputFile , int ibbox[6] )
{

  real bbox[6];

  bboxInRealUnits( bbox , ibbox );

  fprintf( outputFile , "%g %g %g\n\n" , bbox[XLO] , bbox[YLO] , bbox[ZLO] );
      
  return;

}

/* Render a line */
void gnuplotLine( FILE *outputFile , int ibbox[6] , CoordAxis direction )
{

  real bbox[6];

  bboxInRealUnits( bbox , ibbox );

  switch( direction )
  {
  case XDIR:
    fprintf( outputFile , "%g %g %g\n" , bbox[XLO] , bbox[YLO] , bbox[ZLO] );
    fprintf( outputFile , "%g %g %g\n\n" , bbox[XHI] , bbox[YLO] , bbox[ZLO] );
    break;
  case YDIR:
    fprintf( outputFile , "%g %g %g\n" , bbox[XLO] , bbox[YLO] , bbox[ZLO] );
    fprintf( outputFile , "%g %g %g\n\n" , bbox[XLO] , bbox[YHI] , bbox[ZLO] );
    break;
  case ZDIR:
    fprintf( outputFile , "%g %g %g\n" , bbox[XLO] , bbox[YLO] , bbox[ZLO] );
    fprintf( outputFile , "%g %g %g\n\n" , bbox[XLO] , bbox[YLO] , bbox[ZHI] );
    break;
  default:
    assert( 0 );
    break;
  }
        
  return;

}

/* Render a surface. */
void gnuplotSurface( FILE *outputFile , int ibbox[6] , CoordAxis direction )
{

  real bbox[6];
  real ri , rj , rk;

  bboxInRealUnits( bbox , ibbox );

  switch( direction )
  {
  case XDIR:
    for( int i = ibbox[XLO]  ; i <= ibbox[XHI] ; i++ )
    {
      ri = indexInRealUnits( i , XDIR );
      fprintf( outputFile , "%g %g %g\n" , ri , bbox[YLO] , bbox[ZLO] );
      fprintf( outputFile , "%g %g %g\n" , ri , bbox[YLO] , bbox[ZHI] );
      fprintf( outputFile , "%g %g %g\n" , ri , bbox[YHI] , bbox[ZHI] );
      fprintf( outputFile , "%g %g %g\n" , ri , bbox[YHI] , bbox[ZLO] );
      fprintf( outputFile , "%g %g %g\n\n" , ri , bbox[YLO] , bbox[ZLO] );
    }
    break;
  case YDIR:
    for( int j = ibbox[YLO] ; j <= ibbox[YHI] ; j++ )
    {
      rj = indexInRealUnits( j , YDIR );
      fprintf( outputFile , "%g %g %g\n" , bbox[XLO] , rj , bbox[ZLO] );
      fprintf( outputFile , "%g %g %g\n" , bbox[XLO] , rj , bbox[ZHI] );
      fprintf( outputFile , "%g %g %g\n" , bbox[XHI] , rj , bbox[ZHI] );
      fprintf( outputFile , "%g %g %g\n" , bbox[XHI] , rj , bbox[ZLO] );
      fprintf( outputFile , "%g %g %g\n\n" , bbox[XLO], rj , bbox[ZLO] );
    }
    break;
  case ZDIR:
    for( int k = ibbox[ZLO] ; k <= ibbox[ZHI] ; k++ )
    {
      rk = indexInRealUnits( k , ZDIR );
      fprintf( outputFile , "%g %g %g\n" , bbox[XLO] , bbox[YLO] , rk );
      fprintf( outputFile , "%g %g %g\n" , bbox[XLO] , bbox[YHI] , rk );
      fprintf( outputFile , "%g %g %g\n" , bbox[XHI] , bbox[YHI] , rk );
      fprintf( outputFile , "%g %g %g\n" , bbox[XHI] , bbox[YLO] , rk );
      fprintf( outputFile , "%g %g %g\n\n" , bbox[XLO] , bbox[YLO] , rk );
    }
    break;	
  default:
    assert( 0 );
    break;
  }

  return;

}

/* Render a volume. */
void gnuplotVolume( FILE *outputFile , int ibbox[6] )
{

  real bbox[6];

  bboxInRealUnits( bbox , ibbox );

  /* ZLO-1 face. */
  fprintf( outputFile , "%g %g %g\n", bbox[XLO] , bbox[YLO] , bbox[ZLO] );
  fprintf( outputFile , "%g %g %g\n", bbox[XLO] , bbox[YHI] , bbox[ZLO] );
  fprintf( outputFile , "%g %g %g\n", bbox[XHI] , bbox[YHI] , bbox[ZLO] );
  fprintf( outputFile , "%g %g %g\n", bbox[XHI] , bbox[YLO] , bbox[ZLO] );
	
  /* XLO-1 face. */
  fprintf( outputFile , "%g %g %g\n", bbox[XLO] , bbox[YLO] , bbox[ZLO] );
  fprintf( outputFile , "%g %g %g\n", bbox[XLO] , bbox[YHI] , bbox[ZLO] );
  fprintf( outputFile , "%g %g %g\n", bbox[XLO] , bbox[YHI] , bbox[ZHI] );
  fprintf( outputFile , "%g %g %g\n", bbox[XLO] , bbox[YLO] , bbox[ZHI] );
	
  /* YLO-1 face. */
  fprintf( outputFile , "%g %g %g\n", bbox[XLO] , bbox[YLO] , bbox[ZLO] );
  fprintf( outputFile , "%g %g %g\n", bbox[XHI] , bbox[YLO] , bbox[ZLO] );
  fprintf( outputFile , "%g %g %g\n", bbox[XHI] , bbox[YLO] , bbox[ZHI] );
  fprintf( outputFile , "%g %g %g\n", bbox[XLO] , bbox[YLO] , bbox[ZHI] );
  fprintf( outputFile , "%g %g %g\n\n", bbox[XLO] , bbox[YLO] , bbox[ZLO] );
	
  /* ZHI face. */
  fprintf( outputFile , "%g %g %g\n", bbox[XHI] , bbox[YHI] , bbox[ZHI] );
  fprintf( outputFile , "%g %g %g\n", bbox[XLO] , bbox[YHI] , bbox[ZHI] );
  fprintf( outputFile , "%g %g %g\n", bbox[XLO] , bbox[YLO] , bbox[ZHI] );
  fprintf( outputFile , "%g %g %g\n", bbox[XHI] , bbox[YLO] , bbox[ZHI] );
	
  /* XHI face. */
  fprintf( outputFile , "%g %g %g\n", bbox[XHI] , bbox[YHI] , bbox[ZHI] );
  fprintf( outputFile , "%g %g %g\n", bbox[XHI] , bbox[YHI] , bbox[ZLO] );
  fprintf( outputFile , "%g %g %g\n", bbox[XHI] , bbox[YLO] , bbox[ZLO] );
  fprintf( outputFile , "%g %g %g\n", bbox[XHI] , bbox[YLO] , bbox[ZHI] );

  /* YHI face. */
  fprintf( outputFile , "%g %g %g\n", bbox[XHI] , bbox[YHI] , bbox[ZHI] );
  fprintf( outputFile , "%g %g %g\n", bbox[XHI] , bbox[YHI] , bbox[ZLO] );
  fprintf( outputFile , "%g %g %g\n", bbox[XLO] , bbox[YHI] , bbox[ZLO] );
  fprintf( outputFile , "%g %g %g\n", bbox[XLO] , bbox[YHI] , bbox[ZHI] );
  fprintf( outputFile , "%g %g %g\n\n", bbox[XHI] , bbox[YHI] , bbox[ZHI] );
      
  return;

}

/* Render an excitation arrow. */
void gnuplotBoundingBoxArrow( FILE *outputFile , int ibbox[6] , FieldComponent field )
{

  real bbox[6];
  real length = 1.0;
  real halfWidth = 0.5;

  bboxInRealUnits( bbox , ibbox );

  switch( field )
  {
  case EX:
  case HX:
    if( ibbox[XHI] == ibbox[XLO] )
    {
      ibbox[XHI] += 2;
      bboxInRealUnits( bbox , ibbox );
    }
    length = ( bbox[XHI] - bbox[XLO] ) / ( ibbox[XHI] - ibbox[XLO] );
    halfWidth = 0.5 * length;
    fprintf( outputFile , "%g %g %g\n" , bbox[XLO] , ( bbox[YLO] + bbox[YHI] ) * 0.5 , ( bbox[ZLO] + bbox[ZHI] ) * 0.5 );
    fprintf( outputFile , "%g %g %g\n" , bbox[XHI] , ( bbox[YLO] + bbox[YHI] ) * 0.5 , ( bbox[ZLO] + bbox[ZHI] ) * 0.5 );
    fprintf( outputFile , "%g %g %g\n" , bbox[XHI] - length , ( bbox[YLO] + bbox[YHI] ) * 0.5 - halfWidth , ( bbox[ZLO] + bbox[ZHI] ) * 0.5 - halfWidth );
    fprintf( outputFile , "%g %g %g\n" , bbox[XHI] - length , ( bbox[YLO] + bbox[YHI] ) * 0.5 + halfWidth , ( bbox[ZLO] + bbox[ZHI] ) * 0.5 + halfWidth );
    fprintf( outputFile , "%g %g %g\n\n" , bbox[XHI] , ( bbox[YLO] + bbox[YHI] ) * 0.5 , ( bbox[ZLO] + bbox[ZHI] ) * 0.5 );
    break;
  case EY:
  case HY:
    if( ibbox[YHI] == ibbox[YLO] )
    {
      ibbox[YHI] += 2;
      bboxInRealUnits( bbox , ibbox );
    }
    length = ( bbox[YHI] - bbox[YLO] ) / ( ibbox[YHI] - ibbox[YLO] );
    halfWidth = 0.5 * length;
    fprintf( outputFile , "%g %g %g\n" , ( bbox[XLO] + bbox[XHI] ) * 0.5, bbox[YLO] , ( bbox[ZLO] + bbox[ZHI] ) * 0.5 );
    fprintf( outputFile , "%g %g %g\n" , ( bbox[XLO] + bbox[XHI] ) * 0.5, bbox[YHI] , ( bbox[ZLO] + bbox[ZHI] ) * 0.5 );
    fprintf( outputFile , "%g %g %g\n" , ( bbox[XLO] + bbox[XHI] ) * 0.5 - halfWidth , bbox[YHI] - length , ( bbox[ZLO] + bbox[ZHI] ) * 0.5 - halfWidth );
    fprintf( outputFile , "%g %g %g\n" , ( bbox[XLO] + bbox[XHI] ) * 0.5 + halfWidth , bbox[YHI] - length , ( bbox[ZLO] + bbox[ZHI] ) * 0.5 + halfWidth );
    fprintf( outputFile , "%g %g %g\n\n" , ( bbox[XLO] + bbox[XHI] ) * 0.5 , bbox[YHI] , ( bbox[ZLO] + bbox[ZHI] ) * 0.5 );
    break;
  case EZ:
  case HZ:
    if( ibbox[ZHI] == ibbox[ZLO] )
    {
      ibbox[ZHI] += 2;
      bboxInRealUnits( bbox , ibbox );
    }
    length = ( bbox[ZHI] - bbox[ZLO] ) / ( ibbox[ZHI] - ibbox[ZLO] );
    halfWidth = 0.5 * length;
    fprintf( outputFile , "%g %g %g\n" , ( bbox[XLO] + bbox[XHI] ) * 0.5 , ( bbox[YLO] + bbox[YHI] ) * 0.5 , bbox[ZLO] );
    fprintf( outputFile , "%g %g %g\n" , ( bbox[XLO] + bbox[XHI] ) * 0.5 , ( bbox[YLO] + bbox[YHI] ) * 0.5 , bbox[ZHI] );
    fprintf( outputFile , "%g %g %g\n" , ( bbox[XLO] + bbox[XHI] ) * 0.5 - halfWidth , ( bbox[YLO] + bbox[YHI] ) * 0.5 - halfWidth , bbox[ZHI] - length );
    fprintf( outputFile , "%g %g %g\n" , ( bbox[XLO] + bbox[XHI] ) * 0.5 + halfWidth , ( bbox[YLO] + bbox[YHI] ) * 0.5 + halfWidth , bbox[ZHI] - length );
    fprintf( outputFile , "%g %g %g\n\n" , ( bbox[XLO] + bbox[XHI] ) * 0.5 , ( bbox[YLO] + bbox[YHI] ) * 0.5 , bbox[ZHI] );
    break;
  default:
    assert( 0 );
    break;
  }

  return;

}

/* Render a line/arrow defined by its start and end points and a direction for the head. */
void gnuplotArrow( FILE *outputFile , real start[3] , real end[3] , real norm[3] , int headStyle  )
{

  real pstart[3];                 // Start point in physical units.
  real pend[3];                   // End point in physical units.
  real headLengthFraction = 0.25; // Fractional length of head.
  real headLength;                // Length of head.
  real halfWidthFraction = 1.0;   // Fractional width of head.
  real halfWidth;                 // Width of head.

  if( isPhysicalUnits )
  {
    nodeInPhysicalUnits( pstart , start );
    nodeInPhysicalUnits( pend , end );
  }
  else
  {
    pstart[XDIR] = start[XDIR];
    pstart[YDIR] = start[YDIR];
    pstart[ZDIR] = start[ZDIR];
    pend[XDIR] = end[XDIR];
    pend[YDIR] = end[YDIR];
    pend[ZDIR] = end[ZDIR];
  }
  //printf( "*** --> [%.2f,%.2f,%.2f]->[%.2f,%.2f,%.2f]\n" , pstart[XDIR],pstart[YDIR],pstart[ZDIR],pend[XDIR],pend[YDIR],pend[ZDIR]);
  headLength = headLengthFraction * sqrt( ( pend[XDIR] - pstart[XDIR] ) * ( pend[XDIR] - pstart[XDIR] ) +
                                          ( pend[YDIR] - pstart[YDIR] ) * ( pend[YDIR] - pstart[YDIR] ) +
                                          ( pend[ZDIR] - pstart[ZDIR] ) * ( pend[ZDIR] - pstart[ZDIR] ) );
  halfWidth = halfWidthFraction * headLength;

  /* Render the line. */
  fprintf( outputFile , "%g %g %g\n" , pstart[XDIR] , pstart[YDIR] , pstart[ZDIR] );
  fprintf( outputFile , "%g %g %g\n" , pend[XDIR] , pend[YDIR] , pend[ZDIR] );

  /* Render the arrow head if required. */
  switch( headStyle )
  {
  case 0:
    /* Headless. */
    break;
  case 1:
    /* closed head. */
    fprintf( outputFile , "%g %g %g\n" , ( 1 - headLengthFraction ) * pend[XDIR] + headLengthFraction * pstart[XDIR] + halfWidth * norm[XDIR] , 
                                         ( 1 - headLengthFraction ) * pend[YDIR] + headLengthFraction * pstart[YDIR] + halfWidth * norm[YDIR] ,
                                         ( 1 - headLengthFraction ) * pend[ZDIR] + headLengthFraction * pstart[ZDIR] + halfWidth * norm[ZDIR] );
    fprintf( outputFile , "%g %g %g\n" , ( 1 - headLengthFraction ) * pend[XDIR] + headLengthFraction * pstart[XDIR] - halfWidth * norm[XDIR] , 
                                         ( 1 - headLengthFraction ) * pend[YDIR] + headLengthFraction * pstart[YDIR] - halfWidth * norm[YDIR] ,
                                         ( 1 - headLengthFraction ) * pend[ZDIR] + headLengthFraction * pstart[ZDIR] - halfWidth * norm[ZDIR] );
    fprintf( outputFile , "%g %g %g\n" , pend[XDIR] , pend[YDIR] , pend[ZDIR] );
    break;
  case 2:
    /* Open head. */
    fprintf( outputFile , "%g %g %g\n" , ( 1 - headLengthFraction ) * pend[XDIR] + headLengthFraction * pstart[XDIR] + halfWidth * norm[XDIR] , 
                                         ( 1 - headLengthFraction ) * pend[YDIR] + headLengthFraction * pstart[YDIR] + halfWidth * norm[YDIR] ,
                                         ( 1 - headLengthFraction ) * pend[ZDIR] + headLengthFraction * pstart[ZDIR] + halfWidth * norm[ZDIR] );
    fprintf( outputFile , "%g %g %g\n" , pend[XDIR] , pend[YDIR] , pend[ZDIR] );
    fprintf( outputFile , "%g %g %g\n" , ( 1 - headLengthFraction ) * pend[XDIR] + headLengthFraction * pstart[XDIR] - halfWidth * norm[XDIR] , 
                                         ( 1 - headLengthFraction ) * pend[YDIR] + headLengthFraction * pstart[YDIR] - halfWidth * norm[YDIR] ,
                                         ( 1 - headLengthFraction ) * pend[ZDIR] + headLengthFraction * pstart[ZDIR] - halfWidth * norm[ZDIR] );
    break;
  case 3:
    /* Double open head. */
    fprintf( outputFile , "%g %g %g\n" , ( 1 - headLengthFraction ) * pend[XDIR] + headLengthFraction * pstart[XDIR] + halfWidth * norm[XDIR] , 
                                         ( 1 - headLengthFraction ) * pend[YDIR] + headLengthFraction * pstart[YDIR] + halfWidth * norm[YDIR] ,
                                         ( 1 - headLengthFraction ) * pend[ZDIR] + headLengthFraction * pstart[ZDIR] + halfWidth * norm[ZDIR] );
    fprintf( outputFile , "%g %g %g\n" , pend[XDIR] , pend[YDIR] , pend[ZDIR] );
    fprintf( outputFile , "%g %g %g\n" , ( 1 - headLengthFraction ) * pend[XDIR] + headLengthFraction * pstart[XDIR] - halfWidth * norm[XDIR] , 
                                         ( 1 - headLengthFraction ) * pend[YDIR] + headLengthFraction * pstart[YDIR] - halfWidth * norm[YDIR] ,
                                         ( 1 - headLengthFraction ) * pend[ZDIR] + headLengthFraction * pstart[ZDIR] - halfWidth * norm[ZDIR] );
    fprintf( outputFile , "%g %g %g\n" , pend[XDIR] , pend[YDIR] , pend[ZDIR] );
    fprintf( outputFile , "%g %g %g\n" , ( 1 - 2 * headLengthFraction ) * pend[XDIR] + 2 * headLengthFraction * pstart[XDIR] + 0.5 * halfWidth * norm[XDIR] , 
                                         ( 1 - 2 * headLengthFraction ) * pend[YDIR] + 2 * headLengthFraction * pstart[YDIR] + 0.5 * halfWidth * norm[YDIR] ,
                                         ( 1 - 2 * headLengthFraction ) * pend[ZDIR] + 2 * headLengthFraction * pstart[ZDIR] + 0.5 * halfWidth * norm[ZDIR] );
    fprintf( outputFile , "%g %g %g\n" , pend[XDIR] , pend[YDIR] , pend[ZDIR] );
    fprintf( outputFile , "%g %g %g\n" , ( 1 - 2 * headLengthFraction ) * pend[XDIR] + 2 * headLengthFraction * pstart[XDIR] - 0.5 * halfWidth * norm[XDIR] , 
                                         ( 1 - 2 * headLengthFraction ) * pend[YDIR] + 2 * headLengthFraction * pstart[YDIR] - 0.5 * halfWidth * norm[YDIR] ,
                                         ( 1 - 2 * headLengthFraction ) * pend[ZDIR] + 2 * headLengthFraction * pstart[ZDIR] - 0.5 * halfWidth * norm[ZDIR] );
    break;
  default:
    assert( false );
    break;
  }
  
  /* Close the object. */
  fprintf( outputFile , "\n" );
   
  return;

}

/* Convert bounding box to physical units. */
void bboxInRealUnits( real physbbox[6] , int bbox[6] )
{

  if( isPhysicalUnits )
  {
    bboxInPhysicalUnits( physbbox , bbox );
  }
  else
  {
    physbbox[XLO] = (real) bbox[XLO];
    physbbox[XHI] = (real) bbox[XHI];
    physbbox[YLO] = (real) bbox[YLO];
    physbbox[YHI] = (real) bbox[YHI];
    physbbox[ZLO] = (real) bbox[ZLO];
    physbbox[ZHI] = (real) bbox[ZHI];
  }

  return;

}

/* Convert coordinate line index to physical units. */
real indexInRealUnits( int index , CoordAxis dir )
{

  if( isPhysicalUnits )  
    return indexInPhysicalUnits( index , dir );
  else
    return (real) index;

}


/* Determine if longest external edge is within ten percent of smallest external edge */
bool equalBoundaries()
{
  int side[3] = { mbox[XHI] - mbox[XLO] , mbox[YHI] - mbox[YLO] , mbox[ZHI] - mbox[ZLO] };

  int small = 0, large = 0;

  for (int i = 1; i < 3; ++i)
  {
    if (side[i] < side[small])
      small = i;
    if (side[i] > side[large])
      large = i;
  }

  if ( side[large] - side[small] <= side[small] / 10)
    return 1;

  return 0;
}