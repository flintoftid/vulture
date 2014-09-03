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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "fdtd_types.h"
#include "physical.h"
#include "message.h"
#include "alloc_array.h"
#include "mesh.h"
#include "medium.h"
#include "boundary.h"
#include "surface.h"
#include "pml.h"
#include "waveform.h"
#include "source.h"
#include "planewave.h"
#include "observer.h"
#include "gnuplot.h"
#include "gmsh.h"

/* Vulture version. */
static int solverVersion[3] = { 0 , 7 , 0 };

/* Support mesh versions. */
static int meshVersion[2][3] = { { 0 , 0 , 0 } , { 1 , 0 , 0 } };

/* Output formats. */
#define NUM_GRAPHICS_FORMATS 2

/* Graphics formats. */
typedef enum {

  GF_GNUPLOT,
  GF_GMSH,
  GF_UNDEFINED

} GraphicsFormat;

/* Options. */
struct Options_t {

  MessageType logLevel;
  bool isPhysicalUnits; 
  GraphicsFormat format;
  bool isExternalSurfaces;

} options = { MSG_LOG , false , GF_GNUPLOT , true };

/* Private functions. */
void parseOption( int argc , char *argv[] , char meshFileName[] );
void printUsage( void );
void printVersion( void );

/* Main. */
int main ( int argc , char **argv )
{

  /* Mesh file name */
  char meshFileName[PATH_SIZE] = "";

  /* Parse options. */
  parseOption( argc , argv , meshFileName );

  /* Start logging. */
  startMessaging( "gvulture.log" , options.logLevel , "gvulture" , solverVersion[0] , solverVersion[1]  , solverVersion[2] );

  /* Define physical constants. */
  physicalConstants();

  /* Initialise mesh. */
  initMesh();

  /* Read in the mesh. */
  readMesh( meshFileName );

  /* Render mesh. */
  switch( options.format )
  {
    case GF_GNUPLOT:
      gnuplotMesh( options.isPhysicalUnits , options.isExternalSurfaces );
      break;
    case GF_GMSH:
      gmshMesh( options.isPhysicalUnits , options.isExternalSurfaces );
      break;
    default:
      assert( 0 );
      break;
  }
  
  /* Free the mesh. */
  deallocMesh();
 
  /* Tidy up. */
  stopMessaging();

  return 0;

}

/* Parse command line options. */
void parseOption( int argc , char *argv[] , char meshFileName[] )
{

  while ( ( argc > 1 ) && ( argv[1][0] == '-' ) )
  {

    if( strncmp( argv[1] , "-h" , 2 ) == 0  || strncmp( argv[1] , "--help" , 6 ) == 0 )
    {
      printUsage();
      exit( 0 );
    }
    else if( strncmp( argv[1] , "-V" , 2 ) == 0  || strncmp( argv[1] , "--version" , 9 ) == 0 )
    {
      printVersion();
      exit( 0 );
    }
    else if( strncmp( argv[1] , "-v" , 2 ) == 0  || strncmp( argv[1] , "--verbose" , 9 ) == 0 )
    {
      options.logLevel = MSG_DEBUG3;
    }
    else if( strncmp( argv[1] , "-p" , 2 ) == 0  || strncmp( argv[1] , "--physical" , 10 ) == 0 )
    {
      options.isPhysicalUnits = true;
    }
    else if( strncmp( argv[1] , "-e" , 2 ) == 0  || strncmp( argv[1] , "--no-ext-surf" , 13 ) == 0 )
    {
      options.isExternalSurfaces = false;
    }
    else if( strncmp( argv[1] , "-g" , 2 ) == 0  || strncmp( argv[1] , "--gnuplot" , 9 ) == 0 )
    {
      options.format = GF_GNUPLOT;
    }
    else if( strncmp( argv[1] , "-m" , 2 ) == 0  || strncmp( argv[1] , "--gmsh" , 6 ) == 0 )
    {
      options.format = GF_GMSH;
    }
    else
    {
      printf( "\n*** Error: invalid option %s\n" , argv[1] );
      printUsage();
      exit( 1 );
    }

    ++argv;
    --argc;

  } /* while */

  if( argc != 2 )
    printUsage();
  else
    strncpy( meshFileName , argv[1] , PATH_SIZE );

  return;

}

/* Print usage information to standard output. */
void printUsage( void )
{

  printf( "\nUsage:\n\n" );
  printf( "gvulture -h | --help\n" );
  printf( "gvulture -V | --version\n" );
  printf( "gvulture [ option ] <meshFile>\n\n" );
  printf( "\nValid options are:\n\n" );
  printf( "-e, --no-ext-surf\tDo not render mesh external surfaces\n" );  
  printf( "-g, --gnuplot\t\tGenerate gnuplot format output (default)\n" );  
  printf( "-p, --physical\t\tGenerate plot data in physical units (metres)\n" );
  printf( "-v, --verbose\t\tProduce verbose logging information\n\n" );

  return;

}

/* Print solver version information to standard output. */
void printVersion( void )
{

  printf( "\nVulture gnuplot generator version %d.%d.%d\n\n" , solverVersion[0] , solverVersion[1] , solverVersion[2] );

  printf( "  Supported mesh versions %d.%d.%d - %d.%d.%d\n" , meshVersion[0][0] , meshVersion[0][1] , meshVersion[0][2] ,
                                                             meshVersion[1][0] , meshVersion[1][1] , meshVersion[1][2] );
  printf( "\n" );

  return;
 
}

