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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef WITH_OPENMP
#include <omp.h> 
#endif
  
#include "vulture.h"
#include "fdtd_types.h"
#include "physical.h"
#include "message.h"
#include "medium.h"
#include "boundary.h"
#include "surface.h"
#include "waveform.h"
#include "source.h"
#include "planewave.h"
#include "observer.h"
#include "mesh.h"
#include "simulation.h"
#include "block.h"
#include "line.h"
#include "grid.h"
#include "memory.h"


/* Vulture version. */
static int solverVersion[3] = { 0 , 7 , 0 };

/* Support mesh versions. */
static int meshVersion[2][3] = { { 0 , 0 , 0 } , { 1 , 0 , 0 } };

/* Options. */
struct Options_t {

  MessageType logLevel;
  bool readOnly;
  bool preprocessOnly;
  bool dumpGrid;
  int numThread;

} options = { MSG_LOG , false , false , false , -1 };

/* Private functions. */
void parseOption( int argc , char *argv[] , char meshFileName[] );
void printUsage( void );
void printVersion( void );
void printLicence( void );


/* Main. */
int main ( int argc , char **argv )
{

  /* Mesh file name */
  char meshFileName[PATH_SIZE] = "";

  /* Parse options. */
  parseOption( argc , argv , meshFileName );

  /* Start logging. */
  startMessaging( "vulture.log" , options.logLevel , "Vulture" , solverVersion[0] , solverVersion[1]  , solverVersion[2] );

  /* Define physical constants. */
  physicalConstants();

  /* Initialise simulation. */
  initSimulation();

  /* Initialise mesh. */
  initMesh();

  /* Read in the mesh. */
  readMesh( meshFileName );
  if( options.readOnly ) exit( 0 );

  /* Initialise the grid. */
  initGrid();

  /* Initialise boundary types - must be done before internal/external surfaces. */
  initBoundaries();

  /* Initialise media types - must be done before blocks. */
  initMedia();

  /* Initalise grid media to free space. */
  /* Sets all grid, including gibox/gobox surfaces and ghost cells to free-space. */
  initMediaArrays();
  
  /* Apply blocks to grid - must be done before external PML/MUR boundaries. */
  /* This can set medium parameters on the gibox surfaces. */
  initBlocks();

  /* Apply lines to grid - must be done before external PML/MUR boundaries. */
  /* This can set medium parameters on the gibox surfaces. */
  initLines();

  /* Initialise internal surfaces - must be done before external PML/MUR boundaries. */
  /* This can set medium parameters on the gibox surfaces. */
  initInternalSurfaces();
  
  /*  Initialise the external boundaries. */
  initExternalSurfaces();

  /* Initialise the waveforms. */
  initWaveforms();

  /* Initialise the sources. */
  initSources();
  initPlaneWaves();
  
  /* Initialise the observers. */
  initObservers();

  /* Free the mesh. */
  deallocMesh();
 
  /* Log grid report. */
  reportGrid();

  /* Report memory usage. */
  reportMemory();

  /* Limit checking mode. */
#ifdef CHECK_LIMITS
  setNumTimeSteps( 1 );
  message( MSG_LOG , 0 , "\n*** Limit checking run with one timestep ***\n\n" );  
#endif

  if( options.dumpGrid )
  {
    dumpMediaOnGrid( EX );
    dumpMediaOnGrid( EY );
    dumpMediaOnGrid( EZ );
    dumpMediaOnGrid( HX );
    dumpMediaOnGrid( HY );
    dumpMediaOnGrid( HZ );
  }

  /* Set number of threads if given.*/
#ifdef WITH_OPENMP
  if( options.numThread > 0 )  
    omp_set_num_threads( options.numThread ); 
#endif
  
  /* Step the fields. */
  if( !options.preprocessOnly )
    propagate();

#ifdef CHECK_LIMITS
  checkGrid();
  checkExternalSurfaces();
#endif

  /* Tidy up. */
  deallocObservers();
  deallocPlaneWaves();
  deallocSources();
  deallocWaveforms();
  deallocExternalSurfaces();
  deallocInternalSurfaces();
  deallocBoundaries();
  deallocLines();
  deallocBlocks();
  deallocMedia();
  deallocGridArrays();
  deallocSimulation();
  stopMessaging();

  return 0;

}

/* Parse command line options. */
void parseOption( int argc , char *argv[] , char meshFileName[] )
{

  char *ptr;
  
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
    else if( strncmp( argv[1] , "-m" , 2 ) == 0  || strncmp( argv[1] , "--read-mesh" , 11 ) == 0 )
    {
      options.readOnly = true;
    }
    else if( strncmp( argv[1] , "-n" , 2 ) == 0  || strncmp( argv[1] , "--num-proc" , 10 ) == 0 )
    {
      if( argc > 0 )
      {
        options.numThread = strtod( argv[2] , &ptr );      
        if( options.numThread == 0 )
        {
          printf( "\n*** Error: invalid value %s for option %s\n" , argv[2] , argv[1] );
          printUsage();
          exit( 1 );         
        }
        ++argv;
        --argc;
      }
      else
      {
        printf( "\n*** Error: no value for option %s\n" , argv[1] );
        printUsage();
        exit( 1 );
      }

    }
    else if( strncmp( argv[1] , "-p" , 2 ) == 0  || strncmp( argv[1] , "--preprocess" , 12 ) == 0 )
    {
      options.preprocessOnly = true;
    }
    else if( strncmp( argv[1] , "-g" , 2 ) == 0  || strncmp( argv[1] , "--dump-grid" , 11 ) == 0 )
    {
      options.dumpGrid = true;
    }
    else if( strncmp( argv[1] , "-l" , 2 ) == 0  || strncmp( argv[1] , "--licence" , 11 ) == 0 )
    {
      printLicence();
      exit( 0 );
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

/* Print solver usage information to standard output. */
void printUsage( void )
{

  printf( "\nUsage:\n\n" );
  printf( "vulture -h | --help\n" );
  printf( "vulture -V | --version\n" );
  printf( "vulture [ option ] <meshFile>\n\n" );
  printf( "Valid options are:\n\n" );
  printf( "-g, --dump-grid\t\t\tWrite out grid in ASCII format\n" );
  printf( "-m, --readmesh\t\t\tRead the mesh only and stop\n" );
  printf( "-n <int>, --numproc <int> \tSet number of threads\n" );
  printf( "-p, --preprocess\t\tPreprocess the mesh only and stop\n" );
  printf( "-v, --verbose\t\t\tProduce verbose logging information\n\n" );

  return;

}

/* Print solver version information to standard output. */
void printVersion( void )
{

  printf( "\nVulture (version %d.%d.%d) Copyright (C) 2011-2016 Ian David Flintoft\n\n" , 
          solverVersion[0] , solverVersion[1] , solverVersion[2] );
  printf( "Vulture comes with ABSOLUTELY NO WARRANTY; for details type `vulture --licence'.\n" );
  printf( "This is free software, and you are welcome to redistribute it\n" );
  printf( "under certain conditions; type `vulture --licence' for details.\n\n" );

  printf( "  Supported mesh versions %d.%d.%d - %d.%d.%d\n" , meshVersion[0][0] , meshVersion[0][1] , meshVersion[0][2] ,
                                                             meshVersion[1][0] , meshVersion[1][1] , meshVersion[1][2] );
#ifdef WITH_OPENMP
  printf( "  Built with OpenMP parallelisation support.\n" );
#endif
#ifdef USE_SCALED_FIELDS
  printf( "  Using scaled fields.\n" );
#else
  printf( "  Using un-scaled fields.\n" );  
#endif
  printf( "  Field arrays are %d-bytes.\n" , (int) sizeof( real ) );
#ifdef USE_INDEXED_MEDIA
  printf( "  Using indexed media.\n" );
  printf( "  Medium index is %d-bytes.\n" , (int) sizeof( MediumIndex ) );
#else
  printf( "  Using unindexed media.\n" );
  printf( "  Media arrays are %d-bytes.\n" , (int) sizeof( real ) );
#endif  
#ifdef USE_AVERAGED_MEDIA
  printf( "  Using averaged media.\n" );
#else
  printf( "  Using un-averaged media.\n" );  
#endif  
  
  printf( "\n" );

  return;
 
}

/* Print licence information to standard output. */
void printLicence( void )
{
 
  printf( "\nVulture finite-difference time-domain electromagnetic solver.\n" );
  printf( "Copyright (C) 2011-2016 Ian David Flintoft\n" );
  printf( "\n" );
  printf( "This program is free software; you can redistribute it and/or modify\n" );
  printf( "it under the terms of the GNU General Public License as published by\n" );
  printf( "the Free Software Foundation; either version 3 of the License, or\n" );
  printf( "(at your option) any later version.\n" );
  printf( "\n" );
  printf( "This program is distributed in the hope that it will be useful,\n" );
  printf( "but WITHOUT ANY WARRANTY; without even the implied warranty of\n" );
  printf( "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n" );
  printf( "GNU General Public License for more details.\n" );
  printf( "\n" );
  printf( "You should have received a copy of the GNU General Public License\n" );
  printf( "along with this program; if not, write to the Free Software Foundation,\n" );
  printf( "Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA\n" );
  printf( "or go to the web-site http://gnu.org/licenses/gpl.html.\n\n" );

  return;
 
}

