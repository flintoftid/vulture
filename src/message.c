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

#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include "message.h"

static FILE *logFile;
static char *programName;
static MessageType logLevel;


/* Start logging messages. */
/* Messages of level less than minimumLevel are ignored. */
void startMessaging( char *logFileName , MessageType minimumLevel , char *progName , int versionMajor , int versionMinor , int versionPatch )
{

  programName = (char *) malloc( ( strlen( progName ) + 1 ) * sizeof( char ) );
  if( !programName )
  {
    fprintf( stderr , "%s: ERROR: Failed to allocate program name string\n" , progName );
    fflush( stderr );
    exit( 1 );
  }
  strncpy( programName , progName , strlen( progName ) );
  programName[strlen( progName )] = '\0';
  
  logLevel = minimumLevel;

  logFile = fopen( logFileName , "w" );
  if( !logFile )
  {
    fprintf( stderr , "%s: ERROR: Failed to open log file %s\n" , programName , logFileName );
    fflush( stderr );
    exit( 1 );
  }

  fprintf( logFile , "\n *** %s version %d.%d.%d *** \n\n" , programName , versionMajor , versionMinor , versionPatch );

  fflush( logFile );

}

/* Send messgage to logger. */
void message(  MessageType status , int errnum , const char *message , ... )
{

  va_list ap;
  va_list ap2;
  
  if( status < logLevel ) 
    return;

  if( status >= MSG_WARN )
    fprintf ( stderr, "%s: " , programName );

  va_start( ap , message );
  if( status >= MSG_WARN )
    vfprintf ( stderr , message , ap );
  va_end ( ap );
  
  va_start( ap2 , message );  
  vfprintf ( logFile , message , ap2 );
  va_end ( ap2 );

  if ( errnum )
  {
    if( status >= MSG_WARN )
    {
      fprintf ( stderr , ": %s" , strerror ( errnum ) );
    }

    fprintf ( logFile , ": %s" , strerror ( errnum ) );
  }

  fflush( stderr );
  fflush( logFile );

  if ( status >= MSG_ERROR )
     exit( MSG_ERROR );

}

/* Stop logging mesaages. */
void stopMessaging( void )
{

  free( programName );
  fflush( logFile );
  fclose( logFile );

  return;

}
