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
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include "mesh.h"
#include "message.h"
#include "waveform.h"
#include "source.h"
#include "planewave.h"
#include "simulation.h"
#include "observer.h"
#include "medium.h"
#include "block.h"
#include "boundary.h"
#include "wire.h"
#include "line.h"
#include "surface.h"
#include "grid.h"


/* Number of sections in mesh file. */
#define NUM_SECTIONS        4

/* Number of cards. */
#define NUM_CARDS          27

/* Length of card tag (including '\0'). */
#define CARD_STRING_LENGTH  3

/* Card types - values used to index arrays! */
typedef enum {

 /* Section 0 */
 CT_VM = 0, 
 /* Section 1 */ 
 CT_CE , CT_DM , CT_BR , CT_BE , CT_GS , 
 /* Section 2 */ 
 CT_MT , CT_BT , CT_WT , CT_MB , CT_TB , CT_TW , CT_WF , CT_EX , CT_PW , CT_OP , CT_FF , CT_GE , 
 /* Section 3 */ 
 CT_NT , CT_OT , CT_OF , CT_CN , CT_MS , CT_XL , CT_YL , CT_ZL , CT_EN ,
 /* Error condition. */
 CT_ERROR 

} CardType;

/* FILE pointer to mesh file. */
static FILE *meshFilePointer;

/* String for each card indexed by card type above. */
static char CT_STRING[NUM_CARDS][CARD_STRING_LENGTH] =
{ 
  /* Section 0 */
  "VM" ,
  /* Section 1 */ 
  "CE" , "DM" , "BR" , "BE" , "GS" , 
  /* Section 2 */  
  "MT" , "BT" , "WT" , "MB" , "TB" , "TW" , "WF" , "EX" , "PW" , "OP" , "FF" , "GE" , 
  /* Section 3 */
  "NT" , "OT" , "OF" , "CN" , "MS" , "XL" , "YL" , "ZL" , "EN" 
};

/* Card limits for each section. */
CardType sectionLimits[NUM_SECTIONS][2] = {
  { CT_VM , CT_VM } ,
  { CT_CE , CT_GS } ,
  { CT_MT , CT_GE } ,
  { CT_NT , CT_EN } ,
};

/* Which cards are valid in each section */
/* Redundant - can be determined from above. */
static bool isValidInSection[NUM_SECTIONS][NUM_CARDS] =
{
  { /* Section 0 */
    true , 
    false , false , false , false , false , 
    false , false , false , false , false , false , false , false , false , false , false , false ,
    false , false , false , false , false , false , false , false , false
  } ,
  { /* Section 1 */
    false , 
    true  , true  , true  , true  , true  , 
    false , false , false , false , false , false , false , false , false , false , false , false ,
    false , false , false , false , false , false , false , false , false
  } ,
  { /* Section 2 */
    false , 
    false , false , false , false , false , 
    true  , true  , true  , true  , true  , true  , true  , true  , true  , true  , true  , true  ,
    false , false , false , false , false , false , false , false , false
  } ,
  { /* Section 3 */
    false , 
    false , false , false , false , false , 
    false , false , false , false , false , false , false , false , false , false , false , false ,
    true  , true , true  , true  , true  , true  , true  , true  , true
  }
};

/* Which cards are required. */
static bool isRequired[NUM_CARDS] =
{
  /* Section 0 */
  true  , 
  /* Section 1 */
  false , true  , false  , false , true  , 
  /* Section 2 */
  false , false , false , false , false , false , false , false , false , false , false , true  ,
  /* Section 3 */
  true  , false , false , false , false , false , false , false , true
};

/* Which cards can only be given once. */
static bool isSingleton[NUM_CARDS] =
{
  /* Section 0 */
  true  , 
  /* Section 1 */
  false , true  , true  , true  , true  , 
  /* Section 2 */
  false , false , false , false , false , false , false , false , false , false , false , true  ,
  /* Section 3 */
  true  , true , true  , true  , true  , true  , true  , true  , true
};

/* Which cards mark ends of sections. */
/* Redundant - can be determined from above. */
static bool isEndOfSection[NUM_CARDS] =
{
  /* Section 0 */
  true  , 
  /* Section 1 */
  false , false , false , false , true  , 
  /* Section 2 */
  false , false , false , false , false , false , false , false , false , false , false , true  ,
  /* Section 3 */
  false , false , false , false , false , false , false , false , true
};

/* Which cards have been found. */
static bool isFound[NUM_CARDS] =
{
  /* Section 0 */
  false , 
  /* Section 1 */
  false , false , false , false , false , 
  /* Section 2 */
  false , false , false , false , false , false , false , false , false , false , false , false ,
  /* Section 3 */
  false , false , false , false , false , false , false , false , false
};

/* Pointers to card handling functions. */
static bool (*cardFunc[NUM_CARDS])( char * ) = 
{
  /* Section 0 */
  &parseVM , 
  /* Section 1 */
  &parseCE , &parseDM , &parseBR , &parseBE , NULL , 
  /* Section 2 */
  &parseMT , &parseBT , &parseWT , &parseMB , &parseTB , &parseTW , &parseWF , &parseEX , &parsePW , &parseOP , &parseFF , NULL ,
  /* Section 3 */
  &parseNT , &parseOT , &parseOF , &parseCN , &parseMS , &parseXL , &parseYL , &parseZL , NULL 
};

/* Size of line buffer. */
#define BUFFER_LENGTH 1024

static int meshVersionMajor;             // Major version number. 
static int meshVersionMinor;             // Minor number.
static int meshVersionPatch;             // Patch number
static unsigned long meshVersion;        // Mesh version number hash.
static char comment[COMMENT_SIZE];       // Mesh title comment.

/* Private functions.*/
void removeComments( char *string );
void strtrim( char *string );
CardType getCardType( char *line , int lineNumber );
unsigned long meshVersionNumber( int versionMajor , int versionMinor , int versionPatch );
void setComment( char *str );


/* Parse mesh version. */
bool parseVM( char *line )
{

  if( sscanf( line , "%d.%d.%d" , &meshVersionMajor , &meshVersionMinor , &meshVersionPatch ) != 3 )
    return false;

  meshVersion = meshVersionNumber( meshVersionMajor , meshVersionMinor , meshVersionPatch );

  return true;

}

/* Simple hash function for mesh version number. */
unsigned long meshVersionNumber( int versionMajor , int versionMinor , int versionPatch )
{

  if( versionMajor >= 2000 || versionMinor >= 1000 || versionPatch >= 1000 )
    message( MSG_ERROR , 0 , "Mesh version number (%d.%d.%d) too large for hashing function!\n" , 
             versionMajor , versionMinor , versionPatch );

  return 1000000 * versionMajor + 1000 * versionMinor + versionPatch;

}

/* Parse mesh comment. */
bool parseCE( char *line )
{

  setComment( line );
 
  return true;

}

/* Initialise the mesh. */
void initMesh( void )
{

  message( MSG_LOG , 0 , "\nInitialising mesh...\n\n" );

  /* Predefine PEC and free space boundaries. */
  addBoundary( "PEC" , BT_PEC , 0 , 0 , 1.0 , -1.0 , 1.0 , "" , NULL , NULL );
  addBoundary( "PMC" , BT_PMC , 0 , 0 , 1.0 , +1.0 , 1.0 , "" , NULL , NULL );
  addBoundary( "FREE_SPACE" , BT_FREE_SPACE , 0 , 0 , 1.0 , 0.0 , 1.0 , "" , NULL , NULL );

  /* Predefine PEC and free space boundaries. */
  addWire( "PEC" , TW_PEC , 0.0 );
  addWire( "FREE_SPACE" , TW_FREE_SPACE , 0.0 );

  /* Predefine free space and PEC media. */
  addMedium( "FREE_SPACE" , MT_SIMPLE , 1.0 , 0.0 , 1.0 , 0 , NULL , NULL , "" ); 
  addMedium( "PEC" , MT_PEC , 1.0 , 1e8 , 1.0 , 0 , NULL , NULL , "" ); 

  return;
 
}

/* Read the input mesh. */
void readMesh( char fileName[] )
{
  /* Generic parser doesn't handle dependencies on group of cards. */
  /* These flags are used to check a source/output is specficied. */
  bool foundSource = false;
  bool foundOutput = false;

  /* Current input line after processing. */
  char line[BUFFER_LENGTH];

  /* Raw input line used for error reporting. */
  char rawLine[BUFFER_LENGTH];

  /* Line number in mesh file. */
  unsigned long lineNumber = 0;

  /* Currnet section number. */
  int sectionNumber = 0;

  /* Type of current card. */
  CardType thisCardType;

  message( MSG_LOG , 0 , "\n  Reading input mesh file %s...\n\n" , fileName );
 
  /* Open file. */
  meshFilePointer = fopen( fileName , "r" );
  if( meshFilePointer == NULL ) 
     message( MSG_ERROR , 0 , "  ***Error: Cannot open input mesh file %s\n" , fileName );

  while( !feof( meshFilePointer ) ) 
  {

    if( !fgets( rawLine , BUFFER_LENGTH - 1 , meshFilePointer ) )
    {
      if( feof( meshFilePointer ) )
        break;
      else
        message( MSG_ERROR , 0 , "*** Error reading %s after line number %d\n" , fileName , lineNumber );
    }

    lineNumber++; 

    /* Make copy for parsing. Keep original for error reporting. */
    strncpy( line , rawLine , BUFFER_LENGTH );

    /* Remove any comments. */
    removeComments( line );

    /* Trim leading and trailing whitespace. */
    strtrim( line );

    /* Empty line. */
    if( strlen( line ) == 0 )
      continue;
 
    /* Get card type from front of line. */
    thisCardType = getCardType( line , lineNumber );
    if( thisCardType == CT_ERROR )
      message( MSG_ERROR , 0 , "*** Error: Invalid card on line %lu:\n  %s" , lineNumber , rawLine ); 

    /* Check card is valid in this section. */
    if( !isValidInSection[sectionNumber][thisCardType] )
      message( MSG_ERROR , 0 , "*** Error: Card type %s on line %lu is invalid in mesh section %d\n" , CT_STRING[thisCardType] , lineNumber , sectionNumber );

    /* Check singleton cards not repeated. */
    if( isFound[thisCardType] && isSingleton[thisCardType] )
      message( MSG_ERROR , 0 , "*** Error: Card type %s on line %lu is has already been found\n" , CT_STRING[thisCardType] , lineNumber );

    /* Mark card as fouend. */
    isFound[thisCardType] = true;

    message( MSG_DEBUG1 , 0 , "  [%d] %2s: %s" , lineNumber , CT_STRING[thisCardType] , rawLine );

    /* Invoke generic card handler. */
    if( (*cardFunc[thisCardType]) )
      if( !(*cardFunc[thisCardType])( line + 2 ) )
        message( MSG_ERROR , 0 , "*** Error: Failed to parse %s card on line %lu:\n  %s" , CT_STRING[thisCardType] , lineNumber , rawLine );

    /* Specific semantics not supported by generic parser. */
    switch( thisCardType )
    {
    case CT_EX:
    case CT_PW:
      foundSource = true;
      break;
    case CT_OP:
    case CT_FF:
      foundOutput = true;
      break;
    case CT_GE:
      if( !foundSource )
        message( MSG_ERROR , 0 , "*** Error: No source found in section 2 of mesh file!\n" );      if( !foundOutput )
        message( MSG_WARN , 0 , "*** Warning: No outputs found in section 2 of mesh file!\n" ); 
      break;
    case CT_MS:
      isFound[CT_MS] = isFound[CT_XL] = isFound[CT_YL] = isFound[CT_ZL] = true;
      break;
    case CT_XL:
      isFound[CT_MS] = isFound[CT_XL] = true; 
      break;
    case CT_YL:
      isFound[CT_MS] = isFound[CT_YL] = true; 
      break;
    case CT_ZL:
      isFound[CT_MS] = isFound[CT_ZL] = true; 
      break;
    case CT_EN:
      if( !isFound[CT_MS] || ! ( isFound[CT_XL] && isFound[CT_YL] && isFound[CT_ZL] ) )
        message( MSG_ERROR , 0 , "*** Error: No mesh lines found in section 3 of mesh file!\n" );
      break;
    default:
      /* Handled generically. */
      break;
    }

   /* End of section processing. */
   if( isEndOfSection[thisCardType] )
   {
     /* Check required cards are found. */
     for( CardType card = sectionLimits[sectionNumber][0] ; card <= sectionLimits[sectionNumber][1] ; card++ )
       if( isRequired[card] && !isFound[card] )
         message( MSG_ERROR , 0 , "*** Error: Required card %s not found in Section %d\n" , CT_STRING[card] , sectionNumber );         
   
     /* Advance to next section.*/
     message( MSG_LOG , 0 , "  Processed Section %d\n\n" , sectionNumber );   
     sectionNumber++; 
   }

  } /* End of switch. */

  fclose(meshFilePointer);

  reportMesh();

  return;

}

/* Remove leading and trailing whitespace from string. */
void strtrim( char *string )
{
  size_t stringLength = 0;
  char *startPointer = string - 1;
  char *endPointer = NULL;

  if( string == NULL )
    return;

  if( string[0] == '\0' )
    return;

  stringLength = strlen( string );
  endPointer = string + stringLength;

  /* Move the start and end pointers to the first 
   *  non-whitespace characters from each end.
   */
  while( isspace( *(++startPointer) ) );
  while( isspace( *(--endPointer) ) && endPointer != startPointer );

  if( string + stringLength - 1 != endPointer )
    *(endPointer + 1) = '\0';
  else if( startPointer != string && endPointer == startPointer )
    *string = '\0';

  /* Shift the string back to the start of the buffer. */
  endPointer = string;
  if( startPointer != string )
  {
    while( *startPointer ) *endPointer++ = *startPointer++;
    *endPointer = '\0';
  }

  return;

}

/* Remove trailing comment from string. */
void removeComments( char *string )
{
  char *charPointer = string - 1;

  if( string == NULL )
    return;

  if( string[0] == '\0' )
    return;

  while( *(++charPointer) )
  {
    if( *charPointer == '#' )
    {
      *charPointer = '\0';
      return;
    }        
  }

  return;

}

/* Get card type from line and remove it if successful. */
CardType getCardType( char *line , int lineNumber )
{

  char cardString[CARD_STRING_LENGTH];

  if( strlen( line ) < CARD_STRING_LENGTH - 1 )
    return CT_ERROR;

  strncpy( cardString , line , CARD_STRING_LENGTH - 1 );
   
  for( CardType card = CT_VM ; card <= CT_EN ; card++ )
    if( strncmp( cardString , CT_STRING[card] , CARD_STRING_LENGTH - 1 ) == 0 )
        return card;

  return CT_ERROR;

}

/* Read array to reals from mesh file. */
bool meshReadRealArray( int size , real *v )
{

  double value;

  for( int k = 0 ; k < size ; k++ )
  {
    if( fscanf( meshFilePointer , "%lf" , &value ) != 1 )
      return false;
    else
      v[k] = (real)value;
  }

  return true;

}

/* Report mesh information. */
void reportMesh( void )
{

  message( MSG_LOG , 0 , "\nMesh characteristics:\n\n" );

  message( MSG_LOG , 0 , "  Mesh version: %d.%d.%d\n" , meshVersionMajor , meshVersionMinor , meshVersionPatch );
  message( MSG_LOG , 0 , "  Mesh title: %s\n" , comment );

  reportSimulation();
  reportBoundaries();
  reportMedia();
  reportWires();
  reportSurfaces();
  reportBlocks();
  reportLines();
  reportWaveforms();
  reportSources();
  reportPlaneWaves();
  reportObservers();

  message( MSG_LOG , 0 , "\n" );

  return;

}

/* Deallocate mesh arrays. */
void deallocMesh( void )
{

  message( MSG_LOG , 0 , "\nDeallocating the mesh...\n\n" );

  return;

}

/* Set comment string. */
void setComment( char *str )
{

  strncpy( comment , str , COMMENT_SIZE );
  
  return;


}

/* Get pointer to comment string. */
char *getCommentReference( void )
{

  return comment;

}

