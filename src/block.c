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

#include <string.h>

#include "block.h"
#include "debye.h"
#include "utlist.h"
#include "medium.h"
#include "grid.h"
#include "alloc_array.h"
#include "message.h"
#include "bounding_box.h"
#include "gnuplot.h"
#include "gmsh.h"
#include "memory.h"

/* 
 * Private data. 
 */

/* Number of blocks. */
static BlockIndex numBlock = 0;

/* Existance flag for blocks of each medium type, including undefined. */
static bool isBlockType[NUM_MEDIUM_TYPES+1] = { false };

/* List of blocks. */
static BlockItem *blockList = NULL;

/* 
 * Private method interfaces. 
 */

void setMediumOnMesh( MediumIndex ***cellArray , int mbbox[6] , MediumIndex medium );

/*
 * Method Implementations.
 */

/* Parse material blocks. */
bool parseMB( char *line )
{

  int numScanned = 0;
  char mediumName[TAG_SIZE] = "";
  int mbbox[6];
  MediumIndex mediumNumber = 0;
  char maskStr[TAG_SIZE] = "111111";
  FaceMask mask = FACE_MASK_ALL;
  
  numScanned = sscanf( line , "%d %d %d %d %d %d %31s %31s" , 
                       &mbbox[XLO] , &mbbox[XHI] , &mbbox[YLO] , &mbbox[YHI] , &mbbox[ZLO] , &mbbox[ZHI] , mediumName , maskStr );

  if( numScanned < 7 )
    return false;  

  /* Validate bounding box. */ 
  if( !bboxIsNormal( mbbox ) )
  {
    message( MSG_LOG , 0 , "  Bounding box is abnormal!\n" );
    return false;
  }
  else if( !bboxIsWithin( mbbox , mbox ) )
  {
    message( MSG_LOG , 0 , "  Bounding box is outside mesh!\n" );
    return false;
  }
  else if( bboxType( mbbox ) != BB_VOLUME )
  {
    message( MSG_LOG , 0 , "  Bounding box is not a volume!\n" );
    return false;
  }

  /* Check medium exists. */
  if( !isMedium( mediumName , &mediumNumber ) )
  {
    message( MSG_LOG , 0 , "  Medium %s not defined in MB card\n" , mediumName );
    return false;
  }

  if( numScanned <= 8 )
  {
    mask = setFaceMaskFromString( maskStr );
    if( mask == FACE_MASK_ERROR )
    {
      message( MSG_LOG , 0 , "  Face mask %s is invalid\n" , maskStr );
      return false;
    }      
  }

  addBlock( mbbox , mediumName , mediumNumber , mask );

  return true;

}

/* Add block to lists. */
void addBlock( int mbbox[6] , char *mediumName , MediumIndex mediumNumber , FaceMask mask )
{

  BlockItem *item = NULL;

  if( numBlock == MAX_BLOCK )
    message( MSG_ERROR , 0 , "*** Error: Maximum number of blocks exceeded!\n" );

  item = (BlockItem *) malloc( sizeof( BlockItem ) );
  if( !item )
    message( MSG_ERROR , 0 , "*** Error: Failed to allocate block.\n" );

  for( int i = XLO ; i <= ZHI ; i++ ) item->mbbox[i] = mbbox[i];
  item->mediumNumber = mediumNumber;
  strncpy( item->mediumName , mediumName , TAG_SIZE );
  item->mask = mask;

  /* Add to list. */ 
  DL_APPEND( blockList , item );
  numBlock++;
  isBlockType[MT_UNDEFINED] = true;
    
  return;

}

/* Initialise blocks. */
/* Depends: initGrid, initMedia, initExternalPecPmcSurfaces */
void initBlocks( void )
{

  BlockItem *item;
  int gbbox[6];
  MediumType mediumType;
  BlockIndex numDebyeBlocks = 0;

  message( MSG_LOG , 0 , "\nInitialising blocks...\n\n" );

  /* Temporary array for block media. */
  #if USE_AVERAGED_MEDIA
    MediumIndex ***blockArray;
    unsigned long bytes;
    message( MSG_DEBUG1 , 0 , "  Allocating mesh block array\n" );
    blockArray = allocArray( &bytes , sizeof( MediumIndex ) , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );

    for ( int i = gobox[XLO] - 1 ; i <= gobox[XHI] ; i++ )
      for ( int j = gobox[YLO] - 1 ; j <= gobox[YHI] ; j++ )
        for ( int k = gobox[ZLO] - 1 ; k <= gobox[ZHI] ; k++ )
          blockArray[i][j][k] = MT_FREE_SPACE;
  #endif
          
  /* Apply simple media using medium coefficients. */
  DL_FOREACH( blockList , item ) 
  {

    if( !mediumTypeByName( item->mediumName , &mediumType ) )
      assert( 0 );

    isBlockType[mediumType] = true;
    
    switch( mediumType )
    {
    case MT_FREE_SPACE:
    case MT_PEC:
    case MT_SIMPLE:
      offsetBoundingBox( gbbox , item->mbbox , gibox );
      message( MSG_DEBUG3 , 0 , "  Setting SIMPLE block medium#%lu on [%d,%d,%d,%d,%d,%d]/[%d,%d,%d,%d,%d,%d]\n" , 
               item->mediumNumber , item->mbbox[XLO] , item->mbbox[XHI] , item->mbbox[YLO] , item->mbbox[YHI] , 
	       item->mbbox[ZLO] , item->mbbox[ZHI] , gbbox[XLO] , gbbox[XHI] , gbbox[YLO] , gbbox[YHI] , 
	       gbbox[ZLO] , gbbox[ZHI] );
      #if USE_AVERAGED_MEDIA
        setMediumOnMesh( blockArray , gbbox , item->mediumNumber );       
      #else
        setMediumOnGrid( gbbox , item->mediumNumber , item->mask );
      #endif
      break;
    case MT_DEBYE:   
      numDebyeBlocks++;
      offsetBoundingBox( gbbox , item->mbbox , gibox );
      message( MSG_DEBUG3 , 0 , "  Setting DEBYE block medium#%lu on [%d,%d,%d,%d,%d,%d]/[%d,%d,%d,%d,%d,%d]\n" , 
               item->mediumNumber , item->mbbox[XLO] , item->mbbox[XHI] , item->mbbox[YLO] , item->mbbox[YHI] , 
               item->mbbox[ZLO] , item->mbbox[ZHI] , gbbox[XLO] , gbbox[XHI] , gbbox[YLO] , gbbox[YHI] , 
               gbbox[ZLO] , gbbox[ZHI] );
      #if USE_AVERAGED_MEDIA
        setMediumOnMesh( blockArray , gbbox , item->mediumNumber );       
      #else
        setMediumOnGrid( gbbox , item->mediumNumber , item->mask );
      #endif      
      break;
    default:
      assert( 0 );
      break;
    }

  }

  #if USE_AVERAGED_MEDIA
    /* Apply blocks to grid with averaging. */
    applyVoxelsToGrid( blockArray );
    /* Free temporary medium array. */
    deallocArray( blockArray , 3 , numCells[XDIR] , numCells[YDIR] , numCells[ZDIR] );
  #endif

  /* Apply other media. */
  initDebyeBlocks( numDebyeBlocks , blockList );
  
  return;

}

/* Set medium on temporary cell array. */
void setMediumOnMesh( MediumIndex ***blockArray , int gbbox[6] , MediumIndex medium )
{

  /* Care - upper index one less than bbox because these are cells! */
  for( int i = gbbox[XLO]  ; i < gbbox[XHI] ; i++ )
    for( int j = gbbox[YLO]  ; j < gbbox[YHI] ; j++ )
      for( int k = gbbox[ZLO]  ; k < gbbox[ZHI] ; k++ )
        blockArray[i][j][k] = medium;

  return;

}

/* Deallocate blocks. */
void deallocBlocks( void )
{

  BlockItem *item , *tmp;

  message( MSG_DEBUG1 , 0 , "Deallocating blocks...\n" );

  deallocDebyeBlocks();
  
  DL_FOREACH_SAFE( blockList , item , tmp ) 
  {
    DL_DELETE( blockList , item );
    free( item );
  }

  return;

}

/* Report blocks. */
void reportBlocks( void )
{

  BlockItem *item;
  BlockIndex counter = 0;
 
  message( MSG_LOG , 0 , "  Number of blocks: %lu\n" , numBlock );

  DL_FOREACH( blockList , item ) 
  {
    message( MSG_DEBUG3 , 0 , "    Block #%lu: Medium=%s Medium#=%lu BBOX=[%d,%d,%d,%d,%d,%d]\n" , 
             (unsigned long) counter , item->mediumName , (unsigned long) item->mediumNumber ,
             item->mbbox[XLO] , item->mbbox[XHI] , item->mbbox[YLO] , 
             item->mbbox[YHI] , item->mbbox[ZLO] , item->mbbox[ZHI] );
    counter++;
  }

  return;

}

/* Return true if there are blocks of given medium type. */
bool thereAreBlocks( MediumType type )
{
  
  return isBlockType[type];

}

/* Block electric field updates. */
void updateBlocksEfield( void )
{
  
  updateDebyeBlocksEfield();

  return;  

  
}

/* Block magnetic field updates. */
void updateBlocksHfield( void )
{
  
  return;

}

/* Output gnuplot compatible data for blocks. */
void gnuplotBlocks( void )
{

  char blockFileName[] = "gnuplot-block.dat";
  FILE *outputFile;
  BlockItem *item;

  outputFile = fopen( blockFileName , "w" );
  if( !outputFile )
    message( MSG_ERROR , 0 , "*** Error: Failed to open block output file %s\n" , blockFileName );

  gnuplotProblemSize( outputFile , mbox );

  DL_FOREACH( blockList , item ) 
  {
    gnuplotBoundingBox( outputFile , item->mbbox );
  }

  fclose( outputFile );

  return;

}

/* Output gmsh compatible data for blocks. */
void gmshBlocks( void )
{

  int step[3] = { 1 , 1 , 1 };
  BlockItem *item;
  unsigned long entityNumber;
  char name[GMSH_NAME_LENGTH];
  
  DL_FOREACH( blockList , item ) 
  {
    entityNumber = gmshGetEntityNumber();
    snprintf( name , GMSH_NAME_LENGTH - 1 , "MT_%s" , getMediumName( item->mediumNumber ) );
    gmshAddEntity( entityNumber , BB_VOLUME , name , item->mbbox , step );
  }

  return;

}
