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

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <limits.h>

#include "gmsh.h"
#include "utlist.h"
#include "alloc_array.h"
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

/* Number of nodes in element types. */
static int ELEMENT_NUM_NODES[NUM_ELEMENT_TYPES] = { 1 , 2 , 4 , 8 };

/* Gmsh element type codes. */
static int ELEM_TYPE_CODE[NUM_ELEMENT_TYPES] = { 15 , 1 , 3 , 5 };

/* Gmsh element type dimensions. */
static int ELEM_TYPE_DIMS[NUM_ELEMENT_TYPES] = { 0 , 1 , 2 , 3 };

/*
 * Element class.
 */
#define MAX_ELEM_NODES 8

typedef struct ElementItem_t {

  ElementType   type;                   // Element type.
  unsigned long nodes[MAX_ELEM_NODES];  // Element node numbers.
  unsigned long groupNumber;            // Element's physical group number
  unsigned long entityNumber;           // Element's mesh entity number.

  /* UT list/hash. */

  struct ElementItem_t *prev;
  struct ElementItem_t *next;

} ElementItem;

/* 
 * Group class. 
 */
typedef struct GroupItem_t {

  char name[GMSH_NAME_LENGTH]; // Group name.
  unsigned long number;        // Group number, > 0.
  int dimension;               // Group dimension - 1, 2 or 3.

  /* UT list/hash. */

  struct GroupItem_t *prev;
  struct GroupItem_t *next;

  UT_hash_handle hh;      // Group name hash.

} GroupItem;

/* 
 * Private data. 
 */

/* Node map array - holds node number of nodes which are used. 
   Unused nodes are set to sentinel value NODE_UNUSED. */
#define NODE_UNUSED UINT_MAX
unsigned long ***nodeMap = NULL;

/* Number of elements. */
#define MAX_ELEMENT UINT_MAX
static unsigned long numElement = 0;

/* List of elements. */
static ElementItem *elementList = NULL;

/* Number of groups. */
#define MAX_GROUP UINT_MAX
static unsigned long numGroup = 0;

/* Hash of groups by name. */
static GroupItem *groupHash = NULL;

/* Next entity number. */
#define MAX_ENTITY UINT_MAX
static unsigned long nextEntityNumber = 1;

/* Number of mesh nodes in each direction. */
static int numMeshNodes[3];

/* Whether to use indices or physical units for node coordinates. */
static bool isPhysicalUnits = true;

/* Private functions. */
unsigned long encodeNodeNumber( int i , int j , int k );
void gmshWrite( void );
bool isGroup( char *name , unsigned long *number , int *dimension );
unsigned long addGroup( char *name , int dimension );
unsigned long addGroupIfNew( char *name , int dimension );
void addElement( ElementType type , unsigned long nodes[MAX_ELEM_NODES] , unsigned long groupNumber , unsigned long entityNumber );


/* Output mesh is gmsh format. */
void gmshMesh( bool isPhysUnits , bool isExternalSurfaces )
{

  GroupItem *groupItem = NULL;
  GroupItem *tmpGroupItem = NULL;
  ElementItem *elemItem = NULL;
  ElementItem *tmpElemItem = NULL;    
  int numMeshCells[3] = { 0 , 0 , 0 };
  unsigned long bytes = 0;
 
  isPhysicalUnits = isPhysUnits;

  /* Get mesh extents. */
  getGridNumCells( numMeshCells );
  
  /* Number of mesh nodes/lines in each direction. */ 
  numMeshNodes[XDIR] = numMeshCells[XDIR] + 1;
  numMeshNodes[YDIR] = numMeshCells[YDIR] + 1;
  numMeshNodes[ZDIR] = numMeshCells[ZDIR] + 1;
  
  /* Create node map and initialise to unused. */
  nodeMap = allocArray( &bytes , sizeof( unsigned long ) , 3 , numMeshNodes[XDIR] , numMeshNodes[YDIR] , numMeshNodes[ZDIR] );
  for( int i = 0 ; i < numMeshNodes[XDIR] ; i++ )
    for( int j = 0 ; j < numMeshNodes[YDIR] ; j++ )
      for( int k = 0 ; k < numMeshNodes[ZDIR] ; k++ )
         nodeMap[i][j][k] = NODE_UNUSED;


  /* Add mesh objects to gmsh mesh. */
  if( isExternalSurfaces ) 
    gmshExternalSurfaces();
  gmshInternalSurfaces();
  gmshBlocks();
  gmshLines();
  gmshSources();
  gmshPlaneWaves();
  gmshObservers();
  /* [FIXME] is this useful? */
  /* gmshGridLines(); */

  /* Write out mesh. */
  gmshWrite();
       
  /* Deallocate element list. */
  DL_FOREACH_SAFE( elementList , elemItem , tmpElemItem ) 
  {
    DL_DELETE( elementList , elemItem );
    free( elemItem );
  }

  /* Deallocate group name hash. */
  HASH_ITER( hh , groupHash , groupItem , tmpGroupItem )
  {
    HASH_DELETE( hh , groupHash , groupItem );
    free( groupItem );
  }
  
  /* Deallocate node map. */
  deallocArray( nodeMap , 3 , numMeshNodes[XDIR] , numMeshNodes[YDIR] , numMeshNodes[ZDIR] );
  
  return;

}

/* Write gmsh file. */
void gmshWrite( void )
{

  char mshFileName[] = "mesh.msh";
  FILE *mshFile= NULL;
  unsigned long numNodes = 0;
  GroupItem *groupItem = NULL;
  GroupItem *tmpGroupItem = NULL;
  ElementItem *elemItem = NULL;
  unsigned long groupIdx = 0;
  unsigned long elementIdx = 0;
  real nodeCoords[3] = { 0.0 , 0.0 , 0.0 };
  int nodeIndices[3] = { 0 , 0 , 0 };

  mshFile = fopen( mshFileName , "w" );
  if( !mshFile )
    message( MSG_ERROR , 0 , "*** Error: Failed to open output msh file %s\n" , mshFileName );
  
  /* MeshFormat. */
  fprintf( mshFile , "$MeshFormat\n" );
  fprintf( mshFile ,"%.1f %d %d\n" , 2.2 , 0 , 8 );
  fprintf( mshFile , "$EndMeshFormat\n" );

  message( MSG_LOG , 0 , "Wrote mesh format" );
  
  /* PhysicalNames. */
  if( numGroup > 0 )
  {
    fprintf( mshFile , "$PhysicalNames\n" );
    fprintf( mshFile , "%lu\n" , numGroup );
    HASH_ITER( hh , groupHash , groupItem , tmpGroupItem )
    {
      fprintf( mshFile ,"%d %lu \"%s\"\n" , groupItem->dimension , groupItem->number ,  groupItem->name );
      groupIdx++;
    }
    fprintf( mshFile , "$EndPhysicalNames\n" );
  }

  /* Defence work. */
  assert( groupIdx == numGroup );

  message( MSG_LOG , 0 , "Wrote physical names" );

  /* Count number of used nodes. */
  for( int i = 0 ; i < numMeshNodes[XDIR] ; i++ )
    for( int j = 0 ; j < numMeshNodes[YDIR] ; j++ )
      for( int k = 0 ; k < numMeshNodes[ZDIR] ; k++ )
         if( nodeMap[i][j][k] != NODE_UNUSED ) numNodes++;
      
  /* Output mesh nodes. */
  fprintf( mshFile , "$Nodes\n" );
  fprintf( mshFile , "%lu\n" , numNodes );
  for( int i = 0 ; i < numMeshNodes[XDIR] ; i++ )
    for( int j = 0 ; j < numMeshNodes[YDIR] ; j++ )
      for( int k = 0 ; k < numMeshNodes[ZDIR] ; k++ )
      {
         if( nodeMap[i][j][k] != NODE_UNUSED )
         {
           nodeIndices[XDIR] = i;
           nodeIndices[YDIR] = j;
           nodeIndices[ZDIR] = k;
           if( isPhysicalUnits )
           {
             getMeshNodeCoords( nodeCoords , nodeIndices );
           }
           else
           {
             nodeCoords[XDIR] = (real)nodeIndices[XDIR];
             nodeCoords[YDIR] = (real)nodeIndices[YDIR];
             nodeCoords[ZDIR] = (real)nodeIndices[ZDIR];
           }
           fprintf( mshFile , "%lu %e %e %e\n" , nodeMap[i][j][k] , 
                    nodeCoords[XDIR] , nodeCoords[YDIR] , nodeCoords[ZDIR] );
         }
      }
  fprintf( mshFile , "$EndNodes\n" );

  message( MSG_LOG , 0 , "Wrote nodes" );

  /* Elements. */
  fprintf( mshFile , "$Elements\n" );
  fprintf( mshFile , "%lu\n" , numElement );
  DL_FOREACH( elementList , elemItem ) 
  {
    /* gmsh element number must be > 0 so add one here! */
    fprintf( mshFile , "%lu %d %d " , elementIdx + 1UL , ELEM_TYPE_CODE[elemItem->type] , 2 );
    fprintf( mshFile , "%lu %lu" ,  elemItem->groupNumber , elemItem->entityNumber );
    for( int idx = 0 ; idx < ELEMENT_NUM_NODES[elemItem->type] ; idx++ )
      fprintf( mshFile , " %lu" , elemItem->nodes[idx] );
    fprintf( mshFile , "\n" );
    elementIdx++;
  }
  fprintf( mshFile , "$EndElements\n" );

  /* Defence work. */
  assert( elementIdx == numElement );
  
  message( MSG_LOG , 0 , "Wrote elements" );
  
  fclose( mshFile );

  message( MSG_LOG , 0 , "Closed file" );
    
  return;

}

/* Get group number from name if its exists. */
bool isGroup( char *name , unsigned long *number , int *dimension )
{

  GroupItem *item;

  HASH_FIND_STR( groupHash , name , item );
  if( item )
  {
    *number = item->number;
    *dimension = item->dimension;
    return true;
  }
  else
  {
    *number = 0;
    return false;
  }
  
}

/* Add group to hash. */
unsigned long addGroup( char *name , int dimension )
{

  GroupItem *item = NULL;

  if( numGroup == MAX_GROUP )
    message( MSG_ERROR , 0 , "Maximum number of groups exceeded!\n" );

  item = (GroupItem *) malloc( sizeof( GroupItem ) );
  if( !item )
    message( MSG_ERROR , 0 , "Failed to allocate group %s\n" , name );

  strncpy( item->name , name , GMSH_NAME_LENGTH );
  item->dimension = dimension;
  /* gmsh physical group number must be > 0 so add one here! */
  item->number = numGroup + 1;

  /* Add to hash. */
  HASH_ADD_STR( groupHash , name , item );
  numGroup++;
  
  return item->number;

}

/* Add group to hash if it is not already there. Return group number and dimension. */
unsigned long addGroupIfNew( char *name , int dimension )
{

  unsigned long groupNumber;
  int groupDimension;

  if( isGroup( name , &groupNumber , &groupDimension ) )
  {
    /* gmsh groups must contain elements of the same dimensionality! */
    if( groupDimension != dimension )
      message( MSG_ERROR , 0 , "incompatible dimension %d for physical group %s" , dimension , name );
  }
  else
  {
    groupNumber = addGroup( name , dimension );
  }  

  return groupNumber;

}

/* Add element to list. */
void addElement( ElementType type , unsigned long nodes[MAX_ELEM_NODES] , unsigned long groupNumber , unsigned long entityNumber )
{

  ElementItem *item = NULL;

  if( numElement == MAX_ELEMENT )
    message( MSG_ERROR , 0 , "*** Error: Maximum number of elements exceeded!\n" );

  item = (ElementItem *) malloc( sizeof( ElementItem ) );
  if( !item )
    message( MSG_ERROR , 0 , "*** Error: Failed to allocate element item" );

  item->type = type;
  item->groupNumber = groupNumber;
  item->entityNumber = entityNumber;  
  
  for( int idx = 0 ; idx < ELEMENT_NUM_NODES[type] ; idx++ ) 
    item->nodes[idx] = nodes[idx];

  /* Add to list. */
  DL_APPEND( elementList , item );
  numElement++;
  
  return;

}

/* Add entity to mesh. */
void gmshAddEntity( unsigned long entityNumber , BoundingBoxType typeToAdd , char *name , int mbbox[6] , int step[3] )
{

  ElementType elementType;
  unsigned long groupNumber = 0;
  BoundingBoxType mboxType;
  CoordAxis mboxDir;
  unsigned long nodes[MAX_ELEM_NODES];
  int i;
  int j;
  int k;

  message( MSG_DEBUG3 , 0 , "  Adding entity: number=%lu type=%s name=%s mbbox=[%d,%d,%d,%d,%d,%d]\n" ,
           entityNumber , BBOX_STR[typeToAdd] , name , mbbox[XLO] , mbbox[XHI] , mbbox[YLO] , mbbox[YHI] , mbbox[ZLO] , mbbox[ZHI] );
   
  /* Find type and direction of bounding box */
  mboxType = bboxType( mbbox );
  mboxDir = bboxDirection( mbbox );

  message( MSG_DEBUG3 , 0 , "    Mesh bbox: type=%s dir=%s\n" ,
           BBOX_STR[mboxType] , AXIS[mboxDir] );

  /* Add compatible mesh elements for each bounding box to mesh. */
  switch( mboxType )
  {
    case BB_POINT:
      if( typeToAdd == BB_POINT )
      {
        /* Just add the point. */
        elementType = ET_NODE1;
        groupNumber = addGroupIfNew( name , ELEM_TYPE_DIMS[elementType] );
        i = mbbox[XLO];
        j = mbbox[YLO];
        k = mbbox[ZLO];
        nodes[0] = nodeMap[i][j][k] = encodeNodeNumber( i , j , k );
        addElement( elementType , nodes , groupNumber , entityNumber );
      }
      else
      {
        assert( 0 );
      }
      break;
    case BB_LINE:
      if( typeToAdd == BB_POINT )
      {
        /* Add all points on line. */
        elementType = ET_NODE1;
        groupNumber = addGroupIfNew( name , ELEM_TYPE_DIMS[elementType] );
        for( i = mbbox[XLO] ; i <= mbbox[XHI]; i += step[XDIR] )
          for( j = mbbox[YLO] ; j <= mbbox[YHI]; j += step[YDIR] )
            for( k = mbbox[ZLO] ; k <= mbbox[ZHI]; k += step[ZDIR] )
            {
              nodes[0] = nodeMap[i][j][k] = encodeNodeNumber( i , j , k );
              addElement( ET_NODE1 , nodes , groupNumber , entityNumber );
            } 
      }
      else if( typeToAdd == BB_LINE )
      {
        /* Add all edges on line. */
        switch( mboxDir )
        {
          case XDIR:
            elementType = ET_BAR2;
            groupNumber = addGroupIfNew( name , ELEM_TYPE_DIMS[elementType] );
            j = mbbox[YLO];
            k = mbbox[ZLO];
            for( i = mbbox[XLO] ; i < mbbox[XHI]; i++ )
            {
              nodes[0] = nodeMap[i][j][k] = encodeNodeNumber( i , j , k );
              nodes[1] = nodeMap[i+1][j][k] = encodeNodeNumber( i + 1 , j , k );
              addElement( elementType , nodes , groupNumber , entityNumber );
            }
            break;
          case YDIR:
            elementType = ET_BAR2;
            groupNumber = addGroupIfNew( name , ELEM_TYPE_DIMS[elementType] );
            i = mbbox[XLO];
            k = mbbox[ZLO];
            for( j = mbbox[YLO] ; j < mbbox[YHI]; j++ )
            {
              nodes[0] = nodeMap[i][j][k] = encodeNodeNumber( i , j , k );
              nodes[1] = nodeMap[i][j+1][k] = encodeNodeNumber( i , j + 1 , k );
              addElement( elementType , nodes , groupNumber , entityNumber );
            }            
            break;
          case ZDIR:
            elementType = ET_BAR2;
            groupNumber = addGroupIfNew( name , ELEM_TYPE_DIMS[elementType] );
            i = mbbox[XLO];
            j = mbbox[YLO];
            for( k = mbbox[ZLO] ; k < mbbox[ZHI]; k++ )
            {
              nodes[0] = nodeMap[i][j][k] = encodeNodeNumber( i , j , k );
              nodes[1] = nodeMap[i][j][k+1] = encodeNodeNumber( i , j , k + 1 );
              addElement( elementType , nodes , groupNumber , entityNumber );
            }            
            break;
          default:
            assert( 0 );
            break;
        }
      }
      else
      {
        assert( 0 );
      }      
      break;      
    case BB_SURFACE:
      if( typeToAdd == BB_POINT )
      {
        /* Add all points on surface. */
        elementType = ET_NODE1;
        groupNumber = addGroupIfNew( name , ELEM_TYPE_DIMS[elementType] );
        for( i = mbbox[XLO] ; i <= mbbox[XHI]; i += step[XDIR] )
          for( j = mbbox[YLO] ; j <= mbbox[YHI]; j += step[YDIR] )
            for( k = mbbox[ZLO] ; k <= mbbox[ZHI]; k += step[ZDIR] )
            {
              nodes[0] = nodeMap[i][j][k] = encodeNodeNumber( i , j , k );
              addElement( elementType , nodes , groupNumber , entityNumber );
            } 
      }
      else if( typeToAdd == BB_SURFACE )
      {
        /* Add all faces on surface. */
        switch( mboxDir )
        {
          case XDIR:
            elementType = ET_QUAD4;
            groupNumber = addGroupIfNew( name , ELEM_TYPE_DIMS[elementType] );
            i = mbbox[XLO];
            for( j = mbbox[YLO] ; j < mbbox[YHI]; j++ )
              for( k = mbbox[ZLO] ; k < mbbox[ZHI]; k++ )
              {
                nodes[0] = nodeMap[i][j][k]     = encodeNodeNumber( i , j     , k     );
                nodes[1] = nodeMap[i][j+1][k]   = encodeNodeNumber( i , j + 1 , k     );
                nodes[2] = nodeMap[i][j+1][k+1] = encodeNodeNumber( i , j + 1 , k + 1 );
                nodes[3] = nodeMap[i][j][k+1]   = encodeNodeNumber( i , j     , k + 1 );
                addElement( elementType , nodes , groupNumber , entityNumber );
              }
            break;
          case YDIR:
            elementType = ET_QUAD4;
            groupNumber = addGroupIfNew( name , ELEM_TYPE_DIMS[elementType] );
            j = mbbox[YLO];
            for( i = mbbox[XLO] ; i < mbbox[XHI]; i++ )
              for( k = mbbox[ZLO] ; k < mbbox[ZHI]; k++ )
              {
                nodes[0] = nodeMap[i][j][k]     = encodeNodeNumber( i     , j , k     );
                nodes[1] = nodeMap[i+1][j][k]   = encodeNodeNumber( i + 1 , j , k     );
                nodes[2] = nodeMap[i+1][j][k+1] = encodeNodeNumber( i + 1 , j , k + 1 );
                nodes[3] = nodeMap[i][j][k+1]   = encodeNodeNumber( i     , j , k + 1 );
                addElement( elementType , nodes , groupNumber , entityNumber );
              }            
            break;
          case ZDIR:
            elementType = ET_QUAD4;
            groupNumber = addGroupIfNew( name , ELEM_TYPE_DIMS[elementType] );
            k = mbbox[ZLO];
            for( i = mbbox[XLO] ; i < mbbox[XHI]; i++ )
              for( j = mbbox[YLO] ; j < mbbox[YHI]; j++ )
              {
                nodes[0] = nodeMap[i][j][k]     = encodeNodeNumber( i     , j     , k );
                nodes[1] = nodeMap[i+1][j][k]   = encodeNodeNumber( i + 1 , j     , k );
                nodes[2] = nodeMap[i+1][j+1][k] = encodeNodeNumber( i + 1 , j + 1 , k );
                nodes[3] = nodeMap[i][j+1][k]   = encodeNodeNumber( i     , j + 1 , k );
                addElement( elementType , nodes , groupNumber , entityNumber );
              }            
            break;
          default:
            assert( 0 );
            break;
        }        
      }
      else
      {
        assert( 0 );
      } 
      break; 
    case BB_VOLUME:
      if( typeToAdd == BB_POINT )
      {
        /* Add all points in volume. */
        elementType = ET_NODE1;
        groupNumber = addGroupIfNew( name , ELEM_TYPE_DIMS[elementType] );
        for( i = mbbox[XLO] ; i <= mbbox[XHI]; i += step[XDIR] )
          for( j = mbbox[YLO] ; j <= mbbox[YHI]; j += step[YDIR] )
            for( k = mbbox[ZLO] ; k <= mbbox[ZHI]; k += step[ZDIR] )
            {
              nodes[0] = nodeMap[i][j][k] = encodeNodeNumber( i , j , k );
              addElement( elementType , nodes , groupNumber , entityNumber );
            } 
      } 
      else if( typeToAdd == BB_VOLUME )
      {
        /* Add all cubes in volume. */
        elementType = ET_HEX8;
        groupNumber = addGroupIfNew( name , ELEM_TYPE_DIMS[elementType] );
        for( i = mbbox[XLO] ; i < mbbox[XHI]; i++ )
          for( j = mbbox[YLO] ; j < mbbox[YHI]; j++ )
            for( k = mbbox[ZLO] ; k < mbbox[ZHI]; k++ )
            {
              nodes[0] = nodeMap[i][j][k] = encodeNodeNumber( i , j , k );
              nodes[1] = nodeMap[i+1][j][k] = encodeNodeNumber( i + 1 , j , k );
              nodes[2] = nodeMap[i+1][j+1][k] = encodeNodeNumber( i + 1 , j + 1 , k );
              nodes[3] = nodeMap[i][j+1][k] = encodeNodeNumber( i , j + 1 , k );
              nodes[4] = nodeMap[i][j][k+1] = encodeNodeNumber( i , j , k + 1 );
              nodes[5] = nodeMap[i+1][j][k+1] = encodeNodeNumber( i + 1 , j , k + 1 );
              nodes[6] = nodeMap[i+1][j+1][k+1] = encodeNodeNumber( i + 1 , j + 1 , k + 1 );
              nodes[7] = nodeMap[i][j+1][k+1] = encodeNodeNumber( i , j + 1 , k + 1 );
              addElement( elementType , nodes , groupNumber , entityNumber );
            } 
      }
      else
      {
        assert( 0 );
      }      
      break; 
    default:
      assert( 0 );
      break;
  }
 
  return;

}

/* Generate implicit node number for node (i,j,k). */ 
unsigned long encodeNodeNumber( int i , int j , int k )
{
  /* Use the AMELET implicit encoding scheme. Gmah node numbers must be > 0. */
  return (i+1) + j * numMeshNodes[XDIR] + k * numMeshNodes[XDIR] * numMeshNodes[YDIR]; 
  //return 1 + i + j * numMeshNodes[XDIR] + k * numMeshNodes[XDIR] * numMeshNodes[ZDIR]; 
  
}

/* Get next available entity number. */ 
unsigned long gmshGetEntityNumber( void )
{
  if( nextEntityNumber == MAX_ENTITY )
    message( MSG_ERROR , 0 , "*** Error: Maximum number of entities exceeded!\n" );

  /* Note this increments after returning value! */
  return  nextEntityNumber++;
}
