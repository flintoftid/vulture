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

#include <assert.h>
#include <math.h>
#include <string.h>

#include "source.h"
#include "utlist.h"
#include "waveform.h"
#include "alloc_array.h"
#include "message.h"
#include "bounding_box.h"
#include "gnuplot.h"
#include "gmsh.h"
#include "grid.h"
#include "medium.h"
#include "physical.h"

/* 
 * Source class. 
 */

typedef struct SourceItem_t {

  SourceIndex number;             // Source number, assigned in order found.

  /* Parameters from mesh. */
  char name[TAG_SIZE];            // Source name.  
  SourceType type;                // Source type.
  int mbbox[6];                   // Bounding box on mesh.
  WaveformIndex waveformNumber;   // Waveform number reference for FDOM.
  CoordAxis polarisation;         // Polarisation.
  CoordAxis normal;               // Surface normal direction.
  FieldComponent field;           // Field component
  bool isSoft;                    // Soft or hard source.
  real size;                      // Amplitude.
  real delay;                     // Delay [s].
  real resistance;                // Internal resistance [ohms].
  
  /* Derived parameters. */
  int gbbox[6];                   // Bounding box on grid.
  int flim[6][6];                 // Field limits for source.

  /* UT list. */

  struct SourceItem_t *prev;
  struct SourceItem_t *next;

  UT_hash_handle hh;              // Source name hash.
  
} SourceItem;

/*
 * Private data. 
 */

/* Source type strings. */
static char SOURCE_TYPE_STR[NUM_SOURCE_TYPES][23] = { "ELECTRIC_FIELD" , 
                                                      "MAGNETIC_FIELD" , 
                                                      "ELEC_CURR_DENSITY" , 
                                                      "MAGN_CURR_DENSITY" , 
                                                      "ELEC_SURF_CURR_DENSITY" , 
                                                      "MAGN_SURF_CURR_DENSITY" , 
                                                      "ELEC_CURRENT" , 
                                                      "MAGN_CURRENT" , 
                                                      "ELEC_CURRENT_MOMENT" , 
                                                      "MAGN_CURRENT_MOMENT" , 
                                                      "VOLTAGE" , 
                                                      "THEVENIN_VOLTAGE" ,
                                                      "NORTON_CURRENT" };

/* Source tags used in mesh file. */
#define NUM_SOURCE_TAGS 84

static char SOURCE_TAG_STR[NUM_SOURCE_TAGS][TAG_SIZE]     
    = { "EX"     , "EY"     , "EZ"     ,
        "HX"     , "HY"     , "HZ"     ,
        "=EX"    , "=EY"    , "=EZ"    ,
        "=HX"    , "=HY"    , "=HZ"    ,
        "JX"     , "JY"     , "JZ"     ,
        "JMX"    , "JMY"    , "JMZ"    ,
        "=JX"    , "=JY"    , "=JZ"    ,
        "=JMX"   , "=JMY"   , "=JMZ"   ,
        "JSXY"   , "JSYZ"   , "JSZX"   ,
        "JSXZ"   , "JSYX"   , "JSZY"   ,
        "JMSXY"  , "JMSYZ"  , "JMSZX"  ,
        "JMSXZ"  , "JMSYX"  , "JMSZY"  ,     
        "=JSXY"  , "=JSYZ"  , "=JSZX"  ,
        "=JSXZ"  , "=JSYX"  , "=JSZY"  ,
        "=JMSXY" , "=JMSYZ" , "=JMSZX" ,
        "=JMSXZ" , "=JMSYX" , "=JMSZY" ,         
        "IX"     , "IY"     , "IZ"     ,
        "IMX"    , "IMY"    , "IMZ"    ,
        "=IX"    , "=IY"    , "=IZ"    ,
        "=IMX"   , "=IMY"   , "=IMZ"   ,
        "IDX"    , "IDY"    , "IDZ"    ,
        "IMDX"   , "IMDY"   , "IMDZ"   ,
        "=IDX"   , "=IDY"   , "=IDZ"   ,
        "=IMDX"  , "=IMDY"  , "=IMDZ"  ,
        "VX"     , "VY"     , "VZ"     ,	
        "=VX"    , "=VY"    , "=VZ"    ,
        "VRX"    , "VRY"    , "VRZ"    ,
        "IGX"    , "IGY"    , "IGZ"    };

/* Source type tag mappings. */
static SourceType sourceTypeMap[NUM_SOURCE_TAGS]          
    = { ST_EFIELD                 , ST_EFIELD                 , ST_EFIELD                 ,
        ST_HFIELD                 , ST_HFIELD                 , ST_HFIELD                 ,
        ST_EFIELD                 , ST_EFIELD                 , ST_EFIELD                 ,
        ST_HFIELD                 , ST_HFIELD                 , ST_HFIELD                 ,
        ST_ELEC_CURR_DENSITY      , ST_ELEC_CURR_DENSITY      , ST_ELEC_CURR_DENSITY      , 
        ST_MAGN_CURR_DENSITY      , ST_MAGN_CURR_DENSITY      , ST_MAGN_CURR_DENSITY      , 
        ST_ELEC_CURR_DENSITY      , ST_ELEC_CURR_DENSITY      , ST_ELEC_CURR_DENSITY      , 
        ST_MAGN_CURR_DENSITY      , ST_MAGN_CURR_DENSITY      , ST_MAGN_CURR_DENSITY      ,
        ST_ELEC_SURF_CURR_DENSITY , ST_ELEC_SURF_CURR_DENSITY , ST_ELEC_SURF_CURR_DENSITY , 
        ST_ELEC_SURF_CURR_DENSITY , ST_ELEC_SURF_CURR_DENSITY , ST_ELEC_SURF_CURR_DENSITY ,         
        ST_MAGN_SURF_CURR_DENSITY , ST_MAGN_SURF_CURR_DENSITY , ST_MAGN_SURF_CURR_DENSITY , 
        ST_MAGN_SURF_CURR_DENSITY , ST_MAGN_SURF_CURR_DENSITY , ST_MAGN_SURF_CURR_DENSITY ,         
        ST_ELEC_SURF_CURR_DENSITY , ST_ELEC_SURF_CURR_DENSITY , ST_ELEC_SURF_CURR_DENSITY , 
        ST_ELEC_SURF_CURR_DENSITY , ST_ELEC_SURF_CURR_DENSITY , ST_ELEC_SURF_CURR_DENSITY ,         
        ST_MAGN_SURF_CURR_DENSITY , ST_MAGN_SURF_CURR_DENSITY , ST_MAGN_SURF_CURR_DENSITY ,
        ST_MAGN_SURF_CURR_DENSITY , ST_MAGN_SURF_CURR_DENSITY , ST_MAGN_SURF_CURR_DENSITY ,         
        ST_ELEC_CURRENT           , ST_ELEC_CURRENT           , ST_ELEC_CURRENT           , 
        ST_MAGN_CURRENT           , ST_MAGN_CURRENT           , ST_MAGN_CURRENT           , 
        ST_ELEC_CURRENT           , ST_ELEC_CURRENT           , ST_ELEC_CURRENT           , 
        ST_MAGN_CURRENT           , ST_MAGN_CURRENT           , ST_MAGN_CURRENT           , 
        ST_ELEC_CURRENT_MOMENT    , ST_ELEC_CURRENT_MOMENT    , ST_ELEC_CURRENT_MOMENT    , 
        ST_MAGN_CURRENT_MOMENT    , ST_MAGN_CURRENT_MOMENT    , ST_MAGN_CURRENT_MOMENT    , 
        ST_ELEC_CURRENT_MOMENT    , ST_ELEC_CURRENT_MOMENT    , ST_ELEC_CURRENT_MOMENT    , 
        ST_MAGN_CURRENT_MOMENT    , ST_MAGN_CURRENT_MOMENT    , ST_MAGN_CURRENT_MOMENT    , 
        ST_VOLTAGE                , ST_VOLTAGE                , ST_VOLTAGE                , 
        ST_VOLTAGE                , ST_VOLTAGE                , ST_VOLTAGE                ,
        ST_THEVENIN_VOLTAGE       , ST_THEVENIN_VOLTAGE       , ST_THEVENIN_VOLTAGE       ,
        ST_NORTON_CURRENT         , ST_NORTON_CURRENT         , ST_NORTON_CURRENT         };

/* Source type tag field component mappings. */
static FieldComponent sourceFieldCompMap[NUM_SOURCE_TAGS] 
    = { EX , EY , EZ ,
        HX , HY , HZ ,
        EX , EY , EZ ,
        HX , HY , HZ ,
        EX , EY , EZ ,
        HX , HY , HZ ,
        EX , EY , EZ ,
        HX , HY , HZ ,
        EX , EY , EZ ,
        EX , EY , EZ ,
        HX , HY , HZ ,
        HX , HY , HZ ,        
        EX , EY , EZ ,
        EX , EY , EZ ,
        HX , HY , HZ ,        
        HX , HY , HZ ,
        EX , EY , EZ ,
        HX , HY , HZ ,      
        EX , EY , EZ ,
        HX , HY , HZ , 
        EX , EY , EZ ,
        HX , HY , HZ ,      
        EX , EY , EZ ,
        HX , HY , HZ ,      
        EX , EY , EZ ,
        EX , EY , EZ ,
        EX , EY , EZ ,
        EX , EY , EZ };                                                            

/* Source type tag polarisation mapping. */
static CoordAxis sourcePolarisationMap[NUM_SOURCE_TAGS]  
    = { XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,      
        XDIR , YDIR , ZDIR ,      
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR , 
        XDIR , YDIR , ZDIR ,      
        XDIR , YDIR , ZDIR ,      
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR ,       
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR , 
        XDIR , YDIR , ZDIR ,
        XDIR , YDIR , ZDIR };

/* Source surface normal mapping. */
static CoordAxis sourceSurfaceNormal[NUM_SOURCE_TAGS]  
    = { CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,
        YDIR         , ZDIR         , XDIR         ,
        ZDIR         , XDIR         , YDIR         ,
        YDIR         , ZDIR         , XDIR         ,
        ZDIR         , XDIR         , YDIR         ,
        YDIR         , ZDIR         , XDIR         ,
        ZDIR         , XDIR         , YDIR         ,
        YDIR         , ZDIR         , XDIR         ,
        ZDIR         , XDIR         , YDIR         ,
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,      
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,      
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED , 
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,      
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,      
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,       
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED , 
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED ,
        CA_UNDEFINED , CA_UNDEFINED , CA_UNDEFINED };
        
/* Source type tag soft/hard mapping. */
static bool sourceSoftnessMap[NUM_SOURCE_TAGS]            
    = { true  , true  , true  ,
        true  , true  , true  ,
        false , false , false ,
        false , false , false ,
        true  , true  , true  ,
        true  , true  , true  ,
        false , false , false ,
        false , false , false ,
        true  , true  , true  ,
        true  , true  , true  ,
        true  , true  , true  ,
        true  , true  , true  ,        
        false , false , false ,
        false , false , false ,
        false , false , false ,
        false , false , false ,        
        true  , true  , true  ,
        true  , true  , true  ,
        false , false , false ,
        false , false , false , 
        true  , true  , true  ,
        true  , true  , true  ,
        false , false , false ,
        false , false , false ,
        true  , true  , true  ,
        false , false , false ,
        true  , true  , true  ,
        true  , true  , true  };

/* Minimum resistance for resistive voltage source. */
#define MIN_RESISTANCE  1e-2

/* Number of sources. */
static SourceIndex numSource = 0;

/* Existance flag for sources of each type, including undefined. */
static bool isSourceType[NUM_SOURCE_TYPES+1] = { false };

/* List of sources. */
static SourceItem *sourceList;

/* Hash of sources using name. */
static SourceItem *sourceHash = NULL;

/* 
 * Private method interfaces. 
 */

void addSource( int mbbox[6] , char name[TAG_SIZE] , SourceType type , CoordAxis pol , CoordAxis norm , FieldComponent field , 
                bool isSoft , real size , real delay , real resist, WaveformIndex waveformNumber );
void getElectricSourceSize( int gbbox[6] , CoordAxis direction , real *length , real *area , real side[3] );
void getMagneticSourceSize( int gbbox[6] , CoordAxis direction , real *length , real *area , real side[3] );
bool isSource( char *name , SourceIndex *number );

/*
 * Method Implementations.
 */

/* Add source to lists. */
void addSource( int mbbox[6] , char name[TAG_SIZE] , SourceType type , CoordAxis pol , CoordAxis norm , FieldComponent field , 
                bool isSoft , real size , real delay , real resist , WaveformIndex waveformNumber )
{

  SourceItem *item = NULL;

  if( numSource == MAX_SOURCE )
    message( MSG_ERROR , 0 , "*** Error: Maximum number of sources exceeded!\n" );

  item = (SourceItem *) malloc( sizeof( SourceItem ) );
  if( !item )
    message( MSG_ERROR , 0 , "*** Error: Failed to allocate source item" );

  strncpy( item->name , name , TAG_SIZE );
  item->type = type;
  item->number = numSource;
  for( int i = XLO ; i <= ZHI ; i++ ) item->mbbox[i] = mbbox[i];
  item->polarisation = pol;
  item->normal = norm;
  item->field = field;
  item->isSoft = isSoft;
  item->resistance = resist;
  item->waveformNumber = waveformNumber;
  item->size = size;
  item->delay = delay;

  /* Add to list. */
  DL_APPEND( sourceList , item );
  HASH_ADD_STR( sourceHash , name , item );
  numSource++;
  isSourceType[type] = true;
  isSourceType[ST_UNDEFINED] = true;
  
  return;

}

/* Parse field excitations. */
bool parseEX( char *line )
{

  int numScanned = 0;
  char name[TAG_SIZE];
  char typeStr[TAG_SIZE] = "";
  char waveformName[TAG_SIZE] = "";
  char resistName[TAG_SIZE] = "";
  int mbbox[6];
  SourceType type = ST_EFIELD;
  bool foundType = false;
  CoordAxis pol = XDIR;
  FieldComponent field = EX;
  double size = 1.0;
  double delay = 0.0;
  double resist = 0.0;
  bool isSoft = true;
  WaveformIndex waveformNumber;
  SourceIndex sourceNumber;
  CoordAxis norm = CA_UNDEFINED;
  
  numScanned = sscanf( line , "%d %d %d %d %d %d %31s %31s %31s" , 
                       &mbbox[XLO] , &mbbox[XHI] , &mbbox[YLO] , &mbbox[YHI] , &mbbox[ZLO] , &mbbox[ZHI] , 
                       name , typeStr , waveformName );

  if( numScanned < 9 )
    return false;  

  /* Check source is not already defined. */
  if( isSource( name , &sourceNumber ) )
  {
    message( MSG_LOG , 0 , "  Source %s already defined\n" , name );
    return false;
  }

  /* Find field component. */
  for( int i = 0 ; i < NUM_SOURCE_TAGS ; i++ )
  {
    if( strncmp( typeStr , SOURCE_TAG_STR[i] , TAG_SIZE ) == 0 )
    {
      type = sourceTypeMap[i];
      pol = sourcePolarisationMap[i];
      norm = sourceSurfaceNormal[i];
      field = sourceFieldCompMap[i]; 
      isSoft = sourceSoftnessMap[i];
      foundType = true;
    }
  }

  if( !foundType )
  {
    message( MSG_LOG , 0 , "  Invalid source type: %s\n" , typeStr );
    return false;
  }

  /* Validate bounding box. */ 
  if( !bboxIsNormal( mbbox ) )
  {
    message( MSG_LOG , 0 , "  Bounding box is abnormal:\n" );
    return false;
  }
  else if( !bboxIsWithin( mbbox , mbox ) )
  {
    message( MSG_LOG , 0 , "  Bounding box is outside mesh:\n" );
    return false;
  }

  /* Type dependent assignment and validation of paramters. */
  switch( type )
  {
  case ST_THEVENIN_VOLTAGE:
  case ST_NORTON_CURRENT:
    
    numScanned = sscanf( line , "%d %d %d %d %d %d %31s %31s %31s %lf %lf %lf" , 
                         &mbbox[XLO] , &mbbox[XHI] , &mbbox[YLO] , &mbbox[YHI] , &mbbox[ZLO] , &mbbox[ZHI] , 
                         name , typeStr , waveformName , &resist , &size , &delay );

    /* Currently lumped source must be a single edge. */
    if( bboxType( mbbox ) != BB_LINE || !bboxIsElemental( mbbox ) )
    {
      message( MSG_LOG , 0 , "  Lumped source bounding box must be an edge!\n" );
      return false;
    }

    if(  numScanned >= 12 && delay < 0.0 )
      message( MSG_WARN , 0 , "Source delay negative:\n" , delay );
    
    // Add placeholder for lumped resistance material. 
    sprintf( resistName , "__VR_RS_%d__" , numSource + 1 );
    addMedium( resistName , MT_SIMPLE , 1.0 , 0.0 , 1.0 , 0 , NULL , NULL , "" );

    break;
  default:
    
    numScanned = sscanf( line , "%d %d %d %d %d %d %31s %31s %31s %lf %lf" , 
                         &mbbox[XLO] , &mbbox[XHI] , &mbbox[YLO] , &mbbox[YHI] , &mbbox[ZLO] , &mbbox[ZHI] , 
                         name , typeStr , waveformName , &size , &delay );

    if(  numScanned >= 11 && delay < 0.0 )
      message( MSG_WARN , 0 , "Source delay negative:\n" );

    break;  
  }

  /* Check waveform exists. */
  if( !isWaveform( waveformName , &waveformNumber ) )
  {
    message( MSG_LOG , 0 , "  Waveform %s not defined in source card\n" , waveformName );
    return false;
  }

  addSource( mbbox, name , type , pol , norm , field , isSoft , (real)size , (real)delay , (real)resist , waveformNumber );

  return true;

}

/* Initialise sources. */
/* Depends: initGrid, initWaveforms, initMedia */
void initSources( void )
{

  SourceItem *item;
  char resistName[TAG_SIZE] = "";
  MediumIndex resistIndex;
  bool includeBoundary[6] = {  true ,  true ,  true ,  true ,  true ,  true };
  real length;
  real area;
  real side[3];
  real sigma;
 
  message( MSG_LOG , 0 , "\nInitialising sources...\n\n" );

  message( MSG_DEBUG1 , 0 , "  Allocating source array\n" );

  DL_FOREACH( sourceList , item ) 
  {
    
    offsetBoundingBox( item->gbbox , item->mbbox , gibox );
    setFieldLimits( item->gbbox , item->flim , includeBoundary );

    switch( item->type )
    {      
    case ST_EFIELD:
    case ST_HFIELD:   
    case ST_ELEC_CURR_DENSITY:
    case ST_MAGN_CURR_DENSITY:
      /* No-op. */
      break;
    case ST_ELEC_SURF_CURR_DENSITY:
      /* Convert surface density current source to volume current density source. */ 
      getElectricSourceSize( item->gbbox , item->polarisation , &length , &area , side );
      item->size = item->size / side[item->normal];
      item->type = ST_ELEC_CURR_DENSITY;
      break;
    case ST_MAGN_SURF_CURR_DENSITY:
      /* Convert surface density current source to volume current density source. */
      getMagneticSourceSize( item->gbbox , item->polarisation , &length , &area , side );
      item->size = item->size / side[item->normal];
      item->type = ST_MAGN_CURR_DENSITY;
      break;
    case ST_ELEC_CURRENT:
      /* Convert current source to current density source. */
      getElectricSourceSize( item->gbbox , item->polarisation , &length , &area , side );
      item->size = item->size / area;
      item->type = ST_ELEC_CURR_DENSITY;
      break;
    case ST_MAGN_CURRENT:
      /* Convert current source to current density source. */
      getMagneticSourceSize( item->gbbox , item->polarisation , &length , &area , side );
      item->size = item->size / area;
      item->type = ST_MAGN_CURR_DENSITY;
      break;
    case ST_ELEC_CURRENT_MOMENT:
      /* Convert current source to current density source. */
      getElectricSourceSize( item->gbbox , item->polarisation , &length , &area , side );
      item->size = item->size / area / length;
      item->type = ST_ELEC_CURR_DENSITY;
      break;
    case ST_MAGN_CURRENT_MOMENT:
      /* Convert current source to current density source. */
      getMagneticSourceSize( item->gbbox , item->polarisation , &length , &area , side );
      item->size = item->size / area / length;
      item->type = ST_MAGN_CURR_DENSITY;
      break;
    case ST_VOLTAGE:
      /* Convert voltage source to field source. */
      getElectricSourceSize( item->gbbox , item->polarisation , &length , &area , side );
      item->size = -item->size / length;
      item->type = ST_EFIELD;       
      break;
    case ST_THEVENIN_VOLTAGE:
    case ST_NORTON_CURRENT:
      if( item->type == ST_NORTON_CURRENT )
      {
        /* Convert Norton source to Thevenin source. */
        item->size = item->resistance * item->size;
        item->type = ST_THEVENIN_VOLTAGE;
      }
      /* Convert voltage source to current density source and medium parameters. */
      sprintf( resistName , "__VR_RS_%d__" , item->number + 1 );
      if( !isMedium( resistName , &resistIndex ) )
        assert( 0 );
      getElectricSourceSize( item->gbbox , item->polarisation , &length , &area , side );
      if( item->resistance > MIN_RESISTANCE )
      {
        /* Resistive voltage source. */
        item->size = item->size / area / item->resistance;      
        sigma = length / ( area * item->resistance );
        updateSimpleMedium( resistIndex , 1.0 , sigma , 1.0 );
        setMediumOnGrid( item->gbbox , resistIndex , FACE_MASK_ALL );
        item->type = ST_ELEC_CURR_DENSITY;
      }
      else
      {
        /* Ideal zero resistance voltage source. */
        item->size = -0.5 * item->size / length;
        setMediumOnGrid( item->gbbox , MT_PEC , FACE_MASK_ALL );
        item->type = ST_EFIELD;       
      }
      break;
    default:
      break;
    }

    message( MSG_DEBUG3 , 0 , "  Setting %s source \"%s\" on [%d,%d,%d,%d,%d,%d]/[%d,%d,%d,%d,%d,%d]: pol=%s, soft=%s, size=%g, delay=%g, resist=%g\n" , 
               SOURCE_TYPE_STR[item->type] , item->name ,
               item->mbbox[XLO] , item->mbbox[XHI] , item->mbbox[YLO] , item->mbbox[YHI] , 
	       item->mbbox[ZLO] , item->mbbox[ZHI] , item->gbbox[XLO] , item->gbbox[XHI] , 
	       item->gbbox[YLO] , item->gbbox[YHI] , item->gbbox[ZLO] , item->gbbox[ZHI] , 
	       AXIS[item->polarisation] , BOOL[item->isSoft] , item->size, item->delay , item->resistance );
	        
    message( MSG_DEBUG3 , 0 , "    EX BBOX=[%d,%d,%d,%d,%d,%d]\n" , item->flim[EX][XLO] , item->flim[EX][XHI] , item->flim[EX][YLO] , 
	                                                            item->flim[EX][YHI] , item->flim[EX][ZLO] , item->flim[EX][ZHI] );
    message( MSG_DEBUG3 , 0 , "    EY BBOX=[%d,%d,%d,%d,%d,%d]\n" , item->flim[EY][XLO] , item->flim[EY][XHI] , item->flim[EY][YLO] , 
	                                                            item->flim[EY][YHI] , item->flim[EY][ZLO] , item->flim[EY][ZHI] );
    message( MSG_DEBUG3 , 0 , "    EZ BBOX=[%d,%d,%d,%d,%d,%d]\n" , item->flim[EZ][XLO] , item->flim[EZ][XHI] , item->flim[EZ][YLO] , 
	                                                            item->flim[EZ][YHI] , item->flim[EZ][ZLO] , item->flim[EZ][ZHI] );
    message( MSG_DEBUG3 , 0 , "    HX BBOX=[%d,%d,%d,%d,%d,%d]\n" , item->flim[HX][XLO] , item->flim[HX][XHI] , item->flim[HX][YLO] , 
	                                                            item->flim[HX][YHI] , item->flim[HX][ZLO] , item->flim[HX][ZHI] );
    message( MSG_DEBUG3 , 0 , "    HY BBOX=[%d,%d,%d,%d,%d,%d]\n" , item->flim[HY][XLO] , item->flim[HY][XHI] , item->flim[HY][YLO] , 
	                                                            item->flim[HY][YHI] , item->flim[HY][ZLO] , item->flim[HY][ZHI] );
    message( MSG_DEBUG3 , 0 , "    HZ BBOX=[%d,%d,%d,%d,%d,%d]\n" , item->flim[HZ][XLO] , item->flim[HZ][XHI] , item->flim[HZ][YLO] , 
	                                                            item->flim[HZ][YHI] , item->flim[HZ][ZLO] , item->flim[HZ][ZHI] );
	       
  }

  return;

}

/* Update electric field and voltage sources. */
void updateSourcesEfield( real timeE )
{

  real source;
  int i , j , k;
  FieldComponent field;
  SourceItem *item;
  real isSoft;

  DL_FOREACH( sourceList , item ) 
  {
    field = item->field;
    source = item->size * getWaveformValue( timeE , item->waveformNumber , item->delay );
    isSoft = (real)item->isSoft;

    switch( item->type )
    {
    case ST_EFIELD:
      switch( field )
      {
      case EX:
        for( i = item->flim[field][XLO] ; i <= item->flim[field][XHI] ; i++ )
          for( j = item->flim[field][YLO] ; j <= item->flim[field][YHI] ; j++ )
            for( k = item->flim[field][ZLO] ; k <= item->flim[field][ZHI] ; k++ )
              Ex[i][j][k] = isSoft * Ex[i][j][k] + SCALE_Ex( source , i );
        break;
      case EY:
        for( i = item->flim[field][XLO] ; i <= item->flim[field][XHI] ; i++ )
          for( j = item->flim[field][YLO] ; j <= item->flim[field][YHI] ; j++ )
            for( k = item->flim[field][ZLO] ; k <= item->flim[field][ZHI] ; k++ )
              Ey[i][j][k] = isSoft * Ey[i][j][k] + SCALE_Ey( source , j );
        break;
      case EZ:
        for( i = item->flim[field][XLO] ; i <= item->flim[field][XHI] ; i++ )
          for( j = item->flim[field][YLO] ; j <= item->flim[field][YHI] ; j++ )
            for( k = item->flim[field][ZLO] ; k <= item->flim[field][ZHI] ; k++ )
              Ez[i][j][k] = isSoft * Ez[i][j][k] + SCALE_Ez( source , k );
        break;
      default:
        break;
      } // switch( field )
      break;
    case ST_ELEC_CURR_DENSITY:
      switch( field )
      {
      case EX:
        for( i = item->flim[field][XLO] ; i <= item->flim[field][XHI] ; i++ )
          for( j = item->flim[field][YLO] ; j <= item->flim[field][YHI] ; j++ )
            for( k = item->flim[field][ZLO] ; k <= item->flim[field][ZHI] ; k++ )
              Ex[i][j][k] = isSoft * Ex[i][j][k] - BETA_EX(i,j,k) * SCALE_Jx( source , i );
        break;
      case EY:
        for( i = item->flim[field][XLO] ; i <= item->flim[field][XHI] ; i++ )
          for( j = item->flim[field][YLO] ; j <= item->flim[field][YHI] ; j++ )
            for( k = item->flim[field][ZLO] ; k <= item->flim[field][ZHI] ; k++ )
              Ey[i][j][k] = isSoft * Ey[i][j][k] - BETA_EY(i,j,k) * SCALE_Jy( source , j );
        break;
      case EZ:
        for( i = item->flim[field][XLO] ; i <= item->flim[field][XHI] ; i++ )
          for( j = item->flim[field][YLO] ; j <= item->flim[field][YHI] ; j++ )
            for( k = item->flim[field][ZLO] ; k <= item->flim[field][ZHI] ; k++ )
              Ez[i][j][k] = isSoft * Ez[i][j][k] - BETA_EZ(i,j,k) * SCALE_Jz( source , k );
        break;
      default:
        break;
      } // switch( field )
      break;  
    default:
      break;
    } // switch( item->type )
  } // DL_FOREACH( sourceList , item ) 

  return;

}

/* Update magnetic field sources. */
void updateSourcesHfield( real timeH )
{

  real source;
  int i , j , k;
  FieldComponent field;
  SourceItem *item;
  real isSoft;

  DL_FOREACH( sourceList , item ) 
  {
    field = item->field;
    source = item->size * getWaveformValue( timeH , item->waveformNumber , item->delay );
    isSoft = (real)item->isSoft;

    switch( item->type )
    {
    case ST_HFIELD:
      switch( field )
      {
      case HX:
        for( i = item->flim[field][XLO] ; i <= item->flim[field][XHI] ; i++ )
          for( j = item->flim[field][YLO] ; j <= item->flim[field][YHI] ; j++ )
            for( k = item->flim[field][ZLO] ; k <= item->flim[field][ZHI] ; k++ )
              Hx[i][j][k] = isSoft * Hx[i][j][k] + SCALE_Hx( source , i );
        break;
      case HY:
        for( i = item->flim[field][XLO] ; i <= item->flim[field][XHI] ; i++ )
          for( j = item->flim[field][YLO] ; j <= item->flim[field][YHI] ; j++ )
            for( k = item->flim[field][ZLO] ; k <= item->flim[field][ZHI] ; k++ )
              Hy[i][j][k] = isSoft * Hy[i][j][k] + SCALE_Hy( source , j );
        break;
      case HZ:
        for( i = item->flim[field][XLO] ; i <= item->flim[field][XHI] ; i++ )
          for( j = item->flim[field][YLO] ; j <= item->flim[field][YHI] ; j++ )
            for( k = item->flim[field][ZLO] ; k <= item->flim[field][ZHI] ; k++ )
              Hz[i][j][k] = isSoft * Hz[i][j][k] + SCALE_Hz( source , k );
        break;
      default:
        break;
      } // switch( field )
      break;
  case ST_MAGN_CURR_DENSITY:
      switch( field )
      {
      case HX:
        for( i = item->flim[field][XLO] ; i <= item->flim[field][XHI] ; i++ )
          for( j = item->flim[field][YLO] ; j <= item->flim[field][YHI] ; j++ )
            for( k = item->flim[field][ZLO] ; k <= item->flim[field][ZHI] ; k++ )
              Hx[i][j][k] = isSoft * Hx[i][j][k] - GAMMA_HX(i,j,k) * SCALE_JMx( source , i );
        break;
      case HY:
        for( i = item->flim[field][XLO] ; i <= item->flim[field][XHI] ; i++ )
          for( j = item->flim[field][YLO] ; j <= item->flim[field][YHI] ; j++ )
            for( k = item->flim[field][ZLO] ; k <= item->flim[field][ZHI] ; k++ )
              Hy[i][j][k] = isSoft * Hy[i][j][k] - GAMMA_HY(i,j,k) * SCALE_JMy( source , j );
        break;
      case HZ:
        for( i = item->flim[field][XLO] ; i <= item->flim[field][XHI] ; i++ )
          for( j = item->flim[field][YLO] ; j <= item->flim[field][YHI] ; j++ )
            for( k = item->flim[field][ZLO] ; k <= item->flim[field][ZHI] ; k++ )
              Hz[i][j][k] = isSoft * Hz[i][j][k] - GAMMA_HZ(i,j,k) * SCALE_JMz( source , k );     
        break;
      default:
        break;
      }  // switch( field )
      break;      
    default:
      break;
    } // switch( item->type )
  } // DL_FOREACH( sourceList , item ) 

  return;

}

/* Report sources. */
void reportSources( void )
{

  SourceItem *item;

  message( MSG_LOG , 0 , "  Number of sources: %lu\n" , numSource );

  DL_FOREACH( sourceList , item ) 
  {
    message( MSG_DEBUG3 , 0 , "    Source \"%s\" (#%lu): Waveform#=%d Type=%s Pol=%s Soft=%s BBOX=[%d,%d,%d,%d,%d,%d] size=%e delay=%e Z=%e\n" , 
             item->name , (unsigned long) item->number , item->waveformNumber , SOURCE_TYPE_STR[item->type] , 
             AXIS[item->polarisation] , BOOL[item->isSoft] ,
             item->mbbox[XLO] , item->mbbox[XHI] , item->mbbox[YLO] , 
             item->mbbox[YHI] , item->mbbox[ZLO] , item->mbbox[ZHI] ,
             item->size , item->delay , item->resistance );
  }

  return;

}

/* Return true if there are sources of given type. */
bool thereAreSources( SourceType type )
{
  
  return  isSourceType[type];

}

/* Deallocate sources. */
void deallocSources( void )
{

  SourceItem *item , *tmp;

  message( MSG_DEBUG1 , 0 , "Deallocating sources...\n" );

  /* Free source name hash and the sources. */
  HASH_ITER( hh , sourceHash , item , tmp )
  {
    HASH_DELETE( hh , sourceHash , item );
    free( item );
  }

  return;

}

/* Get source number from name. */
bool isSource( char *name , SourceIndex *number )
{

  SourceItem *item;

  HASH_FIND_STR( sourceHash , name , item );
  if( item )
  {
    *number = item->number;
    return true;
  }
  else
  {
    *number = 0;
    return false;
  }
  
}

/* Output gnuplot compatiable plot data for sources. */
void gnuplotSources( void )
{

  char sourceFileName[] = "gnuplot-source.dat";
  FILE *outputFile; 
  SourceItem *item;

  outputFile = fopen( sourceFileName , "w" );
  if( !outputFile )
    message( MSG_ERROR , 0 , "*** Error: Failed to open source output file %s\n" , sourceFileName );

  gnuplotProblemSize( outputFile , mbox );

  DL_FOREACH( sourceList , item ) 
  {
    gnuplotBoundingBox( outputFile , item->mbbox );
    gnuplotBoundingBoxArrow( outputFile , item->mbbox , item->field );
  }

  fclose (outputFile);

  return;

}

/* Output gmsh compatiable plot data for sources. */
void gmshSources( void )
{

  //int bbox[6];
  int step[3] = { 1 , 1 , 1 };
  SourceItem *item;
  unsigned long entityNumber;
  char name[GMSH_NAME_LENGTH];
  BoundingBoxType mboxType;
  //CoordAxis mboxDir;
  
  DL_FOREACH( sourceList , item ) 
  {

    snprintf( name , GMSH_NAME_LENGTH - 1 , "EX_%s" , item->name );

    /* Find type and direction of bounding box */
    mboxType = bboxType( item->mbbox );
    //mboxDir = bboxDirection( item->mbbox );
  
    // [FIXME]
    // Simple:
    //   bbox type BB_POINT - render points
    //   bbox type BB_LINE - render edges
    //   bbox type BB_SURFACE- render surface
    //   bbox type BB_VOLUME - render volume or outer surfaces?
    // Better:
    //   field source - render edges according to polarisation
    //   current source - faces according to polarisation
    switch( mboxType )
    {
      case BB_POINT:
        entityNumber = gmshGetEntityNumber();
        gmshAddEntity( entityNumber , mboxType , name , item->mbbox , step );        
        break;
      case BB_LINE:
        entityNumber = gmshGetEntityNumber();
        gmshAddEntity( entityNumber , mboxType , name , item->mbbox , step ); 
        break;
      case BB_SURFACE:
        entityNumber = gmshGetEntityNumber();
        gmshAddEntity( entityNumber , mboxType , name , item->mbbox , step ); 
        break;
      case BB_VOLUME:
        entityNumber = gmshGetEntityNumber();
        gmshAddEntity( entityNumber , mboxType , name , item->mbbox , step ); 
        break;
      default:
        assert( 0 );
        break;
    }
    
    //entityNumber = gmshGetEntityNumber();
    //getFaceOfBoundingBox( bbox , item->mbbox , ZLO );
    //gmshAddEntity( entityNumber , BB_SURFACE , name , bbox , step );

    //entityNumber = gmshGetEntityNumber();
    //getFaceOfBoundingBox( bbox , item->mbbox , XLO );
    //gmshAddEntity( entityNumber , BB_SURFACE , name , bbox , step );

    //entityNumber = gmshGetEntityNumber();
    //getFaceOfBoundingBox( bbox , item->mbbox , YLO );
    //gmshAddEntity( entityNumber , BB_SURFACE , name , bbox , step );

    //entityNumber = gmshGetEntityNumber();
    //getFaceOfBoundingBox( bbox , item->mbbox , ZHI );
    //gmshAddEntity( entityNumber , BB_SURFACE , name , bbox , step );

    //entityNumber = gmshGetEntityNumber();
    //getFaceOfBoundingBox( bbox , item->mbbox , XHI );
    //gmshAddEntity( entityNumber , BB_SURFACE , name , bbox , step );

    //entityNumber = gmshGetEntityNumber();
    //getFaceOfBoundingBox( bbox , item->mbbox , YHI );
    //gmshAddEntity( entityNumber , BB_SURFACE , name , bbox , step );

  }

  return;

}

/* Get length and area of distributed electric source. */
void getElectricSourceSize( int gbbox[6] , CoordAxis direction , real *length , real *area , real side[3] )
{

  *length = 0.0;
  *area = 0.0;
  side[XDIR] = 0.0;
  side[YDIR] = 0.0;
  side[ZDIR] = 0.0;
  
  switch( direction )
  {
  case XDIR:
    for( int i = gbbox[XLO] ; i <= gbbox[XHI] - 1 ; i++ ) *length += dex[i];
    for( int j = gbbox[YLO] ; j <= gbbox[YHI]     ; j++ ) side[YDIR] += dhy[j];
    for( int k = gbbox[ZLO] ; k <= gbbox[ZHI]     ; k++ ) side[ZDIR] += dhz[k];
    *area = side[YDIR] * side[ZDIR];
    break;
  case YDIR:
    for( int i = gbbox[XLO] ; i <= gbbox[XHI]     ; i++ ) side[XDIR] += dhx[i];
    for( int j = gbbox[YLO] ; j <= gbbox[YHI] - 1 ; j++ ) *length += dey[j];
    for( int k = gbbox[ZLO] ; k <= gbbox[ZHI]     ; k++ ) side[ZDIR] += dhz[k];
    *area = side[ZDIR] * side[XDIR];
    break;
  case ZDIR:
    for( int i = gbbox[XLO] ; i <= gbbox[XHI]     ; i++ ) side[XDIR] += dhx[i];
    for( int j = gbbox[YLO] ; j <= gbbox[YHI]     ; j++ ) side[YDIR] += dhy[j];
    for( int k = gbbox[ZLO] ; k <= gbbox[ZHI] - 1 ; k++ ) *length += dez[k];
    *area = side[XDIR] * side[YDIR];
    break;
  default:
    assert( false );
    break;
  }
  
  return;

}

/* Get length and area of distributed magnetic source. */
void getMagneticSourceSize( int gbbox[6] , CoordAxis direction , real *length , real *area , real side[3] )
{

  *length = 0.0;
  *area = 0.0;
  side[XDIR] = 0.0;
  side[YDIR] = 0.0;
  side[ZDIR] = 0.0;

  switch( direction )
  {
  case XDIR:
    for( int i = gbbox[XLO] ; i <= gbbox[XHI]     ; i++ ) *length += dhx[i];
    for( int j = gbbox[YLO] ; j <= gbbox[YHI] - 1 ; j++ ) side[YDIR] += dey[j];
    for( int k = gbbox[ZLO] ; k <= gbbox[ZHI] - 1 ; k++ ) side[ZDIR] += dez[k];
    *area = side[YDIR] * side[ZDIR];
    break;
  case YDIR:
    for( int i = gbbox[XLO] ; i <= gbbox[XHI] - 1 ; i++ ) side[XDIR] += dex[i];
    for( int j = gbbox[YLO] ; j <= gbbox[YHI]     ; j++ ) *length += dhy[j];
    for( int k = gbbox[ZLO] ; k <= gbbox[ZHI] - 1 ; k++ ) side[ZDIR] += dez[k];
    *area = side[ZDIR] * side[XDIR];
    break;
  case ZDIR:
    for( int i = gbbox[XLO] ; i <= gbbox[XHI] - 1 ; i++ ) side[XDIR] += dex[i];
    for( int j = gbbox[YLO] ; j <= gbbox[YHI] - 1 ; j++ ) side[YDIR] += dey[j];
    for( int k = gbbox[ZLO] ; k <= gbbox[ZHI]     ; k++ ) *length += dhz[k];
    *area = side[XDIR] * side[YDIR];
    break;
  default:
    assert( false );
    break;
  }

  return;

}
