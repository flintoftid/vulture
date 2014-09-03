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

#include "medium.h"
#include "grid.h"
#include "physical.h"
#include "alloc_array.h"
#include "message.h"
#include "bounding_box.h"
#include "memory.h"

/* 
 * Global data.
 */ 

/* Array of pointer to media indexed by medium number/index. */
MediumItem **mediumArray = NULL;

/* 
 * Private data. 
 */

/* Medium type strings. */
char MEDIUM_TYPE_STR[NUM_MEDIUM_TYPES][TAG_SIZE] = { "FREE_SPACE" , "PEC" , "SIMPLE" , "DEBYE" };

/* Number of medium types. */
static MediumIndex numMedium = 0;

/* Existance flag for media of each type, including undefined. */
static bool isMediumType[NUM_MEDIUM_TYPES+1] = { false };

/* List of medium types. */
static MediumItem *mediumList = NULL;

/* Hash of medium type by name.*/
static MediumItem *mediumHash = NULL;

/* Hash of medium type by number.*/
static MediumItem *mediumNumberHash = NULL;

/* 
 * Private method interfaces. 
 */
void readDebyeParameters( char fileName[PATH_SIZE] , real *eps_inf , real *sigma , real *mu_r , 
                          int *numPoles , double complex **residues , double complex **poles );

/*
 * Method Implementations.
 */

/* Parse material type. */
bool parseMT( char *line )
{

  int numScanned = 0;
  char name[TAG_SIZE] = "";
  char typeStr[TAG_SIZE] = "";
  double eps_r = 1.0;
  double sigma = 0.0;
  double mu_r = 1.0;
  MediumType type = MT_FREE_SPACE;
  MediumIndex number;
  double residues[3] = { 0.0 };
  double poles[3] = { 0.0 };
  int numPoles = 0;
  char fileName[PATH_SIZE] = "";
  bool foundType = false;

  numScanned = sscanf( line , "%31s %31s" , name , typeStr );
  if( numScanned < 2 )
    return false;  

  /* Check medium not already defined. */
  if( isMedium( name , &number ) )
  {
    message( MSG_LOG , 0 , "  Medium %s already defined\n" , name );
    return false;
  }

  /* Find type. */
  for( int medium = 0 ; medium < NUM_MEDIUM_TYPES ; medium++ )
    if( strncmp( typeStr , MEDIUM_TYPE_STR[medium] , TAG_SIZE ) == 0 )
    {
      type = (MediumType)medium;      
      foundType = true;
    }

  if( !foundType )
  {
    message( MSG_LOG , 0 , "  Invalid medium: %s\n" , type );
    return false;
  }

  /* Validate parameters. */
  switch( type )
  {
  case MT_FREE_SPACE:
  case MT_PEC:
    break;
  case MT_SIMPLE:
    numScanned = sscanf( line , "%31s %31s %lf %lf %lf" , name , typeStr , &eps_r , &sigma , &mu_r );
    if( numScanned >= 3 && eps_r < 1.0 )
    {
      message( MSG_LOG , 0 , "  Relative permittivity must be >= 1.0\n" );      
      return false;
    }
    if( numScanned >= 4 && sigma < 0.0 )
    {
      message( MSG_LOG , 0 , "  Relative permittivity must be >= 0.0\n" );      
      return false;
    }
    if( numScanned >= 5 && mu_r < 1.0 )
    {
      message( MSG_LOG , 0 , "  Relative permeability must be >= 0.0\n" );      
      return false;
    }
    break;
  case MT_DEBYE:
    /* First try for parameters. */
    numScanned = sscanf( line , "%31s %31s %lf %lf %lf %lf %lf %lf %lf %lf %lf" , name , typeStr , &eps_r , &sigma , &mu_r , 
                         &residues[0] , &poles[0] , &residues[1] , &poles[1] , &residues[2] , &poles[2] );
    if( numScanned >= 3 && eps_r < 1.0 )
    {
      message( MSG_LOG , 0 , "  Relative permittivity must be >= 1.0\n" );      
      return false;
    }
    if( numScanned >= 4 && sigma < 0.0 )
    {
      message( MSG_LOG , 0 , "  Relative permittivity must be >= 0.0\n" );      
      return false;
    }
    if( numScanned >= 5 && mu_r < 1.0 )
    {
      message( MSG_LOG , 0 , "  Relative permeability must be >= 0.0\n" );      
      return false;
    }    
    if( numScanned == 7 || numScanned == 9 || numScanned == 11 )
    {
      numPoles = ( numScanned - 5 ) / 2; 
    }
    else
    {
      /* See if file name given. */
      numScanned = sscanf( line , "%31s %31s \"%1023[^\"]\"" , name , typeStr , fileName );
      if( numScanned != 3 )
      {
        message( MSG_LOG , 0 , "  Unable to parse MT directive for DEBYE type\n" );      
        return false;
      }        
    }
    break;
  default:
    assert( 0 );
    break;
  }
  
  addMedium( name , type , eps_r , sigma , mu_r , numPoles , residues , poles , fileName );

  return true;

}

/* Get medium number from name. */
bool isMedium( char *name , MediumIndex *number )
{
  
  MediumItem *item = NULL;

  HASH_FIND_STR( mediumHash , name , item );
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

/* Get medium type from name. */
bool mediumTypeByName( char *name , MediumType *type )
{
  
  MediumItem *item = NULL;

  HASH_FIND_STR( mediumHash , name , item );
  if( item )
  {
    *type = item->type;
    return true;
  }
  else
  {
    *type = MT_UNDEFINED;
    return false;
  }
  
}

/* Get medium type by number. */
MediumType getMediumType( MediumIndex number )
{

  MediumItem *item;

  HASH_FIND( hhint , mediumNumberHash , &number , sizeof( number ) , item );
  if( !item)
    assert( 0 );

  return item->type;

}

/* Get medium pointer. */
MediumItem *getMedium( MediumIndex number )
{

  MediumItem *item;

  HASH_FIND( hhint , mediumNumberHash , &number , sizeof( number ) , item );
  if( !item)
    assert( 0 );

  return item;

}

/* Add medium to lists. */
void addMedium( char *name , MediumType type , real eps_r , real sigma , real mu_r , 
                int numPoles , double residues[3] , double poles[3] , char fileName[PATH_SIZE] )
{

  MediumItem *item = NULL;

  if( numMedium == MAX_MEDIA )
    message( MSG_ERROR , 0 , "*** Error: Maximum number of media exceeded!\n" );
 
  item = (MediumItem *) malloc( sizeof( MediumItem ) );
  if( !item )
    message( MSG_ERROR , 0 , "*** Error: Failed to allocate medium!\n" );

  strncpy( item->name , name , TAG_SIZE );
  item->number = numMedium;
  item->type = type;
  item->eps_r = eps_r;
  item->sigma = sigma;
  item->mu_r = mu_r;
  strncpy( item->fileName , fileName , PATH_SIZE );     
      
  if( type == MT_DEBYE )
  {
    if( strlen( fileName ) == 0 )
    {
      item->numPoles = numPoles;
      item->residues = malloc( numPoles * sizeof( double complex ) );
      item->poles = malloc( numPoles * sizeof( double complex ) );
      for( int poleIdx = 0 ; poleIdx < numPoles ; poleIdx++ )
      {
        item->residues[poleIdx] = (double complex)residues[poleIdx];
        item->poles[poleIdx] = (double complex)poles[poleIdx];
      }
    }
    else
    {
      readDebyeParameters( item->fileName , &item->eps_r , &item->sigma , &item->mu_r , &item->numPoles , &item->residues , &item->poles );
    }
  }
  else
  {
    item->numPoles = 0;    
    item->residues = NULL;
    item->poles = NULL;
  }
  
  /* Add to list and hash. */
  DL_APPEND( mediumList , item );
  HASH_ADD_STR( mediumHash , name , item );
  HASH_ADD( hhint , mediumNumberHash , number , sizeof( numMedium ) , item );
  numMedium++;
  isMediumType[type] = true;
  isMediumType[MT_UNDEFINED] = true; 

  return;

}

/* Read Debye parameters from external file. */
void readDebyeParameters( char fileName[PATH_SIZE] , real *eps_inf , real *sigma , real *mu_r , 
                          int *numPoles , double complex **residues , double complex **poles )
{
  
  FILE *filePointer = NULL;
  double tmp1 = 0.0;
  double tmp2 = 0.0;
  double tmp3 = 0.0;
  double tmp4 = 0.0;
  int numScanned = 0;

  message( MSG_LOG , 0 , "    Reading Debye parameters from file %s\n" , fileName );
       
  /* Open file. */
  filePointer = fopen( fileName , "r" );
  if( filePointer == NULL ) 
     message( MSG_ERROR , 0 , "  ***Error: Cannot open Debye parameter file %s\n" , fileName );
  
  numScanned = fscanf( filePointer , "%d %lf %lf %lf" , numPoles , &tmp1 , &tmp2 , &tmp3 );
  if( numScanned != 4 )
    message( MSG_ERROR , 0 , "  *** Error: Failed to read first line of Debye parameter file %s\n" , fileName );
  
  if( *numPoles < 0 )
    message( MSG_ERROR , 0 , "  *** Error: Number of poles (%d) in file %s must be >=0\n" , *numPoles , fileName );

  if( tmp1 < 1.0 )
    message( MSG_ERROR , 0 , "  *** Error: High frequency rel. permittivity (%e) in file %s must be >=1\n" , tmp1 , fileName );
  else
    *eps_inf = (real)tmp1;

  if( tmp2 < 0.0 )
    message( MSG_ERROR , 0 , "  *** Error: Conductivity (%e) in file %s must be >=0\n" , tmp2 , fileName );    
  else
    *sigma = (real)tmp2;
  
  if( tmp3 < 1.0 )
    message( MSG_ERROR , 0 , "  *** Error: Rel. permeability (%e) in file %s must be >=0\n" , tmp3 , fileName );    
  else
    *mu_r = (real)tmp3;

  *residues = malloc( *numPoles * sizeof( double complex ) );
  if( !*residues )
    message( MSG_ERROR , 0 , "  *** Failed to allociate memory for residues array for %s\n" , fileName ); 
  *poles = malloc( *numPoles * sizeof( double complex ) );
  if( !*poles )
    message( MSG_ERROR , 0 , "  *** Failed to allocate memory for poles array for %s\n" , fileName ); 

  for( int poleIdx = 0 ; poleIdx < *numPoles ; poleIdx++ )
  {
    numScanned = fscanf( filePointer , "%lf %lf %lf %lf" , &tmp1 , &tmp2 , &tmp3 , &tmp4 );
    if( numScanned != 4 )
      message( MSG_ERROR , 0 , "  *** Error: Failed to read parameters from Debye parameter file %s\n" , fileName ); 
    (*residues)[poleIdx] = tmp1 + I * tmp2;
    if( tmp3 > 0.0 )
      message( MSG_ERROR , 0 , "  *** Error: Unstable pole (%e+j%e) in file %s\n" , tmp3, tmp4 , fileName );
    else
      (*poles)[poleIdx] = tmp3 + I * tmp4;
  }

  fclose( filePointer );
  
  return;

}

/* Initialise the material coefficient arrays. */
/* Depends: initGrid */
void initMedia( void )
{

  MediumIndex mediumNumber = 0;
  MediumItem *item = NULL;
  real dt = 0.0;
  unsigned long bytes = 0;

  message( MSG_LOG , 0 , "\nInitialising media...\n\n" );

  dt = getGridTimeStep();

  /* Allocate material parameter arrays. */
  message( MSG_DEBUG1 , 0 , "  Allocating media array\n" );
  mediumArray = allocArray( &bytes , sizeof( MediumItem * ) , 1 , numMedium );
  memory.media += bytes;

  /* Determine update coefficient for each material. */
  mediumNumber = 0;
  DL_FOREACH( mediumList , item ) 
  {
    mediumArray[mediumNumber] = item;
    if( item->type == MT_DEBYE )
    {
      item->dalpha = malloc( item->numPoles * sizeof( double complex ) );
      item->dbeta = malloc( item->numPoles * sizeof( double complex ) ); 
    }
    calcCoeffFromParam( &(item->alpha) , &(item->beta) , &(item->gamma) , item->dalpha , item->dbeta ,
                        dt , item->eps_r , item->sigma , item->mu_r , 
                        item->numPoles , item->residues , item->poles );
    message( MSG_DEBUG3 , 0 , "    Medium#=%lu: eps_r=%g, sigma=%g, mu_r=%g, alpha=%g, beta=%g, gamma=%g, npole=%d\n" , mediumNumber,
	     item->eps_r , item->sigma , item->mu_r ,  item->alpha , item->beta , item->gamma , item->numPoles );
    for( int poleIdx = 0 ; poleIdx < item->numPoles ; poleIdx++ )
      message( MSG_DEBUG3 , 0 , "      pole#=%d Re(dalpha)=%e Im(dalpha)=%e Re(dbeta)=%e Im(dbeta)=%e\n" ,
               poleIdx , creal( item->dalpha[poleIdx] ) , cimag( item->dalpha[poleIdx] ) , 
               creal( item->dbeta[poleIdx] ) , cimag( item->dbeta[poleIdx] ) );    
    mediumNumber++;
  }

  /* Force medium zero to free-space. */
  mediumArray[MT_FREE_SPACE]->alpha = 1.0;
  mediumArray[MT_FREE_SPACE]->beta  = dt / eps0;
  mediumArray[MT_FREE_SPACE]->gamma = dt / mu0;
  mediumArray[MT_FREE_SPACE]->dalpha = NULL;
  mediumArray[MT_FREE_SPACE]->dbeta = NULL;  
  message( MSG_DEBUG3 , 0 , "    Forcing medium#=%lu: alpha=%g, beta=%g, gamma=%g\n" , MT_FREE_SPACE ,
	   mediumArray[MT_FREE_SPACE]->alpha , mediumArray[MT_FREE_SPACE]->beta , mediumArray[MT_FREE_SPACE]->gamma );
	     
  /* Force medium one to PEC. */
  /* Important that beta is exactly zero for PEC! */
  mediumArray[MT_PEC]->alpha = -1.0;
  mediumArray[MT_PEC]->beta  = 0.0;
  mediumArray[MT_PEC]->gamma = dt / mu0;
  mediumArray[MT_FREE_SPACE]->dalpha = NULL;
  mediumArray[MT_FREE_SPACE]->dbeta = NULL;  
  message( MSG_DEBUG3 , 0 , "    Forcing medium#=%lu: alpha=%g, beta=%g, gamma=%g\n" , MT_PEC ,
	   mediumArray[MT_PEC]->alpha , mediumArray[MT_PEC]->beta , mediumArray[MT_PEC]->gamma );
	   
  return;

}

/* Get update ceofficients of a medium. */
void getSimpleMediumCoefficients( real *alpha , real *beta , real *gamma , MediumIndex medium )
{

  *alpha = mediumArray[medium]->alpha;
  *beta = mediumArray[medium]->beta;
  *gamma = mediumArray[medium]->gamma;

  return;

}

/* Get medium paramters. */
void getSimpleMediumParameters( real *eps_r , real *sigma , real *mu_r , MediumIndex medium )
{

  *eps_r = mediumArray[medium]->eps_r;
  *sigma = mediumArray[medium]->sigma;
  *mu_r = mediumArray[medium]->mu_r;

  return;

}

/* Get pointer to medium name by number. */
char *getMediumName( MediumIndex number )
{

  MediumItem *item;

  HASH_FIND( hhint , mediumNumberHash , &number , sizeof( number ) , item );
  if( !item)
    assert( 0 );

  return item->name;

}

/* Determine Debye medium coefficients from parameters. */
void calcCoeffFromParam( real *alpha , real *beta , real *gamma , double complex *dalpha , double complex *dbeta ,
                         real dt , real eps_r , real sigma , real mu_r , 
                         int numPoles , double complex residues[] , double complex poles[] )
{

  double sum = 0.0;

  for( int poleIdx = 0 ; poleIdx < numPoles ; poleIdx++ )
  {
    dalpha[poleIdx] = ( 1 + 0.5 * dt * poles[poleIdx] ) / ( 1 - 0.5 * dt * poles[poleIdx] );
    dbeta[poleIdx] = eps0 * residues[poleIdx] / ( 1 - 0.5 * dt * poles[poleIdx] );
    sum = sum + creal( dbeta[poleIdx] );
  }

  *alpha = ( 2.0 * eps_r * eps0 + 2.0 * sum * dt - dt * sigma ) / 
           ( 2.0 * eps_r * eps0 + 2.0 * sum * dt + dt * sigma );
  *beta  = ( 2.0 * dt ) / 
           ( 2.0 * eps_r * eps0 + 2.0 * sum * dt + dt * sigma );
  *gamma = dt / ( mu_r * mu0 );

  return;

}

/* Deallocate media. */
void deallocMedia( void )
{

  MediumItem *item , *tmp;

  message( MSG_DEBUG1 , 0 , "Deallocating media...\n" );
 
  deallocArray( mediumArray , 1 , numMedium );

  /* Free the number hash. */
  HASH_ITER( hhint , mediumNumberHash , item , tmp )
  {
    HASH_DELETE( hhint , mediumNumberHash , item );
  }
  
  HASH_ITER( hh , mediumHash , item , tmp )
  {
    if( item->type == MT_DEBYE )
    {
      free( item->residues );
      free( item->poles );
      free( item->dalpha );
      free( item->dbeta ); 
    }
    HASH_DEL(  mediumHash , item );
    free( item );
  }

  return;

}

/* Report media. */
void reportMedia( void )
{

  MediumItem *item;

  message( MSG_LOG , 0 , "  Number of media: %lu\n" , (unsigned long) numMedium );

  DL_FOREACH( mediumList , item ) 
  {
    message( MSG_DEBUG3 , 0 , "    Medium #%lu: Name=%s Type=%s eps_r=%e sigma=%e mu_r=%e npoles=%d\n" , 
             (unsigned long) item->number , item->name , MEDIUM_TYPE_STR[item->type] , 
             item->eps_r , item->sigma , item->mu_r , item->numPoles );
    for( int poleIdx = 0 ; poleIdx < item->numPoles ; poleIdx++ )
      message( MSG_DEBUG3 , 0 , "      pole#=%d Re(residue)=%e Im(residue)=%e Re(pole)=%e Im(pole)=%e\n" ,
               poleIdx , creal( item->residues[poleIdx] ) , cimag( item->residues[poleIdx] ) , 
               creal( item->poles[poleIdx] ) , cimag( item->poles[poleIdx] ) );
  }

  return;

}

/* Return true if there are medium types of given type. */
bool thereAreMedia( MediumType type )
{
    
  return isMediumType[type];

}

/* Update medium parameters for simple mdeium. */
void updateSimpleMedium( MediumIndex index , real eps_r , real sigma , real mu_r )
{

  real dt;
  
  /* Only valid after medium array is initialised. */
  assert( mediumArray );
  
  mediumArray[index]->eps_r = eps_r;
  mediumArray[index]->sigma = sigma;
  mediumArray[index]->mu_r = mu_r;
  
  dt = getGridTimeStep();
  
  calcCoeffFromParam( &(mediumArray[index]->alpha) , &(mediumArray[index]->beta) , &(mediumArray[index]->gamma) , 
                      NULL , NULL , dt , eps_r , sigma , mu_r , 0 , NULL , NULL );

  return; 
  
}
