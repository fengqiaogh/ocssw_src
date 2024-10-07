#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "genutils.h"
#include "sensorInfo.h"

void get_cld_tbl_( int *inst_id, int *fname_id, char *fname, int *status, 
  long len1 )
 /************************************************************************
  get_cld_tbl_ - get the cloud table name that is asked for

  Returns - none

  int *    inst_id     I      instrument id, see oel_util/libpiutils/
                                sensorDefs.h
  int *    fname_id    I      table file name id #, see fid_names below
  char *   fname       O      table file name found
  int *    status      O      a returned status - 0 is good
  int      len1        I~     length of the fname string 

  The file with the names is in each cloud sensor subdirectory with the 
  name cloud_table_names.txt

  W. Robinson, SAIC, 14 Sep 2021
  ************************************************************************/
  {
  #define NFILES 25
  static int get_entered = 0, files_present[NFILES];
  static char file_store[NFILES][563];
  
  char *fid_names[] = { "FAST_TRANS_COEFF_DRY", "FAST_TRANS_COEFF_OZO",
    "FAST_TRANS_COEFF_WCO", "FAST_TRANS_COEFF_WTL", "FAST_TRANS_COEFF_WTS",
    "SNOW_ALBEDO", "ECOLOGY_MAP", "GLOBAL_EMISS", "ALBEDO_CLIM",
    "TRANSMITTANCE_COEFFS", "SINGLE_SCAT", 
    "SCAT_ICE", "SCAT_ICE_03_M_S", "SCAT_ICE_07_M_S", "SCAT_ICE_15_M_S", 
    "SCAT_ICE_SD_03_M_S", "SCAT_ICE_SD_07_M_S", "SCAT_ICE_SD_15_M_S", 
    "SCAT_WATER", "SCAT_WATER_03_M_S", "SCAT_WATER_07_M_S", 
    "SCAT_WATER_15_M_S", "SCAT_WATER_SD_03_M_S", "SCAT_WATER_SD_07_M_S", 
    "SCAT_WATER_SD_15_M_S" };
  char str[563];
  char *id, *floc;
  char file[500], *filep;
  FILE *rd_stream;
  int ilin, i, ifound;

  *status = 0;
  
  for( ilin = 0; ilin < len1; ilin++ )
    fname[ilin] = ' ';
  // if we haven't gotten the file names, get those
  if( get_entered == 0 )
    {
    get_entered = 1;
   // initialize the file name storage and present flag
    for( ilin = 0; ilin < NFILES; ilin++ )
      {
      files_present[ilin] = 0;
      strcpy( file_store[ilin], "NONE" );
      }
    //  Make name of table file list file, open the file, read each line 
    //  and process it
    if( ( filep = getenv("OCDATAROOT") ) == NULL)
      {
      printf( "E: %s, %d: unable to get cloud table namen",
        __FILE__, __LINE__ );
      *status = 1;
      return;
      }
    strcpy( file, filep );
    strcat( file, "/cloud/" );
    strcat( file, sensorId2SensorDir( *inst_id ) );
    strcat( file, "/cloud_table_names.txt" );

    if( ( rd_stream = fopen( file, "r" ) ) == NULL )
      {
      printf( 
        "E: %s, %d: Open error for cloud table file: %s\n",
        __FILE__, __LINE__, file );
      *status = 1;
      return;
      }
    // get all inputs (name_id=name) from the stream
    while( fgets( str, 562, rd_stream ) != NULL )
      {
      if( *str != '#' )
        {
        floc = strchr( str, '=' );
        if( floc == NULL )
          {
          printf( "E: %s, %d: A data line in file: %s\n", 
            __FILE__, __LINE__, file );
          printf( "       Does not have the form NAME_ID=NAME\n" );
          printf( "  Line is: %s\n", str );
          }
        *floc = 0;
        floc++;
        id = str;
        //  find the name and store the file name
        ifound = 0;
        for( ilin = 0; ilin < NFILES; ilin++ )
          {
          if( strcmp( upcase( id ), fid_names[ilin] ) == 0 )
            {
            files_present[ilin] = 1;
            strcpy( file_store[ilin], floc );
            ifound = 1;
            break;
            }
          }
        if( ifound != 1 )
          {
          printf( "I: %s, %d: No match found for cloud table file ID: %s\n",
            __FILE__, __LINE__, id );
          }
        }
      }
    fclose( rd_stream );
   /*
    *  Be sure all the table files are assigned
    */
    for( ilin = 0; ilin < NFILES; ilin++ )
      {
      if( files_present[ilin] == 0 )
        {
        printf( "E: %s, %d: Cloud table: %s  was not found\n", 
          __FILE__, __LINE__, fid_names[ilin] );
        *status = 1;
        }
      }
    }
 /*
  *  now, here we send back the requested file name
  */
  if( *status != 0 ) return;

  // if the file starts with a '/', assume it is an entire file spec
  if( file_store[*fname_id][0] == '/' )
    {
    strcpy( fname, file_store[*fname_id] );
    }
  else
    {
    // Otherwise, assume file is in $OCDATAROOT/cloud/(common)|<sensor path>
    if( ( filep = getenv("OCDATAROOT") ) == NULL)
      {
      printf( "E: %s, %d: unable to get cloud table namen",
        __FILE__, __LINE__ );
      *status = 1;
      }
    else
      {
      strcpy( fname, filep );  
      strcat( fname, "/cloud/" );
      if( *fname_id < 9 )
        strcat( fname, "common/" );
      else
        {
        strcat( fname, sensorId2SensorDir( *inst_id ) );
        strcat( fname, "/" );
        }
      strcat( fname, file_store[*fname_id] );
      }
    }
  // see if we can pad the string with blanks to end
  ilin = strlen( fname );
  for( i = ilin - 1; i < len1; i++ )
    fname[i] = ' ';
 /*
  printf( "%s, %d: temp found file is: %s\n", __FILE__, __LINE__, fname );
  printf( "from inst: %d, sens: %d, file_store: %s\n", *inst_id, *fname_id, file_store[*fname_id] );
  */
  return;
  }
