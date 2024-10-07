/*
 *  mng_ms.cpp - take over the management of the large multi-scattering 
 *    reflectance arrays for the cloud processing
 *  W Robinson, SAIC, 17 Dec 2021
 */
#include <string.h>
#include <netcdf.h>
//#include "mfhdf.h"
#include "dim_mgr.hpp"
#include "passthebuck.h"

struct ms_finfo_struc_def
  {
  int fid_ms;  /* fid of MS table file */
  int fid_ms_std;  /* fid of MS std file (for files over water */
  int sds_id;  /* the SDS id refl table */
  int sds_id_std;  /* the SDS id refl SD table */
  int32_t *dim_siz;  /* array of dim sizes */
  double **dim_vals;  /* arrays of each dim value set */
  };
typedef struct ms_finfo_struc_def ms_finfo_struc;

ms_finfo_struc wtrcld_landsfc_info, icecld_landsfc_info;
//  for the water surface
ms_finfo_struc wtrcld_wtrsfc_info[3], icecld_wtrsfc_info[3];

dim_mgr land_mgr( 3, 83 );  /* manage separately for over water, land */
dim_mgr wtr_mgr( 3, 83 );

//  storage for the size of each reflectance / sd table size
//  the water has 18 radii and the ice has 12
int32_t wtrcld_nelts, icecld_nelts;  // size of water, ice cloud arrays

/*
 * mng_ms_init - start the management of the large tables
 */
/*******************************************************************
   mng_ms_init_  initialize the management of the large cloud multi-scatering 
        reflectance arrays

   Returns:
     int32_t - nothing at the moment

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            fil_wtrcld_landsfc I    file name for water cloud
                                                land sfc arrays
      char *            fil_icecld_landsfc I    file name for ice cloud 
                                                land sfc arrays
      int               len1             I      length of 1st char *
      int               len2             I      length of 2nd char *
      The above is for 2 land surface files and their lengths
      The rest are the files for the water surface and the lengths after
      them all (a fortran convention).  below, I'll just list the file names
      and their lengths same types as above
      file name                   length        Description
      -----------------------     ----------    -----------------------------
      fil_wtrcld_landsfc          len1          water cloud, land sfc
      fil_icecld_landsfc          len2          ice cloud, land sfc
      fil_wtrcld_wtrsfc1          len_ww1       water cloud, water sfc wind1
      fil_wtrcld_wtrsfc1_sd       len_ww1sd     water cloud, water sfc wind1 sd
      fil_icecld_wtrsfc1          len_iw1       ice cloud, water sfc wind1
      fil_icecld_wtrsfc1_sd       len_iw1sd     ice cloud, water sfc wind1 sd
      fil_wtrcld_wtrsfc2          len_ww2       water cloud, water sfc wind2
      fil_wtrcld_wtrsfc2_sd       len_ww2sd     water cloud, water sfc wind2 sd
      fil_icecld_wtrsfc2          len_iw2       ice cloud, water sfc wind2
      fil_icecld_wtrsfc2_sd       len_iw2sd     ice cloud, water sfc wind sd
      fil_wtrcld_wtrsfc3          len_ww3       water cloud, water sfc wind3
      fil_wtrcld_wtrsfc3_sd       len_ww3sd     water cloud, water sfc wind3 sd
      fil_icecld_wtrsfc3          len_iw3       ice cloud, water sfc wind3
      fil_icecld_wtrsfc3_sd       len_iw3sd     ice cloud, water sfc wind3 sd

********************************************************************/
extern"C" int32_t mng_ms_init_( char *fil_wtrcld_landsfc, 
  char *fil_icecld_landsfc, 
  char *fil_wtrcld_wtrsfc1, char *fil_wtrcld_wtrsfc1_sd,
  char *fil_icecld_wtrsfc1, char *fil_icecld_wtrsfc1_sd,
  char *fil_wtrcld_wtrsfc2, char *fil_wtrcld_wtrsfc2_sd,
  char *fil_icecld_wtrsfc2, char *fil_icecld_wtrsfc2_sd,
  char *fil_wtrcld_wtrsfc3, char *fil_wtrcld_wtrsfc3_sd,
  char *fil_icecld_wtrsfc3, char *fil_icecld_wtrsfc3_sd,
  int len1, int len2, int len_ww1, int len_ww1sd, int len_iw1, int len_iw1sd,
  int len_ww2, int len_ww2sd, int len_iw2, int len_iw2sd, 
  int len_ww3, int len_ww3sd, int len_iw3, int len_iw3sd )
  {
  char *fn, *fn2;
  char *ftn_str_cond( char *, int );
  int32_t status;
  int32_t set_ms_file( char *, char *,  int32_t , ms_finfo_struc * );

 /* set up the water and ice cloud over land file table info */
  fn = ftn_str_cond( fil_wtrcld_landsfc, len1 );
  set_ms_file( fn, NULL, 0, &wtrcld_landsfc_info );
  free( fn );

  fn = ftn_str_cond( fil_icecld_landsfc, len2 );
  set_ms_file( fn, NULL, 0, &icecld_landsfc_info );
  free( fn );

 /* make sure the first 5 dims, less radius, are equal 
    Also, set up the # elements for ice, water in the last 3 dims */
  status = 0;
  wtrcld_nelts = wtrcld_landsfc_info.dim_siz[3] + 1;  // for tau
  icecld_nelts = icecld_landsfc_info.dim_siz[3] + 1;
  for( int idim = 0; idim < 6; idim++ )
    {
    if( ( idim != 5 ) && ( wtrcld_landsfc_info.dim_siz[idim] != 
      icecld_landsfc_info.dim_siz[idim] ) )
      {
      status = 1;
      printf( "%s, %d, E: Water and ice cloud over land dimension mismatch\n",
        __FILE__, __LINE__ );
      printf( "   file_water: %s\n", fil_wtrcld_landsfc );
      printf( "   file_ice  : %s\n", fil_icecld_landsfc );
      }
    //  set up the size of the water and ice # elements,
    //  that will go in the data blob
    if( idim > 3 ) // for band and radius
      {
      wtrcld_nelts *= wtrcld_landsfc_info.dim_siz[idim];
      icecld_nelts *= icecld_landsfc_info.dim_siz[idim];
      }
    }
  if( status != 0 ) exit(27);

/*
 *  setup the water and ice cloud over water surface tables
 *  loop through all 3 water surface types (3, 7, 15 m/s), get the SDS IDs
 *  and check consostency of dims with over-land info
 */
  status = 0;
  for( int32_t isfc = 0; isfc < 3; isfc++ )
    {
    if( isfc == 0 ) {
      fn = ftn_str_cond( fil_wtrcld_wtrsfc1, len_ww1 );
      fn2 = ftn_str_cond( fil_wtrcld_wtrsfc1_sd, len_ww1sd );
    } else if( isfc == 1 ) {
      fn = ftn_str_cond( fil_wtrcld_wtrsfc2, len_ww2 );
      fn2 = ftn_str_cond( fil_wtrcld_wtrsfc2_sd, len_ww2sd );
    } else if( isfc == 2 ) {
      fn = ftn_str_cond( fil_wtrcld_wtrsfc3, len_ww3 );
      fn2 = ftn_str_cond( fil_wtrcld_wtrsfc3_sd, len_ww3sd );
    }
    set_ms_file( fn, fn2, 1, &wtrcld_wtrsfc_info[isfc] );
    free( fn );
    free( fn2 );

    if( isfc == 0 ) {
      fn = ftn_str_cond( fil_icecld_wtrsfc1, len_iw1 );
      fn2 = ftn_str_cond( fil_icecld_wtrsfc1_sd, len_iw1sd );
    } else if( isfc == 1 ) { 
      fn = ftn_str_cond( fil_icecld_wtrsfc2, len_iw2 );
      fn2 = ftn_str_cond( fil_icecld_wtrsfc2_sd, len_iw2sd );
    } else if( isfc == 2 ) {
      fn = ftn_str_cond( fil_icecld_wtrsfc3, len_iw3 );
      fn2 = ftn_str_cond( fil_icecld_wtrsfc3_sd, len_iw3sd );
    }
    set_ms_file( fn, fn2, 1, &icecld_wtrsfc_info[isfc] );
    free( fn );
    free( fn2 );

    for( int idim = 0; idim < 6; idim++ )
      {
      if( ( ( idim == 3 ) && ( wtrcld_wtrsfc_info[isfc].dim_siz[idim] != 
         wtrcld_landsfc_info.dim_siz[idim] + 1 ) ) 
        || ( ( idim != 3 ) && (wtrcld_wtrsfc_info[isfc].dim_siz[idim] !=
        wtrcld_landsfc_info.dim_siz[idim] ) ) )
        {
        status = 1;
        printf( 
          "%s, %d, E: water cloud over-water MS file dimension mismatch\n", 
          __FILE__, __LINE__ );
        printf( "     dim # %d, water sfc # %d\n", idim, isfc );
        }
      if( status != 0 ) exit(27);

      if( ( ( idim == 3 ) && ( icecld_wtrsfc_info[isfc].dim_siz[idim] !=
         icecld_landsfc_info.dim_siz[idim] + 1 ) ) 
        || ( ( idim != 3 ) && (icecld_wtrsfc_info[isfc].dim_siz[idim] !=
        icecld_landsfc_info.dim_siz[idim] ) ) )
        {
        status = 1;
        printf( 
          "%s, %d, E: ice cloud over-water MS file dimension mismatch\n", 
          __FILE__, __LINE__ );
        printf( "     dim # %d, water sfc # %d\n", idim, isfc );
        }
     if( status != 0 ) exit(27);
      }
    }
 /* set up the dimension manager for points over land */
  //dim_mgr land_mgr( 3, 80 );
  for( int idim = 0; idim < 3; idim++ )
    {
    land_mgr.init_dim( idim, wtrcld_landsfc_info.dim_siz[idim], 
      wtrcld_landsfc_info.dim_vals[idim] );
    wtr_mgr.init_dim( idim, wtrcld_wtrsfc_info[0].dim_siz[idim],
      wtrcld_wtrsfc_info[0].dim_vals[idim] );
    }

  return 0;
  }

/********************************************************************
   ftn_str_cond - Just condition a string to be null terminated

   Returns:
     char *  the good null terminated astring
  Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            str              I      initial fortran string
      int               len              I      string length
********************************************************************/
extern"C" char *ftn_str_cond( char *str, int len )
  {
  char *fn;

  fn = (char *) malloc( ( len + 1 ) * sizeof( char ) );
  strncpy( fn, str, len );
  fn[ len ] = (char) 0;
  return fn;
  }

extern"C" int32_t set_ms_file( char *fil_ms, char *fil_ms_std, int32_t fil_typ, 
  ms_finfo_struc *finfo )
/*******************************************************************

   set_ms_file - set up all the info for the multi-scattering file table

   Returns type: int - size of dimension or -1 if problem

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            fil_ms           I      name of multi-scatering file
      char *            fil_ms_std       I      name of multi-scatering std 
                                                file
      int32_t           fil_typ          I      0 for a land surface file 
                            (has the MS and MSstd tables, 1 for the 
                            water surface files (each file has the MS or MSstd)
      ms_finfo_struc *  finfo            O      all information of the file(s)

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       20 Dec 2017     original development

*******************************************************************/
  {
  int rank, dim_id_lst[7], nattr, ret1, vtyp;
  size_t dim_lens[7];
  int sds_id;
  size_t start[7], end[7];
  int32_t trg_rank = 6;
  float *fstore;

  char *dim_sds_names[] = { "ReflectanceSensorZenith", "ReflectanceSolarZenith",
    "ReflectanceRelativeAzimuth", "OpticalThickness", "Wavelengths", 
    "ParticleRadius" };

 /* 
  *  open the file and get the refl array SDS ID
  */
  // finfo->fid_ms = SDstart( fil_ms, DFACC_READ );
  DPTB( nc_open( fil_ms, NC_NOWRITE, &(finfo->fid_ms) ) );
  //finfo->sds_id = SDselect( finfo->fid_ms, SDnametoindex( finfo->fid_ms,       
  //  "MultiScatBDReflectance" ) );
  DPTB( nc_inq_varid( finfo->fid_ms, "MultiScatBDReflectance", 
    &(finfo->sds_id) ) );

  //SDgetinfo( finfo->sds_id, bl_nam, &rank, dim_lens, &vtyp, &nattr );
  DPTB( nc_inq_var( finfo->fid_ms, finfo->sds_id, NULL, &vtyp, &rank, 
    dim_id_lst, &nattr ) );
  
  if( rank != trg_rank )
    {
    printf( "E: %s, %d: MS refl in cloud table file has wrong size\n",
      __FILE__, __LINE__ );
    printf( "      file: %s", fil_ms );
    printf( "      SDS: MultiScatBDReflectance\n" );
    exit(27);
    }
 /*
  *  get dim lengths
  */
  for( int i = 0; i < rank; i++ )
    DPTB( nc_inq_dimlen( finfo->fid_ms, dim_id_lst[i], dim_lens + i ) );
 /*
  *  for the wind-dependent tables, the SD is in separate file
  */
  if( fil_typ == 0 )
    finfo->fid_ms_std = finfo->fid_ms;
  else
    DPTB( nc_open( fil_ms_std, NC_NOWRITE, &(finfo->fid_ms_std) ) );
    //finfo->fid_ms_std = SDstart( fil_ms_std, DFACC_READ ); 
 /*
  *  get the Std Dev SDS ID 
  */
  //finfo->sds_id_std = SDselect( finfo->fid_ms_std, 
  //   SDnametoindex( finfo->fid_ms_std, "StdDevMultiScatBDReflectance" ) );
  DPTB( nc_inq_varid( finfo->fid_ms_std, "StdDevMultiScatBDReflectance",
    &(finfo->sds_id_std) ) );

  //SDgetinfo( finfo->sds_id_std, bl_nam, &rank, dim_lens2, &vtyp, &nattr );
  DPTB( nc_inq_var( finfo->fid_ms_std, finfo->sds_id_std, NULL, &vtyp, &rank,
    dim_id_lst, &nattr ) );

  if( rank != trg_rank ) 
    {
    printf( "E: %s, %d: MS refl std in cloud table file has wrong size\n",
      __FILE__, __LINE__ );
    printf( "      file: %s", fil_ms );
    printf( "      SDS: StdMultiScatBDReflectance\n" );
    exit(27);
    }
 /*  Note that the data in the file is float and we are putting it into 
     double, so some translation will be needed */
  finfo->dim_siz = (int32_t *) malloc( trg_rank * sizeof( int32_t ) );
  finfo->dim_vals = ( double **) malloc( trg_rank * sizeof( double * ) );

 /*
  *  get information for all 6 dimensions
  */
  for( int32_t isds = 0; isds < 6; isds++ )
    {
    if( nc_inq_varid( finfo->fid_ms, dim_sds_names[isds], &sds_id )
      != NC_NOERR )
      {
      printf( "E: %s, %d: Can't read table SDS %s\n", __FILE__, __LINE__,
        dim_sds_names[isds] );
      exit(26);
      }
    //SDgetinfo( sds_id, bl_nam, &rank, dim_lens2, &vtyp, &nattr );
    DPTB( nc_inq_var( finfo->fid_ms, sds_id, NULL, &vtyp, &rank,
      dim_id_lst, &nattr ) );
    DPTB( nc_inq_dimlen( finfo->fid_ms, dim_id_lst[0], dim_lens ) );

    finfo->dim_siz[isds] = dim_lens[0];
    finfo->dim_vals[isds] = (double *) malloc( dim_lens[0] * sizeof(double) );
    fstore = (float *) malloc( dim_lens[0] * sizeof(float) );
    start[0] = 0;
    end[0] = dim_lens[0];
    //ret1 = SDreaddata( sds_id, start, NULL, end, (float *)fstore );
    ret1 = nc_get_vara_float( finfo->fid_ms, sds_id, start, end, 
      (float *)fstore );
    if( ret1 != NC_NOERR )
      {
      printf( "E: %s, %d: Can't read table SDS %s\n", __FILE__, __LINE__,
        dim_sds_names[isds] );
      exit(26);
      }
    for( size_t ipx = 0; ipx < end[0]; ipx++ )
      finfo->dim_vals[isds][ipx] = *(fstore + ipx );
    free( fstore );
    //ret2 = SDendaccess( sds_id );
    }
  return 0;
  }

/* the next is the retrieval of data side.  It will have a part to get the 
   structure defining the set of data grid points and the weights to apply
   And, a part to make the particular interpolated array requested
*/

extern"C" int mng_ms_get_lib_( float *sol_ang, float *sen_ang, 
  float *relaz_ang, int *sfc_typ, int *meas_typ, int *scan, float *wtr_int, 
  float *ice_int, int *stat )
/*******************************************************************

   mng_ms_get_lib - get the portion of the MS reflectance library 
       interpolated to the solar, sensor and relative azimuth angles

   Returns type: int - 0 good, otherwise an error
      Note that in th ecase where the point is outside the table, a 2 status 
      is returned ans will be handled by making the final table all 0, 
      hopefully that will just make for fill values in cloud quantities
      dealing with the calling code is horrendous.

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float *           sol_ang          I      solar angle
      float *           sen_ang          I      sensor angle
      float *           relaz_ang        I      relative azimuth angle
      int *             sfc_typ          I      underlying surface type:
                                                0 - land, 1 - water, 3 m/s 
                                                wind, 2 - water, 7 m/s, 
                                                3 - water, 15 m/s
      int *             meas_typ         I      measurement type to return:
                                                0 - reflectance, 1 - Std Dev
      int *             scan             I      current scan line, use as 
                                                access id
      float *           wtr_int          O      interpolated water array
      float *           ice_int          O      interpolated ice array
      int *             stat             O      status - 0 is good, 5 is out
                                                of table bounds

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       23 Dec 2021     original development

NOTE - organization of data blobs
position
  1 - reflectance wtrcld   size wtrcld_nelts
  2 - SD wtrcld            size wtrcld_nelts
  3 - reflectance icecld   size icecld_nelts
  2 - SD icecld            size icecld_nelts
  This repeats 2 more for the water surface arrays

*******************************************************************/
  {
  /*  follow the test program for dim_mgr */
  /*   */
  int32_t access_id;
  //static pt_info_struc *pt_info;
  pt_info_struc *pt_info;
  double pt[3];
  int32_t status, ndat;
  int32_t nsfc;  // nsfc will be for the # of surface types (1 land, 3 water)
  int32_t narr = 4;  // # of arrays to get out from the wtrcld, wtrcld_sd,
                     // icecld, icecld_sd SDSes
  int32_t ndim = 3, off_wtr, off_ice;
  int32_t ntau, ntau1, nwav, nrad_wtr, nrad_ice, nrad;
  static int32_t prune_called = 0;

  //  ntau1 represents ntau + 1 the size in tau of the final storage
  //  the over-water data has ntau+1 while the over-land has ntau
  //  its just the (odd) way the cloud code, files were set up

  // WDR temp to report the points managed when the access_id changes
  /*  leave out normally  */
//printf( "%s, %d, SCAN: %d\n", __FILE__, __LINE__, *scan );
  //if( ( *scan % 20 ) == 19 )
/*
  if( *scan == 300 )
    {
    double way = 0.;
    printf( "%s, %d, LAND MANAGER DUMP\n", __FILE__, __LINE__ );
    land_mgr.dump_mgr( &way );
    //printf( "%s, %d, WATER MANAGER DUMP\n", __FILE__, __LINE__ );
    // no need wtr_mgr.dump_mgr( &way );
    exit( 0 );
    }
  */
  // WDR end of the probe of points

  //  WDR add the prune, every 200 lines prune all 100 less than current line
  if( ( *scan % 200 ) == 199 )
  //if( ( *scan % 300 ) == 299 )
    {
    if( prune_called == 0 )
      {
      //double way = 0.;
      //printf( "LAND MANAGER DUMP PRE, scan %d\n", *scan );
      //land_mgr.dump_mgr( &way );
      //printf( "WATER MANAGER DUMP PRE, scan %d\n", *scan );
      //wtr_mgr.dump_mgr( &way );
      land_mgr.prune( *scan - 100 );
      wtr_mgr.prune( *scan - 100 );
      prune_called = 1;
      //printf( "LAND MANAGER DUMP POST, scan %d\n", *scan );
      //land_mgr.dump_mgr( &way );
      }
    }
  else
    prune_called = 0;

  access_id = *scan;
  ndat = pow( 2, ndim );
  //  Note that the organization is senz, solz, relaz
  pt[0] = *sen_ang;
  pt[1] = *sol_ang;
  pt[2] = *relaz_ang;

  //  do over-land and over-water separately - no need to read land table 
  //  over water or visa-versa
  if( *sfc_typ == 0 )  // a request for over-land table data
    {
    pt_info = land_mgr.mng_pt( pt, access_id, &status );
    nsfc = 1;
    }
  else
    {
    pt_info = wtr_mgr.mng_pt( pt, access_id, &status );
    nsfc = 3;
    }

  if( status != 0 )
    {
    //  trying to deal with exceptions in the calling code is a mess
    //  so the error 2 will end up returning an unfilled table
    if( status == 2 )  // out of table bounds
      {
      *stat = 5;
      return 5;
      }
    printf( "%s, %d: mng_ms_get_lib error found of %d\n", __FILE__,
      __LINE__, status );
    exit( 27 );
    }
  else
    {
    ntau = wtrcld_landsfc_info.dim_siz[3];
    ntau1 = ntau + 1;
    nwav = wtrcld_landsfc_info.dim_siz[4];
    nrad_wtr = wtrcld_landsfc_info.dim_siz[5];
    nrad_ice = icecld_landsfc_info.dim_siz[5];

    if( pt_info->interval_needs_data == 1 )
      {
     /* we need to read data for some of the grid points in the interval */
      int32_t n_new_blobs = 0;
      int32_t offset[3], arr_off, tau_off;
      int tsds_id, tfil_id;
     /* set-up for the use of MS arrays symmetric in senz / solz */
      double *senz_set = &( wtrcld_landsfc_info.dim_vals[0][0] );
      double *solz_set = &( wtrcld_landsfc_info.dim_vals[1][0] );
      int32_t n_senz = wtrcld_landsfc_info.dim_siz[0];
      int32_t n_solz = wtrcld_landsfc_info.dim_siz[1];
      int32_t i;
      double grid_senz, grid_solz;
     /*  end */
      size_t istart[6] = { 0, 0, 0, 0, 0, 0 };
      size_t icount[6] = { 1, 1, 1, (size_t)ntau, (size_t)nwav, 
        (size_t)nrad_wtr };
      if( *sfc_typ != 0 )  icount[3] = ntau1;

      for( int32_t idat = 0; idat < ndat; idat++ ) // loop interval corners
        {
        if( pt_info->pt_status[idat] == -1 )
          {
          n_new_blobs++;
         /*  do the read, reorg, store in struc and deposit in 
             a data blob pointer, use calloc to avoid needing to zero 
             the array for a particular tau index that changes */
          float *dat_blob = (float *)calloc( ( wtrcld_nelts * 2 +
            icecld_nelts * 2 ) * nsfc, sizeof(float) );
          linear_to_offset( ndim, idat, offset );
          for( int32_t idim = 0; idim < 3; idim++ )
            istart[idim] = pt_info->pt_base_loc[idim] + offset[idim];

          // to use only part of MS tables where solz <= senz
          grid_senz = senz_set[ istart[0] ];
          grid_solz = solz_set[ istart[1] ];
          if( grid_senz < grid_solz )
          // for removing the flip if( 0 )
            {
            for( i = 0; i < n_solz; i++ )
              {
              if( solz_set[i] == grid_senz )
                {
                istart[1] = i;
                break;
                }
              }
            for( i = 0; i < n_senz; i++ )
              {
              if( senz_set[i] == grid_solz )
                {
                istart[0] = i;
                break;
                }
              }
            }
         else {
           //printf( "not a flip\n" );
         }
//  just re-set the istart to see if intervening code has any effect
          //for( int32_t idim = 0; idim < 3; idim++ )
          //  istart[idim] = pt_info->pt_base_loc[idim] + offset[idim];

          // loop thru the nsfc surface types (3 for water sfc, 1 for land
          for( int32_t isfc = 0; isfc < nsfc; isfc++ )
            {
            // loop thru the arrays to fill the data blob, do - 4 of them
            // per surface: wtrcld ref, sd, icecld ref, sd (etc for more sfc)
            for( int32_t iarr = 0; iarr < narr; iarr++ )
              {
              nrad = nrad_wtr;
              tau_off = 0;  /* for land sfc: is 0 to transfer data to tau 
                                 dim of blob starting at 0, 1 to start at 
                                 1 for the SD (to get expected arrays) 
                               for wtr sfc: is always 0 */
              if( *sfc_typ == 0 )  // the land sfc tables are 1 way
                {
                switch (iarr) 
                  {
                  case 0:  // water cloud refl
                    tsds_id = wtrcld_landsfc_info.sds_id;
                    tfil_id = wtrcld_landsfc_info.fid_ms;
                    arr_off = 0;
                    break;
                  case 1:  //  water cloud SD
                    tsds_id = wtrcld_landsfc_info.sds_id_std;
                    tfil_id = wtrcld_landsfc_info.fid_ms_std;
                    arr_off = wtrcld_nelts;
                    tau_off = 1;
                    break;
                  case 2:  //  ice cloud refl
                    tsds_id = icecld_landsfc_info.sds_id;
                    tfil_id = icecld_landsfc_info.fid_ms;
                    arr_off = 2 * wtrcld_nelts;
                    nrad = nrad_ice;
                    break;
                  case 3:  //  ice cloud SD
                    tsds_id = icecld_landsfc_info.sds_id_std;
                    tfil_id = icecld_landsfc_info.fid_ms_std;
                    arr_off = 2 * wtrcld_nelts + icecld_nelts;
                    nrad = nrad_ice;
                    tau_off = 1;
                    break;
                  }
                }
              else   // the water surface tables are another way
                {
                tau_off = 0;
                switch (iarr)
                  {
                  case 0:  // water cloud refl
                    tsds_id = wtrcld_wtrsfc_info[isfc].sds_id;
                    tfil_id = wtrcld_wtrsfc_info[isfc].fid_ms;
                    arr_off = isfc * 2 * ( wtrcld_nelts + icecld_nelts );
                    break;
                  case 1:  //  water cloud SD
                    tsds_id = wtrcld_wtrsfc_info[isfc].sds_id_std;
                    tfil_id = wtrcld_wtrsfc_info[isfc].fid_ms_std;
                    arr_off = wtrcld_nelts + 
                      isfc * 2 * ( wtrcld_nelts + icecld_nelts );
                    break;
                  case 2:  //  ice cloud refl
                    tsds_id = icecld_wtrsfc_info[isfc].sds_id;
                    tfil_id = icecld_wtrsfc_info[isfc].fid_ms;
                    arr_off = 2 * wtrcld_nelts +
                      isfc * 2 * ( wtrcld_nelts + icecld_nelts );
                    nrad = nrad_ice;
                    break;
                  case 3:  //  ice cloud SD
                    tsds_id = icecld_wtrsfc_info[isfc].sds_id_std;
                    tfil_id = icecld_wtrsfc_info[isfc].fid_ms_std;
                    arr_off = 2 * wtrcld_nelts + icecld_nelts +
                      isfc * 2 * ( wtrcld_nelts + icecld_nelts );
                    nrad = nrad_ice;
                    break;
                  }
                }
              icount[5] = nrad;
             /* read from the file */
              int32_t ntau_mod = ( *sfc_typ == 0 ) ? ntau : ntau1;
              float *tmp_dat = (float *)malloc( ntau_mod * nwav * nrad *
                sizeof(float) );
              //if( status = ( SDreaddata( tsds_id, istart, NULL, icount, 
              //  (void *)tmp_dat ) ) != SUCCEED )
              if( status = nc_get_vara_float( tfil_id, tsds_id, istart, icount,
                tmp_dat ) != NC_NOERR )
                {
                printf( 
                  "E: %s, %d, Read of cloud table data failed with status: %d\n",
                  __FILE__, __LINE__, status );
                exit(27);
                }
             /* transfer water ref to data blob */
              for( int32_t irad = 0; irad < nrad; irad++ )
                {
                for( int32_t iwav = 0; iwav < nwav; iwav++ )
                  {
                  // for water sfc, they have 1 more tau, so the setting 0 
                  //  not needed
                  if( *sfc_typ == 0 )  // = land, do this stuff
                    {
                    for( int32_t itau = 0; itau < ntau; itau++ )
                      {
                      *( dat_blob + arr_off + ( itau + tau_off ) + ntau1 * 
                        ( iwav + nwav * irad ) ) =
                        *( tmp_dat + irad + nrad * ( iwav + nwav * itau ) );
                      }
                    }
                  else  // - water surface, it has a 0 tau in the sds
                    {
                    for( int32_t itau = 0; itau < ntau + 1; itau++ )
                      {
                      *( dat_blob + arr_off + itau + ntau1 *
                        ( iwav + nwav * irad ) ) =
                        *( tmp_dat + irad + nrad * (iwav + nwav * itau ) );
                      }
                    }
                  }
                }
              free(tmp_dat);
              }
            }
         /* place data blob pointer into the pt_info  */
          pt_info->dat_ptrs[idat] = (void *)dat_blob;
         /* set status of the point to say data is there */
          pt_info->pt_status[idat] = 1;
          }
        }

     /*  update knowlege of point status and data location */
//printf( "%s %d, before update, sfc type: %d\n", __FILE__, __LINE__, *sfc_typ );
      if( *sfc_typ == 0 )    //  over land
        {
        land_mgr.update_new_data( );
        land_mgr.add_pts( n_new_blobs );
        }
      else
        {
        wtr_mgr.update_new_data( );
        wtr_mgr.add_pts( n_new_blobs );
        }
      }

   /*  at this point, go through the interval grid points and weight the 
       proper data and sum to the interpolated array */
    //  NOTE - this assumes the interpolated array is ZERO on input
    int32_t sfc_typ_mult [] = { 0, 0, 1, 2 };
    int32_t sfc_typ_off = 2 * ( wtrcld_nelts + icecld_nelts );

    if( *meas_typ == 0 )  //  refl
      {
      off_wtr = 0 + sfc_typ_mult[*sfc_typ] * sfc_typ_off;
      off_ice = 2 * wtrcld_nelts + sfc_typ_mult[*sfc_typ] * sfc_typ_off;
      }
    else  //  SD
      {
      off_wtr = wtrcld_nelts + sfc_typ_mult[*sfc_typ] * sfc_typ_off;
      off_ice = 2 * wtrcld_nelts + icecld_nelts + 
        sfc_typ_mult[*sfc_typ] * sfc_typ_off;
      }
    for( int idat = 0; idat < ndat; idat++ )
      {
      // this is just for doc, may go soon
      //printf( "%s, %d, for interval at point %d, weight: %f\n", __FILE__, 
      //  __LINE__, idat, pt_info->wt_pt[idat] );
      for( int32_t itau = 0; itau < ntau1; itau++ )
        {
        for( int32_t iwav = 0; iwav < nwav; iwav++ )
          {
          for( int32_t irad = 0; irad < nrad_wtr; irad++ )
            {
            *( wtr_int + itau + ntau1 * ( iwav + nwav * irad ) ) +=
              pt_info->wt_pt[idat] *
              *( (float *)pt_info->dat_ptrs[idat] + off_wtr + itau + ntau1 * 
              ( iwav + nwav * irad ) );
            }

          for( int32_t irad = 0; irad < nrad_ice; irad++ )
            {
            *( ice_int + itau + ntau1 * ( iwav + nwav * irad ) ) += 
              pt_info->wt_pt[idat] *
              *( (float *)pt_info->dat_ptrs[idat] + off_ice + itau + ntau1 *
              ( iwav + nwav * irad ) );
            }
          }
        }
      }
    }
 /* at the end... */ 
  //printf( "I %s, %d, interval loc: %d, %d, %d\n", __FILE__, __LINE__, pt_info->pt_base_loc[0], pt_info->pt_base_loc[1], pt_info->pt_base_loc[2] );
  return 0;
  }
