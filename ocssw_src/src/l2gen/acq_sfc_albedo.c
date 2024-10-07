#include "l12_proto.h"
#include <string.h>

// Sizes of the reflectance and uncertainty arrays read with acq_cth_albedo()
#define CTH_NWAVE  21
#define CTH_NMON   12
#define CTH_NLAT   1440
#define CTH_NLON   2880

//  Note that in May, 2023, we added reading of the albedo needed
//  for doing the cloud top height calculation, init_cld_dat set up 
//  2 more arrays and routine acq_cth_albedo() will read data for the 
//  line of pixels

int init_cld_dat( l1str *l1rec )
/*******************************************************************

   init_cld_dat

   purpose: initialize the cloud structure portion of the l1rec

   Returns type: int32_t - status 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      l1str *           l1rec           I/O     all L1 line info

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       3 Apr 2020     Original development
      W. Robinson, SAIC 1 Jun 2023     set up the cth arrays: cth_alb_init,
                                       cth_alb_unc_init

*******************************************************************/ {
    int32_t npix;

    npix = l1rec->npix;

    if( ( l1rec->cld_dat = malloc( sizeof( cld_struc ) ) ) == NULL )
      {
      printf( "-E- %s, %d: Failure to allocate cloud albedo data structure\n",
        __FILE__, __LINE__ );
      return -1;
      }
    if( ( ( l1rec->cld_dat->sfc_albedo_659 = (float *) 
            malloc( npix * sizeof( float ) ) ) == NULL ) || 
        ( ( l1rec->cld_dat->sfc_albedo_858 = (float *)
            malloc( npix * sizeof( float ) ) ) == NULL ) ||
        ( ( l1rec->cld_dat->sfc_albedo_1240 = (float *)
            malloc( npix * sizeof( float ) ) ) == NULL ) ||
        ( ( l1rec->cld_dat->sfc_albedo_1640 = (float *)
            malloc( npix * sizeof( float ) ) ) == NULL ) ||
        ( ( l1rec->cld_dat->sfc_albedo_2130 = (float *)
            malloc( npix * sizeof( float ) ) ) == NULL ) ||
        ( ( l1rec->cld_dat->cth_alb_init = (float *)
            malloc( npix * sizeof( float ) ) ) == NULL ) ||
        ( ( l1rec->cld_dat->cth_alb_unc_init = (float *)
            malloc( npix * sizeof( float ) ) ) == NULL ) )
      {
      printf( "-E- %s, %d: Failure to allocate cloud albedo data structures\n",
        __FILE__, __LINE__ );
      return -1;
      }
    return 0;
    }

int read_albedo( int32_t d_np, int32_t d_nl, int32_t st_lin, int32_t ix_clim, 
  int ncid, int *d_id, unsigned char ***alb_dat )
/*******************************************************************

   read_albedo

   purpose: get the sfc albedo read from the albedo file

   Returns type: int32_t - status 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32_t           d_np             I      albedo data # pixels
      int32_t           d_nl             I      albedo data # lines
      int32_t           st_lin           I      start line to readfrom albedo
                                                file
      int32_t           ix_clim          I      time to read -1 for no time 
                                                (not climatology) else time
                                                index
      int               ncid             I      albedo file id
      int *             d_id             I      ids of albedo band datasets
      unsigned char *** alb_dat          O      final data buffer

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       30 Mar 2020     Original development
 *******************************************************************/ {
  int32_t iwav, i_np, i_nl;
  size_t start[] = {0,0,0}, count[] = {0,0,0};
  int32_t ilat;
  unsigned char *xfr_arr, **alb_lcl;

  alb_lcl = *alb_dat;
 /*  free permanent storage if needed */
  if( alb_lcl[0] != NULL )
    {
    for( iwav = 0; iwav < 5; iwav++)
      free( alb_lcl[iwav] );
    }
 /*  allocate the transfer array and the data arrays */
  if( ( xfr_arr =  (unsigned char *) malloc( d_np * d_nl * 
     sizeof( unsigned char ) ) ) == NULL )
     {
     printf( "-E- %s, %d: Error allocating albedo read storage\n", 
       __FILE__, __LINE__ );
     return -1;
     }

  i_np = d_np + 1;
  i_nl = d_nl + 1;

 /* get the data from the file to the albedo data array */
  for( iwav = 0; iwav < 5; iwav++)
    {
    if( ( alb_lcl[iwav] = (unsigned char *)
      malloc( i_np * i_nl * sizeof( unsigned char ) ) ) == NULL )
      {
      printf( "-E- %s, %d: Error allocating albedo read storage\n", 
        __FILE__, __LINE__ );
      return -1;
      }
   /*  Read data to the transfer array */
    if( ix_clim < 0 )
      {
     /* not a climate dataset */
      start[0] = st_lin;
      start[1] = 0;
      count[0] = d_nl;
      count[1] = d_np;
      }
    else
      {
     /* climate dataset */
      start[0] = ix_clim;
      start[1] = st_lin;
      start[2] = 0;
      count[0] = 1;
      count[1] = d_nl;
      count[2] = d_np;
      }
    if( nc_get_vara_uchar( ncid, d_id[iwav], start, count, xfr_arr ) != 
      NC_NOERR )
      {
      printf( "-E- %s, %d: Error reading albedo\n",
       __FILE__, __LINE__ );
      return -1;
      }
 
   /*  place it in the final data array */
    for( ilat = 0; ilat < d_nl; ilat++ )
      {
      memcpy( ( alb_lcl[iwav] + ilat * i_np + 1 ),
        ( xfr_arr + ilat * d_np ), d_np );
     /* note that first lon will be copied from the end to get to -180 */
      *( alb_lcl[iwav] + ilat * i_np ) =
        *( xfr_arr + ilat * d_np + ( d_np - 1 ) );
      }
   /* repeat the last line - only of interest at 90 */
    memcpy( ( alb_lcl[iwav] + ( i_nl - 1 ) * i_np ), 
      ( alb_lcl[iwav] + ( i_nl - 2 ) * i_np ), i_np );
    }
  free( xfr_arr );
  return 0;
  }

int acq_sfc_albedo(l1str *l1rec)
/*******************************************************************

   acq__sfc_albedo

   purpose: retrieve the surface albedo for the scan line

   Returns type: int32_t - status 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      l1str *           l1rec           I/O     all L1 line info

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       26 Mar 2020     Original development

   Note that this differs from the anc_acq... 1: only 1 file is read for the 
   data, and 2: the albedo is so large that provision is made to only read 
   in a portion of the albedo data.  Also, at this time, all longitudes are 
   read to limit routine complexity.

 *******************************************************************/ {
  static int32_t firstcall = 0, is_clim;
  static int32_t subset_read = 0, scan_npix;
  static float bdy_siz = 5.;
  static float scale[5];
  static float lat_min, lon_min, res1, res2;
  static double lat_res, lon_res, alb_st_lat, alb_lat_rng[2] = { 900., -900. };
  static unsigned char **alb_dat, fillv[5];
  static int32_t ix_clim = -1, sub_nl;
  static size_t d_np, d_nl, d_nt, i_np;
  int32_t iwav, ipx, sub_st, good_nav, ix_lat, ix_lon;
  int32_t n_done, ip, il;
  float dat_sub[2][2], fin_val, sum;
  float lat_rng[2] = { 900., -900. };
  int16_t year, doy;
  double sec, frac_lat, frac_lon, lat, lon;
  char *alb_file;
  char tit_clim[] = "MODIS/TERRA+Aqua BRDF/Albedo Gap-Filled Snow-Free 8-Day climatology";
  char tit_daily[] = "MODIS/TERRA+Aqua BRDF/Albedo Gap-Filled Snow-Free data";
  char str_tit[FILENAME_MAX];
  static int ncid, d_id[5], dim_id, dim_id2;
  nc_type att_typ;
  char attr_buf[1024];
  char *alb_wav_names[5] = { "Albedo_Map_0.659", "Albedo_Map_0.858", 
    "Albedo_Map_1.24", "Albedo_Map_1.64", "Albedo_Map_2.13" };
  int nc1, nc2;

  if( firstcall == 0 )
    {
    firstcall = 1;
    alb_file = input->sfc_albedo;
    //open the file
    if (nc_open(alb_file, 0, &ncid) != NC_NOERR) {
      fprintf(stderr,
        "-E- %s %d: file: %s is not netcdf, not acceptable albedo file\n",
                    __FILE__, __LINE__, alb_file);
      return -1;
      }
    //read the 'title', one of 2 possibilities:
    nc_get_att_text(ncid, NC_GLOBAL, "title", str_tit);
    if( strncmp( str_tit, tit_daily, strlen(tit_daily) ) == 0 )
      {
      printf( "%s, %d: TEMP: we found the daily file: %s\n", 
        __FILE__, __LINE__, str_tit );
      is_clim = 0;
      }
    else if( strncmp( str_tit, tit_clim, strlen(tit_clim) ) == 0 )
      {
      printf( "%s, %d: TEMP: we found the clim file: %s\n", __FILE__, 
        __LINE__, tit_clim );
      is_clim = 1;
     /* get the data day -> index to correct time in file */
      unix2yds( l1rec->scantime, &year, &doy, &sec);
      ix_clim = ( doy - 1 ) / 8;  // day 1 - 8 is index 0 etc
      printf( "-I- %s, %d: Using climate doy # %d, index: %d\n", __FILE__, 
        __LINE__, doy, ix_clim );
      printf( "-I- From file: %s\n", alb_file );
      }
    else
      {
      printf( "-E- %s, %d: Incorrect title in albedo dataset: %s\n",
        __FILE__, __LINE__, alb_file);
      return -1;
      }
  
   /*  read the geospatial_lat_min etc to get data reolution, start */
    if( ( nc_get_att_float(ncid, NC_GLOBAL, "geospatial_lat_min", &lat_min ) 
        != NC_NOERR ) || 
        ( nc_get_att_float(ncid, NC_GLOBAL, "geospatial_lon_min", &lon_min )  
        != NC_NOERR ) ) 
      {
      printf( "-E- %s, %d: Unable to read geolocation information for albedo file: %s\n",
        __FILE__, __LINE__, alb_file);
      return -1;
      }

   /*  break out the geospatial_lat_resolution, geospatial_lon_resolution
       to have flexability for formats */
    nc_inq_atttype( ncid, NC_GLOBAL, "geospatial_lat_resolution", &att_typ );
    if( att_typ == NC_FLOAT ) {
      nc1 = nc_get_att_float(ncid, NC_GLOBAL, "geospatial_lat_resolution",
        &res1 );
    } else {
      nc1 = nc_get_att(ncid, NC_GLOBAL, "geospatial_lat_resolution",
        attr_buf );
      res1 = str2resolution(attr_buf);
      res1 /= 111000;  // as the above gives meters
    }

    nc_inq_atttype( ncid, NC_GLOBAL, "geospatial_lon_resolution", &att_typ );
    if( att_typ == NC_FLOAT ) {
      nc2 = nc_get_att_float(ncid, NC_GLOBAL, "geospatial_lon_resolution",
        &res2 );
    } else {
      nc1 = nc_get_att(ncid, NC_GLOBAL, "geospatial_lon_resolution",
        attr_buf );
      res2 = str2resolution(attr_buf);
      res2 /= 111321;  // as the above gives meters, NOTE assume km at equator!
    }
    if( ( nc1 != NC_NOERR ) || ( nc2 != NC_NOERR ) )
      {
      printf( "-E- %s, %d: Unable to read geolocation information for albedo file: %s\n",
        __FILE__, __LINE__, alb_file);
      return -1;
      }

    lat_res = res1;
    lon_res = res2;

    scan_npix = l1rec->npix;
    // get the dimension sizes
    if ((nc_inq_dimid(ncid, "lat", &dim_id) != NC_NOERR) ||
        (nc_inq_dimid(ncid, "lon", &dim_id2) != NC_NOERR) )
      {
      printf( "-E- %s, %d: Error retrieving the lat/lon size\n", 
        __FILE__, __LINE__ );
      return -1;
      }
    if ( ( nc_inq_dimlen(ncid, dim_id, &d_nl ) != NC_NOERR) ||
         ( nc_inq_dimlen(ncid, dim_id2, &d_np ) != NC_NOERR) )
      {
      printf( "-E- %s, %d: Error retrieving the lat/lon size\n", 
        __FILE__, __LINE__ );
      return -1;
      }
    if( is_clim == 1 )
      {
      if (nc_inq_dimid(ncid, "doy", &dim_id) != NC_NOERR)
        {
        printf( "-E- %s, %d: Error retrieving the doy size\n", 
          __FILE__, __LINE__ );
        return -1;
        } 
      if ( nc_inq_dimlen(ncid, dim_id, &d_nt ) != NC_NOERR)
        {
        printf( "-E- %s, %d: Error retrieving the doy size\n", 
          __FILE__, __LINE__ );
        return -1;
        }
      }
     i_np = d_np + 1;
    /* allocate the 5 wavelength slots in alb_dat */
     alb_dat = (unsigned char **) malloc( 5 * sizeof( unsigned char * ) );
     alb_dat[0] = NULL;
    /* if there are not too many values to read, allocate all the 
       data space here
    */ 
    subset_read = 1;
    if( (long)d_np * (long)d_nl < 4000000 ) subset_read = 0;

    /* loop through the 5 bands  */
    for( iwav = 0; iwav < 5; iwav++ )
      {

     /* get the data id for the band */
      if( nc_inq_varid(ncid, alb_wav_names[iwav], ( d_id + iwav ) ) != NC_NOERR)
        {
        printf( "-E- %s, %d: Error retrieving the band id\n", 
          __FILE__, __LINE__ );
        return -1;
        }
     /* get the fill and scale values */
      if( ( nc_get_att_uchar( ncid, d_id[iwav], "_FillValue", 
        ( fillv + iwav ) ) != NC_NOERR ) ||
          ( nc_get_att_float( ncid, d_id[iwav], "scale_factor",
        ( scale + iwav ) ) != NC_NOERR ) )
        {
        printf( "-E- %s, %d: Error retrieving the fill or scale  value\n", 
         __FILE__, __LINE__ );
        return -1;
        }
      }
    if( subset_read == 0 )
      {
      if( read_albedo( d_np, d_nl, 0, ix_clim, ncid, d_id, &alb_dat ) != 0 )
        return -1;
      }
      alb_st_lat = -90.;
    } 
 /* initialization is complete and all data is read in if possible */

 /* get the cld_dat structure set if needed */
  if( l1rec->cld_dat == NULL )
    {
    //printf( "\n\n\n-T- %s, %d: setting up the cld_dat in l1rec again\n\n\n",
    //__FILE__, __LINE__ );
    // allocate the l1rec albedo struct
    if( init_cld_dat( l1rec ) != 0 ) return -1;
    }
 /* in the case where only a portion of the data is read in based on the 
    latitude range, compute that range and read in that portion if needed */

  if( subset_read == 1 ) 
    {
   /* find the lat range of the scan */
    for( ipx = 0; ipx < scan_npix; ipx++ )
      {
      lat = l1rec->lat[ipx];
      if( ( lat >= -90. ) && ( lat <= 90. ) )
        {
        if( lat < lat_rng[0] ) lat_rng[0] = lat;
        if( lat > lat_rng[1] ) lat_rng[1] = lat;
        }
      }
   /* see if we need to read a new range */
    if( ( alb_lat_rng[0] > lat_rng[0] ) || ( alb_lat_rng[1] < lat_rng[1] ) )
      {
     /* read in with added boundary */
      alb_lat_rng[0] = lat_rng[0] - bdy_siz;
      if( alb_lat_rng[0] < -90. ) 
        alb_lat_rng[0] = -90.;
      else{
        // make sure the range falls on the grid points
        sub_st = ( alb_lat_rng[0] + 90. ) / lat_res;
        alb_lat_rng[0] = ( sub_st * lat_res ) - 90.;
      }
      alb_lat_rng[1] = lat_rng[1] + bdy_siz;
      if( alb_lat_rng[1] > 90. ) 
        alb_lat_rng[1] = 90.;
      else{
        sub_st = ( alb_lat_rng[1] + 90. ) / lat_res;
        alb_lat_rng[1] = ( sub_st * lat_res ) - 90.;
      }
      alb_st_lat = alb_lat_rng[0];

     // printf( "\n\n\n-T- %s, %d: Reading the albedo data for new range:\n",
     //   __FILE__, __LINE__ );
      printf( "Latitude: %f - %f\n\n\n",
        alb_lat_rng[0], alb_lat_rng[1] );
      sub_nl = ( alb_lat_rng[1] - alb_lat_rng[0] + lat_res/2. ) / lat_res;
      sub_st = ( alb_st_lat + 90. ) / lat_res;
      if( read_albedo( d_np, sub_nl, sub_st, ix_clim, ncid, d_id, 
        &alb_dat ) != 0 )
        return -1;
      }
    }
 /* the grid is ready, for each scan point, get the interpolated albedo */
  for( ipx = 0; ipx < scan_npix; ipx++ )
    {
    fin_val = BAD_FLT;
    lat = l1rec->lat[ipx];
    lon = l1rec->lon[ipx];
    good_nav = 0;
   /* find grid box and position inside the box */
    if( ( lat >= -90. ) && ( lat <= 90. ) 
      && ( lon >= -180. ) && (lon <= 180. ) )
      {
      good_nav = 1;
      frac_lat = ( lat - alb_st_lat ) / lat_res;
      frac_lon = ( lon + 180 ) / lon_res;
      ix_lat = (int32_t) frac_lat;
      ix_lon = (int32_t) frac_lon;
      frac_lat = frac_lat - (float) ix_lat; /* get remainder */
      frac_lon = frac_lon - (float) ix_lon;
      }
   /* for the wavelengths, get the interpolated value */
    for( iwav = 0; iwav < 5; iwav++ )
      {
      if( good_nav == 1 )
        {
       /* get bounding box values */
        n_done = 0;
        sum = 0;
        for( ip = 0; ip < 2; ip++ )
          {
          for( il = 0; il < 2; il++ )
            {
            if( ix_lat + il >= sub_nl )
              dat_sub[ip][il] = fillv[iwav];
            else
              dat_sub[ip][il] = *( alb_dat[iwav] + ( ix_lon + ip ) +
                i_np * ( ix_lat + il ) );
            if( dat_sub[ip][il] != fillv[iwav] )
              {
              sum += dat_sub[ip][il] * scale[iwav];
              n_done++;
              }
            }
          }
        if( n_done == 0 )
          {
          fin_val = BAD_FLT;
          }
        else
          {
          sum /= n_done;
          for( ip = 0; ip < 2; ip++ )
            for( il = 0; il < 2; il++ )
              {
              if( dat_sub[ip][il] == fillv[iwav] )
                dat_sub[ip][il] = sum;
              else
                dat_sub[ip][il] *= scale[iwav];
              }
         /* get the weighted average */
          fin_val = ( 1 - frac_lon ) * ( ( dat_sub[0][0] * ( 1 - frac_lat ) ) +
                                     ( dat_sub[0][1] * frac_lat ) ) +
                          frac_lon *   ( ( dat_sub[1][0] * ( 1 - frac_lat ) ) +
                                     ( dat_sub[1][1] * frac_lat  ) );
          }
        }
     /*  place the final value */
      switch( iwav )
        {
        case 0: l1rec->cld_dat->sfc_albedo_659[ipx] = fin_val;
          break;
        case 1: l1rec->cld_dat->sfc_albedo_858[ipx] = fin_val;
          break;
        case 2: l1rec->cld_dat->sfc_albedo_1240[ipx] = fin_val;
          break;
        case 3: l1rec->cld_dat->sfc_albedo_1640[ipx] = fin_val;
          break;
        case 4: l1rec->cld_dat->sfc_albedo_2130[ipx] = fin_val;
        }
      }
    }
 /* and finish the line */
  return 0;
  }

int acq_cth_albedo(l1str *l1rec)
/*******************************************************************

   acq_cth_albedo

   purpose: retrieve the surface albedo for the scan line for use to make
     cloud top height - This does what A. Sayer did in his IDL code: 
     get data from the month of the granule and for the nearest grid point 
     to the pixel lat, lon

   Returns type: int32_t - status 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      l1str *           l1rec           I/O     all L1 line info

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       31 May 2023     Original development

 *******************************************************************/ {
  static size_t start[] = { 0,0,0,0}, count[] = { 0,0,0,0};
  static float *alb_arr, *alb_unc_arr;
  static int32_t firstcall = 0;
  float *arr_ptr, talb;
  size_t nwave, nmon, nlat, nlon, ix_clim;
  int16_t year, mon, day, hh, mm;
  int32_t ipix, ix, iy;
  double sec;
  char *alb_file, str_tit[5000];
  char tit_clim[] = "Surface Lambertian-equivalent reflectivity (LER) observed by TROPOMI";
  char *arr_nms[] = { "minimum_LER_clear", "uncertainty_clear" };
  int ncid, d_id, dim_id[4];
  int status;

  if( firstcall == 0 )
    {
    firstcall = 1;
    alb_file = input->cth_albedo;
    //open the file
    printf( "%s, %d: I - Loading Cloud Top Height albedo file: %s\n", 
      __FILE__, __LINE__, alb_file );
    if (nc_open(alb_file, 0, &ncid) != NC_NOERR) {
      fprintf(stderr,
        "-E- %s %d: file: %s is not netcdf, not acceptable cth albedo file\n",
                    __FILE__, __LINE__, alb_file);
      return -1;
      }
    /* WDR quick, for the reduced TROPOMI file it has normal strings */
    // read the 'title', and check
    // This can't be done as the attrib is a string array so repl with 7 lines
    status = nc_get_att_text(ncid, NC_GLOBAL, "title", str_tit);
    //size_t varid = 0;
    //nc_inq_attlen( ncid, NC_GLOBAL, "title", &varid);
    //size_t attlen = 0;
    //nc_inq_attlen(ncid, NC_GLOBAL, "title", &attlen);
    //char **string_attr = (char**)malloc(attlen * sizeof(char*));
    //memset(string_attr, 0, attlen * sizeof(char*));
    //nc_get_att_string(ncid, NC_GLOBAL, "title", string_attr);
    
    //if( strncmp( string_attr[0], tit_clim, strlen(tit_clim) ) == 0 )
    if( strncmp( str_tit, tit_clim, strlen(tit_clim) ) == 0 )
      {
      printf( "%s, %d: TEMP: we found the cth albedo clim file: %s\n", 
        __FILE__, __LINE__, tit_clim );
     /* get the month index to correct time in file */
      unix2ymdhms( l1rec->scantime, &year, &mon, &day, &hh, &mm, &sec);
      ix_clim = mon - 1;
      printf( "-I- %s, %d: Using climate month # %d, index: %ld\n", __FILE__, 
        __LINE__, mon, ix_clim );
      printf( "-I- From file: %s\n", alb_file );
      }
    else
      {
      printf( "-E- %s, %d: Incorrect title in albedo dataset: %s\n",
        __FILE__, __LINE__, alb_file);
      return -1;
      }

    // get the dim sizes - these should be as expected
    if( ( nc_inq_dimid( ncid, "wavelength", dim_id ) != NC_NOERR) ||
        ( nc_inq_dimid( ncid, "month", ( dim_id + 1 ) ) != NC_NOERR) ||
        ( nc_inq_dimid( ncid, "latitude", ( dim_id + 2 ) ) != NC_NOERR) ||
        ( nc_inq_dimid( ncid, "longitude", ( dim_id + 3 ) ) != NC_NOERR) )
      {
      printf( "-E- %s, %d: Error retrieving the wave, month, lat, lon size\n", 
        __FILE__, __LINE__ );
      return -1;
      }
    if ( ( nc_inq_dimlen(ncid, dim_id[0], &nwave ) != NC_NOERR) ||
         ( nc_inq_dimlen(ncid, dim_id[1], &nmon ) != NC_NOERR) ||
         ( nc_inq_dimlen(ncid, dim_id[2], &nlat ) != NC_NOERR) ||
         ( nc_inq_dimlen(ncid, dim_id[3], &nlon ) != NC_NOERR) )
      {
      printf( "-E- %s, %d: Error retrieving the wave, month, lat, lon size\n", 
        __FILE__, __LINE__ );
      return -1;
      }
    // check for correct # wave, month, lat, lon
    if( ( nwave != CTH_NWAVE ) || ( nmon != CTH_NMON ) ||
        ( nlat != CTH_NLAT ) || ( nlon != CTH_NLON ) )
      {
      printf( "-E- %s, %d: cth albedo file dimension sizes, wave, month, lat, lon are not as expected\n",
        __FILE__, __LINE__ );
      return -1;
      }
    // Read the 18th wavelength, and the ix_clim month
    // for both albedo and uncertainty
    start[0] = ix_clim;  //month
    start[1] = 18;       // wave
    count[0] = 1;
    count[1] = 1;
    count[2] = CTH_NLON;  // lon
    count[3] = CTH_NLAT;  // lat
    // get the slice for that wavelength and month from the minimum_LER_clear 
    // and uncertainty_clear datasets
    if( ( ( alb_arr = (float *) malloc( nlat * nlon * sizeof( float ) ) ) 
       == NULL ) ||
        ( ( alb_unc_arr = (float *) malloc( nlat * nlon * sizeof( float ) ) )
       == NULL ) )
      {
      printf( "-E- %s, %d: Unable to allocate cth albedo storage\n", 
        __FILE__, __LINE__ );
      return -1;
      }
    // get each dataset
    for( ipix = 0; ipix < 2; ipix++ )
      {
      if( nc_inq_varid( ncid, arr_nms[ipix], &d_id ) != NC_NOERR )
        {
        printf( "-E- %s, %d: Unable to find dataset %s in cth albedo file\n",
          __FILE__, __LINE__, arr_nms[ipix] );
        return -1;
        }
      arr_ptr = ( ipix == 0 ) ? alb_arr : alb_unc_arr;
      if( ( status = nc_get_vara_float( ncid, d_id, start, count, arr_ptr ) ) 
        != NC_NOERR )
        {
        printf( "-E- %s, %d: Unable to read dataset %s in cth albedo file\n",
          __FILE__, __LINE__, arr_nms[0] );
        return -1;
        }
      }
    nc_close( ncid );
    }  // end of initialization
 /*
  loop through the scan and find the closest value in the albedo and unc
  */
  for( ipix = 0; ipix < l1rec->npix; ipix++ )
    {
    iy = ( l1rec->lat[ipix] + 90. ) / 0.125;
    ix = ( l1rec->lon[ipix] + 180. ) / 0.125;
    iy = ( iy < 0 )? 0 : iy;
    iy = ( iy > 1439 ) ? 1439 : iy;
    ix = ( ix < 0 ) ? 0 : ix;
    ix = ( ix > 2879 ) ? 2879 : ix;
    //  make sure albedo is < 0.99 and its unc is > 0.02
    talb = *( alb_arr + iy + 1440 * ix );
    l1rec->cld_dat->cth_alb_init[ipix] = ( talb > 0.99 ) ? 0.99 : talb;
    talb = *( alb_unc_arr + iy + 1440 * ix );
    l1rec->cld_dat->cth_alb_unc_init[ipix] = ( talb < 0.02 ) ? 0.02 : talb;
    }
 /* and finish the line */
  return 0;
  }
