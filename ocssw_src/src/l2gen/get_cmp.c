# include <stdlib.h>
# include <stdio.h>
#include "l12_proto.h"
#include "l2prod.h"
#include <netcdf.h>
#include "met_cvt.h"
#include "scene_meta.h"

extern void ch_cld_sci_( float *, int *, unsigned char *, int *,
  int32_t *, int *, int *, char * );
/*
   routine get_cmp.c - Made from tst_cld_sci.c - this will be a stub to
   read the small chimaera data set and call the chimaera code to run
   a small set of pixels through the science code.

   This will first demonstrate the union of l2gen and the chimaera code

   This is the c portion, which reads in the data from a netcdf file
   wc_file.nc and set up the data to pass to the f90 code: ch_cld_sci.f90

   I have also put in 2 netcdf convenience routines based on the ncio_grab_f_ds
   in the l2gen src area file ncio.c
*/

/* this must be global or the f90 won't see it */
struct
  {
  int npix;  /* array sizes */
  int nlin;
  int nbnd;
  int nbnd_albedo;
  int nlvl_model;
  int scan;      /* mainly used for reading same data */
  int st_samp;   /*   from MODIS L2 cloud file */
  int g_year;   /* the granule year */
  int g_day;  /* the granule day-of-year */ 
  } dim_ctl_;

  static int32_t cur_cmp_rec = -1;

int get_cmp( l2str *l2rec, int prodnum, float prod[])
/*
get_cmp will just initiate the derivation of cloud microphysical properties 
(CMP) for the line if we haven't done so already.  Then, it will pass back the 
proper property array for that line

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     l2str *           l2rec            I      l2 record of data
     int               prodnum          I      Product ID to return
     float []          prod             O      returned data - this storage is                                                 already allocated

W. Robinson, SAIC, 4 Jan 2019

  W. Robinson, 14 Jun 2019  Add the cirrus reflectance data to the radiance set
*/
  {
  int compute_cmp( l2str * );
  void get_cmp_prod_( int *, float *, int * );
  int ipix, n_prd;
  /* if proc_cloud not set, error out   */
  if( input->proc_cloud == 0 ) {
    printf( "%s, %d: Cloud products require the proc_cloud=1 set\n", __FILE__,
      __LINE__ );
    exit(1);
  }
  /* WDR this is just to force the call all the time */
  if( l2rec->l1rec->iscan != cur_cmp_rec ) 
    {
    cur_cmp_rec = l2rec->l1rec->iscan;
    if( compute_cmp( l2rec ) != 0 )
      {
      printf( "%s, %d: E - CMP computation failure\n", __FILE__, __LINE__ );
      exit(1);
      }
    }
  /*
   *  call get_cmp_prod to get the cmp desired
   */
  /* printf( "%s, %d: calling the get_cmp_prod\n", __FILE__, __LINE__ ); */
  /* WDR this is awkward */

  n_prd = 1;
  if( ( prodnum >= 480 ) && ( prodnum <= 482 ) ) n_prd = 10;
  if( prodnum == 492 ) n_prd = 10;
  get_cmp_prod_( &prodnum, prod, &n_prd );
  /*
   *  if any product value is < -900., set prodfail
   *  For now, limit to XXX_2100 prods - some better treatment for the 
   *  different cmp products is needed going forward though
   */
   if( ( prodnum == CAT_CER_2100 ) || ( prodnum == CAT_COT_2100 ) ||
       ( prodnum == CAT_COT_2100 ) )
     {
     for( ipix = 0; ipix < l2rec->l1rec->npix; ipix++ )
       if( *( prod + ipix ) < -900 ) 
         l2rec->l1rec->flags[ipix] |= PRODFAIL;
     }
  /*
   *  do the same for 1600 -> PRODWARN and 1621 -> NAVFAIL
   *  2200 -> ATMFAIL
   */
   if( ( prodnum == CAT_CER_1600 ) || ( prodnum == CAT_COT_1600 ) ||
       ( prodnum == CAT_COT_1600 ) )
     {
     for( ipix = 0; ipix < l2rec->l1rec->npix; ipix++ )
       if( *( prod + ipix ) < -900 )  
         l2rec->l1rec->flags[ipix] |= PRODWARN;
     }
   if( ( prodnum == CAT_CER_1621 ) || ( prodnum == CAT_COT_1621 ) ||
       ( prodnum == CAT_COT_1621 ) )
     {
     for( ipix = 0; ipix < l2rec->l1rec->npix; ipix++ )
       if( *( prod + ipix ) < -900 )
         l2rec->l1rec->flags[ipix] |= NAVFAIL;
     }
   if( ( prodnum == CAT_CER_2200 ) || ( prodnum == CAT_COT_2200 ) ||
       ( prodnum == CAT_COT_2200 ) )
     {
     for( ipix = 0; ipix < l2rec->l1rec->npix; ipix++ )
       if( *( prod + ipix ) < -900 )
         l2rec->l1rec->flags[ipix] |= ATMFAIL;
     }
  return(0);
  }

unsigned char *get_cmp_byt( l2str *l2rec, int prodnum )
/*
get_cmp_byt is like get_cmp but will provide the byte flag values

   Returns unsigned char array of the flag valueson the line

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     l2str *           l2rec            I      l2 record of data
     int               prodnum          I      Product ID to return

   W. Robinson, SAIC, 5 Aug 2019

*/
  {
  int compute_cmp( l2str * );
  static unsigned char *bprod = NULL;
  void get_cmp_byt_( int *, unsigned char * );
  /* WDR this is just to force the call all the time */
  if( l2rec->l1rec->iscan != cur_cmp_rec )
    {
    cur_cmp_rec = l2rec->l1rec->iscan;
    if( compute_cmp( l2rec ) != 0 )
      {
      printf( "%s, %d: E - CMP computation failure\n", __FILE__, __LINE__ );
      exit(1);
      }
    }
  
  if( bprod == NULL )
    {
    if( ( bprod = (unsigned char *) 
      malloc( l2rec->l1rec->npix * sizeof(unsigned char) ) ) == NULL )
      {
      printf( "%s, %d: E - unable to alocate byte storage\n", 
        __FILE__, __LINE__ );
      exit(1);
      }
    }
  /*
   *  call get_cmp_prod to get the cmp desired
   */
  /* printf( "%s, %d: calling the get_cmp_prod\n", __FILE__, __LINE__ ); */
  get_cmp_byt_( &prodnum, bprod );

  return( bprod );
  }

int compute_cmp( l2str *l2rec )
/*
  compute_cmp will derive the cloud microphysical properties 

*/
  {
 /*
  *  set the inputs based on a netcdf file for the test
  */
  int set_cmp( l2str *, float **, int32_t *, int32_t **, int32_t *,
    unsigned char **, int32_t * );
  void get_cmp_byt_( int *, unsigned char * );
  unsigned char *lcl_prd;
  int32_t nfloat, nint32, nubyte;
  static float *tdat;
  static unsigned char *ubdat;
  static int32_t *i32dat;
  int sensor_id = l2rec->l1rec->l1file->sensorID;
  int ip, npix = l2rec->l1rec->npix;
  int cat_val = CAT_Cld_Phase_2100;

 /*
  *  assure that the cloud top height parmeters are made for this line
  */
  get_ctht_lin( l2rec, -1 );
  //printf( "%s, %d: Just created the ctht prods for this scan\n", __FILE__, __LINE__ );
 /*
  *  Fill with data from the l1que
  */
  if( set_cmp( l2rec, &tdat, &nfloat, &i32dat, &nint32, &ubdat, &nubyte ) != 0 )
    return(-1);
 /*
  *  call the chimaera code
  */
   /*  OK, call the f90 routine */
  ch_cld_sci_( tdat, &nfloat, ubdat, &nubyte, i32dat, &nint32, &sensor_id,
    input->cloud_hgt_file );
/*
 *  OK, the cloud phase in 2.1 will be used to set the l2_flags as such:
 *  TURBIDW if phase is 2 or 4 (water) and COCCOLITH if phase is 3 (ice)
 */
  if( ( lcl_prd = (unsigned char *)malloc( npix * sizeof( char ) ) ) == NULL )
    {
    printf( "%s, %d: Unable to allocate cloud phase buffer\n", __FILE__, __LINE__ );
    exit(1);
    }
 /*
  *  lastly, to set the phase at the pixels, do the following.  The cloud 
  *  phase 2100 applies as cloud phase for all the absorbing bands
  */
  get_cmp_byt_( &cat_val, lcl_prd );

  for( ip = 0; ip < npix; ip++ )
    {
    if( *( lcl_prd + ip ) == 3 )
      l2rec->l1rec->flags[ip] |= COCCOLITH;
    if( ( *( lcl_prd + ip ) == 2 ) || ( *( lcl_prd + ip ) == 4 ) )
      l2rec->l1rec->flags[ip] |= TURBIDW;
    }

  free(lcl_prd);

  return(0);
  }

int set_cmp( l2str *l2rec, float **tdat, int32_t *nfloat, int32_t **i32dat, 
  int32_t *nint32, unsigned char **ubdat, int32_t *nubyte )
/*
  set_cmp  will set up all the data required by the chimaera code from
  the contents of the l1rec
  30 Jun 2023, W Robinson, add the cloud top height arrays in (ctht)

*/
  {
  extern l1qstr l1que;
  extern ctht_lins_str ctht_lins;
  int mk_cmp_prof( float *, float *, float *, float *, float, float, float, 
    float, unsigned char *, unsigned char *, double *, double *, double *, 
    double * );
/* WDR **** temp to dummy the profile interp */
  int mk_cmp_prof_dum( float *, float *, float *, float *, float, float, float, 
    float, unsigned char *, unsigned char *, double *, double *, double *, 
    double * );
  int make_profile_101_( double * );
  int32_t npix, nlin, nbnd_albedo, ctr_que,ibnd, iln, qln, ipx;
  l1str *l1rec = l2rec->l1rec;
  int32_t nbnd_ref, nbnd_emis, foff, foff2, uoff, uoff2, nbnd, sensor_id;
  float *F0, rad, solz, out;
  int32_t dbg_print, ix_bnd, force_oci, cld_missed;
  int32_t ncmp_bnd = 15;
 /* WDR a temp print routine for debug */
  int profile_print( float *, unsigned char *, int32_t *, int32_t, int32_t, 
    int32_t, int32_t );

 /* some look-ups to id the type / location  of the l2 rads relative to the 
    chimaera rads */
 /* source for the chimaera rads: 1 reflective (Lt), 0 emmissive (Ltir)
    and cirrus (rho_cirrus) */
  static int32_t ref_for_cmb[] = { 1, 1, 1, 1, 1, 1, 0, 2, 0, 0, 1, 0, 0, 1, 1 };
 /* cloud band is in l1rec if 1, else 0 */
  static char cmp_there[15];
 /* index in Lt (ref) or Ltir (emis) (after # vis bands subtracted) */
  static int32_t l1ix[15];
 /* cloud band ideal wavelengths for MODIS (and any other inst to try) */
  int32_t cld_wav_mod[] = { 645, 859, 1240, 1640, 2130, 930, 3750, 1375, 
    8550, 11000, 443, 12000, 13600, 469, 555 };
 /* cloud bands for OCI (only) */
  // int32_t cld_wav_oci[] = { 645, 859, 1250, 1615, 2130, 930, 3750, 1375, 
  //   8550, 11000, 443, 12000, 13600, 469, 555 };
 /* 10 Sep 21 to match the waves for the Kerry Meyer table */
 /* 29 aug 22 the 8.8 um not in OCI, replace with the 2.26 um */
 /* WDR pre-2.2 band use in OCI
  int32_t cld_wav_oci[] = { 665, 865, 1250, 1616, 2130, 930, 3750, 1375,
    2260, 11000, 443, 12000, 13600, 469, 555 };
  */
 /* WDR 18 Sep 23 set cld_wav_oci as cld_wav_ocis and make cld_wav_oci
    the actual instrument wavelengths
  int32_t cld_wav_oci[] = { 665, 865, 1250, 1616, 2130, 930, 2260, 1375,
    8550, 11000, 443, 12000, 13600, 469, 555 };
  */



  //  for inst OCI
  int32_t cld_wav_oci[] = { 665, 865, 1249, 1618, 2131, 940, 2258, 1378,
    8550, 11000, 442, 12000, 13600, 470, 555 };
  int32_t cld_wav[15];
 /* OCI bands available (for simulation) */
  int32_t oci_bnd_avail[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1 };
 /* minimum bands required for cloud processing */
  static int32_t cld_min_bnd[] = { 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1 };
 /* The following emperically derived factors will correct the l2gen radiances
    and reflectances to those that chimaera expects - probably the difference 
    between OBPG and EOSDIS calibration */
 //  *** I have temporarily removed the can change for a test
  //float chim_rad_corr[] = { 0.965407, 0.973168, 1.00678, 1.00244, 0.963511,
  float chim_rad_corr[] = { 1.0,      1.0,      1.0,     1.0,     1.0,
  //  WDR Jan 2021 try 5% low in 1st 5 bands
  //float chim_rad_corr[] = { 0.5,    0.5,      0.5,    0.5,   0.5,
  /* bands                  645       859       1240     1640     2130  */
    //-1.,      0.999840, -1.,     0.999979, 0.999794, 1.01814, 0.999760,
    -1.,      1.0,      1.0,     1.0,      1.0,      1.0,     1.0,
 /*  930 approx   3750   1375 approx  8550   11000     412      12000 */
    //-1.,           0.986288, 0.994909 };
    -1.,           1.0,      1.0      };
 /*  13600 approx  469      555   */

 /*  Guard the transfer array creation to only 1 per run */
  static int32_t firstcall = 0;
 /* 
  *  for the met data,
  *  set the merra profile heights in reverse order to match the order in the
  *  profiles used by chimaera
  */
  int32_t nlvl = 42, nlvl_cmp = 101, ilvl, is_land, cld_flg;
  int32_t bad_rad, ice_flg, glint_flg, n_sfc_albedo;
  int32_t npix_scan, lst, len, pst, pen, apx, aln, itot, nland;
  unsigned char sfc_lvl, trop_lvl;
  float merra_p_set[] = { 0.1, 0.3, 0.4, 0.5, 0.7, 1., 2., 3., 4., 5., 7., 10., 
      20., 30., 40., 50., 70., 100., 150., 200., 250., 300., 350., 400., 
      450., 500., 550., 600., 650., 700., 725., 750., 775., 800., 825., 
      850., 875., 900., 925., 950., 975., 1000. };
  float merra_t[nlvl], merra_p[nlvl], merra_q[nlvl], merra_h[nlvl];
  double cmp_t[nlvl_cmp], cmp_mixr[nlvl_cmp], cmp_h[nlvl_cmp];
  double cmp_p[nlvl_cmp], cmp_p_set[nlvl_cmp];
  scnstr *meta;
 /*
  get the # pixels
  and for now, make sure they match the # coming in
  */
  ctr_que = l1que.nq / 2;  /* the center of the queue = our primary line */

  npix = l1rec->npix;
  // WDR test nlin = 3;
  nlin = 3;
  dbg_print = 0;

  dim_ctl_.npix = npix;
  dim_ctl_.nbnd_albedo = 6;  /* the MODIS # 'albedo bands' */
  nbnd_albedo = dim_ctl_.nbnd_albedo;
  dim_ctl_.nlin = nlin;
  dim_ctl_.nbnd = ncmp_bnd;  /* a bnd # used from MODIS */
  nbnd = ncmp_bnd;
  dim_ctl_.nlvl_model = nlvl_cmp;  /* the set # interpolated to for chimaera */
  dim_ctl_.scan = l1que.r[ctr_que].iscan;
  dim_ctl_.st_samp = 1; /* probably not used, but ftn 1-origin start samp */
  meta = scene_meta_get();
  /* get the time info down to the day - all that is used */
  dim_ctl_.g_year = meta->start_year;
  dim_ctl_.g_day = meta->start_day;

  sensor_id = l1rec->l1file->sensorID;
  n_sfc_albedo = 5;
  if( sensor_id == OCI ) n_sfc_albedo = 6;

 /*  Allocate arrays to pass to the f90 routine with all the information
    together: all real arrays into a real array etc  */
  *nfloat = npix * nlin * nbnd +        /* reflectance */
    npix * nlin * nbnd_albedo +  /* band uncertainty */
    2 * nbnd_albedo +   /* uncertainty info */
    7 * npix * nlin +   /* geom info */
    4 * npix * nlin * nlvl_cmp +  /* met profiles */
    8 * npix * nlin +   /* 2D anc data */
    n_sfc_albedo * npix * nlin +   /* surface albedo for 5/6 bands */
    3 * npix * nlin;   /* for the 3 arrays of cloud top height, pressure, and 
                          temperature   */

  nbnd_ref = l1rec->l1file->nbands;
  nbnd_emis = l1rec->l1file->nbandsir;
  if( firstcall == 0 )
    {
    firstcall = 1;
    if( ( *tdat = (float *) malloc( *nfloat * sizeof(float) ) ) == NULL )
      {
      printf( "%s - %d: Allocation of real storage failed\n", 
        __FILE__, __LINE__ );
      exit(3);
      }
   /*
    * : all byte arrays
    */
    *nubyte = 2 * npix * nlin +    /*  for the ancillary byte stuff */
              23 * npix * nlin +    /* for the cloud mask stuff */
              ncmp_bnd;            /*  for the band availability array */
    if( ( *ubdat = (unsigned char *) 
      malloc( *nubyte * sizeof(unsigned char) ) ) == NULL )
      {
      printf( "%s - %d: Allocation of byte storage failed\n",
        __FILE__, __LINE__ );
      exit(3);
      }
   /*
    *  : all int32 arrays
    */
    *nint32 = npix * nlin;
    if( ( *i32dat = (int32_t *) malloc( *nint32 * sizeof(int32_t) ) ) == NULL )
      {
      printf( "%s - %d: Allocation of i32 storage failed\n",
        __FILE__, __LINE__ );
      exit(3);
      }
   /*
    *  Set up the cloud band presence and L1 index here
    */
    for( ibnd = 0; ibnd < ncmp_bnd; ibnd++ )
      {
      if( sensor_id == OCI )
        cld_wav[ibnd] = cld_wav_oci[ibnd];
      else
        cld_wav[ibnd] = cld_wav_mod[ibnd];
      }
   /*  if OCI, make the 2.2 reflective band  */
    if( sensor_id == OCI ) 
      {
      ref_for_cmb[6] = 1;
      cld_min_bnd[6] = 1;
      }
    for( ibnd = 0; ibnd < ncmp_bnd; ibnd++ )
      {
      ix_bnd = bindex_get( cld_wav[ibnd] );
      if( ( ix_bnd < 0 ) && ( ref_for_cmb[ibnd] != 2 ) )  // if no matching, 
                                                // non-cirrus band was found
        {
        cmp_there[ibnd] = 0;
        l1ix[ibnd] = -1;
        }
      else
        {
        cmp_there[ibnd] = 1;
        if( ref_for_cmb[ibnd] == 0 )  // emissive
          l1ix[ibnd] = ix_bnd;
        else if( ref_for_cmb[ibnd] == 1 )  // reflective
          l1ix[ibnd] = ix_bnd;
        else  //cirrus
          l1ix[ibnd] = 0;
        }
      }
   /*
    *  This code will optionally force the available bands down to the OCI bands
    */
    force_oci = 0;
    printf( "OCI band forcing is set to %d\n", force_oci );
    if( force_oci == 1 )
      {
      for( ibnd = 0; ibnd < ncmp_bnd; ibnd++ )
        cmp_there[ibnd] *= oci_bnd_avail[ibnd];
      }
   /*
    *  If the minimum band set for 2.1, (2.2, OCI) and 1.6 nm cloud products 
    *  (cld_min_bnd) is not in this instrument (cmp_there), report it and 
    *  stop processing here
    */
    printf( "Talley of current instrument's wavelengths\n" );
    printf( "Cloud        wave is     wave is\n" );
    printf( "wavelength   available   required\n" );
    cld_missed = 0;
    for( ibnd = 0; ibnd < ncmp_bnd; ibnd++ )
      {
      if( ( cld_min_bnd[ibnd] == 1 ) && ( cmp_there[ibnd] == 0 ) )
        cld_missed = 1;
      printf( "%10d   %9d   %8d\n", cld_wav[ibnd], cmp_there[ibnd],
        cld_min_bnd[ibnd] );   
      }
    if( cld_missed == 1 )
      {
      printf( "%s, %d: The required wavelengths for cloud processing were not found, Exiting\n", __FILE__, __LINE__ );
      exit(FATAL_ERROR);
      }
    }

  F0 = l1rec->Fo;

 /* go through the chimaera refl (refl) rad (emiss) bands */
  for( ibnd = 0; ibnd < ncmp_bnd; ibnd++ )
    {
    for( iln = 0; iln < nlin; iln++ )
      {
      qln = iln + ctr_que - nlin / 2;
      for( ipx = 0; ipx < npix; ipx++ )
        {
       /* only write the center line value */
        if( ( iln == 1) && ( dbg_print == 1 ) )
          {
          printf( "%s, %d: bnd %d, lin %d, pix %d, rad %f\n", 
            __FILE__, __LINE__, ibnd, iln, ipx, 
            (*tdat)[ ipx + npix * ( ibnd + ncmp_bnd * iln ) ] );
          }
       /* only substitute if the band exists (for now) */
        if( cmp_there[ibnd] == 1 )
          {
         /* for the reflective */
          if( ref_for_cmb[ibnd] == 1 )
            {
            rad = l1que.r[qln].Lt[ l1ix[ ibnd ] + nbnd_ref * ipx ];
            solz = l1que.r[qln].solz[ipx ];
            out = rad * OEL_PI / ( cos( solz * OEL_PI / 180. ) * F0[ l1ix[ ibnd ] ] );
            (*tdat)[ ipx + npix * ( ibnd + ncmp_bnd * iln ) ] = 
               chim_rad_corr[ibnd] * out;
            }
          else if( ref_for_cmb[ibnd] == 2 )  /* for the cirrus */
            {
            solz = l1que.r[qln].solz[ipx ];
            rad = l1que.r[qln].rho_cirrus[ipx];
            (*tdat)[ ipx + npix * ( ibnd + ncmp_bnd * iln ) ] =
              chim_rad_corr[ibnd] * rad * cos( solz * OEL_PI / 180. );
            }
          else  /* for emissive */
            {
            (*tdat)[ ipx + npix * ( ibnd + ncmp_bnd * iln ) ] = 
              chim_rad_corr[ibnd] * 
              10. * l1que.r[qln].Ltir[ ( l1ix[ ibnd ] - nbnd_ref ) + 
              nbnd_emis * ipx ];
            }
          }
        /* for the missing... try -999.0 */
        else
          {
          (*tdat)[ ipx + npix * ( ibnd + ncmp_bnd * iln ) ] = -999.0;
          }
       /*  repeat write of the center line value */
         if( ( iln == 1 ) && (dbg_print == 1 ) )
          {
          printf( "%s, %d: bnd %d, lin %d, pix %d, rad %f\n", 
            __FILE__, __LINE__, ibnd, iln, ipx,
            (*tdat)[ ipx + npix * ( ibnd + ncmp_bnd * iln ) ] );
          }
        }
      }
    }
 /*
  *  the uncertainty with 6 band_albedo - I have a little on this:
  *  There is a equation in modis_science_module.f90 that shows how the
  *  info is used:
  *  refl_unc = spec_uncertain * exp(band_uncertainty) / uncertain_sf
  *  This is a kind of odd exponential relation.
  *  It also shows that I might be better using 1. for the dummy uncertain_sf
  *  value
  *
  *  so 0 values are OK for now except the uncertain_sf (used in denominator)
  *
  *   *** NOTE *** that the band_uncertainty in the chimaera code is 
  *   actually a unsigned byte or i*1.  Now, assigning the real data to that
  *   in the cloud code should work just fine - just some overkill in the 
  *   value description.  However, making this data a char array would be
  *   work for very little gain and could make problems if done wrong.  It 
  *   could also hide the problem we're looking for too.  So, I'm leaving
  *   this as-is for now.
  */
  foff = npix * nlin * ncmp_bnd; /* offset to the output data */
  for( ibnd = 0; ibnd < nbnd_albedo; ibnd++ )
    {
    for( iln = 0; iln < nlin; iln++ )
      {
      qln = iln + ctr_que - nlin / 2;
      for( ipx = 0; ipx < npix; ipx++ )
        {
       /* only write the center line value */
        if( ( iln == 1 ) && (dbg_print == 1 ) )
          {
          printf( "%s, %d: bnd %d, lin %d, pix %d, UNC %f\n", 
            __FILE__, __LINE__, ibnd, iln, ipx, 
            (*tdat)[ foff + ipx + npix * ( ibnd + nbnd_albedo * iln ) ] );
          }
       /* substitute a 0 for all */
        (*tdat)[ foff + ipx + npix * ( ibnd + nbnd_albedo * iln ) ] = 0.;

       /*  repeat write of the center line value */
        if( ( iln == 1 ) && (dbg_print == 1 ) )
          {
          printf( "%s, %d: bnd %d, lin %d, pix %d, UNC %f\n",
            __FILE__, __LINE__, ibnd, iln, ipx, 
            (*tdat)[ foff + ipx + npix * ( ibnd + nbnd_albedo * iln ) ] );
          }
        }
      }
    }
 /*
  *  also the nbnd_albedo long arrays of spec_uncertain and uncertain_sf see 
  *  above
  */
  foff += npix * nlin * nbnd_albedo;
  foff2 = foff + nbnd_albedo;
  for( ibnd = 0; ibnd < nbnd_albedo; ibnd++ )
    {
    (*tdat)[ foff + ibnd ] = 0.;
    (*tdat)[ foff2 + ibnd ] = 1.;
    }
 /*
  *  transfer the geolocation and view angle information over
  */
  foff = foff2 + nbnd_albedo;
  foff2 = npix * nlin;
  /* lat and lon at the l1rec values */
  for( iln = 0; iln < nlin; iln++ )
    {
    qln = iln + ctr_que - nlin / 2;
    for( ipx = 0; ipx < npix; ipx++ )
      {
      if( ( iln == 1 ) && ( dbg_print == 1 ) )
        {
        printf( "%s, %d: lin %d, pix %d, INITIAL lon %f, lat %f, senz %f, sena %f, solz %f, sola %f, relaz %f\n", 
          __FILE__, __LINE__, iln, ipx, 
          (*tdat)[ foff + ipx + npix * iln ], (*tdat)[ foff + foff2 + ipx + npix * iln ],
          (*tdat)[ foff + 2 * foff2 + ipx + npix * iln ], 
          (*tdat)[ foff + 3 * foff2 + ipx + npix * iln ],
          (*tdat)[ foff + 4 * foff2 + ipx + npix * iln ], 
          (*tdat)[ foff + 5 * foff2 + ipx + npix * iln ], 
          (*tdat)[ foff + 6 * foff2 + ipx + npix * iln ] );
        }
     /* lat and lon */
      out = l1que.r[qln].lat[ipx];
      (*tdat)[ foff + ipx + npix * iln ] = out;
      out = l1que.r[qln].lon[ipx];
      (*tdat)[ foff + foff2 + ipx + npix * iln ] = out;
     /* senz, sena, sola, solz, relaz */
      out = l1que.r[qln].senz[ipx];
      (*tdat)[ foff + 2 * foff2 + ipx + npix * iln ] = out;
      out = l1que.r[qln].sena[ipx];
      (*tdat)[ foff + 3 * foff2 + ipx + npix * iln ] = out;
      out = l1que.r[qln].sola[ipx];
      (*tdat)[ foff + 4 * foff2 + ipx + npix * iln ] = out;
      out = l1que.r[qln].solz[ipx];
      (*tdat)[ foff + 5 * foff2 + ipx + npix * iln ] = out;
      out = l1que.r[qln].delphi[ipx];
      (*tdat)[ foff + 6 * foff2 + ipx + npix * iln ] = -out;  /* their relaz is (-) ours */

      if( ( iln == 1 ) && ( dbg_print == 1 ) )
        {
        printf( "%s, %d: lin %d, pix %d, FINAL lon %f, lat %f, senz %f, sena %f, solz %f, sola %f, relaz %f\n",
          __FILE__, __LINE__, iln, ipx,
          (*tdat)[ foff + ipx + npix * iln ], (*tdat)[ foff + foff2 + ipx + npix * iln ], 
          (*tdat)[ foff + 2 * foff2 + ipx + npix * iln ], 
          (*tdat)[ foff + 3 * foff2 + ipx + npix * iln ],
          (*tdat)[ foff + 4 * foff2 + ipx + npix * iln ], 
          (*tdat)[ foff + 5 * foff2 + ipx + npix * iln ],
          (*tdat)[ foff + 6 * foff2 + ipx + npix * iln ] );
        }
      }
    }
 /*
  *  next is the met info
  */
  foff += 7 * foff2;  /* start of the profiles */
 /* set up the clean hi res profile */
  make_profile_101_( cmp_p_set );

 /*  check that profiles were read in */
  if( l1que.r[ctr_que].anc_add == NULL )
    {
    fprintf( stderr, 
      "-E- %s %d: Ancillary profile required for cloud processing is missing\n",
      __FILE__, __LINE__);
    exit(1);
    }
  uoff = 0;
  for( iln = 0; iln < nlin; iln++ )
    {
    qln = iln + ctr_que - nlin / 2;
    for( ipx = 0; ipx < npix; ipx++ )
      {
      memcpy( merra_p, merra_p_set, nlvl * sizeof(float) );
      for( ilvl = 0; ilvl < nlvl; ilvl++ )
        {
       /*  note that the profiles for chimaera are reverse of the MERRA,
           so in making the merra arrays from the l1rec stuff, we reverse it */
        merra_t[ nlvl - 1 - ilvl ] = l1que.r[qln].anc_add->prof_temp[ ilvl + nlvl * ipx ];
        merra_q[ nlvl - 1 - ilvl ] = l1que.r[qln].anc_add->prof_q[ ilvl + nlvl * ipx ];
        merra_h[ nlvl - 1 - ilvl ] = l1que.r[qln].anc_add->prof_height[ ilvl + nlvl * ipx ];
        }
/* WDR trad test point  tpix = 2, tlin = 1; */
int tpix = 2, tlin = 1;
if( (ipx == tpix) && ( iln == tlin ) && ( dbg_print == 1 ) )
  {
      printf( "\n\n%s, %d: merra_p for pix: %d, lin: %d\n", __FILE__, __LINE__, ipx, iln );
      for( ilvl = 0; ilvl < nlvl; ilvl++ )
        printf( "%f ", merra_p[ilvl] );
      printf( "\n" );
  
      printf( "%s, %d: merra_t for pix: %d, lin: %d\n", __FILE__, __LINE__, ipx, iln );
      for( ilvl = 0; ilvl < nlvl; ilvl++ )
        printf( "%f ", merra_t[ilvl] );
      printf( "\n" );
  
      printf( "%s, %d: merra_q for pix: %d, lin: %d\n", __FILE__, __LINE__, ipx, iln );
      for( ilvl = 0; ilvl < nlvl; ilvl++ )
        printf( "%f ", merra_q[ilvl] );
      printf( "\n" );
  
      printf( "%s, %d: merra_h for pix: %d, lin: %d\n", __FILE__, __LINE__, ipx, iln );
      for( ilvl = 0; ilvl < nlvl; ilvl++ )
        printf( "%f ", merra_h[ilvl] );
      printf( "\n" );
  }

     /*
      *  for this profile, make the 101 level chimaera profile
      *  and put the chimaera profiles into the storage
      */
      memcpy( cmp_p, cmp_p_set, nlvl_cmp * sizeof( double) );
/*  WDR **** test branch to get around interpolation of merra */
/* ilvl = -899;   WDR use a pre-computed fixed profile */
ilvl = 0;  /* WDR use the real profiles */
if( ilvl == -899 ) 
  {
  mk_cmp_prof_dum( merra_p, merra_t, merra_q, merra_h, l1que.r[qln].sfcp[ipx],
        l1que.r[qln].sfct[ipx],  l1que.r[qln].sfcrh[ipx],
        l1que.r[qln].lat[ipx], &sfc_lvl, &trop_lvl, cmp_p, cmp_t,
        cmp_mixr, cmp_h );
  }
else
  {
      mk_cmp_prof( merra_p, merra_t, merra_q, merra_h, l1que.r[qln].sfcp[ipx], 
        l1que.r[qln].sfct[ipx],  l1que.r[qln].sfcrh[ipx], 
        l1que.r[qln].lat[ipx], &sfc_lvl, &trop_lvl, cmp_p, cmp_t, 
        cmp_mixr, cmp_h );
  }
     /* as sfc_lvl, trop_lvl go to ftn, make them 1-origin */
      sfc_lvl = ( sfc_lvl == nlvl_cmp ) ? sfc_lvl : sfc_lvl + 1;
      trop_lvl = ( trop_lvl == nlvl_cmp ) ? trop_lvl : trop_lvl + 1;

     /*
      *  WDR Look through the profiles and see if all is OK with them.
      *  out, no need for debug
      *** 5 Mar 2024 I may check this below in mk_cmp_prof
      *** but comment out for now
      */
/*
      for( ilvl = 0; ilvl < nlvl_cmp; ilvl++ )
        {
        if( ( cmp_t[ilvl] > 350. ) || ( cmp_t[ilvl] < 100. ) )
          {
          printf( "PROFILE ODD T: %f, iln: %d, ipx: %d, ilvl: %d\n", 
            cmp_t[ilvl], iln, ipx, ilvl );
          }
        if( ( cmp_mixr[ilvl] > 100. ) || ( cmp_mixr[ilvl] < 0 ) )
          {
          printf( "PROFILE ODD MIXR: %f, iln: %d, ipx: %d, ilvl: %d\n",
            cmp_mixr[ilvl], iln, ipx, ilvl );
          }
        if( ( cmp_h[ilvl] > 1000000. ) || ( cmp_h[ilvl] < -20000. ) )
          {
          printf( "PROFILE ODD height: %f, iln: %d, ipx: %d, ilvl: %d\n",
            cmp_h[ilvl], iln, ipx, ilvl );
          }
        if( ( cmp_p[ilvl] > 1200. ) || ( cmp_p[ilvl] < 0. ) )
          {
          printf( "PROFILE ODD P: %f, iln: %d, ipx: %d, ilvl: %d\n",
            cmp_p[ilvl], iln, ipx, ilvl );
          }
        }
  end commenting */
/*  the tdat at foff starts with mixr 
    it runs all pixels fastest, then lines, themn levels
    so t is at foff + npix * nlin * nlvl_cmp
       p is at foff + 2 * nlin * nlvl_cmp
       h is  at foff + 3 * nlin * nlvl_cmp

       the offset for pix, lin level is the base above plus
       ipx * npix * ( iln + nlin * ilvl )
*/
/* I need to print the profiles before they get re-done and then again after 
   they are updated.  Maybe a routine for that would be usefull */
if( ( ipx == tpix ) && ( iln == tlin ) && ( dbg_print == 1 ) ) 
  {
  printf( "/n%s, %d: mix, t, p, h profiles BEFORE, pix %d, lin %d\n", __FILE__, __LINE__, ipx, iln );
  profile_print( (*tdat) + foff, (*ubdat), (*i32dat), ipx, iln, npix, nlin );
  printf( "OUR sfcp: %f, sfct: %f\n", l1que.r[qln].sfcp[ipx], l1que.r[qln].sfct[ipx] );
   printf( "mk_cmp_prof sfc_lvl: %d, trop_lvl: %d\n", sfc_lvl, trop_lvl );
  }
      for( ilvl = 0; ilvl < nlvl_cmp; ilvl++ )
        {
        (*tdat)[foff + ipx + npix * ( iln + nlin * ilvl ) ] = cmp_mixr[ilvl];
        foff2 = foff + npix * nlin * nlvl_cmp;
        (*tdat)[foff2 + ipx + npix * ( iln + nlin * ilvl ) ] = cmp_t[ilvl];
        foff2 = foff + 2 * npix * nlin * nlvl_cmp;
        (*tdat)[foff2 + ipx + npix * ( iln + nlin * ilvl ) ] = cmp_p[ilvl];
        foff2 = foff + 3 * npix * nlin * nlvl_cmp;
        (*tdat)[foff2 + ipx + npix * ( iln + nlin * ilvl ) ] = cmp_h[ilvl];
        }
     /*
      *  add the single-level met parameters
      *  sfc T
      */
      foff2 = foff + 4 * npix * nlin * nlvl_cmp;
      (*tdat)[foff2 + ipx + npix * iln ] = 273.15 + l1que.r[qln].sfct[ipx];
     /* sfc P */
      foff2 += npix * nlin;
      (*tdat)[foff2 + ipx + npix * iln ] = l1que.r[qln].sfcp[ipx];
     /* wind speed */
      foff2 += npix * nlin;
      (*tdat)[foff2 + ipx + npix * iln ] = l1que.r[qln].ws[ipx];
     /* O3 */
      foff2 += npix * nlin;
      (*tdat)[foff2 + ipx + npix * iln ] = 1.e3 * l1que.r[qln].oz[ipx];
     /* ice frac and snow frac */
      is_land = l1que.r[qln].land[ipx];
/*
if( ( ipx == 2 ) && ( iln == 1 ) )
printf( "flag: %d\n", l1que.r[qln].flags[ipx]);
*/
      foff2 += npix * nlin;
      (*tdat)[foff2 + ipx + npix * iln ] = ( is_land ) ? 0. : 
        l1que.r[qln].icefr[ipx];

      foff2 += npix * nlin;
      (*tdat)[foff2 + ipx + npix * iln ] = ( is_land ) ? 
        l1que.r[qln].icefr[ipx] : 0;
     /* alt o3 - originally an alt o3 source but... */
      foff2 += npix * nlin;
      (*tdat)[foff2 + ipx + npix * iln ] = 1.e3 * l1que.r[qln].oz[ipx];
     /* alt ice conc */
      foff2 += npix * nlin;
     /*  WDR to globalize this value
      (*tdat)[foff2 + ipx + npix * iln ] = ( is_land ) ? 0. :
        l1que.r[qln].icefr[ipx];
     */
      (*tdat)[foff2 + ipx + npix * iln ] = l1que.r[qln].icefr[ipx];
     /*  We will get the surface albedo from l2gen now */
      foff2 += npix * nlin;
      (*tdat)[foff2 + ipx + npix * iln ] = 
        l1que.r[qln].cld_dat->sfc_albedo_659[ipx];
      foff2 += npix * nlin;
      (*tdat)[foff2 + ipx + npix * iln ] =
        l1que.r[qln].cld_dat->sfc_albedo_858[ipx];
      foff2 += npix * nlin;
      (*tdat)[foff2 + ipx + npix * iln ] =
        l1que.r[qln].cld_dat->sfc_albedo_1240[ipx];
      foff2 += npix * nlin;
      (*tdat)[foff2 + ipx + npix * iln ] =
        l1que.r[qln].cld_dat->sfc_albedo_1640[ipx];
      foff2 += npix * nlin;
      (*tdat)[foff2 + ipx + npix * iln ] =
        l1que.r[qln].cld_dat->sfc_albedo_2130[ipx];
     /* WDR for 2.2 as we don't have 2.2, duplicate the 2.1 here */
      if( sensor_id == OCI  )
        {
        foff2 += npix * nlin;
        (*tdat)[foff2 + ipx + npix * iln ] =
          l1que.r[qln].cld_dat->sfc_albedo_2130[ipx];
        }
     /*  add the cloud top height, pressure, and temperature */
      foff2 += npix * nlin;
      (*tdat)[foff2 + ipx + npix * iln ] =
        ctht_lins.ct[qln]->cth[ipx];
      foff2 += npix * nlin;
      (*tdat)[foff2 + ipx + npix * iln ] =
        ctht_lins.ct[qln]->ctp[ipx];
      foff2 += npix * nlin;
      (*tdat)[foff2 + ipx + npix * iln ] =
        ctht_lins.ct[qln]->ctt[ipx];
     /* end of float data set up */

     /* sfc_lvl and trop_lvl gotten with the profile above */
     /* note that this is unsigned char data now */
     /* also note that the ftn expects 1-origin levels so add 1 */
      (*ubdat)[ uoff + ipx + npix * iln ] = sfc_lvl;
      uoff2 = uoff + npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = trop_lvl;

     /* the I32 met info is the next */
     /* just try inserting the ice fraction as percent */
      (*i32dat)[ ipx + npix * iln ] = 100 * l1que.r[qln].icefr[ipx];

     /* the remaining byte data is all from the cloud mask that is
        a MODIS product.  We have one thing: the cloud mask
        and we'll see how that does... */
     /* some set-up and CLD_DET */
      uoff2 += npix * nlin;
      bad_rad = ( l1que.r[qln].Lt[l1ix[0]+nbnd_ref * ipx] == BAD_INT ) ? 1 : 0;
      cld_flg = l1que.r[qln].cloud[ipx];
      ice_flg = l1que.r[qln].ice[ipx];
      glint_flg = l1que.r[qln].glint[ipx];
      npix_scan = l1que.r[qln].npix;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = ( bad_rad == 0 ) ? 1 : 0;
     // We will use the cloud mask file categories if the file was specified
      if( l1_input->cld_msk_file[0]) {
          char cld_category;
          get_sdps_cld_mask( &(l1que.r[qln]), ipx, &cld_category );
          /* for missing cloud, set bad_rad, no cloud detected and make clear */
          if( cld_category == BAD_BYTE ){
              bad_rad = 1;
              cld_category = 3;
              (*ubdat)[ uoff2 + ipx + npix * iln ] = 0;
          }
          /* CONF_CLD, CLR_66, CLR_95, CLR_99 */
          uoff2 += npix * nlin;
          (*ubdat)[ uoff2 + ipx + npix * iln ] = ( cld_category == 0 ) ? 1 : 0;
          uoff2 += npix * nlin;
          (*ubdat)[ uoff2 + ipx + npix * iln ] = ( cld_category == 1 ) ? 1 : 0;
          uoff2 += npix * nlin;
          (*ubdat)[ uoff2 + ipx + npix * iln ] = ( cld_category == 2 ) ? 1 : 0;
          uoff2 += npix * nlin;
          (*ubdat)[ uoff2 + ipx + npix * iln ] = ( cld_category >= 3 ) ? 1 : 0;
      } else {
         // use the binary cloud flag
         /* CONF_CLD, CLR_66, CLR_95, CLR_99 */
          uoff2 += npix * nlin;
          (*ubdat)[ uoff2 + ipx + npix * iln ] = ( cld_flg != 0 ) ? 1 : 0;
          uoff2 += npix * nlin;
          (*ubdat)[ uoff2 + ipx + npix * iln ] = ( cld_flg != 0 ) ? 1 : 0;
          uoff2 += npix * nlin;
          (*ubdat)[ uoff2 + ipx + npix * iln ] = ( cld_flg != 0 ) ? 1 : 0;
          uoff2 += npix * nlin;
          (*ubdat)[ uoff2 + ipx + npix * iln ] = 0;
      }
     /* SNO_SFC */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = 
        ( ( is_land == 1 ) && ( ice_flg == 1 ) ) ? 1 : 0;
     /* WTR_SFC */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = ( is_land == 0 ) ? 1 : 0;
     /* COAST_SFC */
      uoff2 += npix * nlin;
      lst = ( iln < 1 ) ? iln : iln - 1;
      len = ( iln >= ( nlin - 1 ) ) ? iln : iln + 1;
      pst = ( ipx <= 0 ) ? 0 : ipx -1;
      pen = ( ipx >= ( npix_scan - 1 ) ) ? ipx : ipx + 1;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = 0;
      itot = 0; 
      nland = 0;
      for( aln = lst; aln <= len; aln++ )
        for( apx = pst; apx <= pen; apx++ )
          {
          itot++;
          if( l1que.r[aln].land[apx] ) nland++;
          }
      if( ( nland != 0 ) && ( nland != itot ) )
        (*ubdat)[ uoff2 + ipx + npix * iln ] = 1;
     /* DESERT_SFC */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = 0;
     /* LND_SFC */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = ( is_land ) ? 1 : 0;
     /* NIGHT - solz > 90 */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = 
        ( l1que.r[iln].solz[ipx] > 90. ) ? 1 : 0;
     /* GLINT */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = 
        ( glint_flg ) ? 1 : 0;
     /* OCEAN_NO_SNOW */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] =
        ( ( is_land == 0 ) && ( ice_flg == 0 ) ) ? 1 : 0;
     /* OCEAN_SNOW */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = 
        ( ( is_land == 0 ) && ( ice_flg == 1 ) ) ? 1 : 0;
     /* LND_NO_SNOW */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] =
        ( ( is_land == 1 ) && ( ice_flg == 0 ) ) ? 1 : 0;
     /* LND_SNOW */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] =
        ( ( is_land == 1 ) && ( ice_flg == 1 ) ) ? 1 : 0;
     /* TST_H_138 */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = 0;
     /* TST_VIS_REFL */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = 0;
     /* TST_VIS_RATIO */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = 0;
     /* VIS_CLD_250 */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = 0;
     /* APPL_HCLD_138 */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = 0;
     /* APPL_VIS_REFL */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = 0;
     /* APPL_VIS_NIR_RATIO */
      uoff2 += npix * nlin;
      (*ubdat)[ uoff2 + ipx + npix * iln ] = 0;
     /* end all inputs */
      
if( ( ipx == tpix ) && ( iln == tlin ) && ( dbg_print == 1 ) )
  {
  printf( "/n%s, %d: mix, t, p, h profiles AFTER, pix %d, lin %d\n", __FILE__, __LINE__, ipx, iln );
  profile_print( (*tdat) + foff , (*ubdat), (*i32dat), ipx, iln, npix, nlin );
  }
     /* 
      *  note that we also get the trop_lvl, sfc_lvl which go into another 
      *  storage area further on in the work
      */
      }
    }
 /*
  *  add the cloud band presence array to the end of the byte array
  */
  memcpy( ( *ubdat + 25 * npix * nlin ), cmp_there, ncmp_bnd );

  return(0);
  }

int profile_print( float *prof, unsigned char *ubdat, int32_t *i32dat, 
  int32_t ipx, int32_t iln, int32_t npix, int32_t nlin )
 /*******************************************************************

   profile_print

   purpose: print the met profile and 2d data at a point

   Returns type: int - 0 if good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float *           prof             I      pointer to the start of the 
                                                met data
      unsigned char *   ubdat            I      unsigned char data
      int32_t           ipx              I      pixel to look at
      int32_t           iln              I      line to look at
      int32_t           npix             I      # pixels
      int32_t           nlin             I      # lines
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       28-Jan-2019     Original development

*******************************************************************/
  {
  int nlvl_cmp = 101;
  int32_t off, ilvl;
  /* just print each profile */
  printf( "mixr profile\n" );
  for( ilvl = 0; ilvl < nlvl_cmp; ilvl++ )
    {
    printf( "%f ", prof[ ipx + npix * ( iln + nlin * ilvl ) ] );
    }
  printf( "\nT profile\n" );
  off = npix * nlin * nlvl_cmp;
  for( ilvl = 0; ilvl < nlvl_cmp; ilvl++ )
    {
    printf( "%f ", prof[ off + ipx + npix * ( iln + nlin * ilvl ) ] );
    }
  printf( "\nPRESS profile\n" );
  off += npix * nlin * nlvl_cmp;
  for( ilvl = 0; ilvl < nlvl_cmp; ilvl++ )
    {
    printf( "%f ", prof[ off + ipx + npix * ( iln + nlin * ilvl ) ] );
    }
  printf( "\nHEIGHT profile\n" );
  off += npix * nlin * nlvl_cmp;
  for( ilvl = 0; ilvl < nlvl_cmp; ilvl++ )
    {
    printf( "%f ", prof[ off + ipx + npix * ( iln + nlin * ilvl ) ] );
    }
  printf( "\n" );
 /* on to the single-level params */
  off += npix * nlin * nlvl_cmp;
  printf( "TSFC: %f, ", prof[ off + ipx + npix * iln] );
  off += npix * nlin;
  printf( "PSFC: %f, ", prof[ off + ipx + npix * iln] );
  off += npix * nlin;
  printf( "WIND SP: %f, ", prof[ off + ipx + npix * iln] );
  off += npix * nlin;
  printf( "O3: %f\n", prof[ off + ipx + npix * iln] );

  off += npix * nlin;
  printf( "ICE FR:  %f, ", prof[ off + ipx + npix * iln] );
  off += npix * nlin;
  printf( "SNO FR:  %f, ", prof[ off + ipx + npix * iln] );
  off += npix * nlin;
  printf( "ALT O3:  %f, ", prof[ off + ipx + npix * iln] );
  off += npix * nlin;
  printf( "ALT ICEC:  %f\n", prof[ off + ipx + npix * iln] );
 /*  the first 2 met byte quantities - the 101 level profile surface 
     and tropopause levels */
  off = 0;
  printf( "SFC_LVL: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin; 
  printf( "TROP_LVL: %d\n", ubdat[ off + ipx + npix * iln ] );
 /*  the 'alternate snow/ice type (from NSIDC) */
  printf( "ALT ICE: %d\n", i32dat[ ipx + npix * iln ] );
 /*  the remaining byte values are the cloud mask stuff */
  off = 2 * npix * nlin;
  printf( "---- Cloud Mask values ----\n" );
  printf( "CLD_DET: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "CONF_CLD: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "CLR_66: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "CLR_95: %d\n", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "CLR_99: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "SNO_SFC: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "WTR_SFC: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "COAST_SFC: %d\n", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "DESERT_SFC: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "LND_SFC: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "NIGHT: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "GLINT: %d\n", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "OCEAN_NO_SNOW: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "OCEAN_SNOW: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "LND_NO_SNOW: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "LND_SNOW: %d\n", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "TST_H_138: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "TST_VIS_REFL: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "TST_VIS_RATIO: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "VIS_CLD_250: %d\n", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "APPL_HCLD_138: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "APPL_VIS_REFL: %d, ", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
  printf( "APPL_VIS_NIR_RATIO: %d\n", ubdat[ off + ipx + npix * iln ] );
  off += npix * nlin;
 /* end? */
  return 0;
  }
/* WDR **** test end-around to mk_cmp_prof */
int mk_cmp_prof_dum( float *merra_p, float *merra_t, float *merra_q,
  float *merra_h, float sfcp, float sfct, float sfcrh, float lat,
  unsigned char *sfc_lvl, unsigned char *trop_lvl, double *cmp_p,
  double *cmp_t, double *cmp_mixr, double *cmp_h )
  {
int i;
double cmp_t_sav[] = {169.25176566420146, 197.35942572651811, 229.18349451149692, 256.22238210609601, 236.47864107966123, 244.04297143175555, 250.64605396062848, 255.49783596464286, 259.17682096634508, 260.50430632229126, 257.86788008325004, 255.09753268656993, 251.71866532458688, 246.82264824496451, 241.36301970269031, 235.5349625222656, 231.00564672010188, 230.31934481415178, 230.00717490122784, 230.93984650147675, 231.87556038742846, 230.052070835427, 226.94875351064999, 223.97554681909486, 221.12246382681033, 218.38062746953614, 216.05483908715939, 214.30433557126059, 212.6156777870857, 210.98488863951465, 209.57034993410883, 208.26534656680437, 207.00140487891071, 205.76477951859837, 204.56229701519121, 202.84595260841155, 200.1528837094385, 197.5345872033561, 194.98734605240188, 193.13906753041528, 192.9482175288388, 192.76222219178675, 192.58086243374066, 192.4039337826197, 193.53690804615815, 196.34328186967463, 199.08486131023724, 201.76426680177968, 204.38396142841776, 206.94626326403974, 209.4240557622999, 211.67158053966773, 213.87221009140157, 216.02762928535265, 218.13943261959307, 220.23999924343335, 222.60373354556071, 224.92166874653532, 227.195309078017, 229.42608466977327, 232.27719286025973, 235.2490670940424, 238.1670467211311, 241.03276062011722, 243.91780941249539, 246.75240714620259, 249.53800835678237, 251.86502860055904, 253.8587763082154, 255.81919333430849, 257.63971457966039, 259.29891772708743, 260.9312384604865, 262.72494342240276, 264.64753610551918, 266.53988144929855, 268.11359331445857, 269.60080666913444, 271.01286592906087, 272.25331134687139, 273.4751007120372, 274.88788770642338, 276.32402654921054, 277.67617255366838, 278.94743905007527, 280.29684643297759, 281.87085668556011, 283.47072360439103, 285.24727818913414, 286.9975584364218, 288.51282458375113, 289.86087429979796, 290.9986136887851, 292.01521247752953, 292.98262145756013, 294.3041990659554, 296.79569609374886, 298.98862037321277, 299.76233992810609, 300.50664057118684, 301.24118800223761};
double cmp_mixr_sav[] = {0.0030000000000000001, 0.0030000000000000001, 0.0030000000000000001, 0.0030000000000000001, 0.0039293894442529974, 0.0041105013907719468, 0.0042181900744832808, 0.0042143496342834938, 0.0041817655664085482, 0.0040206101090032517, 0.0038733849643674614, 0.0037377625227534746, 0.0036014648248259636, 0.0034516475997961819, 0.0033167026788735537, 0.0031947714963073749, 0.0030799596646610928, 0.0030107622372746654, 0.0029488582454394335, 0.0028982059410255015, 0.0028502657213637195, 0.0027944201318237261, 0.002735964452701974, 0.0026799596042977963, 0.0026262174675935851, 0.0025745708296473057, 0.0025405554258138494, 0.0025323804296621087, 0.0025244942581909551, 0.0025168783380373483, 0.0025258323028017112, 0.0025409217199495382, 0.0025555363514953174, 0.0025643709192777427, 0.0025716416839911219, 0.0025546327064356832, 0.0024931393172969081, 0.002433353276920958, 0.0023751897106692147, 0.0023510380165387293, 0.0024102799348228623, 0.0024680149123855986, 0.0025243109557570816, 0.0025792315352972882, 0.0028074932565275365, 0.0032578001900602956, 0.0036977103048750353, 0.004127644072243016, 0.004547996715297346, 0.0049591401890269129, 0.0066045793017774166, 0.015756209535939904, 0.024716888476006896, 0.033493476685974542, 0.042092466755785318, 0.053504463790107065, 0.094126872282568774, 0.13396219155187214, 0.17303627285925954, 0.21137369397371381, 0.2419968812989492, 0.27022626866094007, 0.2979437192418275, 0.32516470551490928, 0.47386770138225009, 0.6199703171142984, 0.76354752329508391, 1.0340410272474989, 1.3924930017834096, 1.7449525187114237, 1.9710464191297217, 2.0459308727278538, 2.1196020481996749, 2.5849697692954705, 3.3725913913975889, 4.1478216949005233, 4.7278626505249006, 5.2595586368380491, 5.6780053593354278, 5.7887852404941018, 5.897899014248777, 6.0126629733222074, 6.1272686176990323, 6.2612715902934664, 6.4138250668383652, 6.5300489624532512, 6.6239450570569449, 7.0237148492647679, 7.815759173341303, 8.9338062006928265, 10.017175527482264, 10.798574208009653, 11.482688439271284, 12.855657707404751, 14.06097496976534, 15.451677990360071, 8.9087798522868766, 0.078892832367705168, 0.0013115547243546515, 0.0021613451297737749, 0.0030000000000000001};
double cmp_h_sav[] = {80980.567272357701, 75343.876190403083, 70225.982411112418, 65480.628319062816, 61428.924238065607, 58024.80202191414, 54950.191567802227, 52150.37234583289, 49587.069685576469, 47234.584895480606, 45084.856633888245, 43122.927977172032, 41324.952076504, 39676.318479632821, 38164.924587601548, 36777.388481093178, 35497.403939587413, 34300.298902152295, 33167.424793334678, 32088.91680780418, 31057.094888806638, 30073.813431125316, 29143.272510847557, 28263.483356793786, 27430.050001579777, 26639.098949354469, 25886.659507298184, 25168.400851046954, 24481.003799506794, 23822.279970437601, 23190.014293749664, 22582.119678830768, 21996.942353189086, 21433.083441024752, 20889.269993475449, 20365.028820377422, 19861.183794284312, 19377.767922145697, 18913.565376345625, 18466.731479797825, 18033.779798131261, 17612.242600176087, 17201.593473696255, 16801.341547642598, 16409.703419308389, 16023.364531449324, 15640.568137969527, 15261.317752839021, 14885.61129854926, 14513.441856962587, 14144.824125415456, 13779.920194073944, 13418.849581909997, 13061.570770637654, 12708.041363820343, 12358.193358432209, 12011.718940291992, 11668.35429065583, 11328.08453068694, 10990.892967344218, 10656.282389846725, 10323.666245846885, 9992.9598216377744, 9664.201332244118, 9337.3653911268993, 9012.4264031402581, 8689.4162094312233, 8368.613449980001, 8050.4524006411266, 7735.0947019042451, 7422.5961680835417, 7113.1000037295735, 6806.6644285353668, 6503.107856676912, 6202.1379607025819, 5903.6184955729223, 5607.7124190540144, 5314.6091843005188, 5024.3518657619488, 4737.0487534217273, 4452.7644388850267, 4171.3480686550502, 3892.6336221371671, 3616.5978083090913, 3343.2676087364471, 3072.5944718766254, 2804.3890955354977, 2538.46734163782, 2274.6417716957349, 2012.7539832113782, 1752.8699864080911, 1495.1635489650282, 1239.7954748266484, 986.82637292768879, 736.25029352497961, 487.89682088402327, 241.65371459693083, 0, -243.21093547707787, -481.72170087446557, 0};
*sfc_lvl = 97;
*trop_lvl = 44;
for( i = 0; i < 101; i++ )
  {
  *(cmp_t + i) = *(cmp_t_sav + i);
  *(cmp_mixr + i) = *(cmp_mixr_sav + i);
  *(cmp_h + i) = *(cmp_h_sav + i);
  }
  return 0;
  }

int mk_cmp_prof( float *merra_p, float *merra_t, float *merra_q, 
  float *merra_h, float sfcp, float sfct, float sfcrh, float lat, 
  unsigned char *sfc_lvl, unsigned char *trop_lvl, double *cmp_p, 
  double *cmp_t, double *cmp_mixr, double *cmp_h )
  {
  unsigned char get_trop( double *, double *, unsigned char );
  int height_profile_( double *, double *, double *, double *, int *, 
    double * );
  int profile_to_101_( double *, double *, double *, int *, double *,
                   double *, double *, double *, int * );
  int32_t nlvl = 42, nlvl_cmp = 101, ioz = 0, ilvl, kstart;
  int32_t nlvl_d = nlvl + 1;
  float oval, a_factor;
  double p_init[nlvl_d], lcl_p[nlvl_d], lcl_t[nlvl_d];
  double lcl_w[nlvl_d], lcl_h[nlvl_d], lat_d, sfcp_d;
 /*
  *  Some l1 conventions that need to be addressed:
  *  the p, q, mixr, h arrays need a extra bottom level at 1100 mb (beyond
  *  normal seen p) and need to be double.
  *  the merra_q needs to be converted to mixr via:
  *  w = q ( 1 - q )^-1  q = spec hum, w = mixing ratio
  * and the merra_t needs 273.15 added to it to get deg K
  * and the sfc RH needs to be made into a mixing ratio
  * lastly, for levels with (-) t, ppropagate higher level to it
  */
  lat_d = lat;
  sfcp_d = sfcp;

  for( ilvl = 0; ilvl < nlvl; ilvl++ )
    {
    merra_q[ilvl] = 1000. * merra_q[ilvl] / ( 1. - merra_q[ilvl]);
    merra_t[ilvl] += 273.15;
    }
  sfct +=  273.15;
  /*  use cvt_rh_to_q  */
  met_cvt_rh_to_q( 1, &sfcp, MET_UNITS__P_HPA, &sfct, MET_UNITS__T_K,
    &sfcrh, &oval, MET_UNITS__Q_G_KG );
/*  WDR test to see that the conversion worked 
float ex_val;
met_cvt_q_to_rh(1, &sfcp, MET_UNITS__P_HPA, &sfct, MET_UNITS__T_K, 
   &oval, MET_UNITS__Q_KG_KG, &ex_val );
  END TEST */
  sfcrh = oval  / ( 1. - oval );  /* now, this is the mixing ratio */

 /* transfer the 42 level profile info to the local double 43 level arrays */
  for( ilvl = 0; ilvl < nlvl; ilvl++ )
    {
    lcl_p[ilvl] = merra_p[ilvl];
    lcl_t[ilvl] = merra_t[ilvl];
    lcl_w[ilvl] = merra_q[ilvl];
    lcl_h[ilvl] = merra_h[ilvl];

    if( lcl_t[ilvl] < 0. )
      {
      lcl_t[ilvl] = lcl_t[ ilvl - 1 ];
      lcl_w[ilvl] = lcl_w[ ilvl - 1 ];
      lcl_h[ilvl] = lcl_h[ ilvl - 1 ];
      }
    }
  lcl_p[nlvl] = 1100.;
  lcl_t[nlvl] = lcl_t[ nlvl - 1 ];
  lcl_w[nlvl] = lcl_w[ nlvl - 1 ];
  lcl_h[nlvl] = lcl_h[ nlvl - 1 ];

/*  we will just duplicate the method used in ancillary_module.f90
    to set up the incoming profile and placing the surface values 
    and then interpolation even though I'd probably do it differently
*/
  memcpy( p_init, lcl_p, nlvl_d * sizeof(double) );

 /* first, if sfc p is > next lowest p level, set the p level as the last */
  if( ( sfcp > 0 ) && ( lcl_p[ nlvl_d - 2 ] > sfcp ) )
    *sfc_lvl = nlvl_d - 1;
  else
    *sfc_lvl = 0;

 /* look through the lower half of the levels and find the sfc_lvl for 
    the rest */
  kstart = nlvl_d / 2;

  for( ilvl = kstart; ilvl < nlvl_d; ilvl++ )
    {
    if( ( sfcp > 0 ) && ( lcl_p[ilvl] > sfcp ) )
      {
      if( *sfc_lvl == 0 )
        {
        *sfc_lvl = ilvl;
        lcl_t[ilvl] = sfct;
       /* now this next part is ver-batum from chimaera code but I don't 
          understand the reason they did it this way */
        if( ( sfcp - lcl_p[ ilvl - 1 ] < 5. ) || 
            ( lcl_p[ilvl] - sfcp < 5. ) )
          {
          lcl_p[ilvl] = ( lcl_p[ilvl] + lcl_p[ ilvl - 1 ] ) / 2.;
          }
        else
          {
          lcl_p[ilvl] = sfcp;
          }
        lcl_w[ilvl] = sfcrh;
        }
      else
        {
        lcl_t[ilvl] = sfct;
        lcl_w[ilvl] = sfcrh;
        lcl_p[ilvl] = p_init[ ilvl - 1 ];
        }
      }
    }
 /* more dubvious mods: add the surface level into the coarse profile */
  if( *sfc_lvl != nlvl_d )
    merra_p[ nlvl_d - 1 ] =  p_init[ ilvl - 1 ];
  else
    lcl_p[ nlvl - 1 ] = sfcp;
  lcl_t[ nlvl - 1 ] = sfct;
  lcl_w[ nlvl - 1 ] = sfcrh;

 /* on to getting the lowest level of the 101 level profile */
  kstart = nlvl_cmp / 2;
  for( ilvl = kstart; ilvl < nlvl_cmp; ilvl++ )
    {
    if( cmp_p[ilvl] >= sfcp )
      {
      *sfc_lvl = ilvl;
      break;
      }
    }

  profile_to_101_( lcl_p, lcl_t, lcl_w, &nlvl, &lat_d, cmp_p, cmp_t, 
    cmp_mixr, &ioz );

 /*  instead of assigning the sfct, sfcrh (now mixr) to the sfc_lvl, it
     interpolates in p coordinates the resulting T, mixr */
  a_factor = ( sfcp - cmp_p[ *sfc_lvl - 1 ] ) / 
             ( cmp_p[*sfc_lvl] - cmp_p[ *sfc_lvl - 1 ] );
  cmp_t[*sfc_lvl] = cmp_t[ *sfc_lvl - 1 ] + a_factor * 
    ( cmp_t[*sfc_lvl] - cmp_t[ *sfc_lvl - 1 ] );
  cmp_mixr[*sfc_lvl] = cmp_mixr[ *sfc_lvl - 1 ] + a_factor *
    ( cmp_mixr[*sfc_lvl] - cmp_mixr[ *sfc_lvl - 1 ] );

  cmp_p[*sfc_lvl] = sfcp;

 /* and then follow with the height profile */
  height_profile_( cmp_p, cmp_t, cmp_mixr, cmp_h, &nlvl_cmp, &sfcp_d );
  cmp_h[nlvl_cmp - 1 ] = 0.;

 /* also find the tropopause level - call a chimaera-based algorithm */
  *trop_lvl = get_trop( cmp_t, cmp_p, *sfc_lvl );
  
  return 0;
  }

unsigned char get_trop( double *cmp_t, double *cmp_p, unsigned char sfc_lvl )
  {
  unsigned char trop_lvl;
  int32_t ilev, imin;
  float xmin, ptop = 100.;
  
  xmin = 999999.;
  imin = 1;

  for( ilev = 0; ilev < sfc_lvl - 5; ilev++ )
    {
    if( ( cmp_t[ilev] < xmin ) && ( cmp_p[ilev] > ptop ) )
      {
      xmin = cmp_t[ilev];
      imin = ilev;
      }
    }

 /* don't allow trop height > 400 mb (level 71) */
  trop_lvl = 200;
  for( ilev = imin; ilev < 71; ilev++ )
    {
    if( ( cmp_t[ ilev - 1 ] >= cmp_t[ilev] ) && 
        ( cmp_t[ ilev + 1 ] > cmp_t[ilev] ) )
      {
      trop_lvl = ilev;
      break;
      }
    }

  if( trop_lvl == 200 ) trop_lvl = imin;

  return trop_lvl;
  }

int ncio_grab_ub_ds(int ncid, char *ds_name, unsigned char *data)
/*******************************************************************

   ncio_grab_ub_ds

   purpose: grab dataset and return it in unsigned char format

   Returns type: int -  0 if OK, -1 if problem

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               ncid             I      netcdf id of file
      char *            ds_name          I      name of dataset to read
      unsigned char *   data             O      pointer to a pre-allocated
                                                array to receive the data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       26 Nov 2018     original development

 *******************************************************************/ {
    int status;
    int var_id;
    /*
     *  get ID of the dataset
     */
    if ((status = nc_inq_varid(ncid, ds_name, &var_id)) != NC_NOERR) {
        printf("%s, %d: nc_inq_varid returned error %d\n", __FILE__, __LINE__,
                status);
        return -1;
    }
    /*
     *  read the dataset as a unsigned char
     */
    if ((status = nc_get_var_uchar(ncid, var_id, data)) != NC_NOERR) {
        printf("%s, %d: nc_get_var_uchar returned error %d\n", 
        __FILE__, __LINE__, status);
        return -1;
    }
    return 0;
}

int ncio_grab_i32_ds(int ncid, char *ds_name, int32_t *data)
/*******************************************************************

   ncio_grab_i32_ds

   purpose: grab dataset and return it in nteger format

   Returns type: int -  0 if OK, -1 if problem

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               ncid             I      netcdf id of file
      char *            ds_name          I      name of dataset to read
      int32_t *         data             O      pointer to a pre-allocated
                                                array to receive the data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       26 Nov 2018     original development

 *******************************************************************/ {
    int status;
    int var_id;
    /*
     *  get ID of the dataset
     */
    if ((status = nc_inq_varid(ncid, ds_name, &var_id)) != NC_NOERR) {
        printf("%s, %d: nc_inq_varid returned error %d\n", __FILE__, __LINE__,
                status);
        return -1;
    }
    /*
     *  read the dataset as a unsigned char
     */
    if ((status = nc_get_var_int(ncid, var_id, (int32_t *)data)) != NC_NOERR) {
        printf("%s, %d: nc_get_var_int returned error %d\n", 
        __FILE__, __LINE__, status);
        return -1;
    }
    return 0;
}
