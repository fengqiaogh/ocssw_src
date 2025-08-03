#include <stdlib.h>
#include <stdio.h>
#include "l12_proto.h"
#include "l2prod.h"
#include "l1.h"
#include  "get_ctht.h"
#include "sensorInfo.h"
#define C_TO_K 273.15

ctht_lins_str ctht_lins;
/* 
    file get_ctht.c contains most of the cloud top height source code
*/
float *get_ctht_lin( l2str *l2rec, int prodnum )
/*
get_ctht_lin - get the cloud top height info for a current scan
   Returns - float * - the record of the scan's product data
   
   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     l2str *           l2rec            I      l2 record of data
     int               prodnum          I      Product ID to return
                                               Note that if prodnum = -1
                                               we won't return any product
                                               This is when the code is called 
                                               to prepare for getting cmp 
                                               cloud values

   W. Robinson, SAIC, 22 Jun 2023
*/
  { 
  static int32_t cur_ctht_scan = -1;
                             // being managed: 0 - not filled, 1 filled
                             // (may need status for '-1 line or nlin + 1)
  static int32_t firstin = 0;
  static float *outbuf;
  static int32_t *outint;   //  so ints will get treated right
  int32_t cur_scan = l2rec->l1rec->iscan;
  int32_t npix = l2rec->l1rec->npix;
  int32_t clin, ipix, ntyp = 2, phase_v;
  ctht_lins.nrec = 3;

  outint = (int32_t *)outbuf;

  if( firstin == 0 )
    {
    firstin = 1;
    if( ( outbuf = malloc( npix * ntyp * sizeof( float ) ) ) == NULL )
      {
      printf( "%s, %d: E - Unable to allocate the ctht output buffer\n",
        __FILE__, __LINE__ );
      exit(1);
      }
    }

  clin = ctht_lins.nrec / 2;
  if( cur_ctht_scan != cur_scan )
    {
    cur_ctht_scan = cur_scan;
    if( comp_ctht( l2rec ) != 0 )
      {
      printf( "%s, %d: E - CTH computation failure\n", __FILE__, __LINE__ );
      exit(1);
      }
    }
 /* output the specified data */
   // NOTE this should be in the center line of the 'queue'
  switch( prodnum )
    {
    case CAT_Cld_p:
      memcpy( outbuf, ctht_lins.ct[clin]->ctp, npix * sizeof( float ) );
      return outbuf;
      break;
    case CAT_Cld_t:
      memcpy( outbuf, ctht_lins.ct[clin]->ctt, npix * sizeof( float ) );
      return outbuf;
      break;
    case CAT_Cld_h:
      memcpy( outbuf, ctht_lins.ct[clin]->cth, npix * sizeof( float ) );
      return outbuf;
      break;
    case CAT_cth_cod:
      for( ipix = 0; ipix < npix; ipix++ )
        {
        phase_v = ctht_lins.ct[clin]->ct_phase[ipix];
        outbuf[ipix] = ( phase_v == BAD_UBYTE ) ? BAD_FLT : 
          ctht_lins.ct[clin]->cth_cod_all[ phase_v + ntyp * ipix ];
        }
      return outbuf;
      break; 
    case CAT_cth_cod_all:
      memcpy( outbuf, ctht_lins.ct[clin]->cth_cod_all,
        npix  * ntyp * sizeof( float ) );
      return outbuf;
      break;
    case CAT_cth_alb_all:
      memcpy( outbuf, ctht_lins.ct[clin]->cth_alb_all,
        npix * ntyp * sizeof( float ) );
      return outbuf;
      break;
    case CAT_cth_alb:
      for( ipix = 0; ipix < npix; ipix++ )
        {
        phase_v = ctht_lins.ct[clin]->ct_phase[ipix];
        outbuf[ipix] = ( phase_v == BAD_UBYTE ) ? BAD_FLT :
          ctht_lins.ct[clin]->cth_alb_all[ phase_v + ntyp * ipix ];
        }
      return outbuf;
      break;
    case CAT_cth_dcod:
      for( ipix = 0; ipix < npix; ipix++ )
        {
        phase_v = ctht_lins.ct[clin]->ct_phase[ipix];
        outbuf[ipix] = ( phase_v == BAD_UBYTE ) ? BAD_FLT :
          ctht_lins.ct[clin]->dcod[ phase_v + ntyp * ipix ];
        }
      return outbuf;
      break;
    case CAT_cth_dlcod:
      for( ipix = 0; ipix < npix; ipix++ )
        {
        phase_v = ctht_lins.ct[clin]->ct_phase[ipix];
        outbuf[ipix] = ( phase_v == BAD_UBYTE ) ? BAD_FLT :
          ctht_lins.ct[clin]->dlcod[ phase_v + ntyp * ipix ];
        }
      return outbuf;
      break;
    case CAT_cth_dcth:
      for( ipix = 0; ipix < npix; ipix++ )
        {
        phase_v = ctht_lins.ct[clin]->ct_phase[ipix];
        outbuf[ipix] = ( phase_v == BAD_UBYTE ) ? BAD_FLT :
          ctht_lins.ct[clin]->dcth[ phase_v + ntyp * ipix ];
        }
      return outbuf;
      break;
    case CAT_cth_dctp:
      for( ipix = 0; ipix < npix; ipix++ )
        {
        phase_v = ctht_lins.ct[clin]->ct_phase[ipix];
        outbuf[ipix] = ( phase_v == BAD_UBYTE ) ? BAD_FLT :
          ctht_lins.ct[clin]->dctp[ phase_v + ntyp * ipix ];
        }
      return outbuf;
      break;
    case CAT_cth_dctt:
      for( ipix = 0; ipix < npix; ipix++ )
        {
        phase_v = ctht_lins.ct[clin]->ct_phase[ipix];
        outbuf[ipix] = ( phase_v == BAD_UBYTE ) ? BAD_FLT :
          ctht_lins.ct[clin]->dctt[ phase_v + ntyp * ipix ];
        }
      return outbuf;
      break;
    case CAT_cth_dalb:
      for( ipix = 0; ipix < npix; ipix++ )
        {
        phase_v = ctht_lins.ct[clin]->ct_phase[ipix];
        outbuf[ipix] = ( phase_v == BAD_UBYTE ) ? BAD_FLT :
          ctht_lins.ct[clin]->dalb[ phase_v + ntyp * ipix ];
        }
      return outbuf;
      break;
    case CAT_cth_phase:
      memcpy( outbuf, ctht_lins.ct[clin]->ct_phase, 
        npix  * sizeof( unsigned char ) );
      return outbuf;
      break;
    case CAT_cth_lcod_all:
      memcpy( outbuf, ctht_lins.ct[clin]->cth_lcod_all,
        npix * ntyp * sizeof( float ) );
      return outbuf;
      break;
    case CAT_cth_cth_all:
      memcpy( outbuf, ctht_lins.ct[clin]->cth_all,
        npix * ntyp * sizeof( float ) );
      return outbuf;
      break;
    case CAT_cth_ctp_all:
      memcpy( outbuf, ctht_lins.ct[clin]->ctp_all,
        npix * ntyp * sizeof( float ) );
      return outbuf;
      break;
    case CAT_cth_cost_all:
      memcpy( outbuf, ctht_lins.ct[clin]->cost_ss,
        npix * ntyp * sizeof( float ) );
      return outbuf;
      break;
    case CAT_cth_acost_all:
      memcpy( outbuf, ctht_lins.ct[clin]->acost,
        npix * ntyp * sizeof( float ) );
      return outbuf;
      break;
    case CAT_cth_iter_all:
      memcpy( outint, ctht_lins.ct[clin]->nitr,
        npix * ntyp * sizeof( int32_t) );
      return outbuf;
      break;
    case CAT_cth_cth_raw_all:
      memcpy( outbuf, ctht_lins.ct[clin]->cth_raw,
        npix * ntyp * sizeof( float ) );
      return outbuf;
      break;
    case CAT_cth_ctt_all:
      memcpy( outbuf, ctht_lins.ct[clin]->ctt_all,
        npix * ntyp * sizeof( float ) );
      return outbuf;
      break;
    //  back to all 2-d arrays that get extracted from the 3-d
    case CAT_cth_lcod:
      for( ipix = 0; ipix < npix; ipix++ )
        {
        phase_v = ctht_lins.ct[clin]->ct_phase[ipix];
        outbuf[ipix] = ( phase_v == BAD_UBYTE ) ? BAD_FLT :
          ctht_lins.ct[clin]->cth_lcod_all[ phase_v + ntyp * ipix ];
        }
      return outbuf;
      break;
    case CAT_cth_cost:
      for( ipix = 0; ipix < npix; ipix++ )
        {
        phase_v = ctht_lins.ct[clin]->ct_phase[ipix];
        outbuf[ipix] = ( phase_v == BAD_UBYTE ) ? BAD_FLT :
          ctht_lins.ct[clin]->cost_ss[ phase_v + ntyp * ipix ];
        }
      return outbuf;
      break;
    case CAT_cth_acost:
      for( ipix = 0; ipix < npix; ipix++ )
        {
        phase_v = ctht_lins.ct[clin]->ct_phase[ipix];
        outbuf[ipix] = ( phase_v == BAD_UBYTE ) ? BAD_FLT :
          ctht_lins.ct[clin]->acost[ phase_v + ntyp * ipix ];
        }
      return outbuf;
      break;
    case CAT_cth_iter:
      for( ipix = 0; ipix < npix; ipix++ )
        {
        phase_v = ctht_lins.ct[clin]->ct_phase[ipix];
        outint[ipix] = ( phase_v == BAD_UBYTE ) ? BAD_INT :
          ctht_lins.ct[clin]->nitr[ phase_v + ntyp * ipix ];
        }
      return outbuf;
      break;
    case CAT_cth_cth_raw:
      for( ipix = 0; ipix < npix; ipix++ )
        {
        phase_v = ctht_lins.ct[clin]->ct_phase[ipix];
        outbuf[ipix] = ( phase_v == BAD_UBYTE ) ? BAD_FLT :
          ctht_lins.ct[clin]->cth_raw[ phase_v + ntyp * ipix ];
        }
      return outbuf;
      break;
    case CAT_cth_dcod_all:
      memcpy( outbuf, ctht_lins.ct[clin]->dcod,
        npix * ntyp * sizeof( float ) );
      return outbuf;
      break;
    case CAT_cth_dlcod_all:
      memcpy( outbuf, ctht_lins.ct[clin]->dlcod,
        npix * ntyp * sizeof( float ) );
      return outbuf;
      break;
    case CAT_cth_dcth_all:
      memcpy( outbuf, ctht_lins.ct[clin]->dcth,
        npix * ntyp * sizeof( float ) );
      return outbuf;
      break;
    case CAT_cth_dctp_all:
      memcpy( outbuf, ctht_lins.ct[clin]->dctp,
        npix * ntyp * sizeof( float ) );
      return outbuf;
      break;
    case CAT_cth_dctt_all:
      memcpy( outbuf, ctht_lins.ct[clin]->dctt,
       npix * ntyp * sizeof( float ) );
      return outbuf;
      break;
    case CAT_cth_dalb_all:
      memcpy( outbuf, ctht_lins.ct[clin]->dalb,
       npix * ntyp * sizeof( float ) );
      return outbuf;
      break;
    case -1:
      return NULL;
    }
 /*  if this is reached, we need another CAT */
  printf( "%s, %d: E - CTH computation could not find the catalog # %d\n", 
    __FILE__, __LINE__, prodnum );
  exit(1);
  }

int32_t comp_ctht( l2str *l2rec )
 /*
  comp_ctht will run the cloud top height algorithm on the lines of data in 
  the l1que that have not been processed yet and keep a local copy of the lines

  Returns - int32_t 0 if all  good, no error conditions yet

  Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     l2str *           l2rec            I      l2 record of data

   W. Robinson, SAIC, 22 Jun 2023
  */
  {
  int32_t ctht_nrec = 3;  // # lines in our ctht records, matches those 
                            // in the cmp, see get_cmp.c
  static int32_t *ctht_stat;
  static int32_t ctht_init = 0;
  int32_t ilin, il2, lin_in_que, scn_in_que, ntyp = 2;
  extern l1qstr l1que;
  static ctht_lin_str *ctht_sav;
  int32_t npix = l2rec->l1rec->npix;

  if( input->proc_cloud != 1 ) {
    printf( "E, %s, %d: For cloud processing, proc_cloud must be set to 1\n",
      __FILE__, __LINE__ );
    exit(1);
  }

  if( l2rec->l1rec->anc_add == NULL ) {
    printf( "E, %s, %d: For cloud processing, anc_profile_1, _2, _3 must be specified\n", __FILE__, __LINE__ );
    exit(1);
  }

  //  until climatology sfcp, sfct is cleared up, we can't use climatology
  if( strstr( input->met1, "climatology" ) != NULL ) {
    printf( "E, %s, %d: Ancillary meteorological data from a climatology file cannot be used for cloud processing\n", __FILE__, __LINE__ );
    printf( "Use daily files in MET1, 2, 3 l2gen inputs\n" );
    exit(1);
  }

 /*  initialize the ctht line storage and lines */
  if( ctht_init == 0 )
    {
    ctht_init = 1;
   /* set up the ctht record status */
    if( ( ( ctht_lins.ct = (ctht_lin_str **)malloc( 
        ctht_nrec * sizeof( ctht_lin_str * ) ) ) == NULL ) ||
        ( ( ctht_stat = (int32_t *)malloc(
        ctht_nrec * sizeof( int32_t ) ) ) == NULL ) )
      {
      printf( "%s - %d: Allocation of ctht record space failed\n",
        __FILE__, __LINE__ );
      exit(1);
      }
    for( ilin = 0; ilin < ctht_nrec; ilin++ )
      {
      if( ( ctht_lins.ct[ilin] = 
        (ctht_lin_str *) malloc( sizeof( ctht_lin_str) ) ) == NULL )
        {
        printf( "%s - %d: Allocation of ctht record space failed\n",
          __FILE__, __LINE__ );
        exit(1);
        }
      ctht_lins.ct[ilin]->iscan = -1;
      ctht_lins.ct[ilin]->npix = npix;
      if( ( ( ctht_lins.ct[ilin]->cth = 
            (float *) malloc(npix * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->ctp = 
            (float *) malloc(npix * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->ctt = 
           (float *) malloc(npix * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->ct_phase =
            (unsigned char *) malloc(npix * sizeof(unsigned char) ) ) == NULL )
          || ( ( ctht_lins.ct[ilin]->cth_all =
           (float *) malloc(npix * ntyp * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->cth_raw =
           (float *) malloc(npix * ntyp * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->ctp_all = 
           (float *) malloc(npix * ntyp * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->ctt_all  =
           (float *) malloc(npix * ntyp * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->cth_lcod_all =    
           (float *) malloc(npix * ntyp * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->cth_alb_all =
           (float *) malloc(npix * ntyp * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->cost_ss =
           (float *) malloc(npix * ntyp * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->acost =
           (float *) malloc(npix * ntyp * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->nitr =
           (int32_t *) malloc(npix * ntyp * sizeof(int32_t) ) ) == NULL ) )
        {
        printf( "%s - %d: Allocation of ctht record space failed\n",
          __FILE__, __LINE__ );
        exit(1);
        }
      if( ( ( ctht_lins.ct[ilin]->cth_cod_all =
            (float *) malloc(npix * ntyp * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->oe_akdiag =
            (float *) malloc(npix * 3 * ntyp * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->dlcod =
            (float *) malloc(npix * ntyp * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->dcod =
            (float *) malloc(npix * ntyp * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->dcth =
            (float *) malloc(npix * ntyp * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->dctp =
            (float *) malloc(npix * ntyp * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->dctt =
            (float *) malloc(npix * ntyp * sizeof(float) ) ) == NULL ) ||
          ( ( ctht_lins.ct[ilin]->dalb =
            (float *) malloc(npix * ntyp * sizeof(float) ) ) == NULL ) )
        {
        printf( "%s - %d: Allocation of ctht record space failed\n",
          __FILE__, __LINE__ );
        exit(1); 
        }
      }
    }
 /*  OK, arrange the existant ctht records that we have and process any 
     records we don't have in comp_ctht_lin 
  */

  //  go thru ctht lines
  for( ilin = 0; ilin < ctht_nrec; ilin++ )
    {
    ctht_stat[ilin] = 0;
    }
  for( ilin = 0; ilin < ctht_nrec; ilin++ )
    {
    lin_in_que = ilin + l1que.nq / 2 - ctht_nrec / 2;
    // if lin_in_que < 0 we're in trouble...
    scn_in_que = l1que.r[lin_in_que].iscan;
    //  search for the l1rec match line in ctht
    for( il2 = 0; il2 < ctht_nrec; il2++ )
      {
      if( ctht_lins.ct[il2]->iscan == scn_in_que )
        {
        // switch addresses
        ctht_sav = ctht_lins.ct[ilin];
        ctht_lins.ct[ilin] = ctht_lins.ct[il2];
        ctht_lins.ct[il2] = ctht_sav;
        ctht_stat[ilin] = 1;
        break;
        }
      }
    }
  // for any lines not transferred, process ctht
  for( ilin = 0; ilin < ctht_nrec; ilin++ )
    {
    if( ctht_stat[ilin] == 0 )
      {
      lin_in_que = ilin + l1que.nq / 2 - ctht_nrec / 2;
      scn_in_que = l1que.r[lin_in_que].iscan;
      //if( comp_ctht_lin( l2rec->l1rec, ctht_lins.ct[ilin] ) != 0 )
      if( comp_ctht_lin( &(l1que.r[lin_in_que]), ctht_lins.ct[ilin] ) != 0 )
        {
        printf( "%s - %d: Problem found running comp_ctht_lin\n",
          __FILE__, __LINE__ );
        exit(1);
        }
      }
    }
  return 0;
  }

int32_t comp_ctht_lin( l1str *l1rec, ctht_lin_str *ctht_lin )
  /*
   this will process the cloud top height for a line of satellite data

   Returns -  int32_t  - 0 for all good , only status for now

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
   l1str *      l1rec    input l1 record to use
   ctht_lin_str *ctht_lin  the cloud top height record to fill

   W. Robinson, SAIC, 22 Jun 2023, based on A. Sayer's IDL code
   W. Robinson, SAIC, 3 May 2024, improve the merra 2 profile treatment
  */
  {
  int32_t n_ctht_bnd = 8, n_prd, ny;
  static int32_t ctht_init = 0;
  static int32_t l1_bnd_ix[8];
  static double *sy, *tmp_lut;
  static ctht_lut_str full_lut_wtr, full_lut_ice;
  static ctht_unc_str unc;
  static oe_info_str oe_info;
  static ctht_lut_str *lut;
  static float *lcl_prof_t, *lcl_prof_h, *lcl_p_set;
  static int32_t lcl_nlvl, lcl_lvl_ptr, fill_lvl;
  float solz, senz, relaz, windsp, glint_ref;
  int32_t ipix, ityp, npix, nbnd_ref, do_const;
  int32_t ntyp = 2, foff, stype;
  int32_t lut_dims[4];
  float err1, err2, ev1, ev2;
  double rhot[8], xa[3];
  float *sxax[4], p_log_lo, p_log_hi;
  int32_t snax[4];
  int64_t n_block;
  int32_t ibnd, ibnd2, icod, icth, ialb, icor, ib1, ib2, ilev;
  // old: 
  //int32_t ctht_wav[] = { 755, 757, 760, 762, 765, 767, 770, 772 };
  // new 14 May 2024
  int32_t ctht_wav[] = { 754, 757, 759, 762, 764, 767, 769, 772 };
  int32_t ctht_all_bnd = 0, close_ind;
  int64_t ar_off[16], tmp_off, tmp_off2, off_cost;
  float ar_wt[16], orig_wt[4];
  float dp, found_hgt, lo_cost, spress_out;
  int32_t  bad_lt;
  // for gsl vector, matrix work
  double min_cost, yres[8], xres[8];
  double fg_cost_t1, fg_cost_t2;
  double int_pls, int_neg, int_ctr, coord, pctr;
  int32_t min_cost_loc[3];
  static double *fg_cost;
  static gsl_matrix *gm_sa, *gm_sy, *gm_sy_inv, *gm_sa_inv;
  static gsl_vector *gv_yres, *gv_sy_inv_yres, *gv_xres, *gv_sa_inv_xres;
  //  The final retrieved data to hang on to, yes, maybe we'll need both
  //  (or all) phases NOTE that ret_<name> will have all the product for 
  //  phase 0 followed by phase 1
  //  phase 0 for ipix: ret_<name>[ ipix + npix * ityp ]
  static float *ret_lcod;
  //  MERRA2 p levels
  float merra_p_set[] = { 1000., 975., 950., 925., 900., 875., 850., 825., 
      800., 775., 750., 725., 700., 650., 600., 550., 500., 450., 400.,
      350., 300., 250., 200., 150., 100., 70., 50., 40., 30., 20., 10., 
      7., 5., 4., 3., 2., 1., 0.7, 0.5, 0.4, 0.3, 0.1 };
  int32_t nlvl_merra2 = 42, found_phase;

  //  set some errors
  unc.noise_err = 0.005;
  unc.cal_err = 0.02;
  float lut_dalb = 0.1;  // Temporary till lut read in and albedo delta can be measured

  npix = l1rec->npix;
  do_const = 2;
  if( do_const != 2 )
    {
    comp_ctht_tst( l1rec, do_const, ctht_lin );
    }
  else
    {
    n_prd = 3;
    ny = 8;
   /*
    *  initialize the routine
    */
    if( ctht_init == 0 )
      {
      // get the location of ref_true (rhot) or bands to use for ctht
      ctht_init = 1;
      ctht_all_bnd = 1; // all bands there
      for( ibnd = 0; ibnd < n_ctht_bnd; ibnd++ )
        {
        l1_bnd_ix[ibnd] = bindex_get( ctht_wav[ibnd] );
        if( l1_bnd_ix[ibnd] < 0 )
          {
          ctht_all_bnd = 0;
          printf( "%s, %d, I: CTHT wave: %d not found\n", __FILE__, __LINE__,
            ctht_wav[ibnd] );
          }
        else
          {
          printf( "%s, %d, I: Got CTHT  band wavelength %d as index %d\n", 
            __FILE__, __LINE__, ctht_wav[ibnd], l1_bnd_ix[ibnd] ); 
          }
        }
      if( ctht_all_bnd == 0 )
        {
        printf( "%s, %d, E: Not all bands available for CTHT work\n",
          __FILE__, __LINE__ );
        exit(1);
        }
      // read in the cloud height tables uncertainty file parameters
      if( ctht_tbl_init( n_ctht_bnd, l1rec, &unc, &full_lut_wtr, &full_lut_ice ) != 0 )
        {
        printf( "%s, %d, E: Failure reading unc or lut\n", __FILE__, __LINE__ );
        exit(1);
        }
      //  check the prods and bands (ny) vs expected and 1 table against 
      //  the other
      if( ( full_lut_wtr.lut_nwave != full_lut_ice.lut_nwave ) ||
          ( full_lut_wtr.lut_nwave != ny ) )
        {
        printf( "%s, %d, E: expected # bands is not in the table\n",
          __FILE__,__LINE__ );
        exit(1);
        }
      if( ( full_lut_wtr.lut_nsolz != full_lut_ice.lut_nsolz ) ||
          ( full_lut_wtr.lut_nsenz != full_lut_ice.lut_nsenz ) ||
          ( full_lut_wtr.lut_nrelaz != full_lut_ice.lut_nrelaz ) ||
          ( full_lut_wtr.lut_nlev != full_lut_ice.lut_nlev ) ||
          ( full_lut_wtr.lut_npress != full_lut_ice.lut_npress ) ||
          ( full_lut_wtr.lut_nalb != full_lut_ice.lut_nalb ) ||
          ( full_lut_wtr.lut_ncod != full_lut_ice.lut_ncod ) ||
          ( full_lut_wtr.lut_naod != full_lut_ice.lut_naod ) ||
          ( full_lut_wtr.lut_ncth != full_lut_ice.lut_ncth ) )
        {
        printf( "%s, %d, E: LUTs (ice and water) have different array sizes\n",
          __FILE__,__LINE__ );
        exit(1);
        }
      //  set some arrays used below
      sy = (double *)malloc( n_ctht_bnd * n_ctht_bnd * sizeof(double) );

      if( ( fg_cost = (double *)
        malloc( full_lut_wtr.lut_ncod * full_lut_wtr.lut_ncth * 
        full_lut_wtr.lut_nalb * sizeof(double) ) ) == NULL )
        {
        printf( "%s, %d, E: unable to allocate the cost array\n",
          __FILE__, __LINE__ );
        exit(1);
        }
      //  set up the temp lut once
      lut = &full_lut_wtr;
      if( ( tmp_lut = (double *) malloc( lut->lut_nwave * lut->lut_ncod * 
        lut->lut_ncth * lut->lut_nalb * sizeof(double) ) ) == NULL )
        {
        printf( "%s - %d: E: Failure to allocate the pixel-specific LUT\n",
          __FILE__, __LINE__ );
        exit(1);
        }
      //
      gm_sa = gsl_matrix_alloc( n_prd, n_prd );
      gm_sy = gsl_matrix_alloc( ny, ny );
      gv_yres = gsl_vector_alloc( ny );
      gv_sy_inv_yres = gsl_vector_alloc( ny );
      gv_xres = gsl_vector_alloc(n_prd);
      gv_sa_inv_xres = gsl_vector_alloc(n_prd);

      //  set up the outputs
      //  OK we could do an if( all mallocs == NULL err
      if( ( ret_lcod = (float *) malloc( npix * ntyp * sizeof(float) ) ) == NULL )
        {
        printf( 
          "%s, %d, E: ctht was unable to allocate all the storage needed\n", 
          __FILE__, __LINE__ );
        exit(1);
        }
      // set the local profiles that are specific to each pixel
      if( ( ( lcl_prof_t = (float *) malloc( ( nlvl_merra2 + 1 ) * sizeof(float) ) ) == NULL )
       || ( ( lcl_prof_h = (float *) malloc( ( nlvl_merra2 + 1 ) * sizeof(float) ) ) == NULL )
       || ( ( lcl_p_set = (float *) malloc( ( nlvl_merra2 + 1 ) * sizeof(float) ) ) == NULL ) )
        {
        printf( "%s, %d, E: Unable to allocate local profile storage\n",
          __FILE__, __LINE__ );
        exit(1);
        }
      }  // end general preparation

    npix = l1rec->npix;
    nbnd_ref = l1rec->l1file->nbands;
   /*
    *  loop through all the pixels
    */
    for( ipix = 0; ipix < npix; ipix++ )
      {
     /*
      *  do only the cloudy pixels, all others will be BAD_FLT
      */
      init_ctht_parms( ctht_lin, ipix );
      //  replace below with true eval of cloudy condition
      if( l1rec->cloud[ipix] == 1 )
//printf( "%s,%d T: only testing pix 166\n", __FILE__, __LINE__ );
//if( ipix == 166 )  //  WDR for test of odd point
        {
        solz = l1rec->solz[ipix];
        senz = l1rec->senz[ipix];
        relaz = l1rec->delphi[ipix];
        windsp = l1rec->ws[ipix];

        //  looks like relaz can be outside 0:180, so fix it
        relaz = fmod( fabs( relaz ), 360. );
        relaz = ( relaz <= 180. ) ? relaz : 360. - relaz;
        // get toa reflectance for ctht bands
        bad_lt = 0;
        for( ibnd = 0; ibnd < n_ctht_bnd; ibnd++ )
          {
          foff = l1_bnd_ix[ibnd] + nbnd_ref * ipix;
          if( !( isfinite( l1rec->Lt[foff] ) ) || l1rec->Lt[foff] == BAD_FLT )
            bad_lt = 1;
          rhot[ibnd] = l1rec->Lt[foff] * OEL_PI / ( cos(solz *
            OEL_PI / 180. ) * l1rec->Fo[l1_bnd_ix[ibnd]] );
          }
        if( bad_lt ) break;  // just skip this pixel
        //  get the pixel type: 0 water, 1 land, 2 snow/ice
        stype = 0;  // water
        if( l1rec->land[ipix] ) stype = 1;  // land
        if( l1rec->ice[ipix] ) stype = 2;  // snow, ice has precedence

//  Looks like the xa computation may be able to be done here  and applied to both types

        // for this pixel, create the local profiles of p, t, h
        // first level is the surface
        lcl_prof_t[0] = l1rec->sfct[ipix];
        lcl_prof_h[0] = l1rec->height[ipix];
        lcl_p_set[0] = l1rec->sfcp[ipix];

        lcl_lvl_ptr = 0;
        fill_lvl = 0;
        // fill the profiles
        for( int ilvl = 0; ilvl < nlvl_merra2; ilvl ++ )
          {
          // find the 1st level that is not fill and has p > the level p
          // and the h < level h
          if( lcl_lvl_ptr == 0 )
            {
            if( ( l1rec->anc_add->prof_height[ ipix * nlvl_merra2 + ilvl ] 
              != BAD_FLT ) && 
              ( l1rec->anc_add->prof_temp[ ipix * nlvl_merra2 + ilvl ] 
              != BAD_FLT ) && 
              ( merra_p_set[ilvl] < l1rec->sfcp[ipix] ) && 
              ( l1rec->anc_add->prof_height[ ipix * nlvl_merra2 + ilvl ] 
              > l1rec->height[ipix] ) )
              {
              fill_lvl = 1;
              lcl_lvl_ptr++;
              }
            }

          if( fill_lvl == 1 )
            {
            lcl_prof_t[lcl_lvl_ptr] = l1rec->anc_add->prof_temp[ ipix * nlvl_merra2 + ilvl ];
            lcl_prof_h[lcl_lvl_ptr] = l1rec->anc_add->prof_height[ ipix * nlvl_merra2 + ilvl ];
            lcl_p_set[lcl_lvl_ptr] = merra_p_set[ilvl];
            lcl_lvl_ptr++;
            }
          }
        lcl_nlvl = lcl_lvl_ptr;
        // End local profile set-up
       /*
        *  process for each cloud top particle type (phase)
        */
        for( ityp = 0; ityp < ntyp; ityp++ )
          {
         /*
          *  Set up the inputs for optimization
          */
          //  Sy the uncerainty
          for( ib1= 0; ib1 < n_ctht_bnd; ib1++ )
            {
            for( ib2 = 0; ib2 < n_ctht_bnd; ib2++ )
              {
              *( sy + ib1 + n_ctht_bnd * ib2 ) = 0.;
              if( ib1 == ib2 )
                {
                // add random error
                *( sy + ib1 + n_ctht_bnd * ib2 ) += unc.noise_err * 
                  unc.noise_err * rhot[ib1] * rhot[ib2];
                // add calib error
                *( sy + ib1 + n_ctht_bnd * ib2 ) += unc.cal_err * unc.cal_err * 
                  rhot[ib1] * rhot[ib2];
                }
              ev1 = unc.icpts[ ityp + ntyp * ( stype + unc.nsurf * ib1 ) ] +
                unc.grads[ ityp + ntyp * ( stype + unc.nsurf * ib1 ) ] *
                rhot[ib1];
              ev2 = unc.fit_rms[ ityp + ntyp * ( stype + unc.nsurf * ib1 ) ];
              err1 = ( ev1 > ev2 ) ? ev1 : ev2;

              ev1 = unc.icpts[ ityp + ntyp * ( stype + unc.nsurf * ib2 ) ] +
                unc.grads[ ityp + ntyp * ( stype + unc.nsurf * ib2 ) ] *
                rhot[ib2];
              ev2 = unc.fit_rms[ ityp + ntyp * ( stype + unc.nsurf * ib2 ) ];
              err2 = ( ev1 > ev2 ) ? ev1 : ev2;

              *( sy + ib1 + n_ctht_bnd * ib2 ) += err1 * err2 *
                unc.err_corrs_category[ityp + ntyp * 
                ( stype + unc.nsurf * ( ib1 + n_ctht_bnd * ib2 ) ) ];
              }
            }
         /*
          *  get the glint (Andy's version of Cox, Munk)
          */
          glint_ref = ctht_glint( solz, senz, relaz, windsp );
         /*
          *  compute the xa - a-priori state vector [ hgt, p, albedo ]
          */
          double max_glint_for_table;
          max_glint_for_table = lut->lut_alb[lut->lut_nalb-1];
          xa[0] = 3.;  xa[1] = 3.;
          if( stype == 0 ) //  water
            {
            xa[2] = ( 0.02 + glint_ref ) / lut_dalb;
            }
          else if( stype == 1 )  // land
            {
            xa[2] = l1rec->cld_dat->cth_alb_init[ipix] / lut_dalb;
            }
          else  // snow/ice
            {
            xa[2] = 0.9 / lut_dalb;
            }
          if( xa[2] >= max_glint_for_table / lut_dalb )
            {
            //WDR avoid these as inst can view direct sun glint
            //printf( 
            //  "%s, %d: W: a-priori albedo exceeds table range, set to %f\n", 
            //  __FILE__, __LINE__, max_glint_for_table );
            xa[2] = max_glint_for_table / lut_dalb;
            }
          if( xa[2] < 0 )
            {
            //printf( 
            //  "%s, %d: W: a-priori albedo below table range, set to 0\n", 
            //  __FILE__, __LINE__ );
            xa[2] = 0.;
            }
         /*
          *  make the covariance matrix sa
          */
          gsl_matrix_set_zero( gm_sa );
          // Very weak prior on COD
          gsl_matrix_set( gm_sa, 0, 0, pow( 100., 2 ) ); 
          // Very weak prior on CTH
          gsl_matrix_set( gm_sa, 1, 1, pow( 100., 2 ) ); 

          if( stype == 0 ) //  water
            {
            // Stronger prior on albedo for water. This is +/-0.01. + 
            // 20% of the glint strength
            gsl_matrix_set( gm_sa, 2, 2, 
              pow( ( 0.01 + 0.2 * glint_ref ) / lut_dalb, 2. ) );
            }
          else if( stype == 1 )  // land
            {
            gsl_matrix_set( gm_sa, 2, 2,
              pow( ( l1rec->cld_dat->cth_alb_unc_init[ipix] / lut_dalb ), 2. ) );
            // Prior comes from climatology
            }
          else
            {
            gsl_matrix_set( gm_sa, 2, 2, pow( 0.05 / lut_dalb, 2. ) );
            // Overwrite if snow/ice cover > 50%,
            // assume 0.05 uncertainty on it
            }
            
         /*
          *  interpolate the lut to the current senz, solz, relaz, sfc press
          */
          lut = ( ityp == 0 ) ? &full_lut_wtr : &full_lut_ice;

          sxax[0] = lut->lut_press;
          snax[0] = lut->lut_npress;
          sxax[1] = lut->lut_solz;
          snax[1] = lut->lut_nsolz;
          sxax[2] = lut->lut_senz;
          snax[2] = lut->lut_nsenz;
          sxax[3] = lut->lut_relaz;
          snax[3] = lut->lut_nrelaz;
          n_block = lut->lut_ncod * lut->lut_ncth * lut->lut_nalb * 
            lut->lut_nwave;

          // For the idl match up work, I've been using l1rec->pr[ipix]
          // as sfc pressure - seems to match what Andy has?  but Andy
          // did not want pr as it uses height and sea lvl p to get the sfc p
          // and so wanted to use sfcp - ? a bit of an odd match.  Just the 
          // same, the surface p is used in at least 3 parts below.  I'll 
          // just declare spress_out here as one or the other
          // spress_out = l1rec->pr[ipix]; // used for idl match work
          spress_out = l1rec->sfcp[ipix];
          float xpt[4] = { spress_out, solz, senz, relaz };

          //  looks like the delphi can be outside range 0 - 180, so:
          xpt[3] = fabs( fmod( xpt[3], 180. ) );
          //  This is a quick fix to get things running, the interpolation
          //  that got this bad p should get checked
          xpt[0] = ( xpt[0] > 1050. ) ? 1050. : xpt[0];
          xpt[0] = ( xpt[0] < 600. ) ? 600. : xpt[0];

          //  set up. interpolate the LUT to the current point's 3 view 
          //  angles and the surface p
          //  get the weights and offsets first
//  **** this can be moved outside the 2 LUT types
          if( int_4d( sxax, snax, n_block, xpt, ar_off, ar_wt, orig_wt ) != 0 )
            {
            //printf( "%s - %d: E: Failure in int_4d\n", __FILE__, __LINE__ );
            //for( ibnd = 0; ibnd < 4; ibnd++ )
            //  printf( "x index: %d, value: %f, lo: %f, hi: %f\n", ibnd, 
            //    xpt[ibnd], sxax[ibnd][0], sxax[ibnd][snax[ibnd]-1] );
            break;  // leave the prod arrays as bad for this point
            }
          //  apply the weights to the LUT and reduce it from 9 to 4 dimensions
          //  [ band, cod, cth, alb ]
          //  apply the weights (and re-organize to have bands slowest)
          for( ibnd = 0; ibnd < lut->lut_nwave; ibnd++ )
            {
            for( icod = 0; icod < lut->lut_ncod; icod++ )
              {
              for( icth = 0; icth < lut->lut_ncth; icth++ )
                {
                for( ialb = 0; ialb < lut->lut_nalb; ialb++ )
                  {
                 tmp_off = icod + lut->lut_ncod * ( icth + lut->lut_ncth *
                    ( ialb + lut->lut_nalb * ibnd ) );
                 tmp_off2 = ialb + lut->lut_nalb * ( icth + lut->lut_ncth *
                    ( icod + lut->lut_ncod * ibnd ) );
                  *( tmp_lut + tmp_off ) = 0;
                  for( icor = 0; icor < 16; icor++ )
                    {
                    *( tmp_lut + tmp_off ) += ar_wt[icor] * 
                      *(lut->lut + ar_off[icor] + tmp_off2 );
                    }
                  }
                }
              }
            }
         /*
          *  get the first guess - the last piece in the puzzle
          */
          min_cost = 1.e11;  // This should go down in checks below
          for( icod = 0; icod < lut->lut_ncod; icod++ )
            {
            for( icth = 0; icth < lut->lut_ncth; icth++ )
              {
              //  set the fg_cost to 1.e10, may shift to BAD later
              for( ialb = 0; ialb < lut->lut_nalb; ialb++ )
                {
                // As orig: off_cost = ialb + lut->lut_nalb * 
                //  ( icth + lut->lut_ncth * icod );
                off_cost = icod + lut->lut_ncod *
                  ( icth + lut->lut_ncth * ialb );
                *( fg_cost + off_cost ) = 1.e10;
                }
              ialb = (int32_t) ( xa[2] + .5 );  // a-priori albedo
              // I orig thought off_cost = ialb + lut->lut_nalb * 
              // ( icth + lut->lut_ncth * icod );
              //  but below is correct
              off_cost = icod + lut->lut_ncod * ( icth + lut->lut_ncth * ialb );
              for( ibnd = 0; ibnd < lut->lut_nwave; ibnd++ )
                {
                tmp_off = icod + lut->lut_ncod * ( icth + lut->lut_ncth *
                    ( ialb + lut->lut_nalb * ibnd ) );
                yres[ibnd] = *( tmp_lut + tmp_off ) - rhot[ibnd];
                //  ABOVE looks a lot like it could go in loop further up!!!
                }
              xres[0] = icod - xa[0];
              xres[1] = icth - xa[1];
              xres[2] = ialb - xa[2];

              //  to replace Andy's code line:
              //  fg_cost[aa,bb,cc]=transpose(yres)#invert(Sy)#yres + $
              //  xres#invert(Sa)#xres
              //  make the matrix for sy and the vector for yres
              for( ibnd = 0; ibnd < lut->lut_nwave; ibnd++ )
                {
                gsl_vector_set( gv_yres, ibnd, yres[ibnd] );
                for( ibnd2= 0; ibnd2 < lut->lut_nwave; ibnd2++ )
                  {
                  gsl_matrix_set( gm_sy, ibnd, ibnd2, 
                    *( sy + ibnd + 8 * ibnd2 ) );
                  }
                }

              //  term 1 of fg_cost:  transpose(yres)#invert(Sy)#yres portion
              gm_sy_inv = invert_a_matrix( gm_sy );
              gsl_blas_dgemv( CblasNoTrans, 1., gm_sy_inv, gv_yres, 0., gv_sy_inv_yres );
              gsl_matrix_free( gm_sy_inv );
              gsl_blas_ddot( gv_sy_inv_yres, gv_yres, &fg_cost_t1 );
              //  
              for( ibnd = 0; ibnd < 3; ibnd++ )
                {
                gsl_vector_set( gv_xres, ibnd, xres[ibnd] );
                }

              //  term 2 of fg_cost: xres#invert(Sa)#xres portion
              gm_sa_inv = invert_a_matrix( gm_sa );
              gsl_blas_dgemv( CblasNoTrans, 1., gm_sa_inv, gv_xres, 0., 
                gv_sa_inv_xres );
              gsl_matrix_free( gm_sa_inv );
              gsl_blas_ddot( gv_sa_inv_xres, gv_xres, &fg_cost_t2 );
              // sum the 2 values
              // fg_cost[ialb,icth,icod] = fg_cost_t1 + fg_cost_t2;
              *( fg_cost + off_cost ) = fg_cost_t1 + fg_cost_t2;

              //  end Andy's code line:
             /*
              *  find the locatuion and value of the fg_cost minimum
              */
              if( *( fg_cost + off_cost ) < min_cost )
                {
                min_cost = *( fg_cost + off_cost );
                min_cost_loc[0] = icod;
                min_cost_loc[1] = icth;
                min_cost_loc[2] = ialb;
                }
              //
              }
            }
         /*
          *  check the min cost so it is < 1.e10
          */
          if( min_cost >= 1.e10 )
            {
            printf( "%s, %d, I: the minimum cost is 1.e10. or higher\n",
              __FILE__, __LINE__ );
            }
         /*
          *  perform the optimization, (ams_oe_inversion.pro part, maybe 
          *  supplied by Laughlin
          */
          // oe_out=ams_oe_inversion(ref_true,Sy,xa,Sa,tmp_lut,fg=inds_min)
          // TO
          lut_dims[0] = lut->lut_ncod;
          lut_dims[1] = lut->lut_ncth;
          lut_dims[2] = lut->lut_nalb;
          lut_dims[3] = lut->lut_nwave;
          ams_oe_inversion( rhot, gm_sy, xa, gm_sa, tmp_lut, 
            lut_dims, min_cost_loc, &oe_info );
         /*
          *  post-optimization work
          */
          //  get the retrieved parameters - xlate from index to value
          //  WE'LL try a linear interp to emulate interpolate as Andy does 
          //  26 calls to interpolate.  Try on the 1st 
          //  of his
          //  Do: ret_lcod[out_inds[0],out_inds[1],t] =
          //           interpolate( lut_lcod, oe_out.x[0] )

          // Retrieved parameters
          // ret_lcod
          ctht_lin->cth_lcod_all[ ityp + ntyp * ipix ] = 
           axis_interp( lut->lut_lcod, lut->lut_ncod, oe_info.x_prd_state[0] );
          // ret_cod
          ctht_lin->cth_cod_all[ ityp + ntyp * ipix ] = 
            pow( 10., ctht_lin->cth_lcod_all[ ityp + ntyp * ipix ] );
          //  for height, Use l1rec 'height' 
          // ret_cth
          ctht_lin->cth_all[ ityp + ntyp * ipix ] =
            axis_interp( lut->lut_cth, lut->lut_ncth, oe_info.x_prd_state[1] )
            + l1rec->height[ipix] / 1000.;
          // ret_ctp
          ctht_lin->ctp_all[ ityp + ntyp * ipix ] =
            axis_interp( lut->lut_pp0, lut->lut_ncth, oe_info.x_prd_state[1] ) *
            spress_out;
          // ret_alb
          ctht_lin->cth_alb_all[ ityp + ntyp * ipix ] =
            axis_interp( lut->lut_alb, lut->lut_nalb, oe_info.x_prd_state[2] );
          // Diagnostics
          // ret_ss[npix,nlin,nphase]  - cost
          ctht_lin->cost_ss[ ityp + ntyp * ipix ] = oe_info.cost;
          // ret_acost[npix,nlin,nphase]
          ctht_lin->acost[ ityp + ntyp * ipix ] = oe_info.acost;
          // ret_iter[npix,nlin,nphase]
          ctht_lin->nitr[ ityp + ntyp * ipix ] = oe_info.nitr;

          if (oe_info.nitr == 1 && oe_info.cost > 30){
            l1rec->flags[ipix] |= CHLFAIL;
          }

          // Avg kernel - HOLD outputting this, need yet another 3rd dim size
          // oe_akdiag[npix,nlin,nphase,3]
          int nxx = 3;
          for( int xx = 0; xx < nxx; xx++ )
            ctht_lin->oe_akdiag[ xx + nxx * ( ityp + ntyp * ipix ) ] = 
              gsl_matrix_get( oe_info.gm_ak, xx, xx );

          // Uncertainty estimates
          // ret_dlcod[npix,nlin,nphase]
          /*  initial code:
          ul=10.^interpolate(lut_lcod,oe_out.x[0]+sqrt(oe_out.sx[0,0]) < $
            (lut_ncod-1)) - 10.^interpolate(lut_lcod,oe_out.x[0])
          ll=10.^interpolate(lut_lcod,oe_out.x[0]-sqrt(oe_out.sx[0,0]) > 0 ) - $
            10.^interpolate(lut_lcod,oe_out.x[0])
          ret_dcod[out_inds[0],out_inds[1],t]=max(abs([ul,ll]))

          ul=interpolate(lut_lcod,oe_out.x[0]+sqrt(oe_out.sx[0,0]) < $
            (lut_ncod-1)) - interpolate(lut_lcod,oe_out.x[0])
          ll=interpolate(lut_lcod,oe_out.x[0]-sqrt(oe_out.sx[0,0]) > 0 ) - $
            interpolate(lut_lcod,oe_out.x[0])
          ret_dlcod[out_inds[0],out_inds[1],t]=max(abs([ul,ll]))
          */
          // ret_dcod[npix,nlin,nphase]
          coord = oe_info.x_prd_state[0] + 
            sqrt( gsl_matrix_get( oe_info.gm_sx, 0, 0 ) );
          if( lut->lut_ncod - 1 < coord ) coord = lut->lut_ncod - 1;
          int_ctr = axis_interp( lut->lut_lcod, lut->lut_ncod, 
            oe_info.x_prd_state[0] );
          int_pls = axis_interp( lut->lut_lcod, lut->lut_ncod, coord );

          coord = oe_info.x_prd_state[0] -
            sqrt( gsl_matrix_get( oe_info.gm_sx, 0, 0 ) );
          if( coord < 0 ) coord = 0;
          int_neg = axis_interp( lut->lut_lcod, lut->lut_ncod, coord ); 

          ctht_lin->dlcod[ ityp + ntyp * ipix ] = 
            fmax( fabs( int_pls - int_ctr ), fabs( int_neg -int_ctr ) );
          //  ret_dcod[npix,nlin,nphase]
          pctr = pow( 10., int_ctr );
          int_pls = pow( 10., int_pls );
          int_neg = pow( 10., int_neg );
          ctht_lin->dcod[ ityp + ntyp * ipix ] =
            fmax( fabs( int_pls - pctr ), fabs( int_neg - pctr ) );

          // ret_dcth[npix,nlin,nphase]
          coord = oe_info.x_prd_state[1] +
            sqrt( gsl_matrix_get( oe_info.gm_sx, 1, 1 ) );
          if( lut->lut_ncth - 1 < coord ) coord = lut->lut_ncth - 1;
          int_ctr = axis_interp( lut->lut_cth, lut->lut_ncth,
            oe_info.x_prd_state[1] );
          int_pls = axis_interp( lut->lut_cth, lut->lut_ncth, coord );

          coord = oe_info.x_prd_state[1] -
            sqrt( gsl_matrix_get( oe_info.gm_sx, 1, 1 ) );
          if( coord < 0 ) coord = 0;
          int_neg = axis_interp( lut->lut_cth, lut->lut_ncth, coord );

          ctht_lin->dcth[ ityp + ntyp * ipix ] =
            fmax( fabs( int_pls - int_ctr ), fabs( int_neg - int_ctr ) );

          // Retrieved parameters (more which require MERRA2 profiles)
          // ret_cth_raw[nphase, npix,nlin]
          // the cth_raw will be the non-MERRA2 -derived height, in case 
          // first, save the cth_all
          ctht_lin->cth_raw[ ityp + ntyp * ipix ] =
            ctht_lin->cth_all[ ityp + ntyp * ipix ];
            //  find start of interval in MERRA2 levels containg the ctp
          close_ind = -1;
          for( ilev = 0; ilev < lcl_nlvl; ilev++ )
            {
            if( lcl_p_set[ilev] < ctht_lin->ctp_all[ ityp + ntyp * ipix ] )
              {
              close_ind = ( ilev == 0 ) ? 0 : ilev -1;
              break;
              }
            }
          if( ( close_ind == -1 ) || ( close_ind == lcl_nlvl-1 ) )
            {
            p_log_lo = log10( lcl_p_set[lcl_nlvl-2] );
            p_log_hi = log10( lcl_p_set[lcl_nlvl-1] );
            }
          else
            {
            p_log_lo = log10( lcl_p_set[close_ind] );
            p_log_hi = log10( lcl_p_set[ close_ind + 1 ] );
            }
          dp = ( log10( ctht_lin->ctp_all[ ityp + ntyp * ipix ] ) - p_log_lo )
            / ( p_log_hi - p_log_lo );

          // ret_cth (another calc using MERRA2)[npix,nlin,nphase] above
          // Note that we make the height into km as Andy has
          found_hgt = axis_interp( lcl_prof_h, lcl_nlvl, close_ind + dp );

          if( found_hgt < l1rec->height[ipix] ) found_hgt = 
            l1rec->height[ipix]  + 10.;
          ctht_lin->cth_all[ ityp + ntyp * ipix ] = found_hgt / 1000.;

          // ret_ctt (calc using MERRA2)[npix,nlin,nphase]
          ctht_lin->ctt_all[ ityp + ntyp * ipix ] = C_TO_K + 
            axis_interp( lcl_prof_t, lcl_nlvl, close_ind + dp );
          // <a cth sanity check>
          if( ctht_lin->cth_all[ ityp + ntyp * ipix ] > 100. )
            {
            printf( "%s, %d, W: crazy CTH, profile issue, will use non MERRA\n",
              __FILE__, __LINE__ );
            ctht_lin->cth_all[ ityp + ntyp * ipix ] = 
              ctht_lin->cth_raw[ ityp + ntyp * ipix ];
            }
          // <check the ctt, if > 350, use sfc t (10 m T )>
          if( ctht_lin->ctt_all[ ityp + ntyp * ipix ] >= 350. )
            ctht_lin->ctt_all[ ityp + ntyp * ipix ] = 
            l1rec->sfct[ipix] + C_TO_K;

          // Uncertainty estimates, round 2
          // ret_dctp[npix,nlin,nphase]
          // coord: oe_out.x[1]+sqrt(oe_out.sx[1,1])
          // or lut_ncth-1, which ever lower
          coord = oe_info.x_prd_state[1] +
            sqrt( gsl_matrix_get( oe_info.gm_sx, 1, 1 ) );
          if( lut->lut_ncth - 1 < coord ) coord = lut->lut_ncth - 1;
          int_ctr = axis_interp( lut->lut_pp0, lut->lut_ncth, 
            oe_info.x_prd_state[1] ) * spress_out;
          int_pls = axis_interp( lut->lut_pp0, lut->lut_ncth, coord ) *
            spress_out;
          
          coord = oe_info.x_prd_state[1] -
            sqrt( gsl_matrix_get( oe_info.gm_sx, 1, 1 ) );
          if( coord < 0 ) coord = 0;
          int_neg = axis_interp( lut->lut_pp0, lut->lut_ncth, coord ) *
            spress_out;

          ctht_lin->dctp[ ityp + ntyp * ipix ] =
            fmax( fabs( int_pls - int_ctr ), fabs( int_neg - int_ctr ) );

          // ret_dctt[npix,nlin,nphase] a bit easier, Likewise no super 
          // easy way to make an equivalent for CTT uncertainty.
          // Since average lapse rate in troposphere is 6 K per km,
          // assume that holds, and scale based on CTH uncertainty.
          ctht_lin->dctt[ ityp + ntyp * ipix ] = 
            6. * ctht_lin->dcth[ ityp + ntyp * ipix ];

          // ret_dalb[npix,nlin,nphase]
          coord = oe_info.x_prd_state[2] +
            sqrt( gsl_matrix_get( oe_info.gm_sx, 2, 2 ) );
          if( lut->lut_nalb - 1 < coord ) coord = lut->lut_nalb - 1;
          int_ctr = axis_interp( lut->lut_alb, lut->lut_nalb, 
            oe_info.x_prd_state[2] );
          int_pls = axis_interp( lut->lut_alb, lut->lut_nalb, coord );

          coord = oe_info.x_prd_state[2] -
            sqrt( gsl_matrix_get( oe_info.gm_sx, 2, 2 ) );
          int_neg = axis_interp( lut->lut_alb, lut->lut_nalb, coord );

          ctht_lin->dalb[ ityp + ntyp * ipix ] =
            fmax( fabs( int_pls - int_ctr ), fabs( int_neg - int_ctr ) );

          }   // end cloud top particle type
        //  determine the 'best' phase  - the phase will control what is
        //  output
        lo_cost = 1.e20;
        found_phase = BAD_UBYTE;
        for( ityp = 0; ityp < ntyp; ityp++ )
          {
          if( ctht_lin->cost_ss[ ityp + ntyp * ipix ] == BAD_FLT )
            {
            found_phase = BAD_UBYTE;
            break;
            }
          else
            {
            if( ctht_lin->cost_ss[ ityp + ntyp * ipix ] < lo_cost )
              {
              found_phase = ityp;
              lo_cost = ctht_lin->cost_ss[ ityp + ntyp * ipix ];
              }
            }
          }
        ctht_lin->ct_phase[ipix] = found_phase;
        //   set the main 3 products that go into the CMP code
        if( found_phase != BAD_UBYTE )
          {
          ctht_lin->cth[ipix] = ctht_lin->cth_all[ found_phase + ntyp * ipix ];
          ctht_lin->ctp[ipix] = ctht_lin->ctp_all[ found_phase + ntyp * ipix ];
          ctht_lin->ctt[ipix] = ctht_lin->ctt_all[ found_phase + ntyp * ipix ];
          }
        else
          {
          ctht_lin->cth[ipix] = BAD_FLT;
          ctht_lin->ctp[ipix] = BAD_FLT;
          ctht_lin->ctt[ipix] = BAD_FLT;
          }
        }  // end cloudy pixel work,  may need to set some BADs here too 
      }  //  end pixel loop
// WDR temp, stop after 1st line
// printf( "WDR - get_ctht early end\n" );
// exit(8);

    //  set scan and # pix
    ctht_lin->iscan = l1rec->iscan;
    ctht_lin->npix = l1rec->npix;
    }
  return 0;
  }
int32_t init_ctht_parms( ctht_lin_str *ctht_lin, int32_t ipix )
 /*
  *  just initialize the cloud top height values
  *
  *  Return              0 - all well
  *  ctht_lin_str *     ctht_lin     structure of all ctht info
  *  int32_t            ipix         pixell to set values at
  */
  {
  int32_t ntyp = 2;

  //  basic values needed by chimaera code
  ctht_lin->cth[ipix] = BAD_FLT;
  ctht_lin->ctp[ipix] = BAD_FLT;
  ctht_lin->ctt[ipix] = BAD_FLT;
  //  other cloud height values
  ctht_lin->ct_phase[ipix] = BAD_UBYTE;
  
  //  'all' from all phases
  for( int i = 0; i < ntyp; i++ ) 
    {
    ctht_lin->cth_all[ i + ntyp * ipix ] = BAD_FLT;
    ctht_lin->ctp_all[ i + ntyp * ipix ] = BAD_FLT;
    ctht_lin->ctt_all[ i + ntyp * ipix ] = BAD_FLT;
    ctht_lin->cth_cod_all[ i + ntyp * ipix ] = BAD_FLT;
    ctht_lin->cth_lcod_all[ i + ntyp * ipix ] = BAD_FLT;
    ctht_lin->cth_raw[ i + ntyp * ipix ] = BAD_FLT;
    ctht_lin->cth_alb_all[ i + ntyp * ipix ] = BAD_FLT;
    ctht_lin->cost_ss[ i + ntyp * ipix ] = BAD_FLT;
    ctht_lin->acost[ i + ntyp * ipix ] = BAD_FLT;
    ctht_lin->nitr[ i + ntyp * ipix ] = BAD_INT;
    for( int inx = 0; inx < 3; inx++ )
      ctht_lin->oe_akdiag[ inx * ( i + ntyp * ipix ) ] = BAD_FLT;
    ctht_lin->dlcod[ i + ntyp * ipix ] = BAD_FLT;
    ctht_lin->dcod[ i + ntyp * ipix ] = BAD_FLT;
    ctht_lin->dcth[ i + ntyp * ipix ] = BAD_FLT;
    ctht_lin->dctp[ i + ntyp * ipix ] = BAD_FLT;
    ctht_lin->dctt[ i + ntyp * ipix ] = BAD_FLT;
    ctht_lin->dalb[ i + ntyp * ipix ] = BAD_FLT;
    }
  return 0;
  }
int32_t comp_ctht_tst(l1str *l1rec, int32_t proc_flg, ctht_lin_str *ctht_lin )
  /*
   comp_ctht_tst will supply test ctht data if called for a line of 
   satellite data

   Returns -  int32_t  - 0 for all good , only status for now

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
   l1str *      l1rec    input l1 record to use
   int32_t      proc_flg   either 0 to return constants, or 1 to return 
                           values based on scan mod 3
   ctht_lin_str *ctht_lin  the cloud top height record to fill

   W. Robinson, SAIC, 22 Jun 2023
  */
  {
  int32_t ipix, npix;
  float cth, ctp, ctt;
  
  npix = l1rec->npix;
  if( proc_flg == 0 )
    {
    cth = 3.;
    ctp = 480.;
    ctt = 250.;
    }
  else if( proc_flg== 1 )
    {
    switch ( l1rec->iscan % 3 )
      {
      case 0:
        cth = 2.;
        ctp = 600.;
        ctt = 270.;
        break;
      case 1:
        cth = 3.;
        ctp = 480.;
        ctt = 250.;
        break;
      case 2:
        cth = 4.;
        ctp = 400.;
        ctt = 230.;
        break;
      }
    }
  for( ipix = 0; ipix < npix; ipix++ )
    {
    ctht_lin->cth[ipix] = cth;
    ctht_lin->ctp[ipix] = ctp;
    ctht_lin->ctt[ipix] = ctt;
    }
  ctht_lin->iscan = l1rec->iscan;
  ctht_lin->npix = l1rec->npix;
  return 0;
  }

int32_t ctht_tbl_init( int32_t nband, l1str *l1rec, ctht_unc_str *unc, 
  ctht_lut_str *full_lut_wtr, ctht_lut_str *full_lut_ice )
 /*
  *  ctht_tbl_init - initialize the cloud height tables
  *
  *  returns int32_t  status: 0 good else bad
  *  int32_t   nband   input expected # bands
  *  l1str *   l1rec, the L1 record
  *  ctht_unc_str *unc  output filled uncertainty information
  *  ctht_lut_str *full_lut_wtr  output filled LUT information for water
  *  ctht_lut_str *full_lut_ice  output filled LUT information for ice
  */
  {
  int ncid, dim_loc, varid;
  char *env_arr;
  char unc_file[FILENAME_MAX], lut_wtr_file[FILENAME_MAX], lut_ice_file[FILENAME_MAX];
  char tit_unc[] = { "Forward model uncertainty coefficients for cloud top pressure retrieval from oci" };
  char tit_lut[] = { "LUT for cloud top pressure retrieval from PACE OCI" };
  ctht_lut_str *lut_ptr;
  char *lut_fil;
  int32_t ilut, ilev, ihgt, close_ind;
  float pstep, frac_ind;

 /*  set the path to the tables and table names */
  if ((env_arr = getenv("OCDATAROOT")) == NULL) {
      printf("-E- %s, %d:  Error looking up environmental variable OCDATAROOT\n", __FILE__, __LINE__ );
      exit(1);
  }
  strcpy( unc_file, env_arr );
  strcat( unc_file, "/cloud/" );
  strcat( unc_file, sensorId2SensorDir(l1rec->l1file->sensorID));
  strcat( unc_file, "/" );
  strcpy( lut_wtr_file, unc_file );
  strcpy( lut_ice_file, unc_file );
  //  make individual file names here
  strcat( unc_file, "oci_fm_err_fits_lut_v17_sim_v5.nc" );
  strcat( lut_wtr_file, "oci_o2_cth_lut_v17_water.nc" );
  strcat( lut_ice_file, "oci_o2_cth_lut_v17_column_8elements.nc" );

 /*  open the uncertainty file */
  printf( "%s, %d: I: Reading cloud top uncertainty file: %s\n", 
    __FILE__, __LINE__, unc_file );
  DPTB( nc_open( unc_file, 0, &ncid ) );
 /*  check the title for a match */
 /*  must do the below for a string attr, not text */
  size_t attlen = 0, dim_len;
  DPTB( nc_inq_attlen(ncid, NC_GLOBAL, "description", &attlen) );
  char **string_attr = (char**)malloc(attlen * sizeof(char*));
  memset(string_attr, 0, attlen * sizeof(char*));
  DPTB( nc_get_att_string(ncid, NC_GLOBAL, "description", string_attr) );

  if( strncmp( string_attr[0], tit_unc, strlen(tit_unc) ) != 0 )
    {
    printf( "%s, %d: E - Description of ctht uncertainty file is bad\n", 
      __FILE__, __LINE__ );
    exit(1);
    }
 /*  read sizes */
  DPTB( nc_inq_dimid( ncid, "nphase", &dim_loc ) );
  DPTB( nc_inq_dimlen( ncid, dim_loc, &dim_len ) );
  unc->nphase = dim_len;
  
  DPTB( nc_inq_dimid( ncid, "nsurf", &dim_loc ) );
  DPTB( nc_inq_dimlen( ncid, dim_loc, &dim_len ) );
  unc->nsurf = dim_len;

  DPTB( nc_inq_dimid( ncid, "nbands", &dim_loc ) );
  DPTB( nc_inq_dimlen( ncid, dim_loc, &dim_len ) );
  unc->nbands = dim_len;
 /*  read axis coordinates */
 /*  in this case the wavelengths 'better' match ctht_wav, we have 2 phases
     water and ice and we have 3 surfaces, 0 water, 1 land, 2 snowice */
 /*  read the uncertainty arrays */
  DPTB( nc_inq_varid( ncid, "icpts", &varid ) );
  unc->icpts = (float *)malloc( unc->nphase * unc->nsurf * unc->nbands * 
    sizeof(float) );
  DPTB( nc_get_var_float( ncid, varid, unc->icpts ) );

  DPTB( nc_inq_varid( ncid, "grads", &varid ) );
  unc->grads = (float *)malloc( unc->nphase * unc->nsurf * unc->nbands *
    sizeof(float) );
  DPTB( nc_get_var_float( ncid, varid, unc->grads ) );

  DPTB( nc_inq_varid( ncid, "fit_rms", &varid ) );
  unc->fit_rms = (float *)malloc( unc->nphase * unc->nsurf * unc->nbands *
    sizeof(float) );
  DPTB( nc_get_var_float( ncid, varid, unc->fit_rms ) );

  DPTB( nc_inq_varid( ncid, "err_corrs_category", &varid ) );
  unc->err_corrs_category = (float *)malloc( unc->nphase * unc->nsurf * 
    unc->nbands * unc->nbands * sizeof(float) );
  DPTB( nc_get_var_float( ncid, varid, unc->err_corrs_category ) );
  DPTB( nc_close( ncid ) );
 /*  
  *  open the LUT file for water and ice
  */ 
  for( ilut = 0; ilut < 2; ilut++ )
    {
    lut_ptr = ( ilut == 0 ) ? full_lut_wtr : full_lut_ice;
    lut_fil = ( ilut == 0 ) ? lut_wtr_file : lut_ice_file;

    printf( "%s, %d: I: Reading cloud top phase LUT file: %s\n",
    __FILE__, __LINE__, lut_fil );
    DPTB( nc_open( lut_fil, 0, &ncid ) );
   /*  check the title for a match */
    size_t attlen = 0, dim_len;
    DPTB( nc_inq_attlen(ncid, NC_GLOBAL, "description", &attlen) );
    char **string_attr = (char**)malloc(attlen * sizeof(char*));
    memset(string_attr, 0, attlen * sizeof(char*));
    DPTB( nc_get_att_string(ncid, NC_GLOBAL, "description", string_attr) );

    if( strncmp( string_attr[0], tit_lut, strlen(tit_lut) ) != 0 )
      {
      printf( "%s, %d: E - Description of ctht lut file is bad\n", 
        __FILE__, __LINE__ );
      printf( "For file: %s\n", lut_fil );
      exit(1);
      }
   /*  read sizes */
    DPTB( nc_inq_dimid( ncid, "nsza", &dim_loc ) );
    DPTB( nc_inq_dimlen( ncid, dim_loc, &dim_len ) );
    lut_ptr->lut_nsolz = dim_len;

    DPTB( nc_inq_dimid( ncid, "nvza", &dim_loc ) );
    DPTB( nc_inq_dimlen( ncid, dim_loc, &dim_len ) );
    lut_ptr->lut_nsenz = dim_len;

    DPTB( nc_inq_dimid( ncid, "nazi", &dim_loc ) );
    DPTB( nc_inq_dimlen( ncid, dim_loc, &dim_len ) );
    lut_ptr->lut_nrelaz = dim_len;

    DPTB( nc_inq_dimid( ncid, "nlevs", &dim_loc ) );
    DPTB( nc_inq_dimlen( ncid, dim_loc, &dim_len ) );
    lut_ptr->lut_nlev = dim_len;

    DPTB( nc_inq_dimid( ncid, "npress", &dim_loc ) );
    DPTB( nc_inq_dimlen( ncid, dim_loc, &dim_len ) );
    lut_ptr->lut_npress = dim_len;

    DPTB( nc_inq_dimid( ncid, "nalb", &dim_loc ) );
    DPTB( nc_inq_dimlen( ncid, dim_loc, &dim_len ) );
    lut_ptr->lut_nalb = dim_len;

    DPTB( nc_inq_dimid( ncid, "ncod", &dim_loc ) );
    DPTB( nc_inq_dimlen( ncid, dim_loc, &dim_len ) );
    lut_ptr->lut_ncod = dim_len;

    DPTB( nc_inq_dimid( ncid, "naod", &dim_loc ) );
    DPTB( nc_inq_dimlen( ncid, dim_loc, &dim_len ) );
    lut_ptr->lut_naod = dim_len;

    DPTB( nc_inq_dimid( ncid, "nbands", &dim_loc ) );
    DPTB( nc_inq_dimlen( ncid, dim_loc, &dim_len ) );
    lut_ptr->lut_nwave = dim_len;

    DPTB( nc_inq_dimid( ncid, "ncth", &dim_loc ) );
    DPTB( nc_inq_dimlen( ncid, dim_loc, &dim_len ) );
    lut_ptr->lut_ncth = dim_len;

  //printf( "%s, %d: TEMP: nsolz: %d, nsenz: %d, nrelaz: %d\n", __FILE__, __LINE__, lut_ptr->lut_nsolz, lut_ptr->lut_nsenz, lut_ptr->lut_nrelaz );
  //printf( "nlev: %d, npress: %d, nalb: %d, ncod: %d\n", lut_ptr->lut_nlev, lut_ptr->lut_npress, lut_ptr->lut_nalb, lut_ptr->lut_ncod );
  //printf( "naod: %d, nbands: %d, ncth: %d\n", lut_ptr->lut_naod, lut_ptr->lut_nwave, lut_ptr->lut_ncth );

    if( lut_ptr->lut_nwave != nband )
      {
      printf( "%s, %d: E: CTHT LUT # bands  != expected #: %d\n", __FILE__, __LINE__, nband );
      exit(1);
      }
   /*  read axis coordinates */
    lut_ptr->lut_solz = (float *) malloc( lut_ptr->lut_nsolz * sizeof( float) );
    lut_ptr->lut_senz = (float *) malloc( lut_ptr->lut_nsenz * sizeof( float) );
    lut_ptr->lut_relaz = (float *) malloc( lut_ptr->lut_nrelaz * sizeof( float) );
    if( ( lut_ptr->lut_solz == NULL ) || ( lut_ptr->lut_senz == NULL ) || ( lut_ptr->lut_relaz == NULL ) )
      {
      printf( "%s, %d: E: unable to allocate lut_solz, lut_senz, or lut_relaz\n", __FILE__, __LINE__ );
      exit(1);
      }
    DPTB( nc_inq_varid( ncid, "sza", &varid ) );
    DPTB( nc_get_var_float( ncid, varid, lut_ptr->lut_solz ) );

    DPTB( nc_inq_varid( ncid, "vza", &varid ) );
    DPTB( nc_get_var_float( ncid, varid, lut_ptr->lut_senz) );

    DPTB( nc_inq_varid( ncid, "azi", &varid ) );
    DPTB( nc_get_var_float( ncid, varid, lut_ptr->lut_relaz ) );

    lut_ptr->lut_press = (float *) malloc( lut_ptr->lut_npress * sizeof( float) );
    lut_ptr->lut_alb = (float *) malloc( lut_ptr->lut_nalb * sizeof( float) );
    lut_ptr->lut_cod = (float *) malloc( lut_ptr->lut_ncod * sizeof( float) );
    lut_ptr->lut_lcod = (float *) malloc( lut_ptr->lut_ncod * sizeof( float) );
    if( ( lut_ptr->lut_press == NULL ) || ( lut_ptr->lut_alb == NULL ) || ( lut_ptr->lut_cod == NULL ) )
      {
      printf( "%s, %d: E: unable to allocate lut_press, lut_alb, or lut_cod\n", __FILE__, __LINE__ );
      exit(1);
      }
    DPTB( nc_inq_varid( ncid, "spress", &varid ) );
    DPTB( nc_get_var_float( ncid, varid, lut_ptr->lut_press ) );

    DPTB( nc_inq_varid( ncid, "cod", &varid ) );
    DPTB( nc_get_var_float( ncid, varid, lut_ptr->lut_cod ) );

    DPTB( nc_inq_varid( ncid, "lcod", &varid ) );
    DPTB( nc_get_var_float( ncid, varid, lut_ptr->lut_lcod ) );

    DPTB( nc_inq_varid( ncid, "salb", &varid ) );
    DPTB( nc_get_var_float( ncid, varid, lut_ptr->lut_alb ) );

    lut_ptr->lut_aod = (float *) malloc( lut_ptr->lut_naod * sizeof( float) );
    lut_ptr->lut_wave = (float *) malloc( lut_ptr->lut_nwave * sizeof( float) );
    lut_ptr->lut_cth = (float *) malloc( lut_ptr->lut_ncth * sizeof( float) );
    if( ( lut_ptr->lut_aod == NULL ) || ( lut_ptr->lut_wave == NULL ) || ( lut_ptr->lut_cth == NULL ) )
      {
      printf( "%s, %d: E: unable to allocate lut_aod, lut_wave, or lut_cth\n", __FILE__, __LINE__ );
      exit(1);
      }
    DPTB( nc_inq_varid( ncid, "aod", &varid ) );
    DPTB( nc_get_var_float( ncid, varid, lut_ptr->lut_aod ) );

    DPTB( nc_inq_varid( ncid, "band_wavelengths", &varid ) );
    DPTB( nc_get_var_float( ncid, varid, lut_ptr->lut_wave ) );

    DPTB( nc_inq_varid( ncid, "cth", &varid ) );
    DPTB( nc_get_var_float( ncid, varid, lut_ptr->lut_cth ) );

   /*  there is a state_inds array, but avoid this for now */
   /*  read the LUT array and height rel press */
    lut_ptr->z_and_pp0 = (float *) malloc( lut_ptr->lut_nlev * 2 * sizeof( float) );
    lut_ptr->lut_pp0 = (float *) malloc( lut_ptr->lut_ncth * sizeof( float ) );
    lut_ptr->lut = (float *) malloc( lut_ptr->lut_nrelaz * lut_ptr->lut_nsenz 
      * lut_ptr->lut_nsolz * lut_ptr->lut_npress * lut_ptr->lut_nwave * 
      lut_ptr->lut_ncod * lut_ptr->lut_ncth * lut_ptr->lut_nalb * 
      sizeof( float) );
    if( ( lut_ptr->z_and_pp0 == NULL ) || ( lut_ptr->lut == NULL ) ||
        ( lut_ptr->lut_pp0 == NULL ) )
      {
      printf( "%s, %d: E: unable to allocate lut_z_and_pp0, lut_pp0, or lut\n", 
        __FILE__, __LINE__ );
      exit(1);
      }
     DPTB( nc_inq_varid( ncid, "z_and_pp0", &varid ) );
    DPTB( nc_get_var_float( ncid, varid, lut_ptr->z_and_pp0 ) );

    DPTB( nc_inq_varid( ncid, "lut", &varid ) );
    DPTB( nc_get_var_float( ncid, varid, lut_ptr->lut ) );
    //  that should do it
    DPTB( nc_close( ncid ) );
   /*
    * The lut_pp0 is created from the z_and_pp0 - do that here
    */
    float *and_pp0 = (float *) malloc( lut_ptr->lut_nlev * sizeof(float) );
    for( ihgt = 0; ihgt < lut_ptr->lut_nlev; ihgt++ )
      and_pp0[ihgt] = lut_ptr->z_and_pp0[ 1 + 2 * ihgt ];
    for( ilev = 0; ilev < lut_ptr->lut_ncth; ilev++ )
      {
      // get the pressure for all heights in the height axis
      lut_ptr->lut_pp0[ilev] = BAD_FLT;
      close_ind = -1;
      for( ihgt = 0; ihgt < lut_ptr->lut_nlev; ihgt++ )
        {
        if( lut_ptr->z_and_pp0[ ihgt * 2 ] < lut_ptr->lut_cth[ilev] )
          {
          close_ind = ( ihgt == 0 ) ? 0 : ihgt - 1;
          break;
          }
        }
      if( close_ind < lut_ptr->lut_nlev - 1 )
        pstep = lut_ptr->z_and_pp0[ ( close_ind + 1 ) * 2  ] - 
          lut_ptr->z_and_pp0[ close_ind * 2  ];
      else
        pstep = lut_ptr->z_and_pp0[ close_ind * 2  ] -
          lut_ptr->z_and_pp0[ ( close_ind - 1 ) * 2 ];
      frac_ind = close_ind + ( lut_ptr->lut_cth[ilev] - 
        lut_ptr->z_and_pp0[ close_ind * 2 ] ) / pstep;
      lut_ptr->lut_pp0[ilev] = axis_interp( and_pp0, lut_ptr->lut_nlev, 
        frac_ind );
      }
    free(and_pp0);
    }
 /*
  *  Both tables read
  */
  return 0;
  }

float ctht_glint( float solz, float senz, float relaz, float windsp )
 /*
  *  ctht_glint will vget the glint for a pixel and wind speed, adapted from 
  *    calc_sunglint.pro f A. Sayer
  *
  *  return float     glint value as reflectance
  *  float    solz    solar zenith, in deg
  *  float    senz    sensor zenith, degrees
  *  float    relaz   relative azimuth
  *  float    windsp  wind speed, m/s
  *
  *  follows the Cox and Munk 1954 paper
  *
  *  24 Jan 2024, WDR make more variables as double and avoid NaN in 
  *  R calculation
  *
  */
  {
  double cos_two_omega, cos_beta;
  double a1, b1, c1, d1, rgl, wprime, ref_coef;
  float kai, nr;
  double zxprime, zyprime, sigx, sigy, w;
  double zeta, eta, p_gram, zx, zy;
  double dsenz, dsolz, drelaz;

  int32_t do_gordon = 1;

  dsenz = (double) senz * OEL_DEGRAD;
  dsolz = (double) solz * OEL_DEGRAD;
  drelaz = ( do_gordon ) ? ( 180. - (double)relaz ) * OEL_DEGRAD : 
                           (double)relaz * OEL_DEGRAD;

  nr = 1.341;

  //  Specular reflectance from ruffled sea
  kai = 0.;
  rgl = 0.0;
  kai = kai * OEL_DEGRAD;

  //  Defines surface slopes
  zx = ( -1.0 * sin(dsenz) * sin(drelaz) ) / ( cos(dsolz) + cos(dsenz) );
  zy = ( sin(dsolz) + sin(dsenz) * cos( drelaz ) ) / 
    ( cos(dsolz) + cos(dsenz) );

  //  Use these lines to make independent of wind direction
  zxprime = zx;
  zyprime = zy;

  //  Coefficients from Cox and Munk, 1954 (Statistics of...).
  //  Wind isotropic slope as function of find speed.
  sigx = sqrt( 0.003 + 0.00192 * windsp );
  sigy = sqrt( 0.00316 * windsp );

  zeta = zxprime / sigx;
  eta = zyprime / sigy;

  p_gram = ( 1. / ( 2. * OEL_PI * sigx * sigy ) ) * 
    exp( -0.5 * ( pow( zeta, 2. ) + pow( eta, 2. ) ) );

  //  Cox and Munk (1954) geometry (Measurements of...)
  //  2 omega = angle between incident ray (sun) and instrument
  //  with the surface as the sloping sea

  cos_two_omega = cos(dsenz) * cos(dsolz) + 
    sin(dsenz) * sin( dsolz ) * cos(drelaz);
  cos_beta = ( cos(dsolz) + cos(dsenz) ) / ( sqrt( 2 + 2 * cos_two_omega ) );
  w = ( cos_two_omega >= 1. ) ? 0. : 0.5 * acos( cos_two_omega );

  // Fresnel component
  wprime = asin( 1.00029 * sin(w) / nr );
  a1 = sin( w - wprime );
  b1 = sin( w + wprime );
  c1 = tan( w - wprime );
  d1 = tan( w + wprime );
  // make glint refl 0 if b1 or d1 is 0
  if( ( b1 == 0 ) || ( d1 == 0 ) )
    rgl = 0.;
  else {
    ref_coef = 0.5 * ( ( a1 * a1 ) / ( b1 * b1 ) + ( c1 * c1 ) / ( d1 * d1 ) );
    //  compute the glint reflectance
    rgl = OEL_PI * p_gram * ref_coef / ( 4. * cos( dsolz ) * cos( dsenz ) * 
      ( pow( cos_beta, 4. ) ) );
  }
  return rgl;
  }
// the matrix invert:

gsl_matrix *invert_a_matrix(gsl_matrix *matrix)
 /**************************************************
  *  invert_a_matrix was copied as a way to do an inverse of a gsl matrix
  *  The only thing was that it affects the input matris.  So, we copy it 
  *  and work on that
  *************************************************/
{
    size_t size = matrix->size1;
    gsl_permutation *p = gsl_permutation_alloc(size);
    gsl_matrix *gm_work;
    int s;

    //  copy the input matris to work on
    gm_work = gsl_matrix_alloc( size, size );
    gsl_matrix_memcpy( gm_work, matrix );

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(gm_work, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(size, size);
    gsl_linalg_LU_invert(gm_work, p, inv);

    gsl_permutation_free(p);
    gsl_matrix_free( gm_work );

    return inv;
}

//  a diagnostic routines to allow printing matricise and vectors
//  (the gsl_matrix_fprintf is  bit vague)
//  the gsl_vector_fprintf has no such confusion
void
print_mat_contents(gsl_matrix *matrix, int nrow, int ncol )
{
    double element;
    int i, j;
    printf( "%d rows, %d columns\n", nrow, ncol );
    for (i = 0; i < nrow; ++i) {
        for (j = 0; j < ncol; ++j) {
            element = gsl_matrix_get(matrix, i, j);
            printf("%f ", element);
        }
        printf("\n");
    }
}
