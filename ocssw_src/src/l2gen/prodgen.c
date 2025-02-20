/* =========================================================== */
/* Module prodgen.c                                            */
/*                                                             */
/* Function prodgen returns a pointer to the geophysical       */
/* product associated with the input product catalog entry and */
/* the input level-2 record.  The specific product will be     */
/* extracted from the level-2 record or computed using the     */
/* information in the level-2 record and knowledge of the      */
/* required algorithm to be called.  The output is stored in   */
/* a static array, and a pointer is returned to the caller.    */
/*                                                             */
/* Written By:                                                 */
/*     Bryan A. Franz, NASA/OBPG, August 2008.                 */
/* =========================================================== */

#include "l12_proto.h"

#include <stdint.h>
#include <inttypes.h>

static int32 numScans;
static int32 numPixels;
static int32 numBands;
static int32 numBandsIR;


/* ----------------------------------------------------------- */
/* extract_band() - extracts a product from a BIL array        */
/* ----------------------------------------------------------- */
VOIDP extract_band(float *fbuf, l2prodstr *p, int32 nbands) {
    static float32 *fbuf2 = NULL;

    int32 band = MIN(MAX(p->prod_ix, 0), nbands - 1);
    static int32 npix = 0;
    int32 ip;

    if (fbuf2 == NULL) {
        npix = p->dim[1];
        if ((fbuf2 = calloc(npix, sizeof (float32))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    } else if (npix != p->dim[1]) {
        npix = p->dim[1];
        free(fbuf2);
        // this is put here for l3gen which varies the number of pixels on each line

        if ((fbuf2 = (float32 *) calloc(npix, sizeof (float32)))
                == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }

    for (ip = 0; ip < npix && ip < numPixels; ip++) {
        fbuf2[ip] = fbuf[ip * nbands + band];
    }

    return ((VOIDP) fbuf2);
}

/* ----------------------------------------------------------- */
/* extract_band() - extracts a product from a BIL array        */
/* ----------------------------------------------------------- */
float *extract_band_double(double *in_buf, l2prodstr *p, int32 nbands) {
    static float32 *fbuf2 = NULL;

    int32 band = MIN(MAX(p->prod_ix, 0), nbands - 1);
    static int32 npix = 0;
    int32 ip;

    if (fbuf2 == NULL) {
        npix = p->dim[1];
        if ((fbuf2 = calloc(npix, sizeof (float32))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    } else if (npix != p->dim[1]) {
        npix = p->dim[1];
        free(fbuf2);
        // this is put here for l3gen which varies the number of pixels on each line

        if ((fbuf2 = (float32 *) calloc(npix, sizeof (float32)))
                == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }

    for (ip = 0; ip < npix && ip < numPixels; ip++) {
        fbuf2[ip] = (float)in_buf[ip * nbands + band];
    }

    return (fbuf2);
}

float *extract_band_3d(float *in_buf, float *out_buf) {
    float *out_ptr = out_buf;
    for (int pix = 0; pix < numPixels; pix++) {
        float *in_ptr = in_buf + pix * numBands;
        for (int band_3d = 0; band_3d < input->nwavelengths_3d; band_3d++) {
            int band = input->wavelength_3d_index[band_3d];
            *out_ptr = in_ptr[band];
            out_ptr++;
        }
    }
    return out_buf;
}

float *extract_band_3d_double(double *in_buf, float *out_buf) {
    float *out_ptr = out_buf;
    for (int pix = 0; pix < numPixels; pix++) {
        double *in_ptr = in_buf + pix * numBands;
        for (int band_3d = 0; band_3d < input->nwavelengths_3d; band_3d++) {
            int band = input->wavelength_3d_index[band_3d];
            *out_ptr = (float)in_ptr[band];
            out_ptr++;
        }
    }
    return out_buf;
}

/**
 * loop through a float array applying a multiplier
 *
 * @param in input array
 * @param out output array
 * @param count size of array
 * @param multiplier factor to apply
 */
void applyMultiplier(float* in, float* out, int count, float multiplier) {
    int i;
    for (i = 0; i < count; i++) {
        if (in[i] == BAD_FLT)
            out[i] = BAD_FLT;
        else
            out[i] = in[i] * multiplier;
    }
}

/* ----------------------------------------------------------- */
/* prodgen() - returns pointer the the requested product       */

/* ----------------------------------------------------------- */
VOIDP prodgen(l2prodstr *p, l2str *l2rec) {
    static int firstCall = 1;

    static float32 *fbuf = NULL;

    VOIDP pbuf = NULL;

    if (firstCall) {
        firstCall = 0;
        numPixels = l2rec->l1rec->npix;
        numScans = l2rec->l1rec->l1file->nscan;
        numBands = l2rec->l1rec->l1file->nbands;
        numBandsIR = l2rec->l1rec->l1file->nbandsir;
        if ((fbuf = (float32 *) calloc(numBands*numPixels,sizeof(float32))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    } else if (l2rec->l1rec->npix != numPixels) {

        // this is put here for l3gen which varies the number of pixels on each line
        if (l2rec->l1rec->npix > numPixels) {

            free(fbuf);
            if ((fbuf = (float32 *) calloc(numBands*l2rec->l1rec->npix, sizeof (float32)))
                    == NULL) {
                fprintf(stderr,
                        "-E- %s line %d: Unable to allocate buffer space.\n",
                        __FILE__, __LINE__);
                exit(1);
            }

        }
        numPixels = l2rec->l1rec->npix;
    }

    l1str* l1rec = l2rec->l1rec;
    switch (p->cat_ix) {

        //
        // Band-dependent, precomputed products
        //
    case CAT_nLw:
        if (p->rank == 3)
            applyMultiplier(extract_band_3d(l2rec->nLw, fbuf), fbuf, numPixels * input->nwavelengths_3d, 10.0);
        else
            applyMultiplier(extract_band(l2rec->nLw, p, numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
    case CAT_Lw:
        if (p->rank == 3)
            applyMultiplier(extract_band_3d(l2rec->Lw, fbuf), fbuf, numPixels * input->nwavelengths_3d, 10.0);
        else
            applyMultiplier(extract_band(l2rec->Lw, p, numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
    case CAT_Rrs:
        pbuf = p->rank == 3 ? extract_band_3d(l2rec->Rrs, fbuf) : extract_band(l2rec->Rrs, p, numBands);
        break;
    case CAT_Rrs_unc:
        if(!input->proc_uncertainty) {
            printf("-E- Can not produce product Rrs_unc without setting proc_uncertainty\n");
            exit(EXIT_FAILURE);
        }
        pbuf = p->rank == 3 ? extract_band_3d(l2rec->Rrs_unc, fbuf) : extract_band(l2rec->Rrs_unc, p, numBands);
        break;
    case CAT_Taua:
        pbuf = p->rank == 3 ? extract_band_3d(l2rec->taua, fbuf) : extract_band(l2rec->taua, p, numBands);
        break;
    case CAT_Lr:
        if (p->rank == 3)
            applyMultiplier(extract_band_3d(l1rec->Lr, fbuf), fbuf, numPixels * input->nwavelengths_3d, 10.0);
        else
            applyMultiplier(extract_band(l1rec->Lr, p, numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
    case CAT_L_q:
        if (p->rank == 3)
            applyMultiplier(extract_band_3d(l1rec->L_q, fbuf), fbuf, numPixels * input->nwavelengths_3d, 10.0);
        else
            applyMultiplier(extract_band(l1rec->L_q, p, numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
    case CAT_L_u:
        if (p->rank == 3)
            applyMultiplier(extract_band_3d(l1rec->L_u, fbuf), fbuf, numPixels * input->nwavelengths_3d, 10.0);
        else
            applyMultiplier(extract_band(l1rec->L_u, p, numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
    case CAT_polcor:
        pbuf = p->rank == 3 ? extract_band_3d(l1rec->polcor, fbuf) : extract_band(l1rec->polcor, p, numBands);
        break;
    case CAT_La:
        if (p->rank == 3)
            applyMultiplier(extract_band_3d(l2rec->La, fbuf), fbuf, numPixels * input->nwavelengths_3d, 10.0);
        else
            applyMultiplier(extract_band(l2rec->La, p, numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
    case CAT_TLg:
        if (p->rank == 3)
            applyMultiplier(extract_band_3d(l1rec->TLg, fbuf), fbuf, numPixels * input->nwavelengths_3d, 10.0);
        else
            applyMultiplier(extract_band(l1rec->TLg, p, numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
    case CAT_tLf:
        if (p->rank == 3)
            applyMultiplier(extract_band_3d(l1rec->tLf, fbuf), fbuf, numPixels * input->nwavelengths_3d, 10.0);
        else
            applyMultiplier(extract_band(l1rec->tLf, p, numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
    case CAT_brdf:
        pbuf = p->rank == 3 ? extract_band_3d(l2rec->brdf, fbuf) : extract_band(l2rec->brdf, p, numBands);
        break;
    case CAT_Lt:
        if (p->rank == 3) {
            applyMultiplier(extract_band_3d(l1rec->Lt, fbuf), fbuf, numPixels * input->nwavelengths_3d, 10.0);
            pbuf = fbuf;
        } else {
            if (p->prod_ix < numBands) {
                applyMultiplier(extract_band(l2rec->l1rec->Lt, p, numBands), fbuf, numPixels, 10.0);
                pbuf = fbuf;
            } else {
                // first subtract the num of visible bands
                p->prod_ix -= numBands;
                applyMultiplier(extract_band(l1rec->Ltir, p, NBANDSIR), fbuf, numPixels, 10.0);
                pbuf = fbuf;
                p->prod_ix += numBands;
            }
        }
        break;
    case CAT_BT:
        // first subtract the num of visible bands
        p->prod_ix -= numBands;
        pbuf = extract_band(l1rec->Bt, p, NBANDSIR);
        p->prod_ix += numBands;
        break;
    case CAT_rhos:
        pbuf = p->rank == 3 ? extract_band_3d(l1rec->rhos, fbuf) : extract_band(l1rec->rhos, p, numBands);
        break;
    case CAT_nw:
        pbuf = extract_band(l1rec->sw_n, p, numBands);
        break;
    case CAT_aw:
        pbuf = extract_band(l1rec->sw_a, p, numBands);
        break;
    case CAT_bbw:
        pbuf = extract_band(l1rec->sw_bb, p, numBands);
        break;
    case CAT_bbws:
        get_bbws(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_a:
        pbuf = p->rank == 3 ? extract_band_3d(l2rec->a, fbuf) : extract_band(l2rec->a, p, numBands);
        break;
    case CAT_bb:
        pbuf = p->rank == 3 ? extract_band_3d(l2rec->bb, fbuf) : extract_band(l2rec->bb, p, numBands);
        break;
    case CAT_t_sol:
        pbuf = p->rank == 3 ? extract_band_3d(l1rec->t_sol, fbuf) : extract_band(l1rec->t_sol, p, numBands);
        break;
    case CAT_t_sen:
        pbuf = p->rank == 3 ? extract_band_3d(l1rec->t_sen, fbuf) : extract_band(l1rec->t_sen, p, numBands);
        break;
    case CAT_tg_sol:
        pbuf = p->rank == 3 ? extract_band_3d_double(l1rec->tg_sol, fbuf) : extract_band_double(l1rec->tg_sol, p, numBands);
        break;
    case CAT_tg_sen:
        pbuf = p->rank == 3 ? extract_band_3d_double(l1rec->tg_sen, fbuf) : extract_band_double(l1rec->tg_sen, p, numBands);
        break;
    case CAT_t_h2o:
       pbuf = p->rank == 3 ? extract_band_3d(l1rec->t_h2o, fbuf) : extract_band(l1rec->t_h2o, p, numBands);
        break;
    case CAT_t_o2:
        pbuf = p->rank == 3 ? extract_band_3d(l1rec->t_o2, fbuf) : extract_band(l1rec->t_o2, p, numBands);
        break;
    case CAT_dpol:
        pbuf = p->rank == 3 ? extract_band_3d(l1rec->dpol, fbuf) : extract_band(l1rec->dpol, p, numBands);
        break;
    case CAT_BT_39:
        pbuf = extract_band(l1rec->Bt, p, NBANDSIR);
        break;
    case CAT_BT_40:
        pbuf = extract_band(l1rec->Bt, p, NBANDSIR);
        break;
    case CAT_BT_11:
        pbuf = extract_band(l1rec->Bt, p, NBANDSIR);
        break;
    case CAT_BT_12:
        pbuf = extract_band(l1rec->Bt, p, NBANDSIR);
        break;
        //
        // Band independent or dependent, precomputed products
        // for the view angles and related values
        //
    case CAT_senz:
    case CAT_solz:
    case CAT_sola:
    case CAT_sena:
    case CAT_relaz:
    case CAT_scattang:
        if ((p->prod_ix >= 0) && (l1rec->geom_per_band == NULL)) {
            fprintf(stderr,
                    "-E- %s, %d: Geometry-dependent view angle information",
                    __FILE__, __LINE__);
            fprintf(stderr, " is unavailable for product catalog ID %d.\n",
                    p->cat_ix);
            fprintf(stderr, "or, geometry-dependent view angle information was specifically not requested\n");
            exit(1);
        }
        switch (p->cat_ix) {
        case CAT_senz:
            if (p->prod_ix < 0) {
                pbuf = (VOIDP) l1rec->senz;
            } else {
                pbuf = extract_band(l1rec->geom_per_band->senz, p, numBands);
            }
            break;
        case CAT_solz:
            if (p->prod_ix < 0) {
                pbuf = (VOIDP) l1rec->solz;
            } else {
                pbuf = extract_band(l1rec->geom_per_band->solz, p, numBands);
            }
            break;
        case CAT_sena:
            if (p->prod_ix < 0) {
                pbuf = (VOIDP) l1rec->sena;
            } else {
                pbuf = extract_band(l1rec->geom_per_band->sena, p, numBands);
            }
            break;
        case CAT_sola:
            if (p->prod_ix < 0) {
                pbuf = (VOIDP) l1rec->sola;
            } else {
                pbuf = extract_band(l1rec->geom_per_band->sola, p, numBands);
            }
            break;
        case CAT_relaz:
            if (p->prod_ix < 0) {
                pbuf = (VOIDP) l1rec->delphi;
            } else {
                pbuf = extract_band(l1rec->geom_per_band->delphi, p, numBands);
            }
            break;
        case CAT_scattang:
            if (p->prod_ix < 0) {
                pbuf = (VOIDP) l1rec->scattang;
            } else {
                pbuf = extract_band(l1rec->geom_per_band->scattang, p, numBands);
            }
            break;
        }
        break;

        //
        // Band independent, precomputed products
        //
    case CAT_epsilon:
        pbuf = (VOIDP) l2rec->eps;
        break;
    case CAT_chisqr_mbac:
        pbuf = (VOIDP) l2rec->chi2;
        break;
    case CAT_alpha:
        pbuf = (VOIDP) l1rec->alpha;
        break;
    case CAT_ozone:
        pbuf = (VOIDP) l1rec->oz;
        break;
    case CAT_no2_tropo:
        pbuf = (VOIDP) l1rec->no2_tropo;
        break;
    case CAT_no2_strat:
        pbuf = (VOIDP) l1rec->no2_strat;
        break;
    case CAT_no2_frac:
        pbuf = (VOIDP) l1rec->no2_frac;
        break;
    case CAT_windspeed:
        pbuf = (VOIDP) l1rec->ws;
        break;
    case CAT_windangle:
        pbuf = (VOIDP) l1rec->wd;
        break;
    case CAT_zwind:
        pbuf = (VOIDP) l1rec->zw;
        break;
    case CAT_mwind:
        pbuf = (VOIDP) l1rec->mw;
        break;
    case CAT_pressure:
        pbuf = (VOIDP) l1rec->pr;
        break;
    case CAT_water_vapor:
        pbuf = (VOIDP) l1rec->wv;
        break;
    case CAT_humidity:
        pbuf = (VOIDP) l1rec->rh;
        break;
    case CAT_sfc_pressure:
        pbuf = (VOIDP) l1rec->sfcp;
        break;
    case CAT_sfc_humidity:
        pbuf = (VOIDP) l1rec->sfcrh;
        break;
    case CAT_sfc_temp:
        pbuf = (VOIDP) l1rec->sfct;
        break;
    case CAT_T_prof:
    case CAT_RH_prof:
    case CAT_HGT_prof:
    case CAT_Q_prof:
    case CAT_O3_prof:
        if( p->rank == 3 ) {
          switch (p->cat_ix) {
          case CAT_T_prof:
            pbuf = (VOIDP) l1rec->anc_add->prof_temp;
            break;
          case CAT_RH_prof:
            pbuf = (VOIDP) l1rec->anc_add->prof_rh;
            break;
          case CAT_HGT_prof:
            pbuf = (VOIDP) l1rec->anc_add->prof_height;
            break;
          case CAT_Q_prof:
            pbuf = (VOIDP) l1rec->anc_add->prof_q;
            break;
          case CAT_O3_prof:
            pbuf = (VOIDP) l1rec->anc_add->prof_o3;
            break;
          }
        } else {
          fprintf( stderr,
          "-I- %s %d: Single-level profile output not implemented yet\n",
          __FILE__, __LINE__ );
        }
        break;
    case CAT_height:
        pbuf = (VOIDP) l1rec->height;
        break;
    case CAT_sstref:
        pbuf = (VOIDP) l1rec->sstref;
        break;
    case CAT_sssref:
        pbuf = (VOIDP) l1rec->sssref;
        break;
    case CAT_glint_coef:
        pbuf = (VOIDP) l1rec->glint_coef;
        break;
    case CAT_cloud_albedo:
        pbuf = (VOIDP) l1rec->cloud_albedo;
        break;
    case CAT_rho_cirrus:
        pbuf = (VOIDP) l1rec->rho_cirrus;
        break;
    case CAT_aerindex:
        pbuf = (VOIDP) l2rec->aerindex;
        break;
    case CAT_aer_ratio:
        if (p->prod_ix == 1)
            pbuf = (VOIDP) l2rec->aerratio;
        else
            pbuf = (VOIDP) l2rec->aerratio2;
        break;

        //
        // Integer, precomputed  products
        //
    case CAT_l2_flags:
        pbuf = (VOIDP) l2rec->l1rec->flags;
        break;
    case CAT_num_iter:
        pbuf = (VOIDP) l2rec->num_iter;
        break;
    case CAT_slot:
        pbuf = (VOIDP) l2rec->l1rec->slot;
        break;
    case CAT_aer_model:
        switch (p->prod_ix) {
        case 1:
            pbuf = l2rec->aermodmin;
            break;
        case 2:
            pbuf = l2rec->aermodmax;
            break;
        case 3:
            pbuf = l2rec->aermodmin2;
            break;
        case 4:
            pbuf = l2rec->aermodmax2;
            break;
        }
        break;

        //
        // 1-dimensional, precomputed products
        //
    case CAT_fsol:
        pbuf = (VOIDP) & l1rec->fsol;
        break;
    case CAT_pixnum:
        pbuf = (VOIDP) l1rec->pixnum;
        break;
    case CAT_detnum:
        pbuf = (VOIDP) & l1rec->detnum;
        break;
    case CAT_mside:
        pbuf = (VOIDP) & l1rec->mside;
        break;

        //
        // Derived products
        //
    case CAT_angstrom:
        get_angstrom(l2rec, p->prod_ix, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_ms_epsilon:
        get_ms_epsilon(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_Es:
        get_es(l2rec, p->prod_ix, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_rhot:
        get_toa_refl(l2rec, p->prod_ix, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_rhom:
        get_rho_mumm(l2rec, -1, p->prod_ix, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_depth_class:
        get_depth_classification(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_fqy:
        get_fqy(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_fsat:
        get_fsat(l2rec, fbuf);
        applyMultiplier(fbuf, fbuf, numPixels, 10.0);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_ipar:
        get_ipar(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_par_scalar:
        get_par_scalar(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_mu_0:
        get_mu_cosine(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_par_below_surface:
        get_par_below_surface(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_BSi:
        get_bsi(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_Kd_mueller:
    case CAT_Kd_532:
    case CAT_Kd_obpg:
    case CAT_Kd_lee:
    case CAT_Kd_unc_lee:
    case CAT_Kd_morel:
    case CAT_Kd_KD2:
    case CAT_Kd_unc_KD2:
    case CAT_KPAR_morel:
    case CAT_KPAR_lee:
    case CAT_Zhl_morel:
    case CAT_Kd_jamet:
    case CAT_Kd_rhos:
    case CAT_nKd_lin:
        get_Kd(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_Zphotic_lee:
    case CAT_Zeu_morel:
    case CAT_Zsd_morel:
    case CAT_Zsd_gbr:
        get_photic_depth(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_tindx_morel:
        tindx_morel(l2rec, -1, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_tindx_shi:
        tindx_shi(l2rec, -1, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_iCDOM_morel:
    case CAT_pCDOM_morel:
    case CAT_adg_morel:
    case CAT_chl_morel:
    case CAT_chl_cdomcorr_morel:
        get_cdom_morel(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_owt:
    case CAT_owtn:
    case CAT_chl_owterr:
        optical_water_type(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_owtd:
        optical_water_type(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_ndvi:
    case CAT_evi:
    case CAT_evi2:
    case CAT_evi3:
        get_ndvi_evi(l1rec, p->cat_ix, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_smoke:
        get_smoke(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_par:
        get_par(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_par2:
        get_par2(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_taucld:
        get_taucld(l2rec,fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_clfr:
        get_clfr(l2rec,fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_ipar2:
        get_ipar2(l2rec,fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_ipar_below_surface:
        get_ipar_below_surface(l2rec,fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_ipar_scalar:
        get_ipar_scalar(l2rec,fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_chl_oc2:
    case CAT_chl_oc3:
    case CAT_chl_oc3c:
    case CAT_chl_oc4:
    case CAT_chl_hu:
    case CAT_chl_oci:
    case CAT_chl_oci2:
    case CAT_chl_cdr:
    case CAT_chl_abi:
        get_chl(l2rec, p->cat_ix, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_chl_unc:
        pbuf=l2rec->chl_unc;
        break;
    case CAT_chl_mgiop:
    case CAT_bbp_mgiop:
    case CAT_adg_mgiop:
    case CAT_aph_mgiop:
    case CAT_npix_mgiop:
    case CAT_crat_mgiop:
    case CAT_fitpar_mgiop:
        get_mgiop(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_a_pml:
    case CAT_aph_pml:
    case CAT_adg_pml:
    case CAT_bb_pml:
    case CAT_bbp_pml:
        get_pml(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_a_qaa:
    case CAT_aph_qaa:
    case CAT_adg_qaa:
    case CAT_bb_qaa:
    case CAT_bbp_qaa:
    case CAT_mod_rrs_qaa:
        get_qaa(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_flags_qaa:
        pbuf = (VOIDP) get_flags_qaa(l2rec);
        break;
    case CAT_chl_carder:
    case CAT_a_carder:
    case CAT_bb_carder:
    case CAT_aph_carder:
    case CAT_adg_carder:
    case CAT_bbp_carder:
        get_carder(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_chl_carder_emp:
        chl_carder_empirical(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_flags_carder:
        pbuf = (VOIDP) get_flags_carder(l2rec);
        break;
    case CAT_flags_giop:
        pbuf = (VOIDP) get_flags_giop(l2rec);
        break;
    case CAT_iter_giop:
        pbuf = (VOIDP) get_iter_giop(l2rec);
        break;
    case CAT_chl_giop:
    case CAT_a_giop:
    case CAT_bb_giop:
    case CAT_aph_giop:
    case CAT_adg_giop:
    case CAT_bbp_giop:
    case CAT_chl_unc_giop:
    case CAT_a_unc_giop:
    case CAT_bb_unc_giop:
    case CAT_aph_unc_giop:
    case CAT_adg_unc_giop:
    case CAT_bbp_unc_giop:
    case CAT_aphs_giop:
    case CAT_adgs_giop:
    case CAT_bbps_giop:
    case CAT_rrsdiff_giop:
    case CAT_mRrs_giop:
    case CAT_chisqr_giop:
    case CAT_fitpar_giop:
    case CAT_acdom_giop:
    case CAT_anap_giop:
    case CAT_bbph_giop:
    case CAT_bbnap_giop:
    case CAT_acdom_unc_giop:
    case CAT_anap_unc_giop:
    case CAT_bbph_unc_giop:
    case CAT_bbnap_unc_giop:
    case CAT_opt_siop_giop:
        get_giop(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_iter_gsm:
        pbuf = (VOIDP) get_iter_gsm(l2rec);
        break;
    case CAT_chl_gsm:
    case CAT_a_gsm:
    case CAT_bb_gsm:
    case CAT_aph_gsm:
    case CAT_adg_gsm:
    case CAT_bbp_gsm:
        get_gsm(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_a_las:
    case CAT_b_las:
    case CAT_c_las:
    case CAT_bb_las:
    case CAT_bbp_las:
    case CAT_bbps_las:
        get_las(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_a_niwa:
    case CAT_bb_niwa:
        get_niwa(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_flags_niwa:
        pbuf = (VOIDP) get_flags_niwa(l2rec);
        break;
//    case CAT_chl_soa:
//    case CAT_adg_soa:
//    case CAT_bbp_soa:
//    case CAT_pcentcdm_soa:
//    case CAT_w0_soa:
//    case CAT_v_soa:
//        get_soa(l2rec, p->cat_ix, fbuf);
//        pbuf = (VOIDP) fbuf;
//        break;
//    case CAT_chl_sma:
//    case CAT_adg_sma:
//    case CAT_bbp_sma:
//    case CAT_w0_sma:
//    case CAT_dom_sma:
//        get_sma(l2rec, p->cat_ix, fbuf);
//        pbuf = (VOIDP) fbuf;
//        break;
    case CAT_poc_stramski_443:
    case CAT_poc_stramski_490:
    case CAT_poc_stramski_hybrid:
        get_poc(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_ag_412_mlrc:
    case CAT_Sg_275_295_mlrc:
    case CAT_Sg_300_600_mlrc:
        cdom_mannino(l2rec, p->cat_ix, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_calcite:
    case CAT_calcite_2b:
    case CAT_calcite_3b:
    case CAT_calcite_ci2:
    case CAT_calcite_ci748:
    case CAT_calcite_ci869:
        calcite(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_Rrs_vc:
        if (input->band_shift_opt == 1)
            bioOptBandShift(l2rec, p, fbuf);
        else
            virtual_constellation(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_chl_vc:
        virtual_constellation(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_sst:
        pbuf = (VOIDP) get_sst(l2rec);
        break;
    case CAT_sst4:
        pbuf = (VOIDP) get_sst4(l2rec);
        break;
    case CAT_sst3:
        pbuf = (VOIDP) get_sst_triple(l2rec);
        break;
    case CAT_dsdi:
        pbuf = (VOIDP) get_sst_dust_correction(l2rec);
        break;
    case CAT_bias_sst:
        pbuf = (VOIDP) get_bias_sst(l2rec);
        break;
    case CAT_bias_sst4:
        pbuf = (VOIDP) get_bias_sst4(l2rec);
        break;
    case CAT_bias_sst3:
        pbuf = (VOIDP) get_bias_sst_triple(l2rec);
        break;
    case CAT_bias_mean_sst:
        pbuf = (VOIDP) get_bias_mean_sst(l2rec);
        break;
    case CAT_bias_mean_sst4:
        pbuf = (VOIDP) get_bias_mean_sst4(l2rec);
        break;
    case CAT_bias_mean_sst3:
        pbuf = (VOIDP) get_bias_mean_sst_triple(l2rec);
        break;
    case CAT_counts_sst:
        pbuf = (VOIDP) get_counts_sst(l2rec);
        break;
    case CAT_counts_sst4:
        pbuf = (VOIDP) get_counts_sst4(l2rec);
        break;
    case CAT_counts_sst3:
        pbuf = (VOIDP) get_counts_sst_triple(l2rec);
        break;
    case CAT_stdv_sst:
        pbuf = (VOIDP) get_stdv_sst(l2rec);
        break;
    case CAT_stdv_sst4:
        pbuf = (VOIDP) get_stdv_sst4(l2rec);
        break;
    case CAT_stdv_sst3:
        pbuf = (VOIDP) get_stdv_sst_triple(l2rec);
        break;
    case CAT_flags_sst:
        pbuf = (VOIDP) get_flags_sst(l2rec);
        break;
    case CAT_flags_sst4:
        pbuf = (VOIDP) get_flags_sst4(l2rec);
        break;
    case CAT_flags_sst3:
        pbuf = (VOIDP) get_flags_sst_triple(l2rec);
        break;
    case CAT_qual_sst:
        pbuf = (VOIDP) get_qual_sst(l2rec);
        break;
    case CAT_qual_sst4:
        pbuf = (VOIDP) get_qual_sst4(l2rec);
        break;
    case CAT_qual_sst3:
        pbuf = (VOIDP) get_qual_sst_triple(l2rec);
        break;
    case CAT_vLt:
    case CAT_vtLw:
    case CAT_vLw:
    case CAT_vnLw:
        vcal(l2rec, p, fbuf);
        applyMultiplier(fbuf, fbuf, numPixels, 10.0);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_vgain:
    case CAT_vbsat:
    case CAT_vbtgt:
        vcal(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_ice_frac:
        get_ice_frac(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_class_ward_owmc:
        pbuf = get_class_ward_owmc(l2rec);
        break;
    case CAT_class_k_owmc:
        pbuf = get_class_k_owmc(l2rec);
        break;
    case CAT_class_34k_w_owmc:
        pbuf = get_class_34k_w_owmc(l2rec);
        break;
    case CAT_a_swim:
    case CAT_bb_swim:
    case CAT_adg_swim:
    case CAT_aph_swim:
    case CAT_bbp_swim:
        get_swim(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_elev:
        pbuf = (VOIDP) l1rec->dem;
        break;
    case CAT_microplankton_hirata:
    case CAT_diatoms_hirata:
    case CAT_greenalgae_hirata:
    case CAT_picoplankton_hirata:
    case CAT_prokaryotes_hirata:
    case CAT_prochlorococcus_hirata:
    case CAT_dinoflagellates_hirata:
    case CAT_nanoplankton_hirata:
    case CAT_picoeukaryotes_hirata:
    case CAT_prymnesiophytes_hirata:
        get_pft_hirata(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_microplankton_uitz:
    case CAT_nanoplankton_uitz:
    case CAT_picoplankton_uitz:
        get_pft_uitz(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_npp_vgpm:
    case CAT_npp_eppley:
    case CAT_npp_cbpm2:
    case CAT_npp_mld:
    case CAT_npp_zno3:
    case CAT_npp_bbp:
    case CAT_npp_par:
    case CAT_npp_cafe:
        get_npp(l2rec, p->cat_ix, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_CI_stumpf:
    case CAT_CI_cyano:
    case CAT_CI_noncyano:
    case CAT_MCI_stumpf:
        get_habs_ci(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_MPH_chl:
        get_habs_mph(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_flags_habs_mph:
        pbuf = (VOIDP) get_flags_habs_mph(l2rec);
        break;
    case CAT_flags_habs:
        pbuf = (VOIDP) get_flags_habs(l2rec);
        break;

    case CAT_microplankton_abundanceksm:
    case CAT_nanoplankton_abundanceksm:
    case CAT_picoplankton_abundanceksm:
    case CAT_microplankton_volumeksm:
    case CAT_nanoplankton_volumeksm:
    case CAT_picoplankton_volumeksm:
    case CAT_microplankton_ratioksm:
    case CAT_nanoplankton_ratioksm:
    case CAT_picoplankton_ratioksm:
        get_psd_ksm(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;

    case CAT_nitrate:
    	get_nitrate(l2rec,p,fbuf);
        pbuf = (VOIDP) fbuf;
    	break;
    case CAT_avw:
        get_avw(l2rec,fbuf);
        pbuf=(VOIDP) fbuf;
        break;
    case CAT_Rrs_brightness:
        get_Rrs_brightness(l2rec,fbuf);
        pbuf=(VOIDP) fbuf;
        break;
    case CAT_lambda_max:
        get_lambda_max(l2rec,fbuf);
        pbuf=(VOIDP) fbuf;
        break;
    case CAT_Cphyt:
        get_Cphyt(l2rec,fbuf);
        pbuf=(VOIDP) fbuf;
        break;
    case CAT_Cphyt_unc:
        get_Cphyt_unc(l2rec,fbuf);
        pbuf=(VOIDP) fbuf;
        break;
    case CAT_prochlorococcus:
    case CAT_synechococcus:
    case CAT_autotrophic_picoeukaryotes:
        get_Cpicophyt(l2rec,p,fbuf);
        pbuf=(VOIDP)fbuf;
        break;
    
    case CAT_iparb:
    case CAT_parb:
    	get_bpar(l2rec,p,fbuf);
        pbuf = (VOIDP) fbuf;
    	break;
#ifdef BUILD_CLOUD        
   /* Chimaera cloud products */
    case CAT_CER_2100:
    case CAT_CER_2200:
    case CAT_CER_1600:
    case CAT_COT_2100:
    case CAT_COT_2200:
    case CAT_COT_1600:
    case CAT_CER_1621:
    case CAT_COT_1621:
    case CAT_CWP_2100:
    case CAT_CWP_2200:
    case CAT_CWP_1621:
    case CAT_CWP_1600:
    case CAT_Cld_Top_Refl_650:
    case CAT_Cld_Top_Refl_860:
    case CAT_Cld_Top_Refl_1200:
    case CAT_Cld_Top_Refl_1600:
    case CAT_Cld_Top_Refl_2100:
    case CAT_Cld_Top_Refl_2200:

    case CAT_COT_fail_2100:
    case CAT_COT_fail_2200:
    case CAT_COT_fail_1600:
    case CAT_COT_fail_1621:
    case CAT_CER_fail_2100:
    case CAT_CER_fail_2200:
    case CAT_CER_fail_1600:
    case CAT_CER_fail_1621:
    case CAT_CMP_fail_pct_2100:
    case CAT_CMP_fail_pct_2200:
    case CAT_CMP_fail_pct_1600:
    case CAT_CMP_fail_pct_1621:
    case CAT_refl_loc_1600:
    case CAT_refl_loc_2100:
    case CAT_refl_loc_2200:
    case CAT_refl_loc_1621:
/*   initial albedo */
    case CAT_Surface_Albedo_650:
    case CAT_Surface_Albedo_860:
    case CAT_Surface_Albedo_1200:
    case CAT_Surface_Albedo_1600:
    case CAT_Surface_Albedo_2100:
    case CAT_Surface_Albedo_2200:
        get_cmp(l2rec, p->cat_ix, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_Cld_Sfc_Type:
    case CAT_Cld_Phase_2100:
    case CAT_Cld_Phase_2200:
    case CAT_Cld_Non_Abs_Band:
    case CAT_Cld_Phase_1600:
    case CAT_Cld_Phase_1621:
    case CAT_Cld_water_cloud:
    case CAT_Cld_ice_cloud:
        pbuf = (VOIDP)get_cmp_byt( l2rec, p->cat_ix );
        break;
#endif    
/*  WDR cloud top height parameters, we just get the values from the record */
    case CAT_cth_alb_init:
        pbuf = (VOIDP) l1rec->cld_dat->cth_alb_init;
        break;
    case CAT_cth_alb_unc_init:
        pbuf = (VOIDP) l1rec->cld_dat->cth_alb_unc_init;
        break;
/*  WDR cloud top height parameters to be computed */
    case CAT_Cld_p:
    case CAT_Cld_t:
    case CAT_Cld_h:
    case CAT_cth_cod:
    case CAT_cth_alb:
    case CAT_cth_phase:
    case CAT_cth_cod_all:
    case CAT_cth_alb_all:
    case CAT_cth_lcod_all:
    case CAT_cth_cth_all:
    case CAT_cth_ctp_all:
    case CAT_cth_cost_all:
    case CAT_cth_acost_all:
    case CAT_cth_iter_all:
    case CAT_cth_cth_raw_all:
    case CAT_cth_ctt_all:
    case CAT_cth_lcod:
    case CAT_cth_cost:
    case CAT_cth_acost:
    case CAT_cth_iter:
    case CAT_cth_cth_raw:
    case CAT_cth_dcod_all:
    case CAT_cth_dlcod_all:
    case CAT_cth_dcth_all:
    case CAT_cth_dctp_all:
    case CAT_cth_dctt_all:
    case CAT_cth_dalb_all:
    case CAT_cth_dcod:
    case CAT_cth_dlcod:
    case CAT_cth_dcth:
    case CAT_cth_dctp:
    case CAT_cth_dctt:
    case CAT_cth_dalb:
        pbuf = (VOIDP) get_ctht_lin( l2rec, p->cat_ix );
        break;
    case CAT_ndii:
        get_ndii( l1rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_ndwi:
        get_ndwi( l1rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_ndsi:
        get_ndsi( l1rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    case CAT_cci :
        get_cci ( l1rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
    // Hyperspectral vegetation indices
    case CAT_pri: // Photochemical Reflectance Index
    case CAT_cire: // Chlorophyll Index Red Edge
    case CAT_car: // Chlorophyll-Carotenoid Index
    case CAT_ari: // Anthocyanin Rreflectance Index
        get_hyper_vi ( l1rec, p->cat_ix, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
	
    //SDP pigments products
    case CAT_tchl_sdp:
    case CAT_zea_sdp:
    case CAT_dvchla_sdp:
    case CAT_butfuco_sdp:
    case CAT_hexfuco_sdp:
    case CAT_allo_sdp:
    case CAT_mvchlb_sdp: 
    case CAT_neo_sdp:
    case CAT_viola_sdp:
    case CAT_fuco_sdp:
    case CAT_chlc12_sdp:
    case CAT_chlc3_sdp:
    case CAT_perid_sdp:
    case CAT_flags_sdp:
    case CAT_mRrs_sdp:
    case CAT_mRrs_diff_sdp:
    case CAT_d2Rrs_diff_sdp: 
        get_sdp(l2rec,p,fbuf);
        pbuf=(VOIDP) fbuf;
        break;

    case CAT_tchl_gpig:
    case CAT_chlc12_gpig:
    case CAT_tchlb_gpig:
    case CAT_ppc_gpig:
    case CAT_flags_gpig:
        get_gpig(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;  

    default:
        fprintf(stderr, "-E- %s Line %d: Unknown product catalogue ID %d.\n",
                __FILE__, __LINE__, p->cat_ix);
        exit(1);
    }

    return (pbuf);
}

