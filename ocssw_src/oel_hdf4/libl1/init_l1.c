
#include "l1.h"

#include <stdio.h>
#include <stdlib.h>

/* --------------------------------------------------------- */
/* init_l1() - initialize a Level-1 record                   */

/* --------------------------------------------------------- */
void init_l1(l1str *l1rec) {
    int32_t nbands = l1rec->l1file->nbands;
    int32_t nbir = NBANDSIR;
    int32_t npix = l1rec->npix = l1rec->l1file->npix;

    int32_t ip, ib, ipb;

    for (ip = 0; ip < npix; ip++) {
        l1rec->lon [ip] = BAD_FLT;
        l1rec->lat [ip] = BAD_FLT;
        l1rec->solz [ip] = BAD_FLT;
        l1rec->sola [ip] = BAD_FLT;
        l1rec->senz [ip] = BAD_FLT;
        l1rec->sena [ip] = BAD_FLT;
        l1rec->delphi [ip] = BAD_FLT;
        l1rec->csolz [ip] = BAD_FLT;
        l1rec->csenz [ip] = BAD_FLT;
        l1rec->alpha [ip] = BAD_FLT;
        l1rec->scattang [ip] = BAD_FLT;

        l1rec->ws [ip] = BAD_FLT;
        l1rec->wd [ip] = BAD_FLT;
        l1rec->mw [ip] = BAD_FLT;
        l1rec->zw [ip] = BAD_FLT;
        l1rec->pr [ip] = BAD_FLT;
        l1rec->oz [ip] = BAD_FLT;
        l1rec->wv [ip] = BAD_FLT;
        l1rec->rh [ip] = BAD_FLT;
        l1rec->sfcp [ip] = BAD_FLT;
        l1rec->sfcrh [ip] = BAD_FLT;
        l1rec->sfct [ip] = BAD_FLT;
        l1rec->icefr [ip] = BAD_FLT;
        l1rec->no2_tropo [ip] = BAD_FLT;
        l1rec->no2_strat [ip] = BAD_FLT;
        l1rec->no2_frac [ip] = BAD_FLT;
        l1rec->height [ip] = BAD_FLT;
        l1rec->glint_coef [ip] = BAD_FLT;
        l1rec->cloud_albedo[ip] = BAD_FLT;
        l1rec->aerindex [ip] = BAD_FLT;
        l1rec->sstref [ip] = BAD_FLT;
        l1rec->sssref [ip] = BAD_FLT;
        l1rec->rho_cirrus [ip] = BAD_FLT;

        l1rec->ssttype[ip] = 0;
        l1rec->pixnum [ip] = 0;
        l1rec->slot [ip] = 0;
        l1rec->flags [ip] = 0;
        l1rec->dem [ip] = 0;
        l1rec->ancqc [ip] = 0;
        l1rec->mask [ip] = 0;
        l1rec->hilt [ip] = 0;
        l1rec->cloud [ip] = 0;
        l1rec->glint [ip] = 0;
        l1rec->land [ip] = 0;
        l1rec->swater [ip] = 0;
        l1rec->ice [ip] = 0;
        l1rec->solzmax[ip] = 0;
        l1rec->senzmax[ip] = 0;
        l1rec->stlight[ip] = 0;
        l1rec->absaer [ip] = 0;
        l1rec->navfail[ip] = 0;
        l1rec->navwarn[ip] = 0;
        l1rec->filter [ip] = 0;
        l1rec->pixdet [ip] = 0;
        l1rec->nobs [ip] = 1;

        for (ib = 0; ib < nbands; ib++) {
            ipb = ip * nbands + ib;
            l1rec->Lt [ipb] = BAD_FLT;
            l1rec->t_h2o [ipb] = BAD_FLT;
            l1rec->t_o2 [ipb] = BAD_FLT;
            l1rec->tg_sol [ipb] = BAD_FLT;
            l1rec->tg_sen [ipb] = BAD_FLT;
            l1rec->t_sol [ipb] = BAD_FLT;
            l1rec->t_sen [ipb] = BAD_FLT;
            l1rec->rhof [ipb] = BAD_FLT;
            l1rec->tLf [ipb] = BAD_FLT;
            l1rec->Lr [ipb] = BAD_FLT;
            l1rec->L_q [ipb] = BAD_FLT;
            l1rec->L_u [ipb] = BAD_FLT;
            l1rec->polcor [ipb] = BAD_FLT;
            l1rec->dpol [ipb] = BAD_FLT;
            l1rec->TLg [ipb] = BAD_FLT;
            l1rec->rhos [ipb] = BAD_FLT;
            l1rec->sw_a [ipb] = BAD_FLT;
            l1rec->sw_bb [ipb] = BAD_FLT;
            l1rec->sw_a_avg [ipb] = BAD_FLT;
            l1rec->sw_bb_avg[ipb] = BAD_FLT;
            l1rec->sw_n [ipb] = BAD_FLT;
            l1rec->radcor [ipb] = 0.0;

        }

        for (ib = 0; ib < nbir; ib++) {
            ipb = ip * nbir + ib;
            l1rec->Ltir[ipb] = BAD_FLT;
            l1rec->Bt [ipb] = BAD_FLT;
        }

    } // ip

    l1rec->alt = BAD_FLT;
    l1rec->scn_fmt = 0;
    l1rec->is_l2 = 0;

    for (ib = 0; ib < nbands; ib++) {
        l1rec->Fo[ib] = BAD_FLT;
    }
}



