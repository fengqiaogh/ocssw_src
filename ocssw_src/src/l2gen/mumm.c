/* =================================================================== */
/* MUMM module - turbid water correction for Gordon & Wang atmosphere  */
/*                                                                     */
/* Reference:                                                          */
/*                                                                     */
/* Ruddick, K., F.Ovidio & M.Rijkeboer (2000). Atmospheric correction  */
/* of SeaWiFS imagery for turbid coastal and inland waters,            */
/* Applied Optics, 39(6), pp897-912.                                   */
/*                                                                     */
/* Written By:                                                         */
/*                                                                     */
/* B. Franz, NASA/OBPG, 03 November 2006 (based on implementation from */
/* http://www.mumm.ac.be/OceanColour/Products/Software/index.php)      */
/* W. Robinson, SAIC, 26 Apr 2017  add band dependence for solz        */
/*                                                                     */
/* =================================================================== */

#include "l12_proto.h"

/* ------------------------------------------------------------------- */
/* get_rho_mumm(): compute quasi-surface reflectance preferred by MUMM.*/

/* ------------------------------------------------------------------- */

void get_rho_mumm(l2str *l2rec, int32_t ipix, int32_t iw, float *rhom) {
    int32_t ip, ip1, ip2, ipb;
    float Ltemp;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;
    int32_t nbands = l1file->nbands;
    float solz;

    if (ipix < 0) {
        ip1 = 0;
        ip2 = l1rec->npix - 1;
    } else {
        ip1 = ipix;
        ip2 = ipix;
    }

    for (ip = ip1; ip <= ip2; ip++) {

        if (l1rec->mask[ip])
            *rhom++ = 0.0;

        else {

            ipb = ip * nbands + iw;
            solz = (l1rec->geom_per_band == NULL) ? l1rec->solz[ip] :
                    l1rec->geom_per_band->solz[ipb];

            Ltemp = (l2rec->l1rec->Lt[ipb]
                    / l1rec->tg_sol[ipb] / l1rec->tg_sen[ipb]
                    / l1rec->polcor[ipb]
                    - l1rec->tLf[ipb]
                    - l1rec->Lr[ipb]
                    ) / l1rec->t_o2[ipb]
                    - l1rec->TLg[ipb];

            *rhom++ = OEL_PI * Ltemp / l1rec->Fo[iw] / cos(solz / OEL_RADEG);

        }
    }
}


/* ------------------------------------------------------------------- */
/* get_rhown_mumm(): compute the normalized water-leaving reflectance  */
/* contribution in the NIR using MUMM algorithm.                       */

/* ------------------------------------------------------------------- */

void get_rhown_mumm(l2str *l2rec, int32_t ip, int32_t nir_s, int32_t nir_l, float rhown[]) {
    float alpha = input->mumm_alpha;
    float gamma = input->mumm_gamma;
    float epsilon = input->mumm_epsilon;

    float rhom_s, rhom_l;
    float rhoa_s, rhoa_l;

    get_rho_mumm(l2rec, ip, nir_s, &rhom_s);
    get_rho_mumm(l2rec, ip, nir_l, &rhom_l);

    rhoa_l = MAX(MIN((rhom_l * gamma * alpha - rhom_s) / (gamma * alpha - epsilon), rhom_l), 0.0);
    rhoa_s = MIN(epsilon*rhoa_l, rhom_s);

    rhown[nir_s] = rhom_s - rhoa_s;
    rhown[nir_l] = rhom_l - rhoa_l;

    return;
}

