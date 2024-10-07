/* -------------------------------------------------------------------- */
/* l1subpix() - sub-samples and/or crops a level-1 record in place      */
/*                                                                      */
/* Returns 0 on succss, 1 on error.                                     */
/*                                                                      */
/* Notes: Only record fields which are filled by the L1 read routines   */
/*        must be updated here.  All other fields will be filled later  */
/*        using the sub-sampled geometry and radiances.                 */
/*                                                                      */
/* Written By: B. A. Franz                                              */
/* W. Robinson, SAIC 15Feb2017  add geom_per_band shifting              */
/* -------------------------------------------------------------------- */

#include <string.h>
#include "l1.h"

int l1subpix(filehandle *l1file, l1str *l1rec) {
    int32_t sp = l1file->spix;
    int32_t ep = l1file->epix;
    int32_t dp = l1_input->dpixl;
    int32_t nbands = l1file->nbands;

    if ((sp == 0) && (ep == l1file->npix - 1) && dp == 1)
        return (0);

    if (sp > ep || dp < 1 || dp > (ep - sp + 1)) {
        fprintf(stderr,
                "-E- %s Line %d: subpixel specification error (sp=%d,ep=%d,dp=%d).\n",
                __FILE__, __LINE__, sp, ep, dp);
        return (1);
    }

    if (dp == 1) {

        int32_t length;

        l1rec->npix = ep - sp + 1;

        length = l1rec->npix * sizeof (float);

        memmove(l1rec->pixnum, &l1rec->pixnum[sp], l1rec->npix * sizeof (int32_t));
        memmove(l1rec->slot, &l1rec->slot [sp], l1rec->npix * sizeof (unsigned char));
        memmove(l1rec->nobs, &l1rec->nobs [sp], l1rec->npix * sizeof (int32_t));

        memmove(l1rec->lon, &l1rec->lon [sp], length);
        memmove(l1rec->lat, &l1rec->lat [sp], length);
        memmove(l1rec->solz, &l1rec->solz[sp], length);
        memmove(l1rec->sola, &l1rec->sola[sp], length);
        memmove(l1rec->senz, &l1rec->senz[sp], length);
        memmove(l1rec->sena, &l1rec->sena[sp], length);

        memmove(l1rec->alpha, &l1rec->alpha[sp], length);
        memmove(l1rec->height, &l1rec->height[sp], length);
        memmove(l1rec->flags, &l1rec->flags[sp], l1rec->npix * sizeof (int32_t));
        memmove(l1rec->hilt, &l1rec->hilt[sp], l1rec->npix * sizeof (char));
        memmove(l1rec->stlight, &l1rec->stlight[sp], l1rec->npix * sizeof (char));
        memmove(l1rec->navfail, &l1rec->navfail[sp], l1rec->npix * sizeof (char));
        memmove(l1rec->navwarn, &l1rec->navwarn[sp], l1rec->npix * sizeof (char));

        memmove(l1rec->Lt, &l1rec->Lt [sp * nbands], length * nbands);

        memmove(l1rec->sw_n, &l1rec->sw_n [sp * nbands], length * nbands);
        memmove(l1rec->sw_a, &l1rec->sw_a [sp * nbands], length * nbands);
        memmove(l1rec->sw_bb, &l1rec->sw_bb [sp * nbands], length * nbands);
        memmove(l1rec->sw_a_avg, &l1rec->sw_a_avg [sp * nbands], length * nbands);
        memmove(l1rec->sw_bb_avg, &l1rec->sw_bb_avg [sp * nbands], length * nbands);

        memmove(l1rec->Ltir, &l1rec->Ltir [sp * NBANDSIR], length * NBANDSIR);
        memmove(l1rec->Bt, &l1rec->Bt [sp * NBANDSIR], length * NBANDSIR);

        memmove(l1rec->rho_cirrus, &l1rec->rho_cirrus[sp], length);

        if (l1rec->geom_per_band != NULL) {
            memmove(l1rec->geom_per_band->senz,
                    &l1rec->geom_per_band->senz[ sp * nbands ], length * nbands);
            memmove(l1rec->geom_per_band->sena,
                    &l1rec->geom_per_band->sena[ sp * nbands ], length * nbands);
            memmove(l1rec->geom_per_band->solz,
                    &l1rec->geom_per_band->solz[ sp * nbands ], length * nbands);
            memmove(l1rec->geom_per_band->sola,
                    &l1rec->geom_per_band->sola[ sp * nbands ], length * nbands);
        }

        /* MERIS */
        memmove(l1rec->pixdet, &l1rec->pixdet[sp], l1rec->npix * sizeof (int32_t));
        memmove(l1rec->radcor, &l1rec->radcor[sp * nbands], length * nbands);

    } else {

        int32_t i;

        l1rec->npix = (ep - sp) / dp + 1;

        for (i = 0; i < l1rec->npix; i++) {

            l1rec->pixnum[i] = l1rec->pixnum[i * dp + sp];
            l1rec->slot[i] = l1rec->slot[i * dp + sp];
            l1rec->nobs[i] = l1rec->nobs[i * dp + sp];

            l1rec->lon [i] = l1rec->lon [i * dp + sp];
            l1rec->lat [i] = l1rec->lat [i * dp + sp];
            l1rec->solz[i] = l1rec->solz[i * dp + sp];
            l1rec->sola[i] = l1rec->sola[i * dp + sp];
            l1rec->senz[i] = l1rec->senz[i * dp + sp];
            l1rec->sena[i] = l1rec->sena[i * dp + sp];

            l1rec->alpha[i] = l1rec->alpha[i * dp + sp];
            l1rec->height[i] = l1rec->height[i * dp + sp];
            l1rec->flags[i] = l1rec->flags[i * dp + sp];

            l1rec->hilt [i] = l1rec->hilt [i * dp + sp];
            l1rec->stlight[i] = l1rec->stlight[i * dp + sp];
            l1rec->navfail[i] = l1rec->navfail[i * dp + sp];
            l1rec->navwarn[i] = l1rec->navwarn[i * dp + sp];

            l1rec->rho_cirrus[i] = l1rec->rho_cirrus[i * dp + sp];

            memmove(&l1rec->Lt [i * nbands],
                    &l1rec->Lt [(i * dp + sp) * nbands], nbands * sizeof (float));

            memmove(&l1rec->sw_n [i * nbands],
                    &l1rec->sw_n [(i * dp + sp) * nbands], nbands * sizeof (float));
            memmove(&l1rec->sw_a [i * nbands],
                    &l1rec->sw_a [(i * dp + sp) * nbands], nbands * sizeof (float));
            memmove(&l1rec->sw_bb[i * nbands],
                    &l1rec->sw_bb[(i * dp + sp) * nbands], nbands * sizeof (float));
            memmove(&l1rec->sw_a_avg [i * nbands],
                    &l1rec->sw_a_avg [(i * dp + sp) * nbands], nbands * sizeof (float));
            memmove(&l1rec->sw_bb_avg[i * nbands],
                    &l1rec->sw_bb_avg[(i * dp + sp) * nbands], nbands * sizeof (float));

            memmove(&l1rec->Ltir [i * NBANDSIR],
                    &l1rec->Ltir [(i * dp + sp) * NBANDSIR], NBANDSIR * sizeof (float));

            memmove(&l1rec->Bt [i * NBANDSIR],
                    &l1rec->Bt [(i * dp + sp) * NBANDSIR], NBANDSIR * sizeof (float));

            if (l1rec->geom_per_band != NULL) {
                memmove(&l1rec->geom_per_band->senz[ i * nbands ],
                        &l1rec->geom_per_band->senz[ (i * dp + sp) * nbands ],
                        nbands * sizeof (float));
                memmove(&l1rec->geom_per_band->sena[ i * nbands ],
                        &l1rec->geom_per_band->sena[ (i * dp + sp) * nbands ],
                        nbands * sizeof (float));
                memmove(&l1rec->geom_per_band->solz[ i * nbands ],
                        &l1rec->geom_per_band->solz[ (i * dp + sp) * nbands ],
                        nbands * sizeof (float));
                memmove(&l1rec->geom_per_band->sola[ i * nbands ],
                        &l1rec->geom_per_band->sola[ (i * dp + sp) * nbands ],
                        nbands * sizeof (float));
            }

            /* MERIS */
            l1rec->pixdet[i] = l1rec->pixdet[i * dp + sp];
            memmove(&l1rec->radcor [i * nbands],
                    &l1rec->radcor [(i * dp + sp) * nbands], nbands * sizeof (float));

        }
    }

    return (0);
}

