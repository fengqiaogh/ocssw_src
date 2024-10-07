/*---------------------------------------------------------------------*/
/* ipar.c -  functions to compute instantaneous PAR                    */
/*---------------------------------------------------------------------*/

#include "l12_proto.h"
#include "l2prod.h"

static float unitc = 119.625e8; /* conversion to einsteins/m^2/s */

#define PARW1 400
#define PARW2 700
#define PARWN (PARW2-PARW1)+1

/*---------------------------------------------------------------------*/
/* get_ipar - computes ipar for each pixel in L2 record                */

/*---------------------------------------------------------------------*/
void get_ipar(l2str *l2rec, float ipar[]) {
    static float badval = BAD_FLT;
    static int firstCall = 1;
    static int32_t nwave;
    static float *wave;
    static float F0vis[PARWN];
    static float *ta;

    float Ed0p;
    int32_t ip, ipb, iw, ib;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;

    if (firstCall) {

        firstCall = 0;
        if ((wave = (float *) calloc(l1file->nbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for ipar:get_ipar.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        for (iw = 0; iw < l1file->nbands; iw++) {
            wave[iw] = l1file->fwave[iw];
        }

        // instantaneous solar irradiance at 1-nm intervals
        for (iw = PARW1; iw <= PARW2; iw++) {
            ib = iw - PARW1;
            get_f0_thuillier_ext(iw, 1, &F0vis[ib]);
            F0vis[ib] *= l1rec->fsol;
        }

        nwave = MIN(windex(PARW2, wave, l1file->nbands) + 1, l1file->nbands);

        if ((ta = (float *) calloc(nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for ipar:get_ipar.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

    }

    for (ip = 0; ip < l1rec->npix; ip++) {

        ipar[ip] = 0.;

        if (!l1rec->mask[ip]) {

            for (iw = 0; iw < nwave; iw++) {

                ipb = ip * l1file->nbands + iw;

                // atmospheric transmittance (sun to surface)
                ta[iw] = l1rec->tg_sol[ipb] * l1rec->t_sol[ipb];
            }

            for (iw = PARW1; iw <= PARW2; iw++) {
                ib = iw - PARW1;
                Ed0p = F0vis[ib] * l1rec->csolz[ip] * linterp(wave, ta, nwave, (float) iw);
                ipar[ip] += ((float) iw) * Ed0p / unitc;
            }
            ipar[ip]*=1e6;
            if (!isfinite(ipar[ip])) {
                ipar[ip] = badval;
                l1rec->flags[ip] |= PRODFAIL;
            }

        } else {
            ipar[ip] = badval;
            l1rec->flags[ip] |= PRODFAIL;
        }
    }
}

/*
        Energy of one photon = hc/lamda.
        So the number of photons generating Ed (watts/m2/nm) is

        Ed / (hc/lamba) = Ed*lamda /(hc)

        = Ed*lamda/ (6.626e-34 J Sec. x 2.998e8 m/Sec()	
        = Ed*lamda/ (1.9865e-25 J m)	

        So using lamda in nm

        (J/Sec/m2/nm) * nm * 1m/1e9nm / (1.9865e-25 J m x 6.022e23 photons / Ein)

        = (Ein/Sec/m2/nm) / 1.19625e8.	
 */
