/*
 *  file geom_per_band.c - generic operations needed for the 
 *  band-dependent geometry capability
 */
#include "l1.h"

int init_geom_per_band(l1str *l1rec)
/*-----------------------------------------------------------------------
 init_geom_per_band

 purpose: set up the band-dependent geometry storage
  
 Returns 0 if all checks are OK
 Parameters: (in calling order)
 Type              Name            I/O     Description
 ----              ----            ---     -----------
 l1str *           l1rec            I      record containing the band-
                                        dependent geometry fields and 
                                        the sizes in pixels, bands

 Modification history:
 Programmer        Date            Description of change
 ----------        ----            ---------------------
 Wayne Robinson    2 Dec 2016     Original development

 -----------------------------------------------------------------------*/ {
    int32_t nbands, npix;

    nbands = l1rec->l1file->nbands;
    npix = l1rec->npix;
    l1rec->geom_per_band = (geom_struc *) malloc(sizeof (geom_struc));
    l1rec->geom_per_band->senz =
            (float *) malloc(npix * nbands * sizeof (float));
    l1rec->geom_per_band->sena =
            (float *) malloc(npix * nbands * sizeof (float));
    l1rec->geom_per_band->csenz =
            (float *) malloc(npix * nbands * sizeof (float));
    l1rec->geom_per_band->solz =
            (float *) malloc(npix * nbands * sizeof (float));
    l1rec->geom_per_band->sola =
            (float *) malloc(npix * nbands * sizeof (float));
    l1rec->geom_per_band->csolz =
            (float *) malloc(npix * nbands * sizeof (float));
    l1rec->geom_per_band->delphi =
            (float *) malloc(npix * nbands * sizeof (float));
    l1rec->geom_per_band->scattang =
            (float *) malloc(npix * nbands * sizeof (float));
    return 0;
}

int geom_per_band_deriv(l1str *l1rec)
/*-----------------------------------------------------------------------
 geom_per_band_deriv

 purpose: derive the other band-dependent geometry quantities

 Returns 0 if all checks are OK
 Parameters: (in calling order)
 Type              Name            I/O     Description
 ----              ----            ---     -----------
 l1str *           l1rec            I      record containing the band-
                                        dependent geometry fields and
                                        the sizes in pixels, bands

 Modification history:
 Programmer        Date            Description of change
 ----------        ----            ---------------------
 Wayne Robinson   10 Jan 2017      Original development

 -----------------------------------------------------------------------*/ {
    static double radeg = RADEG;
    double temp;
    int32_t nbands, npix, ip, ib, ipb;
    geom_struc *geom = l1rec->geom_per_band;

    nbands = l1rec->l1file->nbands;
    npix = l1rec->npix;
    /*  This just creates the csolz, csola, delphi, and scattang for the
        band-dependent geometry values */
    for (ib = 0; ib < nbands; ib++) {
        for (ip = 0; ip < npix; ip++) {
            ipb = ip * nbands + ib;

            /* relative azimuth */
            geom->delphi[ipb] = geom->sena[ipb] - 180.0 - geom->sola[ipb];
            if (geom->delphi[ipb] < -180.) geom->delphi[ipb] += 360.0;
            else if (geom->delphi[ipb] > 180.0)
                geom->delphi[ipb] -= 360.0;

            /* frequently used trig relations */
            geom->csolz[ipb] = cos(geom->solz[ipb] / radeg);
            geom->csenz[ipb] = cos(geom->senz[ipb] / radeg);

            /* Scattering angle */
            temp = sqrt((1.0 - geom->csenz[ipb] * geom->csenz[ipb])*
                    (1.0 - geom->csolz[ipb] * geom->csolz[ipb]))
                    * cos(geom->delphi[ipb] / radeg);
            geom->scattang[ipb] = acos(MAX(-geom->csenz[ipb] * geom->csolz[ipb] +
                    temp, -1.0)) * radeg;
        }
    }
    return 0;
}
