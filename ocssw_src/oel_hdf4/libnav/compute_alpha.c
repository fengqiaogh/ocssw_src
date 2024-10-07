#include <math.h>

/* ----------------------------------------------------------------------------- */
/* compute_alpha() - computes polarization correction frame rotation angle       */
/*                                                                               */
/* algorithm provided by: F. S. Patt, SAIC, February 2003.                       */
/* ----------------------------------------------------------------------------- */
void compute_alpha(float lon[],
                   float lat[],
                   float senz[],
                   float sena[],
                   double mnorm[3],
                   int npix,
                   float alpha[]) {
    static double radeg = 57.29577951;

    int ip, i;

    double slon, clon;
    double slat, clat;
    double szen, czen;
    double sazi, cazi;

    double e[3], n[3], v[3], r[3], s[3];
    double sdvcr, vdr, sdr, sdv;

    /* invert mirror normal */
    for (i = 0; i < 3; i++)
        r[i] = -mnorm[i];

    for (ip = 0; ip < npix; ip++) {
        slon = sin(lon[ip] / radeg);
        clon = cos(lon[ip] / radeg);
        slat = sin(lat[ip] / radeg);
        clat = cos(lat[ip] / radeg);
        szen = sin(senz[ip] / radeg);
        czen = cos(senz[ip] / radeg);
        sazi = sin(sena[ip] / radeg);
        cazi = cos(sena[ip] / radeg);

        /* pixel coordinate system (north, east, vertical) in ECR */
        e[0] = -slon;
        e[1] = clon;
        e[2] = 0.0;

        n[0] = -slat * clon;
        n[1] = -slat * slon;
        n[2] = clat;

        v[0] = clat * clon;
        v[1] = clat * slon;
        v[2] = slat;

        /* sensor view vector in ECR */
        for (i = 0; i < 3; i++)
            s[i] = e[i] * szen * sazi + n[i] * szen * cazi + v[i] * czen;

        /* compute rotation angle (alpha) from pixel normal (v) to mirror */
        /* normal (r) about sensor view vector (s)  (Wertz, p. 757)       */
        sdvcr = s[0] * (v[1] * r[2] - v[2] * r[1])
            + s[1] * (v[2] * r[0] - v[0] * r[2])
            + s[2] * (v[0] * r[1] - v[1] * r[0]);

        vdr = v[0] * r[0] + v[1] * r[1] + v[2] * r[2];
        sdr = s[0] * r[0] + s[1] * r[1] + s[2] * r[2];
        sdv = v[0] * s[0] + v[1] * s[1] + v[2] * s[2];

        /* negated to be consistent with Gordon et al. */
        alpha[ip] = -(radeg * atan2(sdvcr, (vdr - sdr * sdv)) - 90.0);
    }
}
