#include "viirs_sim_sdr.h"
#include <math.h>
#include <stdlib.h>

int scan_cvt(in_rec_struc *in_rec, out_rec_struc *out_rec)
/*-----------------------------------------------------------------------------
   scan_cvt

   purpose:  make proper output record from the input record 
     based on the scan formats for both

   Returns type: int - 0 if good

   Parameters (in calling order):
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      in_rec_struc *    in_rec          I/O  controls for input record reading
      out_rec_struc *   out_rec         I/O  controls for output file writing

   Modification history:

   W. Robinson, SAIC  29 Apr 2010  Original development
   W. Robinson, SAIC  18 Nov 2010  create the QF1_VIIRSMBANDSDR value in 
      out_rec array qual1_m from saturation and gain info in in_rec

----------------------------------------------------------------------------*/ {
    static int cvt_mode = -1; /* conversion mode for this run: -1 not started,
                                0 - no conversion, 1 - make aggregated,
                                2 - make unaggregated */
    static int npixin, npixout;
    static float *cvt_lat, *cvt_lon, *cvt_senz, *cvt_sena, *cvt_solz, *cvt_sola,
            *cvt_bnd_lt[MAX_BND];
    int ibnd, ilin, olin, opix, nxfr, nxfrb, ipix_st, ag_rng, num, iag, ipix;
    float theta1, theta2, theta, phi1, phi2, phi, sum;
    static unsigned char *cvt_bnd_q[MAX_BND];
    unsigned char uc_val, ua_gain_vals[] = {4, 0};
    unsigned char ua_sat_vals[] = {0, 8};
    int v_unsat, v_sat, ag_sat_vals[] = {8, 0, 8, 4};

    /*
     *  if in_rec->lat is null, this is signal to for end, so free local storage
     */
    if (in_rec->lat == NULL) {
        free(cvt_lat);
        free(cvt_lon);
        free(cvt_senz);
        free(cvt_sena);
        free(cvt_solz);
        free(cvt_sola);
        for (ibnd = 0; ibnd < out_rec->nbnd; ibnd++) {
            free(cvt_bnd_lt[ibnd]);
            free(cvt_bnd_q[ibnd]);
        }
        return 0;
    }
    /*
     *  initialization phase
     */
    if (cvt_mode == -1) {
        /*
         *  determine the conversion mode based on the input scan format and 
         *  requested output scan format change
         */
        if ((in_rec->scn_fmt == 1) && (out_rec->scn_fmt == 0))
            cvt_mode = 1;
        else if ((in_rec->scn_fmt == 2) && (out_rec->scn_fmt == 0))
            cvt_mode = 1;
        else if ((in_rec->scn_fmt == 2) && (out_rec->scn_fmt == 1))
            cvt_mode = 2;
        else {
            cvt_mode = 0;
        }
        if (cvt_mode != 0) {
            /*
             *  A conversion is needed for all scans.  allocate space for
             *  the converted data and re-assign it to the out_rec locations
             *  (from default done in init_sdr with in_rec)
             */
            npixin = in_rec->npix;
            npixout = out_rec->npix;
            if ((cvt_lat = (float *)
                    malloc(out_rec->npix * out_rec->ndet_scan * sizeof (float)))
                    == NULL) {
                printf("%s, %d: Unable to allocate storage for lat conversion\n",
                        __FILE__, __LINE__);
                return 1;
            }
            out_rec->lat = cvt_lat;
            if ((cvt_lon = (float *)
                    malloc(out_rec->npix * out_rec->ndet_scan * sizeof (float)))
                    == NULL) {
                printf("%s, %d: Unable to allocate storage for lon conversion\n",
                        __FILE__, __LINE__);
                return 1;
            }
            out_rec->lon = cvt_lon;
            if ((cvt_senz = (float *)
                    malloc(out_rec->npix * out_rec->ndet_scan * sizeof (float)))
                    == NULL) {
                printf("%s, %d: Unable to allocate storage for senz conversion\n",
                        __FILE__, __LINE__);
                return 1;
            }
            out_rec->senz = cvt_senz;
            if ((cvt_sena = (float *)
                    malloc(out_rec->npix * out_rec->ndet_scan * sizeof (float)))
                    == NULL) {
                printf("%s, %d: Unable to allocate storage for sena conversion\n",
                        __FILE__, __LINE__);
                return 1;
            }
            out_rec->sena = cvt_sena;
            if ((cvt_solz = (float *)
                    malloc(out_rec->npix * out_rec->ndet_scan * sizeof (float)))
                    == NULL) {
                printf("%s, %d: Unable to allocate storage for solz conversion\n",
                        __FILE__, __LINE__);
                return 1;
            }
            out_rec->solz = cvt_solz;
            if ((cvt_sola = (float *)
                    malloc(out_rec->npix * out_rec->ndet_scan * sizeof (float)))
                    == NULL) {
                printf("%s, %d: Unable to allocate storage for sola conversion\n",
                        __FILE__, __LINE__);
                return 1;
            }
            out_rec->sola = cvt_sola;
            for (ibnd = 0; ibnd < out_rec->nbnd; ibnd++) {
                if ((cvt_bnd_lt[ibnd] = (float *)
                        malloc(out_rec->npix * out_rec->ndet_scan * sizeof (float)))
                        == NULL) {
                    printf(
                            "%s, %d: Unable to allocate storage for bnd_lt[%d] conversion\n",
                            __FILE__, __LINE__, ibnd);
                    return 1;
                }
                out_rec->bnd_lt[ibnd] = cvt_bnd_lt[ibnd];
                /*
                 */
                if ((cvt_bnd_q[ibnd] = (unsigned char *)
                        malloc(out_rec->npix * out_rec->ndet_scan *
                        sizeof (unsigned char))) == NULL) {
                    printf(
                            "%s, %d: Unable to allocate storage for bnd_q[%d] conversion\n",
                            __FILE__, __LINE__, ibnd);
                    return 1;
                }
                out_rec->bnd_q[ibnd] = cvt_bnd_q[ibnd];
            }
        }
    }
    /*
     *  work on different conversion modes
     */
    if (cvt_mode == 0) {
        /*
         *  if the data is aggregated, set up qual1_m from the dn_sat * 4
         *  (most likely, dn_sat always 0 in this case, but if not, we'll
         *  assume the 3 state format for the standard saturation field)
         */
        if (out_rec->scn_fmt == 0) {
            for (ibnd = 0; ibnd < in_rec->nbnd; ibnd++) {
                for (ilin = 0; ilin < in_rec->ndet_scan; ilin++) {
                    for (ipix = 0; ipix < in_rec->npix; ipix++) {
                        *(out_rec->qual1_m[ibnd] + ipix + in_rec->npix * ilin) =
                                *(in_rec->dn_sat[ibnd] + ipix + in_rec->npix * ilin) * 4;
                    }
                }
            }
        }            /*
    *  if the data is un-aggregated, fold the gain and sat into qual1_m
    *  bit 2 = gain 0 hi, 1 lo, bit 3 = sat 0 not, 1 saturated
    */
        else {
            for (ibnd = 0; ibnd < in_rec->nbnd; ibnd++) {
                for (ilin = 0; ilin < in_rec->ndet_scan; ilin++) {
                    for (ipix = 0; ipix < in_rec->npix; ipix++) {
                        uc_val = ua_gain_vals[ *(in_rec->gain_bit[ibnd] + ipix +
                                in_rec->npix * ilin) ];
                        uc_val = uc_val | ua_sat_vals[ (int) *(in_rec->dn_sat[ibnd] + ipix +
                                in_rec->npix * ilin) ];
                        *(out_rec->qual1_m[ibnd] + ipix + in_rec->npix * ilin) = uc_val;
                    }
                }
            }
        }
        return 0;
    } else if (cvt_mode == 2) {
        /*
         *  mode 2 will go from margin to unaggregated, so move to proper storage
         */
        nxfr = npixout * sizeof (float);
        nxfrb = npixout * sizeof (unsigned char);
        for (olin = 0, ilin = in_rec->margin[0]; olin < NDET; olin++, ilin++) {
            memcpy((void *) (cvt_lat + npixout * olin),
                    (void *) (in_rec->lat + npixin * ilin + in_rec->margin[1]), nxfr);
            memcpy((void *) (cvt_lon + npixout * olin),
                    (void *) (in_rec->lon + npixin * ilin + in_rec->margin[1]), nxfr);
            memcpy((void *) (cvt_sena + npixout * olin),
                    (void *) (in_rec->sena + npixin * ilin + in_rec->margin[1]), nxfr);
            memcpy((void *) (cvt_senz + npixout * olin),
                    (void *) (in_rec->senz + npixin * ilin + in_rec->margin[1]), nxfr);
            memcpy((void *) (cvt_sola + npixout * olin),
                    (void *) (in_rec->sola + npixin * ilin + in_rec->margin[1]), nxfr);
            memcpy((void *) (cvt_solz + npixout * olin),
                    (void *) (in_rec->solz + npixin * ilin + in_rec->margin[1]), nxfr);
            for (ibnd = 0; ibnd < out_rec->nbnd; ibnd++) {
                memcpy((void *) (cvt_bnd_lt[ibnd] + npixout * olin),
                        (void *) (in_rec->bnd_lt[ibnd] + npixin * ilin + in_rec->margin[1]),
                        nxfr);
                memcpy((void *) (cvt_bnd_q[ibnd] + npixout * olin),
                        (void *) (in_rec->bnd_q[ibnd] + npixin * ilin + in_rec->margin[1]),
                        nxfrb);
                /*
                 *  handle the conversion and transfer of saturation and gain
                 */
                for (ipix = in_rec->margin[1], opix = 0;
                        ipix < in_rec->npix - in_rec->margin[1]; ipix++, opix++) {
                    uc_val = ua_gain_vals[ *(in_rec->gain_bit[ibnd] + ipix +
                            in_rec->npix * ilin) ];
                    uc_val = uc_val | ua_sat_vals[ (int) *(in_rec->dn_sat[ibnd] + ipix +
                            in_rec->npix * ilin) ];
                    *(out_rec->qual1_m[ibnd] + opix + npixout * olin) = uc_val;
                }
            }
        }
    } else if (cvt_mode == 1) {
        /*
         *  mode 1 will involve aggregating and possibly trimming off the margins
         *  Method for aggregating depends on the agg zone:
         *                     ZONE (scan location)
         *  PARM          1:1 agg (edge)     2:1 agg      3:1 agg (center)
         *  -------
         *  lat, lon,         Just           'Angle          Use center
         *  view angles     Transfer        Average'           angle
         * 
         *  Lt data          values         Average          Average
         *                                  2 values         3 values
         */
        for (olin = 0, ilin = in_rec->margin[0]; olin < out_rec->ndet_scan;
                olin++, ilin++) {
            for (opix = 0; opix < npixout; opix++) {
                if ((opix >= 0) && (opix <= 639)) {
                    ipix_st = opix + in_rec->margin[1];
                    ag_rng = 1;
                } else if ((opix >= 640) && (opix <= 1007)) {
                    ipix_st = 640 + in_rec->margin[1] + (opix - 640) * 2;
                    ag_rng = 2;
                } else if ((opix >= 1008) && (opix <= 2191)) {
                    ipix_st = 1376 + in_rec->margin[1] + (opix - 1008) * 3;
                    ag_rng = 3;
                } else if ((opix >= 2192) && (opix <= 2559)) {
                    ipix_st = 4928 + in_rec->margin[1] + (opix - 2192) * 2;
                    ag_rng = 2;
                } else if ((opix >= 2560) && (opix <= 3199)) {
                    ipix_st = 5664 + in_rec->margin[1] + opix - 2560;
                    ag_rng = 1;
                }
                /*
                 *  Handle each agg range for the parameters
                 */
                switch (ag_rng) {
                case 1:
                    /*
                     *  at 1:1 zones, just transfer (to proper locations)
                     */
                    *(cvt_lat + npixout * olin + opix) =
                            *(in_rec->lat + ipix_st + npixin * ilin);
                    *(cvt_lon + npixout * olin + opix) =
                            *(in_rec->lon + ipix_st + npixin * ilin);
                    *(cvt_senz + npixout * olin + opix) =
                            *(in_rec->senz + ipix_st + npixin * ilin);
                    *(cvt_sena + npixout * olin + opix) =
                            *(in_rec->sena + ipix_st + npixin * ilin);
                    *(cvt_solz + npixout * olin + opix) =
                            *(in_rec->solz + ipix_st + npixin * ilin);
                    *(cvt_sola + npixout * olin + opix) =
                            *(in_rec->sola + ipix_st + npixin * ilin);
                    for (ibnd = 0; ibnd < out_rec->nbnd; ibnd++) {
                        *(cvt_bnd_lt[ibnd] + npixout * olin + opix) =
                                *(in_rec->bnd_lt[ibnd] + ipix_st + npixin * ilin);
                        *(cvt_bnd_q[ibnd] + npixout * olin + opix) =
                                *(in_rec->bnd_q[ibnd] + ipix_st + npixin * ilin);
                        /*
                         *  saturation aggregation - the unagg saturation is either 0 or 1
                         *  it xlates to 0 or 2 (all sat) and at bit location it 
                         *  translates to 0 or 8 (x4)
                         */
                        *(out_rec->qual1_m[ibnd] + opix + npixout * olin) =
                                *(in_rec->dn_sat[ibnd] + ipix_st + npixin * ilin);
                    }
                    break;
                case 2:
                    /*
                     *  at 2:1 zones, average the radiances and treat lat, lon and 
                     *  view angles as vectors to find the average of
                     *
                     *  Lt
                     */
                    for (ibnd = 0; ibnd < out_rec->nbnd; ibnd++) {
                        sum = 0.;
                        num = 0;
                        v_unsat = 0;
                        v_sat = 0;
                        for (iag = 0; iag < ag_rng; iag++) {
                            if (
                                    (*(in_rec->bnd_lt[ibnd] + ilin * npixin + ipix_st + iag) > 0)
                                    && (*(in_rec->bnd_q[ibnd] + ilin * npixin + ipix_st + iag) == 0)) {
                                sum +=
                                        *(in_rec->bnd_lt[ibnd] + ilin * npixin + ipix_st + iag);
                                num++;
                                /*
                                 * see if saturated or unsaturated or both - the v_unsat and
                                 * v_sat work together with ag_sat_vals to set proper value
                                 */
                                if (*(in_rec->dn_sat[ibnd] + ilin + npixin + ipix_st + iag)
                                        == 0) v_unsat = 1;
                                if (*(in_rec->dn_sat[ibnd] + ilin + npixin + ipix_st + iag)
                                        == 1) v_sat = 2;
                            }
                        }
                        if (num > 0) {
                            *(cvt_bnd_lt[ibnd] + npixout * olin + opix) = sum / num;
                            *(cvt_bnd_q[ibnd] + npixout * olin + opix) = 0;
                        } else {
                            *(cvt_bnd_lt[ibnd] + npixout * olin + opix) = -32767.;
                            *(cvt_bnd_q[ibnd] + npixout * olin + opix) = 2;
                        }
                        /*  for saturation setting  */
                        *(out_rec->qual1_m[ibnd] + npixout * olin + opix) =
                                ag_sat_vals[ v_sat + v_unsat ];
                    }
                    /*
                     *  lat, lon, and view angles
                     */
                    theta1 = *(in_rec->lat + ilin * npixin + ipix_st);
                    theta2 = *(in_rec->lat + ilin * npixin + ipix_st + 1);
                    phi1 = *(in_rec->lon + ilin * npixin + ipix_st);
                    phi2 = *(in_rec->lon + ilin * npixin + ipix_st + 1);
                    ang_avg(theta1, phi1, theta2, phi2, &theta, &phi);
                    *(cvt_lat + npixout * olin + opix) = theta;
                    *(cvt_lon + npixout * olin + opix) = phi;

                    theta1 = *(in_rec->senz + ilin * npixin + ipix_st);
                    theta2 = *(in_rec->senz + ilin * npixin + ipix_st + 1);
                    phi1 = *(in_rec->sena + ilin * npixin + ipix_st);
                    phi2 = *(in_rec->sena + ilin * npixin + ipix_st + 1);
                    ang_avg(theta1, phi1, theta2, phi2, &theta, &phi);
                    *(cvt_senz + npixout * olin + opix) = theta;
                    *(cvt_sena + npixout * olin + opix) = phi;

                    theta1 = *(in_rec->solz + ilin * npixin + ipix_st);
                    theta2 = *(in_rec->solz + ilin * npixin + ipix_st + 1);
                    phi1 = *(in_rec->sola + ilin * npixin + ipix_st);
                    phi2 = *(in_rec->sola + ilin * npixin + ipix_st + 1);
                    ang_avg(theta1, phi1, theta2, phi2, &theta, &phi);
                    *(cvt_solz + npixout * olin + opix) = theta;
                    *(cvt_sola + npixout * olin + opix) = phi;
                    break;
                case 3:
                    /*
                     *  at 3:1 zones, average 3 lt and take center angle
                     */
                    for (ibnd = 0; ibnd < out_rec->nbnd; ibnd++) {
                        sum = 0.;
                        num = 0;
                        v_unsat = 0;
                        v_sat = 0;
                        for (iag = 0; iag < ag_rng; iag++) {
                            if (
                                    (*(in_rec->bnd_lt[ibnd] + ilin * npixin + ipix_st + iag) > 0)
                                    && (*(in_rec->bnd_q[ibnd] + ilin * npixin + ipix_st + iag) == 0)) {
                                sum +=
                                        *(in_rec->bnd_lt[ibnd] + ilin * npixin + ipix_st + iag);
                                num++;
                                /*
                                 * see if saturated or unsaturated or both - the v_unsat and
                                 * v_sat work together with ag_sat_vals to set proper value
                                 */
                                if (*(in_rec->dn_sat[ibnd] + ilin + npixin + ipix_st + iag)
                                        == 0) v_unsat = 1;
                                if (*(in_rec->dn_sat[ibnd] + ilin + npixin + ipix_st + iag)
                                        == 1) v_sat = 2;
                            }
                        }
                        if (num > 0) {
                            *(cvt_bnd_lt[ibnd] + npixout * olin + opix) = sum / num;
                            *(cvt_bnd_q[ibnd] + npixout * olin + opix) = 0;
                        } else {
                            *(cvt_bnd_lt[ibnd] + npixout * olin + opix) = -32767.;
                            *(cvt_bnd_q[ibnd] + npixout * olin + opix) = 2;
                        }
                        /*  for saturation setting  */
                        *(out_rec->qual1_m[ibnd] + npixout * olin + opix) =
                                ag_sat_vals[ v_sat + v_unsat ];
                    }

                    *(cvt_lat + npixout * olin + opix) =
                            *(in_rec->lat + ilin * npixin + ipix_st + 1);
                    *(cvt_lon + npixout * olin + opix) =
                            *(in_rec->lon + ilin * npixin + ipix_st + 1);
                    *(cvt_senz + npixout * olin + opix) =
                            *(in_rec->senz + ilin * npixin + ipix_st + 1);
                    *(cvt_sena + npixout * olin + opix) =
                            *(in_rec->sena + ilin * npixin + ipix_st + 1);
                    *(cvt_solz + npixout * olin + opix) =
                            *(in_rec->solz + ilin * npixin + ipix_st + 1);
                    *(cvt_sola + npixout * olin + opix) =
                            *(in_rec->sola + ilin * npixin + ipix_st + 1);
                    break;
                }
            }
        }
    }
    return 0;
}

int ang_avg(float theta1, float phi1, float theta2, float phi2,
        float *theta, float *phi)
/*-----------------------------------------------------------------------------
   ang_avg

   purpose:  get resultant vector that is the average of 2 vectors

   Returns type: int - 0 if good

   Parameters (in calling order):
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float             theta1           I      first zenith or lat angle in 
                                                  degrees
      float             phi1             I      first azimuth or lon angle in 
                                                  degrees
      float             theta2           I      second zenith or lat angle in
                                                  degrees
      float             phi2             I      second azimuth or lon angle in
                                                  degrees
      float *           theta            O      resultant zenith or lat angle in
                                                  degrees
      float *           phi              O      result azimuth or lon angle in
                                                  degrees

   Modification history:

   W. Robinson, SAIC  4 May 2010  Original development

----------------------------------------------------------------------------*/ {
    float v1[2], v2[2], v[2];
    float rad2deg = 180. / M_PI;
    /*
     *  convert to just the x, y parts of a vector - this assumes all are
     *  unit length
     *** don't know why but sinf, cosf... do not work, so use less compatible
            with floats sin, cos (gcc?)
     */
    v1[0] = sin(theta1 / rad2deg) * cos(phi1 / rad2deg);
    v1[1] = sin(theta1 / rad2deg) * sin(phi1 / rad2deg);

    v2[0] = sin(theta2 / rad2deg) * cos(phi2 / rad2deg);
    v2[1] = sin(theta2 / rad2deg) * sin(phi2 / rad2deg);

    /* average the vectors and return to unit length */
    v[0] = (v1[0] + v2[0]) / 2.;
    v[1] = (v1[1] + v2[1]) / 2.;

    /* go back to angles */
    *theta = asin(sqrt(v[0] * v[0] + v[1] * v[1])) * rad2deg;
    *phi = atan2(v[1], v[0]) * rad2deg;

    return 0;
}
