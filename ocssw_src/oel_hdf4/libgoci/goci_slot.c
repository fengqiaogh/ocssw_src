
#include "hdf5.h"
#include "hdf5_hl.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <proj.h>
#include "goci.h"

#define NFIELDS  (hsize_t)     28
#define NRECORDS (hsize_t)     128
#define NAV_GRP  "HDFEOS/POINTS/Navigation for GOCI/Data"
#define TABLE_NAME             "Navigation for GOCI"

/*  quick number sorter for qsort */
static int numcomp(const void *p1, const void *p2) {
    return *(char *) p1 - *(char *) p2;
}

int32_t goci_slot_init(hid_t file_id, hsize_t *dims, float *slot_rel_time,
        unsigned char *slot_asg, int32_t *slot_nav_avail)
/*******************************************************************

   goci_slot_init

   purpose: set up tables to find the time for a pixel in a GOCI scene

   Returns type: int - 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      hid_t             file_id          I      ID of GOCI L1B file
      int *             dims             I      scene size in [ lines, pixels ]
      float *           slot_rel_time    O      table of relative time for the
      unsigned char *   slot_asg         O      size [# scene pixels, # scene
                                                lines] esitmate of the slot 
                                                for each scene pixel
      int32_t *         slot_nav_avail   O      flag indicating if the slot 
                                                navigation is valid from this 
                                                L1B (a 1 if valid, 0 if not)

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 26 Nov 2014     original development

 *******************************************************************/ {
    slot_nav_str *slot_nav;
    /*
     *  allocate storage for the table of slot nav info
     *  to be read from the data table: "Navigation for GOCI"
     */
    if ((slot_nav = (slot_nav_str *)
            malloc(NRECORDS * sizeof ( slot_nav_str))) == NULL) {
        printf("%s,%d-E Unable to allocate the slot navigation structure\n",
                __FILE__, __LINE__);
        return 1;
    }

    /* Calculate the size and the offsets of our struct members in memory */
    size_t nav_size = sizeof ( slot_nav_str);
    size_t nav_offset[NFIELDS] = {HOFFSET(slot_nav_str, band_num),
        HOFFSET(slot_nav_str, slot_num),
        HOFFSET(slot_nav_str, rel_time),
        HOFFSET(slot_nav_str, sc_att),
        HOFFSET(slot_nav_str, xo),
        HOFFSET(slot_nav_str, yo),
        HOFFSET(slot_nav_str, xs),
        HOFFSET(slot_nav_str, ys),
        HOFFSET(slot_nav_str, xpo),
        HOFFSET(slot_nav_str, ypo),
        HOFFSET(slot_nav_str, xps),
        HOFFSET(slot_nav_str, yps),
        HOFFSET(slot_nav_str, num_a_parm),
        HOFFSET(slot_nav_str, a_parm),
        HOFFSET(slot_nav_str, num_b_parm),
        HOFFSET(slot_nav_str, b_parm),
        HOFFSET(slot_nav_str, num_c_parm),
        HOFFSET(slot_nav_str, c_parm),
        HOFFSET(slot_nav_str, num_d_parm),
        HOFFSET(slot_nav_str, d_parm),
        HOFFSET(slot_nav_str, num_ap_parm),
        HOFFSET(slot_nav_str, ap_parm),
        HOFFSET(slot_nav_str, num_bp_parm),
        HOFFSET(slot_nav_str, bp_parm),
        HOFFSET(slot_nav_str, num_cp_parm),
        HOFFSET(slot_nav_str, cp_parm),
        HOFFSET(slot_nav_str, num_dp_parm),
        HOFFSET(slot_nav_str, dp_parm)};

    size_t nav_sizes[NFIELDS] = {sizeof ( slot_nav[0].band_num),
        sizeof ( slot_nav[0].slot_num),
        sizeof ( slot_nav[0].rel_time),
        sizeof ( slot_nav[0].sc_att),
        sizeof ( slot_nav[0].xo),
        sizeof ( slot_nav[0].yo),
        sizeof ( slot_nav[0].xs),
        sizeof ( slot_nav[0].ys),
        sizeof ( slot_nav[0].xpo),
        sizeof ( slot_nav[0].ypo),
        sizeof ( slot_nav[0].xps),
        sizeof ( slot_nav[0].yps),
        sizeof ( slot_nav[0].num_a_parm),
        sizeof ( slot_nav[0].a_parm),
        sizeof ( slot_nav[0].num_b_parm),
        sizeof ( slot_nav[0].b_parm),
        sizeof ( slot_nav[0].num_c_parm),
        sizeof ( slot_nav[0].c_parm),
        sizeof ( slot_nav[0].num_d_parm),
        sizeof ( slot_nav[0].d_parm),
        sizeof ( slot_nav[0].num_ap_parm),
        sizeof ( slot_nav[0].ap_parm),
        sizeof ( slot_nav[0].num_bp_parm),
        sizeof ( slot_nav[0].bp_parm),
        sizeof ( slot_nav[0].num_cp_parm),
        sizeof ( slot_nav[0].cp_parm),
        sizeof ( slot_nav[0].num_dp_parm),
        sizeof ( slot_nav[0].dp_parm)};
    /*
     *  goci slot navigation info end
     */

    hid_t grp_id;
    int32_t npix, nlin, step, nsx, nsy, trg_bnd, ix, iy;
    int32_t itile, ipix, ilin, lin_st, lin_en, pix_st, pix_en;
    hsize_t nfields, nrecords, nbnd, nslot;
    unsigned char *slot_asg_sml, *bnd_tile_lut, curtil;
    unsigned char box_pts[4];
    float minrad, crad;
    float min_t, max_t, rel_t;
    int32_t ibnd, ilut;

    npix = dims[1];
    nlin = dims[0];
    nbnd = 8;
    nslot = 16;
    bnd_tile_lut = NULL;

    /*
     *  set to the group
     */
    if ((grp_id = H5Gopen1(file_id, NAV_GRP)) < 0) {
        printf("%s,%d:E Unable to open group: %s\n", __FILE__, __LINE__, NAV_GRP);
        return 1;
    }
    /*
     *  see if it has the fields, records expected
     */
    *slot_nav_avail = 0;
    if (H5TBget_table_info(grp_id, TABLE_NAME, &nfields, &nrecords) < 0) {
        printf("%s,%d:E Unable to get table info for: %s\n", __FILE__,
                __LINE__, TABLE_NAME);
        return 1;
    }
    printf("# fields: %d, # records: %d\n", (int) nfields, (int) nrecords);
    if ((nfields != NFIELDS) || (nrecords != NRECORDS)) {
        *slot_nav_avail = 0;
        printf("%s,%d:W L1B GOCI input file\n does not have %d fields or %d records\n",
                __FILE__, __LINE__, (int) NFIELDS, (int) NRECORDS);
    } else {
        *slot_nav_avail = 1;
        /* read the table */
        if (H5TBread_table(grp_id, TABLE_NAME, nav_size, nav_offset,
                nav_sizes, slot_nav) < 0) {
            printf("%s,%d:E Unable to read table info for: %s\n", __FILE__,
                    __LINE__, TABLE_NAME);
            return 1;
        }

        if (*slot_nav_avail == 1) {
            /*
             *  next, make a super-grid generally defining the slot assignment
             */
            printf("Begin GOCI slot assignment\n");
            step = 18; /* reduces calls to goci_slot_nav in super and full grid 
                     steps for a 5k x 5k scene */
            trg_bnd = 7;
            nsx = 2 + npix / step;
            nsy = 2 + nlin / step;
            if ((slot_asg_sml = (unsigned char *)
                    malloc(nsx * nsy * sizeof ( unsigned char))) == NULL) {
                printf("%s,%d:E Unable to allocate space for slot_asg_sml array\n",
                        __FILE__, __LINE__);
                return 1;
                H5Gclose(grp_id);
            }
            if ((bnd_tile_lut = (unsigned char *)
                    malloc(nbnd * nslot * sizeof ( unsigned char))) == NULL) {
                printf("%s,%d:E Unable to allocate space for bnd_tile_lut array\n",
                        __FILE__, __LINE__);
                H5Gclose(grp_id);
                return 1;
                ;
            }
            *bnd_tile_lut = 254;
            for (iy = 0; iy < nsy; iy++) {
                ilin = iy * step;
                for (ix = 0; ix < nsx; ix++) {
                    ipix = ix * step;
                    minrad = 200.;
                    for (itile = 0; itile < 16; itile++) {
                        goci_slot_nav(ipix, ilin, trg_bnd, itile, slot_nav, nbnd, nslot,
                                bnd_tile_lut, &crad);
                        if (crad < minrad) {
                            *(slot_asg_sml + ix + nsx * iy) = itile;
                            minrad = crad;
                        }
                    }
                }
            }
            printf("GOCI supergrid made, starting full slot assignment\n");
            /*
             *  fill the full size slot assignment array
             */
            for (iy = 0; iy < (nsy - 1); iy++) {
                lin_st = iy * step;
                lin_en = (iy + 1) * step;
                if (lin_en > nlin) lin_en = nlin;
                for (ix = 0; ix < (nsx - 1); ix++) {
                    pix_st = ix * step;
                    pix_en = (ix + 1) * step;
                    if (pix_en > npix) pix_en = npix;
                    /*
                     *  if the supergrid box has all the same slot, assign that value 
                     *  to the pixel, line range
                     */
                    box_pts[0] = *(slot_asg_sml + ix + nsx * iy);
                    box_pts[1] = *(slot_asg_sml + (ix + 1) + nsx * iy);
                    box_pts[2] = *(slot_asg_sml + ix + nsx * (iy + 1));
                    box_pts[3] = *(slot_asg_sml + (ix + 1) + nsx * (iy + 1));
                    qsort(box_pts, 4, sizeof ( char), numcomp);

                    if (box_pts[0] == box_pts[3]) {
                        /*
                         *  assign one tile type
                         */
                        for (ilin = lin_st; ilin < lin_en; ilin++)
                            for (ipix = pix_st; ipix < pix_en; ipix++)
                                *(slot_asg + ipix + npix * ilin) = box_pts[0];
                    } else {
                        /*
                         *  step through each point and determine best slot
                         *  using up to 4 candidate slots
                         */
                        for (ilin = lin_st; ilin < lin_en; ilin++) {
                            for (ipix = pix_st; ipix < pix_en; ipix++) {
                                minrad = 200.;
                                curtil = -1;
                                for (itile = 0; itile < 4; itile++) {
                                    if (box_pts[itile] != curtil) {
                                        curtil = box_pts[itile];
                                        goci_slot_nav(ipix, ilin, trg_bnd, curtil, slot_nav,
                                                nbnd, nslot, bnd_tile_lut, &crad);
                                        if (crad < minrad) {
                                            *(slot_asg + ipix + npix * ilin) = curtil;
                                            minrad = crad;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    /*  for checkout, put in original supergrid point * 10 */
                    /*
                     *( slot_asg + pix_st + npix * lin_st ) = 
                     *( slot_asg_sml + ix + nsx * iy ) * 10;
                     */
                }
            }
        }
        /*
         *  set up the time offsets per slot.  As each band has a time offset, but
         *  we have no way to address individual bands, we'll use the mean of the
         *  band time range
         */
        /*
         *  loop thru slots and get each mean time from bands
         */
        for (itile = 0; itile < nslot; itile++) {
            min_t = 5000.;
            max_t = -5000.;
            for (ibnd = 0; ibnd < nbnd; ibnd++) {
                ilut = *(bnd_tile_lut + ibnd + nbnd * itile);
                rel_t = slot_nav[ilut].rel_time;
                if (rel_t > max_t) max_t = rel_t;
                if (rel_t < min_t) min_t = rel_t;
            }
            slot_rel_time[itile] = (min_t + max_t) / 2.;
        }
        printf("GOCI slot, time assignments completed\n");
        /* free space used and close the group */
        free(slot_asg_sml);
        free(bnd_tile_lut);
        free(slot_nav);
    }
    H5Gclose(grp_id);

    return 0;
}

int goci_slot_nav(int32_t ipix, int32_t ilin, int32_t bnd, int32_t itile,
        slot_nav_str *slot_nav, int32_t nbnd, int32_t nslot,
        unsigned char *bnd_tile_lut, float *lindist)
/*******************************************************************

   goci_slot_nav

   purpose:  transform a goci scene point into the point on a specific tile
     and return the largest linear normalized distance from the center

   Returns type: int - 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32_t           ipix             I      Pixel to transform
      int32_t           ilin             I      Line to transform
      int32_t           bnd              I      Band number
      int32_t           itile            I      Tile or slot number of GOCI
      slot_nav_str     *slot_nav         I      structure with transform 
                                        coefficients and normalization values
      int32_t           nbnd             I      number of bands
      int32_t           nslot            I      number of slots
      unsigned char *   bnd_tile_lut    I/O     storage for a look-up for 
                                             proper element in slot_nav
      float *           lindist          O      The largest normalized 
                                                linear distance from tile 
                                                center

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 1 Dec 2014      original development
      W. Robinson, SAIC 22 May 2017     switch criteria for slot choice from
                                        minimum distance from center to minimum
                                        linear distance from center

 *******************************************************************/ {

    int32_t ilut, nlut;
    int i, num_a_parm, num_b_parm, num_c_parm, num_d_parm;
    float xo, yo, xs, ys, *a_parm, *b_parm, *c_parm, *d_parm;
    float dist;
    double xn, yn, vec[16], xpn, ypn, numer, denom;

    /*
     *  set up the look-up to find proper slot, band
     */
    nlut = nbnd * nslot;
    if (*bnd_tile_lut == 254) {
        for (ilut = 0; ilut < nlut; ilut++)
            *(bnd_tile_lut + slot_nav[ilut].band_num +
                nbnd * slot_nav[ilut].slot_num) = ilut;
    }
    /*
     *  get normalization coefficients for scene to normalized scene
     */
    ilut = *(bnd_tile_lut + bnd + nbnd * itile);
    xo = slot_nav[ilut].xo;
    yo = slot_nav[ilut].yo;
    xs = slot_nav[ilut].xs;
    ys = slot_nav[ilut].ys;
    /*
     *  for the transform: scene -> tile normalized
     */
    num_a_parm = slot_nav[ilut].num_a_parm;
    a_parm = slot_nav[ilut].a_parm;
    num_b_parm = slot_nav[ilut].num_b_parm;
    b_parm = slot_nav[ilut].b_parm;
    num_c_parm = slot_nav[ilut].num_c_parm;
    c_parm = slot_nav[ilut].c_parm;
    num_d_parm = slot_nav[ilut].num_d_parm;
    d_parm = slot_nav[ilut].d_parm;
    /*
     *  make normalized scene location
     */
    xn = ((double) ipix - (double) xo) / (double) xs;
    yn = ((double) ilin - (double) yo) / (double) ys;
    /*
     *  transform from scene to tile normalized
     */
    for (i = 0; i < nslot; i++)
        vec[i] = 0.;
    vec[0] = 1.;
    vec[1] = xn;
    vec[2] = yn;
    vec[3] = xn * yn;
    vec[4] = pow(xn, 2.);
    vec[5] = pow(yn, 2.);
    vec[6] = pow(xn, 2.) * yn;
    vec[7] = pow(yn, 2.) * xn;
    vec[8] = pow(yn, 2.) * pow(xn, 2.);
    /*  X  */
    for (i = 0, numer = 0.; i < num_a_parm; i++)
        numer += vec[i] * a_parm[i];
    for (i = 0, denom = 1.; i < num_b_parm; i++)
        denom += vec[i + 1] * b_parm[i];

    xpn = numer / denom;

    /*  Y  */
    for (i = 0, numer = 0.; i < num_c_parm; i++)
        numer += vec[i] * c_parm[i];
    for (i = 0, denom = 1.; i < num_d_parm; i++)
        denom += vec[i + 1] * d_parm[i];

    ypn = numer / denom;
    /*
     *  find the normalized radius  *** OLD WAY
     *nradsq = (float) ( pow( xpn, 2. ) + pow( ypn, 2. ) );
     */
    /*
     *  find the smallest distance of this point from the 4 edges and 
     *  return 1 - this as the longest linear distance from the benter
     */
    *lindist = 200.;
    /* +X edge */
    dist = 1. - xpn;
    *lindist = (dist < *lindist) ? dist : *lindist;
    /* -X edge */
    dist = 1. + xpn;
    *lindist = (dist < *lindist) ? dist : *lindist;
    /* +Y edge */
    dist = 1. - ypn;
    *lindist = (dist < *lindist) ? dist : *lindist;
    /* -Y edge */
    dist = 1. + ypn;
    *lindist = (dist < *lindist) ? dist : *lindist;

    /* Make it a distance from center */
    *lindist = 1. - *lindist;
    return 0;
}

unsigned char goci_slot_time(int32_t ipix, int32_t ilin, goci_l1b_t *goci_l1b,
        float *rel_sec)
/*******************************************************************

   goci_slot_time

   purpose:  return the scene start relative time given a GOCI pixel, line

   Returns type: int - 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32_t           ipix             I      Pixel to transform
      int32_t           ilin             I      Line to transform
      goci_l1b_t *      goci_l1b         I      GOCI information structure
      float *           rel_sec          O      mean time of the pixel 
                                                relative to scene start time

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 10 Dec 2014     original development

 *******************************************************************/ {
    unsigned char p_slot;

    p_slot = *(goci_l1b->slot_asg + ipix + goci_l1b->npixels * ilin);

    *rel_sec = *(goci_l1b->slot_rel_time + p_slot);
    return p_slot;
}
