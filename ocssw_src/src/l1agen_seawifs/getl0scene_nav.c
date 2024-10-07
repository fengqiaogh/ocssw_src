/* -------------------------------------------------------------------- */
/* getl0scene_nav() - adds coordinate info to the scene records.        */
/*                                                                      */
/* B. A. Franz, GSC, November 1997                                      */
/* ---------------------------------------------------------------------*/

#include <stdio.h>
#include <time.h>
#include <math.h>
#include "swl0_proto.h"
#include <timeutils.h>

FLOAT32 westernmost(FLOAT32 lon1, FLOAT32 lon2);
FLOAT32 easternmost(FLOAT32 lon1, FLOAT32 lon2);
int getGeoNav(navblk_sType *navblk,
        INT32 spix,
        INT32 ipix,
        INT32 npix,
        FLOAT32 lat[],
        FLOAT32 lon[],
        FLOAT32 solz[],
        FLOAT32 sola[],
        FLOAT32 senz[],
        FLOAT32 sena[]);

/* ---------------------------------------------------------------------*/

/* ---------------------------------------------------------------------*/
int getl0scene_nav(FLOAT32 xnodel,
        INT32 tnode,
        navblk_sType navblk[],
        tilt_states_sType *tiltblk,
        swl0scene *scene) {
    int i, j;
    int nbad = 0;
    INT32 first = -1;
    INT32 first_plus_1 = -1;
    INT32 last = -1;
    INT32 center = -1;
    INT16 srec;
    INT16 erec;
    INT16 ascdsc;
    char node[2][11] = {"Ascending", "Descending"};

    INT32 spix;
    INT32 ipix;
    INT32 npix;
    INT32 cpix;
    FLOAT32 lat [NPIXLAC];
    FLOAT32 lon [NPIXLAC];
    FLOAT32 solz[NPIXLAC];
    FLOAT32 sola[NPIXLAC];
    FLOAT32 senz[NPIXLAC];
    FLOAT32 sena[NPIXLAC];

    FLOAT32 lastLat;
    FLOAT32 lastLon0;
    FLOAT32 lastLonN;

    FLOAT64 sec;
    INT16 day;
    INT16 year;

    BYTE navFlag[MAXFRAMES];


    /*                                                                */
    /* Set any datatype specific parameters                           */
    /*                                                                */
    if (scene->mnftype == GACTYPE) {
        spix = SPIXGAC;
        ipix = IPIXGAC;
        npix = NPIXGAC;
    } else {
        spix = SPIXLAC;
        ipix = IPIXLAC;
        npix = NPIXLAC;
    }
    cpix = npix / 2;


    /*                                                                */
    /* Break-down time for later use                                  */
    /*                                                                */
    unix2yds(scene->stime, &year, &day, &sec);

    /*                                                                */
    /* Copy node and tilt info                                        */
    /*                                                                */
    scene->node_lon = xnodel;
    scene->node_time = yds2unix(year, day, (FLOAT64) tnode / 1000.);
    scene->ntilts = tiltblk->ntilts;
    memcpy(scene->tilt_flags, tiltblk->tilt_flags,
            sizeof (scene->tilt_flags));
    memcpy(scene->tilt_ranges, tiltblk->tilt_ranges,
            sizeof (scene->tilt_ranges));

    /*                                                                */
    /* Generate navigation quality index for this scene               */
    /*                                                                */
    for (i = 0; i < scene->nscan; i++) {
        if (getGeoNav(&navblk[i], spix, ipix, npix, lat, lon, solz, sola, senz, sena))
            navFlag[i] = 1; /* Good nav */
        else {
            navFlag[i] = 0; /* Bad  nav */
            nbad++;
            if ((scene->mnftype != LUNTYPE) && (scene->mnftype != SOLTYPE) && nbad < 10)
                printf("-W- %s: bad navigation at scan %d\n", __FILE__, i);
        }
        /*                                                            */
        /* Update nflag[0] with results of more stringent test        */
        /*                                                            */
        if ((navblk[i].nflag[0] == 0) && (navFlag[i] == 0)) {
            navblk[i].nflag[0] = 1;
            navblk[i].nflag[1] = 1;
        }
    }

    /*                                                                */
    /* Require at least 3 navigatable scans, else set nav vals to 0   */
    /* Do the same for lunar and solar calibration, but return OK.    */
    /*                                                                */
    if ((nbad > scene->nscan - 3) || (scene->mnftype == LUNTYPE)
            || (scene->mnftype == SOLTYPE)) {
        scene->upper_left_lon = -999.0;
        scene->upper_left_lat = -99.0;
        scene->upper_right_lon = -999.0;
        scene->upper_right_lat = -99.0;
        scene->lower_left_lon = -999.0;
        scene->lower_left_lat = -99.0;
        scene->lower_right_lon = -999.0;
        scene->lower_right_lat = -99.0;
        scene->center_lon = -999.0;
        scene->center_lat = -99.0;
        scene->center_solz = -999.0;
        scene->start_center_lon = -999.0;
        scene->start_center_lat = -99.0;
        scene->end_center_lon = -999.0;
        scene->end_center_lat = -99.0;
        scene->start_center_lon = -999.0;
        scene->northern_lat = -99.0;
        scene->southern_lat = -99.0;
        scene->eastern_lon = -999.0;
        scene->western_lon = -999.0;

        strncpy(scene->start_node, node[0], 11);
        strncpy(scene->end_node, node[0], 11);

        if (scene->type != HRPT)
            return (0);
        else
            return (1);
    }

    /*                                                                */
    /* Find first and last scan with good navigation                  */
    /*                                                                */
    for (i = 0; i < scene->nscan; i++) {
        if (navFlag[i]) {
            if (first < 0)
                first = i;
            else if (first_plus_1 < 0)
                first_plus_1 = i;
            last = i;
        }
    }

    /*                                                                */
    /* Find center scan with good navigation                          */
    /*                                                                */
    center = (last + first) / 2;
    while ((center > first) && (!navFlag[center]))
        center--;
    scene->center_scan_line = center + 1;


    /*                                                                */
    /* Geolocate first, central, and last navigatable scans, and      */
    /* determine scene limits and ascending, descending state         */
    /*                                                                */

    scene->northern_lat = -90.0;
    scene->southern_lat = 90.0;
    scene->eastern_lon = 0.0;
    scene->western_lon = 0.0;
    lastLat = 90.0;

    for (i = 0; i < scene->nscan; i++) {

        if (navFlag[i]) {

            geonav_(navblk[i].orb_vec,
                    navblk[i].sen_mat,
                    navblk[i].scan_ell,
                    navblk[i].sun_ref,
                    &spix, &ipix, &npix,
                    lat, lon, solz, sola, senz, sena);

            /*                                                        */
            /* Ascending or Descending ?                              */
            /*                                                        */
            if (lat[cpix] > lastLat)
                ascdsc = 0;
            else
                ascdsc = 1;

            /*                                                        */
            /* Update first, last, center info                        */
            /*                                                        */
            if (i == first) {
                scene->upper_left_lon = lon[0];
                scene->upper_left_lat = lat[0];
                scene->upper_right_lon = lon[npix - 1];
                scene->upper_right_lat = lat[npix - 1];
                scene->start_center_lon = lon[cpix];
                scene->start_center_lat = lat[cpix];

            } else if (i == first_plus_1) {
                strncpy(scene->start_node, node[ascdsc], 11);

                /* now that we know which direction we are going, */
                /* reinitialize eastermost and westermost         */
                if (ascdsc == 1) {
                    scene->western_lon = westernmost(lastLon0, lon[0]);
                    scene->eastern_lon = easternmost(lastLonN, lon[npix - 1]);
                } else {
                    scene->western_lon = westernmost(lastLonN, lon[npix - 1]);
                    scene->eastern_lon = easternmost(lastLon0, lon[0]);
                }

            } else if (i == center) {
                scene->center_lon = lon[cpix];
                scene->center_lat = lat[cpix];
                scene->center_solz = solz[cpix];

            } else if (i == last) {
                scene->lower_left_lon = lon[0];
                scene->lower_left_lat = lat[0];
                scene->lower_right_lon = lon[npix - 1];
                scene->lower_right_lat = lat[npix - 1];
                scene->end_center_lon = lon[cpix];
                scene->end_center_lat = lat[cpix];
                strncpy(scene->end_node, node[ascdsc], 11);
            }

            /*                                                        */
            /* Update extrema                                         */
            /*                                                        */
            for (j = 0; j < npix; j++) {
                scene->northern_lat = MAX(scene->northern_lat, lat[j]);
                scene->southern_lat = MIN(scene->southern_lat, lat[j]);
            }

            if (ascdsc == 1) {
                /* Descending */
                scene->western_lon = westernmost(scene->western_lon,
                        lon[0]);
                scene->eastern_lon = easternmost(scene->eastern_lon,
                        lon[npix - 1]);
            } else {
                /*Ascending */
                scene->western_lon = westernmost(scene->western_lon,
                        lon[npix - 1]);
                scene->eastern_lon = easternmost(scene->eastern_lon,
                        lon[0]);
            }

            lastLat = lat[cpix];
            lastLon0 = lon[0];
            lastLonN = lon[npix - 1];
        }
    }


    /*                                                                */
    /* Geolocate scans at start and end of each tilt period. Try to   */
    /* compensate for nav errors, but if all lines are bad, set to 0  */
    /*                                                                */
    for (i = 0; i < tiltblk->ntilts; i++) {

        srec = tiltblk->tilt_ranges[i][0] - 1;
        erec = tiltblk->tilt_ranges[i][1] - 1;

        while (!navFlag[srec]) {
            if (srec == erec) break;
            srec++;
        }

        while (!navFlag[erec]) {
            if (erec == srec) break;
            erec--;
        }

        if (navFlag[srec]) {

            geonav_(navblk[srec].orb_vec,
                    navblk[srec].sen_mat,
                    navblk[srec].scan_ell,
                    navblk[srec].sun_ref,
                    &spix, &ipix, &npix,
                    lat, lon, solz, sola, senz, sena);

            scene->tilt_lons[i][0][0] = lon[0];
            scene->tilt_lons[i][0][1] = lon[npix - 1];
            scene->tilt_lats[i][0][0] = lat[0];
            scene->tilt_lats[i][0][1] = lat[npix - 1];

        } else {
            scene->tilt_lons[i][0][0] = 0.0;
            scene->tilt_lons[i][0][1] = 0.0;
            scene->tilt_lats[i][0][0] = 0.0;
            scene->tilt_lats[i][0][1] = 0.0;
        }
        if (navFlag[erec]) {

            geonav_(navblk[erec].orb_vec,
                    navblk[erec].sen_mat,
                    navblk[erec].scan_ell,
                    navblk[erec].sun_ref,
                    &spix, &ipix, &npix,
                    lat, lon, solz, sola, senz, sena);

            scene->tilt_lons[i][1][0] = lon[0];
            scene->tilt_lons[i][1][1] = lon[npix - 1];
            scene->tilt_lats[i][1][0] = lat[0];
            scene->tilt_lats[i][1][1] = lat[npix - 1];

        } else {
            scene->tilt_lons[i][1][0] = 0.0;
            scene->tilt_lons[i][1][1] = 0.0;
            scene->tilt_lats[i][1][0] = 0.0;
            scene->tilt_lats[i][1][1] = 0.0;
        }

    }

    return (0);
}

/*                                                          */
/* getGeoNav() - runs geonav and makes additional checks on */
/* navigation quality.  Returns 1 if OK, 0 if bad.          */

/*                                                          */
int getGeoNav(navblk_sType *navblk,
        INT32 spix,
        INT32 ipix,
        INT32 npix,
        FLOAT32 lat[],
        FLOAT32 lon[],
        FLOAT32 solz[],
        FLOAT32 sola[],
        FLOAT32 senz[],
        FLOAT32 sena[]) {
    if (navblk->nflag[0] != 0)
        return (0);

    else {

        geonav_(navblk->orb_vec,
                navblk->sen_mat,
                navblk->scan_ell,
                navblk->sun_ref,
                &spix, &ipix, &npix,
                lat, lon, solz, sola, senz, sena);

        if (isnan(lat[0]) || lat[0] == 999.0)
            return (0);

        if (isnan(lat[npix / 2]) || lat[npix / 2] == 999.0)
            return (0);

        if (isnan(lat[npix - 1]) || lat[npix - 1] == 999.0)
            return (0);
    }

    return (1);
}

FLOAT32 westernmost(FLOAT32 lon1, FLOAT32 lon2) {
    if (fabs(lon1 - lon2) < 190.0)
        return ( MIN(lon1, lon2));
    else
        return ( MAX(lon1, lon2));
}

FLOAT32 easternmost(FLOAT32 lon1, FLOAT32 lon2) {
    if (fabs(lon1 - lon2) < 190.0)
        return ( MAX(lon1, lon2));
    else
        return ( MIN(lon1, lon2));
}



