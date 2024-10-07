#include "scene_meta.h"

static scnstr meta;

static int first = 1;
static int second = 1;
static int centerfound = 0;
static int start_node = UNKNOWNNODE;
static int end_node = UNKNOWNNODE;
static int ascdsc = UNKNOWNNODE;
static int daynight = UNKNOWNSCENE;
static float lastLat = -999;
static float lastLon0 = -999;
static float lastLonN = -999;
static int lastScan = -1;
static float northern_lat = -90.0;
static float southern_lat = +90.0;
static float eastern_lon = -180.0;
static float western_lon = +180.0;

static float lat_rec[10] ={-999., -999., -999., -999., -999., -999., -999., -999., -999., -999,};
static int lat_rec_ptr = 0;

static char *nodestr [3] = {"Ascending", "Descending", "Unknown"};
static char *daynightstr[4] = {"Day", "Night", "Mixed", "Unknown"};

float westernmost(float lon1, float lon2) {
    if (lon1 == -999.) lon1 = lon2;

    if (fabs(lon1 - lon2) < 190.0)
        return ( MIN(lon1, lon2));
    else
        return ( MAX(lon1, lon2));
}

float easternmost(float lon1, float lon2) {
    if (lon1 == -999.) lon1 = lon2;

    if (fabs(lon1 - lon2) < 190.0)
        return ( MAX(lon1, lon2));
    else
        return ( MIN(lon1, lon2));
}

void scene_meta_init(void) {
    strcpy(meta.daynight, daynightstr[UNKNOWNSCENE]);
    strcpy(meta.start_node, nodestr[UNKNOWNNODE]);
    strcpy(meta.end_node, nodestr[UNKNOWNNODE]);
    meta.start_time = 0;
    meta.center_time = 0;
    meta.end_time = 0;

    meta.start_year = 0;
    meta.start_day = 0;
    meta.start_msec = 0;

    meta.end_year = 0;
    meta.end_day = 0;
    meta.end_msec = 0;

    meta.northern_lat = 90.0;
    meta.southern_lat = -90.0;
    meta.eastern_lon = 180.0;
    meta.western_lon = -180.0;

    meta.start_center_lon = 0.0;
    meta.start_center_lat = 0.0;
    meta.end_center_lon = 0.0;
    meta.end_center_lat = 0.0;
}

void scene_meta_put(l1str *l1rec) {
    static int firstCall = 1;
    static int32_t cscan;
    int16_t year, day;
    double msec;

    int32_t cpix = l1rec->npix / 2;
    int32_t npix = l1rec->npix;
    int32_t epix = -1;
    int32_t spix = -1;
    int32_t ip;

    if (firstCall) {
        firstCall = 0;
        scene_meta_init();
        if (l1_input->eline < 1)
            cscan = ((l1rec->l1file->nscan - 1) - MAX(l1_input->sline - 1, 0)) / 2;
        else
            cscan = ((l1_input->eline - 1) - MAX(l1_input->sline - 1, 0)) / 2;
    }

    /* ignore buffering scans */
    if (l1rec->iscan == lastScan) return;
    if (l1rec->iscan < MAX(l1_input->sline - 1, 0)) return;
    if (l1_input->eline < 1) {
        if (l1rec->iscan > (l1rec->l1file->nscan - 1))
            return;
    } else {
        if (l1rec->iscan > (l1_input->eline - 1))
            return;
    }

    for (ip = 0; ip < cpix; ip++) {
        if ((l1rec->flags[ip] & NAVFAIL) == 0) {
            spix = ip;
            break;
        }
    }
    for (ip = l1rec->npix - 1; ip > cpix; ip--) {
        if ((l1rec->flags[ip] & NAVFAIL) == 0) {
            epix = ip;
            break;
        }
    }

    if (spix != -1 &&
            epix != -1 &&
            (l1rec->flags[cpix] & NAVFAIL) == 0) {

        if (first) {

            first = 0;
            daynight = UNKNOWNSCENE;

            meta.upperleft_lon = l1rec->lon[0];
            meta.upperright_lon = l1rec->lon[npix - 1];
            meta.upperleft_lat = l1rec->lat[0];
            meta.upperright_lat = l1rec->lat[npix - 1];
            meta.start_center_lon = l1rec->lon[cpix];
            meta.start_center_lat = l1rec->lat[cpix];
            unix2yds(l1rec->scantime, &year, &day, &msec);
            meta.start_year = (int32_t) year;
            meta.start_day = (int32_t) day;
            meta.start_msec = (int32_t) (msec * 1.e3);
            meta.start_time = l1rec->scantime;

        } else if (second) {

            second = 0;
            if (l1rec->lat[cpix] > lastLat)
                ascdsc = 0; /* Ascending  */
            else
                ascdsc = 1; /* Descending */
            start_node = ascdsc;
            end_node = ascdsc;

            ascdsc = (l1rec->l1file->sensorID == CZCS || l1rec->l1file->sensorID == OCI) ? 1 - ascdsc : ascdsc;
            if (ascdsc == 1) { /* Descending */
                western_lon = westernmost(lastLon0, l1rec->lon[spix]);
                eastern_lon = easternmost(lastLonN, l1rec->lon[epix]);
            } else { /* Ascending */
                western_lon = westernmost(lastLonN, l1rec->lon[epix]);
                eastern_lon = easternmost(lastLon0, l1rec->lon[spix]);
            }

        } else {

            if (l1rec->lat[cpix] > lastLat)
                ascdsc = 0; /* Ascending  */
            else
                ascdsc = 1; /* Descending */
            end_node = ascdsc;

            ascdsc = (l1rec->l1file->sensorID == CZCS || l1rec->l1file->sensorID == OCI) ? 1 - ascdsc : ascdsc;
            if (ascdsc == 1) { /* Descending */
                western_lon = westernmost(western_lon, l1rec->lon[spix]);
                eastern_lon = easternmost(eastern_lon, l1rec->lon[epix]);
            } else { /*Ascending */
                western_lon = westernmost(western_lon, l1rec->lon[epix]);
                eastern_lon = easternmost(eastern_lon, l1rec->lon[spix]);
            }

            meta.lowerleft_lon = l1rec->lon[0];
            meta.lowerright_lon = l1rec->lon[npix - 1];
            meta.lowerleft_lat = l1rec->lat[0];
            meta.lowerright_lat = l1rec->lat[npix - 1];
            meta.end_center_lon = l1rec->lon[cpix];
            meta.end_center_lat = l1rec->lat[cpix];
            unix2yds(l1rec->scantime, &year, &day, &msec);
            meta.end_year = (int32_t) year;
            meta.end_day = (int32_t) day;
            meta.end_msec = (int32_t) (msec * 1.e3);
            meta.end_time = l1rec->scantime;
        }

        if (l1rec->iscan >= cscan && !centerfound) {
            meta.scene_center_lon = l1rec->lon[cpix];
            meta.scene_center_lat = l1rec->lat[cpix];
            meta.scene_center_solz = l1rec->solz[cpix];
            meta.center_time = l1rec->scantime;
            centerfound = 1;
        }

        for (ip = spix; ip <= epix; ip++) {
            if ((l1rec->flags[ip] & NAVFAIL) == 0) {
                northern_lat = MAX(northern_lat, l1rec->lat[ip]);
                southern_lat = MIN(southern_lat, l1rec->lat[ip]);
                if (daynight != DAYANDNIGHT) {
                    if (l1rec->solz[ip] >= SOLZNIGHT) {
                        if (daynight == DAYSCENE)
                            daynight = DAYANDNIGHT;
                        else
                            daynight = NIGHTSCENE;
                    } else {
                        if (daynight == NIGHTSCENE)
                            daynight = DAYANDNIGHT;
                        else
                            daynight = DAYSCENE;
                    }
                }
            }
        }

        if (lat_rec[lat_rec_ptr] > -100.)
            lastLat = lat_rec[lat_rec_ptr];
        else
            lastLat = l1rec->lat[cpix];
        lat_rec[lat_rec_ptr++] = l1rec->lat[cpix];
        lat_rec_ptr = lat_rec_ptr % 10;
        lastLon0 = l1rec->lon[spix];
        lastLonN = l1rec->lon[epix];
        lastScan = l1rec->iscan;

        strcpy(meta.start_node, nodestr[start_node]);
        strcpy(meta.end_node, nodestr[end_node]);
        strcpy(meta.daynight, daynightstr[daynight]);

        meta.northern_lat = northern_lat;
        meta.southern_lat = southern_lat;
        meta.eastern_lon = eastern_lon;
        meta.western_lon = western_lon;
    }

}

scnstr *scene_meta_get(void) {
    return (&meta);
}

void scene_meta_write(idDS ds_id) {
    scnstr *m = &meta;

    if (ds_id.fftype == DS_NCDF) {

        // 1994-11-05T13:15:30Z
        SetChrGA(ds_id, "time_coverage_start", unix2isodate(m->start_time, 'G'));

        SetChrGA(ds_id, "time_coverage_end", unix2isodate(m->end_time, 'G'));

        SetF32GA(ds_id, "start_center_longitude", m->start_center_lon);
        SetF32GA(ds_id, "start_center_latitude", m->start_center_lat);
        SetF32GA(ds_id, "end_center_longitude", m->end_center_lon);
        SetF32GA(ds_id, "end_center_latitude", m->end_center_lat);
        SetF32GA(ds_id, "northernmost_latitude", m->northern_lat);
        SetF32GA(ds_id, "southernmost_latitude", m->southern_lat);
        SetF32GA(ds_id, "easternmost_longitude", m->eastern_lon);
        SetF32GA(ds_id, "westernmost_longitude", m->western_lon);
        SetChrGA(ds_id, "geospatial_lat_units", "degrees_north");
        SetChrGA(ds_id, "geospatial_lon_units", "degrees_east");
        SetF32GA(ds_id, "geospatial_lat_max", m->northern_lat);
        SetF32GA(ds_id, "geospatial_lat_min", m->southern_lat);
        SetF32GA(ds_id, "geospatial_lon_max", m->eastern_lon);
        SetF32GA(ds_id, "geospatial_lon_min", m->western_lon);
        SetChrGA(ds_id, "startDirection", m->start_node);
        SetChrGA(ds_id, "endDirection", m->end_node);
        SetChrGA(ds_id, "day_night_flag", m->daynight);

    } else {

        // write to an HDF4 file
        SetChrGA(ds_id, "Start Time", unix2ydhmsf(m->start_time, 'G'));
        SetI16GA(ds_id, "Start Year", m->start_year);
        SetI16GA(ds_id, "Start Day", m->start_day);
        SetI32GA(ds_id, "Start Millisec", m->start_msec);

        SetF32GA(ds_id, "Upper Left Longitude", m->upperleft_lon);
        SetF32GA(ds_id, "Upper Right Longitude", m->upperright_lon);
        SetF32GA(ds_id, "Upper Left Latitude", m->upperleft_lat);
        SetF32GA(ds_id, "Upper Right Latitude", m->upperright_lat);

        SetF32GA(ds_id, "Start Center Longitude", m->start_center_lon);
        SetF32GA(ds_id, "Start Center Latitude", m->start_center_lat);

        SetChrGA(ds_id, "Latitude Units", "degrees North");
        SetChrGA(ds_id, "Longitude Units", "degrees East");

        SetChrGA(ds_id, "Scene Center Time", unix2ydhmsf(m->center_time, 'G'));
        SetF32GA(ds_id, "Scene Center Longitude", m->scene_center_lon);
        SetF32GA(ds_id, "Scene Center Latitude", m->scene_center_lat);
        SetF32GA(ds_id, "Scene Center Solar Zenith", m->scene_center_solz);

        SetChrGA(ds_id, "End Time", unix2ydhmsf(m->end_time, 'G'));
        SetI16GA(ds_id, "End Year", m->end_year);
        SetI16GA(ds_id, "End Day", m->end_day);
        SetI32GA(ds_id, "End Millisec", m->end_msec);

        SetF32GA(ds_id, "Lower Left Longitude", m->lowerleft_lon);
        SetF32GA(ds_id, "Lower Right Longitude", m->lowerright_lon);
        SetF32GA(ds_id, "Lower Left Latitude", m->lowerleft_lat);
        SetF32GA(ds_id, "Lower Right Latitude", m->lowerright_lat);

        SetF32GA(ds_id, "End Center Longitude", m->end_center_lon);
        SetF32GA(ds_id, "End Center Latitude", m->end_center_lat);

        SetF32GA(ds_id, "Northernmost Latitude", m->northern_lat);
        SetF32GA(ds_id, "Southernmost Latitude", m->southern_lat);
        SetF32GA(ds_id, "Easternmost Longitude", m->eastern_lon);
        SetF32GA(ds_id, "Westernmost Longitude", m->western_lon);

        SetChrGA(ds_id, "Start Node", m->start_node);
        SetChrGA(ds_id, "End Node", m->end_node);
        SetChrGA(ds_id, "Day or Night", m->daynight);

    }

}
