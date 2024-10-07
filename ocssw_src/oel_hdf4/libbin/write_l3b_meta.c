/*
 * write_l3b_meta.c
 *
 *  Created on: Jun 12, 2014
 *      Author: dshea
 */

#include <netcdf.h>
#include <meta_l3b.h>
#include <timeutils.h>
#include <sensorInfo.h>
#include <genutils.h>

#include <hdf.h>
#include <mfhdf.h>


void calculate_temporal_range(meta_l3bType *meta_l3b) {
    double timediff;
    double days;
    int nday;
    char *period;

    period = (char *) malloc(SM_ATTRSZ * sizeof (char));

    timediff = meta_l3b->endTime - meta_l3b->startTime;
    days = timediff / 86400.;

    if (days < 0.75) {
        nday = floor(days * 24.) + 1;
        sprintf(period, "%d-hour", nday);
        strcpy(meta_l3b->prod_type, period);
    } else if (days >= 0.75 && days < 1.5) {
        strcpy(meta_l3b->prod_type, "day");
    } else if (days >= 1.5 && days < 27.5) {
        nday = round(days);
        sprintf(period, "%d-day", nday);
        strcpy(meta_l3b->prod_type, period);
    } else if (days >= 27.5 && days < 31.5) {
        strcpy(meta_l3b->prod_type, "month");
    } else if (days > 31.5 && days < 32.75) {
        // special case for rolling 32D composites
        strcpy(meta_l3b->prod_type, "32-day");
    } else if (days >= 32.75 && days < 364) {
        nday = round(days / 30.);
        sprintf(period, "%d-month", nday);
        strcpy(meta_l3b->prod_type, period);
    } else if (days > 364 && days < 367) {
        strcpy(meta_l3b->prod_type, "year");
    } else if (days >= 367) {
        nday = round(days / 365.);
        sprintf(period, "%d-year", nday);
        strcpy(meta_l3b->prod_type, period);
    }

    free(period);
}

int write_l3b_meta_hdf4(int32_t sd_id, meta_l3bType *meta_l3b) {
    /* Write Global Attributes */
    /* ----------------------- */
    char buf[LG_ATTRSZ];

    strcpy(buf, meta_l3b->product_name);
    SDsetattr(sd_id, "Product Name", DFNT_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->sensor_name);
    strcat(buf, " Level-3 Binned Data");
    SDsetattr(sd_id, "Title", DFNT_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->sensor_name);
    SDsetattr(sd_id, "Sensor Name", DFNT_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->mission);
    SDsetattr(sd_id, "Mission", DFNT_CHAR, strlen(buf) + 1, buf);

    calculate_temporal_range(meta_l3b);
    strcpy(buf, meta_l3b->prod_type);
    SDsetattr(sd_id, "Product Type", DFNT_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->pversion);
    SDsetattr(sd_id, "Processing Version", DFNT_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->soft_name);
    SDsetattr(sd_id, "Software Name", DFNT_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->soft_ver);
    SDsetattr(sd_id, "Software Version", DFNT_CHAR, strlen(buf) + 1, buf);

    SDsetattr(sd_id, "Start Orbit", DFNT_INT32, 1, &meta_l3b->start_orb);
    SDsetattr(sd_id, "End Orbit", DFNT_INT32, 1, &meta_l3b->end_orb);

    strcpy(buf, meta_l3b->ptime);
    SDsetattr(sd_id, "Processing Time", DFNT_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->proc_con);
    SDsetattr(sd_id, "Processing Control", DFNT_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->input_parms);
    SDsetattr(sd_id, "Input Parameters", DFNT_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->infiles);
    SDsetattr(sd_id, "Input Files", DFNT_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->flag_names);
    SDsetattr(sd_id, "L2 Flag Names", DFNT_CHAR, strlen(buf) + 1, buf);

    // if bin_syear isn't set - assume (yes, I know) the start/end year/day need to be set
    int16_t syear, sday, eyear, eday;
    double msec;
    int32_t smsec, emsec;
    unix2yds(meta_l3b->startTime, &syear, &sday, &msec);
    smsec = (int32_t) (msec * 1000.);
    unix2yds(meta_l3b->endTime, &eyear, &eday, &msec);
    emsec = (int32_t) (msec * 1000.);

    //    if (meta_l3b->bin_syear == 0){
    //        meta_l3b->bin_syear = syear;
    //        meta_l3b->bin_sday = sday;
    //        meta_l3b->bin_eyear = eyear;
    //        meta_l3b->bin_eday = eday;
    //    }
    SDsetattr(sd_id, "Period Start Year", DFNT_INT16, 1, &syear);
    SDsetattr(sd_id, "Period Start Day", DFNT_INT16, 1, &sday);
    SDsetattr(sd_id, "Period End Year", DFNT_INT16, 1, &eyear);
    SDsetattr(sd_id, "Period End Day", DFNT_INT16, 1, &eday);

    strcpy(buf, ydhmsf(meta_l3b->startTime, 'G'));
    SDsetattr(sd_id, "Start Time", DFNT_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, ydhmsf(meta_l3b->endTime, 'G'));
    SDsetattr(sd_id, "End Time", DFNT_CHAR, strlen(buf) + 1, buf);

    SDsetattr(sd_id, "Start Year", DFNT_INT16, 1, &syear);
    SDsetattr(sd_id, "Start Day", DFNT_INT16, 1, &sday);
    SDsetattr(sd_id, "Start Millisec", DFNT_INT32, 1, &smsec);
    SDsetattr(sd_id, "End Year", DFNT_INT16, 1, &eyear);
    SDsetattr(sd_id, "End Day", DFNT_INT16, 1, &eday);
    SDsetattr(sd_id, "End Millisec", DFNT_INT32, 1, &emsec);

    strcpy(buf, "degrees North");
    SDsetattr(sd_id, "Latitude Units", DFNT_CHAR, strlen(buf) + 1, buf);
    strcpy(buf, "degrees East");
    SDsetattr(sd_id, "Longitude Units", DFNT_CHAR, strlen(buf) + 1, buf);

    SDsetattr(sd_id, "Northernmost Latitude", DFNT_FLOAT32, 1, &meta_l3b->north);
    SDsetattr(sd_id, "Southernmost Latitude", DFNT_FLOAT32, 1, &meta_l3b->south);
    SDsetattr(sd_id, "Easternmost Longitude", DFNT_FLOAT32, 1, &meta_l3b->east);
    SDsetattr(sd_id, "Westernmost Longitude", DFNT_FLOAT32, 1, &meta_l3b->west);

    int32 tmpInt = meta_l3b->data_bins;
    SDsetattr(sd_id, "Data Bins", DFNT_INT32, 1, &tmpInt);
    SDsetattr(sd_id, "Percent Data Bins", DFNT_FLOAT32, 1,
            &meta_l3b->pct_databins);

    strcpy(buf, meta_l3b->units);
    SDsetattr(sd_id, "Units", DFNT_CHAR, strlen(buf) + 1, buf);
    strcpy(buf, resolution2str(meta_l3b->resolution));
    SDsetattr(sd_id, "Bin Resolution", DFNT_CHAR, strlen(buf) + 1, buf);

    return 0;
}

int write_l3b_meta_netcdf4(idDS ds_id, meta_l3bType *meta_l3b, int write64bit) {
    // Write Global Attributes
    // -----------------------
    char buf[LG_ATTRSZ];

    strcpy(buf, meta_l3b->product_name);
    setAttr(ds_id, "product_name", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->sensor);
    strcat(buf, " Level-3 Binned Data");
    setAttr(ds_id, "title", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->sensor);
    setAttr(ds_id, "instrument", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->mission);
    if (strcmp(meta_l3b->mission, "") != 0)
        setAttr(ds_id, "platform", NC_CHAR, strlen(buf) + 1, buf);

    calculate_temporal_range(meta_l3b);
    strcpy(buf, meta_l3b->prod_type);
    setAttr(ds_id, "temporal_range", NC_CHAR, strlen(buf) + 1, buf);

    setAttr(ds_id, "start_orbit_number", NC_INT, 1, &meta_l3b->start_orb);
    setAttr(ds_id, "end_orbit_number", NC_INT, 1, &meta_l3b->end_orb);

    strcpy(buf, meta_l3b->ptime);
    setAttr(ds_id, "date_created", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->pversion);
    setAttr(ds_id, "processing_version", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, "satellite observations from ");
    strncat(buf, meta_l3b->sensor,sizeof(buf) - strlen(meta_l3b->sensor) - 1);
    if (strcmp(meta_l3b->mission, "") != 0){
        strcat(buf, "-");
        strncat(buf, meta_l3b->mission,sizeof(buf) - strlen(meta_l3b->mission) - 1);
    }
    setAttr(ds_id, "source", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->proc_con);
    setAttr(ds_id, "history", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, unix2isodate(meta_l3b->startTime, 'G'));
    setAttr(ds_id, "time_coverage_start", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, unix2isodate(meta_l3b->endTime, 'G'));
    setAttr(ds_id, "time_coverage_end", NC_CHAR, strlen(buf) + 1, buf);

    setAttr(ds_id, "northernmost_latitude", NC_FLOAT, 1,
            (VOIDP) & meta_l3b->north);
    setAttr(ds_id, "southernmost_latitude", NC_FLOAT, 1,
            (VOIDP) & meta_l3b->south);
    setAttr(ds_id, "easternmost_longitude", NC_FLOAT, 1,
            (VOIDP) & meta_l3b->east);
    setAttr(ds_id, "westernmost_longitude", NC_FLOAT, 1,
            (VOIDP) & meta_l3b->west);
    double val = (double) meta_l3b->north;
    setAttr(ds_id, "geospatial_lat_max", NC_DOUBLE, 1, (VOIDP) & val);
    val = (double) meta_l3b->south;
    setAttr(ds_id, "geospatial_lat_min", NC_DOUBLE, 1, (VOIDP) & val);
    val = (double) meta_l3b->east;
    setAttr(ds_id, "geospatial_lon_max", NC_DOUBLE, 1, (VOIDP) & val);
    val = (double) meta_l3b->west;
    setAttr(ds_id, "geospatial_lon_min", NC_DOUBLE, 1, (VOIDP) & val);

    strcpy(buf, "degrees_north");
    setAttr(ds_id, "geospatial_lat_units", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, "degrees_east");
    setAttr(ds_id, "geospatial_lon_units", NC_CHAR, strlen(buf) + 1, buf);
    const char* resString = resolution2str(meta_l3b->resolution);
    int len = strlen(resString) + 1;
    setAttr(ds_id, "geospatial_lat_resolution", NC_CHAR, len , (VOIDP)resString);
    setAttr(ds_id, "geospatial_lon_resolution", NC_CHAR, len , (VOIDP)resString);
    setAttr(ds_id, "spatialResolution", NC_CHAR, len , (VOIDP)resString);

    if(write64bit)
         setAttr(ds_id, "data_bins", NC_INT64, 1, (VOIDP) &meta_l3b->data_bins);
    else {
        int tmpInt = meta_l3b->data_bins;
        setAttr(ds_id, "data_bins", NC_INT, 1, (VOIDP) &tmpInt);
    }

    setAttr(ds_id, "percent_data_bins", NC_FLOAT, 1,
            (VOIDP) & meta_l3b->pct_databins);

    strcpy(buf, meta_l3b->units);
    setAttr(ds_id, "units", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->binning_scheme);
    setAttr(ds_id, "binning_scheme", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, "Ocean Biology Processing Group (NASA/GSFC/OBPG)");
    setAttr(ds_id, "project", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf,
            "NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group");
    setAttr(ds_id, "institution", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, "CF Standard Name Table v36");
    setAttr(ds_id, "standard_name_vocabulary", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, "CF-1.6 ACDD-1.3");
    setAttr(ds_id, "Conventions", NC_CHAR, strlen(buf) + 1, buf);

    char const *naming_authority = "gov.nasa.gsfc.sci.oceandata";
    strcpy(buf, naming_authority);
    setAttr(ds_id, "naming_authority", NC_CHAR, strlen(buf) + 1, buf);
    // create id
    strcpy(buf, meta_l3b->pversion);
    if (strcmp(meta_l3b->pversion, "Unspecified") != 0) {
        strcpy(buf, meta_l3b->product_name);
        strcat(buf, "/L3/");

    } else {
        strcpy(buf, "L3/");
    }
    strcat(buf, meta_l3b->product_name);
    setAttr(ds_id, "id", NC_CHAR, strlen(buf) + 1, buf);

    char const *license =
            "https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/";
    strcpy(buf, license);
    setAttr(ds_id, "license", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, "NASA/GSFC/OBPG");
    setAttr(ds_id, "creator_name", NC_CHAR, strlen(buf) + 1, buf);
    setAttr(ds_id, "publisher_name", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, "data@oceancolor.gsfc.nasa.gov");
    setAttr(ds_id, "creator_email", NC_CHAR, strlen(buf) + 1, buf);
    setAttr(ds_id, "publisher_email", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, "https://oceandata.sci.gsfc.nasa.gov");
    setAttr(ds_id, "creator_url", NC_CHAR, strlen(buf) + 1, buf);
    setAttr(ds_id, "publisher_url", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, "L3 Binned");
    setAttr(ds_id, "processing_level", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, "point");
    setAttr(ds_id, "cdm_data_type", NC_CHAR, strlen(buf) + 1, buf);

    if(meta_l3b->doi[0]) {
        strcpy(buf, "https://dx.doi.org");
        setAttr(ds_id, "identifier_product_doi_authority",
                NC_CHAR, strlen(buf) + 1, buf);
        setAttr(ds_id, "identifier_product_doi", NC_CHAR, strlen(meta_l3b->doi) + 1, meta_l3b->doi);
    }

    if(meta_l3b->keywords[0]) {
        strcpy(buf, "NASA Global Change Master Directory (GCMD) Science Keywords");
        setAttr(ds_id, "keywords_vocabulary", NC_CHAR, strlen(buf) + 1, buf);
        setAttr(ds_id, "keywords", NC_CHAR, strlen(meta_l3b->keywords) + 1, meta_l3b->keywords);
    }

    nc_def_grp(ds_id.fid, "processing_control", &ds_id.fid);

    strcpy(buf, meta_l3b->soft_name);
    setAttr(ds_id, "software_name", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->soft_ver);
    setAttr(ds_id, "software_version", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->infiles);
    setAttr(ds_id, "input_sources", NC_CHAR, strlen(buf) + 1, buf);

    strcpy(buf, meta_l3b->flag_names);
    setAttr(ds_id, "l2_flag_names", NC_CHAR, strlen(buf) + 1, buf);

    nc_def_grp(ds_id.fid, "input_parameters", &ds_id.fid);
    char *cleanInputParams = replace_ocroots(meta_l3b->input_parms);
    char *end_str;
    char *token = strtok_r(cleanInputParams, "\n", &end_str);
    char tmp_buf[16384];
    while (token != NULL) {
        char *end_token;
        strcpy(tmp_buf, token);
        char *name = strtok_r(token, "=", &end_token);
        trimBlanks(name);
        strcpy(tmp_buf, strtok_r(NULL, "\n", &end_token));
        trimBlanks(tmp_buf);
        if (name[0] != '#')
            PTB(SetChrGA(ds_id, name, tmp_buf));
        token = strtok_r(NULL, "\n", &end_str);
    }
    free(cleanInputParams);
    
    return 0;
}
