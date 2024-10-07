#include "viirs_sim_sdr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <time.h>
/*  this is used mainly in setting up the grougs in the output SDRs */
char *core_g_nm[] = {"VIIRS-MOD-GEO-TC", "VIIRS-M1-SDR",
    "VIIRS-M2-SDR", "VIIRS-M3-SDR", "VIIRS-M4-SDR", "VIIRS-M5-SDR",
    "VIIRS-M6-SDR", "VIIRS-M7-SDR", "VIIRS-M8-SDR", "VIIRS-M9-SDR",
    "VIIRS-M10-SDR", "VIIRS-M11-SDR", "VIIRS-M12-SDR", "VIIRS-M13-SDR",
    "VIIRS-M14-SDR", "VIIRS-M15-SDR", "VIIRS-M16-SDR"};

int init_sdr(ctl_struc *ctl, sdr_info_struc *sdr_info, in_rec_struc *in_rec,
        out_rec_struc *out_rec)
/*-----------------------------------------------------------------------------
    Program:   init_sdr.c

    Description:  create the geolocation and band files with everything but 
      data lines

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        ctl_struc *  ctl        I    input controls
        sdr_info_struc * sdr_info I  general SDR information
        in_rec_struc * in_rec    I/O  controls for input record reading
        out_rec_struc *  out_rec  I/O  controls for output file writing

    Modification history:

    W. Robinson, SAIC  15 Oct 2008  Original development

----------------------------------------------------------------------------*/ {
    h5io_str dat1_g_id, dat2_g_id;
    char out_file[500], grp_nam[50], *sdr_base;
    int isdr;
    char lcl_out_bnd_typ[] = {0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0};
    char lcl_meas_typ[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
    /*
     *  set the radiance storage and type flag
     */
    for (isdr = 0; isdr < MAX_BND; isdr++) {
        out_rec->out_bnd_typ[isdr] = lcl_out_bnd_typ[isdr];
        out_rec->meas_typ[isdr] = lcl_meas_typ[isdr];
    }
    /*
     *  initialize the lat, lon limit values
     */
    in_rec->ll_lims[0] = 99.;
    in_rec->ll_lims[1] = -99.;
    in_rec->ll_lims[2] = 200.;
    in_rec->ll_lims[3] = -200.;
    in_rec->ll_lims[4] = 200.;
    in_rec->ll_lims[5] = -200.;
    /*
     *  based on the controls and the input scan format, set controls for
     *  output scan format information
     *  base state is to leave as-is (out_scn_fmt = 2)
     */
    out_rec->nbnd = in_rec->nbnd;
    out_rec->npix = in_rec->npix;
    out_rec->nlin = in_rec->nlin;
    out_rec->nscan = in_rec->nscan;
    out_rec->scn_fmt = in_rec->scn_fmt;
    out_rec->margin[0] = in_rec->margin[0];
    out_rec->margin[1] = in_rec->margin[1];
    out_rec->ndet_scan = in_rec->ndet_scan;

    if (ctl->out_scn_fmt == 0) /* make aggregated */ {
        out_rec->scn_fmt = 0;
        out_rec->npix = 3200;
        out_rec->nlin = in_rec->nscan * NDET;
        out_rec->margin[0] = 0;
        out_rec->margin[1] = 0;
        out_rec->ndet_scan = NDET;
    } else if (ctl->out_scn_fmt == 1) /* make unaggregated with no margin */ {
        out_rec->margin[0] = 0;
        out_rec->margin[1] = 0;
        out_rec->ndet_scan = NDET;
        out_rec->nlin = in_rec->nscan * NDET;

        if (in_rec->scn_fmt == 0) {
            printf(
                    "%s, %d: Warning - requested output format (%d) of unaggregated\n",
                    __FILE__, __LINE__, in_rec->scn_fmt);
            printf("          will default to input scan format of aggregated\n");
            out_rec->npix = 3200;
            out_rec->scn_fmt = 0;
        } else {
            out_rec->npix = 6304;
            out_rec->scn_fmt = 1;
        }
    }
    /*
     *  go through the SDR files and initialize each one
     */
    for (isdr = 0; isdr < in_rec->nbnd + 1; isdr++) {
        /*
         *  divide the work among creating the major parts of the file
         *  first, open output file and set up the top attributes
         */
        if (gen_sdr_fname(isdr, ctl->out_loc, sdr_info, ctl->fname_opt,
                out_file) != 0) {
            printf("%s, %d: Unable to generate output file for sdr index %d\n\n",
                    __FILE__, __LINE__, isdr);
            return 1;
        }
        sdr_base = basename(out_file);
        strcpy(sdr_info->sdr_files[isdr], sdr_base);
        if (h5io_openw(out_file, ctl->sdr_overwrite,
                &(out_rec->sdr_fid[isdr])) != 0) {
            printf("%s, %d: Unable to open output geo file: \n%s\n\n",
                    __FILE__, __LINE__, out_file);
            return 1;
        }
        if (init_sdr_top(isdr, sdr_info, out_rec) != 0) return 1;
        /*
         *  and then make the data group and write the data
         */
        if (h5io_mk_grp(&(out_rec->sdr_fid[isdr]), "All_Data",
                &(out_rec->sdr_dat_gid[0][isdr])) != 0) {
            printf("%s, %d - could not create data group: All_Data\n",
                    __FILE__, __LINE__);
            return 1;
        }
        sprintf(grp_nam, "%s_All", core_g_nm[isdr]);
        if (h5io_mk_grp(&(out_rec->sdr_dat_gid[0][isdr]), grp_nam,
                &(out_rec->sdr_dat_gid[1][isdr])) != 0) {
            printf("%s, %d - could not create data group: %s\n",
                    __FILE__, __LINE__, grp_nam);
            return 1;
        }
        if (isdr == 0) {
            if (init_geo_data(sdr_info, in_rec, out_rec) != 0) return 1;
        } else {
            /*
             *  allocate the dn storage here (depends on ctl) but do rest in 
             *  call to init_bnd_data
             */
            if (ctl->count_cal_opt != 0) {
                if ((in_rec->dn[ isdr - 1 ] = (float *)
                        malloc(in_rec->ndet_scan * in_rec->npix * sizeof (float)))
                        == NULL) {
                    printf("%s, %d: Error, allocation of dn storage failed\n",
                            __FILE__, __LINE__);
                    return 1;
                }
            }
            if ((in_rec->gain_bit[ isdr - 1 ] = (unsigned char *)
                    calloc(in_rec->ndet_scan * in_rec->npix, sizeof (unsigned char)))
                    == NULL) {
                printf(
                        "%s, %d: Error, allocation of count gain bit storage failed\n",
                        __FILE__, __LINE__);
                return 1;
            }
            if ((in_rec->dn_sat[ isdr - 1 ] = (char *)
                    calloc(in_rec->ndet_scan * in_rec->npix, sizeof (char)))
                    == NULL) {
                printf(
                        "%s, %d: Error, allocation of dn saturation storage failed\n",
                        __FILE__, __LINE__);
                return 1;
            }
            if (init_bnd_data(isdr - 1, sdr_info, in_rec, out_rec) != 0) return 1;
        }
        /*
         *  make the other group pair for data products and fill
         */
        if (h5io_mk_grp(&(out_rec->sdr_fid[isdr]),
                "Data_Products", &dat1_g_id) != 0) {
            printf("%s, %d - could not create data group: Data_Products\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if (h5io_mk_grp(&dat1_g_id, core_g_nm[isdr], &dat2_g_id) != 0) {
            printf("%s, %d - could not create data group: %s\n",
                    __FILE__, __LINE__, core_g_nm[isdr]);
            return 1;
        }
        /*
         *  set up attributes in the group, the aggregate and granule datasets
         *  and their attributes
         */
        if (init_sdr_dpattr(isdr, &dat2_g_id, sdr_info) != 0) return 1;
        if (init_sdr_agg(isdr, &dat2_g_id, sdr_info) != 0) return 1;
        if (init_sdr_gran(isdr, &dat2_g_id, sdr_info, out_rec) != 0) return 1;
        /*
         *  and finish up, closing the group ids
         */
        if (h5io_close(&dat2_g_id) != 0) {
            printf("%s, %d - could not close dat2_g_id\n", __FILE__, __LINE__);
            return 1;
        }
        if (h5io_close(&dat1_g_id) != 0) {
            printf("%s, %d - could not close dat1_g_id\n", __FILE__, __LINE__);
            return 1;
        }
    }
    /*
     *  we still need to place the data lines
     */
    return 0;
}

int init_sdr_top(int isdr, sdr_info_struc *sdr_info, out_rec_struc *out_rec)
/*-----------------------------------------------------------------------------
    Routine:   init_sdr_top

    Description:  make the top-level attributes for the geolocation or 
      band SDR file

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       isdr          I    SDR file to work on: 0 geo, > 0 band isdr
        sdr_info_struc * sdr_info I  general SDR information
        out_rec_struc * out_rec I    output dataset information

    Modification history:

    W. Robinson, SAIC  15 Oct 2008  Original development
    W. Robinson, SAIC  18 Mar 2010  place non-standard attributes here for the 
       non-aggregated file types

----------------------------------------------------------------------------*/ {
    int i, n_attr = 13, dims_1_1[] = {1, 1}, len_geo;
    int dims_1[] = {1}, dims_2[] = {2};
    char geo_name[150], *bloc = ".";
    int16_t lcl_scn_fmt, lcl_margin[2], lcl_ndet_scan;
    int32_t lcl_npix, lcl_nlin;

    gen_sdr_fname(0, bloc, sdr_info, 0, geo_name);
    len_geo = strlen(geo_name);

    h5attr_struc attrs[] = {
        { 1, 1, "Distributor", H5T_NATIVE_CHAR, 4, 2, dims_1_1,
            (void *) (sdr_info->origin)},
        { 1, 1, "Instrument_Short_Name", H5T_NATIVE_CHAR, 6, 2, dims_1_1,
            (void *) "VIIRS"},
        { 1, 1, "Mission_Name", H5T_NATIVE_CHAR, 4, 2, dims_1_1,
            (void *) "NPP"},
        { 1, 1, "N_Dataset_Source", H5T_NATIVE_CHAR, 5, 2, dims_1_1,
            (void *) "OBPG"},
        { 1, 1, "N_GEO_Ref", H5T_NATIVE_CHAR, len_geo, 2, dims_1_1,
            (void *) geo_name},
        { 1, 1, "N_HDF_Creation_Date", H5T_NATIVE_CHAR, 9, 2, dims_1_1,
            (void *) sdr_info->cre_date},
        { 1, 1, "N_HDF_Creation_Time", H5T_NATIVE_CHAR, 14, 2, dims_1_1,
            (void *) sdr_info->cre_time},
        { 1, 1, "Platform_Short_Name", H5T_NATIVE_CHAR, 4, 2, dims_1_1,
            (void *) "NPP"},
        { 0, 0, "Data Scan Format", H5T_STD_I16BE, 0, 1, dims_1,
            (void *) &lcl_scn_fmt},
        { 0, 0, "Scan Margin (track, scan)", H5T_STD_I16BE, 0, 1, dims_2,
            (void *) lcl_margin},
        { 0, 0, "Pixels per Scan Line", H5T_STD_I32BE, 0, 1, dims_1,
            (void *) &lcl_npix},
        { 0, 0, "Number of Scan Lines", H5T_STD_I32BE, 0, 1, dims_1,
            (void *) &lcl_nlin},
        { 0, 0, "Number of Detectors per Scan", H5T_STD_I16BE, 0, 1, dims_1,
            (void *) &lcl_ndet_scan}
    };
    /*
     *  for the geo file, don't output the N_GEO_Ref attrib
     */
    if (isdr == 0) attrs[4].express = 0;
    /*
     *  for non-aggregated file formats, add the descriptors
     */
    if (out_rec->scn_fmt != 0) {
        for (i = 0; i < 5; i++)
            attrs[ i + 8 ].express = 1;
        lcl_scn_fmt = (int16_t) out_rec->scn_fmt;
        lcl_margin[0] = (int16_t) out_rec->margin[0];
        lcl_margin[1] = (int16_t) out_rec->margin[1];
        lcl_npix = (int32_t) out_rec->npix;
        lcl_nlin = (int32_t) out_rec->nlin;
        lcl_ndet_scan = (int16_t) out_rec->ndet_scan;
    }
    /*
     *  output the attributes to the location
     */
    if (wr_attr_seq(&(out_rec->sdr_fid[isdr]), n_attr, attrs) != 0) {
        printf("%s, %d: Could not write top attributes\n", __FILE__, __LINE__);
        return 1;
    }
    return 0;
}

int init_geo_data(sdr_info_struc *sdr_info, in_rec_struc *in_rec,
        out_rec_struc *out_rec)
/*-----------------------------------------------------------------------------
    Routine:   init_geo_data

    Description:  make the data arrays for the geo file and fill some
      Also, get some data start, end times

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        sdr_info_struc * sdr_info I/O  general SDR information
        in_rec_struc * in_rec  I/O  input file information
        out_rec_struc * out_rec I/O  output file information

    Modification history:

    W. Robinson, SAIC  21 Oct 2008  Original development
    W. Robinson, SAIC  02 Oct 2009  use geoloc file attitude, position, 
                                    velocity in the SDR

----------------------------------------------------------------------------*/ {
    h5io_str ds_id, in_ds_id, *gid;
    int dim_siz[2], dim_siz2[2], *arr_int, i;
    float *flt_data;
    double *dbl_data;
    int64 *llon_data;
    unsigned char *uchar_data;
    /*
     *  define gid for convenience
     */
    gid = &(out_rec->sdr_dat_gid[1][0]);
    /*
     *  datasets are:
     *  Height -- may be derivable but not necessary
     */
    dim_siz[0] = out_rec->nlin;
    dim_siz[1] = out_rec->npix;

    if ((flt_data = (float *) malloc(out_rec->npix * out_rec->nlin *
            sizeof (float))) == NULL) {
        printf("%s, %d: Could not allocate local float buffer\n",
                __FILE__, __LINE__);
        return 1;
    }

    for (i = 0; i < out_rec->nlin * out_rec->npix; i++)
        *(flt_data + i) = 0;
    if (h5io_mk_ds(gid, "Height", H5T_IEEE_F32BE, 2, dim_siz,
            &ds_id) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for Height\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) flt_data) != 0) {
        printf("%s, %d: Could not write to Height\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close Height\n", __FILE__, __LINE__);
        return 1;
    }
    /*
     *  Latitude
     *  The lat, lon, senz, sena, solz, sola arrays will be filled scae-by-scan
     *  later.  just create these datasets here (others will be filled here)
     */
    if (h5io_set_ds(&(in_rec->geo_fid), "latitude",
            &(in_rec->geo_dat_id[0])) != 0) {
        printf("%s, %d: Could not do h5io_set_ds for latitude\n",
                __FILE__, __LINE__);
        return 1;
    }

    dim_siz[0] = out_rec->nlin;
    dim_siz[1] = out_rec->npix;
    if (h5io_mk_ds(gid, "Latitude", H5T_IEEE_F32BE, 2, dim_siz,
            &(out_rec->geo_dat_id[0])) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for Latitude\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  allocate the data transfer buffer and make it the same for output
     **** good for now as long as out size = in size but will need change 
            if scn_fnt is changed
     */
    if ((in_rec->lat = (float *) malloc(in_rec->npix * in_rec->ndet_scan *
            sizeof (float))) == NULL) {
        printf("%s, %d: Could not allocate space for lat data transfer buffer\n",
                __FILE__, __LINE__);
        return 1;
    }
    out_rec->lat = in_rec->lat;
    /*
     *  Longitude
     */
    if (h5io_set_ds(&(in_rec->geo_fid), "longitude",
            &(in_rec->geo_dat_id[1])) != 0) {
        printf("%s, %d: Could not do h5io_set_ds for longitude\n",
                __FILE__, __LINE__);
        return 1;
    }

    if (h5io_mk_ds(gid, "Longitude", H5T_IEEE_F32BE, 2, dim_siz,
            &(out_rec->geo_dat_id[1])) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for Longitude\n",
                __FILE__, __LINE__);
        return 1;
    }

    if ((in_rec->lon = (float *) malloc(in_rec->npix * in_rec->ndet_scan *
            sizeof (float))) == NULL) {
        printf("%s, %d: Could not allocate space for lon data transfer buffer\n",
                __FILE__, __LINE__);
        return 1;
    }
    out_rec->lon = in_rec->lon;
    /*
     *  MidTime - an 8 byte int of microsecs past 1/1958 -- for now,
     *  fill the ScanStartTime but leave this 0
     */
    llon_data = (int64 *) calloc(out_rec->nscan, sizeof (int64));
    dim_siz2[0] = out_rec->nscan;
    if (h5io_mk_ds(gid, "MidTime", H5T_STD_I64BE, 1, dim_siz2, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for MidTime\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) llon_data) != 0) {
        printf("%s, %d: Could not write to MidTime\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close MidTime\n", __FILE__, __LINE__);
        return 1;
    }
    free(llon_data);
    /*
     *  ModeGran 1 value, 1 is day, so fill with 1  
     */
    uchar_data = (unsigned char *) malloc(out_rec->nscan *
            sizeof ( unsigned char));
    *uchar_data = 1;
    dim_siz2[0] = 1;
    if (h5io_mk_ds(gid, "ModeGran", H5T_STD_U8BE, 1, dim_siz2, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for ModeGran\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) uchar_data) != 0) {
        printf("%s, %d: Could not write to ModeGran\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close ModeGran\n", __FILE__, __LINE__);
        return 1;
    }
    free(uchar_data);
    /*
     *  ModeScan 48 bytes, 1 is day, so fill with 1 
     */
    uchar_data = (unsigned char *) malloc(out_rec->nscan *
            sizeof ( unsigned char));
    for (i = 0; i < out_rec->nscan; i++)
        *(uchar_data + i) = 1;
    dim_siz2[0] = out_rec->nscan;
    if (h5io_mk_ds(gid, "ModeScan", H5T_STD_U8BE, 1, dim_siz2, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for ModeScan\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) uchar_data) != 0) {
        printf("%s, %d: Could not write to ModeScan\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close ModeScan\n", __FILE__, __LINE__);
        return 1;
    }
    free(uchar_data);
    /*
     *  NumberOfScans
     */
    arr_int = (int *) malloc(sizeof ( int));
    *arr_int = out_rec->nscan;
    dim_siz2[0] = 1;
    if (h5io_mk_ds(gid, "NumberOfScans", H5T_STD_I32BE, 1, dim_siz2, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for NumberOfScans\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) arr_int) != 0) {
        printf("%s, %d: Could not write to NumberOfScans\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close NumberOfScans\n", __FILE__, __LINE__);
        return 1;
    }
    free(arr_int);
    /*
     *  PadByte1 ( 3 values of 0)  -- just fill with 0 for the 3 values
     */
    uchar_data = (unsigned char *) malloc(3 * sizeof ( unsigned char));
    for (i = 0; i < 3; i++)
        *(uchar_data + i) = 0;
    dim_siz2[0] = 3;
    if (h5io_mk_ds(gid, "PadByte1", H5T_STD_U8BE, 1, dim_siz2, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for PadByte\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) uchar_data) != 0) {
        printf("%s, %d: Could not write to PadByte\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close PadByte\n", __FILE__, __LINE__);
        return 1;
    }
    free(uchar_data);
    /*
     *  QF1_SCAN_VIIRSSDRGEO - contains info on HAM encoder qual and attitude,
     *  ephem availability good is  48 * 0
     */
    uchar_data = (unsigned char *) malloc(out_rec->nscan *
            sizeof ( unsigned char));
    for (i = 0; i < out_rec->nscan; i++)
        *(uchar_data + i) = 0;
    dim_siz2[0] = out_rec->nscan;
    if (h5io_mk_ds(gid, "QF1_SCAN_VIIRSSDRGEO", H5T_STD_U8BE, 1,
            dim_siz2, &ds_id) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for QF1_SCAN_VIIRSSDRGEO\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) uchar_data) != 0) {
        printf("%s, %d: Could not write to QF1_SCAN_VIIRSSDRGEO\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close QF1_SCAN_VIIRSSDRGEO\n",
                __FILE__, __LINE__);
        return 1;
    }
    free(uchar_data);
    /*
     *  QF2_VIIRSSDRGEO - valid / invalid state for 4 items: input data,
     *  pointing, terrain and solar angles - all 0 for all pixels is all good
     */
    uchar_data = (unsigned char *) calloc(out_rec->npix * out_rec->nlin,
            sizeof ( unsigned char));
    if (h5io_mk_ds(gid, "QF2_VIIRSSDRGEO", H5T_STD_U8BE, 2, dim_siz, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for QF2_VIIRSSDRGEO\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) uchar_data) != 0) {
        printf("%s, %d: Could not write to QF2_VIIRSSDRGEO\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close QF2_VIIRSSDRGEO\n", __FILE__, __LINE__);
        return 1;
    }
    free(uchar_data);
    /*
     *  SCAttitude attitude information 48 X 3 from geoloc dataset
     *  convert from deg to arcsec
     */
    for (i = 0; i < out_rec->nscan * 3; i++)
        *(flt_data + i) = *(sdr_info->geo_att + i) * 3600.;
    /*  free the space as it is no longer in use */
    free(sdr_info->geo_att);
    dim_siz2[0] = out_rec->nscan;
    dim_siz2[1] = 3;
    if (h5io_mk_ds(gid, "SCAttitude", H5T_IEEE_F32BE, 2, dim_siz2,
            &ds_id) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for SCAttitude\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) flt_data) != 0) {
        printf("%s, %d: Could not write to SCAttitude\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close SCAttitude\n", __FILE__, __LINE__);
        return 1;
    }
    /*
     *  SCPosition SC position, ECR coords., from geoloc dataset
     */
    for (i = 0; i < out_rec->nscan * 3; i++)
        *(flt_data + i) = *(sdr_info->geo_pos + i) * 1000.;
    /*  free the space as it is no longer in use */
    free(sdr_info->geo_pos);
    dim_siz2[0] = out_rec->nscan;
    dim_siz2[1] = 3;
    if (h5io_mk_ds(gid, "SCPosition", H5T_IEEE_F32BE, 2, dim_siz2,
            &ds_id) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for SCPosition\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) flt_data) != 0) {
        printf("%s, %d: Could not write to SCPosition\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close SCPosition\n", __FILE__, __LINE__);
        return 1;
    }
    /*
     *  SCSolarAzimuthAngle - 48 array of az angle on solar diffuser
     */
    for (i = 0; i < out_rec->nscan; i++)
        *(flt_data + i) = 0.;
    dim_siz2[0] = out_rec->nscan;
    if (h5io_mk_ds(gid, "SCSolarAzimuthAngle", H5T_IEEE_F32BE, 1, dim_siz2,
            &ds_id) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for SCSolarAzimuthAngle\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) flt_data) != 0) {
        printf("%s, %d: Could not write to SCSolarAzimuthAngle\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close SCSolarAzimuthAngle\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  SCSolarZenithAngle - a 48 array of zen angle on solar diffuser
     */
    for (i = 0; i < out_rec->nscan; i++)
        *(flt_data + i) = 0.;
    dim_siz2[0] = out_rec->nscan;
    if (h5io_mk_ds(gid, "SCSolarZenithAngle", H5T_IEEE_F32BE, 1, dim_siz2,
            &ds_id) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for SCSolarZenithAngle\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) flt_data) != 0) {
        printf("%s, %d: Could not write to SCSolarZenithAngle\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close SCSolarZenithAngle\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  SCVelocity -- get from geoloc dataset
     */
    dim_siz2[0] = out_rec->nscan;
    dim_siz2[1] = 3;
    for (i = 0; i < out_rec->nscan * 3; i++)
        *(flt_data + i) = *(sdr_info->geo_vel + i) * 1000.;
    /*  free the space as it is no longer in use */
    free(sdr_info->geo_vel);
    if (h5io_mk_ds(gid, "SCVelocity", H5T_IEEE_F32BE, 2, dim_siz2,
            &ds_id) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for SCVelocity\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) flt_data) != 0) {
        printf("%s, %d: Could not write to SCVelocity\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close SCVelocity\n", __FILE__, __LINE__);
        return 1;
    }
    /*
     *  SatelliteAzimuthAngle
     */
    if (h5io_set_ds(&(in_rec->geo_fid), "sena",
            &(in_rec->geo_dat_id[2])) != 0) {
        printf("%s, %d: Could not do h5io_set_ds for sena\n", __FILE__, __LINE__);
        return 1;
    }

    if (h5io_mk_ds(gid, "SatelliteAzimuthAngle", H5T_IEEE_F32BE, 2,
            dim_siz, &(out_rec->geo_dat_id[2])) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for SatelliteAzimuthAngle\n",
                __FILE__, __LINE__);
        return 1;
    }

    if ((in_rec->sena = (float *) malloc(in_rec->npix * in_rec->ndet_scan *
            sizeof (float))) == NULL) {
        printf(
                "%s, %d: Could not allocate transfer buffer for SatelliteAzimuthAngle\n",
                __FILE__, __LINE__);
        return 1;
    }
    out_rec->sena = in_rec->sena;
    /*
     *  SatelliteRange -- We may be able to get, but not necessary now
     *    Just zero out
     */
    for (i = 0; i < out_rec->npix * out_rec->nlin; i++)
        *(flt_data + i) = 0.;
    if (h5io_mk_ds(gid, "SatelliteRange", H5T_IEEE_F32BE, 2, dim_siz, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for SatelliteRange\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) flt_data) != 0) {
        printf("%s, %d: Could not write to SatelliteRange\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close SatelliteRange\n", __FILE__, __LINE__);
        return 1;
    }
    /*
     *  SatelliteZenithAngle
     */
    if (h5io_set_ds(&(in_rec->geo_fid), "senz",
            &(in_rec->geo_dat_id[3])) != 0) {
        printf("%s, %d: Could not do h5io_set_ds for senz\n", __FILE__, __LINE__);
        return 1;
    }

    if (h5io_mk_ds(gid, "SatelliteZenithAngle", H5T_IEEE_F32BE, 2,
            dim_siz, &(out_rec->geo_dat_id[3])) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for SatelliteZenithAngle\n",
                __FILE__, __LINE__);
        return 1;
    }

    if ((in_rec->senz = (float *) malloc(in_rec->npix * in_rec->ndet_scan *
            sizeof (float))) == NULL) {
        printf(
                "%s, %d: Could not allocate transfer buffer for SatelliteZenithAngle\n",
                __FILE__, __LINE__);
        return 1;
    }
    out_rec->senz = in_rec->senz;
    /*
     *  SolarAzimuthAngle
     */
    if (h5io_set_ds(&(in_rec->geo_fid), "sola",
            &(in_rec->geo_dat_id[4])) != 0) {
        printf("%s, %d: Could not do h5io_set_ds on sola\n", __FILE__, __LINE__);
        return 1;
    }

    if (h5io_mk_ds(gid, "SolarAzimuthAngle", H5T_IEEE_F32BE, 2, dim_siz,
            &(out_rec->geo_dat_id[4])) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for SolarAzimuthAngle\n",
                __FILE__, __LINE__);
        return 1;
    }

    if ((in_rec->solz = (float *) malloc(in_rec->npix * in_rec->ndet_scan *
            sizeof (float))) == NULL) {
        printf(
                "%s, %d: Could not allocate transfer buffer for SolarAzimuthAngle\n",
                __FILE__, __LINE__);
        return 1;
    }
    out_rec->solz = in_rec->solz;
    /*
     *  SolarZenithAngle
     */
    if (h5io_set_ds(&(in_rec->geo_fid), "solz",
            &(in_rec->geo_dat_id[5])) != 0) {
        printf("%s, %d: Could not do h5io_set_ds on solz\n", __FILE__, __LINE__);
        return 1;
    }

    if (h5io_mk_ds(gid, "SolarZenithAngle", H5T_IEEE_F32BE, 2, dim_siz,
            &(out_rec->geo_dat_id[5])) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for SolarZenithAngle\n",
                __FILE__, __LINE__);
        return 1;
    }

    if ((in_rec->sola = (float *) malloc(in_rec->npix * in_rec->ndet_scan *
            sizeof (float))) == NULL) {
        printf("%s, %d: Could not allocate transfer buffer for SolarZenithAngle\n",
                __FILE__, __LINE__);
        return 1;
    }
    out_rec->sola = in_rec->sola;
    /*
     *  StartTime, again the microsecs
     */
    dbl_data = (double *) malloc(out_rec->nscan * sizeof (double));
    if (h5io_set_ds(&(in_rec->geo_fid), "time", &in_ds_id) != 0) {
        printf("%s, %d: Could not do h5io_set_ds for time\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_rd_ds(&in_ds_id, (void *) dbl_data) != 0) {
        printf("%s, %d: Could not read time\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&in_ds_id) != 0) {
        printf("%s, %d: Could not close time\n", __FILE__, __LINE__);
        return 1;
    }

    llon_data = (int64 *) malloc(out_rec->nscan * sizeof (int64));
    for (i = 0; i < out_rec->nscan; i++)
        *(llon_data + i) = sdr_info->t58_day +
            (int64) (*(dbl_data + i) * 1e6);
    dim_siz2[0] = out_rec->nscan;
    if (h5io_mk_ds(gid, "StartTime", H5T_STD_I64BE, 1, dim_siz2, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for StartTime\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) llon_data) != 0) {
        printf("%s, %d: Could not write StartTime\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close StartTime\n", __FILE__, __LINE__);
        return 1;
    }

    free(flt_data);
    return 0;
}

int init_bnd_data(int ibnd, sdr_info_struc *sdr_info, in_rec_struc *in_rec,
        out_rec_struc *out_rec)
/*-----------------------------------------------------------------------------
    Routine:   init_bnd_data

    Description:  make the data arrays for the band files and fill some

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       ibnd          I    band index 0 is band M1...
        sdr_info_struc * sdr_info I/O  general SDR information
        in_rec_struc * in_rec  I/O  input file information
        out_rec_struc * out_rec I/O  output file information

    Modification history:

    W. Robinson, SAIC  9 Dec 2008  Original development
    W. Robinson, SAIC  7 Jul 2010  add ability to make all 16 M bands
    W. Robinson, SAIC  24 Feb 2011 switch scale, offset values to those used
                          by IDPS (keep old method commented just in case)

----------------------------------------------------------------------------*/ {
    h5io_str ds_id, *gid;
    int dim_siz[2], dim_siz2[2], *arr_int, i, isdr, iret;
    static unsigned char m_side = 0;
    unsigned char *uchar_data;
    unsigned int *uint32_data;
    float rad_fac[2], refl_fac[2]; /* versions of scale, offset for sdr
      scaled values */
    /*  float lam_um;  --- not used currently, but may come back */
    /*  these values are scale, offset used/found in IDPS SDR files.  Some 
        match my F0 assumptions well, but just use IDPS to avoid confusion 
        Used for RadianceFactors, ReflectanceFactors and 
        BrightnessTemperatureFactors */
    static float lcl_scale[] = {0.0112657, 0.0125841, 1., 1., 1., 0.000752209,
        1., 0.00302196, 0.00141331, 0.00130450, 0.000582661, 5.17344e-05,
        1., 0.000321547, 0.000313153, 0.000265539};
    static float lcl_offset[] = {-0.21, -0.2, 0., 0., 0., -0.09, 0., -0.14,
        -0.09, -0.04, -0.02, 0., 0., -0.03, -0.02, -0.02};
    static float lcl_scale_ref[] = {2.28913e-05, 2.28913e-05, 2.28913e-05,
        2.28913e-05, 2.28913e-05, 2.28913e-05, 2.28913e-05, 2.28913e-05,
        2.28913e-05, 2.28913e-05, 2.28913e-05, 0.00251805, 1., 0.00373892,
        0.00412044, 0.00425779};
    static float lcl_offset_ref[] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 203., 0., 120., 111., 103.};
    /*
     *  define gid for convenience
     */
    m_side = 0;
    isdr = ibnd + 1; /*  isdr is also the commonly known band # M1, ...  */
    gid = &(out_rec->sdr_dat_gid[1][isdr]);
    /*
     *  datasets are:
     *  ModeGran 1 value, 1 is day, so fill with 1
     */
    uchar_data = (unsigned char *) malloc(out_rec->nscan *
            sizeof ( unsigned char));
    *uchar_data = 1;
    dim_siz2[0] = 1;
    if (h5io_mk_ds(gid, "ModeGran", H5T_STD_U8BE, 1, dim_siz2, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for ModeGran\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) uchar_data) != 0) {
        printf("%s, %d: Could not write to ModeGran\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close ModeGran\n", __FILE__, __LINE__);
        return 1;
    }
    free(uchar_data);
    /*
     *  ModeScan 48 bytes, 1 is day, so fill with 1
     */
    uchar_data = (unsigned char *) malloc(out_rec->nscan *
            sizeof ( unsigned char));
    for (i = 0; i < out_rec->nscan; i++)
        *(uchar_data + i) = 1;
    dim_siz2[0] = out_rec->nscan;
    if (h5io_mk_ds(gid, "ModeScan", H5T_STD_U8BE, 1, dim_siz2, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for ModeScan\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) uchar_data) != 0) {
        printf("%s, %d: Could not write to ModeScan\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close ModeScan\n", __FILE__, __LINE__);
        return 1;
    }
    free(uchar_data);
    /*
     *  NumberOfBadChecksums
     */
    arr_int = (int *) malloc(out_rec->nscan * sizeof ( int));
    for (i = 0; i < out_rec->nscan; i++)
        *(arr_int + i) = 1;
    dim_siz2[0] = out_rec->nscan;
    if (h5io_mk_ds(gid, "NumberOfBadChecksums", H5T_STD_I32BE, 1,
            dim_siz2, &ds_id) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for NumberOfBadChecksums\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) arr_int) != 0) {
        printf("%s, %d: Could not write to NumberOfBadChecksums\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close NumberOfBadChecksums\n",
                __FILE__, __LINE__);
        return 1;
    }
    free(arr_int);
    /*
     *  NumberOfDiscardedPackets
     */
    arr_int = (int *) malloc(out_rec->nscan * sizeof ( int));
    for (i = 0; i < out_rec->nscan; i++)
        *(arr_int + i) = 1;
    dim_siz2[0] = out_rec->nscan;
    if (h5io_mk_ds(gid, "NumberOfDiscardedPackets", H5T_STD_I32BE, 1,
            dim_siz2, &ds_id) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for NumberOfDiscardedPackets\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) arr_int) != 0) {
        printf("%s, %d: Could not write to NumberOfDiscardedPackets\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close NumberOfDiscardedPackets\n",
                __FILE__, __LINE__);
        return 1;
    }
    free(arr_int);
    /*
     *  NumberOfMissingPkts
     */
    arr_int = (int *) malloc(out_rec->nscan * sizeof ( int));
    for (i = 0; i < out_rec->nscan; i++)
        *(arr_int + i) = 1;
    dim_siz2[0] = out_rec->nscan;
    if (h5io_mk_ds(gid, "NumberOfMissingPkts", H5T_STD_I32BE, 1,
            dim_siz2, &ds_id) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for NumberOfMissingPkts\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) arr_int) != 0) {
        printf("%s, %d: Could not write to NumberOfMissingPkts\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close NumberOfMissingPkts\n",
                __FILE__, __LINE__);
        return 1;
    }
    free(arr_int);
    /*
     *  NumberOfScans
     */
    arr_int = (int *) malloc(sizeof ( int));
    *arr_int = out_rec->nscan;
    dim_siz[0] = 1;
    if (h5io_mk_ds(gid, "NumberOfScans", H5T_STD_I32BE, 1, dim_siz, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for NumberOfScans\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) arr_int) != 0) {
        printf("%s, %d: Could not write to NumberOfScans\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close NumberOfScans\n", __FILE__, __LINE__);
        return 1;
    }
    free(arr_int);
    /*
     *  PadByte1 ( 3 values of 0)  -- just fill with 0 for the 3 values
     */
    uchar_data = (unsigned char *) malloc(3 * sizeof ( unsigned char));
    for (i = 0; i < 3; i++)
        *(uchar_data + i) = 0;
    dim_siz[0] = 3;
    if (h5io_mk_ds(gid, "PadByte1", H5T_STD_U8BE, 1, dim_siz, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for PadByte1\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) uchar_data) != 0) {
        printf("%s, %d: Could not write to PadByte1\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close PadByte1\n", __FILE__, __LINE__);
        return 1;
    }
    free(uchar_data);
    /*
     *  QF1_VIIRSMBANDSDR  Qual flag for all pixels
     *  set up dataset ids for each band
     */
    dim_siz[0] = out_rec->nlin;
    dim_siz[1] = out_rec->npix;
    if (h5io_mk_ds(gid, "QF1_VIIRSMBANDSDR", H5T_STD_U8BE, 2, dim_siz,
            &(out_rec->qual1_m_id[ibnd])) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for QF1_VIIRSMBANDSDR\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  and set a 0 filled default array
     */
    if ((out_rec->qual1_m[ibnd] = (unsigned char *)
            calloc(out_rec->npix * out_rec->ndet_scan, sizeof ( unsigned char)))
            == NULL) {
        printf("%s, %d: Could not allocate space for qual1_m output for band %d\n",
                __FILE__, __LINE__, ibnd);
        return 1;
    }
    /*
     *  QF2_SCAN_SDR has more qual plus HAM mirror side in low bit
     *  also set it in the sdr_info ham_side array
     */
    uchar_data = (unsigned char *) calloc(out_rec->nscan,
            sizeof ( unsigned char));
    for (i = 0; i < out_rec->nscan; i++) {
        *(uchar_data + i) = m_side;
        m_side = 1 - m_side;
        *(sdr_info->ham_side + i) = m_side;
    }
    dim_siz[0] = out_rec->nscan;
    if (h5io_mk_ds(gid, "QF2_SCAN_SDR", H5T_STD_U8BE, 1, dim_siz, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for QF2_SCAN_SDR\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) uchar_data) != 0) {
        printf("%s, %d: Could not write to QF2_SCAN_SDR\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close QF2_SCAN_SDR\n", __FILE__, __LINE__);
        return 1;
    }
    free(uchar_data);
    /*
     *  QF3_SCAN_RDR  a nscan x 4 value array
     */
    dim_siz[0] = out_rec->nscan;
    dim_siz[1] = 4;
    uint32_data = (unsigned int *)
            calloc(dim_siz[0] * dim_siz[1], sizeof ( unsigned int));
    if (h5io_mk_ds(gid, "QF3_SCAN_RDR", H5T_STD_U32BE, 2, dim_siz, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for QF3_SCAN_RDR\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) uint32_data) != 0) {
        printf("%s, %d: Could not write to QF3_SCAN_RDR\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close QF3_SCAN_RDR\n", __FILE__, __LINE__);
        return 1;
    }
    free(uint32_data);
    /*
     * QF4_SCAN_SDR 768 values reduced line quality # lines in size
     */
    uchar_data = (unsigned char *) calloc(out_rec->nlin,
            sizeof ( unsigned char));
    dim_siz2[0] = out_rec->nlin;
    if (h5io_mk_ds(gid, "QF4_SCAN_SDR", H5T_STD_U8BE, 1, dim_siz2, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for QF4_SCAN_SDR\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) uchar_data) != 0) {
        printf("%s, %d: Could not write to QF4_SCAN_SDR\n", __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close QF4_SCAN_SDR\n", __FILE__, __LINE__);
        return 1;
    }
    free(uchar_data);
    /*
     * QF5_GRAN_BADDETECTOR 16 values 0 if not any bad detectors
     */
    uchar_data = (unsigned char *) calloc(out_rec->ndet_scan,
            sizeof ( unsigned char));
    dim_siz2[0] = out_rec->ndet_scan;
    if (h5io_mk_ds(gid, "QF5_GRAN_BADDETECTOR", H5T_STD_U8BE, 1, dim_siz2,
            &ds_id) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for QF5_GRAN_BADDETECTOR\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) uchar_data) != 0) {
        printf("%s, %d: Could not write to QF5_GRAN_BADDETECTOR\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close QF5_GRAN_BADDETECTOR\n",
                __FILE__, __LINE__);
        return 1;
    }
    free(uchar_data);
    /*
     *  Radiance This is a scan-by-scan output item
     */
    dim_siz[0] = out_rec->nlin;
    dim_siz[1] = out_rec->npix;
    /*
     *  There will have to be an initialization for the input radiance 
     *  dataset dependent on if we do a dummy or take from a L2. If we 
     *  follow the model for init_geo_data, it is put here.  However, input
     *  data for the TOA radiances either are dummy values or come from
     *  the L2 made by running l2gen in reverse.  So, the initialization of 
     *  input rads is best done in rd_l2_init under rd_sim_init. 
     */
    /*
     *  The radiance dataset is either unsigned short or float depending on
     *  the band number (argh), so set up accordingly
     */
    if (out_rec->out_bnd_typ[ibnd] == 0)
        iret = h5io_mk_ds(gid, "Radiance", H5T_STD_U16BE, 2, dim_siz,
            &(out_rec->bnd_dat_id[0][ibnd]));
    else
        iret = h5io_mk_ds(gid, "Radiance", H5T_IEEE_F32BE, 2, dim_siz,
            &(out_rec->bnd_dat_id[0][ibnd]));
    if (iret != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for Radiance\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  allocate the data transfer buffer and make it the same for output
     *  Note that dn storage allocation (if needed) is done above in init_sdr
     */
    if ((in_rec->bnd_lt[ibnd] = (float *)
            malloc(in_rec->npix * in_rec->ndet_scan * sizeof (float)))
            == NULL) {
        printf(
                "%s, %d: Could not allocate space for radiance data transfer buffer\n",
                __FILE__, __LINE__);
        return 1;
    }

    out_rec->bnd_lt[ibnd] = in_rec->bnd_lt[ibnd];
    /*
     *  set up the flags for the radiances - for missing or bad value output
     */
    if ((in_rec->bnd_q[ibnd] = (unsigned char *)
            malloc(in_rec->npix * in_rec->ndet_scan * sizeof (unsigned char)))
            == NULL) {
        printf(
                "%s, %d: Could not allocate space for output radiance flag buffer\n",
                __FILE__, __LINE__);
        return 1;
    }
    out_rec->bnd_q[ibnd] = in_rec->bnd_q[ibnd];
    /*
     *  RadianceFactors - 2 floats of scale and offset:
     *  only needed for scaled bands
     *  Lt = <scale> * count + <offset>
     *  short of finding min, max in granule, could do
     *  <offset> = 0, <scale> = 1.3 * f0(lambda) / ( PI * ( max_short - 10 ) )
     *  the f0 is in MKS already and the 1.3 is to allow some room
     *
     *  24 Feb 2011 WDR This and reflectance set-up does suprisingly well at 
     *  matching the IDPS scale... except for offset and in bands that saturate
     *  KEEP the old methods for setup and add the aray list set
     */
    dim_siz[0] = 2;
    /*
     *  The scaling is only set different from [1., 0.] for the scaled radiances
     */
    if (out_rec->out_bnd_typ[ibnd] == 0) {
        /*  KEEP, initial way of computing using F0 or top BBT of 500K
         *** ALSO note that if this is used again, the scale, offset needs to 
           be for radiance all the way through, not changing to BBT for M12 - 16
        if( out_rec->meas_typ[ibnd] == 0 )
          {
          out_rec->scale[ibnd] = ( 1.3 * out_rec->f0[ibnd] ) / 
            ( M_PI * ( float ) ( SOUB_UINT16_FILL - 1 ) );
          out_rec->offset[ibnd] = 0.;
          }
        else
          {
          lam_um = (float) out_rec->lam_band[ibnd] /1000.;
          out_rec->scale[ibnd] = bbt_2_rad( 500., lam_um ) /
            (float) ( SOUB_UINT16_FILL - 1 );
          out_rec->offset[ibnd] = 0.;
          }
         */
        /*  Use IDPS scale offset for radiance scaling */
        out_rec->scale[ibnd] = lcl_scale[ibnd];
        out_rec->offset[ibnd] = lcl_offset[ibnd];
    } else {
        out_rec->scale[ibnd] = 1.;
        out_rec->offset[ibnd] = 0.;
    }
    rad_fac[0] = out_rec->scale[ibnd];
    rad_fac[1] = out_rec->offset[ibnd];
    if (h5io_mk_ds(gid, "RadianceFactors", H5T_IEEE_F32BE, 1, dim_siz, &ds_id)
            != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for RadianceFactors\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) rad_fac) != 0) {
        printf("%s, %d: Could not write to RadianceFactors\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close RadianceFactors\n", __FILE__, __LINE__);
        return 1;
    }
    if (out_rec->meas_typ[ibnd] == 0) {
        /*
         *  Reflectance - just set up the dataset for output
         */
        dim_siz[0] = out_rec->nlin;
        dim_siz[1] = out_rec->npix;
        if (h5io_mk_ds(gid, "Reflectance", H5T_STD_U16BE, 2, dim_siz,
                &(out_rec->bnd_dat_id[1][ibnd])) != 0) {
            printf("%s, %d: Could not do h5io_mk_ds for Latitude\n",
                    __FILE__, __LINE__);
            return 1;
        }
        /*
         *  ReflectanceFactors - 2 floats of scale and offset for all bands
         *  ref = <scale> * count + <offset>
         *  short of finding min, max in granule, could do
         *  <offset> = 0, <scale> = 1.5 / ( max_short - 10 )
         */
        dim_siz[0] = 2;
        /*  AGAIN as in RadianceFactors above, KEEP this, same caveats
        out_rec->refl_scale[ibnd] = 1.5 / ( float) ( SOUB_UINT16_FILL - 1 );
        out_rec->refl_offset[ibnd] = 0.;
         */
        out_rec->refl_scale[ibnd] = lcl_scale_ref[ibnd];
        out_rec->refl_offset[ibnd] = lcl_offset_ref[ibnd];
        refl_fac[0] = out_rec->refl_scale[ibnd];
        refl_fac[1] = out_rec->refl_offset[ibnd];
        if (h5io_mk_ds(gid, "ReflectanceFactors", H5T_IEEE_F32BE, 1,
                dim_siz, &ds_id) != 0) {
            printf("%s, %d: Could not do h5io_mk_ds for ReflectanceFactors\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if (h5io_wr_ds(&ds_id, (void *) refl_fac) != 0) {
            printf("%s, %d: Could not write to ReflectanceFactors\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if (h5io_close(&ds_id) != 0) {
            printf("%s, %d: Could not close ReflectanceFactors\n",
                    __FILE__, __LINE__);
            return 1;
        }
    } else {
        /*
         *  BBT, This is a scan-by-scan output item
         */
        dim_siz[0] = out_rec->nlin;
        dim_siz[1] = out_rec->npix;
        /*
         *  The BBT dataset is either unsigned short or float depending on
         *  the band number, so set up accordingly
         */
        if (out_rec->out_bnd_typ[ibnd] == 0)
            iret = h5io_mk_ds(gid, "BrightnessTemperature", H5T_STD_U16BE, 2,
                dim_siz, &(out_rec->bnd_dat_id[1][ibnd]));
        else
            iret = h5io_mk_ds(gid, "BrightnessTemperature", H5T_IEEE_F32BE, 2,
                dim_siz, &(out_rec->bnd_dat_id[1][ibnd]));
        if (iret != 0) {
            printf("%s, %d: Could not do h5io_mk_ds for BrightnessTemperature\n",
                    __FILE__, __LINE__);
            return 1;
        }
        /*
         *  BrightnessTemperatureFactors - 2 floats of scale and offset:
         *  only needed for scaled bands
         *  BBT = <scale> * count + <offset>
         *  short of finding min, max in granule, could do
         *  <offset> = 0, <scale> = 500. / ( max_short - 10 )
         *
         *  use the refl factor storage for the BBT scaling
         */
        dim_siz[0] = 2;
        /*
         *  The scaling is only set different from [1., 0.] for the scaled radiances
         */
        if (out_rec->out_bnd_typ[ibnd] == 0) {
            /*  AGAIN as in RadianceFactors above, KEEP this, same caveats
            out_rec->refl_scale[ibnd] = 500. / ( ( float ) ( SOUB_UINT16_FILL - 1 ) );
            out_rec->refl_offset[ibnd] = 0.;
             */
            out_rec->refl_scale[ibnd] = lcl_scale_ref[ibnd];
            out_rec->refl_offset[ibnd] = lcl_offset_ref[ibnd];
        } else {
            out_rec->refl_scale[ibnd] = 1.;
            out_rec->refl_offset[ibnd] = 0.;
        }
        rad_fac[0] = out_rec->refl_scale[ibnd];
        rad_fac[1] = out_rec->refl_offset[ibnd];
        if (h5io_mk_ds(gid, "BrightnessTemperatureFactors", H5T_IEEE_F32BE,
                1, dim_siz, &ds_id) != 0) {
            printf(
                    "%s, %d: Could not do h5io_mk_ds for BrightnessTemperatureFactors\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if (h5io_wr_ds(&ds_id, (void *) rad_fac) != 0) {
            printf("%s, %d: Could not write to BrightnessTemperatureFactors\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if (h5io_close(&ds_id) != 0) {
            printf("%s, %d: Could not close BrightnessTemperatureFactors\n",
                    __FILE__, __LINE__);
            return 1;
        }
    }
    return 0;
}

int init_sdr_dpattr(int isdr, h5io_str *dat2_g_id, sdr_info_struc *sdr_info)
/*-----------------------------------------------------------------------------
    Routine:   init_sdr_dpattr

    Description:  make the data product attributes for the sdr files

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       isdr          I    SDR file #
        h5io_str *  dat2_g_id   I    id for group to place attributes in
        sdr_info_struc*  sdr_info  I  structure to hold geolocation
                                     information

    Modification history:

    W. Robinson, SAIC  21 Oct 2008  Original development

----------------------------------------------------------------------------*/ {
    int n_attr = 7, dims_1_1[] = {1, 1};
    int fsw_v = 0;
    char op_mode[] = "NPP Normal Operations, VIIRS Operational";
    int len_op_mode, len_domain;

    len_op_mode = strlen(op_mode);
    len_domain = strlen(sdr_info->domain);
    /* note that the attribute info is set up for the geo file and changes 
       are made for any band differences */
    h5attr_struc attrs[] = {
        { 1, 1, "Instrument_Short_name", H5T_NATIVE_CHAR, 6, 2, dims_1_1,
            (void *) "VIIRS"},
        { 1, 1, "N_Anc_Type_Tasked", H5T_NATIVE_CHAR, 9, 2, dims_1_1,
            (void *) "Official"},
        { 1, 1, "N_Collection_Short_Name", H5T_NATIVE_CHAR, 17, 2, dims_1_1,
            (void *) core_g_nm[isdr]},
        { 1, 1, "N_Dataset_Type_Tag", H5T_NATIVE_CHAR, 4, 2, dims_1_1,
            (void *) "GEO"},
        { 0, 0, "N_Instrument_Flight_SW_Version", H5T_STD_I32BE, 0, 2, dims_1_1,
            (void *) &fsw_v},
        { 1, 1, "N_Processing_Domain", H5T_NATIVE_CHAR, len_domain, 2, dims_1_1,
            (void *) sdr_info->domain},
        { 1, 1, "Operation_Mode", H5T_NATIVE_CHAR, len_op_mode, 2, dims_1_1,
            (void *) op_mode}
    };
    /*
     *  ATTRIB NOTES - hard coded except for domain and band specific for 
     *  N_Collection_Short_Name.  Some geo / band specific attrs
     */
    if (isdr > 0) {
        /* the N_Anc_Type_Tasked needs to be repressed in the SDR */
        attrs[1].express = 0;
        /* the N_Dataset_Type_Tag is SDR for the band files */
        attrs[3].data = (void *) "SDR";
        /* the N_Instrument_Flight_SW_Version needs to be expressed in the SDR */
        attrs[4].express = 1;
    }
    /*
     *  output the attributes to the location
     */
    if (wr_attr_seq(dat2_g_id, n_attr, attrs) != 0) {
        printf("%s, %d: Could not write data prod attributes\n",
                __FILE__, __LINE__);
        return 1;
    }
    return 0;
}

int init_sdr_agg(int isdr, h5io_str *gid, sdr_info_struc *sdr_info)
/*-----------------------------------------------------------------------------
    Routine:   init_sdr_agg

    Description:  create the aggregation dataset and attributes 
      in the data products group for a SDR file

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       isdr          I    SDR file number
        h5io_str *  gid         I    id for group to place dataset and 
                                     attributes 
        sdr_info_struc *  sdr_info  I  structure to hold geolocation
                                     information

    Modification history:

    W. Robinson, SAIC  21 Oct 2008  Original development

----------------------------------------------------------------------------*/ {
    int n_attr = 9, dims_1_1[] = {1, 1};
    char grp_name[100];
    uint64 abon = 14;
    int ang = 1;
    h5attr_struc attrs[] = {
        { 1, 1, "AggregateBeginningDate", H5T_NATIVE_CHAR, 9, 2, dims_1_1,
            (void *) sdr_info->st_date},
        { 1, 1, "AggregateBeginningGranuleID", H5T_NATIVE_CHAR, 16, 2, dims_1_1,
            (void *) "NPP001212012151"},
        { 1, 0, "AggregateBeginningOrbitNumber", H5T_STD_U64BE, 0, 2, dims_1_1,
            (void *) &abon},
        { 1, 1, "AggregateBeginningTime", H5T_NATIVE_CHAR, 15, 2, dims_1_1,
            (void *) sdr_info->st_time},
        { 1, 1, "AggregateEndingDate", H5T_NATIVE_CHAR, 9, 2, dims_1_1,
            (void *) sdr_info->en_date},
        { 1, 1, "AggregateEndingGranuleID", H5T_NATIVE_CHAR, 16, 2, dims_1_1,
            (void *) "NPP001212012151"},
        { 1, 0, "AggregateEndingOrbitNumber", H5T_STD_U64BE, 0, 2, dims_1_1,
            (void *) &abon},
        { 1, 1, "AggregateEndingTime", H5T_NATIVE_CHAR, 15, 2, dims_1_1,
            (void *) sdr_info->en_time},
        { 1, 0, "AggregateNumberGranules", H5T_STD_I32BE, 0, 2, dims_1_1,
            (void *) &ang}
    };
    /*
     *  ATTRIB NOTES - some values are hard-coded now
     */
    h5io_str ds_id;
    int dim_siz[2], *arr_int;
    /*
     *  the dataset is a reference dataset, so we'll forgo and just place an int
     */
    arr_int = (int *) malloc(sizeof ( int));
    *arr_int = 666;
    dim_siz[0] = 1;
    sprintf(grp_name, "%s_Aggr", core_g_nm[isdr]);
    if (h5io_mk_ds(gid, grp_name, H5T_NATIVE_INT, 1,
            dim_siz, &ds_id) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for grp VIIRS-MOD-GEO-TC_Aggr\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) arr_int) != 0) {
        printf("%s, %d: Could not write to VIIRS-MOD-GEO-TC_Aggr\n",
                __FILE__, __LINE__);
        return 1;
    }
    free(arr_int);
    /*
     *  output the attributes to the location
     */
    if (wr_attr_seq(&ds_id, n_attr, attrs) != 0) {
        printf(
                "%s, %d: Could not write top attributes for VIIRS-MOD-GEO-TC_Aggr\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s, %d: Could not close VIIRS-MOD-GEO-TC_Aggr\n",
                __FILE__, __LINE__);
        return 1;
    }
    return 0;
}

int init_sdr_gran(int isdr, h5io_str *gid, sdr_info_struc *sdr_info, out_rec_struc *out_rec)
/*-----------------------------------------------------------------------------
    Routine:   init_geo_gran

    Description:  make the granule dataset and attributes under the data 
      products groups

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       isdr          I    SDR file number
        h5io_str *  dat2_g_id   I    id for group to place attributes in
        sdr_info_struc *  sdr_info  I  structure to hold geolocation
                                     information
        out_rec_struc * out_rec I    output record information

    Modification history:

    W. Robinson, SAIC  21 Oct 2008  Original development
    W. Robinson, SAIC  13 Oct 2011  change N_Granule_ID to be VYYYYMMDDHHMMSS

----------------------------------------------------------------------------*/ {
    int n_attr = 49, dims_1_1[] = {1, 1}, len_str, i;
    int dims_8_1[] = {8, 1}, dims_9_1[] = {9, 1};
    int dims_2_1[] = {2, 1}, dims_anc[] = {19, 1}, dims_aux[] = {5, 1};
    int qual_sum_val[2];
    float g_ring_lat[8], g_ring_lon[8], gran_bdy[4], nadir_bdy[4];
    float fill = 0.;
    uint64 beg_orb = 14;
    char char33_9[33 * 9], band_id[4], char_100_49[100 * 49 ], grp_nam[100];
    char gran_vers[10] = "A1M", doc_ref[10] = "N/A";
    char *sw_vers = "viirs_sim_sdr_v0.0", gran_id[16];
    char qual_sum_name[ 26 * 2 ], char_anc[ 19 * 106 ], char_aux[ 5 * 106];
    unsigned char asc_dec = 0;

    len_str = strlen(sw_vers);
    sprintf(gran_id, "V%8.8s%6.6s", sdr_info->st_date, sdr_info->st_time);
    h5attr_struc attrs[] = {
        { 1, 0, "Ascending/Descending_Indicator", H5T_NATIVE_UCHAR, 0, 2, dims_1_1,
            (void *) &asc_dec},
        { 0, 1, "Band_ID", H5T_NATIVE_CHAR, 2, 2, dims_1_1,
            (void *) band_id},
        { 1, 1, "Beginning_Date", H5T_NATIVE_CHAR, 9, 2, dims_1_1,
            (void *) sdr_info->st_date},
        { 1, 1, "Beginning_Time", H5T_NATIVE_CHAR, 15, 2, dims_1_1,
            (void *) sdr_info->st_time},
        { 0, 0, "East_Bounding_Coordinate", H5T_IEEE_F32BE, 0, 2, dims_1_1,
            (void *) gran_bdy},
        { 1, 1, "Ending_Date", H5T_NATIVE_CHAR, 9, 2, dims_1_1,
            (void *) sdr_info->en_date},
        { 1, 1, "Ending_Time", H5T_NATIVE_CHAR, 15, 2, dims_1_1,
            (void *) sdr_info->en_time},
        { 1, 0, "G-Ring-Latitude", H5T_IEEE_F32BE, 0, 2, dims_8_1,
            (void *) g_ring_lat},
        { 1, 0, "G-Ring-Longitude", H5T_IEEE_F32BE, 0, 2, dims_8_1,
            (void *) g_ring_lon},
        { 1, 1, "N_Algorithm_Version", H5T_NATIVE_CHAR, 8, 2, dims_1_1,
            (void *) "Fill"},
        { 1, 1, "N_Anc_Filename", H5T_NATIVE_CHAR, 106, 2, dims_anc,
            (void *) char_anc},
        { 1, 1, "N_Aux_Filename", H5T_NATIVE_CHAR, 106, 2, dims_aux,
            (void *) char_aux},
        { 1, 0, "N_Beginning_Orbit_Number", H5T_STD_U64BE, 0, 2, dims_1_1,
            (void *) &beg_orb},
        { 1, 0, "N_Beginning_Time_IET", H5T_STD_U64BE, 0, 1, dims_1_1,
            (void *) &(sdr_info->st_58_t)},
        { 1, 1, "N_Creation_Date", H5T_NATIVE_CHAR, 9, 2, dims_1_1,
            (void *) sdr_info->cre_date},
        { 1, 1, "N_Creation_Time", H5T_NATIVE_CHAR, 15, 2, dims_1_1,
            (void *) sdr_info->cre_time},
        { 1, 0, "N_Ending_Time_IET", H5T_STD_U64BE, 0, 1, dims_1_1,
            (void *) &(sdr_info->en_58_t)},
        { 0, 1, "N_Graceful_Degradation", H5T_NATIVE_CHAR, 3, 2, dims_1_1,
            (void *) "No"},
        { 1, 1, "N_Granule_ID", H5T_NATIVE_CHAR, 15, 2, dims_1_1,
            (void *) gran_id},
        /*  old val inserted: (void *)"NPP001212012151" }, */
        { 1, 1, "N_Granule_Status", H5T_NATIVE_CHAR, 5, 2, dims_1_1,
            (void *) "Fill"},
        { 1, 1, "N_Granule_Version", H5T_NATIVE_CHAR, 4, 2, dims_1_1,
            (void *) gran_vers},
        { 1, 1, "N_Input_Prod", H5T_NATIVE_CHAR, 33, 2, dims_9_1,
            (void *) char33_9},
        { 1, 1, "N_LEOA_Flag", H5T_NATIVE_CHAR, 4, 2, dims_1_1,
            (void *) "Off"},
        { 1, 0, "N_Nadir_Latitude_Max", H5T_IEEE_F32BE, 0, 2, dims_1_1,
            (void *) (nadir_bdy + 2)},
        { 1, 0, "N_Nadir_Latitude_Min", H5T_IEEE_F32BE, 0, 2, dims_1_1,
            (void *) (nadir_bdy + 3)},
        { 1, 0, "N_Nadir_Longitude_Max", H5T_IEEE_F32BE, 0, 2, dims_1_1,
            (void *) nadir_bdy},
        { 1, 0, "N_Nadir_Longitude_Min", H5T_IEEE_F32BE, 0, 2, dims_1_1,
            (void *) (nadir_bdy + 1)},
        { 1, 1, "N_NPOES_Document_Ref", H5T_NATIVE_CHAR, 4, 2, dims_1_1,
            (void *) doc_ref},
        { 1, 0, "N_Number_Of_Scans", H5T_STD_I32BE, 0, 2, dims_1_1,
            (void *) &(out_rec->nscan)},
        { 0, 0, "N_Percent_Erroneous_Data", H5T_IEEE_F32BE, 0, 2, dims_1_1,
            (void *) &fill},
        { 0, 0, "N_Percent_Missing_Data", H5T_IEEE_F32BE, 0, 2, dims_1_1,
            (void *) &fill},
        { 0, 0, "N_Percent_Not_Applicable_Data", H5T_IEEE_F32BE, 0, 2, dims_1_1,
            (void *) &fill},
        { 1, 1, "N_Quality_Summary_Names", H5T_NATIVE_CHAR, 26, 2, dims_2_1,
            (void *) qual_sum_name},
        { 1, 0, "N_Quality_Summary_Values", H5T_STD_I32BE, 0, 2, dims_2_1,
            (void *) qual_sum_val},
        { 1, 1, "N_Reference_ID", H5T_NATIVE_CHAR, 10, 2, dims_1_1,
            (void *) "Who_Cares"},
        { 0, 0, "N_Satellite/Local_Azimuth_Angle_Max", H5T_IEEE_F32BE, 0, 2,
            dims_1_1, (void *) &fill},
        { 0, 0, "N_Satellite/Local_Azimuth_Angle_Min", H5T_IEEE_F32BE, 0, 2,
            dims_1_1, (void *) &fill},
        { 0, 0, "N_Satellite/Local_Zenith_Angle_Max", H5T_IEEE_F32BE, 0, 2,
            dims_1_1, (void *) &fill},
        { 0, 0, "N_Satellite/Local_Zenith_Angle_Min", H5T_IEEE_F32BE, 0, 2,
            dims_1_1, (void *) &fill},
        { 1, 1, "N_Software_Version", H5T_NATIVE_CHAR, len_str, 2, dims_1_1,
            (void *) sw_vers},
        { 0, 0, "N_Solar_Azimuth_Angle_Max", H5T_IEEE_F32BE, 0, 2,
            dims_1_1, (void *) &fill},
        { 0, 0, "N_Solar_Azimuth_Angle_Min", H5T_IEEE_F32BE, 0, 2,
            dims_1_1, (void *) &fill},
        { 0, 0, "N_Solar_Zenith_Angle_Max", H5T_IEEE_F32BE, 0, 2,
            dims_1_1, (void *) &fill},
        { 0, 0, "N_Solar_Zenith_Angle_Min", H5T_IEEE_F32BE, 0, 2,
            dims_1_1, (void *) &fill},
        { 1, 1, "N_Spacecraft_Maneuver", H5T_NATIVE_CHAR, 18, 2, dims_1_1,
            (void *) "Normal Operations"},
        { 0, 0, "North_Bounding_Coordinate", H5T_IEEE_F32BE, 0, 2, dims_1_1,
            (void *) (gran_bdy + 2)},
        { 0, 0, "South_Bounding_Coordinate", H5T_IEEE_F32BE, 0, 2, dims_1_1,
            (void *) (gran_bdy + 3)},
        { 0, 0, "West_Bounding_Coordinate", H5T_IEEE_F32BE, 0, 2, dims_1_1,
            (void *) (gran_bdy + 1)},
        { 1, 1, "N_Day_Night_Flag", H5T_NATIVE_CHAR, 4, 2, dims_1_1,
            (void *) "Day"}
    };
    /*
     *  ATTRIB NOTES - many attrs are hard-coded now
     */
    h5io_str ds_id;
    int dim_siz[2], *arr_int;
    /*
     *  set the g_ring values here - coordinates of granule perimeter)
     */
    for (i = 0; i < 8; i++) {
        *(g_ring_lat + i) = (float) i * 2.5;
        *(g_ring_lon + i) = (float) i * 25;
    }
    /*
     *  set something in the input prod strings demo
     */
    for (i = 0; i < 9; i++)
        sprintf((char33_9 + (i * 33)), "Input_prod_#_%d", i);
    /*
     *  Modifications for the band files
     */
    if (isdr > 0) {
        for (i = 0; i < n_attr; i++)
            attrs[i].express = 1;
        sprintf(band_id, "M%d", isdr);
    }
    for (i = 0; i < 49; i++)
        sprintf((char_100_49 + (i * 100)), "Aux file #%d", i);
    for (i = 0; i < 4; i++) {
        *(gran_bdy + i) = 0.;
        *(nadir_bdy + i) = 0.;
    }
    strcpy(gran_vers, "A1");
    strcpy(doc_ref, "TBD");
    /*
     *  quality summary names fill
     */
    for (i = 0; i < 2; i++) {
        sprintf((qual_sum_name + (i * 26)), "qual sum name %d", i);
        *(qual_sum_val + i) = i;
    }
    /*
     * anc and aux fill
     */
    for (i = 0; i < 19; i++)
        sprintf((char_anc + (i * 106)), "Anc file # %d", i);
    for (i = 0; i < 5; i++)
        sprintf((char_aux + (i * 106)), "Aux file # %d", i);
    /*
     *  the dataset is a reference dataset, so we'll forgo and just place an int
     */
    arr_int = (int *) malloc(sizeof ( int));
    *arr_int = 666;
    dim_siz[0] = 1;
    sprintf(grp_nam, "%s_Gran_0", core_g_nm[isdr]);
    if (h5io_mk_ds(gid, grp_nam, H5T_NATIVE_INT, 1,
            dim_siz, &ds_id) != 0) {
        printf("%s, %d: Could not do h5io_mk_ds for VIIRS-MOD-GEO-TC_Gran_0\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_wr_ds(&ds_id, (void *) arr_int) != 0) {
        printf("%s, %d: Could not write to VIIRS-MOD-GEO-TC_Gran_0\n",
                __FILE__, __LINE__);
        return 1;
    }
    free(arr_int);
    /*
     *  output the attributes to the location
     */
    if (wr_attr_seq(&ds_id, n_attr, attrs) != 0) {
        printf(
                "%s %d: Could not write top attributes for VIIRS-MOD-GEO-TC_Gran_0\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s %d: Could not close VIIRS-MOD-GEO-TC_Gran_0\n",
                __FILE__, __LINE__);
        return 1;
    }
    return 0;
}
