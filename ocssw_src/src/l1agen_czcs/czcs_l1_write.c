#include <stdio.h>
#include <libgen.h>
#include "hdf.h"
#include "mfhdf.h"

#include "l1czcs.h"       /* local include file */

int czcs_l1_write(char *ofile, l1_data_struc l1_data, gattr_struc gattr)
/*******************************************************************

   czcs_l1_write

   purpose: Write the CZCS L1A dataset

   Returns type: int 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            ofile            I      output file name
      l1_data_struc     l1_data          I      data counts, lat, lon
      gattr_struc       gattr            I      file attributes

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson        9-Mar-2004     adapted from seadas version
                                        by removing global attributes: 
                                          Area Code, Circle Parameters 
                                        adding global attributes: 
                                          Slope, Intercept, Roll, Pitch, Yaw

 *******************************************************************/ {
    int32 sdfid, fid, raw_vid, nav_vid, sl_attr_vid;
    int iret, tot_line, tot_pixel;
    char *file;

    /*
     *  get count array size from structure
     */
    tot_line = gattr.scan_lines;
    tot_pixel = gattr.pix_per_scan;

    /*
     *  open output file
     */
    if ((sdfid = SDstart(ofile, DFACC_CREATE)) < 0) {
        fprintf(stderr, "[czcs_l1_write] Error: SDstart failed on %s\n", ofile);
        return FAIL;
    }

    if ((fid = Hopen(ofile, DFACC_RDWR, 0)) < 0) {
        fprintf(stderr, "[czcs_l1_write] Error: Hopen failed on %s\n", ofile);
        return FAIL;
    }

    /*  create global attributes  */

    file = basename(ofile);
    iret = create_global_attribute(file, sdfid, gattr);

    Vstart(fid);

    /* create all the vgroups */

    sl_attr_vid = Vattach(fid, -1, "w");
    Vsetname(sl_attr_vid, "Scan-Line Attributes");

    raw_vid = Vattach(fid, -1, "w");
    Vsetname(raw_vid, "Raw CZCS Data");

    nav_vid = Vattach(fid, -1, "w");
    Vsetname(nav_vid, "Navigation");

    /* create all bands SDSs */

    iret = FAIL;
    if (create_band_sds(sdfid, raw_vid, l1_data.counts, tot_pixel,
            tot_line) >= 0) {
        /* add the quality per line data */
        if (wrt_czcs_qual(sdfid, raw_vid, tot_line, l1_data)
                == SUCCEED) {
            /* create navigation SDS */
            if (set_czcs_ctl_data(sdfid, nav_vid, gattr, l1_data)
                    >= 0) {
                if (wrt_czcs_sla(sdfid, sl_attr_vid, tot_line, l1_data)
                        == SUCCEED)
                    iret = SUCCEED;
            }
        }
    }
    Vdetach(raw_vid);
    Vdetach(nav_vid);
    Vdetach(sl_attr_vid);
    Vend(fid);
    SDend(sdfid);
    Hclose(fid);

    return iret;
}

int create_global_attribute(char *file, int sdfid, gattr_struc gattr) {
#define  TITLE  "Title"
#define  TITLE_VAL "CZCS Level-1A Data"
#define  SENSOR_VAL     "Coastal Zone Color Scanner (CZCS)"
#define  MISSION        "Mission"
#define  MISSION_VAL    "Nimbus CZCS"
#define  MISSIONCHAR    "Mission Characteristics"
#define  MISSIONCHAR_VAL "Nominal orbit: inclination = 99.3 (Sun-Synchronous); node = 11:52 AM local (ascending); eccentricity =< 0.0009; altitude = 955km; ground speed = 6.4km/sec"
#define  SENSOR         "Sensor"
#define  SENSOR_VAL     "Coastal Zone Color Scanner (CZCS)"
#define  SENSORCHAR     "Sensor Characteristics"
#define  SENSORCHAR_VAL "Number of bands = 6; number of active bands = 6; wavelengths per band (nm) = 443, 520, 550, 670, 750, 11500; bits per pixel = 8; instantaneous field-of-view = .865 mrad; pixels per scan = 1968; scan rate = 8.08/sec"
#define  REPL_FLG       "ORIGINAL"

    int cl;

    /***  write Product Name and "Title"  */
    if ((SDsetattr(sdfid, "Product Name", DFNT_CHAR, strlen(file) + 1,
            (VOIDP) file)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, TITLE, DFNT_CHAR, strlen(TITLE_VAL) + 1,
            (VOIDP) TITLE_VAL)) < 0)
        return FAIL;

    /***  write "Data type"  */
    if ((SDsetattr(sdfid, "Data Type", DFNT_CHAR, strlen(gattr.datatype) + 1,
            (VOIDP) gattr.datatype)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Replacement Flag", DFNT_CHAR, strlen(REPL_FLG) + 1,
            (VOIDP) REPL_FLG)) < 0)
        return FAIL;

    /***  write "Data Center" */
    if ((SDsetattr(sdfid, "Data Center", DFNT_CHAR, strlen(gattr.datacenter) + 1,
            (VOIDP) gattr.datacenter)) < 0)
        return FAIL;

    /***  write "Station Name" */
    if ((SDsetattr(sdfid, "Station Name", DFNT_CHAR, strlen(gattr.stn_name) + 1,
            (VOIDP) gattr.stn_name)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Station Latitude", DFNT_FLOAT32, 1,
            (VOIDP) & gattr.stn_lat)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Station Longitude", DFNT_FLOAT32, 1,
            (VOIDP) & gattr.stn_lon)) < 0)
        return FAIL;

    /***  write "Mission" */
    if ((SDsetattr(sdfid, MISSION, DFNT_CHAR, strlen(MISSION_VAL) + 1,
            (VOIDP) MISSION_VAL)) < 0)
        return FAIL;

    /***  write "Mission Characteristics" */
    if ((SDsetattr(sdfid, MISSIONCHAR, DFNT_CHAR, strlen(MISSIONCHAR_VAL) + 1,
            (VOIDP) MISSIONCHAR_VAL)) < 0)
        return FAIL;

    /***  write "Sensor" */
    if ((SDsetattr(sdfid, SENSOR, DFNT_CHAR, strlen(SENSOR_VAL) + 1,
            (VOIDP) SENSOR_VAL)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, SENSORCHAR, DFNT_CHAR, strlen(SENSORCHAR_VAL) + 1,
            (VOIDP) SENSORCHAR_VAL)) < 0)
        return FAIL;

    /***  write "Software ID" */
    if ((SDsetattr(sdfid, "Software ID", DFNT_CHAR,
            strlen(gattr.soft_id) + 1, (VOIDP) gattr.soft_id)) < 0)
        return FAIL;

    /***  write "Processing Time" */
    if ((SDsetattr(sdfid, "Processing Time", DFNT_CHAR,
            strlen(gattr.process_time) + 1, (VOIDP) gattr.process_time)) < 0)
        return FAIL;

    /***  write "Input Files" */
    if ((SDsetattr(sdfid, "Input Files", DFNT_CHAR,
            strlen(gattr.input_files) + 1, (VOIDP) gattr.input_files)) < 0)
        return FAIL;

    /***  write "Processing Control" */
    if ((SDsetattr(sdfid, "Processing Control", DFNT_CHAR,
            strlen(gattr.proc_ctl) + 1, (VOIDP) gattr.proc_ctl)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Pixels per Scan Line", DFNT_INT32, 1,
            (VOIDP) & gattr.pix_per_scan)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Number of Scan Lines", DFNT_INT32, 1,
            &gattr.scan_lines)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Number of Pixel Control Points", DFNT_INT32, 1,
            &gattr.n_ctl_pt)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Number of Scan Control Points", DFNT_INT32, 1,
            &gattr.scan_lines)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "LAC Pixel Start Number", DFNT_INT32, 1,
            &gattr.lac_pixl_start_no)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "LAC Pixel Subsampling", DFNT_INT32, 1,
            &gattr.lac_pixl_subsample)) < 0)
        return FAIL;

    cl = gattr.scan_lines / 2 + 1;
    if ((SDsetattr(sdfid, "Scene Center Scan Line", DFNT_INT32, 1, &cl)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Filled Scan Lines", DFNT_INT32, 1,
            &gattr.scan_lines)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Start Time", DFNT_CHAR, strlen(gattr.start_time) + 1,
            &gattr.start_time)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "End Time", DFNT_CHAR, strlen(gattr.end_time) + 1,
            &gattr.end_time)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Scene Center Time", DFNT_CHAR,
            strlen(gattr.center_time) + 1, &gattr.center_time)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Start Year", DFNT_INT16, 1,
            &gattr.start_year)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Start Day", DFNT_INT16, 1,
            &gattr.start_day)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Start Millisec", DFNT_INT32, 1,
            &gattr.start_msec)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "End Year", DFNT_INT16, 1,
            &gattr.end_year)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "End Day", DFNT_INT16, 1,
            &gattr.end_day)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "End Millisec", DFNT_INT32, 1,
            &gattr.end_msec)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Start Node", DFNT_CHAR,
            strlen("Ascending") + 1, "Ascending")) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "End Node", DFNT_CHAR,
            strlen("Ascending") + 1, "Ascending")) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Orbit Number", DFNT_INT32, 1,
            &gattr.orbit)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Sensor Tilt", DFNT_FLOAT32, 1,
            &gattr.tilt)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Scene Gain", DFNT_INT32, 1,
            &gattr.gain)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Thresh", DFNT_INT32, 1,
            &gattr.thresh)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Latitude Units", DFNT_CHAR,
            strlen("degrees North") + 1, "degrees North")) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Longitude Units", DFNT_CHAR,
            strlen("degrees East") + 1, "degrees East")) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Scene Center Latitude",
            DFNT_FLOAT32, 1, (VOIDP) & gattr.center_lat)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Scene Center Longitude",
            DFNT_FLOAT32, 1, (VOIDP) & gattr.center_lon)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Scene Center Solar Zenith",
            DFNT_FLOAT32, 1, (VOIDP) & gattr.cntr_sol_zen)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Upper Left Latitude",
            DFNT_FLOAT32, 1, (VOIDP) & gattr.up_lft_lat)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Upper Left Longitude",
            DFNT_FLOAT32, 1, (VOIDP) & gattr.up_lft_lon)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Upper Right Latitude",
            DFNT_FLOAT32, 1, (VOIDP) & gattr.up_rgt_lat)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Upper Right Longitude",
            DFNT_FLOAT32, 1, (VOIDP) & gattr.up_rgt_lon)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Lower Left Latitude",
            DFNT_FLOAT32, 1, (VOIDP) & gattr.lo_lft_lat)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Lower Left Longitude",
            DFNT_FLOAT32, 1, (VOIDP) & gattr.lo_lft_lon)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Lower Right Latitude",
            DFNT_FLOAT32, 1, (VOIDP) & gattr.lo_rgt_lat)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Lower Right Longitude",
            DFNT_FLOAT32, 1, (VOIDP) & gattr.lo_rgt_lon)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Start Center Latitude",
            DFNT_FLOAT32, 1, (VOIDP) & gattr.start_cntr_lat)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Start Center Longitude",
            DFNT_FLOAT32, 1, (VOIDP) & gattr.start_cntr_lon)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "End Center Latitude",
            DFNT_FLOAT32, 1, (VOIDP) & gattr.end_cntr_lat)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "End Center Longitude",
            DFNT_FLOAT32, 1, (VOIDP) & gattr.end_cntr_lon)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Northernmost Latitude", DFNT_FLOAT32, 1,
            (VOIDP) & gattr.limits[0])) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Southernmost Latitude", DFNT_FLOAT32, 1,
            (VOIDP) & gattr.limits[1])) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Westernmost Longitude", DFNT_FLOAT32, 1,
            (VOIDP) & gattr.limits[2])) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Easternmost Longitude", DFNT_FLOAT32, 1,
            (VOIDP) & gattr.limits[3])) < 0)
        return FAIL;

    /*
     *  also output the slope, intercept, roll, pitch and yaw
     */
    if ((SDsetattr(sdfid, "Calibration Slope", DFNT_FLOAT32, 6,
            (VOIDP) & gattr.slope)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Calibration Intercept", DFNT_FLOAT32, 6,
            (VOIDP) & gattr.intercept)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Center Roll", DFNT_FLOAT32, 1, (VOIDP) & gattr.roll)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Center Pitch", DFNT_FLOAT32, 1,
            (VOIDP) & gattr.pitch)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Center Yaw", DFNT_FLOAT32, 1, (VOIDP) & gattr.yaw)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "ILT Flags", DFNT_UINT8, 1,
            (VOIDP) & gattr.ilt_flags)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Parameter Presence Code", DFNT_UINT8, 1,
            (VOIDP) & gattr.parm_presence)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Number of Missing Scan Lines", DFNT_INT16, 1,
            (VOIDP) & gattr.n_miss_scans)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Number of Scans with Missing Channels",
            DFNT_INT16, 6, (VOIDP) & gattr.n_scan_mis_chan)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Number of HDT Sync Losses", DFNT_INT16, 1,
            (VOIDP) & gattr.n_hdt_sync_loss)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Number of HDT Parity Errors", DFNT_INT16, 1,
            (VOIDP) & gattr.n_hdt_parity_err)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Number of WBVT Sync Losses", DFNT_INT16, 1,
            (VOIDP) & gattr.n_wbvt_sync_loss)) < 0)
        return FAIL;

    if ((SDsetattr(sdfid, "Number of WBVT Slip Occurrences", DFNT_INT16, 1,
            (VOIDP) & gattr.n_wbvt_slips)) < 0)
        return FAIL;

    return SUCCEED;
}

int create_band_sds(int sdfid, int raw_vid, unsigned char *counts[],
        int tot_pixel, int tot_line) {
    int32 create_sds(int32, char *, int32, int32, int32 *, int32, VOIDP *);
    int32 set_dim_names(int32, char *, char *, char *);
    int dims[3] = {0, 0, 0}, sdsid;
    char *band_name[6] = {"band1", "band2", "band3", "band4", "band5", "band6"};
    unsigned char *tmp_ptr;
    char sds_attr[128];
    int i, iret;

    dims[0] = tot_line;
    dims[1] = tot_pixel;

    for (i = 0; i < 6; i++) {
        tmp_ptr = (unsigned char *) (counts[i]);
        if ((sdsid = create_sds(sdfid, band_name[i], DFNT_UINT8, 2,
                (int32 *) dims, raw_vid, (VOIDP) tmp_ptr)) < 0) {
            fprintf(stderr, "create_band_sds: create_sds failed on %s\n",
                    band_name[i]);
            return FAIL;
        }


        if ((iret =
                set_dim_names(sdsid, "Number of Scan Lines", "Pixels per Scan Line", NULL))
                < 0) {
            fprintf(stderr, "create_band_sds: set_dim_names failed on %s\n", band_name[i]);
            return FAIL;
        }

        sprintf(sds_attr, "Level-1A %s data", band_name[i]);
        if ((iret = SDsetattr(sdsid, "long_name", DFNT_CHAR, strlen(sds_attr) + 1,
                (VOIDP) sds_attr)) < 0) {
            fprintf(stderr, "create_band_sds: SDsetattr failed on %s\n", band_name[i]);
            return FAIL;
        };
        strcpy(sds_attr, "radiance counts");
        if ((iret = SDsetattr(sdsid, "units", DFNT_CHAR, strlen(sds_attr) + 1,
                (VOIDP) sds_attr)) < 0) {
            fprintf(stderr, "create_band_sds: SDsetattr failed on %s\n", band_name[i]);
            return FAIL;
        };

    }

    return SUCCEED;
}
