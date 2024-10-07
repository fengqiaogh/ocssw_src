#include "l1czcs.h"
#include "hdfio.h"

int cz_l1_read(char *file, int r_mode, gattr_struc *cz_attr,
        l1_data_struc *cz_dat)
/*******************************************************************

   cz_l1_read

   purpose: read the attributes and data from a CZCS L1 file

   Returns type: int - return status: 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      name of file to open
      int               r_mode           I      read mode: 0 - read all the 
                                                data into the structures
                                                1 - read time and quality data
                                                only
      gattr_struc *     cz_attr          O      structure of global attributes
      l1_data_struc *   cz_dat           O      structure of SDS data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 6 Aug 2004      Original development
      W. Robinson, SAIC 20 Dec 2005     add the pos_err or position error 
                                        sds to the reading

 *******************************************************************/
 {
    hdfio_struc hdfid;
    int nscan, n_ctl_pt;
    unsigned char *tmp_ptr;
    /*
     *  open the file
     */
    if (hdfio_open(file, &hdfid) != 0) return -1;
    /*
     *  read the basic attributes first
     *  # scans
     */
    if (hdfio_rd_gattr(hdfid, "Number of Scan Lines", DFNT_INT32,
            1, (void *) &cz_attr->scan_lines) != 0) return -1;
    nscan = cz_attr->scan_lines;

    /* parm presence */
    if (hdfio_rd_gattr(hdfid, "Parameter Presence Code", DFNT_UINT8, 1,
            (void *) &cz_attr->parm_presence) != 0) return -1;

    /* # pixel ctl points  */
    if (hdfio_rd_gattr(hdfid, "Number of Pixel Control Points", DFNT_INT32,
            1, (void *) &cz_attr->n_ctl_pt) != 0) return -1;
    n_ctl_pt = cz_attr->n_ctl_pt;

    /*  start day */
    if (hdfio_rd_gattr(hdfid, "Start Day", DFNT_INT16,
            1, (void *) &cz_attr->start_day) != 0) return -1;

    /*
     *  set up the data storage 
     */
    cz_dat_alloc(nscan, n_ctl_pt, r_mode, cz_dat);
    /*
     *  read in the arrays
     */
    if (hdfio_rd_sd(hdfid, "msec", (void *) cz_dat->msec) != 0) return -1;
    if (hdfio_rd_sd(hdfid, "cal_sum", (void *) cz_dat->cal_sum) != 0)
        return -1;
    if (hdfio_rd_sd(hdfid, "cal_scan", (void *) cz_dat->cal_scan) != 0)
        return -1;
    /*
     *  if mode is set for reading all the data, read the rest
     */
    if (r_mode == 0) {
        /*
         *  first, read the rest of the attributes
         */
        if (hdfio_rd_gattr(hdfid, "Product Name", DFNT_CHAR,
                40, (void *) cz_attr->prod_name) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Title", DFNT_CHAR,
                40, (void *) cz_attr->f_title) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Data Center", DFNT_CHAR,
                256, (void *) cz_attr->datacenter) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Station Name", DFNT_CHAR,
                256, (void *) cz_attr->stn_name) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Station Latitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->stn_lat) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Station Longitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->stn_lon) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Mission", DFNT_CHAR,
                40, (void *) cz_attr->mission) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Mission Characteristics", DFNT_CHAR,
                256, (void *) cz_attr->mission_char) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Sensor", DFNT_CHAR,
                40, (void *) cz_attr->sensor) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Sensor Characteristics", DFNT_CHAR,
                256, (void *) cz_attr->sensor_char) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Data Type", DFNT_CHAR,
                16, (void *) cz_attr->datatype) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Replacement Flag", DFNT_CHAR,
                40, (void *) cz_attr->repl_flg) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Software ID", DFNT_CHAR,
                128, (void *) cz_attr->soft_id) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Processing Time", DFNT_CHAR,
                17, (void *) cz_attr->process_time) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Input Files", DFNT_CHAR,
                1280, (void *) cz_attr->input_files) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Processing Control", DFNT_CHAR,
                1024, (void *) cz_attr->proc_ctl) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Start Time", DFNT_CHAR,
                17, (void *) cz_attr->start_time) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "End Time", DFNT_CHAR,
                17, (void *) cz_attr->end_time) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Scene Center Time", DFNT_CHAR,
                17, (void *) cz_attr->center_time) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Start Year", DFNT_INT16,
                1, (void *) &cz_attr->start_year) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Start Millisec", DFNT_INT32,
                1, (void *) &cz_attr->start_msec) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "End Year", DFNT_INT16,
                1, (void *) &cz_attr->end_year) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "End Day", DFNT_INT16,
                1, (void *) &cz_attr->end_day) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "End Millisec", DFNT_INT32,
                1, (void *) &cz_attr->end_msec) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Start Node", DFNT_CHAR,
                11, (void *) cz_attr->start_node) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "End Node", DFNT_CHAR,
                11, (void *) cz_attr->end_node) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Orbit Number", DFNT_INT32,
                1, (void *) &cz_attr->orbit) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Pixels per Scan Line", DFNT_INT32,
                1, (void *) &cz_attr->pix_per_scan) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Number of Scan Lines", DFNT_INT32,
                1, (void *) &cz_attr->scan_lines) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Number of Scan Control Points", DFNT_INT32,
                1, (void *) &cz_attr->n_ctl_lin) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "LAC Pixel Start Number", DFNT_INT32,
                1, (void *) &cz_attr->lac_pixl_start_no) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "LAC Pixel Subsampling", DFNT_INT32,
                1, (void *) &cz_attr->lac_pixl_subsample) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Scene Center Scan Line", DFNT_INT32,
                1, (void *) &cz_attr->cntr_scn_line) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Filled Scan Lines", DFNT_INT32,
                1, (void *) &cz_attr->filled_lines) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Sensor Tilt", DFNT_FLOAT32,
                1, (void *) &cz_attr->tilt) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Scene Gain", DFNT_INT32,
                1, (void *) &cz_attr->gain) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Thresh", DFNT_INT32,
                1, (void *) &cz_attr->thresh) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Calibration Slope", DFNT_FLOAT32,
                6, (void *) cz_attr->slope) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Calibration Intercept", DFNT_FLOAT32,
                6, (void *) cz_attr->intercept) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Center Roll", DFNT_FLOAT32,
                1, (void *) &cz_attr->roll) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Center Pitch", DFNT_FLOAT32,
                1, (void *) &cz_attr->pitch) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Center Yaw", DFNT_FLOAT32,
                1, (void *) &cz_attr->yaw) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "ILT Flags", DFNT_UINT8,
                1, (void *) &cz_attr->ilt_flags) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Number of Missing Scan Lines", DFNT_INT16,
                1, (void *) &cz_attr->n_miss_scans) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Number of Scans with Missing Channels",
                DFNT_INT16,
                6, (void *) cz_attr->n_scan_mis_chan) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Number of HDT Sync Losses", DFNT_INT16,
                1, (void *) &cz_attr->n_hdt_sync_loss) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Number of HDT Parity Errors", DFNT_INT16,
                1, (void *) &cz_attr->n_hdt_parity_err) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Number of WBVT Sync Losses", DFNT_INT16,
                1, (void *) &cz_attr->n_wbvt_sync_loss) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Number of WBVT Slip Occurrences", DFNT_INT16,
                1, (void *) &cz_attr->n_wbvt_slips) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Latitude Units", DFNT_CHAR,
                15, (void *) cz_attr->lat_unit) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Longitude Units", DFNT_CHAR,
                15, (void *) cz_attr->lon_unit) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Scene Center Latitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->center_lat) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Scene Center Longitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->center_lon) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Scene Center Solar Zenith", DFNT_FLOAT32,
                1, (void *) &cz_attr->cntr_sol_zen) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Upper Left Latitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->up_lft_lat) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Upper Left Longitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->up_lft_lon) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Upper Right Latitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->up_rgt_lat) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Upper Right Longitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->up_rgt_lon) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Lower Left Latitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->lo_lft_lat) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Lower Left Longitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->lo_lft_lon) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Lower Right Latitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->lo_rgt_lat) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Lower Right Longitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->lo_rgt_lon) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Northernmost Latitude", DFNT_FLOAT32,
                1, (void *) cz_attr->limits) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Southernmost Latitude", DFNT_FLOAT32,
                1, (void *) cz_attr->limits + 1) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Westernmost Longitude", DFNT_FLOAT32,
                1, (void *) cz_attr->limits + 2) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Easternmost Longitude", DFNT_FLOAT32,
                1, (void *) cz_attr->limits + 3) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Start Center Latitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->start_cntr_lat) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "Start Center Longitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->start_cntr_lon) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "End Center Latitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->end_cntr_lat) != 0) return -1;

        if (hdfio_rd_gattr(hdfid, "End Center Longitude", DFNT_FLOAT32,
                1, (void *) &cz_attr->end_cntr_lon) != 0) return -1;

        /*
         *  read in the SDSes
         */
        if (hdfio_rd_sd(hdfid, "slat", (void *) cz_dat->slat) != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "slon", (void *) cz_dat->slon) != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "clat", (void *) cz_dat->clat) != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "clon", (void *) cz_dat->clon) != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "elat", (void *) cz_dat->elat) != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "elon", (void *) cz_dat->elon) != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "tilt", (void *) cz_dat->tilt) != 0)
            return -1;

        tmp_ptr = (unsigned char *) (cz_dat->counts[0]);
        if (hdfio_rd_sd(hdfid, "band1", (void *) tmp_ptr) != 0)
            return -1;

        tmp_ptr = (unsigned char *) (cz_dat->counts[1]);
        if (hdfio_rd_sd(hdfid, "band2", (void *) tmp_ptr) != 0)
            return -1;

        tmp_ptr = (unsigned char *) (cz_dat->counts[2]);
        if (hdfio_rd_sd(hdfid, "band3", (void *) tmp_ptr) != 0)
            return -1;

        tmp_ptr = (unsigned char *) (cz_dat->counts[3]);
        if (hdfio_rd_sd(hdfid, "band4", (void *) tmp_ptr) != 0)
            return -1;

        tmp_ptr = (unsigned char *) (cz_dat->counts[4]);
        if (hdfio_rd_sd(hdfid, "band5", (void *) tmp_ptr) != 0)
            return -1;

        tmp_ptr = (unsigned char *) (cz_dat->counts[5]);
        if (hdfio_rd_sd(hdfid, "band6", (void *) tmp_ptr) != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "orb_vec", (void *) cz_dat->orb_vec) != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "att_ang", (void *) cz_dat->att_ang) != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "pos_err", (void *) cz_dat->pos_err) != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "cntl_pt_cols", (void *) cz_dat->ctl_pt_cols)
                != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "cntl_pt_rows", (void *) cz_dat->ctl_pt_rows)
                != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "latitude", (void *) cz_dat->ctl_pt_lat) != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "longitude", (void *) cz_dat->ctl_pt_lon) != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "gain", (void *) cz_dat->gain) != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "slope", (void *) cz_dat->slope) != 0)
            return -1;

        if (hdfio_rd_sd(hdfid, "intercept", (void *) cz_dat->intercept) != 0)
            return -1;

    }
    /*
     *  close the file and end
     */
    hdfio_close(hdfid);

    return 0;
}
