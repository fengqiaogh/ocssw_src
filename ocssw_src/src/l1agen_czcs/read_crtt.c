#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "l1czcs.h"
#include "navigation.h"

#define CTL_PT_INC 8  /* control point incriment (fixed for now) */
#define N_ANCHOR 77

int read_crtt(char *crtt_file, gattr_struc *gattr, l1_data_struc *l1_data)
/*******************************************************************

   read_crtt

   purpose: read the needed data from the CRTT raw data file

   Returns type: int - 0 no problems, else file open, read problems

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            crtt_file        I      name of the raw file
      gattr_struc *     gattr            O      structure with final 
                                                attributes
      l1_data_struc *   l1_data          O      arrays data counts, lat, lons

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       10-Mar-2004     Original development, based on
                                        read_l1_czcs program
      W. Robinson       2-Dec-2005      add capability to use Nimbus 7 orbit 
                                        data, calling fill_orb_dat

 *******************************************************************/ {
    char *get_file_specifier();
    char *get_file_name();
    char *file_name;
    int yr, doy;
    int i = 0, j, n_ctl, n_ctl_pt;
    FILE *fptr;
    /* def of the cal_czcs */
    void cal_czcs(int, gattr_struc *, l1_data_struc *);

    HEADER1_TYPE header1;
    HEADER2_TYPE header2;
    DATA_REC_TYPE data;

    short header2_size, /* size of header 2 record in blocks */
            data_rec_size, /* size of data record in blocks */
            offset_header, /* offset to first header record in blocks */
            offset_data, /* offset to first data record in blocks */
            num_of_data_recs;

    /* get file name: */
    file_name = basename(crtt_file);

    sscanf(file_name, "%2d%3d", &yr, &doy);
    yr = yr + 1900;
    fptr = fopen(crtt_file, "r");
    if (!(fptr = fopen(crtt_file, "r"))) {
        printf("File \"%s\" not found.\n", crtt_file);
        exit(-1);
    }

    /* read VAX header and get record info from it: */

    if (fread(&header1, sizeof ( header1), 1, fptr))
        get_record_info(header1, &header2_size, &data_rec_size,
            &offset_header, &offset_data, &num_of_data_recs);
    else {
        printf("Could not read VAX header.\nProgram terminated.");
        return -1;
    }

    /*
     *  seek, read and print first document record: 
     *  Nothing is used from here, but keep for now
     */

    fseek(fptr, offset_header * 1024, 0); /* ie, look for document  */
    /*  rec at a displacement of    */
    /* offset_header - 1 * 512 bytes*/
    /* from begining of file        */
    if (fread(&header2, header2_size * 512, 1, fptr)) {
#ifdef __DEC_MACHINE
        convert_header_rec_to_dec(&header2);
#endif
        hdr_2_gattr(header2, gattr);
    } else {
        printf("\nNote:  could not read first documentation record.\n");
        return -1;
    }

    /* seek, read and print trailing document record: */

    fseek(fptr, (offset_data * 512) +
            (num_of_data_recs * data_rec_size * 512), 0);
    if (fread(&header2, header2_size * 512, 1, fptr)) {
#ifdef __DEC_MACHINE
        convert_header_rec_to_dec(&header2);
#endif
        hdr_2_gattr(header2, gattr);
    } else {
        printf("\nError: Could not read trailing documentation record.\n");
        return -1;
    }
    /*
     *  compute the # of control points
     */
    n_ctl = 1 + (NCZCS_PIX - 1) / CTL_PT_INC;
    if ((n_ctl - 1) * CTL_PT_INC == NCZCS_PIX + 1)
        n_ctl_pt = n_ctl;
    else
        n_ctl_pt = n_ctl + 1;

    gattr->n_ctl_pt = n_ctl_pt;
    gattr->ctl_pt_incr = CTL_PT_INC;
    /*
     *  The gattr->scan_lines can be wrong and != num_of_data_recs
     *  So, for now, if not =, note and use num_of_data_recs
     *  (if num_of_data_recs wrong, it couldn't get trailing header and would die)
     */
    if (num_of_data_recs != gattr->scan_lines) {
        printf("\nWarning, # scan lines in trailing header incorrect\nUsing VAX header value\n");
        gattr->scan_lines = num_of_data_recs;
    }

    /*
     *  set up space for all the data arrays and set the control point 
     *  column values
     */
    cz_dat_alloc(num_of_data_recs, n_ctl_pt, 0, l1_data);

    for (j = 0; j < n_ctl_pt; j++) {
        l1_data->ctl_pt_cols[ j ] = j * CTL_PT_INC + 1;
        if (j == (n_ctl_pt - 1)) l1_data->ctl_pt_cols[ j ] = NCZCS_PIX;
    }
    /*
     *  initialize some nav values needed in ftn calls in czcs_ctl_pt
     */
    cdata_();

    for (i = 0; i < num_of_data_recs; i++) {
        l1_data->ctl_pt_rows[i] = i + 1;
        fseek(fptr, (offset_data + (data_rec_size * i)) * 512, 0);
        if (fread(&data, data_rec_size * 512, 1, fptr)) {
#ifdef __DEC_MACHINE
            convert_data_rec_to_dec(&data);
#endif
            /*
             * deposit the line of counts
             * this will make pixel, line, band organization with pixel the fastest
             */
            for (j = 0; j < 6; j++) {
                memcpy(l1_data->counts[j] + NCZCS_PIX * i,
                        *(data.radiance_cnts) + j * NCZCS_PIX, NCZCS_PIX);
            }
            /*
             * set up the times needed from the data records
             * (needs to be used in czcs_ctl_pt)
             */
            l1_data->msec[i] = data.msec_of_day;
            /*
             * take the anchor points from the file, correct them and 
             * make an evenly spaced control point grid
             */
            czcs_ctl_pt(data, gattr, i, l1_data);
            /*
             * set up the quality data
             */
            l1_data->cal_sum[ i * 5 ] = data.info.qstnble_ephemeris;
            l1_data->cal_sum[ i * 5 + 1 ] = data.info.qstnble_attitude;
            l1_data->cal_sum[ i * 5 + 2 ] = data.info.channel_not_present;
            l1_data->cal_sum[ i * 5 + 3 ] = data.info.cal_outside_range;
            l1_data->cal_sum[ i * 5 + 4 ] = data.info.volt_outside_range;

            for (j = 0; j < 6; j++)
                l1_data->cal_scan[ i * 6 + j ] = data.cal_qual[j];

            /*
             *  if enabled, calibrate the counts
             */
#ifdef GEOM_CAL
            cal_czcs(i, gattr, l1_data);
#endif
        } else {
            printf("Error:  Could not read data record # %d.\n", i);
            return -1;
        }
    }
    /*
     * the last data record provides the end time
     */
    time_str(data.year, data.doy, data.msec_of_day, gattr->end_time);
    gattr->end_year = data.year;
    gattr->end_day = data.doy;
    gattr->end_msec = data.msec_of_day;
    /*
     *  Set up some lat, lon info using improved control points
     *  and set line-by-line info from ctl pts and global attrs
     *  Also, get orbit info in using SMMR-derived Nimbus 7 orbit
     */
    cz_ll_upd(l1_data, gattr);
    cz_sd_set(l1_data, gattr);
    if (fill_orb_dat(l1_data, gattr) != 0) return -1;
    return 0;
}

float fixed_pt_2_floating_pt(int num, int position)

/* returns the floating point value of a fixed pt #.  The MSB indicates 
the sign, and the parameter "position" is the position of the decimal pt 
(0 means to the right of the LSB, 1 means between the LSB and second to 
the least most significant bit). 
 */
 {

    int divisor = 1;

    divisor = divisor << position;

    return ((float) /*(1 - (2 * ((num & 0x80000000) != 0))) **/ num / divisor);
}

void get_record_info(HEADER1_TYPE header1, short *hdrsz, short *datarecsz,
        short *offsethdr, short *offsetdata, short *numrecs)

/*
Get the following quantities from the VAX header recorder:
     document record size (hdrsz)
     data record size
     offset to the first document record (offsethdr)
     offset to the first data record
     number of data records (numrecs).
If __DEC_MACHINE is not defined, reverse the byte order of the 2-byte
words.
 */
 {
#ifdef __DEC_MACHINE
    *hdrsz = (header1.doc_rec_len / 512) +
            ((header1.doc_rec_len % 512) != 0);

    *datarecsz = (header1.data_rec_len / 512) +
            ((header1.data_rec_len % 512) != 0);

    *offsethdr = header1.header_rec_offset;
    *offsetdata = header1.offset_1st_data_rec;
    *numrecs = header1.num_of_data_recs;
#else             /* reverse byte order */
    short temp;
    temp = header1.doc_rec_len;
    temp = reverse_short_int(temp);
    *hdrsz = (temp / 512) + ((temp % 512) != 0);

    temp = header1.data_rec_len;
    temp = reverse_short_int(temp);
    *datarecsz = (temp / 512) + ((temp % 512) != 0);

    *offsethdr = reverse_short_int(header1.header_rec_offset);
    *offsetdata = reverse_short_int(header1.offset_1st_data_rec);
    *numrecs = reverse_short_int(header1.num_of_data_recs);
#endif
}

short reverse_short_int(short num)
/* returns the parameter passed with the byte positions reversed. */ {

    union {
        unsigned short i;
        char c[2];
    } short_struct;
    char temp;

    short_struct.i = num;
    temp = short_struct.c[0];
    short_struct.c[0] = short_struct.c[1];
    short_struct.c[1] = temp;

    return (short_struct.i);
}

int32_t reverse_long_int(int32_t num)
/* returns the parameter passed with the byte positions reversed. */ {

    union {
        int i;
        char c[4];
    } long_struct;
    char temp;

    long_struct.i = num;
    temp = long_struct.c[0];
    long_struct.c[0] = long_struct.c[3];
    long_struct.c[3] = temp;
    temp = long_struct.c[1];
    long_struct.c[1] = long_struct.c[2];
    long_struct.c[2] = temp;

    return (long_struct.i);
}


#ifdef __DEC_MACHINE

void convert_data_rec_to_dec(DATA_REC_TYPE *rec)
/* converts integers in data record to dec format--ie, reverses
byte order.
 */
 {
    int i;

    union {
        BIT_FLD_2_TYPE bit_struct;
        unsigned short words[2];
        int i;
    } u;

    u.bit_struct = rec->info;
    u.words[0] = reverse_short_int(u.words[0]);
    rec->info.p_rec_num = u.words[0] >> 4;
    rec->info.sp1 = u.words[0] & 0x00f;

    rec->seq_no = reverse_short_int(rec->seq_no);
    rec->year = reverse_short_int(rec->year);
    rec->doy = reverse_short_int(rec->doy);
    rec->subcommuted_data_val_cnt =
            reverse_short_int(rec->subcommuted_data_val_cnt);
    for (i = 0; i < 6; i++)
        rec->cal_lamp_rad[i] = reverse_short_int(rec->cal_lamp_rad[i]);
    rec->blk_bdy_tmp_cnt = reverse_short_int(rec->blk_bdy_tmp_cnt);
    rec->wbvt_bit_slips_summary =
            reverse_short_int(rec->wbvt_bit_slips_summary);
    rec->hdt_sync_losses = reverse_short_int(rec->hdt_sync_losses);
    rec->hdt_parity_errs = reverse_short_int(rec->hdt_parity_errs);
    rec->wbvt_bit_slips = reverse_short_int(rec->wbvt_bit_slips);
    rec->pxl_no_nadir = reverse_short_int(rec->pxl_no_nadir);

    rec->msec_of_day = reverse_long_int(rec->msec_of_day);
    for (i = 0; i < 77; i++)
        rec->lat_anchr_pts[i] = reverse_long_int(rec->lat_anchr_pts[i]);
    for (i = 0; i < 77; i++)
        rec->long_anchr_pts[i] = reverse_long_int(rec->long_anchr_pts[i]);
    for (i = 0; i < 28; i++)
        rec->spare3[i] = reverse_long_int(rec->spare3[i]);

}
#endif


#ifdef __DEC_MACHINE

void convert_header_rec_to_dec(HEADER2_TYPE *rec)
/*
 *  Revision history:
 *   W. Robinson, SAIC, 7 Jun 2004  fix setting of the no_of_missing_scans
 *        value so {1] - [5] is correct
 */
/* converts integers in header record to dec format--ie, reverses
byte order.
 */
 {
    int i;

    union {
        BIT_FLD_1_TYPE bit_struct;
        unsigned short words[2];
        int i;
    } u;

    u.bit_struct = rec->file_info;
    u.words[0] = reverse_short_int(u.words[0]);
    rec->file_info.p_rec_num = u.words[0] >> 4;
    rec->file_info.sp1 = u.words[0] & 0x00f;

    rec->tape_seq_no = reverse_long_int(rec->tape_seq_no);
    rec->film_frame_no = reverse_long_int(rec->film_frame_no);
    rec->starting_yr = reverse_short_int(rec->starting_yr);
    rec->starting_day = reverse_short_int(rec->starting_day);

    rec->start_msec = reverse_long_int(rec->start_msec);
    rec->inc_in_msec = reverse_long_int(rec->inc_in_msec);

    rec->orbit = reverse_short_int(rec->orbit);
    rec->no_of_scans = reverse_short_int(rec->no_of_scans);
    rec->lat_cntr = reverse_short_int(rec->lat_cntr);
    rec->long_cntr = reverse_short_int(rec->long_cntr);
    rec->lat_crnr_fitl = reverse_short_int(rec->lat_crnr_fitl);
    rec->long_crnr_fitl = reverse_short_int(rec->long_crnr_fitl);
    rec->lat_crnr_fitr = reverse_short_int(rec->lat_crnr_fitr);
    rec->long_crnr_fitr = reverse_short_int(rec->long_crnr_fitr);
    rec->lat_crnr_litl = reverse_short_int(rec->lat_crnr_litl);
    rec->long_crnr_litl = reverse_short_int(rec->long_crnr_litl);
    rec->lat_crnr_litr = reverse_short_int(rec->lat_crnr_litr);
    rec->long_crnr_litr = reverse_short_int(rec->long_crnr_litr);
    rec->no_of_missing_scans_all =
            reverse_short_int(rec->no_of_missing_scans_all);

    for (i = 0; i < 6; i++)
        rec->no_of_missing_scans[i] =
            reverse_short_int(rec->no_of_missing_scans[i]);
    rec->decom_run_no = reverse_long_int(rec->decom_run_no);
    rec->decom_reel_no = reverse_long_int(rec->decom_reel_no);

    rec->no_of_hdt_sync_losses = reverse_short_int(rec->no_of_hdt_sync_losses);
    rec->no_of_hdt_parity_errs = reverse_short_int(rec->no_of_hdt_parity_errs);
    rec->no_of_wbvt_sync_losses = reverse_short_int(rec->no_of_wbvt_sync_losses);
    rec->no_of_wbvt_bit_slips = reverse_short_int(rec->no_of_wbvt_bit_slips);
    for (i = 0; i < 32; i++)
        rec->ave_subcomputed_data[i] =
            reverse_short_int(rec->ave_subcomputed_data[i]);
    rec->baseplate_temp = reverse_short_int(rec->baseplate_temp);
    rec->tilt = reverse_short_int(rec->tilt);

    for (i = 0; i < 134; i++)
        rec->spare4[i] = reverse_long_int(rec->spare4[i]);

    rec->scene_cntr_yr = reverse_short_int(rec->scene_cntr_yr);
    rec->scene_cntr_doy = reverse_short_int(rec->scene_cntr_doy);
    rec->scene_cntr_msec = reverse_long_int(rec->scene_cntr_msec);

    rec->solar_el_cntr = reverse_short_int(rec->solar_el_cntr);
    rec->solar_az_cntr = reverse_short_int(rec->solar_az_cntr);
    rec->cntr_roll = reverse_short_int(rec->cntr_roll);
    rec->cntr_pitch = reverse_short_int(rec->cntr_pitch);
    rec->cntr_yaw = reverse_short_int(rec->cntr_yaw);
    rec->top_lft_tick_label = reverse_short_int(rec->top_lft_tick_label);
    rec->top_rgt_tick_label = reverse_short_int(rec->top_rgt_tick_label);
    rec->bot_lft_tick_label = reverse_short_int(rec->bot_lft_tick_label);
    rec->bot_rgt_tick_label = reverse_short_int(rec->bot_rgt_tick_label);
    rec->lft_top_tick_label = reverse_short_int(rec->lft_top_tick_label);
    rec->lft_bot_tick_label = reverse_short_int(rec->lft_bot_tick_label);
    rec->rgt_top_tick_label = reverse_short_int(rec->rgt_top_tick_label);
    rec->rgt_bot_tick_label = reverse_short_int(rec->rgt_bot_tick_label);

    for (i = 0; i < 27; i++)
        rec->top_tick_loc[i] = reverse_short_int(rec->top_tick_loc[i]);
    for (i = 0; i < 27; i++)
        rec->bot_tick_loc[i] = reverse_short_int(rec->bot_tick_loc[i]);
    for (i = 0; i < 27; i++)
        rec->lft_tick_loc[i] = reverse_short_int(rec->lft_tick_loc[i]);
    for (i = 0; i < 27; i++)
        rec->rgt_tick_loc[i] = reverse_short_int(rec->rgt_tick_loc[i]);

    for (i = 0; i < 6; i++) {
        rec->chan_line[i].slope = reverse_long_int(rec->chan_line[i].slope);
        rec->chan_line[i].intercept =
                reverse_long_int(rec->chan_line[i].intercept);
    }

    for (i = 0; i < 256; i++)
        rec->conversion_tab_chan6[i] =
            reverse_short_int(rec->conversion_tab_chan6[i]);

    for (i = 0; i < 6; i++) {
        rec->enhancement_eqs[i].slope =
                reverse_short_int(rec->enhancement_eqs[i].slope);
        rec->enhancement_eqs[i].intercept =
                reverse_short_int(rec->enhancement_eqs[i].intercept);
    }
    for (i = 0; i < 2; i++)
        rec->spare5[i] = reverse_long_int(rec->spare5[i]);
    for (i = 0; i < 945; i++)
        rec->czcs_1lt[i] = reverse_long_int(rec->czcs_1lt[i]);

}
#endif
