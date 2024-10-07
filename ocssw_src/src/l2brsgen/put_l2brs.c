#include "l2brsgen.h"

#include <get_product_table.h>

#include <math.h>
#include <strings.h>
#include <png.h>
#include <setupflags.h>
#include <hdf.h>

static float32 base = BASE_VAL, slope, intercept;
static char sc_eqn[128], sc_type[64], param_val[128], units_val[32];

/*-----------------------------------------------------------------------------

    Function: put_l2brs

    Returns: int32(status)
        The return code is FAIL (-1) if an error occurs, SUCCEED (0)
        otherwise.

    Description:
        The function put_l2brs creates a level 2 Browse data file with
        the file name given by l2brs_path.  Level 2  file meta data, palette,
        and data will be written.  See "INTERFACE SPECIFICATIONS FOR SeaWiFS
        OPERATIONAL PRODUCT INPUT AND OUTPUT SOFTWARE" (By Fred Patt,
        Jim Firestone, and Mike Darzi).

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      char *    l2brs_path  I     directory path & file name for browse product
      char *    replaces    I     filename of previously generated product
                                  that current product is intended to replace
      char *    ptime       I     processing start time
      char *    infiles     I     name of input product
      int32     px_start    I     scan data start pixel (1 based)
      int32     px_end      I     scan data end pixel (1 based)
      int32     px_subsamp  I     pixel subsampling rate
      int32     sc_start    I     start scan line (1 based)
      int32     sc_end      I     end scan line (1 based)
      int32     sc_subsamp  I     scan subsampling rate
      char *    l2brs_name  I     name of the geophysical parameter to which
                                  l2brs_data corresponds
      float32 * l2brs_data  I     image data for parameter l2brs_name
      int16 *   l2brs_flags I     level2 data flags
      char  *   flag_names  I     list of algorithm names
      char  *   mskflg      I     names of all flag algorithms used as masks
      usigned char *palette I     RGB wts for each of 256 gray-levels of the
                                  l2brs_data byte values
      int32     px_ll_num   I     number of lat/lon pairs in px_ll_first and
                                        px_ll_last
      float32 * px_ll_first I     lat/lon pairs of pixels along the first scan
      float32 * px_ll_last  I     lat/lon pairs of pixels along the last scan
      int32     sc_ll_num   I     number of lat/lon pairs in sc_ll_first and
                                   sc_ll_last
      float32 * sc_ll_first I     lat/lon pairs of scan line start points
      float32 * sc_ll_last  I     lat/lon pairs of scal line end points
      char *    proc_con    I     processing control input
      char *    proc_log    I     processing log output
      int16     syear       I     yr of data start time; as ret by get_l2_open
      int16     sday        I     day-of-year for data start time; as ret by
                                        get_l2_open
      int32     smsec       I     millisecs of day for data start time; as ret
                                        by get_l2_open
      int16     eyear       I     yr data end time; as ret by get_l2_open
      int16     eday        I     day-of-year for data end time; as ret by
                                        get_l2_open
      int16     emsec       I     millisecs of day for data end time; as ret by
                                        get_l2_open
      char *    dtype       I     data type flag; as ret by get_l2_oepn
      int32     nrec        I     no. of scan lines; as ret by get_l2_open
      int32     nsamp       I     pixels/scan line; as ret by get_l2_open
      int32     ntilts      I     Number of tilt states in the scene; as ret by
                                        get_l2_open
      uint8 *   tilt_flags  I     tilt flags corres. to eh tilt state in the
                                        scene; as ret by get_l2_open
      int16 *   tilt_ranges I     first and last scan line nos. corrsp. to eh.
                                  tilt state in the scee; as ret by get_l2_open
      float32 * tilt_lats   I     lats of the end pixels for the scan lines of
                                  tilt_ranges; as ret by get_l2_open
      float32 * tilt_lons   I     lons of the end pixels for the scan lines of
                                  tilt_ranges; as ret by get_l2_open
      meta_l2Type *meta_l2  I     L2 product metadata read by get_l2_open

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/23/94    Original development
        Lakshmi Kumar    Hughes STX      06/16/94    Modified according to
                                                     Interface specs v3.1
        Lakshmi Kumar    Hughes STX      10/06/94    Modified according to
                                                     Interface specs v3.3
        Lakshmi Kumar 	 Hughes STX      01/24/94    Corrected ctime/ntime
                                                     datatype to character type
        Lakshmi Kumar	 HITC		 05/19/95    Added code to write
                                                     Navigation vgp, its data
                                                     sets and tilt vgroup
        Lakshmi Kumar	 Hughes STX	 07/25/95    Added code to define gray
                                                     levels between 250 & 255
        Lakshmi Kumar	 Hughes STX	 11/01/95    Added legend, start and
                                                     end center coordinates as
                                                     global attributes (V4.4)
                                                     Sets l2brs pix to 0 if
                                                     mskflg indicates some flag
                                                     & masking bit is turned on
                                                     for that pixel
        Lakshmi Kumar    Hughes STX	 02/15/96    Added some comments and
                                                     deleted a redundent if
                                                     check
        Lakshmi Kumar    Hughes STX      02/23/96    Added gray level 252 for
                                                     representing invalid data
        Lakshmi Kumar    Hughes STX      06/18/96    Changed defn. of MAX to
                                                     MAXVAL inorder to remove
                                                     compile time warning
        Lakshmi Kumar	 Hughes STX	 11/01/96    Gray levels from 251-255
                                                     are reset as shown below:
                                                     255 - missing data
                                                     253 - land
                                                     254 - cloud & ice
                                                     252 - glint
                                                     251 - invalid data
                                                     (Ref. V3.0 product specs.)
        Lakshmi Kumar	 Hughes STX	 08/08/97    ff_mising and sdps_missing
                                                     flags needed casting as
                                                     their input values were
                                                     bytes and output values
                                                     were 4bytes int.
        Lakshmi Kumar	 Hughes STX	 10/16/97    updated to add upper left,
                                                     right and lower left,
                                                     right lat/lon attributes
        Norman Kuring    NASA/GSFC	 11/12/97    Mask stray-light pixels
        Norman Kuring	 NASA/GSFC	 11/14/97    More mask fiddling
-----------------------------------------------------------------------------*/
int32 put_l2brs(char *l2brs_path, char *replaces, char *ptime, char *infiles,
        int32 px_start, int32 px_end, int32 px_subsamp, int32 brs_nsamp,
        int32 sc_start, int32 sc_end, int32 sc_subsamp, int32 brs_nrec,
        char *l2brs_name, float32 *l2brs_data, int32 *l2brs_flags,
        char *flag_names,int32_t * l2_bits, char *flaguse, unsigned char *palette,
        float32 *px_ll_first, float32 *px_ll_last, float32 *sc_ll_first,
        float32 *sc_ll_last, char *proc_con, int16 syear, int16 sday,
        int32 smsec, int16 eyear, int16 eday, int32 emsec, char *dtype,
        int32 nrec, int32 nsamp, int32 ntilts, short *tilt_flags,
        int16 *tilt_ranges, int16 *cntl_pt_lat, int16 *cntl_pt_lon,
        meta_l2Type *meta_l2, product_table_t *ptable_rec, const char* oformat, int32 apply_pal) {

    int32 i, fid, sdfid, sdsid, rank = 2, dimsizes[3], start[3];

    uint32 mask_253, mask_254;
    double scaled_value;
    unsigned char *brs_data;
    uint32 required_mask = 0;


    FILE *outfp = NULL;

    if (strcasecmp(ptable_rec->scaling, "logarithmic") == 0) {
        intercept = log10(ptable_rec->min);
        slope = (log10(ptable_rec->max) - intercept) / 250.0;
        strcpy(sc_type, "logarithmic");
        strcpy(sc_eqn, "Base**((Slope*brs_data) + Intercept) = ");
        strcat(sc_eqn, ptable_rec->description);
    } else {
        intercept = ptable_rec->min;
        slope = (ptable_rec->max - ptable_rec->min) / 250.0;
        strcpy(sc_type, "linear");
        strcpy(sc_eqn, "Slope*brs_data + Intercept = ");
        strcat(sc_eqn, ptable_rec->description);
    }
    strcpy(param_val, ptable_rec->description);
    strcpy(units_val, ptable_rec->units);

    brs_data = (unsigned char *) malloc(sizeof (unsigned char) * brs_nsamp * brs_nrec);
    if (brs_data == NULL) {
        fprintf(stderr, "\nput_l2brs: Unable to allocate space for scaling brs data\n");
        return FAIL;
    }

    // land mask
    mask_253 = 1 << 1;
    mask_254 = 0;
    required_mask = 0;

    /*  Check masks against their names and set the masking fld "xmask" */
    if (flag_names == NULL) {
        fprintf(stderr, "\nError: put_l2brs: flag_names = null string");
        fprintf(stderr, "\n\t masking will not be applied ");
    } else {
        if (flaguse != NULL && flaguse[0] != '\0') {
            int status;
            setupflags(flag_names, flaguse, &mask_254, &required_mask, &status,l2_bits);
        }
    }

    for (i = 0; i < brs_nsamp * brs_nrec; i++) {
        if (l2brs_flags[i] & mask_253)
            brs_data[i] = 253;
        else if (l2brs_flags[i] & mask_254)
            brs_data[i] = 254;
        else if ((l2brs_flags[i] & required_mask) != required_mask)
            brs_data[i] = 254;
        else if (l2brs_data[i] < ptable_rec->min)
            brs_data[i] = 0;
        else {
            if (strcasecmp(ptable_rec->scaling, "logarithmic") == 0)
                scaled_value = floor((log10(l2brs_data[i]) - intercept) / slope + 0.5);
            else
                scaled_value = floor((l2brs_data[i] - intercept) / slope + 0.5);

            if (scaled_value > 250.0)
                brs_data[i] = 250;
            else
                brs_data[i] = (unsigned char) scaled_value;
        }
    }

    if (strcmp(oformat, "HDF4") == 0) { // hdf
        sdfid = SDstart(l2brs_path, DFACC_CREATE);
        fid = Hopen(l2brs_path, DFACC_RDWR, 0);
        Vstart(fid);

        write_attrs(sdfid, l2brs_path, replaces, ptime, infiles, px_start, px_end,
                px_subsamp, brs_nsamp, sc_start, sc_end, sc_subsamp, brs_nrec,
                l2brs_name, proc_con, syear, sday, smsec, eyear, eday, emsec, dtype,
                nrec, nsamp, meta_l2);

        write_image(l2brs_path, brs_data, brs_nsamp, brs_nrec, palette);

        free(brs_data);

        /*   write px_ll_first SDS */
        start[0] = start[1] = start[2] = 0;
        dimsizes[0] = brs_nsamp;
        dimsizes[1] = 2;
        if ((sdsid = write_SDS(sdfid, PX_LL_FIRST, DFNT_FLOAT32, rank, dimsizes,
                start, px_ll_first)) < 0)
            return FAIL;
        SDsetattr(sdsid, LONGNAME, DFNT_CHAR, strlen(PX_LL_FST_ATTR) + 1,
                (VOIDP) PX_LL_FST_ATTR);

        SDsetdimname(SDgetdimid(sdsid, 0), "Pixels per Scan");
        SDsetdimname(SDgetdimid(sdsid, 1), "Lat/Lon");

        /*   write px_ll_last SDS */

        dimsizes[0] = brs_nsamp;
        dimsizes[1] = 2;
        if ((sdsid = write_SDS(sdfid, PX_LL_LAST, DFNT_FLOAT32, rank, dimsizes,
                start, px_ll_last)) < 0)
            return FAIL;
        SDsetattr(sdsid, LONGNAME, DFNT_CHAR, strlen(PX_LL_LST_ATTR) + 1,
                (VOIDP) PX_LL_LST_ATTR);

        SDsetdimname(SDgetdimid(sdsid, 0), "Pixels per Scan");
        SDsetdimname(SDgetdimid(sdsid, 1), "Lat/Lon");

        /*   write sc_ll_first SDS */

        dimsizes[0] = brs_nrec;
        dimsizes[1] = 2;
        if ((sdsid = write_SDS(sdfid, SC_LL_FIRST, DFNT_FLOAT32, rank, dimsizes,
                start, sc_ll_first)) < 0)
            return FAIL;
        SDsetattr(sdsid, LONGNAME, DFNT_CHAR, strlen(SC_LL_FST_ATTR) + 1,
                (VOIDP) SC_LL_FST_ATTR);

        SDsetdimname(SDgetdimid(sdsid, 0), "Number of Scans");
        SDsetdimname(SDgetdimid(sdsid, 1), "Lat/Lon");

        /*   write sc_ll_last SDS */

        dimsizes[0] = brs_nrec;
        dimsizes[1] = 2;
        if ((sdsid = write_SDS(sdfid, SC_LL_LAST, DFNT_FLOAT32, rank, dimsizes,
                start, sc_ll_last)) < 0)
            return FAIL;
        SDsetattr(sdsid, LONGNAME, DFNT_CHAR, strlen(SC_LL_LST_ATTR) + 1,
                (VOIDP) SC_LL_LST_ATTR);

        SDsetdimname(SDgetdimid(sdsid, 0), "Number of Scans");
        SDsetdimname(SDgetdimid(sdsid, 1), "Lat/Lon");

        /*** create tilt vgroup and its datasets */
        if (ntilts > 0) {
            if ((write_tilt_sets(fid, sdfid, ntilts, tilt_flags, tilt_ranges)) < 0)
                return FAIL;
        }

        /*** create navigation vgroup and its datasets */
        if ((write_nav_sets(fid, sdfid, brs_nrec, brs_nsamp, cntl_pt_lat, cntl_pt_lon)) < 0)
            return FAIL;

        SDend(sdfid);
        Vend(fid);
        Hclose(fid);

        // end outmode = 0 = hdf

    } else if (strcmp(oformat, "PPM") == 0) {
        if ((outfp = fopen(l2brs_path, "w")) == NULL) {
            fprintf(stderr, "put_l2brs: Error: Unable to open %s for writing.\n", l2brs_path);
            exit(1);
        }

        if (apply_pal) {

            /* Write ppm header */
            fprintf(outfp, "P6\n");
            fprintf(outfp, "%d %d\n", brs_nsamp, brs_nrec);
            fprintf(outfp, "255\n");
            for (i = 0; i < brs_nsamp * brs_nrec; i++) {
                fwrite(palette + brs_data[i]*3, 1, 3, outfp);
            }

        } else {

            /* Write pgm header */
            fprintf(outfp, "P5\n");
            fprintf(outfp, "%d %d\n", brs_nsamp, brs_nrec);
            fprintf(outfp, "255\n");

            for (i = 0; i < brs_nsamp * brs_nrec; i++) {
                fwrite(brs_data + i, 1, 1, outfp);
            }
        }

        fclose(outfp);

        //   outmode = 1 = ppm or pgm 

    } else if (strcmp(oformat, "PNG") == 0) { // png

        outfp = fopen(l2brs_path, "w");
        if (!outfp) {
            fprintf(stderr, "put_l2brs: Error: Unable to open %s for writing.\n", l2brs_path);
            exit(1);
        }

        png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                NULL, NULL, NULL);
        if (!png_ptr) {
            fprintf(stderr, "put_l2brs: Error: Unable to create PNG write structure.\n");
            exit(1);
        }

        png_infop info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr) {
            png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
            fprintf(stderr, "put_l2brs: Error: Unable to create PNG info structure.\n");
            exit(1);
        }
        if (setjmp(png_jmpbuf(png_ptr))) {
            png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
            fprintf(stderr, "put_l2brs: Error: Unable to call PNG setjmp().\n");
            exit(1);
        }
        png_init_io(png_ptr, outfp);

        if (apply_pal) {
            // color
            png_set_IHDR(png_ptr, info_ptr, brs_nsamp, brs_nrec,
                    8, PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE,
                    PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
            png_set_PLTE(png_ptr, info_ptr, (png_const_colorp) palette, 256);

        } else {
            // Grayscale
            png_set_IHDR(png_ptr, info_ptr, brs_nsamp, brs_nrec,
                    8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
                    PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
        }

        uint8 * row_pointers[brs_nrec];
        for (i = 0; i < brs_nrec; i++)
            row_pointers[i] = brs_data + (i * brs_nsamp);
        png_set_rows(png_ptr, info_ptr, row_pointers);

        png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

        /* clean up after the write, and free any memory allocated */
        png_destroy_write_struct(&png_ptr, (png_infopp) NULL);

        fclose(outfp);

        // end 2 = png

    } else {
        fprintf(stderr, "put_l2brs: Error: oformat=%s is not supported.\n", oformat);
        exit(1);
    }

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function: write_attrs

    Returns: void

    Description:
        The function write_attrs writes metadata (global attributes).

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     sdfid       I     file ID
      char *    l2brs_path  I     directory path & file name for browse product
      char *    replaces    I     filename of previously generated product
                                  that current product is intended to replace
      char *    ptime       I     processing start time
      char *    infiles     I     name of input product
      int32     px_start    I     scan data start pixel
      int32     px_subsamp  I     pixel subsampling rate
      int32     px_num      I     no. of pixels in the browse image
      int32     sc_start    I     start scan line
      int32     sc_subsamp  I     scan subsampling rate
      int32     sc_num      I     no. of scan lines in the browse image
      char *    l2brs_name  I     name of the geophysical parameter to which
                                  l2brs_data corresponds
      int32     px_ll_num   I     number of lat/lon pairs in px_ll_first and
                                        px_ll_last
      int32     sc_ll_num   I     number of lat/lon pairs in sc_ll_first and
                                   sc_ll_last
      char *    proc_con    I     processing control input
      char *    proc_log    I     processing log output
      int16     syear       I     yr of data start time; as ret by get_l2_open
      int16     sday        I     day-of-year for data start time; as ret by
                                        get_l2_open
      int32     smsec       I     millisecs of day for data start time; as ret
                                        by get_l2_open
      int16     eyear       I     yr data end time
      int16     eday        I     day-of-year for data end time
      int16     emsec       I     millisecs of day for data end time
      char *    dtype       I     data type flag
      int32     nrec        I     no. of scan lines
      int32     nsamp       I     pixels/scan line
      meta_l2Type *meta_l2  I     L2 product metadata read by get_l2_open

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/23/94    Original development
        Lakshmi Kumar    Hughes STX      06/16/94    Made changes to agree with
                                                     Interface specs v3.1
        Lakshmi Kumar    Hughes STX      10/06/94    Made changes to agree with
                                                     Interface specs v3.3
        Lakshmi Kumar 	 Hughes STX      01/24/94    Corrected ctime/ntime
                                                     datatype to character type
        Lakshmi Kumar	 Hughes STX	 10/24/95    Added legend, start and
                                                     end center coordinates as
                                                     global attributes (V4.4)
        Lakshmi Kumar	 Hughes STX	 10/16/97    Added sc_ll_first and
                                                     sc_ll_last as input params.
        Norman Kuring    NASA/GSFC	 11/26/97    Removed hard-coded data-
                                                     type from legend attribute
        Norman Kuring	 NASA/GSFC	 02/09/97    Removed sc_ll_first and
                                                     sc_ll_last as inputs; get
                                                     corner coordinates from
                                                     meta_l2 structure instead.
-----------------------------------------------------------------------------*/

void write_attrs(int32 sdfid, char *l2brs_path, char *replaces, char *ptime,
        char *infiles, int32 px_start, int32 px_end, int32 px_subsamp,
        int32 brs_nsamp, int32 sc_start, int32 sc_end, int32 sc_subsamp,
        int32 brs_nrec, char *l2brs_name, char *proc_con, int16 syear,
        int16 sday, int32 smsec, int16 eyear, int16 eday, int32 emsec,
        char *dtype, int32 nrec, int32 nsamp, meta_l2Type *meta_l2) {
    div_t quot1, quot2, quot3;
    char string[20], pname[2048], *p, *pinfile;
    int32 lac_px_st, lac_px_subsamp;
    char legend[132];

    char *soft_name = "l2brsgen";
    char *soft_ver = "5.1";

    /*  write out global attributes  */
    p = strrchr(l2brs_path, '/');
    if (p == NULL)
        strcpy(pname, l2brs_path);
    else
        strcpy(pname, ++p);

    SDsetattr(sdfid, L2B_PNAME, DFNT_CHAR, strlen(pname) + 1, (VOIDP) pname);

    strcpy(pname, meta_l2->sensor_name);
    strcat(pname, TITLE_VAL);
    SDsetattr(sdfid, L2B_TITLE, DFNT_CHAR, strlen(pname) + 1, pname);
    SDsetattr(sdfid, SENSOR_NAME, DFNT_CHAR, strlen(meta_l2->sensor_name) + 1,
            meta_l2->sensor_name);

    sprintf(legend,
            "%s Level-2 %s %s browse data, day %3.3d, %4.4d",
            meta_l2->sensor_name, dtype, param_val, sday, syear);
    SDsetattr(sdfid, LEGEND, DFNT_CHAR, strlen(legend) + 1, (VOIDP) legend);

    SDsetattr(sdfid, DCENTER, DFNT_CHAR, strlen(DCENTER_VAL) + 1, DCENTER_VAL);
    SDsetattr(sdfid, SOFT_NAME, DFNT_CHAR, strlen(soft_name) + 1,
            (VOIDP) soft_name);
    SDsetattr(sdfid, SOFT_VER, DFNT_CHAR, strlen(soft_ver) + 1,
            (VOIDP) soft_ver);
    SDsetattr(sdfid, REPLACES, DFNT_CHAR, strlen(replaces) + 1, (VOIDP) replaces);

    if (meta_l2->mission != 0x0) {
        SDsetattr(sdfid, L2BRS_MISSION, DFNT_CHAR, strlen(meta_l2->mission) + 1,
                (VOIDP) meta_l2->mission);
    }

    if (meta_l2->mission_char != 0x0) {
        SDsetattr(sdfid, MSNCHAR, DFNT_CHAR, strlen(meta_l2->mission_char) + 1,
                (VOIDP) meta_l2->mission_char);
    }

    if (meta_l2->sensor != 0x0) {
        SDsetattr(sdfid, SENSOR, DFNT_CHAR, strlen(meta_l2->sensor) + 1,
                (VOIDP) meta_l2->sensor);
    }

    if (meta_l2->sensor_char != 0x0) {
        SDsetattr(sdfid, SNSCHAR, DFNT_CHAR, strlen(meta_l2->sensor_char) + 1,
                (VOIDP) meta_l2->sensor_char);
    }

    SDsetattr(sdfid, PTIME, DFNT_CHAR, strlen(ptime) + 1, (VOIDP) ptime);
    SDsetattr(sdfid, L2BRS_INFILES, DFNT_CHAR,
            strlen(infiles) + 1, (VOIDP) infiles);
    SDsetattr(sdfid, PROC_CON, DFNT_CHAR, strlen(proc_con) + 1, (VOIDP) proc_con);

    if(meta_l2->infiles) {
        pinfile = strrchr(meta_l2->infiles, '/');
        if (pinfile == NULL)
            pinfile = meta_l2->infiles;
        else
            pinfile++;
    } else {
        pinfile = "unknown";
    }
    SDsetattr(sdfid, PINFILES, DFNT_CHAR, strlen(pinfile) + 1,
            (VOIDP) pinfile);

    SDsetattr(sdfid, DTYPE, DFNT_CHAR, strlen(dtype) + 1, (VOIDP) dtype);
    SDsetattr(sdfid, NSAMP, DFNT_INT32, 1, (VOIDP) & nsamp);

    SDsetattr(sdfid, NREC, DFNT_INT32, 1, (VOIDP) & nrec);
    SDsetattr(sdfid, SNCNTR, DFNT_INT32, 1, (VOIDP) & meta_l2->ncrec);
    SDsetattr(sdfid, L2BRS_PCTFLAG, DFNT_FLOAT32, 32, (VOIDP) meta_l2->flags_pc);

    quot1 = div(smsec, MSECHOUR);
    quot2 = div(quot1.rem, MSECMIN);
    quot3 = div(quot2.rem, MSECSEC);
    sprintf(string, "%4.4d%3.3d%2.2d%2.2d%2.2d%3.3d", syear, sday,
            quot1.quot, quot2.quot, quot3.quot, quot3.rem);
    SDsetattr(sdfid, STIME, DFNT_CHAR, strlen(string) + 1, (VOIDP) string);
    quot1 = div(emsec, MSECHOUR);
    quot2 = div(quot1.rem, MSECMIN);
    quot3 = div(quot2.rem, MSECSEC);
    sprintf(string, "%4.4d%3.3d%2.2d%2.2d%2.2d%3.3d", eyear, eday,
            quot1.quot, quot2.quot, quot3.quot, quot3.rem);
    SDsetattr(sdfid, END_TIME, DFNT_CHAR,
            strlen(string) + 1, (VOIDP) string);

    if (meta_l2->ctime != 0x0)
        SDsetattr(sdfid, CTIME, DFNT_CHAR, strlen(meta_l2->ctime) + 1,
            (VOIDP) meta_l2->ctime);
    if (meta_l2->ntime != 0x0)
        SDsetattr(sdfid, NTIME, DFNT_CHAR, strlen(meta_l2->ntime) + 1,
            (VOIDP) meta_l2->ntime);

    SDsetattr(sdfid, SYEAR, DFNT_INT16, 1, (VOIDP) & syear);
    SDsetattr(sdfid, SDAY, DFNT_INT16, 1, (VOIDP) & sday);
    SDsetattr(sdfid, SMSEC, DFNT_INT32, 1, (VOIDP) & smsec);
    SDsetattr(sdfid, EYEAR, DFNT_INT16, 1, (VOIDP) & eyear);
    SDsetattr(sdfid, EDAY, DFNT_INT16, 1, (VOIDP) & eday);
    SDsetattr(sdfid, EMSEC, DFNT_INT32, 1, (VOIDP) & emsec);

    if (meta_l2->snode != 0x0)
        SDsetattr(sdfid, SNODE, DFNT_CHAR, strlen(meta_l2->snode) + 1,
            (VOIDP) meta_l2->snode);
    if (meta_l2->enode != 0x0)
        SDsetattr(sdfid, ENODE, DFNT_CHAR, strlen(meta_l2->enode) + 1,
            (VOIDP) meta_l2->enode);

    SDsetattr(sdfid, ORBNUM, DFNT_INT32, 1, (VOIDP) & meta_l2->orbnum);
    SDsetattr(sdfid, LATUNITS, DFNT_CHAR, strlen(LATUNITS_VAL) + 1,
            (VOIDP) LATUNITS_VAL);
    SDsetattr(sdfid, LONUNITS, DFNT_CHAR, strlen(LONUNITS_VAL) + 1,
            (VOIDP) LONUNITS_VAL);
    SDsetattr(sdfid, NLAT, DFNT_FLOAT32, 1, (VOIDP) & meta_l2->northlat);
    SDsetattr(sdfid, SLAT, DFNT_FLOAT32, 1, (VOIDP) & meta_l2->southlat);
    SDsetattr(sdfid, WLON, DFNT_FLOAT32, 1, (VOIDP) & meta_l2->westlon);
    SDsetattr(sdfid, ELON, DFNT_FLOAT32, 1, (VOIDP) & meta_l2->eastlon);
    SDsetattr(sdfid, STCLAT, DFNT_FLOAT32, 1, (VOIDP) & meta_l2->startclat);
    SDsetattr(sdfid, STCLON, DFNT_FLOAT32, 1, (VOIDP) & meta_l2->startclon);
    SDsetattr(sdfid, ENDCLAT, DFNT_FLOAT32, 1, (VOIDP) & meta_l2->endclat);
    SDsetattr(sdfid, ENDCLON, DFNT_FLOAT32, 1, (VOIDP) & meta_l2->endclon);

    SDsetattr(sdfid, NODEL, DFNT_FLOAT32, 1, (VOIDP) & meta_l2->nodel);

    SDsetattr(sdfid, PARAM, DFNT_CHAR, strlen(param_val) + 1, (VOIDP) param_val);
    SDsetattr(sdfid, UNITS, DFNT_CHAR, strlen(units_val) + 1, (VOIDP) units_val);
    SDsetattr(sdfid, PX_START, DFNT_INT32, 1, (VOIDP) & px_start);
    SDsetattr(sdfid, PX_END, DFNT_INT32, 1, (VOIDP) & px_end);
    lac_px_st = meta_l2->pix_start + ((px_start - 1) * meta_l2->pix_sub);
    SDsetattr(sdfid, LAC_PX_ST, DFNT_INT32, 1, (VOIDP) & lac_px_st);
    SDsetattr(sdfid, PX_SUBSAMP, DFNT_INT32, 1, (VOIDP) & px_subsamp);
    lac_px_subsamp = meta_l2->pix_sub * px_subsamp;
    SDsetattr(sdfid, LAC_PX_SUBSAMP, DFNT_INT32, 1, (VOIDP) & lac_px_subsamp);
    SDsetattr(sdfid, PX_NUM, DFNT_INT32, 1, (VOIDP) & brs_nsamp);
    SDsetattr(sdfid, SC_START, DFNT_INT32, 1, (VOIDP) & sc_start);
    SDsetattr(sdfid, SC_END, DFNT_INT32, 1, (VOIDP) & sc_end);
    SDsetattr(sdfid, SC_SUBSAMP, DFNT_INT32, 1, (VOIDP) & sc_subsamp);
    SDsetattr(sdfid, SC_NUM, DFNT_INT32, 1, (VOIDP) & brs_nrec);
    SDsetattr(sdfid, PX_LL_NUM, DFNT_INT32, 1, (VOIDP) & brs_nsamp);
    SDsetattr(sdfid, SC_LL_NUM, DFNT_INT32, 1, (VOIDP) & brs_nrec);
    SDsetattr(sdfid, SC_TYPE, DFNT_CHAR, strlen(sc_type) + 1, (VOIDP) sc_type);
    SDsetattr(sdfid, SC_EQN, DFNT_CHAR, strlen(sc_eqn) + 1, (VOIDP) sc_eqn);
    SDsetattr(sdfid, BASE, DFNT_FLOAT32, 1, (VOIDP) & base);
    SDsetattr(sdfid, SLOPE, DFNT_FLOAT32, 1, (VOIDP) & slope);
    SDsetattr(sdfid, INTERCEPT, DFNT_FLOAT32, 1, (VOIDP) & intercept);

}

/*-----------------------------------------------------------------------------
    Function: write_image

    Returns: int32 (status)
        The return code is FAIL (-1) if an error occurs, SUCCEED (0)
        otherwise.

    Description:
        The function write_image writes the image in to the product file

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      char *    l2brs_path  I     directory path & file name for browse product
      char      l2brs_data  I     image data
      int32     px_num      I     no. of pixels in the browse image
      int32     sc_num      I     no. of scan lines in the browse image
      uchar *   palette     I     RGB wts for each of 256 gray-levels of the
                                  l2brs_data byte values
    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/23/94    Original development
        Lakshmi Kumar    Hughes STX      06/16/94    Made changes to agree with
                                                     Interface specs v3.1
-----------------------------------------------------------------------------*/
int32
write_image(char *l2brs_path, unsigned char *l2brs_data, int32 brs_nsamp, int32 brs_nrec,
        unsigned char *palette) {
    int32 img_ref;

    DFR8setpalette(palette);

    if ((DFR8addimage(l2brs_path, (VOIDP) l2brs_data, brs_nsamp, brs_nrec, 0)) < 0) {
        fprintf(stderr, "\n Error writing browse image\n");
        return FAIL;
    }

    if ((img_ref = DFR8lastref()) > 0) {
        if ((DFANputlabel(l2brs_path, DFTAG_RIG, img_ref, "brs_data")) < 0) {
            fprintf(stderr, "\nwrite_image: Error writing - brs_data label\n");
            fprintf(stderr, "\n No label is written to the raster image\n");
        }
    } else {
        fprintf(stderr, "\nwrite_image: Error reading image ref. no. \n");
        fprintf(stderr, "\n No label is written to the raster image\n");
    }

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function: write_SDS

    Returns: int32 (status)
        The return code is FAIL (-1) if an error occurs, SUCCEED (0)
        otherwise.

    Description:
        The function write_image writes the image in to the product file

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     sdfid       I     file ID
      char *    label       I     Name of the dataset
      int32     ntype       I     Data type
      int32     rank        I     Number of dimensions of the dataset
      int32 *   dimsizes    I     Dimensions of the data
      int32 *   start       I     start dimensions of the data
      void  *   buf         I     data buffer

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/23/94    Original development
        Lakshmi Kumar    Hughes STX      06/16/94    Made changes to agree with
                                                     Interface specs v3.1
-----------------------------------------------------------------------------*/
int32
write_SDS(int32 sdfid, char *label, int32 ntype, int32 rank, int32 *dimsizes,
        int32 *start, void *buf) {
    int32 sdsid, ret;

    if ((sdsid = SDcreate(sdfid, label, ntype, rank, dimsizes)) < 0)
        return FAIL;

    if ((ret = SDwritedata(sdsid, start, NULL, dimsizes, (VOIDP) buf)) < 0)
        return FAIL;

    SDendaccess(sdsid);
    return sdsid;
}

/*-----------------------------------------------------------------------------
    Function: write_nav_sets

    Returns: int32 (status)
        The return code is FAIL (-1) if an error occurs, SUCCEED (0)
        otherwise.

    Description:
        The function write_nav_sets creates navigation vgroup and its
        datasets.

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     fid         I     file ID
      int32     sdfid       I     SD interface file ID
      float32   *orb_vec    I     Orbit position vector data buffer
      float32   *l_vert     I     local vertical vector in ECEF frame databuf
      float32   *sun_ref    I     referecnce sun vector data buffer
      float32   *att_ang    I     computed yaw, roll, pitch data buffer
      float32   *sen_mat    I     sensor-frame matrix data buffer
      float32   *scan_ell   I     scan-track ellipse coefficients data buffer
      int32     *nflag      I     navigation flags buffer

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      05/12/95    Original development
-----------------------------------------------------------------------------*/
int32
write_nav_sets(int32 fid, int32 sdfid, int32 brs_nrec, int32 brs_nsamp,
        int16 *cntl_pt_lat, int16 *cntl_pt_lon) {

    int32 vid, sdsid, sdsref;
    int32 rank, start[3] = {0, 0, 0}, dimsizes[3];
    char FUNC[] = "write_nav_sets:";
    float32 f32;

    if ((vid = Vattach(fid, -1, "w")) < 0) {
        sprintf(ERR_MSG, "\n%s Vattach failed for navigation vgp", FUNC);
        return FAIL;
    }

    Vsetname(vid, "Navigation");

    /****   write cntl_pt_lat SDS and link it to navigation vgroup */
    dimsizes[0] = brs_nrec;
    dimsizes[1] = brs_nsamp;
    dimsizes[2] = 0;
    if ((sdsid = write_SDS(sdfid, CNTLPTLAT, DFNT_INT16, rank = 2,
            dimsizes, start, cntl_pt_lat)) < 0) return FAIL;

    if ((SDsetattr(sdsid, LONGNAME, DFNT_CHAR, strlen(CNTLPTLAT_NAME) + 1,
            (VOIDP) CNTLPTLAT_NAME)) < 0) {
        sprintf(ERR_MSG,
                "\n%s Error writing attribute %s for cnlt_pt_lat", FUNC, LONGNAME);
        return FAIL;
    }

    if ((sdsref = SDidtoref(sdsid)) < 0) {
        sprintf(ERR_MSG, "\nError SDidtoref failed for cntl_pt_lat dataset");
        return FAIL;
    }

    if ((Vaddtagref(vid, DFTAG_NDG, sdsref)) < 0) {
        sprintf(ERR_MSG, "\nError Vaddtagref failed for cntl_pt_lat");
        return FAIL;
    }

    SDsetdimname(SDgetdimid(sdsid, 0), "Number of Scans");
    SDsetdimname(SDgetdimid(sdsid, 1), "Pixels per Scan");

    f32 = 1. / 360;
    SDsetattr(sdsid, "slope", DFNT_FLOAT32, 1, &f32);
    f32 = 0.0;
    SDsetattr(sdsid, "intercept", DFNT_FLOAT32, 1, &f32);


    /****   write cntl_pt_lon SDS and link it to navigation vgroup */
    dimsizes[0] = brs_nrec;
    dimsizes[1] = brs_nsamp;
    dimsizes[2] = 0;
    if ((sdsid = write_SDS(sdfid, CNTLPTLON, DFNT_INT16, rank = 2,
            dimsizes, start, cntl_pt_lon)) < 0) return FAIL;

    if ((SDsetattr(sdsid, LONGNAME, DFNT_CHAR, strlen(CNTLPTLON_NAME) + 1,
            (VOIDP) CNTLPTLON_NAME)) < 0) {
        sprintf(ERR_MSG,
                "\n%s Error writing attribute %s for cnlt_pt_lon", FUNC, LONGNAME);
        return FAIL;
    }

    if ((sdsref = SDidtoref(sdsid)) < 0) {
        sprintf(ERR_MSG, "\nError SDidtoref failed for cntl_pt_lon dataset");
        return FAIL;
    }

    if ((Vaddtagref(vid, DFTAG_NDG, sdsref)) < 0) {
        sprintf(ERR_MSG, "\nError Vaddtagref failed for cntl_pt_lon");
        return FAIL;
    }

    SDsetdimname(SDgetdimid(sdsid, 0), "Number of Scans");
    SDsetdimname(SDgetdimid(sdsid, 1), "Pixels per Scan");

    f32 = 1.0 / 180.0;
    SDsetattr(sdsid, "slope", DFNT_FLOAT32, 1, &f32);
    f32 = 0.0;
    SDsetattr(sdsid, "intercept", DFNT_FLOAT32, 1, &f32);

    Vdetach(vid);

    return SUCCEED;

}

/*-----------------------------------------------------------------------------
    Function: write_tilt_sets

    Returns: int32 (status)
        The return code is FAIL (-1) if an error occurs, SUCCEED (0)
        otherwise.

    Description:
        The function write_tilt_sets creates tilt vgroup and its
        datasets.

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     fid         I     file ID
      int32     sdfid       I     SD interface file ID
      int32     ntilts	    I     Number of tilts
      short	tilt_flags  I     tilt flags corresponding to each tilt state
      short     tilt_ranges I     first and last scan line nos. corresponding
                                        to each tilt state
      float32   tilt_lats   I     latitudes of the end pixels for the scan
                                        lines of tilt_ranges
      float32   tilt_lons   I     longitudes of the end pixels for the scan
                                        lines of tilt_ranges

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      05/12/95    Original development
-----------------------------------------------------------------------------*/
int32
write_tilt_sets(int32 fid, int32 sdfid, int32 ntilts, short *tilt_flags,
        int16 *tilt_ranges) {
    int32 vid, sdsid, sdsref;
    int32 rank, start[3] = {0, 0, 0}, dimsizes[3] = {0, 0, 0};
    int16 int_range[2];
    char FUNC[] = "write_tilt_sets";


    if ((vid = Vattach(fid, -1, "w")) < 0) {
        sprintf(ERR_MSG, "\n%s Vattach failed for tilt vgp", FUNC);
        return FAIL;
    }

    Vsetname(vid, "Sensor Tilt");

    /*   write ntilts SDS */

    dimsizes[0] = 1;
    if ((sdsid = write_SDS(sdfid, NTILTS, DFNT_INT32, rank = 1,
            dimsizes, start, &ntilts)) < 0)
        return FAIL;
    SDsetattr(sdsid, LONGNAME, DFNT_CHAR, strlen(NTILTS_NAME) + 1,
            (VOIDP) NTILTS_NAME);

    if ((sdsref = SDidtoref(sdsid)) < 0) {
        sprintf(ERR_MSG, "\nError SDidtoref failed for ntilts dataset");
        return FAIL;
    }

    if ((Vaddtagref(vid, DFTAG_NDG, sdsref)) < 0) {
        sprintf(ERR_MSG, "\nError Vaddtagref failed for ntilts");
        return FAIL;
    }

    /*   write tilt_flags SDS */

    dimsizes[0] = MAXTILTS;
    if ((sdsid = write_SDS(sdfid, T_FLAGS, DFNT_INT16, rank = 1, dimsizes,
            start, tilt_flags)) < 0)
        return FAIL;
    SDsetattr(sdsid, LONGNAME, DFNT_CHAR, strlen(T_FLAGS_NAME) + 1,
            (VOIDP) T_FLAGS_NAME);

    int_range[0] = -1;
    int_range[1] = 3;
    SDsetattr(sdsid, RANGE, DFNT_INT16, 2, (VOIDP) int_range);

    if ((sdsref = SDidtoref(sdsid)) < 0) {
        sprintf(ERR_MSG, "\nError SDidtoref failed for tilt_flags dataset");
        return FAIL;
    }

    if ((Vaddtagref(vid, DFTAG_NDG, sdsref)) < 0) {
        sprintf(ERR_MSG, "\nError Vaddtagref failed for tilt_flags");
        return FAIL;
    }


    /*   write tilt_ranges SDS */

    dimsizes[0] = MAXTILTS;
    dimsizes[1] = 2;
    if ((sdsid = write_SDS(sdfid, T_RANGES, DFNT_INT16, rank = 2,
            dimsizes, start, tilt_ranges)) < 0)
        return FAIL;
    SDsetattr(sdsid, LONGNAME, DFNT_CHAR, strlen(T_RANGES_NAME) + 1,
            (VOIDP) T_RANGES_NAME);

    if ((sdsref = SDidtoref(sdsid)) < 0) {
        sprintf(ERR_MSG, "\nError SDidtoref failed for tilt_ranges dataset");
        return FAIL;
    }

    if ((Vaddtagref(vid, DFTAG_NDG, sdsref)) < 0) {
        sprintf(ERR_MSG, "\nError Vaddtagref failed for tilt_ranges");
        return FAIL;
    }

    Vdetach(vid);

    return SUCCEED;

}


