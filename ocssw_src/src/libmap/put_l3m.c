#include <ctype.h>
#include "map.h"
#include "mapattr.h"
#include "meta_l3b.h"
#include "mapproto_o.h"
#include "maplists.h"
#include <palette.h>

#include <math.h>

#define VERSION "1.0"

#define max(A, B) ((A) > (B) ?  (A) : (B))
#define min(A, B) ((A) > (B) ?  (B) : (A))

/*-----------------------------------------------------------------------------

    Function: put_l3m

    Returns: intn (status)
        The return code is 
                FAIL (-1) 	 If an I/O error occurs.
                NOMATCH_ERR (-2) The routine checks l3m_name to find out for
                                 which parameter the data is being written. 
                                 It tries to match a substring of the paramter
                                 name with given l3m_name.  If no match occurs
                                 it returns NOMATCH_ERR (-2).
                SUCCEED (0)	 If no errors occur.

    Description:
        The function put_l3m.c creates a level 3 Standard Map Image (SMI)
        file with the file name given by l3m_path.  Level 3 SMI file meta
        data, palette, and data will be written.  See "INTERFACE SPECIFICATIONS
        FOR SeaWiFS OPERATIONAL PRODUCT INPUT AND OUTPUT AND SENSOR CALIBRATION
        SOFTWARE", v3.3, dated 12th July 94, (By Fred Patt, Jim Firestone, 
        and Mike Darzi) and "SeaWiFS OPERATIONAL ARCHIVE PRODUCT 
        SPECIFICATIONS, v1.0 (Incomplete DRAFT version), dated 22 July 1994 
        (By F. Patt, J. Firestone, B. Schieber, L. Kumar, D. Ilg, and M. Darzi)
        for details.

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      char *    l3m_path    I     directory path & file name for SMI product
      char *    replaces    I     filename of previously generated product 
                                  that current product is intended to replace
      int16 	bin_syear   I     year of start of binning period of the parent
                                  product 
      int16     bin_sday    I     GMT day-of-year of start of binning period
                                  of the parent product
      int16     bin_eyear   I     year of end of binning period of the parent
                                  product
      int16	bin_eday    I     GMT day-of-year of end of binning period
                                  of the parent product
      int16     syear       I     data start time year (from l3 bin product)
      int16     sday        I     data start time day-of-year
      int32     smsec       I     data start time milliseconds-of-day
      int16     eyear       I     data end time year; read from l3 bin product
      int16     eday        I     data end day-of-year
      int16     emsec       I     data end time milliseconds-of-day
      float32  *lat_range   I     lats of the outside edges of the northernmost
                                  and easternmost columns for l3mdata
      float32  *lon_range   I     lons of the outside edges of the westernmost
                                  and easternmost columns for l3m_data
      int32     lines	    I     number of points in the vertical direction
      int32     columns	    I     number of points in horizontal direction
      char     *flag_use    I     specifies use of parent l2 product's l2_flags
      uint8    *eng_q_use   I     specifies use of l2 product's eng_qua values
      char *    ptime       I     processing start time
      char *    infiles     I     name of input product
      char *    prod_type   I     binning period description(scene, day.....)
      int32     nbins       I     Number of bins containing data in the parent
      char *    l3m_name    I     name of the geophysical parameter to which
                                  l3m_data corresponds
      void *    l3m_data    I     image data for parameter l3m_name
                                  product
      char  *   measure	    I     indicates whether it is mean, median, mode,
                                  standard deviation, scenes or pixels   
      char *    proc_con    I     processing control input
      char *    proc_log    I     processing log output

    Notes: (see V4.2 I/O & V2.6 Product Specs for updated changes) 

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      12/16/93    Original development
        Lakshmi Kumar    Hughes STX      08/29/94    Revised version     
        Lakshmi Kumar    Hughes STX      03/01/95    changed l3m data type to
                                                     float & added code to
                                                     scale the data to byte
                                                     Added some global attrs
                                                     Changed l3m_name "epsilon"
                                                     to "eps_68"
        Lakshmi Kumar    Hughes STX      03/14/95    Added label "palette" to
                                                     palette.
        Lakshmi Kumar	 Hughes STX	 05/08/95    Removed longname attribute
                                                     from l3m_data SDS
        Lakshmi Kumar	 HITC		 05/18/95    Changed datatype of
                                                     flag_use global attribute	
        Lakshmi Kumar	 Hughes STX	 05/27/95    Removed palette fr put_l3m
                                                     interface call & is 
                                                     defined in an include file 
        Lakshmi Kumar	 Hughes STX	 09/28/95    Added orbit & flag_names
                                                     as input args and written
                                                     them as global attributes
                                                     Added SW point lat/lon as
                                                     global attributes 
                                                     (V2.8 prod & V4.4 I/O 
                                                        specs)
        Lakshmi Kumar	 Hughes STX	 10/27/95    Added "Start Orbit" and 
                                                     "End Orbit" global attrs.
                                                     moved the prototype defn.
                                                     of this fun. to an .h file
        Lakshmi Kumar	 Hughes STX	 11/08/95    Changed eps_68 to eps_78
        Lakshmi Kumar	 Hughes STX	 02/23/96    Added rounding factor to
                                                     the scaling equation
                                                     changed local variable
                                                     names from min & max to 
                                                     datamin and datamax.
        Lakshmi Kumar    Hughes STX      06/18/96    Changed defn. of MAX to
                                                     MAXVAL inorder to remove
                                                     compile time warning
        Lakshmi Kumar	 Hughes STX 	 10/16/96    Removed 'median' & 'mode'
                                                     as input for measure.
                                                     Sets NODATA pixels to 255.
                                                     Ref. V5.0 I/O specs.
        Lakshmi Kumar	 Hughes STX	 12/11/96    Removed 'orbit' as input
                                                     (it can be accessed thro'
                                                      meta_l3b structure)
                                                     Also, removed non-ANSI
                                                     declarations
-----------------------------------------------------------------------------*/
int32
put_l3m(char *l3m_path, char *replaces, int16 bin_syear, int16 bin_sday,
        int16 bin_eyear, int16 bin_eday, int16 syear, int16 sday,
        int32 smsec, int16 eyear, int16 eday, int32 emsec,
        float32 *lat_range, float32 *lon_range, int32 lines, int32 columns,
        char *flag_names, char *flag_use, uint8 *eng_q_use, char *ptime,
        char *infiles, char *prod_type, int32 nbins, char *l3m_name,
        void *l3m_data, char *measure, char *proc_con, char *proc_log,
        meta_l3bType *meta_l3b) {

    int32 fid, sdfid, sdsid;
    int32 nrows = lines, ncols = columns;
    int32 i, rank = 2, index, pal_ref;
    int32 start[2] = {0, 0};
    int32 dimsizes[2];
    uint8 *l3m_outdata = NULL; /* buf to hold scaled data */
    const uint8 NODATA_VAL = 255;
    div_t quot1, quot2, quot3;
    char parm_name[MAXVAL];
    char string[MAXVAL], *loc_path, *str, loc_measure[25];
    char *softid = VERSION; /* SW id obtained fr the makefile*/
    float32 max_val, min_val, datamax, datamin;
    float32 lat_step, lon_step, slope, intercept;
    float32 swlat, swlon;
    const float32 NODATA_FLAG = -9E6f;
    VOIDP pal;
    float *l3m_fdata = (float *) l3m_data;
    uint8 *l3m_bdata = (uint8 *) l3m_data;

    uint8 sixteenbit = 0;

    /*** Open/create output file in the given l3m_path */

    sdfid = SDstart(l3m_path, DFACC_CREATE);
    fid = Hopen(l3m_path, DFACC_RDWR, 0);

    if (sdfid < 0 || fid < 0) {
        printf("\n put_l3m: Error in creating file %s ", l3m_path);
        return FAIL;
    }

    /*** Convert l3m_name to lower case */

    for (i = 0; i < strlen(l3m_name); i++)
        parm_name[i] = tolower(l3m_name[i]);

    parm_name[i] = 0;

    /*** Use parameter name as index in to parameter full name list, units list,
            scales list and intercept list  */

    printf("name=%s\n", parm_name);

    if ((strcmp(parm_name, "nlw_412")) == 0)
        index = NLW_412;
    else if ((strcmp(parm_name, "nlw_443")) == 0)
        index = NLW_443;
    else if ((strcmp(parm_name, "nlw_490")) == 0)
        index = NLW_490;
    else if ((strcmp(parm_name, "nlw_510")) == 0)
        index = NLW_510;
    else if ((strcmp(parm_name, "nlw_555")) == 0)
        index = NLW_555;
    else if ((strcmp(parm_name, "nlw_670")) == 0)
        index = NLW_670;
    else if ((strcmp(parm_name, "czcs_pigment")) == 0)
        index = CZCS_PIGMENT;
    else if ((strcmp(parm_name, "chlor_a")) == 0)
        index = CHLOR_A;
    else if ((strcmp(parm_name, "k_490")) == 0)
        index = K_490;
    else if ((strcmp(parm_name, "chlor_a_k_490")) == 0)
        index = CHOR_A_K_490;
    else if ((strcmp(parm_name, "eps_78")) == 0)
        index = EPSILON;
    else if ((strcmp(parm_name, "tau_865")) == 0)
        index = TAU_865;
    else if ((strcmp(parm_name, "angstrom_510")) == 0)
        index = ANGSTROM_510;
    else if ((strcmp(parm_name, "pixels")) == 0)
        index = PIXELS;
    else if ((strcmp(parm_name, "scenes")) == 0)
        index = SCENES;
    else if ((strcmp(parm_name, "ndvi")) == 0)
        index = NDVI;
    else if ((strcmp(parm_name, "biosphere")) == 0)
        index = BIOSPHERE;
    else {
        printf("\nput_l3m: Error: given l3m_name - %s unrecognized\n", l3m_name);
        free(l3m_outdata);
        return FAIL;
    }

    if (nbins < 0) {
        nbins = -nbins;
        sixteenbit = 1;
    }


    /***  Write out global attributes  */

    if ((loc_path = strrchr(l3m_path, '/')) != NULL)
        loc_path++;
    else
        loc_path = l3m_path;
    SDsetattr(sdfid, L3M_PNAME, DFNT_CHAR, strlen(loc_path) + 1, (VOIDP) loc_path);
    SDsetattr(sdfid, L3M_TITLE, DFNT_CHAR, strlen(L3M_TITLE_VAL) + 1, L3M_TITLE_VAL);
    SDsetattr(sdfid, L3M_DCENTER, DFNT_CHAR,
            strlen(L3M_DCENTER_VAL) + 1, L3M_DCENTER_VAL);
    SDsetattr(sdfid, L3M_MISSION, DFNT_CHAR,
            strlen(L3M_MISSION_VAL) + 1, (VOIDP) L3M_MISSION_VAL);
    SDsetattr(sdfid, L3M_SENSOR_NAME, DFNT_CHAR,
            strlen(meta_l3b->sensor_name) + 1, (VOIDP) meta_l3b->sensor_name);
    SDsetattr(sdfid, L3M_PRODTYPE, DFNT_CHAR,
            strlen(prod_type) + 1, (VOIDP) prod_type);
    if (replaces == NULL || (strcmp(replaces, "")) == 0) {
        strcpy(string, "ORIGINAL");
        replaces = string;
    }
    SDsetattr(sdfid, L3M_REPLACE, DFNT_CHAR,
            strlen(replaces) + 1, (VOIDP) replaces);
    SDsetattr(sdfid, L3M_SOFTID, DFNT_CHAR, strlen(softid) + 1, (VOIDP) softid);
    SDsetattr(sdfid, L3M_PTIME, DFNT_CHAR, strlen(ptime) + 1, (VOIDP) ptime);

    if ((str = strrchr(infiles, '/')) != NULL)
        str++;
    else
        str = infiles;
    SDsetattr(sdfid, L3M_INFILES, DFNT_CHAR, strlen(str) + 1, (VOIDP) str);
    SDsetattr(sdfid, L3M_PROCCON, DFNT_CHAR, strlen(proc_con) + 1, (VOIDP) proc_con);
    /*
    SDsetattr(sdfid, L3M_PROCLOG, DFNT_CHAR, strlen(proc_log)+1,(VOIDP)proc_log);
     */
    SDsetattr(sdfid, L3M_FLAG_NAMES, DFNT_CHAR,
            strlen(flag_names) + 1, (VOIDP) flag_names);
    quot1 = div(smsec, MSECHOUR);
    quot2 = div(quot1.rem, MSECMIN);
    quot3 = div(quot2.rem, MSECSEC);
    sprintf(string, "%4.4d%3.3d%2.2d%2.2d%2.2d%3.3d", syear, sday,
            quot1.quot, quot2.quot, quot3.quot, quot3.rem);
    SDsetattr(sdfid, L3M_STIME, DFNT_CHAR, strlen(string) + 1, (VOIDP) string);
    quot1 = div(emsec, MSECHOUR);
    quot2 = div(quot1.rem, MSECMIN);
    quot3 = div(quot2.rem, MSECSEC);
    sprintf(string, "%4.4d%3.3d%2.2d%2.2d%2.2d%3.3d", eyear, eday,
            quot1.quot, quot2.quot, quot3.quot, quot3.rem);
    SDsetattr(sdfid, L3M_ETIME, DFNT_CHAR, strlen(string) + 1, (VOIDP) string);
    SDsetattr(sdfid, L3M_ORBIT, DFNT_INT32, 1, (VOIDP) & meta_l3b->orbit);
    SDsetattr(sdfid, L3M_SORBIT, DFNT_INT32, 1, (VOIDP) & meta_l3b->start_orb);
    SDsetattr(sdfid, L3M_EORBIT, DFNT_INT32, 1, (VOIDP) & meta_l3b->end_orb);
    SDsetattr(sdfid, L3M_MAPPROJ, DFNT_CHAR, strlen(L3M_MAPPROJ_VAL) + 1,
            (VOIDP) L3M_MAPPROJ_VAL);
    SDsetattr(sdfid, L3M_LATUNITS, DFNT_CHAR, strlen(L3M_LATUNITS_VAL) + 1,
            (VOIDP) L3M_LATUNITS_VAL);
    SDsetattr(sdfid, L3M_LONUNITS, DFNT_CHAR, strlen(L3M_LONUNITS_VAL) + 1,
            (VOIDP) L3M_LONUNITS_VAL);

    /*** check validity of given latitude ranges-- 1st value should be greater
            than the second one.  If in error, output error message  	*/
    if (lat_range[0] <= lat_range[1]) {
        printf("\nput_l3m: Given latitude range is in error \n");
        free(l3m_outdata);
        return FAIL;
    }
    SDsetattr(sdfid, L3M_NLAT, DFNT_FLOAT32, 1, (VOIDP) & lat_range[0]);
    SDsetattr(sdfid, L3M_SLAT, DFNT_FLOAT32, 1, (VOIDP) & lat_range[1]);
    SDsetattr(sdfid, L3M_WLON, DFNT_FLOAT32, 1, (VOIDP) & lon_range[0]);
    SDsetattr(sdfid, L3M_ELON, DFNT_FLOAT32, 1, (VOIDP) & lon_range[1]);

    lat_step = (lat_range[0] - lat_range[1]) / lines;
    SDsetattr(sdfid, L3M_LAT_STEP, DFNT_FLOAT32, 1, (VOIDP) & lat_step);

    if (lon_range[1] > lon_range[0])
        lon_step = (lon_range[1] - lon_range[0]) / columns;
    else
        lon_step = (360 - lon_range[0] + lon_range[1]) / columns;
    SDsetattr(sdfid, L3M_LON_STEP, DFNT_FLOAT32, 1, (VOIDP) & lon_step);

    swlat = (lat_range[1] + lat_step / 2.0);
    swlon = (lon_range[0] + lon_step / 2.0);
    SDsetattr(sdfid, L3M_SWLAT, DFNT_FLOAT32, 1, (VOIDP) & swlat);
    SDsetattr(sdfid, L3M_SWLON, DFNT_FLOAT32, 1, (VOIDP) & swlon);

    SDsetattr(sdfid, L3M_DATABINS, DFNT_INT32, 1, (VOIDP) & nbins);
    SDsetattr(sdfid, L3M_NROWS, DFNT_INT32, 1, (VOIDP) & nrows);
    SDsetattr(sdfid, L3M_NCOLS, DFNT_INT32, 1, (VOIDP) & ncols);

    SDsetattr(sdfid, L3M_PARAMETER, DFNT_CHAR, strlen(parmname_list[index]) + 1,
            (VOIDP) parmname_list[index]);

    for (i = 0; i < strlen(measure); i++)
        loc_measure[i] = tolower(measure[i]);
    loc_measure[i] = 0;

    if (strcmp(loc_measure, "mean") == 0)
        strcpy(loc_measure, "Mean");
    else
        if (strcmp(loc_measure, "standard deviation") == 0)
        strcpy(loc_measure, "Standard deviation");
    else
        if (strcmp(loc_measure, "maximum") == 0)
        strcpy(loc_measure, "Max Value");
    else
        if (strcmp(loc_measure, "scenes") == 0)
        strcpy(loc_measure, "Scenes per bin");
    else
        if (strcmp(loc_measure, "pixels") == 0)
        strcpy(loc_measure, "Pixels per bin");
    else {
        printf("\nput_l3m: measure does not contain valid string\n");
        return FAIL;
    }


    SDsetattr(sdfid, L3M_MEASURE, DFNT_CHAR, strlen(loc_measure) + 1,
            (VOIDP) loc_measure);

    SDsetattr(sdfid, L3M_UNITS, DFNT_CHAR, strlen(unit_list[index]) + 1,
            (VOIDP) unit_list[index]);

    if (index == CHLOR_A || index == CZCS_PIGMENT) {
        SDsetattr(sdfid, L3M_SCALING, DFNT_CHAR, strlen(L3M_LOG_SCALE) + 1,
                (VOIDP) L3M_LOG_SCALE);
        SDsetattr(sdfid, L3M_SC_EQN, DFNT_CHAR, strlen(L3M_LOG_EQN) + 1,
                (VOIDP) L3M_LOG_EQN);
        SDsetattr(sdfid, L3M_BASE, DFNT_FLOAT32, 1, (VOIDP) & base);
    } else {
        SDsetattr(sdfid, L3M_SCALING, DFNT_CHAR, strlen(L3M_LINEAR_SCALE) + 1,
                (VOIDP) L3M_LINEAR_SCALE);
        SDsetattr(sdfid, L3M_SC_EQN, DFNT_CHAR, strlen(L3M_LINEAR_EQN) + 1,
                (VOIDP) L3M_LINEAR_EQN);
    }

    SDsetattr(sdfid, L3M_SLOPE, DFNT_FLOAT32, 1, (VOIDP) & slope_list[index]);
    SDsetattr(sdfid, L3M_INTERCEPT, DFNT_FLOAT32, 1,
            (VOIDP) & intercept_list[index]);

    slope = slope_list[index];
    intercept = intercept_list[index];

    datamin = 9E6f;
    datamax = -9E6f;
    dimsizes[0] = lines;
    dimsizes[1] = columns;


    if (index == BIOSPHERE)
        l3m_outdata = l3m_bdata;

    else if (index == CHLOR_A || index == CZCS_PIGMENT) {
        if ((l3m_outdata = (uint8 *) malloc(lines * columns * sizeof (uint8))) == NULL) {
            printf("\n put_l3m: malloc Error:  Memory allocation unsuccessful\n");
            return FAIL;
        }
        max_val = pow(base, slope * 254 + intercept);
        min_val = pow(base, intercept);
        for (i = 0; i < nrows * ncols; i++) {
            if (l3m_fdata[i] > datamax)
                datamax = l3m_fdata[i];
            if (l3m_fdata[i] < datamin && l3m_fdata[i] > NODATA_FLAG)
                datamin = l3m_fdata[i];
            /* add 0.5 at the end to accomplish rounding */
            if (l3m_fdata[i] < NODATA_FLAG)
                l3m_outdata[i] = NODATA_VAL; /* NODATA_VAL = 255 */
            else
                l3m_outdata[i] =
                    ((log10(min((max(l3m_fdata[i], min_val)), max_val))
                    - intercept) / slope) + 0.5;
        }

    } else {
        if ((l3m_outdata = (uint8 *) malloc(lines * columns * sizeof (uint8))) == NULL) {
            printf("\n put_l3m: malloc Error:  Memory allocation unsuccessful\n");
            return FAIL;
        }
        min_val = intercept;
        max_val = slope * 254 + intercept;
        for (i = 0; i < nrows * ncols; i++) {
            if (l3m_fdata[i] > datamax)
                datamax = l3m_fdata[i];
            if (l3m_fdata[i] < datamin && l3m_fdata[i] > NODATA_FLAG)
                datamin = l3m_fdata[i];
            if (l3m_fdata[i] < NODATA_FLAG)
                l3m_outdata[i] = NODATA_VAL;
            else
                l3m_outdata[i] = ((min((max(l3m_fdata[i], min_val)), max_val)
                    - intercept) / slope) + 0.5;
        }
    }

    SDsetattr(sdfid, L3M_MIN, DFNT_FLOAT32, 1, (VOIDP) & datamin);
    SDsetattr(sdfid, L3M_MAX, DFNT_FLOAT32, 1, (VOIDP) & datamax);

    /***  write  SDS */
    if (sixteenbit) {
        if ((sdsid = SDcreate(sdfid, SDS_NAME, DFNT_UINT16, rank, dimsizes)) < 0) {
            printf("\nput_l3m: SDcreate call to create %s is unsuccessful",
                    SDS_NAME);
            free(l3m_outdata);
            return FAIL;
        }

        if ((SDwritedata(sdsid, start, NULL, dimsizes, (VOIDP) l3m_outdata)) < 0) {
            printf("\nput_l3m: SDwritedata unsuccessful\n");
            free(l3m_outdata);
            return FAIL;
        }
    } else {
        if ((sdsid = SDcreate(sdfid, SDS_NAME, DFNT_UINT8, rank, dimsizes)) < 0) {
            printf("\nput_l3m: SDcreate call to create %s is unsuccessful",
                    SDS_NAME);
            free(l3m_outdata);
            return FAIL;
        }

        if ((SDwritedata(sdsid, start, NULL, dimsizes, (VOIDP) l3m_outdata)) < 0) {
            printf("\nput_l3m: SDwritedata unsuccessful\n");
            free(l3m_outdata);
            return FAIL;
        }
    }

    SDendaccess(sdsid);

    /***  write palette  */
    switch (index) {
    case NDVI: pal = (VOIDP) ndvi_palette;
        break;
    case BIOSPHERE: pal = (VOIDP) bios_palette;
        break;
    default: pal = (VOIDP) default_palette;
        break;
    }

    if ((DFPaddpal(l3m_path, pal)) < 0)
        return FAIL;

    if ((pal_ref = DFPlastref()) > 0) {
        if ((DFANputlabel(l3m_path, DFTAG_IP8, pal_ref, "palette")) < 0) {
            printf("\nput_l3m: Error writing - palette label\n");
            printf("\n No label is written to the palette\n");
        }
    }

    /*** free the buffer space that was allocated */
    /*
      free(l3m_outdata);
     */

    /***  close down interfaces and close file  */
    SDend(sdfid);
    Hclose(fid);
    return SUCCEED;
}

