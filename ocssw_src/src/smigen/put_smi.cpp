#include "netcdf.h"  // Needs to be first to define netcdf stuff JMG
#include <stdlib.h>
#include "smi_map.h"
#include "meta_l3b.h"
#include "mapattr.h"
#include "smiinc.h"
#include "smigen_input.h"
#include "hdf5.h"
#include "hdf5utils.h"
#include "productInfo.h"
#include "sensorDefs.h"
#include "sensorInfo.h"
#include <timeutils.h>
#include <nc4utils.h>

#define MAXVAL     255

using namespace Hdf;

intn setAttr(uint8 isHDF5, int32 obj_id, const char *attr_name, int32 data_type,
        int32 count, VOIDP values) {
    if (isHDF5 == 1) {
        hid_t dtype = H5T_NATIVE_CHAR;
        if (data_type == DFNT_CHAR && count > 1) dtype = H5T_STRING;
        else if (data_type == DFNT_CHAR) dtype = H5T_NATIVE_CHAR;
        else if (data_type == DFNT_UINT8) dtype = H5T_NATIVE_UCHAR;
        else if (data_type == DFNT_INT16) dtype = H5T_NATIVE_SHORT;
        else if (data_type == DFNT_UINT16) dtype = H5T_NATIVE_USHORT;
        else if (data_type == DFNT_INT32) dtype = H5T_STD_I32LE;
        else if (data_type == DFNT_UINT32) dtype = H5T_STD_U32LE;
        else if (data_type == DFNT_FLOAT32) dtype = H5T_NATIVE_FLOAT;

        SetScalarH5A((hid_t) obj_id, attr_name, dtype, values);
        return 0;
    } else if (isHDF5 == 2) {
        int dtype = NC_CHAR;

        if (data_type == DFNT_CHAR) dtype = NC_CHAR;
        else if (data_type == DFNT_UINT8) dtype = NC_UBYTE;
        else if (data_type == DFNT_INT16) dtype = NC_SHORT;
        else if (data_type == DFNT_UINT16) dtype = NC_SHORT;
        else if (data_type == DFNT_INT32) dtype = NC_INT;
        else if (data_type == DFNT_UINT32) dtype = NC_INT;
        else if (data_type == DFNT_FLOAT32) dtype = NC_FLOAT;

        int status = nc_put_att(obj_id, NC_GLOBAL, attr_name, dtype, count, values);
        check_err(status, __LINE__, __FILE__);
    } else {
        return SDsetattr(obj_id, attr_name, data_type, count, (VOIDP) values);
    }
    return 0;
}


int32_t put_smi(char *l3m_path,
        char *l3m_name,
        uint8 *l3m_data,
        int32 *dim_sizes,
        float32 *lat_range,
        float32 *lon_range,
        char *measure,
        char *scale_type,
        float32 *si_used,
        float32 *aminmax,
        char *atype,
        char *aopt,
        char *infiles,
        float32 *l3m_dminmax,
        meta_l3bType *meta_l3b,
        unsigned char *map_palette,
        char *softid,
        char *proc_con,
        instr input,
        char *precision,
        /*	       float32 *si8_used,*/
        uint8 *qual_byt,
        uint8 isHDF5,
        VOIDP fill_value);

/*-----------------------------------------------------------------------------
  Function: put_smi
  
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
  The function put_smi.c creates a level 3 Standard Map Image (SMI)
  file with the file name given by l3m_path.  Level 3 SMI file meta
  data, map_palette, and data will be written.  See "INTERFACE SPECIFICATIONS
  FOR SeaWiFS OPERATIONAL PRODUCT INPUT AND OUTPUT AND SENSOR CALIBRATION
  SOFTWARE", v3.3, dated 12th July 94, (By Fred Patt, Jim Firestone, 
  and Mike Darzi) and "SeaWiFS OPERATIONAL ARCHIVE PRODUCT 
  SPECIFICATIONS, v1.0 (Incomplete DRAFT version), dated 22 July 1994 
  (By F. Patt, J. Firestone, B. Schieber, L. Kumar, D. Ilg, and M. Darzi)
  for details.
  
  Arguments: (in calling order)
  Type           Name        I/O  Description
  ----           ----        ---  -----------
  char          *l3m_path    I    directory path & file name for SMI product
  char          *l3m_name    I    name of the geophysical parameter to which
                                  l3m_data corresponds
  float         *l3m_data    I    image data for parameter l3m_name
                                  product
  int32         *dim_sizes   I    #lines, #columns
  float32       *lat_range   I    lats of the outside edges of the northernmost
                                  and easternmost columns for l3mdata
  float32       *lon_range   I    lons of the outside edges of the westernmost
                                  and easternmost columns for l3m_data
  char          *measure     I    indicates whether it is mean,  variance,
                                  standard deviation, scenes or pixels   
  char          *scale_type  I    whether was scaled LOG or LINEAR
  float32       *si_used     I    slope & intercept used to scale data
  char          *infiles     I    name of input product
  float32       *l3m_dminmax I    data min/max found for L3-binned file
  meta_l3bType  *meta_l3b    I    meta data from input L3-binned file
  unsigned char *map_palette     I    color map_palette
  char          *softid      I    SeaDAS software ID
  char          *proc_con    I    SeaDAS processing log - command line
  
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
  Karen Settle    GSC              03/28/95    Generalized for SeaDAS.
                                               Used meta data argument
                                               for passing meta data from
                                               parent, added arguments,
                                               removed redundant arguments,
                                               expanded maplists.h 
                                               include file.
  --------------------------------------------------------------------------*/
int32_t put_smi(char *l3m_path, char*l3m_name, uint8 *l3m_data,
        int32 *dim_sizes,
        float32 *lat_range,
        float32 *lon_range,
        char *measure,
        char *scale_type,
        float32 *si_used,
        float32 *aminmax,
        char *atype,
        char *aopt,
        char *infiles,
        float32 *l3m_dminmax,
        meta_l3bType *meta_l3b,
        uint8 *map_palette,
        char *softid,
        char *proc_con,
        instr input,
        char *precision,
        /*	     float32 *si8_used,*/
        uint8 *qual_byt,
        uint8 isHDF5,
        VOIDP fill_value
        ) {
    /*-----------------------------------------------------------------------*/
    /*-----------------------------------------------------------------------*/
#include "smi_maplists.h"
    int32 fid, sdfid, sdsid;
    int32 i, rank = 2, mix, pal_ref;
    int32 start[2] = {0, 0};

    char attr_buffer[MAXVAL];

    char *loc_path, *str, loc_measure[MAXVAL];
    char data_center[MAXVAL];
    float32 lat_step, lon_step, swlat, swlon;
    hid_t h5fid, grp0;
    int ncid, status;
    float *lonarray;
    float *latarray;


    /*** Open/create output file in the given l3m_path */

    int32 valid_range[2] = {0, 2};

    if (isHDF5 == 1) {
        h5fid = H5Fcreate(l3m_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (h5fid < 0) {
            printf("\nError creating file %s \n", l3m_path);
            exit(EXIT_FAILURE);
        }
        grp0 = H5Gopen1(h5fid, "/");
        sdfid = h5fid;
    } else if (isHDF5 == 2) {
        status = nc_create(l3m_path, NC_NETCDF4, &ncid);
        if (status != NC_NOERR) {
            printf("\nError creating file %s \n", l3m_path);
            exit(EXIT_FAILURE);
        }
        sdfid = ncid;
    } else {
        sdfid = SDstart(l3m_path, DFACC_CREATE);
        fid = Hopen(l3m_path, DFACC_RDWR, 0);

        if (sdfid < 0 || fid < 0) {
            printf("\nError creating file %s \n", l3m_path);
            exit(EXIT_FAILURE);
        }
    }

    /* Get data center info */
    strcpy(data_center, meta_l3b->data_center);

    /*-----------------------------------------------------------------------*/
    /***  Write out global attributes  */

    if ((loc_path = strrchr(l3m_path, '/')) != NULL)
        loc_path++;
    else
        loc_path = l3m_path;

    //  if (isHDF5 < 2) {
    if (isHDF5 < 1) {
        setAttr(isHDF5, sdfid, L3M_PNAME, DFNT_CHAR,
                strlen(loc_path) + 1, (VOIDP) loc_path);

        setAttr(isHDF5, sdfid, L3M_SENSOR_NAME, DFNT_CHAR,
                strlen(meta_l3b->sensor_name) + 1, (VOIDP) meta_l3b->sensor_name);

        setAttr(isHDF5, sdfid, L3M_SENSOR, DFNT_CHAR,
                strlen(meta_l3b->sensor) + 1, (VOIDP) meta_l3b->sensor);

        sprintf(meta_l3b->title, "%s%s", meta_l3b->sensor_name,
                " Level-3 Standard Mapped Image");
        setAttr(isHDF5, sdfid, "Title", DFNT_CHAR,
                strlen(meta_l3b->sensor_name) +
                strlen(" Level-3 Standard Mapped Image") + 1,
                meta_l3b->title);

        strcpy(meta_l3b->data_center, "NASA/GSFC OBPG");
        setAttr(isHDF5, sdfid, L3M_DCENTER, DFNT_CHAR,
                strlen(meta_l3b->data_center) + 1, meta_l3b->data_center);

        //    setAttr(isHDF5, sdfid, L3M_STATION,   DFNT_CHAR,
        //	    strlen(meta_l3b->station) + 1, (VOIDP)meta_l3b->station);

        //    ft_fld = meta_l3b->station_lat;
        //setAttr(isHDF5, sdfid, L3M_STLAT,     DFNT_FLOAT32, 1, (VOIDP)&ft_fld);

        //ft_fld = meta_l3b->station_lon;
        //setAttr(isHDF5, sdfid, L3M_STLON,     DFNT_FLOAT32, 1, (VOIDP)&ft_fld);

        setAttr(isHDF5, sdfid, L3M_MISSION, DFNT_CHAR,
                strlen(meta_l3b->mission) + 1, (VOIDP) meta_l3b->mission);

        setAttr(isHDF5, sdfid, L3M_MSNCHAR, DFNT_CHAR,
                strlen(meta_l3b->mission_char) + 1, (VOIDP) meta_l3b->mission_char);

        setAttr(isHDF5, sdfid, L3M_SNSCHAR, DFNT_CHAR,
                strlen(meta_l3b->sensor_char) + 1, (VOIDP) meta_l3b->sensor_char);

        setAttr(isHDF5, sdfid, L3M_PRODTYPE, DFNT_CHAR,
                strlen(meta_l3b->prod_type) + 1, (VOIDP) meta_l3b->prod_type);

        setAttr(isHDF5, sdfid, L3M_PVERSION, DFNT_CHAR,
                strlen(input.pversion) + 1, (VOIDP) input.pversion);

        setAttr(isHDF5, sdfid, L3M_SOFTNM, DFNT_CHAR,
                strlen(L3M_SOFTNM_VAL) + 1, (VOIDP) L3M_SOFTNM_VAL);

        setAttr(isHDF5, sdfid, L3M_SOFTVER, DFNT_CHAR,
                strlen(L3M_SOFTVER_VAL) + 1, (VOIDP) L3M_SOFTVER_VAL);


        get_time(meta_l3b->ptime);
        setAttr(isHDF5, sdfid, L3M_PTIME, DFNT_CHAR,
                strlen(meta_l3b->ptime) + 1, (VOIDP) meta_l3b->ptime);

        if ((str = strrchr(infiles, '/')) != NULL)
            str++;
        else
            str = infiles;

        setAttr(isHDF5, sdfid, L3M_INFILES, DFNT_CHAR,
                strlen(str) + 1, (VOIDP) str);

        setAttr(isHDF5, sdfid, L3M_PROCCON, DFNT_CHAR,
                strlen(proc_con) + 1, (VOIDP) proc_con);

        setAttr(isHDF5, sdfid, L3M_INPARMS, DFNT_CHAR,
                strlen(input.parms) + 1, (VOIDP) input.parms);

        setAttr(isHDF5, sdfid, L3M_FLAG_NAMES, DFNT_CHAR,
                strlen(meta_l3b->flag_names) + 1, (VOIDP) meta_l3b->flag_names);

        int16 syear, sday, eyear, eday;
        double msec;
        int32 smsec, emsec;
        unix2yds(meta_l3b->startTime, &syear, &sday, &msec);
        smsec = (int32) (msec * 1000.);
        unix2yds(meta_l3b->endTime, &eyear, &eday, &msec);
        emsec = (int32) (msec * 1000.);

        setAttr(isHDF5, sdfid, L3M_PSYEAR, DFNT_INT16,
                1, (VOIDP) & syear);

        setAttr(isHDF5, sdfid, L3M_PSDAY, DFNT_INT16,
                1, (VOIDP) & sday);

        setAttr(isHDF5, sdfid, L3M_PEYEAR, DFNT_INT16,
                1, (VOIDP) & eyear);

        setAttr(isHDF5, sdfid, L3M_PEDAY, DFNT_INT16,
                1, (VOIDP) & eday);

        strcpy(attr_buffer, ydhmsf(meta_l3b->startTime, 'G'));
        setAttr(isHDF5, sdfid, L3M_STIME, DFNT_CHAR, strlen(attr_buffer) + 1, (VOIDP) attr_buffer);

        strcpy(attr_buffer, ydhmsf(meta_l3b->endTime, 'G'));
        setAttr(isHDF5, sdfid, L3M_ETIME, DFNT_CHAR, strlen(attr_buffer) + 1, (VOIDP) attr_buffer);


        setAttr(isHDF5, sdfid, L3M_SYEAR, DFNT_INT16, 1, (VOIDP) & syear);
        setAttr(isHDF5, sdfid, L3M_SDAY, DFNT_INT16, 1, (VOIDP) & sday);
        setAttr(isHDF5, sdfid, L3M_SMSEC, DFNT_INT32, 1, (VOIDP) & smsec);
        setAttr(isHDF5, sdfid, L3M_EYEAR, DFNT_INT16, 1, (VOIDP) & eyear);
        setAttr(isHDF5, sdfid, L3M_EDAY, DFNT_INT16, 1, (VOIDP) & eday);
        setAttr(isHDF5, sdfid, L3M_EMSEC, DFNT_INT32, 1, (VOIDP) & emsec);
        setAttr(isHDF5, sdfid, L3M_SORBIT, DFNT_INT32, 1, (VOIDP) & meta_l3b->start_orb);
        setAttr(isHDF5, sdfid, L3M_EORBIT, DFNT_INT32, 1, (VOIDP) & meta_l3b->end_orb);
        //    setAttr(isHDF5, sdfid, L3M_ORBIT,     DFNT_INT32, 1, (VOIDP)&meta_l3b->orbit);
        setAttr(isHDF5, sdfid, L3M_MAPPROJ, DFNT_CHAR, strlen(L3M_MAPPROJ_VAL) + 1,
                (VOIDP) L3M_MAPPROJ_VAL);
        setAttr(isHDF5, sdfid, L3M_LATUNITS, DFNT_CHAR, strlen(L3M_LATUNITS_VAL) + 1,
                (VOIDP) L3M_LATUNITS_VAL);
        setAttr(isHDF5, sdfid, L3M_LONUNITS, DFNT_CHAR, strlen(L3M_LONUNITS_VAL) + 1,
                (VOIDP) L3M_LONUNITS_VAL);

        /*** check validity of given latitude ranges-- 1st value should be greater
             than the second one.  If in error, output error message  	*/
        if (lat_range[0] <= lat_range[1]) {
            printf("\nput_smi: Given latitude range is in error \n");
            return FAIL;
        }
        setAttr(isHDF5, sdfid, L3M_NLAT, DFNT_FLOAT32, 1, (VOIDP) & lat_range[0]);
        setAttr(isHDF5, sdfid, L3M_SLAT, DFNT_FLOAT32, 1, (VOIDP) & lat_range[1]);
        setAttr(isHDF5, sdfid, L3M_WLON, DFNT_FLOAT32, 1, (VOIDP) & lon_range[0]);
        setAttr(isHDF5, sdfid, L3M_ELON, DFNT_FLOAT32, 1, (VOIDP) & lon_range[1]);

        lat_step = (lat_range[0] - lat_range[1]) / dim_sizes[0];
        setAttr(isHDF5, sdfid, L3M_LAT_STEP, DFNT_FLOAT32, 1, (VOIDP) & lat_step);

        if (lon_range[1] > lon_range[0])
            lon_step = (lon_range[1] - lon_range[0]) / dim_sizes[1];
        else
            lon_step = (360 - lon_range[0] + lon_range[1]) / dim_sizes[1];
        setAttr(isHDF5, sdfid, L3M_LON_STEP, DFNT_FLOAT32, 1, (VOIDP) & lon_step);

        swlat = (lat_range[1] + lat_step / 2.0);
        swlon = (lon_range[0] + lon_step / 2.0);
        setAttr(isHDF5, sdfid, L3M_SWLAT, DFNT_FLOAT32, 1, (VOIDP) & swlat);
        setAttr(isHDF5, sdfid, L3M_SWLON, DFNT_FLOAT32, 1, (VOIDP) & swlon);

        int32 tmpInt = meta_l3b->data_bins;
        setAttr(isHDF5, sdfid, L3M_DATABINS, DFNT_INT32, 1, (VOIDP) &tmpInt);
        setAttr(isHDF5, sdfid, L3M_NROWS, DFNT_INT32, 1, (VOIDP) & dim_sizes[0]);
        setAttr(isHDF5, sdfid, L3M_NCOLS, DFNT_INT32, 1, (VOIDP) & dim_sizes[1]);

        /*
          setAttr(isHDF5, sdfid, L3M_PARAMETER,DFNT_CHAR, strlen(parmname_list[index])+1,
          (VOIDP)parmname_list[index]);
         */
        setAttr(isHDF5, sdfid, L3M_PARAMETER, DFNT_CHAR, strlen(input.proddesc) + 1,
                (VOIDP) input.proddesc);

        for (i = 0; i < (int32_t) strlen(measure); i++)
            loc_measure[i] = measure[i];
        loc_measure[i] = (char) NULL;

        mix = 0;
        while ((strcmp(loc_measure, measure_list[mix]) != 0) &&
                (mix < L3M_MEASURES)) {
            mix++;
        }
        if (mix >= L3M_MEASURES) {
            printf("\nput_smi: Error: measure - %s unrecognized\n", measure);
            return FAIL;
        }

        setAttr(isHDF5, sdfid, L3M_MEASURE, DFNT_CHAR, strlen(measure_list[mix]) + 1,
                (VOIDP) measure_list[mix]);
        /*
          setAttr(isHDF5, sdfid, L3M_UNITS,      DFNT_CHAR, strlen(unit_list[index]) + 1, 
          (VOIDP)unit_list[index]);
         */
        setAttr(isHDF5, sdfid, L3M_UNITS, DFNT_CHAR, strlen(input.units) + 1,
                (VOIDP) input.units);

        if (strcmp(scale_type, "LOG") == 0) {
            setAttr(isHDF5, sdfid, L3M_SCALING, DFNT_CHAR, strlen(L3M_LOG_SCALE) + 1,
                    (VOIDP) L3M_LOG_SCALE);
            setAttr(isHDF5, sdfid, L3M_SC_EQN, DFNT_CHAR, strlen(L3M_LOG_EQN) + 1,
                    (VOIDP) L3M_LOG_EQN);
            setAttr(isHDF5, sdfid, L3M_BASE, DFNT_FLOAT32, 1, (VOIDP) & base);
        } else {
            setAttr(isHDF5, sdfid, L3M_SCALING, DFNT_CHAR, strlen(L3M_LINEAR_SCALE) + 1,
                    (VOIDP) L3M_LINEAR_SCALE);
            setAttr(isHDF5, sdfid, L3M_SC_EQN, DFNT_CHAR, strlen(L3M_LINEAR_EQN) + 1,
                    (VOIDP) L3M_LINEAR_EQN);
        }

        setAttr(isHDF5, sdfid, L3M_SLOPE, DFNT_FLOAT32, 1, (VOIDP) & si_used[0]);
        setAttr(isHDF5, sdfid, L3M_INTERCEPT, DFNT_FLOAT32, 1,
                (VOIDP) & si_used[1]);

        setAttr(isHDF5, sdfid, L3M_MIN, DFNT_FLOAT32, 1, (VOIDP) & l3m_dminmax[0]);
        setAttr(isHDF5, sdfid, L3M_MAX, DFNT_FLOAT32, 1, (VOIDP) & l3m_dminmax[1]);

        setAttr(isHDF5, sdfid, "Suggested Image Scaling Minimum", DFNT_FLOAT32, 1,
                (VOIDP) & aminmax[0]);
        setAttr(isHDF5, sdfid, "Suggested Image Scaling Maximum", DFNT_FLOAT32, 1,
                (VOIDP) & aminmax[1]);
        setAttr(isHDF5, sdfid, "Suggested Image Scaling Type", DFNT_CHAR, strlen(atype) + 1,
                (VOIDP) atype);
        setAttr(isHDF5, sdfid, "Suggested Image Scaling Applied", DFNT_CHAR, strlen(aopt) + 1,
                (VOIDP) aopt);

        setAttr(isHDF5, sdfid, "_lastModified", DFNT_CHAR,
                strlen(meta_l3b->ptime) + 1, (VOIDP) meta_l3b->ptime);

    } else {
        // NETCDF4
        setAttr(isHDF5, sdfid, "product_name", DFNT_CHAR,
                strlen(loc_path) + 1, (VOIDP) loc_path);

        setAttr(isHDF5, sdfid, "instrument", DFNT_CHAR,
                strlen(meta_l3b->sensor_name) + 1, (VOIDP) meta_l3b->sensor_name);

        sprintf(meta_l3b->title, "%s%s", meta_l3b->sensor_name,
                " Level-3 Standard Mapped Image");
        setAttr(isHDF5, sdfid, "title", DFNT_CHAR,
                strlen(meta_l3b->sensor_name) +
                strlen(" Level-3 Standard Mapped Image") + 1,
                meta_l3b->title);

        char buf[2048];

        strcpy(buf, "Ocean Biology Processing Group (NASA/GSFC/OBPG)");
        status = setAttr(isHDF5, sdfid, "project", DFNT_CHAR, strlen(buf) + 1, buf);

        if (strcmp(meta_l3b->mission, "") != 0)
            setAttr(isHDF5, sdfid, "platform", DFNT_CHAR,
                strlen(meta_l3b->mission) + 1, (VOIDP) meta_l3b->mission);

        setAttr(isHDF5, sdfid, "temporal_range", DFNT_CHAR,
                strlen(meta_l3b->prod_type) + 1, (VOIDP) meta_l3b->prod_type);

        setAttr(isHDF5, sdfid, "processing_version", DFNT_CHAR,
                strlen(input.pversion) + 1, (VOIDP) input.pversion);

        strcpy(meta_l3b->ptime, unix2isodate(time(NULL), 'G'));

        setAttr(isHDF5, sdfid, "date_created", DFNT_CHAR,
                strlen(meta_l3b->ptime) + 1, (VOIDP) meta_l3b->ptime);

        setAttr(isHDF5, sdfid, "history", DFNT_CHAR,
                strlen(proc_con) + 1, (VOIDP) proc_con);

        setAttr(isHDF5, sdfid, "l2_flag_names", DFNT_CHAR,
                strlen(meta_l3b->flag_names) + 1, (VOIDP) meta_l3b->flag_names);

        strcpy(attr_buffer, unix2isodate(meta_l3b->startTime, 'G'));
        setAttr(isHDF5, sdfid, "time_coverage_start", DFNT_CHAR, strlen(attr_buffer) + 1, attr_buffer);

        strcpy(attr_buffer, unix2isodate(meta_l3b->endTime, 'G'));
        setAttr(isHDF5, sdfid, "time_coverage_end", DFNT_CHAR, strlen(attr_buffer) + 1, attr_buffer);

        setAttr(isHDF5, sdfid, "start_orbit_number", DFNT_INT32, 1, (VOIDP) & meta_l3b->start_orb);
        setAttr(isHDF5, sdfid, "end_orbit_number", DFNT_INT32, 1, (VOIDP) & meta_l3b->end_orb);

        setAttr(isHDF5, sdfid, "map_projection", DFNT_CHAR, strlen(L3M_MAPPROJ_VAL) + 1,
                (VOIDP) L3M_MAPPROJ_VAL);

        strcpy(attr_buffer, "degrees_north");
        setAttr(isHDF5, sdfid, "latitude_units", DFNT_CHAR, strlen(attr_buffer) + 1,
                (VOIDP) attr_buffer);

        strcpy(attr_buffer, "degrees_east");
        setAttr(isHDF5, sdfid, "longitude_units", DFNT_CHAR, strlen(attr_buffer) + 1,
                (VOIDP) attr_buffer);

        /*** check validity of given latitude ranges-- 1st value should be greater
             than the second one.  If in error, output error message  	*/
        if (lat_range[0] <= lat_range[1]) {
            printf("\nput_smi: Given latitude range is in error \n");
            return FAIL;
        }

        setAttr(isHDF5, sdfid, "northernmost_latitude", DFNT_FLOAT32, 1, (VOIDP) & lat_range[0]);
        setAttr(isHDF5, sdfid, "southernmost_latitude", DFNT_FLOAT32, 1, (VOIDP) & lat_range[1]);
        setAttr(isHDF5, sdfid, "westernmost_longitude", DFNT_FLOAT32, 1, (VOIDP) & lon_range[0]);
        setAttr(isHDF5, sdfid, "easternmost_longitude", DFNT_FLOAT32, 1, (VOIDP) & lon_range[1]);
        setAttr(isHDF5, sdfid, "geospatial_lat_max", DFNT_FLOAT32, 1, (VOIDP) & lat_range[0]);
        setAttr(isHDF5, sdfid, "geospatial_lat_min", DFNT_FLOAT32, 1, (VOIDP) & lat_range[1]);
        setAttr(isHDF5, sdfid, "geospatial_lon_max", DFNT_FLOAT32, 1, (VOIDP) & lon_range[1]);
        setAttr(isHDF5, sdfid, "geospatial_lon_min", DFNT_FLOAT32, 1, (VOIDP) & lon_range[0]);

        setAttr(isHDF5, sdfid, "grid_mapping_name", DFNT_CHAR, 18, (VOIDP) "latitude_longitude");

        lat_step = (lat_range[0] - lat_range[1]) / dim_sizes[0];
        setAttr(isHDF5, sdfid, "latitude_step", DFNT_FLOAT32, 1, (VOIDP) & lat_step);

        if (lon_range[1] > lon_range[0])
            lon_step = (lon_range[1] - lon_range[0]) / dim_sizes[1];
        else
            lon_step = (360 - lon_range[0] + lon_range[1]) / dim_sizes[1];
        setAttr(isHDF5, sdfid, "longitude_step", DFNT_FLOAT32, 1, (VOIDP) & lon_step);

        if ((lonarray = (float *) malloc(dim_sizes[1] * sizeof (float))) == NULL) {
            printf(
                    "-E- : Error allocating memory for lon array\n");
            exit(EXIT_FAILURE);
        }
        if ((latarray = (float *) malloc(dim_sizes[0] * sizeof (float))) == NULL) {
            printf(
                    "-E- : Error allocating memory for lat array\n");
            exit(EXIT_FAILURE);
        }
        for (i = 0; i < dim_sizes[1]; i++) {
            lonarray[i] = lon_range[0] + lon_step * i + lon_step / 2.0;
        }
        for (i = 0; i < dim_sizes[0]; i++) {
            latarray[i] = lat_range[0] - lat_step * i - lat_step / 2.0;
        }
        swlat = (lat_range[1] + lat_step / 2.0);
        swlon = (lon_range[0] + lon_step / 2.0);
        setAttr(isHDF5, sdfid, "sw_point_latitude", DFNT_FLOAT32, 1, (VOIDP) & swlat);
        setAttr(isHDF5, sdfid, "sw_point_longitude", DFNT_FLOAT32, 1, (VOIDP) & swlon);

        float geo_resolution;
        char geospatial_units[8];
        if (strcmp(input.resolution, "LAND") == 0) {
            geo_resolution = 4.6;
            strcpy(geospatial_units, "km");
        } else if (strcmp(input.resolution, "9km") == 0) {
            geo_resolution = 9.2;
            strcpy(geospatial_units, "km");
        } else if (strcmp(input.resolution, "18km") == 0) {
            geo_resolution = 18.5;
            strcpy(geospatial_units, "km");
        } else if (strcmp(input.resolution, "4km") == 0) {
            geo_resolution = 4.6;
            strcpy(geospatial_units, "km");
        } else if (strcmp(input.resolution, "2km") == 0) {
            geo_resolution = 2.3;
            strcpy(geospatial_units, "km");
        } else if (strcmp(input.resolution, "1km") == 0) {
            geo_resolution = 1.2;
            strcpy(geospatial_units, "km");
        } else if (strcmp(input.resolution, "hkm") == 0) {
            geo_resolution = 575;
            strcpy(geospatial_units, "m");
        } else if (strcmp(input.resolution, "qkm") == 0) {
            geo_resolution = 288;
            strcpy(geospatial_units, "m");
        } else if (strcmp(input.resolution, "36km") == 0) {
            geo_resolution = 36;
            strcpy(geospatial_units, "km");
        } else if (strcmp(input.resolution, "90km") == 0) {
            geo_resolution = 90;
            strcpy(geospatial_units, "km");
        } else if (strcmp(input.resolution, "1deg") == 0) {
            geo_resolution = 1;
            strcpy(geospatial_units, "deg");
        } else if (strcmp(input.resolution, "hdeg") == 0) {
            geo_resolution = 0.5;
            strcpy(geospatial_units, "deg");
        } else if (strcmp(input.resolution, "qdeg") == 0) {
            geo_resolution = 0.25;
            strcpy(geospatial_units, "deg");
        } else if (strcmp(input.resolution, "10deg") == 0) {
            geo_resolution = 10;
            strcpy(geospatial_units, "deg");
        } else if (strncmp(input.resolution, "udeg", 4) == 0) {
            geo_resolution = (float) atof(&input.resolution[5]);
            strcpy(geospatial_units, "deg");
        } else if (strncmp(input.resolution, "ukm", 3) == 0) {
            geo_resolution = (float) atof(&input.resolution[4]);
            strcpy(geospatial_units, "km");
        }

        setAttr(isHDF5, sdfid, "geospatial_lon_resolution", DFNT_FLOAT, 1, &geo_resolution);
        setAttr(isHDF5, sdfid, "geospatial_lat_resolution", DFNT_FLOAT, 1, &geo_resolution);


        setAttr(isHDF5, sdfid, "geospatial_lat_units", DFNT_CHAR,
                strlen(geospatial_units) + 1, geospatial_units);

        setAttr(isHDF5, sdfid, "geospatial_lon_units", DFNT_CHAR,
                strlen(geospatial_units) + 1, geospatial_units);

        sprintf(buf, "%3.2f %s", geo_resolution, geospatial_units);
        setAttr(isHDF5, sdfid, "spatialResolution", DFNT_CHAR,
                strlen(buf) + 1, buf);

        int32 tmpInt = meta_l3b->data_bins;
        setAttr(isHDF5, sdfid, "data_bins", DFNT_INT32, 1, (VOIDP) &tmpInt);
        setAttr(isHDF5, sdfid, "number_of_lines", DFNT_INT32, 1, (VOIDP) & dim_sizes[0]);
        setAttr(isHDF5, sdfid, "number_of_columns", DFNT_INT32, 1, (VOIDP) & dim_sizes[1]);

        for (i = 0; i < (int32_t) strlen(measure); i++)
            loc_measure[i] = measure[i];
        loc_measure[i] = (char) NULL;

        mix = 0;
        while ((strcmp(loc_measure, measure_list[mix]) != 0) &&
                (mix < L3M_MEASURES)) {
            mix++;
        }
        if (mix >= L3M_MEASURES) {
            printf("\nput_smi: Error: measure - %s unrecognized\n", measure);
            return FAIL;
        }

        setAttr(isHDF5, sdfid, "measure", DFNT_CHAR, strlen(measure_list[mix]) + 1,
                (VOIDP) measure_list[mix]);

        setAttr(isHDF5, sdfid, "data_minimum", DFNT_FLOAT32, 1, (VOIDP) & l3m_dminmax[0]);
        setAttr(isHDF5, sdfid, "data_maximum", DFNT_FLOAT32, 1, (VOIDP) & l3m_dminmax[1]);

        setAttr(isHDF5, sdfid, "suggested_image_scaling_minimum", DFNT_FLOAT32, 1,
                (VOIDP) & aminmax[0]);
        setAttr(isHDF5, sdfid, "suggested_image_scaling_maximum", DFNT_FLOAT32, 1,
                (VOIDP) & aminmax[1]);
        setAttr(isHDF5, sdfid, "suggested_image_scaling_type", DFNT_CHAR, strlen(atype) + 1,
                (VOIDP) atype);
        setAttr(isHDF5, sdfid, "suggested_image_scaling_applied", DFNT_CHAR, strlen(aopt) + 1,
                (VOIDP) aopt);

        setAttr(isHDF5, sdfid, "_lastModified", DFNT_CHAR,
                strlen(meta_l3b->ptime) + 1, (VOIDP) meta_l3b->ptime);

        setAttr(isHDF5, sdfid, "Conventions", DFNT_CHAR, strlen("CF-1.6") + 1,
                (VOIDP) "CF-1.6");

        setAttr(isHDF5, sdfid, "institution", DFNT_CHAR,
                strlen("NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group") + 1,
                (VOIDP) "NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group");

        char const *standard_name_vocabulary =
                "NetCDF Climate and Forecast (CF) Metadata Convention";
        strcpy(buf, standard_name_vocabulary);
        status = setAttr(isHDF5, sdfid, "standard_name_vocabulary", DFNT_CHAR,
                strlen(buf) + 1, buf);

        char const *Metadata_Conventions = "Unidata Dataset Discovery v1.0";
        strcpy(buf, Metadata_Conventions);
        status = setAttr(isHDF5, sdfid, "Metadata_Conventions", DFNT_CHAR, strlen(buf) + 1, buf);

        char const *naming_authority = "gov.nasa.gsfc.sci.oceandata";
        strcpy(buf, naming_authority);
        status = setAttr(isHDF5, sdfid, "naming_authority", DFNT_CHAR, strlen(buf) + 1, buf);
        // create id
        strcpy(buf, meta_l3b->pversion);
        if (strcmp(meta_l3b->pversion, "Unspecified") != 0) {
            strcpy(buf, meta_l3b->product_name);
            strcat(buf, "/L3/");

        } else {
            strcpy(buf, "L3/");
        }
        strcat(buf, meta_l3b->product_name);
        status = setAttr(isHDF5, sdfid, "id", DFNT_CHAR, strlen(buf) + 1, buf);

        char const *license =
                "http://science.nasa.gov/earth-science/earth-science-data/data-information-policy/";
        strcpy(buf, license);
        status = setAttr(isHDF5, sdfid, "license", DFNT_CHAR, strlen(buf) + 1, buf);

        strcpy(buf, "NASA/GSFC/OBPG");
        status = setAttr(isHDF5, sdfid, "creator_name", DFNT_CHAR, strlen(buf) + 1, buf);
        status = setAttr(isHDF5, sdfid, "publisher_name", DFNT_CHAR, strlen(buf) + 1, buf);

        strcpy(buf, "data@oceancolor.gsfc.nasa.gov");
        status = setAttr(isHDF5, sdfid, "creator_email", DFNT_CHAR, strlen(buf) + 1, buf);
        status = setAttr(isHDF5, sdfid, "publisher_email", DFNT_CHAR, strlen(buf) + 1, buf);

        strcpy(buf, "http://oceandata.sci.gsfc.nasa.gov");
        status = setAttr(isHDF5, sdfid, "creator_url", DFNT_CHAR, strlen(buf) + 1, buf);
        status = setAttr(isHDF5, sdfid, "publisher_url", DFNT_CHAR, strlen(buf) + 1, buf);

        strcpy(buf, "L3 Mapped");
        status = setAttr(isHDF5, sdfid, "processing_level", DFNT_CHAR, strlen(buf) + 1, buf);

        strcpy(buf, "grid");
        status = setAttr(isHDF5, sdfid, "cdm_data_type", DFNT_CHAR, strlen(buf) + 1, buf);

        strcpy(buf, "http://dx.doi.org");
        status = setAttr(isHDF5, sdfid, "identifier_product_doi_authority", DFNT_CHAR,
                strlen(buf) + 1, buf);

        if (strcmp(meta_l3b->sensor_name, "CZCS") == 0)
            strcpy(buf, "10.5067/NIMBUS-7/CZCS_OC.2014.0");
        if (strcmp(meta_l3b->sensor_name, "OCTS") == 0)
            strcpy(buf, "10.5067/ADEOS/OCTS_OC.2014.0");
        if (strcmp(meta_l3b->sensor_name, "SEAWIFS") == 0)
            strcpy(buf, "10.5067/ORBVIEW-2/SEAWIFS_OC.2014.0");
        if (strstr(meta_l3b->sensor_name, "MODIS") != NULL) {
            if (strstr(meta_l3b->mission, "Aqua") != NULL)
                strcpy(buf, "10.5067/AQUA/MODIS_OC.2014.0");
            if (strstr(meta_l3b->mission, "Terra") != NULL)
                strcpy(buf, "10.5067/TERRA/MODIS_OC.2014.0");
        }


        status = setAttr(isHDF5, sdfid, "identifier_product_doi", DFNT_CHAR,
                strlen(buf) + 1, buf);

        char const *keyword_oc =
                "Oceans > Ocean Chemistry > Chlorophyll; Oceans > Ocean Optics > Ocean Color";
        char const *keyword_sst =
                "Oceans > Ocean Temperature > Sea Surface Temperature";

        char const *keyword_salinity =
                "SURFACE SALINITY, SALINITY, AQUARIUS, Jet Propulsion Laboratory, NASA, http://aquarius.nasa.gov/, AQUARIUS SAC-D, Aquarius Scatterometer, Aquarius Radiometer";

        if (strstr(input.prod, "sst") != NULL)
            strcpy(buf, keyword_sst);
        else
            strcpy(buf, keyword_oc);

        if (strstr(input.prod, "SSS") != NULL)
            strcpy(buf, keyword_salinity);

        status = setAttr(isHDF5, sdfid, "keywords", DFNT_CHAR, strlen(buf) + 1, buf);

        char const *keywords_vocabulary =
                "NASA Global Change Master Directory (GCMD) Science Keywords";
        strcpy(buf, keywords_vocabulary);
        status = setAttr(isHDF5, sdfid, "keywords_vocabulary", DFNT_CHAR,
                strlen(buf) + 1, buf);

        int grpid;
        if (isHDF5 == 2)
            nc_def_grp(sdfid, "processing_control", &grpid);
        else
            grpid = sdfid;

        setAttr(isHDF5, grpid, "software_name", DFNT_CHAR,
                strlen(L3M_SOFTNM_VAL) + 1, (VOIDP) L3M_SOFTNM_VAL);

        setAttr(isHDF5, grpid, "software_version", DFNT_CHAR,
                strlen(L3M_SOFTVER_VAL) + 1, (VOIDP) L3M_SOFTVER_VAL);

        if ((str = strrchr(infiles, '/')) != NULL)
            str++;
        else
            str = infiles;
        setAttr(isHDF5, grpid, "source", DFNT_CHAR,
                strlen(str) + 1, (VOIDP) str);

        //    setAttr(isHDF5, grpid, "l2_flag_names",   DFNT_CHAR,
        //	    strlen(meta_l3b->flag_names) + 1, (VOIDP)meta_l3b->flag_names);


        if (isHDF5 == 2) {
            nc_def_grp(grpid, "input_parameters", &grpid);
            char *end_str;
            char *token = strtok_r(input.parms, "|", &end_str);
            char tmp_buf[2048];
            while (token != NULL) {
                char *end_token;
                strcpy(tmp_buf, token);
                char *name = strtok_r(token, "=", &end_token);
                for (uint32_t i = 0; i < strlen(name); i++) {
                    if (name[i] == ' ') {
                        name[i] = 0;
                        break;
                    }
                }
                strcpy(tmp_buf, strtok_r(NULL, "|", &end_token));

                setAttr(isHDF5, grpid, name, DFNT_CHAR,
                        strlen(tmp_buf), (VOIDP) & tmp_buf[1]);
                token = strtok_r(NULL, "|", &end_str);
            }
        }
    }
    // HDF5
    if (isHDF5 == 1) {
        hid_t dataset, datatype;
        hsize_t dims[2] = {(hsize_t) dim_sizes[0], (hsize_t) dim_sizes[1]};
        hid_t dataspace = H5Screate_simple(rank, dims, NULL);

        if (strcmp(precision, "F") == 0) {
            datatype = H5T_NATIVE_FLOAT;
        } else if (strcmp(precision, "I") == 0) {
            datatype = H5T_NATIVE_USHORT;
        } else {
            datatype = H5T_NATIVE_UCHAR;
        }

        if ((dataset = H5Dcreate1(grp0, SDS_NAME, datatype, dataspace,
                H5P_DEFAULT)) < 0) {
            printf("\nput_smi: H5Dcreate call to create %s is unsuccessful",
                    SDS_NAME);
            return FAIL;
        }

        sdsid = dataset;

        if ((H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                (VOIDP) l3m_data)) < 0) {
            printf("\nput_smi: H5Dwrite unsuccessful\n");
            return FAIL;
        }

        if (strcmp(scale_type, "LOG") == 0) {
            setAttr(isHDF5, sdsid, L3M_SCALING, DFNT_CHAR, strlen(L3M_LOG_SCALE) + 1,
                    (VOIDP) L3M_LOG_SCALE);
            setAttr(isHDF5, sdsid, L3M_SC_EQN, DFNT_CHAR, strlen(L3M_LOG_EQN) + 1,
                    (VOIDP) L3M_LOG_EQN);
            setAttr(isHDF5, sdsid, L3M_BASE, DFNT_FLOAT32, 1, (VOIDP) & base);
        } else {
            setAttr(isHDF5, sdsid, L3M_SCALING, DFNT_CHAR, strlen(L3M_LINEAR_SCALE) + 1,
                    (VOIDP) L3M_LINEAR_SCALE);
            setAttr(isHDF5, sdsid, L3M_SC_EQN, DFNT_CHAR, strlen(L3M_LINEAR_EQN) + 1,
                    (VOIDP) L3M_LINEAR_EQN);
        }

        setAttr(isHDF5, sdsid, L3M_SLOPE, DFNT_FLOAT32, 1, (VOIDP) & si_used[0]);
        setAttr(isHDF5, sdsid, L3M_INTERCEPT, DFNT_FLOAT32, 1,
                (VOIDP) & si_used[1]);

        setAttr(isHDF5, sdsid, "_FillValue", DFNT_FLOAT32, 1, fill_value);

        setAttr(isHDF5, sdsid, "scale_factor", DFNT_FLOAT32, 1, (VOIDP) & si_used[0]);
        setAttr(isHDF5, sdsid, "add_offset", DFNT_FLOAT32, 1,
                (VOIDP) & si_used[1]);

        H5Sclose(dataspace);
        H5Dclose(dataset);

        //  write map_palette
        dims[0] = 3;
        dims[1] = 256;
        dataspace = H5Screate_simple(2, dims, NULL);
        if ((dataset = H5Dcreate1(grp0, "palette", H5T_NATIVE_UCHAR, dataspace,
                H5P_DEFAULT)) < 0) {
            printf("\nput_smi: H5Dcreate call to create %s is unsuccessful",
                    "palette");
            return FAIL;
        }

        if ((H5Dwrite(dataset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                (VOIDP) map_palette)) < 0) {
            printf("\nput_smi: H5Dwrite unsuccessful\n");
            return FAIL;
        }

        //  close down interfaces and close file 
        H5Fclose(h5fid);

    } else if (isHDF5 == 2) {
        // NETCDF4
        int dimids[2];

        status = nc_def_dim(sdfid, "lat", dim_sizes[0], &dimids[0]);
        check_err(status, __LINE__, __FILE__);

        status = nc_def_dim(sdfid, "lon", dim_sizes[1], &dimids[1]);
        check_err(status, __LINE__, __FILE__);

        int datatype;
        if (strcmp(precision, "F") == 0) {
            datatype = NC_FLOAT;
        } else if (strcmp(precision, "I") == 0) {
            datatype = NC_USHORT;
        } else {
            datatype = NC_UBYTE;
        }

        int varid;
        /*
         * removed geophysical_data group to make SMI files more CF compliant...
         *
         *  status = nc_def_grp( sdfid, "geophysical_data", &grpid);
         */
        status = nc_def_var(sdfid, l3m_name, datatype, 2, dimids, &varid);
        check_err(status, __LINE__, __FILE__);

        /* Set compression */
        if (input.deflate > 0) {
            /* First set chunking */
            size_t chunksize[2] = {64, 64};
            status = nc_def_var_chunking(sdfid, varid, NC_CHUNKED, chunksize);
            check_err(status, __LINE__, __FILE__);
            if (status != NC_NOERR) exit(1);

            /* Now we can set compression */
            status = nc_def_var_deflate(sdfid, varid, NC_NOSHUFFLE, 1,
                    input.deflate);
            check_err(status, __LINE__, __FILE__);
            if (status != NC_NOERR) exit(1);
        }

        static productInfo_t *p_info;
        p_info = allocateProductInfo();

        int sensorId = sensorName2SensorId(meta_l3b->sensor_name);
        if (sensorId == -1)
            sensorId = instrumentPlatform2SensorId(meta_l3b->sensor_name, meta_l3b->mission);
        if (sensorId == -1) {
            printf("%s is unknown\n", meta_l3b->sensor_name);
            exit(EXIT_FAILURE);

        }
        if (!findProductInfo(l3m_name, sensorId, p_info)) {
            printf("%s not found in XML product table\n", l3m_name);
            exit(EXIT_FAILURE);
        }

        status = nc_put_att_text(sdfid, varid, "long_name",
                strlen(p_info->description) + 1,
                p_info->description);

        status = nc_put_att_text(sdfid, varid, "units",
                strlen(p_info->units) + 1,
                p_info->units);


        if (p_info->standardName != NULL) {
            status = nc_put_att_text(sdfid, varid, "standard_name",
                    strlen(p_info->standardName) + 1,
                    p_info->standardName);
        }

        if (datatype == NC_SHORT) {
            short fv_i16 = (int) p_info->fillValue;
            status = nc_put_att(sdfid, varid, "_FillValue", NC_SHORT, 1, &fv_i16);

            short int valid;
            valid = p_info->validMin;
            status = nc_put_att_short(sdfid, varid, "valid_min", NC_SHORT, 1,
                    &valid);
            valid = p_info->validMax;
            status = nc_put_att_short(sdfid, varid, "valid_max", NC_SHORT, 1,
                    &valid);


        } else if (datatype == NC_INT) {
            float fv_i = (int) p_info->fillValue;
            status = nc_put_att(sdfid, varid, "_FillValue", NC_INT, 1, &fv_i);

            int valid;
            valid = p_info->validMin;
            status = nc_put_att_int(sdfid, varid, "valid_min", NC_INT, 1,
                    &valid);
            valid = p_info->validMax;
            status = nc_put_att_int(sdfid, varid, "valid_max", NC_INT, 1,
                    &valid);


        } else if (datatype == NC_FLOAT) {
            float fv_f32 = (float) p_info->fillValue;
            status = nc_put_att(sdfid, varid, "_FillValue", NC_FLOAT, 1, &fv_f32);

            float valid;
            valid = p_info->validMin;
            status = nc_put_att_float(sdfid, varid, "valid_min", NC_FLOAT, 1,
                    &valid);
            valid = p_info->validMax;
            status = nc_put_att_float(sdfid, varid, "valid_max", NC_FLOAT, 1,
                    &valid);
        }
        check_err(status, __LINE__, __FILE__);


        status = nc_put_att_text(sdfid, varid, "display_scale",
                strlen(p_info->displayScale) + 1,
                p_info->displayScale);

        status = nc_put_att_double(sdfid, varid, "display_min", NC_DOUBLE, 1,
                &p_info->displayMin);
        status = nc_put_att_double(sdfid, varid, "display_max", NC_DOUBLE, 1,
                &p_info->displayMax);

        status = nc_put_att_float(sdfid, varid, "scale_factor", NC_FLOAT, 1,
                &si_used[0]);
        status = nc_put_att_float(sdfid, varid, "add_offset", NC_FLOAT, 1,
                &si_used[1]);

        if (&p_info->reference[0] != 0)
            status = nc_put_att_text(sdfid, varid, "reference",
                strlen(p_info->reference) + 1,
                p_info->reference);

        if (&p_info->comment[0] != 0)
            status = nc_put_att_text(sdfid, varid, "comment",
                strlen(p_info->comment) + 1,
                p_info->comment);

        status = nc_put_var(sdfid, varid, l3m_data);
        check_err(status, __LINE__, __FILE__);

        if (qual_byt != NULL) {
            int qualid;
            char l3m_qualname[128];
            strcpy(l3m_qualname, "qual_");
            strcat(l3m_qualname, l3m_name);
            status = nc_def_var(sdfid, l3m_qualname, NC_BYTE, 2, dimids, &qualid);
            check_err(status, __LINE__, __FILE__);

            status = nc_put_var(sdfid, qualid, qual_byt);
            check_err(status, __LINE__, __FILE__);
        }

        /*
         * Add lon and lat vectors
         */
        //Lat

        status = nc_def_var(sdfid, "lat", NC_FLOAT, 1, &dimids[0], &varid);
        check_err(status, __LINE__, __FILE__);
        if (!findProductInfo("lat", sensorId, p_info)) {
            printf("%s not found in XML product table\n", l3m_name);
            exit(EXIT_FAILURE);
        }
        status = nc_put_att_text(sdfid, varid, "long_name",
                strlen(p_info->description) + 1,
                p_info->description);

        status = nc_put_att_text(sdfid, varid, "units",
                strlen(p_info->units) + 1,
                p_info->units);


        if (p_info->standardName != NULL) {
            status = nc_put_att_text(sdfid, varid, "standard_name",
                    strlen(p_info->standardName) + 1,
                    p_info->standardName);
        }
        float fv_f32 = (float) p_info->fillValue;
        status = nc_put_att(sdfid, varid, "_FillValue", NC_FLOAT, 1, &fv_f32);
        float valid;
        valid = p_info->validMin;
        status = nc_put_att_float(sdfid, varid, "valid_min", NC_FLOAT, 1,
                &valid);
        valid = p_info->validMax;
        status = nc_put_att_float(sdfid, varid, "valid_max", NC_FLOAT, 1,
                &valid);
        status = nc_put_var(sdfid, varid, latarray);
        check_err(status, __LINE__, __FILE__);

        //Lon
        status = nc_def_var(sdfid, "lon", NC_FLOAT, 1, &dimids[1], &varid);
        check_err(status, __LINE__, __FILE__);
        if (!findProductInfo("lon", sensorId, p_info)) {
            printf("%s not found in XML product table\n", l3m_name);
            exit(EXIT_FAILURE);
        }
        status = nc_put_att_text(sdfid, varid, "long_name",
                strlen(p_info->description) + 1,
                p_info->description);

        status = nc_put_att_text(sdfid, varid, "units",
                strlen(p_info->units) + 1,
                p_info->units);


        if (p_info->standardName != NULL) {
            status = nc_put_att_text(sdfid, varid, "standard_name",
                    strlen(p_info->standardName) + 1,
                    p_info->standardName);
        }
        fv_f32 = (float) p_info->fillValue;
        status = nc_put_att(sdfid, varid, "_FillValue", NC_FLOAT, 1, &fv_f32);
        valid = p_info->validMin;
        status = nc_put_att_float(sdfid, varid, "valid_min", NC_FLOAT, 1,
                &valid);
        valid = p_info->validMax;
        status = nc_put_att_float(sdfid, varid, "valid_max", NC_FLOAT, 1,
                &valid);
        status = nc_put_var(sdfid, varid, lonarray);
        check_err(status, __LINE__, __FILE__);

        free(latarray);
        free(lonarray);

        //  write map_palette
        status = nc_def_dim(sdfid, "rgb", 3, &dimids[0]);
        check_err(status, __LINE__, __FILE__);

        status = nc_def_dim(sdfid, "eightbitcolor", 256, &dimids[1]);
        check_err(status, __LINE__, __FILE__);

        status = nc_def_var(sdfid, "palette", NC_UBYTE, 2, dimids, &varid);
        check_err(status, __LINE__, __FILE__);

        status = nc_put_var(sdfid, varid, map_palette);
        check_err(status, __LINE__, __FILE__);
        // close the file
        status = nc_close(sdfid);
        check_err(status, __LINE__, __FILE__);

    } else {
        // HDF4
        /*-----------------------------------------------------------------------*/
        /***  write  SDS */
        if (strcmp(precision, "F") == 0) {
            if ((sdsid = SDcreate(sdfid, SDS_NAME, DFNT_FLOAT32, rank, dim_sizes)) < 0) {
                printf("\nput_smi: SDcreate call to create %s is unsuccessful",
                        SDS_NAME);
                return FAIL;
            }
            setAttr(isHDF5, sdsid, "Fill", DFNT_FLOAT32, 1, fill_value);
        } else if (strcmp(precision, "I") == 0) {
            if ((sdsid = SDcreate(sdfid, SDS_NAME, DFNT_UINT16, rank, dim_sizes)) < 0) {
                printf("\nput_smi: SDcreate call to create %s is unsuccessful",
                        SDS_NAME);
                return FAIL;
            }
            setAttr(isHDF5, sdsid, "Fill", DFNT_UINT16, 1, fill_value);
        } else {
            if ((sdsid = SDcreate(sdfid, SDS_NAME, DFNT_UINT8, rank, dim_sizes)) < 0) {
                printf("\nput_smi: SDcreate call to create %s is unsuccessful",
                        SDS_NAME);
                return FAIL;
            }
            setAttr(isHDF5, sdsid, "Fill", DFNT_UINT8, 1, fill_value);
        }

        if ((SDwritedata(sdsid, start, NULL, dim_sizes, (VOIDP) l3m_data)) < 0) {
            printf("\nput_smi: SDwritedata unsuccessful\n");
            return FAIL;
        }

        if (strcmp(scale_type, "LOG") == 0) {
            setAttr(isHDF5, sdsid, L3M_SCALING, DFNT_CHAR, strlen(L3M_LOG_SCALE) + 1,
                    (VOIDP) L3M_LOG_SCALE);
            setAttr(isHDF5, sdsid, L3M_SC_EQN, DFNT_CHAR, strlen(L3M_LOG_EQN) + 1,
                    (VOIDP) L3M_LOG_EQN);
            setAttr(isHDF5, sdsid, L3M_BASE, DFNT_FLOAT32, 1, (VOIDP) & base);
        } else {
            setAttr(isHDF5, sdsid, L3M_SCALING, DFNT_CHAR, strlen(L3M_LINEAR_SCALE) + 1,
                    (VOIDP) L3M_LINEAR_SCALE);
            setAttr(isHDF5, sdsid, L3M_SC_EQN, DFNT_CHAR, strlen(L3M_LINEAR_EQN) + 1,
                    (VOIDP) L3M_LINEAR_EQN);
        }

        setAttr(isHDF5, sdsid, L3M_SLOPE, DFNT_FLOAT32, 1, (VOIDP) & si_used[0]);
        setAttr(isHDF5, sdsid, L3M_INTERCEPT, DFNT_FLOAT32, 1,
                (VOIDP) & si_used[1]);

        SDendaccess(sdsid);

        if (qual_byt != NULL) {
            if ((sdsid = SDcreate(sdfid, "l3m_qual", DFNT_UINT8, rank, dim_sizes)) < 0) {
                printf("\nput_smi: SDcreate call to create %s is unsuccessful",
                        SDS_NAME);
                return FAIL;
            }

            if ((SDwritedata(sdsid, start, NULL, dim_sizes, (VOIDP) qual_byt)) < 0) {
                printf("\nput_smi: SDwritedata unsuccessful (qual_byt)\n");
                return FAIL;
            }

            setAttr(isHDF5, sdsid, "valid_range", DFNT_INT32, 2, (VOIDP) & valid_range[0]);

            SDendaccess(sdsid);
        }



        /*-----------------------------------------------------------------------*/
        /*  write map_palette */

        if ((DFPaddpal(l3m_path, (VOIDP) map_palette)) < 0) {
            fprintf(stderr, "put_smi: Error writing map_palette.\n");
            return FAIL;
        }

        if ((pal_ref = DFPlastref()) > 0) {
            if ((DFANputlabel(l3m_path, DFTAG_IP8, pal_ref, (char*) "palette")) < 0) {
                printf("\nput_smi: Error writing palette label\n");
            }
        }

        /*-----------------------------------------------------------------------*/
        /***  close down interfaces and close file  */
        SDend(sdfid);
        Hclose(fid);
    }

    //  free(attr_buffer);
    return SUCCEED;
}

