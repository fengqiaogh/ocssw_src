#include "map.h"
#include "mapproto.h"
#include "mapattr.h"
#include <hdf.h>
#include <mfhdf.h>
#include "hdf5.h"

uint8_t isHDF5;

/*-----------------------------------------------------------------------------

 Function: get_l3m

 Returns: intn (status)
 The return code is FAIL (-1) if an error occurs, SUCCEED (0)
 otherwise.

 Description:
 The function get_l3m.c reads a level 3 Standard Map Image (SMI)
 file with the file name given by l3m_path.  Level 3 SMI file meta
 data, palette, data and calibration information will be retrieved.
 See "INTERFACE SPECIFICATIONS FOR SeaWiFS OPERATIONAL PRODUCT INPUT
 AND OUTPUT SOFTWARE" (By Fred Patt, Jim Firestone, and Mike Darzi)
 and "SeaWiFS Level 3 Standard Map CDL" (By Fred Patt) for details.

 Arguments: (in calling order)
 Type      Name        I/O   Description
 ----      ----        ---   -----------
 char *    l3m_path    I     directory path & file name for SMI product
 int16     syear       O     data start time year (from l3 bin product)
 int16     sday        O     data start time day-of-year
 int32     smsec       O     data start time milliseconds-of-day
 int16     eyear       O     data end time year; read from l3 bin product
 int16     eday        O     data end day-of-year
 int16     emsec       O     data end time milliseconds-of-day
 char *    prod_type   O     binning period description(scene, day.....)
 char *    l3m_name    O     name of the geophysical parameter to which
 l3m_data corresponds
 uint8 *   l3m_data    O     image data for parameter l3m_name
 uint8 *   palette     O     RGB wts for each of 256 gray-levels of the
 l3m_data byte values
 struct *  meta_l3m    O     structure containing station info, product
 name, input files, ptime, proc_con, proc_log,
 parameter info.

 Notes:

 Modification history:
 Programmer     Organization      Date      Description of change
 --------------   ------------    --------    ---------------------
 Lakshmi Kumar    Hughes STX      12/16/93    Original development
 Lakshmi Kumar	 Hughes STX	 08/30/94    Revised version
 Lakshmi Kumar	 Hughes STX	 05/18/95    Modified code for access
 ing flag_use as string
 Lakshmi Kumar	 Hughes STX	 09/28/95    Added orbit as output arg
 (see V4.4 I/O Specs.)
 Lakshmi Kumar	 Hughes STX	 10/27/95    Read "Start Orbit" &
 "End Orbit" global attrs.
 Lakshmi Kumar    Hughes STX      06/18/96    Changed defn. of MAX to
 MAXVAL inorder to remove
 compile time warning
 Lakshmi Kumar	 Hughes STX	 12/11/96    Removed orbit from inter
 face.  It is stored in
 meta_l3m structure
 Also removed non-ANSI
 declarations
 Joel Gales	 Futuretech	 09/27/13    Removed netcdf code
 HDF5 calls can be used
 for NETCDF4 files
 -----------------------------------------------------------------------------*/

int32_t get_l3m(char *l3m_path, int16_t *syear, int16_t *sday, int32_t *smsec,
        int16_t *eyear, int16_t *eday, int32_t *emsec, char *prod_type,
        char *l3m_name, uint8_t *l3m_data, unsigned char *palette,
        meta_struct *meta_l3m) {

    int i;
    int32 fid, sdfid, sdsid, dimsizes[2], start[2];
    int32 rank = 2, nattrs, numbertype, index = 0;
    char name[MAXVAL];
    size_t buf_size = 255;
    hid_t h5fid, dataset, datatype, filespace, grp0, attr, atype, atype_mem;
    hsize_t dims[2];
    hsize_t num_obj;
    H5O_info_t infobuf;

    start[0] = start[1] = 0;

    if (H5Fis_hdf5(l3m_path)) {
        isHDF5 = 1;
        h5fid = H5Fopen(l3m_path, H5F_ACC_RDONLY, H5P_DEFAULT);
        sdfid = h5fid;
    } else {
        isHDF5 = 0;
        sdfid = SDstart(l3m_path, DFACC_RDONLY);
        fid = Hopen(l3m_path, DFACC_RDONLY, 0);
    }

    /*  read global attributes  */
    if (meta_l3m != NULL) {

        if (isHDF5 == 1) {
            H5E_auto_t old_func;
            void *old_client_data;
            H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
            H5Eset_auto(H5E_DEFAULT, NULL, NULL); // Turn off error handling
            H5Gget_num_objs(h5fid, &num_obj);
            index = -1;
            for (i = 0; i < num_obj; i++) {
                H5Gget_objname_by_idx(h5fid, i, name, 80);
                H5Oget_info_by_name(h5fid, name, &infobuf, H5O_INFO_BASIC, H5P_DEFAULT);
                if (infobuf.type == H5O_TYPE_DATASET) {

                    if ((strcmp("lon", name) != 0) && (strcmp("lat", name) != 0)
                            && (strcmp("eightbitcolor", name) != 0)
                            && (strcmp("palette", name) != 0)
                            && (strcmp("rgb", name) != 0)
                            && (strstr(name, "qual") == NULL)
                            ) {
                        index = i;
                        break;
                    }
                }
            }
            // Restore previous error handler
            H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);
            dataset = H5Dopen1(h5fid, name);

            datatype = H5Dget_type(dataset);
            filespace = H5Dget_space(dataset);
            H5Sget_simple_extent_dims(filespace, dims, NULL);
            meta_l3m->nrows = dims[0];
            meta_l3m->ncols = dims[1];
            numbertype = datatype;
            H5Sclose(filespace);
            H5Dclose(dataset);

            grp0 = H5Gopen1(h5fid, "/");
            H5Oget_info(grp0, &infobuf, H5O_INFO_NUM_ATTRS);
            char string_out[80];
            H5T_class_t type_class;
            for (i = 0; i < infobuf.num_attrs; i++) {
                attr = H5Aopen_by_idx(grp0, "/", H5_INDEX_CRT_ORDER, H5_ITER_INC, (hsize_t) i, H5P_DEFAULT, H5P_DEFAULT);
                H5Aget_name(attr, buf_size, name);
                atype = H5Aget_type(attr);
                type_class = H5Tget_class(atype);
                if (type_class == H5T_FLOAT) {
                    atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
                    if ((strcmp("suggested_image_scaling_minimum", name) == 0)
                            || (strcmp("Suggested Image Scaling Minimum", name) == 0)) {
                        H5Aread(attr, atype_mem, &meta_l3m->scaled_data_min);
                    }
                    if ((strcmp("suggested_image_scaling_maximum", name) == 0)
                            || (strcmp("Suggested Image Scaling Maximum", name) == 0)) {
                        H5Aread(attr, atype_mem, &meta_l3m->scaled_data_max);
                    }
                    H5Aread(attr, atype_mem, string_out);
                } else if (type_class == H5T_STRING) {
                    atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
                    if ((strcmp("suggested_image_scaling_type", name) == 0) ||
                            (strcmp("Suggested Image Scaling Type", name) == 0)) {
                        if ((meta_l3m->scaled_data_type = (char *) malloc(
                                sizeof (char) * 7 + 1)) == NULL)
                            return FAIL;
                        H5Aread(attr, atype, meta_l3m->scaled_data_type);
                    }
                }
            }

            H5Aclose(attr);
            H5Gclose(grp0);
        } else {
            // HDF4
            if ((read_meta(sdfid, meta_l3m)) < 0)
                return FAIL;

            *syear = meta_l3m->syear;
            *sday = meta_l3m->sday;
            *smsec = meta_l3m->smsec;
            *eyear = meta_l3m->eyear;
            *eday = meta_l3m->eday;
            *emsec = meta_l3m->emsec;

            if (prod_type != NULL)
                strcpy(prod_type, meta_l3m->prodtype);

            if (l3m_name != NULL)
                strcpy(l3m_name, meta_l3m->parameter);

            if ((sdsid = SDselect(sdfid, 0)) < 0)
                return FAIL;
            if ((SDgetinfo(sdsid, name, &rank, dimsizes, &numbertype, &nattrs))
                    < 0)
                return FAIL;
            SDendaccess(sdsid);
        }

        meta_l3m->bintype = numbertype;
    }

    /*  read   SDS */
    if (l3m_data != NULL) {

        if (strcmp(l3m_name, "l3m_qual") == 0)
            index = 1;

        if (isHDF5 == 1) {
            H5E_auto_t old_func;
            void *old_client_data;
            H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
            H5Eset_auto(H5E_DEFAULT, NULL, NULL); // Turn off error handling
            H5Gget_num_objs(h5fid, &num_obj);
            index = -1;

            if (strcmp(l3m_name, "l3m_qual") != 0) {
                for (i = 0; i < num_obj; i++) {
                    H5Gget_objname_by_idx(h5fid, i, name, 80);
                    H5Oget_info_by_name(h5fid, name, &infobuf, H5O_INFO_BASIC, H5P_DEFAULT);
                    if (infobuf.type == H5O_TYPE_DATASET) {

                        if ((strcmp("lon", name) != 0) && (strcmp("lat", name) != 0)
                                && (strcmp("eightbitcolor", name) != 0)
                                && (strcmp("palette", name) != 0)
                                && (strcmp("rgb", name) != 0)
                                && (strstr(name, "qual") == NULL)) {
                            index = i;
                            break;
                        }
                    }
                }
            } else {
                for (i = 0; i < num_obj; i++) {
                    H5Gget_objname_by_idx(h5fid, i, name, 80);
                    H5Oget_info_by_name(h5fid, name, &infobuf, H5O_INFO_BASIC, H5P_DEFAULT);
                    if (infobuf.type == H5O_TYPE_DATASET) {
                        if ((strstr(name, "qual")) != NULL) {
                            index = i;
                            break;
                        }
                    }
                }
            }
            // Restore previous error handler
            H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

            if (index >= 0) {
                dataset = H5Dopen1(h5fid, name);
            } else {
                H5Fclose(h5fid);
                return 0;
            }

            datatype = H5Dget_type(dataset);
            sdsid = dataset;
            if ((H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    l3m_data)) < 0)
                return FAIL;

            H5Dclose(dataset);
        } else {
            if ((sdsid = SDselect(sdfid, index)) < 0)
                return FAIL;

            if ((SDgetinfo(sdsid, name, &rank, dimsizes, &numbertype, &nattrs))
                    < 0)
                return FAIL;

            if ((SDreaddata(sdsid, start, NULL, dimsizes, (VOIDP) l3m_data))
                    < 0)
                return FAIL;

            SDendaccess(sdsid);
        }
    }

    /*  read  palette  */
    if (isHDF5 && palette != NULL) {
        dataset = H5Dopen1(h5fid, "palette");
        datatype = H5Dget_type(dataset);
        sdsid = dataset;
        if ((H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, palette))
                < 0)
            return FAIL;
        H5Dclose(dataset);
    } else {
        if (palette != NULL) {
            if ((DFPgetpal(l3m_path, (VOIDP) palette)) < 0)
                return FAIL;
        }
    }

    /*  close down interfaces and close file  */
    if (isHDF5 == 1) {
        H5Fclose(h5fid);
    } else {
        SDend(sdfid);
        Hclose(fid);
    }

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
 Function:  read_meta

 Returns:   int32 (status)
 The return code is a negative value if any error occurs, otherwise,
 returns 0.

 Description:
 The function read_meta reads all the global attributes from the given
 map file

 Parameters: (in calling order)
 Type      Name        I/O   Description
 ----      ----        ---   -----------
 int32     sdfid        I    ID req to access HDF SDS interface
 meta_struct *meta_l3m  O    Product meta data passed to the calling
 function
 Modification history:
 Programmer     Organization   Date      Description of change
 -------------- ------------   --------  ---------------------
 Lakshmi Kumar  Hughes STX     11/23/94  Original development
 Lakshmi Kumar  Hughes STX	  02/27/94  Added code to read the new global
 attrs measure, flags, max & min
 (I/O specs V4.2 & product specs 2.6)
 ----------------------------------------------------------------------------*/
int32_t read_meta(int32_t sdfid, meta_struct *meta_l3m) {
    int32_t nt, count;

    /****  read global attribute "product_name"  */
    if ((getattrsz(sdfid, L3M_PNAME, &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->pname = (char *) malloc(sizeof (char) * count + 1)) == NULL)
        return FAIL;
    if ((rdattr(sdfid, L3M_PNAME, (VOIDP *) meta_l3m->pname)) < 0)
        return FAIL;
    meta_l3m->pname[count] = 0;

    /**** read global attribute "title" */
    if ((getattrsz(sdfid, L3M_TITLE, &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->title = (char *) malloc(sizeof (char) * count + 1)) == NULL)
        return FAIL;
    if ((rdattr(sdfid, L3M_TITLE, (VOIDP *) meta_l3m->title)) < 0)
        return FAIL;
    meta_l3m->title[count] = 0;

    /**** read global attribute "project" */
    if ((getattrsz(sdfid, L3M_DCENTER, &nt, &count)) >= 0) {
        if ((meta_l3m->dcenter = (char *) malloc(sizeof (char) * count + 1))
                == NULL)
            return FAIL;
        if ((rdattr(sdfid, L3M_DCENTER, (VOIDP *) meta_l3m->dcenter)) < 0)
            return FAIL;
        meta_l3m->dcenter[count] = 0;
    } else
        meta_l3m->dcenter = 0;

    /**** read global attribute "mission" */
    if ((getattrsz(sdfid, L3M_MISSION, &nt, &count)) >= 0) {
        if ((meta_l3m->mission = (char *) malloc(sizeof (char) * count + 1))
                == NULL)
            return FAIL;
        if ((rdattr(sdfid, L3M_MISSION, (VOIDP *) meta_l3m->mission)) < 0)
            return FAIL;
        meta_l3m->mission[count] = 0;
    } else
        meta_l3m->mission = 0;

    /**** read global attribute "sensor" */
    if ((getattrsz(sdfid, L3M_SENSOR_NAME, &nt, &count)) >= 0) {
        if ((meta_l3m->sensor_name = (char *) malloc(sizeof (char) * count + 1))
                == NULL)
            return FAIL;
        if ((rdattr(sdfid, L3M_SENSOR_NAME, (VOIDP *) meta_l3m->sensor_name))
                < 0)
            return FAIL;
        meta_l3m->sensor_name[count] = 0;
    } else
        meta_l3m->sensor_name = 0;

    /**** read global attribute "Data Type" */
    if ((getattrsz(sdfid, L3M_PRODTYPE, &nt, &count)) >= 0) {
        if ((meta_l3m->prodtype = (char *) malloc(sizeof (char) * count + 1))
                == NULL)
            return FAIL;
        if ((rdattr(sdfid, L3M_PRODTYPE, (VOIDP *) meta_l3m->prodtype)) < 0)
            return FAIL;
        meta_l3m->prodtype[count] = 0;
    } else
        meta_l3m->prodtype = 0;

    /**** read global attribute "Replacement Flag" */
    if ((getattrsz(sdfid, L3M_REPLACE, &nt, &count)) >= 0) {
        if ((meta_l3m->replace = (char *) malloc(sizeof (char) * count + 1))
                == NULL)
            return FAIL;
        if ((rdattr(sdfid, L3M_REPLACE, (VOIDP *) meta_l3m->replace)) < 0)
            return FAIL;
        meta_l3m->replace[count] = 0;
    } else
        meta_l3m->replace = 0;

    /**** read global attribute "Software ID" */
    if ((getattrsz(sdfid, L3M_SOFTID, &nt, &count)) >= 0) {
        if ((meta_l3m->softid = (char *) malloc(sizeof (char) * count + 1))
                == NULL)
            return FAIL;
        if ((rdattr(sdfid, L3M_SOFTID, (VOIDP *) meta_l3m->softid)) < 0)
            return FAIL;
        meta_l3m->softid[count] = 0;
    } else
        meta_l3m->softid = 0;

    /**** read global attribute "Software Name" */
    if ((getattrsz(sdfid, L3M_SOFTNM, &nt, &count)) >= 0) {
        if ((meta_l3m->softnm = (char *) malloc(sizeof (char) * count + 1))
                == NULL)
            return FAIL;
        if ((rdattr(sdfid, L3M_SOFTNM, (VOIDP *) meta_l3m->softnm)) < 0)
            return FAIL;
        meta_l3m->softnm[count] = 0;
    } else
        meta_l3m->softnm = 0;

    /**** read global attribute "Software Version" */
    if ((getattrsz(sdfid, L3M_SOFTVER, &nt, &count)) >= 0) {
        if ((meta_l3m->softver = (char *) malloc(sizeof (char) * count + 1))
                == NULL)
            return FAIL;
        if ((rdattr(sdfid, L3M_SOFTVER, (VOIDP *) meta_l3m->softver)) < 0)
            return FAIL;
        meta_l3m->softver[count] = 0;
    } else
        meta_l3m->softver = 0;

    /**** read global attribute "Processing Time" */
    if ((getattrsz(sdfid, L3M_PTIME, &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->ptime = (char *) malloc(sizeof (char) * count + 1)) == NULL)
        return FAIL;
    if ((rdattr(sdfid, L3M_PTIME, (VOIDP *) meta_l3m->ptime)) < 0)
        return FAIL;
    meta_l3m->ptime[count] = 0;

    /**** read global attribute "Input Files" */
    if ((getattrsz(sdfid, L3M_INFILES, &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->infiles = (char *) malloc(sizeof (char) * count + 1)) == NULL)
        return FAIL;
    if ((rdattr(sdfid, L3M_INFILES, (VOIDP *) meta_l3m->infiles)) < 0)
        return FAIL;
    meta_l3m->infiles[count] = 0;

    /**** read global attribute "Processing Control" */
    if ((getattrsz(sdfid, L3M_PROCCON, &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->proccon = (char *) malloc(sizeof (char) * count + 1)) == NULL)
        return FAIL;
    if ((rdattr(sdfid, L3M_PROCCON, (VOIDP *) meta_l3m->proccon)) < 0)
        return FAIL;
    meta_l3m->proccon[count] = 0;

    /**** read global attribute "Processing Log" */
    /*
     if ((getattrsz(sdfid, L3M_PROCLOG, &nt, &count)) < 0)  return FAIL;
     if ((meta_l3m->proclog = (char *) malloc (sizeof(char) * count + 1)) == NULL)
     return FAIL;
     if ((rdattr(sdfid, L3M_PROCLOG, (VOIDP *)meta_l3m->proclog))<0) return FAIL;
     meta_l3m->proclog[count] = 0;
     */

    /**** read global attribute "L2 Flag Usage" */
    /*
     if ((getattrsz(sdfid, L3M_FLAG_USE, &nt, &count)) < 0) return FAIL;
     if ((meta_l3m->flag_use =(char *)malloc(sizeof(char) * count + 1)) == NULL)
     return FAIL;
     if ((rdattr(sdfid, L3M_FLAG_USE, (VOIDP *)meta_l3m->flag_use))<0)
     return FAIL;
     */

    /**** read global attribute "L2 Engineering Quality Usage" */
    /*
     if ((rdattr(sdfid, L3M_ENG_Q_USE, (VOIDP *)&meta_l3m->eng_q_use)) < 0)
     return FAIL;
     */

    /**** read global attribute "time_coverage_start" */
    if ((getattrsz(sdfid, L3M_STIME, &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->stime = (char *) malloc(sizeof (char) * count + 1)) == NULL)
        return FAIL;
    if ((rdattr(sdfid, L3M_STIME, (VOIDP *) meta_l3m->stime)) < 0)
        return FAIL;
    meta_l3m->stime[count] = 0;

    /**** read global attribute "time_coverage_end" */
    if ((getattrsz(sdfid, L3M_ETIME, &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->etime = (char *) malloc(sizeof (char) * count + 1)) == NULL)
        return FAIL;
    if ((rdattr(sdfid, L3M_ETIME, (VOIDP *) meta_l3m->etime)) < 0)
        return FAIL;
    meta_l3m->etime[count] = 0;

    /**** read global attribute "Orbit" */
    //    if ((rdattr(sdfid, L3M_ORBIT, (VOIDP *) &meta_l3m->orbit)) < 0)
    //    return FAIL;

    /**** read global attribute "Start Orbit" */
    if ((rdattr(sdfid, L3M_SORBIT, (VOIDP *) & meta_l3m->start_orb)) < 0)
        return FAIL;

    /**** read global attribute "End Orbit" */
    if ((rdattr(sdfid, L3M_EORBIT, (VOIDP *) & meta_l3m->end_orb)) < 0)
        return FAIL;

    /**** read global attribute "Map Projection" */
    if ((getattrsz(sdfid, L3M_MAPPROJ, &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->mapproj = (char *) malloc(sizeof (char) * count + 1)) == NULL)
        return FAIL;
    if ((rdattr(sdfid, L3M_MAPPROJ, (VOIDP *) meta_l3m->mapproj)) < 0)
        return FAIL;
    meta_l3m->mapproj[count] = 0;

    /**** read global attribute "geospatial_lat_units" */
    if ((getattrsz(sdfid, L3M_LATUNITS, &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->lat_units = (char *) malloc(sizeof (char) * count + 1))
            == NULL)
        return FAIL;
    if ((rdattr(sdfid, L3M_LATUNITS, (VOIDP *) meta_l3m->lat_units)) < 0)
        return FAIL;
    meta_l3m->lat_units[count] = 0;

    /**** read global attribute "geospatial_lon_units" */
    if ((getattrsz(sdfid, L3M_LONUNITS, &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->lon_units = (char *) malloc(sizeof (char) * count + 1))
            == NULL)
        return FAIL;
    if ((rdattr(sdfid, L3M_LONUNITS, (VOIDP *) meta_l3m->lon_units)) < 0)
        return FAIL;
    meta_l3m->lon_units[count] = 0;

    /**** read global attribute "Nrothernmost Latitude" */
    if ((rdattr(sdfid, L3M_NLAT, (VOIDP *) & meta_l3m->nlat)) < 0)
        return FAIL;

    /**** read global attribute "Southernmost Latitude" */
    if ((rdattr(sdfid, L3M_SLAT, (VOIDP *) & meta_l3m->slat)) < 0)
        return FAIL;

    /**** read global attribute "Westernmost Longitude" */
    if ((rdattr(sdfid, L3M_WLON, (VOIDP *) & meta_l3m->wlon)) < 0)
        return FAIL;

    /**** read global attribute "Easternmost Latitude" */
    if ((rdattr(sdfid, L3M_ELON, (VOIDP *) & meta_l3m->elon)) < 0)
        return FAIL;

    /**** read global attribute "Latitude Step" */
    if ((rdattr(sdfid, L3M_LAT_STEP, (VOIDP *) & meta_l3m->lat_step)) < 0)
        return FAIL;

    /**** read global attribute "Longitude Step" */
    if ((rdattr(sdfid, L3M_LON_STEP, (VOIDP *) & meta_l3m->lon_step)) < 0)
        return FAIL;

    /**** read global attribute "Data Bins" */
    if ((rdattr(sdfid, L3M_DATABINS, (VOIDP *) & meta_l3m->nbins)) < 0)
        return FAIL;

    /**** read global attribute "Number of Lines" */
    if ((rdattr(sdfid, L3M_NROWS, (VOIDP *) & meta_l3m->nrows)) < 0)
        return FAIL;

    /**** read global attribute "Number of Columns" */
    if ((rdattr(sdfid, L3M_NCOLS, (VOIDP *) & meta_l3m->ncols)) < 0)
        return FAIL;

    /**** read global attribute "Parameter" */
    if ((getattrsz(sdfid, L3M_PARAMETER, &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->parameter = (char *) malloc(sizeof (char) * count + 1))
            == NULL)
        return FAIL;
    if ((rdattr(sdfid, L3M_PARAMETER, (VOIDP *) meta_l3m->parameter)) < 0)
        return FAIL;
    meta_l3m->parameter[count] = 0;

    /**** read global attribute "Measure" */
    if ((getattrsz(sdfid, L3M_MEASURE, &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->measure = (char *) malloc(sizeof (char) * count + 1)) == NULL)
        return FAIL;
    if ((rdattr(sdfid, L3M_MEASURE, (VOIDP *) meta_l3m->measure)) < 0)
        return FAIL;
    meta_l3m->measure[count] = 0;

    /**** read global attribute "Units" */
    if ((getattrsz(sdfid, L3M_UNITS, &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->units = (char *) malloc(sizeof (char) * count + 1)) == NULL)
        return FAIL;
    if ((rdattr(sdfid, L3M_UNITS, (VOIDP *) meta_l3m->units)) < 0)
        return FAIL;
    meta_l3m->units[count] = 0;

    /**** read global attribute "Scaling" */
    if ((getattrsz(sdfid, L3M_SCALING, &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->scaling_type = (char *) malloc(sizeof (char) * count + 1))
            == NULL)
        return FAIL;
    if ((rdattr(sdfid, L3M_SCALING, (VOIDP *) meta_l3m->scaling_type)) < 0)
        return FAIL;
    meta_l3m->scaling_type[count] = 0;

    /**** read global attribute "Scaling Equation" */
    if ((getattrsz(sdfid, L3M_SC_EQN, &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->scaling_eqn = (char *) malloc(sizeof (char) * count + 1))
            == NULL)
        return FAIL;
    if ((rdattr(sdfid, L3M_SC_EQN, (VOIDP *) meta_l3m->scaling_eqn)) < 0)
        return FAIL;
    meta_l3m->scaling_eqn[count] = 0;

#if 0
    /**** read global attribute "Slope" */
    if ((rdattr(sdfid, L3M_SLOPE, (VOIDP *) & meta_l3m->slope)) < 0) return FAIL;

    /**** read global attribute "Intercept" */
    if ((rdattr(sdfid, L3M_INTERCEPT, (VOIDP *) & meta_l3m->intercept)) < 0)
        return FAIL;
#endif

    /**** read global attribute "Data Minimum" */
    if ((rdattr(sdfid, L3M_MIN, (VOIDP *) & meta_l3m->data_min)) < 0)
        return FAIL;

    /**** read global attribute "Data Maximum" */
    if ((rdattr(sdfid, L3M_MAX, (VOIDP *) & meta_l3m->data_max)) < 0)
        return FAIL;

    /**** read global attribute "Scaled Data Minimum" */
    if ((rdattr(sdfid, "Suggested Image Scaling Minimum",
            (VOIDP *) & meta_l3m->scaled_data_min)) < 0)
        return FAIL;

    /**** read global attribute "Scaled Data Maximum" */
    if ((rdattr(sdfid, "Suggested Image Scaling Maximum",
            (VOIDP *) & meta_l3m->scaled_data_max)) < 0)
        return FAIL;

    /**** read global attribute "Scaled Data Type" */
    if ((getattrsz(sdfid, "Suggested Image Scaling Type", &nt, &count)) < 0)
        return FAIL;
    if ((meta_l3m->scaled_data_type = (char *) malloc(sizeof (char) * count + 1))
            == NULL)
        return FAIL;
    if ((rdattr(sdfid, "Suggested Image Scaling Type",
            (VOIDP *) meta_l3m->scaled_data_type)) < 0)
        return FAIL;
    meta_l3m->scaled_data_type[count] = 0;

#if 0
    /**** read global attribute "Slope (8 bit)" */
    meta_l3m->slope_8bit = 0.0;
    rdattr(sdfid, "B2F_slope", (VOIDP *) & meta_l3m->slope_8bit);

    /**** read global attribute "Intercept (8 bit)" */
    meta_l3m->intercept_8bit = 0.0;
    rdattr(sdfid, "B2F_intercept", (VOIDP *) & meta_l3m->intercept_8bit);
#endif

    return 0;
}

/*-----------------------------------------------------------------------------
 Function:  getattrsz

 Returns:   int32 (status)
 The return code is a negative value if any error occurs, otherwise,
 returns 0.

 Description:
 The function getattrsz passes the requested global attribute's
 number type (data type) and the number of values

 Parameters: (in calling order)
 Type      Name        I/O   Description
 ----      ----        ---   -----------
 int32     sdfid        I    ID req to access HDF SDS interface
 char  *   attr_name    I    attribute name
 int32 *   nt           O    HDF data type
 int32 *   count        O    number of values in the specified attribute

 Modification history:
 Programmer     Organization   Date      Description of change
 -------------- ------------   --------  ---------------------
 Lakshmi Kumar  Hughes STX     11/07/94  Original development

 ----------------------------------------------------------------------------*/

int32 getattrsz(int32 id, char *attr_name, int32 *nt, int32 *count) {
    int32 attrnum;
    char name[MAXVAL];

    attrnum = SDfindattr(id, attr_name);
    if ((SDattrinfo(id, attrnum, name, nt, count)) < 0) {
        /*
         printf("\ngetattrsz: Error reading attribute size of - %s\n", attr_name);
         */
        return FAIL;
    }
    return SUCCEED;
}

/*-----------------------------------------------------------------------------
 Function:  rdattr

 Returns:   int32 (status)
 The return code is a negative value if any error occurs, otherwise,
 returns 0.

 Description:
 The function rdattr reads the requested global attribute

 Parameters: (in calling order)
 Type      Name        I/O   Description
 ----      ----        ---   -----------
 int32     sdfid        I    ID req to access HDF SDS interface
 char  *   attr_name    I    attribute name
 void  *   buf         I/O   pointer to data buffer

 Modification history:
 Programmer     Organization   Date      Description of change
 -------------- ------------   --------  ---------------------
 Lakshmi Kumar  Hughes STX     11/07/94  Original development

 ----------------------------------------------------------------------------*/

int32 rdattr(int32 sdfid, char *attr_name, void *buf) {
    int32 attrnum;

    attrnum = SDfindattr(sdfid, attr_name);
    if ((SDreadattr(sdfid, attrnum, buf)) < 0) {
        /*
         printf("\nrdattr: Error reading attribute - %s\n", attr_name);
         */
        return FAIL;
    }
    return SUCCEED;
}
