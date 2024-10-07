/*-----------------------------------------------------------------------------
    File:  HDFroutines.c

    Contents:
        openHDF		-  Opens the given HDF file
        rdlatlon	-  Reads the geometry vdata and writes lat/lons
        get_clim_data	-  Reads climatology file for the req month 
        get_NRT_data	-  Reads NRT data from the given file
        get_refs	-  Traverses thro' the NRT file vgroup and outputs
                           the SDS ref number to access the requested data
        closeHDF	-  Closes the HDF file
  
    Other relevant files:
        anc.h         - various #defined constants for ancillary data, also
                        includes hdf.h
        ancproto.h    - prototypes for ancillary data input routines
        getanc.c      - a higher layer of ancillary data input functions
        dataintp.f    - interpolation routine

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      04/02/93    Original development
        Lakshmi Kumar    Hughes STX      12/07/93    Modified to incorporate
                                                     HDF3.3r1 changes
        Lakshmi Kumar	 Hughes STX	 08/11/95    Fixed couple of Unitializ-
                                                     ed memory read problems
        Lakshmi Kumar    Hughes STX      08/02/96    corrected the order of 
                                                     arguments passed in the 
                                                     function 'get_clim_data'
                                                     non-prototype defn. 
                                                     MAX has been redefined as 
                                                     MAXVAL to avoid compile 
                                                     time warnings, also 
                                                     removed unused variables
        Lakshmi Kumar  Hughes STX     	01/15/97     Reversed the order of
                                                     longitudes in lon_buf
        Lakshmi Kumar  Hughes STX     	03/31/97     changed to return the 
                                                     datatype of the attribute 
                                                     read.  Remove non-proto
                                                     function declarations.
-----------------------------------------------------------------------------*/

#include "netcdf.h"
#include <mfhdf.h>
#include "anc.h"
#include "ancproto.h"

/*-----------------------------------------------------------------------------
    Function:  openHDF      

    Returns:   int32 (status)
        The return code is a negative value if any error occurs.

    Description:
        The function openHDF opens the requested HDF file

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        char *    infile       I    input HDF file name  
        int32 *   sdfid        O    ID req to access HDF SDS interface 
        int32 *   fid          O    file ID  

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     12/09/93  Original development
 
----------------------------------------------------------------------------*/
int openHDF(char *infile, int32_t *sdfid, int32_t *fid) {
    int32 lsdfid, lfid;
    char *FUNC = "openHDF";

    if ((lsdfid = SDstart(infile, DFACC_RDONLY)) < 0) {
        sprintf(ERR_MSG, "%s: Cannot open file -- SDstart failed on %s",
                FUNC, infile);
        return FAIL;
    }

    if ((lfid = Hopen(infile, 0, DFACC_RDONLY)) < 0) {
        sprintf(ERR_MSG, "%s: Cannot open file -- Hopen failed on %s",
                FUNC, infile);
        return FAIL;
    }

    Vstart(lfid);

    *sdfid = lsdfid;
    *fid = lfid;

    return SUCCEED;

} /* startHDF */

/*-----------------------------------------------------------------------------
    Function:  rdancattr

    Returns:   int32 (status)
        The return code is a negative value if any error occurs, otherwise,
        returns the data type of the attribute read.

    Description:
        The function rdancattr reads the requested global attribute

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
    Lakshmi Kumar  Hughes STX	  03/31/97  changed to return the datatype of
                                            the attribute read
----------------------------------------------------------------------------*/

int32_t rdancattr(int32_t sdfid, char *attr_name, void *buf) {
    int32 attrnum, nt, count;
    char name[255];
    char *FUNC = "rdancattr";

    if ((attrnum = SDfindattr(sdfid, attr_name)) < 0) {
        sprintf(ERR_MSG, "%s: SDfindattr failed for %s", FUNC, attr_name);
        return FAIL;
    }

    if ((SDattrinfo(sdfid, attrnum, name, &nt, &count)) < 0) {
        sprintf(ERR_MSG, "%s: SDattrinfo failed for %s", FUNC, attr_name);
        return FAIL;
    }

    if ((SDreadattr(sdfid, attrnum, buf)) < 0) {
        sprintf(ERR_MSG, "%s: SDreadattr failed for %s", FUNC, attr_name);
        return FAIL;
    }

    return nt;
}

/*-----------------------------------------------------------------------------
    Function:  rdlatlon

    Returns:   int32 (status)
        The return code is a negative value if any error occurs.

    Description:
        The function rdlatlon reads an appropriate geometry vdata for the
        stepsize and corner coordinates.  Using the stepsize (vsize, hsize)
        and corner coordinates writes latitudes and longitudes in to the 
        lat/lon buffers.

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     fid          I    File ID 
        int32 *   dims         I    dimensions of the data
        float32 * lat_buf      O    latitude buffer
        float32 * lon_buf      O    longitude buffer

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     12/09/93  Original development
    Lakshmi Kumar  Hughes STX     08/11/95  Removed VSdetach call
    Lakshmi Kumar  Hughes STX     01/15/97  Reversed the order of longitudes
----------------------------------------------------------------------------*/
int rdlatlon(int32_t sdfid, int32_t *dims, float *lat_buf, float *lon_buf) {
    int32 i;
    float32 vsize, hsize, sw_lat, sw_lon;

    /**** read global attribute "Latitude Step" */
    if ((rdancattr(sdfid, VSIZE, (VOIDP *) & vsize)) < 0)
        return FAIL;

    /**** read global attribute "Longitude Step" */
    if ((rdancattr(sdfid, HSIZE, (VOIDP *) & hsize)) < 0)
        return FAIL;

    /**** read global attribute "SW Point Latitude" */
    if ((rdancattr(sdfid, SWPT_LAT, (VOIDP *) & sw_lat)) < 0)
        return FAIL;

    /**** read global attribute "SW Point Longitude" */
    if ((rdancattr(sdfid, SWPT_LON, (VOIDP *) & sw_lon)) < 0)
        return FAIL;

    lat_buf[dims[0] - 1] = sw_lat;
    for (i = dims[0] - 2; i >= 0; i--)
        lat_buf[i] = (lat_buf[i + 1] + vsize);

    for (lon_buf[0] = sw_lon, i = 1; i < dims[1]; i++)
        lon_buf[i] = (lon_buf[i - 1] + hsize);

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  closeHDF

    Returns:   none

    Description:
        The function closeHDF closes the given HDF file.

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     sdfid        I    ID req to access SDS interface
        int32     fid          I    File ID

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     12/09/93  Original development

----------------------------------------------------------------------------*/
int closeHDF(int32_t sdfid, int32_t fid) {
    int result;

    if ((result = Vend(fid)) < 0) perror("closing V file");
    if ((result = SDend(sdfid)) < 0) perror("closing SDS");
    if ((result = Hclose(fid)) < 0) perror("closing H file");

    return result;
}

/*-----------------------------------------------------------------------------
    Function:  get_clim_data

    Returns:   int32 (status)
        The return code is a negative value if any error occurs.

    Description:
        The function get_clim_data finds the right month vgroup,
        obtains right parameter SDS reference number and reads data in 
        to the output buffer.  Calls rd_latlon if required.

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     fid          I    File ID
        int32     sdfid        I    ID req to access SDS interface
        int16     parm_flag    I    flag indicating the parameter type
        int16     month        I    requested data month
        int32 *   dims         O    dimensions of the data
        float32 * lat_buf      O    latitude buffer
        float32 * lon_buf      O    longitude buffer
        void    * parm_buf     O    data buffer 

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     12/09/93  Original development
    Lakshmi Kumar  Hughes STX     08/11/95  Passing proper (MAX) vlaue to 
                                            Vgettagrefs 
    Lakshmi Kumar  Hughes STX     08/02/96  corrected the order of arguments
                                            passed in the non-prototype defn.
    W. Robinson, SAIC, 2 Dec 2013  remove the pr_dims and usage
----------------------------------------------------------------------------*/
int get_clim_data(int32_t fid, int32_t sdfid, int16_t parm_flag,
        int16_t month, int32_t *dims, float *lat_buf,
        float *lon_buf, void *parm_buf) {
    int32 i, done = 0;
    int32 vref, vid, index, sdsid, nrefs;
    int32 start[3] = {0, 0, 0}, taglist[MAXVAL], reflist[MAXVAL];
    int32 rank, nattrs, numbertype;
    char name[MAXVAL], *FUNC = "get_clim_data";

    if ((vref = Vfind(fid, clim_vgps[month - 1])) <= 0) {
        sprintf(ERR_MSG, "%s: Vfind failed to find vgroup %s",
                FUNC, clim_vgps[month - 1]);
        return FAIL;
    }

    if ((vid = Vattach(fid, vref, "r")) < 0) {
        sprintf(ERR_MSG, "%s: Vattach failed for vgroup %s",
                FUNC, clim_vgps[month - 1]);
        return FAIL;
    }

    if ((Vgettagrefs(vid, taglist, reflist, MAXVAL)) < 0) {
        sprintf(ERR_MSG, "%s: Vgettagrefs failed for vgroup %s",
                FUNC, clim_vgps[month - 1]);
        return FAIL;
    }

    if ((nrefs = Vntagrefs(vid)) < 0)
        return FAIL;

    /*** get sdsid for the requested parameter sds */
    for (i = 0; i < nrefs && !done; i++) {
        if (taglist[i] == DFTAG_NDG) {
            if ((index = SDreftoindex(sdfid, reflist[i])) < 0)
                return FAIL;
            if ((sdsid = SDselect(sdfid, index)) < 0)
                return FAIL;
            if ((SDgetinfo(sdsid, name, &rank, dims, &numbertype, &nattrs)) < 0)
                return FAIL;
            if (strcmp(name, clim_datasets[parm_flag]) == 0)
                done = 1;
        }
    }

    if (!done) {
        sprintf(ERR_MSG,
                "%s: requested dataset %s not found", FUNC, clim_datasets[parm_flag]);
        return FAIL;
    }

    if ((rdlatlon(sdfid, dims, lat_buf, lon_buf)) < 0)
        return FAIL;

    if ((SDreaddata(sdsid, start, NULL, dims, parm_buf)) < 0)
        return FAIL;

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  get_NRT_data

    Returns:   int32 (status)
        The return code is a negative value if any error occurs.

    Description:
        The function get_NRT_data finds the right parameter vgroup and
        calls get_refs to get the requested SDS reference number and
        dimensions of the data.  Calls rd_latlon if required.
        

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     fid          I    File ID
        int32     sdfid        I    ID req to access SDS interface
        char *    anc_cor_file I    ancillary correction file, if NULL,
                                    apply no correction t(to ozone)
        int16     parm_flag    I    flag indicating the parameter type
        int32 *   dims         O    dimensions of the data
        float32 * lat_buf      O    latitude buffer
        float32 * lon_buf      O    longitude buffer
        void    * parm_buf     O    data buffer
        int8    * qc_buf       O    qc data buffer

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     12/09/93  Original development
    Lakshmi Kumar  HITC		  05/25/95  Modified code to read the re-
                                            structured data files (see product
                                            specs 2.7)
    W. Robinson, SAIC, 2 Dec 2013  remove the pr_dims and usage
    W. Robinson, SAIC, 17 Jan 2014  add capability to read and apply the 
         correction to N7, EP TOMS and AURAOMI
----------------------------------------------------------------------------*/
int get_NRT_data(int32_t fid, int32_t sdfid, char *anc_cor_file, int16_t parm_flag,
        int32_t *dims, float *lat_buf, float *lon_buf,
        void *parm_buf, int8_t *qc_buf) {

    int32 sdsid, index, numbertype;
    int32 rank, start[3] = {0, 0, 0}, nattrs;
    char name[132], *FUNC = "get_NRT_data", dsrc[100];
    int32_t oz_typ;
    int16 st_yr, st_day;
    int corr_oz(char *, int16 *, char *, int32 *, int16, int16, int);

    /*  read   data SDS */
    if ((index = SDnametoindex(sdfid, data_sdsnames[parm_flag])) < 0) {
        sprintf(ERR_MSG, "%s: SDnametoindex failed for %s",
                FUNC, data_sdsnames[parm_flag]);
        return FAIL;
    }

    if ((sdsid = SDselect(sdfid, index)) < 0) {
        sprintf(ERR_MSG, "%s: SDselect failed for %s",
                FUNC, data_sdsnames[parm_flag]);
        return FAIL;
    }

    if ((SDgetinfo(sdsid, name, &rank, dims, &numbertype, &nattrs)) < 0)
        return FAIL;

    if ((SDreaddata(sdsid, start, NULL, dims, (VOIDP) parm_buf)) < 0) {
        sprintf(ERR_MSG, "%s: SDreaddata faild for sds - %s", FUNC, name);
        return FAIL;
    }

    SDendaccess(sdsid);


    /*  read  QC  SDS */
    if ((index = SDnametoindex(sdfid, QC_sdsnames[parm_flag])) < 0) {
        sprintf(ERR_MSG, "%s: SDnametoindex failed for %s ",
                FUNC, QC_sdsnames[parm_flag]);
        return FAIL;
    }

    if ((sdsid = SDselect(sdfid, index)) < 0)
        return FAIL;

    if ((SDgetinfo(sdsid, name, &rank, dims, &numbertype, &nattrs)) < 0)
        return FAIL;

    if ((SDreaddata(sdsid, start, NULL, dims, (VOIDP) qc_buf)) < 0)
        return FAIL;

    SDendaccess(sdsid);

    if ((rdlatlon(sdfid, dims, lat_buf, lon_buf)) < 0)
        return FAIL;

    /*
     *  for the ozone, if it has a source of the TOMS types, see about correcting
     */
    if (parm_flag == 4) {
        if (rdancattr(sdfid, "Data Source", (void *) dsrc) != DFNT_CHAR) {
            printf("%s, %d E: Could not read 'Data Source' from ozone file\n",
                    __FILE__, __LINE__);
            return FAIL;
        }
        oz_typ = -1;
        if (strcmp(dsrc, "N7TOMS") == 0)
            oz_typ = 0;
        else if (strcmp(dsrc, "EPTOMS") == 0)
            oz_typ = 1;
        else if (strcmp(dsrc, "AURAOMI") == 0)
            oz_typ = 2;

        if (oz_typ != -1) {
            /*
             *  for applicable ozones, see if the oz_corr_file is not null
             *  use clo_get... to get the oz_corr_file
             *  initially just set here
             *
             *  at this time, we'd have to get the oz_corr_file into msl12_input.c and 
             *  input structure and then ferry it from setanc - get_ancillary - 
             *  read_NRT - get_NRT_data to use OR use clo_get_String if list was 
             *  global, later, probably do the 1st option
             *  for now, just hard code it
             */
            if (strcmp(anc_cor_file, "") != 0) {
                /*
                 *  now, read the Start Year, Start Day and apply the correction
                 */
                if ((rdancattr(sdfid, "Start Year", (void *) &st_yr) != DFNT_INT16)
                        || (rdancattr(sdfid, "Start Day", (void *) &st_day) != DFNT_INT16)) {
                    printf(
                            "%s, %d E: Could not read 'Start Year' or 'Start Day' from ozone file\n",
                            __FILE__, __LINE__);
                    return FAIL;
                }
                if (corr_oz(anc_cor_file, parm_buf, (char*) qc_buf, dims, st_yr,
                        st_day, oz_typ) != 0)
                    return FAIL;
            }
        }
    }

    return SUCCEED;
}

/*******************************************************************

   corr_oz

   purpose: read out the proper ozone correction to apply to the ozone
     data array

   Returns type: int - 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            oz_corr_file     I      name of ozone correction file
      int16 *           oz_dat          I/O     array of ozone data
      char *            qc_dat           I      array of qc data, 20 indicates 
                                                  bad data
      int32 *           dims             I      array of x, y size of oz_dat
      int16             yr, doy          I      year of the ozone array
      int               oz_typ           I      ozone type: 0 N7, 1 EP, 2 AURA

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC, 17 Jan 2014    original development

 *******************************************************************/
int corr_oz(char *oz_corr_file, int16 *oz_dat, char *qc_dat, int32 *dims,
        int16 yr, int16 doy, int oz_typ) {
    int ncid, gid, status, var_id;
    char *oz_inst_grp[] = {"N7TOMS", "EPTOMS", "AURAOMI", "Extrap_AURAOMI"};
    int16 corr_time_st[4], corr_time_en[4], mon, day;
    int32_t yrs[2], mons[2], npix_oz, nlin_oz, ilin, ipix, ixmon, tot_mon;
    int32_t ixlat;
    float t_wt1, t_wt2, l_wt1, l_wt2, corr_arr[2 * 36], lat, slat, del_lat;
    float pos, corr_fac;
    size_t start[2], count[2];
    /*
     *  open the oz_corr_file
     */
    if (nc_open(oz_corr_file, 0, &ncid) != NC_NOERR) {
        printf("%s, %d E: Unable to open ozone correction file: %s\n",
                __FILE__, __LINE__, oz_corr_file);
        return -1;
    }
    /*
     *  go to the group of corrections for that instruent
     */
    if (nc_inq_grp_ncid(ncid, oz_inst_grp[oz_typ], &gid) != NC_NOERR) {
        printf("%s, %d E: Unable to set group: %s for ozone correction file: %s\n",
                __FILE__, __LINE__, oz_inst_grp[oz_typ], oz_corr_file);
        return -1;
    }
    printf("%s, %d I: Correcting ozone using file:\n", __FILE__, __LINE__);
    printf("          %s\n", oz_corr_file);
    printf("          ozone inst: %s\n", oz_inst_grp[oz_typ]);
    /*
     *  read the time array for the instrument
     */
    if ((status = nc_inq_varid(gid, "start_time", &var_id)) != NC_NOERR) {
        printf("%s, %d E: Unable to read 'start_time' from correction file: %s\n",
                __FILE__, __LINE__, oz_corr_file);
        return -1;
    }
    if ((status = nc_get_var_short(gid, var_id, corr_time_st)) != NC_NOERR) {
        printf("%s, %d E: Unable to read 'start_time' from correction file: %s\n",
                __FILE__, __LINE__, oz_corr_file);
        return -1;
    }
    ;
    if ((status = nc_inq_varid(gid, "end_time", &var_id)) != NC_NOERR) {
        printf("%s, %d E: Unable to read 'end_time' from correction file: %s\n",
                __FILE__, __LINE__, oz_corr_file);
        return -1;
    }
    if ((status = nc_get_var_short(gid, var_id, corr_time_en)) != NC_NOERR) {
        printf("%s, %d E: Unable to read 'end_time' from correction file: %s\n",
                __FILE__, __LINE__, oz_corr_file);
        return -1;
    }
    /*
     *  determine which months to read for the ozone time
     */
    gregor(doy, yr, &mon, &day);
    if (day <= 15) {
        yrs[1] = yr;
        mons[1] = mon;
        if (mon == 1) {
            yrs[0] = yr - 1;
            mons[0] = 12;
        } else {
            yrs[0] = yr;
            mons[0] = mon - 1;
        }
        t_wt1 = (15 - day) / 30.;
        t_wt2 = (15 + day) / 30.;
    } else {
        yrs[0] = yr;
        mons[0] = mon;
        if (mon == 12) {
            yrs[1] = yr + 1;
            mons[1] = 1;
        } else {
            yrs[1] = yr;
            mons[1] = mon + 1;
        }
        t_wt1 = (45 - day) / 30.;
        t_wt2 = (day - 15) / 30.;
    }
    /*
     *  compute the position of the start of the 2 months to interpolate on
     *  and the total # months
     */
    ixmon = (yrs[0] - corr_time_st[0]) * 12 +
            (mons[0] - corr_time_st[1]);
    tot_mon = (corr_time_en[0] - corr_time_st[0]) * 12 +
            (corr_time_en[1] - corr_time_st[1] + 1);
    /*
     *  check for out of bounds and use the end time alone unless beyond the
     *  AURA end time - then, use the extrapolated month set
     */
    if (ixmon < 0) {
        ixmon = 0;
        t_wt1 = 1.;
        t_wt2 = 0.;
    } else if (ixmon >= (tot_mon - 1)) {
        if (oz_typ == 2) /*  it is AURA, we use extrapolation array  */ {
            /*
             *  the extrapolated AURA is set up to start in January and repeat
             *  January so that a pair of months will always be side-by-side.
             *  The time weights are still good from above, so we only need
             *  to set to the Extrapolated_AURA group , set the ixmon and the
             *  following code to read will work fine.
             */
            oz_typ = 3;
            if (nc_inq_grp_ncid(ncid, oz_inst_grp[oz_typ], &gid) != NC_NOERR) {
                printf(
                        "%s, %d E: Unable to set group: %s for ozone correction file: %s\n",
                        __FILE__, __LINE__, oz_inst_grp[oz_typ], oz_corr_file);
                return -1;
            }
            printf("%s, %d I: Extrapolated ozone correction used, name: %s\n",
                    __FILE__, __LINE__, oz_inst_grp[oz_typ]);
            ixmon = mons[0] - 1;
        } else {
            ixmon = tot_mon - 2;
            t_wt1 = 0.;
            t_wt2 = 1.;
        }
    }
    /*
     * read in the 2 month by 36 lat portion of the proper dataset 
     */
    if ((status = nc_inq_varid(gid, "correction_array", &var_id))
            != NC_NOERR) {
        printf("%s, %d E: Unable to read 'correction_array' from correction file: %s\n",
                __FILE__, __LINE__, oz_corr_file);
        return -1;
    }
    start[0] = 0;
    start[1] = ixmon;
    count[0] = 36;
    count[1] = 2;

    if ((status =
            nc_get_vars_float(gid, var_id, start, count, NULL, corr_arr))
            != NC_NOERR) {
        printf("%s, %d E: Unable to read 'correction_array' from correction file: %s\n",
                __FILE__, __LINE__, oz_corr_file);
        return -1;
    }
    npix_oz = dims[1];
    nlin_oz = dims[0];
    /*
     *  close the correction file
     */
    nc_close(ncid);
    /*
     *  NOTE!!! that the ozone array is stored so that the start line is at 89.5
     *  latitude, not -89.5 (SW Latitude in dataset).  So, we need to consider 
     *  that when applying the correction (which is not reversed in this silly way)
     *
     *  The ozone grid is grid box centered and starts at 90 - box siz / 2
     */
    /*  if storage started at -90...
     del_lat = 180. / nlin_oz;
     slat = -90. + del_lat / 2.;
     */
    del_lat = -180. / nlin_oz;
    slat = 90. + del_lat / 2.;
    /*
     *  loop through the lines and pixels of the ozone and correct all the non-bad
     *  values with the interpolated correction value
     *
     */
    for (ilin = 0; ilin < nlin_oz; ilin++) {
        lat = slat + ilin * del_lat;
        pos = (lat + 87.5) / 5.; /* pos is in units of the correction grid */
        ixlat = (int) (pos + 1) - 1;
        if (ixlat < 0) {
            ixlat = 0;
            l_wt1 = 1.;
            l_wt2 = 0.;
        } else if (ixlat > 34) {
            ixlat = 34;
            l_wt1 = 0;
            l_wt2 = 1.;
        } else {
            l_wt1 = (float) (ixlat + 1) - pos;
            l_wt2 = pos - (float) ixlat;
        }
        /*
         *  get the correction factor for the time and latitude
         */
        corr_fac = t_wt1 * (l_wt1 * *(corr_arr + 2 * ixlat) +
                l_wt2 * *(corr_arr + 2 * (ixlat + 1))) +
                t_wt2 * (l_wt1 * *(corr_arr + 1 + 2 * ixlat) +
                l_wt2 * *(corr_arr + 1 + 2 * (ixlat + 1)));
        /*
         *  apply this to all the pixels, the 0.5 is to get closest correction 
         *  to the integer values of ozone (I'd change this in re-write too)
         */
        for (ipix = 0; ipix < npix_oz; ipix++) {
            if (*(qc_dat + ipix + npix_oz * ilin) != 20)
                *(oz_dat + ipix + npix_oz * ilin) += (corr_fac + 0.5);
        }
    }
    /*  and end with ozone modified  */
    return 0;
}

/*-----------------------------------------------------------------------------
    Function:  get_refs

    Returns:   int32 (status)
        The return code is a negative value if any error occurs.

    Description:
        The function get_refs attaches the right parameter vgroup and
        traverses through the vgroup to find the parameter SDS and the
        QC SDS. Outputs the reference numbers of the parameter and the
        QC SDSs.

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     fid          I    File ID
        int32     sdfid        I    ID req to access SDS interface
        int32     vid          I    Parameter Vgroup ID
        char *    parm_label   I    parameter SDS label
        char *    QC_label     I    QC SDS label         
        int32 *   parm_sdsid   O    parameter SDS ID      
        int32 *   QC_sdsid     O    QC SDS ID       
        int32 *   dims         O    dimensions of the data

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     12/09/93  Original development

----------------------------------------------------------------------------*/
int get_refs(int32_t fid, int32_t sdfid, int32_t vid, char *parm_label,
        char *QC_label, int32_t *geom_id, int32_t *parm_sdsid, int32_t *QC_sdsid,
        int32_t *dims) {
    int32 i, j, ret = 0, vg, natts = MAXVAL, nrefs, tagarray[MAXVAL];
    int32 refarray[MAXVAL];
    int32 index, sdsid, rank, numbertype, nattrs, numSDSs;
    int32 sdrefs[25];
    char name[256];

    if ((vg = Vattach(fid, vid, "r")) < 0)
        return FAIL;

    Vgettagrefs(vg, tagarray, refarray, natts);
    nrefs = Vntagrefs(vg);
    for (i = 0, j = 0; i < nrefs; i++) {
        if (tagarray[i] == DFTAG_NDG || tagarray[i] == DFTAG_SD) {
            sdrefs[j++] = refarray[i];
        } else
            if (tagarray[i] == DFTAG_VH || tagarray[i] == DFTAG_VS)
            *geom_id = refarray[i];
    }

    numSDSs = j;

    for (i = 0; i < numSDSs && ret >= 0; i++) {
        index = SDreftoindex(sdfid, sdrefs[i]);
        sdsid = SDselect(sdfid, index);
        if ((SDgetinfo(sdsid, name, &rank, dims, &numbertype, &nattrs)) < 0)
            return FAIL;
        if (strcmp(name, parm_label) == 0) {
            *parm_sdsid = sdsid;
        } else
            if (strcmp(name, QC_label) == 0)
            *QC_sdsid = sdsid;
    }

    Vdetach(vg);

    return SUCCEED;
}
