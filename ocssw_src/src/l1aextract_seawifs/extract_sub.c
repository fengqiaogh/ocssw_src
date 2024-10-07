/*----------------------------------------------------------------------------
    File:  regen.c

    Contents:
        regen                   -  Creates Level1a or level2 data file with
                                        the requested subset of the given
                                        input data file

        set_subsc_params        -  Error checks the given input subscene
                                        parameters

        rdattr                  -  Reads requested attribute

        getset_gattrs           -  Reads all the global attributes from the
                                        given input file and writes them to the
                                        output file

        dupHDF                  -  Reads given input file and creates skeleton
                                        output file with the input file
                                        structure

        set_sds                 -  Reads and writes skeleton SDS in the output
                                        file with all the attribute and
                                        dimension information

        write_data              -  Reads given SDS data, calls subsampling
                                        routine to subsample and writes the
                                        subsampled data to the output SDS

        rdslice                 -  Reads requested slice of data from the
                                        given named dataset

        subsample               -  Calls proper subsampling routine depending
                                        upon whether the subsampling is along
                                        one dimension or two dimensions

        update_sds              -  Updates the specified chunk of data with the
                                        given data of the named SDS

        subsamp_rec             -  Subsamples data along last dimension

        subsamp_2D              -  Subsamples 2D data along two dimensions

        alloc_nav_buffs         -  Allocates data buffers

        get_navdata             -  Reads sensor metrics from the given input
                                        data file

        alloc_geonav_buffs      -  Allocates data buffers for geolocation data

        get_geodata             -  Calls geonav_ for each subsampled scanline
                                        and outputs the geolocation buffers

        geonav_                 -  Computes geolocation information for the
                                        requested pixels of the given scanline

        cdata_                  -  Should be called before geonav_

        write_coords            -  Updates scene corordinates info in the
                                        output file

        get_tiltdata            -  Reads tilt data from the input data file

        set_tiltdata            -  Updates the tilt datasets in the output file
                                        for the requested subset of data

        set_l1adata             -  Processes level 1a data, if the input file
                                        is level 1a

        subsamp_l1adata         -  Subsamples the given 3D data along last
                                        two dimensions ('C' order)

        set_flags               -  Accumulates level 2 flags, if the given
                                        input file is level 2

        set_globalattrs         -  Updates some of the global attributes that
                                        depends on the subset

        free_geonav_buffs       -  frees buffer space

        free_nav_buffs          -  frees buffer space

    Notes:

    Other relevant files:
        regen.h         - various #defined constants for level 1a and level 2
                                data, also #includes hdf.h and mfhdf.h
        regen_proto.h   - prototypes of all the functions defined


    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      11/22/93    Original development
        Lakshmi Kumar    Hughes STX      12/14/95    Fixed a bug in set_l1adata
                                                     & a bug in subsamp_l1adata
        Lakshmi Kumar    Hughes STX      02/21/96    Fixed the above bugfix
                                                     which had some problem
        Lakshmi Kumar    Hughes STX      07/29/96    The variable sun has been
                                                     renamed as sunref and
                                                     PROTOTYPE defs have been
                                                     removed
        Lakshmi Kumar    Hughes STX      02/18/97    Changed variable type of
                                                     flags from int16 to int32
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "hdf.h"
#include "mfhdf.h"
#include "regen.h"
#include "regen_attr.h"
#include "regen_proto.h"
#include "navigation.h"
#include <timeutils.h>

#define  CHUNK_SZ 500
#define  L1ADATA "l1a_data"

char ERR_MSG[255];

/*-----------------------------------------------------------------------------
    Function:  regen

    Returns:   int32 (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        The function regen regenerates a given level 1a or level2 file.
        Allows subsampling and subsetting of the data.

    Parameters: (in calling order)
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        char      *infile       I   input file name
        int32     *ssamp        I   starting pixel
        int32     *esamp        I   ending pixel
        int32     *srec         I   starting scan line
        int32     *erec         I   ending scan line
        int32     px_sub        I   pixel subsampling rate
        int32     sc_sub        I   scan subsampling rate
        char      *parm_list    I   list of parameter names that needs to be
                                    subsampled (Optional)
        char      *outfile      I   output file name

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     03/31/95  Original development
----------------------------------------------------------------------------*/
int32
regen(char *infile, int32 *ssamp, int32 *esamp, int32 *srec, int32 *erec,
        int32 px_sub, int32 sc_sub, char *parm_list, char *outfile) {
    int32 subsc_samps, subsc_recs;
    int32 ifid, ofid, isdfid, osdfid;
    int32 xsub = px_sub, ysub = sc_sub;
    int32 nvgps, nsamp, nrec;
    int32 max_samp_used, max_rec_used;
    int32 pix_start, pix_sub;
    int16 l1aflag;
    int32 flags[16];
    char dtype[MAXVAL], title[MAXVAL];
    NavType nav_rec;
    GeoType geo_rec;
    tilt_Type oldtilt_params, newtilt_params;
    FilemetricsType fm_rec;

    /*** open the given input HDF file */
    if ((isdfid = SDstart(infile, DFACC_RDONLY)) < 0) {
        sprintf(ERR_MSG, "regen: File open unsuccessful");
        return FAIL;
    }
    ifid = Hopen(infile, DFACC_READ, 0);
    Vstart(ifid);

    /*** open the given output HDF file */
    osdfid = SDstart(outfile, DFACC_CREATE);
    ofid = Hopen(outfile, DFACC_RDWR, 0);
    Vstart(ofid);

    /*** find out how many vgroups have been written */
    if ((nvgps = Hnumber(ifid, DFTAG_VG)) < 0) {
        sprintf(ERR_MSG, "\n No vgps found\n");
        return FAIL;
    }

    /***  read global attribute "Data Type"  */
    if ((rdattr(isdfid, DTYPE, dtype)) < 0)
        return FAIL;

    /***  read global attribute "Title"  */
    if ((rdattr(isdfid, "Title", title)) < 0)
        return FAIL;

    /***  read nsamp and nrec */
    if ((rdattr(isdfid, NSAMP, &nsamp)) < 0)
        return FAIL;

    if ((rdattr(isdfid, NREC, &nrec)) < 0)
        return FAIL;

    /*** check validity of input paramters and set subscene parameters */
    set_subsc_params(nsamp, nrec, xsub, ysub, ssamp, esamp, srec, erec,
            &subsc_samps, &subsc_recs, &max_samp_used,
            &max_rec_used);

    /***  set global attributes */
    if ((getset_gattrs(isdfid, osdfid)) < 0)
        return FAIL;

    /*** create vgroups and its objects in the output file */
    if ((dupHDF(ifid, ofid, isdfid, osdfid, ssamp, esamp, srec, erec,
            xsub, ysub, subsc_samps, subsc_recs, nvgps, parm_list)) < 0)
        return FAIL;

    /*** read sensor metrics */
    if ((alloc_nav_buffs(nrec, &nav_rec)) < 0)
        return FAIL;

    if ((get_navdata(isdfid, nrec, &nav_rec)) < 0)
        return FAIL;

    /*** get & update geonavigation datasets and related global attributes */
    if ((alloc_geonav_buffs(subsc_recs, &geo_rec)) < 0)
        return FAIL;

    if ((rdattr(isdfid, PIX_START, &pix_start)) < 0)
        return FAIL;
    if ((rdattr(isdfid, PIX_SUB, &pix_sub)) < 0)
        return FAIL;

    if ((get_geodata(pix_start, pix_sub, *srec, *ssamp, max_rec_used, subsc_samps,
            ysub, xsub, dtype, &nav_rec, &geo_rec)) < 0)
        return FAIL;

    if ((write_coords(osdfid, subsc_recs, &geo_rec)) < 0)
        return FAIL;

    /*** update tilt datasets and related global attributes */
    if ((get_tiltdata(isdfid, &oldtilt_params)) < 0) return FAIL;

    if ((set_tiltdata(osdfid, *srec, max_rec_used, ysub, geo_rec.slat,
            geo_rec.slon, geo_rec.elat, geo_rec.elon,
            &oldtilt_params, &newtilt_params)) < 0)
        return FAIL;

    /*** set pct_flags */
    if ((strcmp(title, "SeaWiFS Level-2 Data") == 0) ||
            (strcmp(title, "SeaWiFS Level-2 Q/C Data") == 0)) {
        if ((set_flags(osdfid, subsc_recs, subsc_samps, flags)) < 0)
            return FAIL;
        l1aflag = 0;
    } else {
        if ((strcmp(title, "SeaWiFS Level-1A Data")) == 0)
            l1aflag = 1;
        else
            return FAIL;
    }
    /*** update l1a_data and related global attributes */
    if (l1aflag) {
        if ((set_l1adata(isdfid, osdfid, nrec, nsamp, subsc_recs, subsc_samps,
                *srec, *ssamp, *erec, *esamp, xsub, ysub, &fm_rec)) < 0)
            return FAIL;
    }

    /*** update required global attributes */
    if ((set_globalattrs(outfile, osdfid, l1aflag, subsc_recs, subsc_samps,
            *ssamp, xsub, flags, &geo_rec, &fm_rec)) < 0)
        return FAIL;

    /*** free all the allocated buffers */
    free_geonav_buffs(&geo_rec);
    free_nav_buffs(&nav_rec);

    Vend(ifid);
    SDend(isdfid);
    Hclose(ifid);

    Vend(ofid);
    SDend(osdfid);
    Hclose(ofid);

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  set_nparams

    Returns:   None

    Description:
        The function set_nparams checks the given input parameters (ssamp,
        esamp, srec, erec).  If any of these parameters are in out of range
        then the rtn sets proper  value to that parameter.  Computes
        subscene dimensions and sets some output paramters

    Parameters: (in calling order)
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int32     nsamp         I   number of pixels in the given input file
        int32     nrec          I   number of records/scans in the input file
        int32     *ssamp        I   starting pixel
        int32     *esamp        I   ending pixel
        int32     *srec         I   starting scan line
        int32     *erec         I   ending scan line
        int32     *subsc_recs   O   subscene records
        int32     *subsc_samps  O   subscene pixels
        int32     *max_samp_used O  max pixel used in subscene image
        int32     *max_rec_used O   max rec/scanline used in subscece image
        char      *proc_log     O   processing log

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     04/03/95  Original development
----------------------------------------------------------------------------*/
void
set_subsc_params(int32 nsamp, int32 nrec, int32 xsub, int32 ysub,
        int32 *ssamp, int32 *esamp, int32 *srec, int32 *erec,
        int32 *subsc_samps, int32 *subsc_recs, int32 *max_samp_used,
        int32 *max_rec_used) {

    (*ssamp)--;
    (*esamp)--;
    (*srec)--;
    (*erec)--;

    /*** TODO: make handling of boundary conditions consistent
         with other extractors; this seems too lenient. ***/

    if (*ssamp < 0 || *ssamp >= nsamp) {
        *ssamp = 0;
        printf("\n***** The given start sample is in error");
        printf("\n\tresetting ssamp to zero\n\n");
    }

    if (*srec < 0 || *srec >= nrec) {
        *srec = 0;
        printf("\n***** The given start record is in error");
        printf("\n\tresetting srec to zero \n\n");
    }

    if (*esamp >= nsamp || *esamp < *ssamp) {
        *esamp = nsamp - 1;
        printf("\n****** The given end sample is in error ");
        printf("\n\tresetting esamp to nsamp = %d\n\n", nsamp);
    }

    if (*erec >= nrec || *erec < *srec) {
        *erec = nrec - 1;
        printf("\n****** The given end record is in error ");
        printf("\n\tresetting erec to nrec = %d \n\n", nrec);
    }

    /*** compute subscene image dimensions */

    *subsc_recs = ((*erec + 1) - (*srec + 1)) / ysub + 1;
    *subsc_samps = ((*esamp + 1) - (*ssamp + 1)) / xsub + 1;


    /*** compute the maximum record and sample numbers used to make the subsampled
            browse output */

    *max_rec_used = (*subsc_recs - 1) * ysub + *srec;
    *max_samp_used = (*subsc_samps - 1) * xsub + *ssamp;
}

/*-----------------------------------------------------------------------------
    Function:  getset_gattrs

    Returns:   int32 (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        The function getset_gattrs reads all the global attributes from
        the given input file and writes them to the output file.

    Parameters: (in calling order)
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int32     isdfid        I    SDinterface ID of the input file
        int32     osdfid        I    SDinterface ID of the output file

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     04/04/95  Original development
----------------------------------------------------------------------------*/
int32 getset_gattrs(int32 isdfid, int32 osdfid) {
    int32 i, numtype[MAXATTRS], cnt[MAXATTRS];
    int32 ndatasets, nglobal_attrs;
    char char_data[MAXVAL], attr_name[MAXATTRS][132];

    char *charbuf;
    int32 i32buf[30];
    int16 i16buf[30];
    int8 i8buf[30];
    uint8 ui8buf[30];
    float32 f32buf[30];
    float64 f64buf[30];

    if ((SDfileinfo(isdfid, &ndatasets, &nglobal_attrs)) < 0)
        return FAIL;

    for (i = 0; i < nglobal_attrs; i++) {
        if ((SDattrinfo(isdfid, i, &attr_name[i][0], &numtype[i], &cnt[i])) < 0)
            return FAIL;

        switch (numtype[i]) {
        case DFNT_CHAR:
            if (cnt[i] > MAXVAL)
                charbuf = (char *) calloc(cnt[i] + 1, sizeof (char));
            else
                charbuf = &char_data[0];
            if ((SDreadattr(isdfid, i, charbuf)) < 0) {
                sprintf(ERR_MSG,
                        "get_globalattrs: SDreadattr not successful for attribute %s", attr_name[i]);
                return FAIL;
            }
            if ((SDsetattr(osdfid, attr_name[i], numtype[i], cnt[i],
                    (VOIDP) charbuf)) < 0) {
                sprintf(ERR_MSG,
                        "get_globalattrs: SDsetattr not successful for attribute %s",
                        attr_name[i]);
                return FAIL;
            }
            if (cnt[i] > MAXVAL)
                free(charbuf);
            break;
        case DFNT_INT32:
            if ((SDreadattr(isdfid, i, i32buf)) < 0) {
                sprintf(ERR_MSG,
                        "get_globalattrs: SDreadattr not successful for attribute %s",
                        attr_name[i]);
                return FAIL;
            }
            if ((SDsetattr(osdfid, attr_name[i], numtype[i], cnt[i],
                    (VOIDP) i32buf)) < 0) {
                sprintf(ERR_MSG,
                        "get_globalattrs: SDsetattr not successful for attribute %s",
                        attr_name[i]);
                return FAIL;
            }
            break;
        case DFNT_INT16:
            if ((SDreadattr(isdfid, i, i16buf)) < 0) {
                sprintf(ERR_MSG,
                        "get_globalattrs:SDreadattr not successful for attribute %s", attr_name[i]);
                return FAIL;
            }
            if ((SDsetattr(osdfid, attr_name[i], numtype[i], cnt[i],
                    (VOIDP) i16buf)) < 0) {
                sprintf(ERR_MSG,
                        "get_globalattrs: SDsetattr not successful for attribute %s", attr_name[i]);

                return FAIL;
            }
            break;
        case DFNT_INT8:
            if ((SDreadattr(isdfid, i, i8buf)) < 0) {
                sprintf(ERR_MSG,
                        "get_globalattrs: SDreadattr not successful for attribute %s",
                        attr_name[i]);
                return FAIL;
            }
            if ((SDsetattr(osdfid, attr_name[i], numtype[i], cnt[i],
                    (VOIDP) i8buf)) < 0) {
                sprintf(ERR_MSG,
                        "get_globalattrs: SDsetattr not successful for attribute %s",
                        attr_name[i]);
                return FAIL;
            }
            break;
        case DFNT_UINT8:
            if ((SDreadattr(isdfid, i, ui8buf)) < 0) {
                sprintf(ERR_MSG,
                        "get_globalattrs:SDreadattr unsuccessful for attribute %s",
                        attr_name[i]);
                return FAIL;
            }
            if ((SDsetattr(osdfid, attr_name[i], numtype[i], cnt[i],
                    (VOIDP) ui8buf)) < 0) {
                sprintf(ERR_MSG,
                        "get_globalattrs:SDsetattr unsuccessful for attribute %s",
                        attr_name[i]);
                return FAIL;
            }
            break;
        case DFNT_FLOAT32:
            if ((SDreadattr(isdfid, i, f32buf)) < 0) {
                sprintf(ERR_MSG,
                        "get_globalattrs: SDreadattr unsuccessful for attr %s",
                        attr_name[i]);
                return FAIL;
            }
            if ((SDsetattr(osdfid, attr_name[i], numtype[i], cnt[i],
                    (VOIDP) f32buf)) < 0) {
                sprintf(ERR_MSG,
                        "get_globalattrs:SDsetattr unsuccessful for attr %s",
                        attr_name[i]);
                return FAIL;
            }
            break;
        case DFNT_FLOAT64:
            if ((SDreadattr(isdfid, i, f64buf)) < 0) {
                sprintf(ERR_MSG,
                        "get_globalattrs: SDreadattr unsuccessful for attr %s",
                        attr_name[i]);
                return FAIL;
            }
            if ((SDsetattr(osdfid, attr_name[i], numtype[i], cnt[i],
                    (VOIDP) f64buf)) < 0) {
                sprintf(ERR_MSG,
                        "get_globalattrs:SDsetattr unsuccessful for attr %s",
                        attr_name[i]);
                return FAIL;
            }
            break;
        default:
            sprintf(ERR_MSG, "getset_gattrs: Unkown data type - %d\n",
                    numtype[i]);
            return FAIL;
        }
    }
    /*
       printf("\n No. of SDSs = %d and global attrs = %d\n", ndatasets,
                            nglobal_attrs);
     */
    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function: dupHDF

    Returns:   int32 (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        The function dupHDF copies the structure of the given input
        data file to the given output file

    Parameters: (in calling order)
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int32     ifid          I   input file ID
        int32     ofid          I   output file ID
        int32     isdfid        I   SDinterface ID of input file
        int32     osdfid        I   SDinterface ID of output file
        int32     *ssamp        I   starting pixel
        int32     *esamp        I   ending pixel
        int32     *srec         I   starting scan line
        int32     *erec         I   ending scan line
        int32     xsub          I   pixel subsampling rate
        int32     ysub          I   scan subsampling rate
        int32     subsc_samps   I   no. of req. subscene samples
        int32     subsc_recs    I   no. of req. subscene records
        int32     nvgps         I   no. of vgroups present in the input file
        char      *parm_list    I   list of parameter names that needs to be
        char      *parm_list    I   list of parameter names that needs to be
                                    processed

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     03/31/95  Original development
----------------------------------------------------------------------------*/
int32 dupHDF(int32 ifid, int32 ofid, int32 isdfid, int32 osdfid, int32 *ssamp,
        int32 *esamp, int32 *srec, int32 *erec, int32 xsub,
        int32 ysub, int32 subsc_samps, int32 subsc_recs, int32 nvgps,
        char *parm_list) {
    int32 i, j, n_entries;
    int32 vref, ivid, ovid;
    int32 isdsid, osdsid;
    int32 tag_list[MAXVAL], ref_list[MAXVAL];
    char vgname[VGNAMELENMAX];
    char vgclass[VGNAMELENMAX];

    vref = -1;
    for (i = 0; i < nvgps; i++) {
        if ((vref = Vgetid(ifid, vref)) < 0)
            break;
        ivid = Vattach(ifid, vref, "r");

        if ((Vinquire(ivid, &n_entries, vgname)) < 0) {
            printf("\nError: Vinquire failed for vid = %d\n", ivid);
            exit(-1);
        }
        Vgetclass(ivid, vgclass);

        if ((strcmp(vgclass, "Dim0.0") != 0) &&
                (strcmp(vgclass, "Var0.0") != 0) &&
                (strcmp(vgclass, "CDF0.0") != 0)) {
            ovid = Vattach(ofid, -1, "w");
            Vsetclass(ovid, vgclass);
            Vsetname(ovid, vgname);
            if ((Vgettagrefs(ivid, tag_list, ref_list, MAXVAL)) < 0)
                exit(-1);
            for (j = 0; j < n_entries; j++) {
                switch (tag_list[j]) {
                case DFTAG_NDG:
                case DFTAG_SDG:
                    if ((set_sds(isdfid, osdfid, vgname, ovid, tag_list[j],
                            ref_list[j], subsc_samps, subsc_recs, parm_list,
                            &isdsid, &osdsid)) < 0) return FAIL;
                    if ((write_data(isdfid, osdfid, isdsid, osdsid, *ssamp, *srec,
                            *esamp, *erec, xsub, ysub)) < 0) return FAIL;

                    break;
                default:
                    printf("\n not processed -- tag = %d", tag_list[j]);
                    break;
                }
            }
            Vdetach(ovid);
        }

        Vdetach(ivid);
    }

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  set_sds

    Returns:   int32 (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        The function set_sds reads all the required information of the
        requested dataset (SDS) from the given input file and creates an
        SDS with no data in the output file with the same same, attributes
        and dimensions

    Parameters: (in calling order)
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int32     isdfid        I    SDinterface ID of the input file
        int32     osdfid        I    SDinterface ID of the output file
        char      *vgname       I    name of the Vgroup for which the
                                        referenced SDS belongs
        int32     ovid          I    output file Vgroup ID for which the new
                                        SDS should be linked
        int32     tag           I    input SDS tag
        int32     ref           I    input SDS reference number
        int32     subsc_samps   I    no. of output samples/pixels
        int32     subsc_recs    I    no. of output records/scans
        char      *parm_list    I    list of l2 parameters that should be
                                        processed
        int32     *in_sdsid     O    SDS id of the input file
        int32     *out_sdsid    O    SDS id of the output file

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     03/31/95  Original development
----------------------------------------------------------------------------*/
int32 set_sds(int32 isdfid, int32 osdfid, char *vgname, int32 ovid,
        int32 tag, int32 ref, int32 subsc_samps, int32 subsc_recs,
        char *parm_list, int32 *in_sdsid, int32 *out_sdsid) {
    int32 i, index, ret;
    int32 isdsid, osdsid;
    int32 sdsref, rank, num_type, nattrs, attr_nt;
    int32 count, idim_id, odim_id, dim_count, dim_num, dim_attrs;
    int32 olddims[3] = {0, 0, 0}, newdims[3] = {0, 0, 0};
    int32 start[3] = {0, 0, 0}, edge[3] = {0, 0, 0};
    int32 found = 0; /* flags */
    char dim_name[3][MAXVAL], buf[100];
    char *parm, attr_name[MAXVAL], sdsname[MAXVAL], str[MAXVAL];

    if ((index = SDreftoindex(isdfid, ref)) < 0)
        return FAIL;
    isdsid = SDselect(isdfid, index);
    if ((ret = SDgetinfo(isdsid, sdsname, &rank, olddims, &num_type, &nattrs)) < 0)
        return FAIL;

    if (((strcmp(vgname, "Geophysical Data")) == 0) &&
            ((strcmp(parm_list, "ALL")) != 0) &&
            ((strcmp(sdsname, "l2_flags")) != 0)) {

        strcpy(str, parm_list);
        parm = strtok(str, " ");
        while (parm != NULL) {
            if ((strcmp(parm, sdsname)) == 0) {
                found = 1;
                break;
            }
            parm = strtok(NULL, " ");
        }

        if (!found) {
            SDendaccess(isdsid);
            return SUCCEED;
        }
    }

    for (i = 0; i < rank; i++) {
        newdims[i] = olddims[i];
        idim_id = SDgetdimid(isdsid, i);
        if ((ret = SDdiminfo(idim_id, dim_name[i], &dim_count, &dim_num,
                &dim_attrs)) < 0) return FAIL;
        /*
            here "rec" and "nsamp" are kept to be consistent with files created
            by old (frank's l1agen)         G. Fu  12/11/97
         */
        if ((strcmp(dim_name[i], "Number of Scan Lines") == 0 ||
                (strcmp(dim_name[i], "rec")) == 0))
            newdims[i] = subsc_recs;
        if ((strcmp(dim_name[i], "Pixels per Scan Line") == 0 ||
                (strcmp(dim_name[i], "nsamp")) == 0))
            newdims[i] = subsc_samps;
    }

    /*** create SDS with given name and dimensions */
    if ((osdsid = SDcreate(osdfid, sdsname, num_type, rank, newdims)) < 0)
        return FAIL;

    if ((ret = SDwritedata(osdsid, start, NULL, edge, (VOIDP) buf)) < 0)
        return FAIL;

    for (i = 0; i < rank; i++) {
        odim_id = SDgetdimid(osdsid, i);
        ret = SDsetdimname(odim_id, dim_name[i]);
    }

    SDendaccess(isdsid);
    SDendaccess(osdsid);

    /*** read attributes of the given sds from the input SeaWiFS file and
            write them to the output SeaWiFS file */
    for (ret = 0, i = 0; ret >= 0; i++) {
        ret = SDattrinfo(isdsid, i, attr_name, &attr_nt, &count);
        if (ret >= 0) {
            if ((rdattr(isdsid, attr_name, (VOIDP) str)) < 0) {
                printf("Error: Unsuccessful reading attribute:%s", attr_name);
                return FAIL;
            }
            if ((SDsetattr(osdsid, attr_name, attr_nt, count, (VOIDP) str)) < 0) {
                printf("Error: Unsuccessful writing attribute %s", attr_name);
                return FAIL;
            }
        }
    }

    /*** link this sds to the given vgruoup ID */
    if ((sdsref = SDidtoref(osdsid)) < 0) return FAIL;
    if ((Vaddtagref(ovid, tag, sdsref)) < 0) return FAIL;

    *in_sdsid = isdsid;
    *out_sdsid = osdsid;

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  write_data

    Returns:   int32 (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        The function write_data reads input SDS, calls subsampling routines
        and writes it to the output SDS

    Parameters: (in calling order)
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int32     isdfid        I   SDinterface ID of input file
        int32     osdfid        I   SDinterface ID of output file
        int32     isdsid        I   input SDS ID
        int32     osdsid        I   output SDS ID
        int32     ssamp         I   starting pixel
        int32     srec          I   starting scan line
        int32     esamp         I   ending pixel
        int32     erec          I   ending scan line
        int32     xsub          I   pixel subsampling rate
        int32     ysub          I   scan subsampling rate

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     03/31/95  Original development
----------------------------------------------------------------------------*/
int32 write_data(int32 isdfid, int32 osdfid, int32 isdsid, int32 osdsid,
        int32 ssamp, int32 srec, int32 esamp, int32 erec,
        int32 xsub, int32 ysub) {

    int32 i;
    int32 rank, num_type, nattrs;
    int32 idims[3] = {0, 0, 0}, odims[3] = {0, 0, 0};
    int32 istart[3] = {0, 0, 0}, iedge[3] = {1, 1, 1};
    int32 ostart[3] = {0, 0, 0}, oedge[3] = {1, 1, 1};
    char sdsname[MAXVAL];
    int32 *ii32p, *oi32p;
    int16 *ii16p, *oi16p;
    int8 *ii8p, *oi8p;
    uint8 *iui8p, *oui8p;
    float32 *if32p, *of32p;
    float64 *if64p, *of64p;
    char *icharp, *ocharp;

    if ((SDgetinfo(isdsid, sdsname, &rank, idims, &num_type, &nattrs)) < 0)
        return FAIL;

    if ((SDgetinfo(osdsid, sdsname, &rank, odims, &num_type, &nattrs)) < 0)
        return FAIL;

    if (strcmp(sdsname, L1ADATA) == 0)
        return SUCCEED;

    for (i = 0; i < rank; i++) {
        iedge[i] = idims[i];
        oedge[i] = odims[i];
    }

    switch (num_type) {
    case DFNT_CHAR:
        if ((icharp = (char *) calloc(iedge[0] * iedge[1] * iedge[2],
                sizeof (char))) == NULL)
            return MALLOC_ERR;
        if ((ocharp = (char *) calloc(oedge[0] * oedge[1] * oedge[2],
                sizeof (char))) == NULL)
            return MALLOC_ERR;

        if ((rdslice(isdfid, sdsname, istart, idims, (VOIDP) icharp)) < 0)
            return FAIL;

        if ((idims[0] != odims[0]) || (idims[1] != odims[1])) {
            subsample(rank, idims, odims, num_type, srec, erec,
                    ssamp, esamp, xsub, ysub, icharp, ocharp);

            if ((update_sds(osdfid, sdsname, ostart, odims, (VOIDP) ocharp)) < 0)
                return FAIL;
        } else
            if ((update_sds(osdfid, sdsname, ostart, odims, (VOIDP) icharp)) < 0)
            return FAIL;

        free(icharp);
        free(ocharp);

        break;

    case DFNT_INT32:
        if ((ii32p = (int32 *) calloc(iedge[0] * iedge[1] * iedge[2],
                sizeof (int32))) == NULL)
            return MALLOC_ERR;
        if ((oi32p = (int32 *) calloc(oedge[0] * oedge[1] * oedge[2],
                sizeof (int32))) == NULL)
            return MALLOC_ERR;

        if ((rdslice(isdfid, sdsname, istart, idims, (VOIDP) ii32p)) < 0)
            return FAIL;

        if ((idims[0] != odims[0]) || (idims[1] != odims[1])) {
            subsample(rank, idims, odims, num_type, srec, erec,
                    ssamp, esamp, xsub, ysub, ii32p, oi32p);

            if ((update_sds(osdfid, sdsname, ostart, odims, (VOIDP) oi32p)) < 0)
                return FAIL;
        } else
            if ((update_sds(osdfid, sdsname, ostart, odims, (VOIDP) ii32p)) < 0)
            return FAIL;

        free(oi32p);
        free(ii32p);

        break;

    case DFNT_INT16:
        if ((ii16p = (int16 *) calloc(iedge[0] * iedge[1] * iedge[2],
                sizeof (int16))) == NULL)
            return MALLOC_ERR;
        if ((oi16p = (int16 *) calloc(oedge[0] * oedge[1] * oedge[2],
                sizeof (int16))) == NULL)
            return MALLOC_ERR;

        if ((rdslice(isdfid, sdsname, istart, idims, (VOIDP) ii16p)) < 0)
            return FAIL;

        if ((idims[0] != odims[0]) || (idims[1] != odims[1])) {
            subsample(rank, idims, odims, num_type, srec, erec,
                    ssamp, esamp, xsub, ysub, ii16p, oi16p);
            if ((update_sds(osdfid, sdsname, ostart, odims, (VOIDP) oi16p)) < 0)
                return FAIL;
        } else
            if ((update_sds(osdfid, sdsname, ostart, odims, (VOIDP) ii16p)) < 0)
            return FAIL;

        free(ii16p);
        free(oi16p);
        break;

    case DFNT_INT8:
        if ((ii8p = (int8 *) calloc(iedge[0] * iedge[1] * iedge[2],
                sizeof (int8))) == NULL)
            return MALLOC_ERR;
        if ((oi8p = (int8 *) calloc(oedge[0] * oedge[1] * oedge[2],
                sizeof (int8))) == NULL)
            return MALLOC_ERR;

        if ((rdslice(isdfid, sdsname, istart, idims, (VOIDP) ii8p)) < 0)
            return FAIL;

        if ((idims[0] != odims[0]) || (idims[1] != odims[1])) {
            subsample(rank, idims, odims, num_type, srec, erec,
                    ssamp, esamp, xsub, ysub, ii8p, oi8p);

            if ((update_sds(osdfid, sdsname, ostart, odims, (VOIDP) oi8p)) < 0)
                return FAIL;
        } else
            if ((update_sds(osdfid, sdsname, ostart, odims, (VOIDP) ii8p)) < 0)
            return FAIL;
        free(ii8p);
        free(oi8p);
        break;

    case DFNT_UINT8:
        if ((iui8p = (uint8 *) calloc(iedge[0] * iedge[1] * iedge[2],
                sizeof (uint8))) == NULL)
            return MALLOC_ERR;
        if ((oui8p = (uint8 *) calloc(oedge[0] * oedge[1] * oedge[2],
                sizeof (uint8))) == NULL)
            return MALLOC_ERR;

        if ((rdslice(isdfid, sdsname, istart, idims, (VOIDP) iui8p)) < 0)
            return FAIL;

        if ((idims[0] != odims[0]) || (idims[1] != odims[1])) {
            subsample(rank, idims, odims, num_type, srec, erec,
                    ssamp, esamp, xsub, ysub, iui8p, oui8p);

            if ((update_sds(osdfid, sdsname, ostart, odims, (VOIDP) oui8p)) < 0)
                return FAIL;
        } else
            if ((update_sds(osdfid, sdsname, ostart, odims, (VOIDP) iui8p)) < 0)
            return FAIL;

        free(iui8p);
        free(oui8p);
        break;

    case DFNT_FLOAT32:
        if ((if32p = (float32 *) calloc(iedge[0] * iedge[1] * iedge[2],
                sizeof (float32))) == NULL)
            return MALLOC_ERR;
        if ((of32p = (float32 *) calloc(oedge[0] * oedge[1] * oedge[2],
                sizeof (float32))) == NULL)
            return MALLOC_ERR;

        if ((rdslice(isdfid, sdsname, istart, idims, (VOIDP) if32p)) < 0)
            return FAIL;

        if ((idims[0] != odims[0]) || (idims[1] != odims[1])) {
            subsample(rank, idims, odims, num_type, srec, erec,
                    ssamp, esamp, xsub, ysub, if32p, of32p);

            if ((update_sds(osdfid, sdsname, ostart, odims, (VOIDP) of32p)) < 0)
                return FAIL;
        } else
            if ((update_sds(osdfid, sdsname, ostart, odims, (VOIDP) if32p)) < 0)
            return FAIL;

        free(if32p);
        free(of32p);

        break;

    case DFNT_FLOAT64:
        if ((if64p = (float64 *) calloc(iedge[0] * iedge[1] * iedge[2],
                sizeof (float64))) == NULL)
            return MALLOC_ERR;
        if ((of64p = (float64 *) calloc(oedge[0] * oedge[1] * oedge[2],
                sizeof (float64))) == NULL)
            return MALLOC_ERR;

        if ((rdslice(isdfid, sdsname, istart, idims, (VOIDP) if64p)) < 0)
            return FAIL;

        if ((idims[0] != odims[0]) || (idims[1] != odims[1])) {
            subsample(rank, idims, odims, num_type, srec, erec,
                    ssamp, esamp, xsub, ysub, if64p, of64p);

            if ((update_sds(osdfid, sdsname, ostart, odims, (VOIDP) of64p)) < 0)
                return FAIL;
        } else
            if ((update_sds(osdfid, sdsname, ostart, odims, (VOIDP) if64p)) < 0)
            return FAIL;

        free(if64p);
        free(of64p);

        break;
    default:
        sprintf(ERR_MSG,
                "write_data: sds %s not processed - unkonwn data type - %d",
                sdsname, num_type);
        return FAIL;
    }

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  subsample

    Returns:   int32 (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        The function subsample calls proper subsampling routine depending
        upon whether the subsampling is along one dimension or two dimensions

    Parameters: (in calling order)
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int32     rank          I   no. of dimensions
        int32     *idims        I   input data dimensions
        int32     *odims        I   output data dimensions
        int32     num_type      I   data type
        int32     srec          I   starting scan line
        int32     erec          I   ending scan line
        int32     ssamp         I   starting pixel
        int32     esamp         I   ending pixel
        int32     xsub          I   pixel subsampling rate
        int32     ysub          I   scan subsampling rate
        void      *ibuf         I   input data buffer
        void      *obuf         O   output data buffer

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     03/31/95  Original development
----------------------------------------------------------------------------*/
int32 subsample(int32 rank, int32 *idims, int32 *odims, int32 num_type,
        int32 srec, int32 erec, int32 ssamp, int32 esamp,
        int32 xsub, int32 ysub, void *ibuf, void *obuf) {

    int32 subsc_recs;
    int32 maxrec;


    /*** compute subscene image dimensions */

    subsc_recs = ((erec + 1) - (srec + 1)) / ysub + 1;

    /*** compute the last record no. that should be processed*/

    maxrec = (subsc_recs - 1) * ysub + srec;

    if (idims[1] == odims[1]) {
        if ((subsamp_rec(srec, maxrec, rank, odims, ysub, num_type, ibuf, obuf))
                < 0) return FAIL;
    } else
        subsamp_2D(ssamp, srec, esamp, erec, idims[1], xsub, ysub, num_type,
            ibuf, obuf);

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
        void  *   buf          O    data buffer

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     11/07/94  Original development
----------------------------------------------------------------------------*/
int32 rdattr(int32 sdfid, char *attr_name, void *buf) {
    int32 attrnum;

    attrnum = SDfindattr(sdfid, attr_name);
    if ((SDreadattr(sdfid, attrnum, buf)) < 0) {
        sprintf(ERR_MSG, "rdattr: SDreadattr not successful for attribute %s",
                attr_name);
        return FAIL;
    }
    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  rdslice

    Returns:   int32 (sdsid)
        The return code is a negative value if any error occurs, otherwise,
        returns sdsid.

    Description:
        The function rdslice reads requested slice of data from the
        given named dataset

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     sdfid         I    ID req to access HDF SDS interface
        char      *name         I   SDS name
        int32     *start        I   start data dimension
        int32     *edge         I   no. of values to be read
        void      *buf          O   SDS data buffer

    NOTE:

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     11/02/94  Original development
----------------------------------------------------------------------------*/
int32 rdslice(int32 sdfid, char *name, int32 *start, int32 *edge, void *buf) {
    int32 index, sdsid, rank, num_type, nattrs;
    char sdsname[MAXVAL];
    /*
      clock_t  val, val2;
     */
    if ((index = SDnametoindex(sdfid, name)) < 0) {
        sprintf(ERR_MSG, "rdslice: Cannot locate sds \"%s\" ", name);
        return FAIL;
    }
    if ((sdsid = SDselect(sdfid, index)) < 0) {
        sprintf(ERR_MSG, "rdslice: SDselect failed for sds \"%s\" ", name);
        return FAIL;
    }

    if (edge[0] == 0 && edge[1] == 0 && edge[2] == 0)
        if ((SDgetinfo(sdsid, sdsname, &rank, edge, &num_type, &nattrs)) < 0)
            return FAIL;
    /*
      val = clock();
     */
    if ((SDreaddata(sdsid, start, NULL, edge, buf)) < 0) {
        sprintf(ERR_MSG,
                "rdslice: SDreaddata error while reading \"%s\" ", name);
        return FAIL;
    }
    /*
      val2 = clock();
      printf("\ntime took for %s = %d\n", name, val2-val);
     */
    SDendaccess(sdsid);
    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  update_sds

    Returns:   int32 (sdsid)
        The return code is a negative value if any error occurs, otherwise,
        returns sdsid.

    Description:
        The function update_sds writes/overwrites the given data of the
        named sds of the given file.  The start and edge specifies the
        starting location and no. of values that should be overwritten

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     sdfid        I    SDinterface ID of the input data file
        char      *name         I   SDS name
        int32     *start        I   start data dimensions
        int32     *edge         I   no. of values to be read alon each dim.
        void      *buf          O   data buffer

    NOTE:

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     11/02/94  Original development
----------------------------------------------------------------------------*/
int32 update_sds(int32 sdfid, char *name, int32 *start, int32 *edge, void *buf) {
    int32 sdsid, index;
    int32 rank, num_type, nattrs;
    char sdsname[MAXVAL];

    if ((index = SDnametoindex(sdfid, name)) < 0) {
        sprintf(ERR_MSG, "update_sds: Cannot locate sds \"%s\" ", name);
        return FAIL;
    }
    if ((sdsid = SDselect(sdfid, index)) < 0) {
        sprintf(ERR_MSG, "update_sds: SDselect failed for sds \"%s\" ", name);
        return FAIL;
    }

    if (edge[0] == 0 && edge[1] == 0 && edge[2] == 0)
        if ((SDgetinfo(sdsid, sdsname, &rank, edge, &num_type, &nattrs)) < 0)
            return FAIL;

    if ((SDwritedata(sdsid, start, NULL, edge, buf)) < 0) {
        sprintf(ERR_MSG,
                "update_sds: SDwritedata error while writing \"%s\" ", name);
        return FAIL;
    }
    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  alloc_geonav_buffs

    Returns:   int32 (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        The function alloc_geonav_buffs allocates buffers for storing
        starting, ending and center scene coordinates of each scan line rec

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     nrec          I   number of records/scanlines
        GeoType  *geo_rec       O   structrue to hold allocated buffers

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     04/03/95  Original development
----------------------------------------------------------------------------*/
int32 alloc_geonav_buffs(int32 nrec, GeoType *geo_rec) {

    /*** allocate space for geodata information */
    if ((geo_rec->slat = (float32 *) calloc(nrec, sizeof (float32))) == NULL)
        return MALLOC_ERR;

    if ((geo_rec->slon = (float32 *) calloc(nrec, sizeof (float32))) == NULL)
        return MALLOC_ERR;

    if ((geo_rec->elat = (float32 *) calloc(nrec, sizeof (float32))) == NULL)
        return MALLOC_ERR;

    if ((geo_rec->elon = (float32 *) calloc(nrec, sizeof (float32))) == NULL)
        return MALLOC_ERR;

    if ((geo_rec->clat = (float32 *) calloc(nrec, sizeof (float32))) == NULL)
        return MALLOC_ERR;

    if ((geo_rec->clon = (float32 *) calloc(nrec, sizeof (float32))) == NULL)
        return MALLOC_ERR;

    if ((geo_rec->csol_z = (float32 *) calloc(nrec, sizeof (float32))) == NULL)
        return MALLOC_ERR;

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  free_geonav_buffs

    Returns:   None

    Description:
        The function free_geonav_buffs frees all the scene coordinate buffers

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        GeoType  *geo_rec       I   structrue containing scene coordinates
                                        buffers

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     04/03/95  Original development
----------------------------------------------------------------------------*/
void free_geonav_buffs(GeoType *geo_rec) {
    free(geo_rec->slat);
    free(geo_rec->slon);
    free(geo_rec->elat);
    free(geo_rec->elon);
    free(geo_rec->clat);
    free(geo_rec->clon);
    free(geo_rec->csol_z);
}

/*-----------------------------------------------------------------------------
    Function:  get_geodata

    Returns:   int32 (sdsid)
        The return code is a negative value if any error occurs, otherwise,
        returns sdsid.

    Description:
        The function get_geodata calls cdata_ and geonav_ (Fortran functions)
        to get scene coordinates information for the starting, center and
        end pixels of each output (subscene) scanline

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     pix_start     I   LAC pixel start number from the input file
        int32     pix_sub       I   LAC pixel subsampling from the input file
        int32     srec          I   start record number
        int32     ssamp         I   start sample number
        int32     max_rec_used  I   last record number to be processed
        int32     subsc_samps   I   output scene samples
        int32     ysub          I   subsampling rate along Y-axis
        int32     xsub          I   subsampling rate along X-axis
        char     *dtype         I   data type (GAC, LAC, etc..)
        NavType  *nav_rec       I   sensor metrics data structure
        GeoType  *geo_rec       O   scene coordinates data structure

    NOTE:

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     11/02/94  Original development
----------------------------------------------------------------------------*/
int32 get_geodata(int32 pix_start, int32 pix_sub, int32 srec, int32 ssamp,
        int32 max_rec_used, int32 subsc_samps, int32 ysub, int32 xsub,
        char *dtype, NavType *nav_rec, GeoType *geo_rec) {

    int32 i, j, k, row, col, rec, scan, st_px, px_sub, npix;
    float32 xlat[1285], xlon[1285], solz[1285], sola[1285];
    float32 senz[1285], sena[1285];
    float32 pos[3], smat[3][3], coef[6], sunref[3];

    if (strcmp(dtype, "GAC") == 0) {
        st_px = pix_start + ssamp * 4; /* updated by G. Fu 12/5/97 */
        px_sub = ((subsc_samps - 1) * xsub * pix_sub) / 2;
    } else {
        st_px = pix_start + ssamp;
        px_sub = ((subsc_samps - 1) * xsub * pix_sub) / 2;
    }

    npix = 3;

    /* call cdata.f to initialize global FORTRAN common block data  */
    cdata_();

    for (scan = 0, rec = srec; rec <= max_rec_used; rec += ysub, scan++) {
        for (j = rec * 3, row = 0; row < 3; row++, j++) {
            pos[row] = nav_rec->orb_vec[j];
            sunref[row] = nav_rec->sun_ref[j];
            for (k = j * 3, col = 0; col < 3; col++, k++)
                smat[row][col] = nav_rec->sen_mat[k];
        }

        for (i = rec * 6, row = 0; row < 6; row++, i++)
            coef[row] = nav_rec->scan_ell[i];

        geonav_(pos, (float*) smat, coef, sunref, &st_px, &px_sub, &npix, xlat, xlon,
                solz, sola, senz, sena);

        geo_rec->slat[scan] = xlat[0];
        geo_rec->slon[scan] = xlon[0];
        geo_rec->clat[scan] = xlat[1];
        geo_rec->clon[scan] = xlon[1];
        geo_rec->elat[scan] = xlat[2];
        geo_rec->elon[scan] = xlon[2];
        geo_rec->csol_z[scan] = solz[1];
    }

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  write_coords

    Returns:   int32 (sdsid)
        The return code is a negative value if any error occurs, otherwise,
        returns sdsid.

    Description:
        The function write_coords

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     sdfid        I    SDinterface ID of output file
        int32     nrec          I   number of records
        GeoType  *geo_rec       I   scene coordinates data structure

    NOTE:

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     11/02/94  Original development
----------------------------------------------------------------------------*/
int32 write_coords(int32 sdfid, int32 nrec, GeoType *geo_rec) {
    int32 start[3] = {0, 0, 0}, edge[3] = {0, 0, 0};

    edge[0] = nrec;

    /*** call getdims only once for getting scene coords dimensions */

    if ((update_sds(sdfid, SLAT, start, edge, (VOIDP) geo_rec->slat)) < 0)
        return FAIL;
    if ((update_sds(sdfid, SLON, start, edge, (VOIDP) geo_rec->slon)) < 0)
        return FAIL;
    if ((update_sds(sdfid, CLAT, start, edge, (VOIDP) geo_rec->clat)) < 0)
        return FAIL;
    if ((update_sds(sdfid, CLON, start, edge, (VOIDP) geo_rec->clon)) < 0)
        return FAIL;
    if ((update_sds(sdfid, ELAT, start, edge, (VOIDP) geo_rec->elat)) < 0)
        return FAIL;
    if ((update_sds(sdfid, ELON, start, edge, (VOIDP) geo_rec->elon)) < 0)
        return FAIL;
    if ((update_sds(sdfid, CSOLZ, start, edge, (VOIDP) geo_rec->csol_z)) < 0)
        return FAIL;
    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  get_tiltdata

    Returns:   int32 (sdsid)
        The return code is a negative value if any error occurs, otherwise,
        returns sdsid.

    Description:
        The function get_tiltdata reads tilt datasets from the given input
        file

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     sdfid        I    ID req to access HDF SDS interface
        tilt_Type *tiltrec      I   structure containing tilt data buffers

    NOTE:

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     11/02/94  Original development
----------------------------------------------------------------------------*/
int32 get_tiltdata(int32 sdfid, tilt_Type *tiltrec) {
    int32 start[3] = {0, 0, 0}, edge[3];

    edge[0] = edge[1] = edge[2] = 0;
    if ((rdslice(sdfid, NTILTS, start, edge, (VOIDP) & tiltrec->ntilts)) < 0)
        return FAIL;

    edge[0] = edge[1] = edge[2] = 0;
    if ((rdslice(sdfid, TILTFLAGS, start, edge, (VOIDP) tiltrec->tilt_flags)) < 0)
        return FAIL;

    edge[0] = edge[1] = edge[2] = 0;
    if ((rdslice(sdfid, TILTRANGES, start, edge, (VOIDP) tiltrec->tilt_ranges)) < 0)
        return FAIL;

    edge[0] = edge[1] = edge[2] = 0;
    if ((rdslice(sdfid, TILTLATS, start, edge, (VOIDP) tiltrec->tilt_lats)) < 0)
        return FAIL;

    edge[0] = edge[1] = edge[2] = 0;
    if ((rdslice(sdfid, TILTLONS, start, edge, (VOIDP) tiltrec->tilt_lons)) < 0)
        return FAIL;

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  set_tiltdata

    Returns:   int32 (sdsid)
        The return code is a negative value if any error occurs, otherwise,
        returns sdsid.

    Description:
        The function set_tiltdata

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     sdfid        I    ID req to access HDF SDS interface
        int32     srec          I   start record
        int32     erec          I   ending record/scanline
        int32     ysub          I   scan subsampling rate
        float32   *slatrec      I   starting latitude points
        float32   *slonrec      I   starting longitude points
        float32   *elatrec      I   ending latitude points
        float32   *elonrec      I   ending longitude points
        tilt_Type *old_tiltrec  I   structure containing input tilt datasets
        tilt_Type *new_tiltrec  O   structure containing output tilt datasets

    NOTE:

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     04/03/95  Original development
----------------------------------------------------------------------------*/
int32 set_tiltdata(int32 sdfid, int32 srec, int32 erec, int32 ysub,
        float32 *slatrec, float32 *slonrec, float32 *elatrec,
        float32 *elonrec, tilt_Type *old_tiltrec,
        tilt_Type *new_tiltrec) {
    int32 i, j, k;
    int16 lsrec, oranges[20][2], nranges[20][2];
    int32 oi, ni, ontilts, erange, subsc_recs;
    int32 slat = 0, elat = 1, sscan = 0, escan = 1;
    int32 slon = 0, elon = 1;
    int32 edge[3] = {0, 0, 0}, start[3] = {0, 0, 0};

    subsc_recs = ((erec + 1) - (srec + 1)) / ysub + 1;

    new_tiltrec->ntilts = 0;
    for (i = 0; i < 20; i++) {
        new_tiltrec->tilt_flags[i] = 0;
        new_tiltrec->tilt_ranges[i][0] = 0;
        new_tiltrec->tilt_ranges[i][1] = 0;
        for (j = 0; j < 2; j++) {
            for (k = 0; k < 2; k++) {
                new_tiltrec->tilt_lats[i][j][k] = 0;
                new_tiltrec->tilt_lons[i][j][k] = 0;
            }
        }
    }

    for (i = 0; i < 20; i++)
        for (j = 0; j < 2; j++)
            oranges[i][j] = nranges[i][j] = 0;

    lsrec = srec + 1;

    ontilts = old_tiltrec->ntilts;

    for (i = 0; i < ontilts; i++)
        for (j = 0; j < 2; j++)
            oranges[i][j] = old_tiltrec->tilt_ranges[i][j];

    for (ni = 0, oi = 0; ((oi < ontilts) && (lsrec > oranges[oi][1])); oi++)
        ;

    nranges[ni][0] = 1;
    if (erec <= oranges[oi][1]) {
        nranges[ni][1] = subsc_recs;
        new_tiltrec->tilt_flags[ni] = old_tiltrec->tilt_flags[oi];
        new_tiltrec->ntilts++;
    } else {
        if (oranges[oi][1] == oranges[oi][0])
            nranges[ni][1] = oranges[oi][1];
        else
            nranges[ni][1] = ((oranges[oi][1] - lsrec) / ysub + 1);
        nranges[ni + 1][0] = nranges[ni][1] + 1;
        new_tiltrec->tilt_flags[ni] = old_tiltrec->tilt_flags[oi];
        new_tiltrec->ntilts++;
        for (++ni, ++oi; oi < ontilts; oi++) {
            nranges[ni][0] = nranges[ni - 1][1] + 1;
            if (erec >= oranges[oi][0]) {
                if (erec <= oranges[oi][1]) {
                    nranges[ni][1] = subsc_recs;
                    new_tiltrec->tilt_flags[ni] = old_tiltrec->tilt_flags[oi];
                    new_tiltrec->ntilts++;
                    ni++;
                } else {
                    erange = (oranges[oi][1] - oranges[oi][0]);
                    if (erange >= ysub) {
                        erange = erange / ysub;
                        nranges[ni][1] = nranges[oi][0] + erange;
                        new_tiltrec->tilt_flags[ni] = old_tiltrec->tilt_flags[oi];
                        new_tiltrec->ntilts++;
                        ni++;
                    }
                }
            }
        }
    }
    for (i = 0; i < new_tiltrec->ntilts; i++)
        for (j = 0; j < 2; j++)
            new_tiltrec->tilt_ranges[i][j] = nranges[i][j];

    for (i = 0; i < new_tiltrec->ntilts; i++) {
        new_tiltrec->tilt_lats[i][sscan][slat] = slatrec[nranges[i][0] - 1];
        new_tiltrec->tilt_lats[i][sscan][elat] = elatrec[nranges[i][0] - 1];

        new_tiltrec->tilt_lats[i][escan][slat] = slatrec[nranges[i][1] - 1];
        new_tiltrec->tilt_lats[i][escan][elat] = elatrec[nranges[i][1] - 1];

        new_tiltrec->tilt_lons[i][sscan][slon] = slonrec[nranges[i][0] - 1];
        new_tiltrec->tilt_lons[i][sscan][elon] = elonrec[nranges[i][0] - 1];

        new_tiltrec->tilt_lons[i][escan][slon] = slonrec[nranges[i][1] - 1];
        new_tiltrec->tilt_lons[i][escan][elon] = elonrec[nranges[i][1] - 1];
    }

    edge[0] = edge[1] = edge[2] = 0;
    if ((update_sds(sdfid, NTILTS, start, edge,
            (VOIDP) & new_tiltrec->ntilts)) < 0) return FAIL;

    edge[0] = edge[1] = edge[2] = 0;
    if ((update_sds(sdfid, TILTFLAGS, start, edge,
            (VOIDP) new_tiltrec->tilt_flags)) < 0) return FAIL;

    edge[0] = edge[1] = edge[2] = 0;
    if ((update_sds(sdfid, TILTRANGES, start, edge,
            (VOIDP) new_tiltrec->tilt_ranges)) < 0) return FAIL;

    edge[0] = edge[1] = edge[2] = 0;
    if ((update_sds(sdfid, TILTLATS, start, edge,
            (VOIDP) new_tiltrec->tilt_lats)) < 0) return FAIL;

    edge[0] = edge[1] = edge[2] = 0;
    if ((update_sds(sdfid, TILTLONS, start, edge,
            (VOIDP) new_tiltrec->tilt_lons)) < 0) return FAIL;

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  get_geodata

    Returns:   int32 (sdsid)
        The return code is a negative value if any error occurs, otherwise,
        returns sdsid.

    Description:
        The function set_flags reads l2_flags dataset for calculating the
        flag percentages

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     sdfid         I   SDinterface ID of output file
        int32     nrec          I   number of records
        int32     nsamp         I   number of samples
        int32     *flags        O   count of flags set for each band

    NOTE:

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     04/03/95  Original development
----------------------------------------------------------------------------*/
int32 set_flags(int32 sdfid, int32 nrec, int32 nsamp, int32 *flags) {
    int32 i, j, start[3] = {0, 0, 0}, edge[3] = {0, 0, 0};
    int16 *l2flags;
    int16 masks[16];

    if ((l2flags = (int16 *) calloc(nrec * nsamp, sizeof (int16))) == NULL)
        return MALLOC_ERR;

    if ((rdslice(sdfid, L2FLAGS, start, edge, (VOIDP) l2flags)) < 0)
        return FAIL;

    for (flags[0] = 0, masks[0] = 1, i = 1; i < 16; i++) {
        masks[i] = (masks[i - 1] * 2);
        flags[i] = 0;
    }

    for (i = 0; i < nrec * nsamp; i++)
        for (j = 0; j < 16; j++) {
            if ((l2flags[i] & masks[j]) == masks[j])
                flags[j]++;
        }

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  alloc_nav_buffs

    Returns:   int32 (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        The function alloc_nav_buffs allocates buffers for navigation data

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     nrec          I   number of records/scanlines
        GeoType  *nav_rec       O   structrue to hold allocated buffers

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     04/03/95  Original development
----------------------------------------------------------------------------*/
int32 alloc_nav_buffs(int32 nrec, NavType *nav_rec) {
    if ((nav_rec->orb_vec = (float32 *) calloc(nrec * 3, sizeof (float32))) == NULL)
        return MALLOC_ERR;

    if ((nav_rec->l_vert = (float32 *) calloc(nrec * 3, sizeof (float32))) == NULL)
        return MALLOC_ERR;

    if ((nav_rec->sun_ref = (float32 *) calloc(nrec * 3, sizeof (float32))) == NULL)
        return MALLOC_ERR;

    if ((nav_rec->att_ang = (float32 *) calloc(nrec * 3, sizeof (float32))) == NULL)
        return MALLOC_ERR;

    if ((nav_rec->sen_mat = (float32 *) calloc(nrec * 3 * 3, sizeof (float32))) == NULL)
        return MALLOC_ERR;

    if ((nav_rec->scan_ell = (float32 *) calloc(nrec * 6, sizeof (float32))) == NULL)
        return MALLOC_ERR;

    if ((nav_rec->nflag = (int32 *) calloc(nrec * 8, sizeof (int32))) == NULL)
        return MALLOC_ERR;

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  get_navdata

    Returns:   int32 (sdsid)
        The return code is a negative value if any error occurs, otherwise,
        returns sdsid.

    Description:
        The function get_navdata reads sensor metrics info from the given
        input data file

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     sdfid         I   SDinterface ID of input file
        int32     nrec          I   number of input data records
        NavType *nav_rec    O   sensor metrics data structure

    NOTE:

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     04/03/95  Original development
----------------------------------------------------------------------------*/
int32 get_navdata(int32 sdfid, int32 nrec, NavType *nav_rec) {
    int32 start[3] = {0, 0, 0}, edge[3] = {0, 0, 0};

    edge[0] = nrec;
    edge[1] = 3;
    edge[2] = 0;
    if ((rdslice(sdfid, ORBVEC, start, edge, (VOIDP) nav_rec->orb_vec)) < 0)
        return FAIL;

    edge[0] = nrec;
    edge[1] = 3;
    edge[2] = 0;
    if ((rdslice(sdfid, LVERT, start, edge, (VOIDP) nav_rec->l_vert)) < 0)
        return FAIL;

    edge[0] = nrec;
    edge[1] = 3;
    edge[2] = 0;
    if ((rdslice(sdfid, SUNREF, start, edge, (VOIDP) nav_rec->sun_ref)) < 0)
        return FAIL;

    edge[0] = nrec;
    edge[1] = 3;
    edge[2] = 0;
    if ((rdslice(sdfid, ATTANG, start, edge, (VOIDP) nav_rec->att_ang)) < 0)
        return FAIL;

    edge[0] = nrec;
    edge[1] = 3;
    edge[2] = 3;
    if ((rdslice(sdfid, SENMAT, start, edge, (VOIDP) nav_rec->sen_mat)) < 0)
        return FAIL;

    edge[0] = nrec;
    edge[1] = 6;
    edge[2] = 0;
    if ((rdslice(sdfid, SCANELL, start, edge, (VOIDP) nav_rec->scan_ell)) < 0)
        return FAIL;

    edge[0] = nrec;
    edge[1] = 8;
    edge[2] = 0;
    if ((rdslice(sdfid, NFLAG, start, edge, (VOIDP) nav_rec->nflag)) < 0)
        return FAIL;

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  free_nav_buffs

    Returns:   None

    Description:
        The function free_nav_buffs frees all the sensor metrics buffers

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        NavType *nav_rec    I   structrue containing sensor metrics
                                        buffers
    NOTE:

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     04/03/95  Original development
----------------------------------------------------------------------------*/
void
#ifdef PRTOTYPE
free_nav_buffs(NavType *nav_rec)
#else
free_nav_buffs(nav_rec)
NavType *nav_rec;
#endif
{
    free(nav_rec->orb_vec);
    free(nav_rec->l_vert);
    free(nav_rec->sun_ref);
    free(nav_rec->att_ang);
    free(nav_rec->sen_mat);
    free(nav_rec->scan_ell);
    free(nav_rec->nflag);
}

/*-----------------------------------------------------------------------------
    Function:  subsamp_rec

    Returns:   int32 (sdsid)
        The return code is a negative value if any error occurs, otherwise,
        returns sdsid.

    Description:
        The function subsamp_rec subsamples the given data along the last
        dimension ('C' order)

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     srec          I   starting record
        int32     maxrec        I   last record that should be processed
        int32     rank          I   number of data dimensions
        int32     *dims         I   dimensions of the data
        int32     ysub          I   scanline subsampling rate
        int32     nt            I   number type of data
        void      *obuf         I   input data buffer
        void      *nbuf         O   output data buffer

    NOTE:

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     04/03/95  Original development
----------------------------------------------------------------------------*/
int32 subsamp_rec(int32 srec, int32 maxrec, int32 rank, int32 *dims,
        int32 ysub, int32 nt, void *obuf, void *nbuf) {
    int32 i, j, nv, ni, oi;
    int32 *i32_obuf, *i32_nbuf;
    int16 *i16_obuf, *i16_nbuf;
    uint8 *ui8_obuf, *ui8_nbuf;
    float32 *f32_obuf, *f32_nbuf;


    for (nv = 1, i = 1; i < rank; i++)
        nv *= dims[i];

    for (ni = 0, oi = nv * srec; oi <= maxrec * nv; oi += ysub * nv) {
        if (nt == DFNT_INT32) {
            i32_obuf = (int32 *) obuf;
            i32_nbuf = (int32 *) nbuf;
            for (j = 0; j < nv; j++)
                i32_nbuf[ni++] = i32_obuf[oi + j];
        } else
            if (nt == DFNT_INT16) {
            i16_obuf = (int16 *) obuf;
            i16_nbuf = (int16 *) nbuf;
            for (j = 0; j < nv; j++)
                i16_nbuf[ni++] = i16_obuf[oi + j];
        } else
            if (nt == DFNT_UINT8) {
            ui8_obuf = (uint8 *) obuf;
            ui8_nbuf = (uint8 *) nbuf;
            for (j = 0; j < nv; j++)
                ui8_nbuf[ni++] = ui8_obuf[oi + j];
        } else
            if (nt == DFNT_FLOAT32) {
            f32_obuf = (float32 *) obuf;
            f32_nbuf = (float32 *) nbuf;
            for (j = 0; j < nv; j++)
                f32_nbuf[ni++] = f32_obuf[oi + j];
        }
    }
    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  subsamp_2D

    Returns:   int32 (sdsid)
        The return code is a negative value if any error occurs, otherwise,
        returns sdsid.

    Description:
        The function subsamp_2D subsamples the given data along x and y

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     spix         I    start pixel position
        int32     sscan        I    start scan line (y pos) position
        int32     epix         I    ending pixel position
        int32     xdim         I    no. of pixels (x dim of inbuf)
        int32     xsub         I    pixel subsampling rate
        int32     ysub         I    scan line subsampling rate
        int32     nt           I    data type (HDF number type) of data
        void     *inbuf        I    input data buffer
        void     *outbuf       O    output data buffer (calling routine should
                                        provide proper sized buffer

    NOTE:

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     04/03/95  Original development
----------------------------------------------------------------------------*/
void
subsamp_2D(int32 spix, int32 sscan, int32 epix, int32 escan, int32 xdim,
        int32 xsub, int32 ysub, int32 nt, void *inbuf, void *outbuf) {
    int32 newrec, rec, px, newpx, newrow, npix = xdim;
    int32 max_rec_used, subsc_samps, subsc_recs;
    int32 *i32_obuf, *i32_ibuf;
    int16 *i16_obuf, *i16_ibuf;
    uint8 *ui8_obuf, *ui8_ibuf;
    float32 *f32_obuf, *f32_ibuf;

    /*** compute subscene image dimensions */

    subsc_recs = ((escan + 1) - (sscan + 1)) / ysub + 1;
    subsc_samps = ((epix + 1) - (spix + 1)) / xsub + 1;

    /*** compute the maximum record and sample numbers used to make the subsampled
            browse output */

    max_rec_used = (subsc_recs - 1) * ysub + sscan;

    /*** subsample the data */

    if (nt == DFNT_INT32) {
        i32_obuf = (int32 *) outbuf;
        i32_ibuf = (int32 *) inbuf;
        for (newrec = 0, rec = sscan; rec <= max_rec_used; rec += ysub, newrec++) {
            newrow = newrec * subsc_samps;
            for (px = rec * npix + (spix), newpx = newrow;
                    newpx < newrow + subsc_samps; px += xsub, newpx++) {
                i32_obuf[newpx] = i32_ibuf[px];
            }
        }
    } else
        if (nt == DFNT_INT16) {
        i16_obuf = (int16 *) outbuf;
        i16_ibuf = (int16 *) inbuf;
        for (newrec = 0, rec = sscan; rec <= max_rec_used; rec += ysub, newrec++) {
            newrow = newrec * subsc_samps;
            for (px = rec * npix + (spix), newpx = newrow;
                    newpx < newrow + subsc_samps; px += xsub, newpx++) {
                i16_obuf[newpx] = i16_ibuf[px];
            }
        }
    } else
        if (nt == DFNT_UINT8) {
        ui8_obuf = (uint8 *) outbuf;
        ui8_ibuf = (uint8 *) inbuf;
        for (newrec = 0, rec = sscan; rec <= max_rec_used;
                rec += ysub, newrec++) {
            newrow = newrec * subsc_samps;
            for (px = rec * npix + (spix), newpx = newrow;
                    newpx < newrow + subsc_samps; px += xsub, newpx++) {
                ui8_obuf[newpx] = ui8_ibuf[px];
            }
        }
    } else
        if (nt == DFNT_FLOAT32) {
        f32_obuf = (float32 *) outbuf;
        f32_ibuf = (float32 *) inbuf;
        for (newrec = 0, rec = sscan; rec <= max_rec_used;
                rec += ysub, newrec++) {
            newrow = newrec * subsc_samps;
            for (px = rec * npix + (spix), newpx = newrow;
                    newpx < newrow + subsc_samps; px += xsub, newpx++) {
                f32_obuf[newpx] = f32_ibuf[px];
            }
        }
    }
}

/*-----------------------------------------------------------------------------
    Function:  set_l1adata

    Returns:   int32 (sdsid)
        The return code is a negative value if any error occurs, otherwise,
        returns sdsid.

    Description:
        The function set_l1adata reads level 1a data, gain, and dark_rest
        from the given input file, calls subsamp_l1adata to subsample the
        data and writes it to the output file

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     isdfid        I    ID req to access HDF SDS interface
        int32     osdfid        I    SDinterface ID of output file
        int32     num_recs      I    number of input records
        int32     num_samps     I    number of input samples
        int32     subsc_recs    I    number of output records
        int32     subsc_samps   I    number of output samples
        int32     srec          I    starting record
        int32     ssamp         I    starting pixel
        int32     erec          I    ending record
        int32     esamp         I    ending sample
        int32     xsub          I    pixel subsampling rate
        int32     ysub          I    scanline subsampling rate
        FilemetricsType *fm_rec I    File metrics

    NOTE:

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     04/03/95  Original development
    Lakshmi Kumar  Hughes STX     12/12/95  Fixed a bug wh was trying to write
                                                more than subsc_recs
    Lakshmi Kumar  Hughes STX     02/21/96  Fixed the above bugfix which had
                                                introduced another bug
----------------------------------------------------------------------------*/
int32
set_l1adata(int32 isdfid, int32 osdfid, int32 num_recs, int32 num_samps,
        int32 subsc_recs, int32 subsc_samps, int32 srec, int32 ssamp,
        int32 erec, int32 esamp, int32 xsub, int32 ysub,
        FilemetricsType *fm_rec) {
    int32 nsamp, nrec, last_rec_used, new_buf_recs;
    int32 old_start[3] = {0, 0, 0}, new_start[3] = {0, 0, 0};
    int32 old_edge[3], new_edge[3], rd_recs, buf_srec, buf_erec;
    int16 nbands = NBANDS, *gain, *dark_rest, *inbuf, *outbuf;
    int16 *s_satp, *s_zerop;
    int32 start[3] = {0, 0, 0}, edge[3];

    nrec = num_recs;
    nsamp = num_samps;

    old_edge[0] = nrec;
    old_edge[1] = nsamp;
    old_edge[2] = nbands;

    new_edge[0] = subsc_recs;
    new_edge[1] = subsc_samps;
    new_edge[2] = nbands;

    if ((gain = (int16 *) calloc(subsc_recs * NBANDS, sizeof (int16))) == NULL)
        return MALLOC_ERR;

    if ((dark_rest = (int16 *) calloc(subsc_recs * NBANDS, sizeof (int16))) == NULL)
        return MALLOC_ERR;

    if ((s_satp = (int16 *) calloc(subsc_recs * NBANDS, sizeof (int16))) == NULL)
        return MALLOC_ERR;

    if ((s_zerop = (int16 *) calloc(subsc_recs * NBANDS, sizeof (int16))) == NULL)
        return MALLOC_ERR;

    /*** read gain and dark_rest datasets */
    edge[0] = edge[1] = edge[2] = 0;
    if ((rdslice(osdfid, GAIN, start, edge, (VOIDP) gain)) < 0)
        return FAIL;

    edge[0] = edge[1] = edge[2] = 0;
    if ((rdslice(osdfid, DARKREST, start, edge, (VOIDP) dark_rest)) < 0)
        return FAIL;

    nrec = erec - srec + 1;
    old_start[0] = srec;
    new_start[0] = -1;
    if (nrec > CHUNK_SZ)
        rd_recs = CHUNK_SZ;
    else
        rd_recs = nrec;

    /****  allocate space for l1a buffers */

    if ((inbuf =
            (int16 *) calloc(CHUNK_SZ * old_edge[1] * old_edge[2], sizeof (int16))) == NULL) {
        sprintf(ERR_MSG, "Calloc error - cannot allocate memory for reading %s",
                L1ADATA);
        return MALLOC_ERR;
    }

    if ((outbuf =
            (int16 *) calloc(CHUNK_SZ * new_edge[1] * new_edge[2],
            sizeof (int16))) == NULL) {
        sprintf(ERR_MSG, "calloc error - cannot allocate memory for writing %s",
                L1ADATA);
        return MALLOC_ERR;
    }

    while (rd_recs > 0) {
        old_edge[0] = rd_recs;
        if ((rdslice(isdfid, L1ADATA, old_start, old_edge, (VOIDP) inbuf)) < 0)
            return FAIL;
        buf_srec = 0;
        buf_erec = rd_recs - 1;
        last_rec_used = subsamp_l1adata(buf_srec, buf_erec, ssamp, esamp, xsub,
                ysub, nsamp, gain, dark_rest, fm_rec, inbuf, outbuf,
                s_satp, s_zerop, &new_buf_recs);
        old_start[0] += ysub + last_rec_used;
        nrec -= last_rec_used + ysub;
        if (nrec < rd_recs)
            rd_recs = nrec;

        if (new_start[0] == -1) {
            new_start[0] = 0;
            new_edge[0] = new_buf_recs;
        } else {
            new_start[0] += new_edge[0];
            new_edge[0] = new_buf_recs;
        }
        if ((update_sds(osdfid, L1ADATA, new_start, new_edge, (VOIDP) outbuf)) < 0)
            return FAIL;
    } /* end while */

    new_start[0] = new_start[1] = new_start[2] = 0;
    new_edge[0] = subsc_recs;
    new_edge[1] = NBANDS;
    new_edge[2] = 0;
    if ((update_sds(osdfid, SSATP, new_start, new_edge, (VOIDP) s_satp)) < 0)
        return FAIL;
    if ((update_sds(osdfid, SZEROP, new_start, new_edge, (VOIDP) s_zerop)) < 0)
        return FAIL;

    free(s_satp);
    free(s_zerop);
    free(gain);
    free(dark_rest);
    free(inbuf);
    free(outbuf);
    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  subsamp_l1adata

    Returns:   int32 (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        The function subsamp_l1adata subsamples the given level 1a 3D data
        and calculates gain1, gain2 saturated and non-saturated pixels, zero
        pixels and mean gain1 and mean gain2 radiances

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     srec          I   starting scan line
        int32     erec          I   ending scan line
        int32     ssamp         I   starting sample
        int32     esamp         I   ending sample
        int32     xsub          I   pixel subsampling rate
        int32     ysub          I   scanline subsampling rate
        int32     nsamp         I   number of samples
        int16    *gain          I   gain dataset
        int16    *dark_rest     I   dark_restore data set
        FilemetricsType *fm_rec I   File metrics
        int16    *inbuf         I   input level1a data
        int16    *outbuf        O   output level1a data
        int32    *new_buf_recs  O   output buffer records

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     04/03/95  Original development
    Lakshmi Kumar  Hughes STX     12/14/95  Fixed the return value (before
                                                it was returning wrong value)
----------------------------------------------------------------------------*/
int32
subsamp_l1adata(int32 srec, int32 erec, int32 ssamp, int32 esamp,
        int32 xsub, int32 ysub, int32 nsamp, int16 *gain,
        int16 *dark_rest, FilemetricsType *fm_rec,
        int16 *inbuf, int16 *outbuf, int16 *s_satp,
        int16 *s_zerop, int32 *new_buf_recs) {
    int32 rec_sz, band, px_index, out_index, rec, rec_index;
    int32 subsc_recs, subsc_samps;
    int16 radiance, i = 0;
    int32 last_rec_used, last_samp_used;
    static int32 gn1px[NBANDS], gn2px[NBANDS];
    static double gainr1[NBANDS], gainr2[NBANDS];
    static int32 first_call = 1, gn_index;

    if (first_call) {
        first_call = 0;
        gn_index = 0;
        for (i = 0; i < NBANDS; i++) {
            fm_rec->satg1[i] = fm_rec->satg2[i] = 0;
            fm_rec->nsatg1[i] = fm_rec->nsatg2[i] = 0;
            fm_rec->zeroes[i] = 0;
            fm_rec->meanr1[i] = fm_rec->meanr2[i] = 0;
            gn1px[i] = gn2px[i] = 0;
            gainr1[i] = gainr2[i] = 0;
        }
    }

    *new_buf_recs = subsc_recs = ((erec + 1) - (srec + 1)) / ysub + 1;
    subsc_samps = ((esamp + 1) - (ssamp + 1)) / xsub + 1;

    last_rec_used = (subsc_recs - 1) * ysub + srec;
    last_samp_used = (subsc_samps - 1) * xsub + ssamp;

    rec_sz = nsamp * NBANDS;
    for (out_index = 0, rec = srec; rec <= last_rec_used;
            rec += ysub, gn_index++) {
        rec_index = rec * rec_sz;
        for (px_index = rec_index + (ssamp * NBANDS);
                px_index <= (rec_index + last_samp_used * NBANDS);
                px_index += (xsub * NBANDS)) {
            for (band = 0; band < NBANDS; band++) {
                radiance = outbuf[out_index++] = inbuf[px_index + band];
                if (gain[gn_index + band] == 0) {
                    gainr1[band] += radiance;
                    gn1px[band]++;
                    if (radiance >= 1023) {
                        s_satp[gn_index * NBANDS + band]++;
                        fm_rec->satg1[band]++;
                    } else
                        fm_rec->nsatg1[band]++;
                } else
                    if (gain[gn_index * NBANDS + band] == 2) {
                    gainr2[band] += radiance;
                    gn2px[band]++;
                    if (radiance >= 1023) {
                        s_satp[gn_index + band]++;
                        fm_rec->satg2[band]++;
                    } else
                        fm_rec->nsatg2[band]++;
                }
                if ((radiance - dark_rest[gn_index * NBANDS + band]) <= 2) {
                    s_zerop[gn_index * NBANDS + band]++;
                    fm_rec->zeroes[band]++;
                }
            }
        }
    }
    for (band = 0; band < NBANDS; band++) {
        if (gn1px[band] > 0)
            fm_rec->meanr1[band] = gainr1[band] / gn1px[band];
        if (gn2px[band] > 0)
            fm_rec->meanr2[band] = gainr2[band] / gn2px[band];
    }

    return last_rec_used;
}

/*-----------------------------------------------------------------------------
    Function:  set_globalattrs

    Returns:   int32 (sdsid)
        The return code is a negative value if any error occurs, otherwise,
        returns sdsid.

    Description:
        The function set_globalattrs updates data dependent global attributes

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        char     *outfile       I   output file name
        int32     sdfid         I   SDinterface ID of output file
        int32     l1aflag       I   flag indicating l1a or l2 file
        int32     nrec          I   number of records
        int32     nsamp         I   number of samples
        int32     ssamp         I   starting sample
        int32     xsub          I   pixle subsampling rate
        char      *dtype        I   data type (GAC, LAC, etc.,)
        int16     *flags        I   level2 flags
        GeoType   *geo_rec      I   structure containing scene coordinate data
        FilemetricsType *fm_rec I   File metrics data structure

    NOTE:

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     11/02/94  Original development
----------------------------------------------------------------------------*/
int32
set_globalattrs(char *outfile, int32 sdfid, int32 l1aflag, int32 nrec,
        int32 nsamp, int32 ssamp, int32 xsub, int32 *flags,
        GeoType *geo_rec, FilemetricsType *fm_rec) {
    uint8 *s_flags;
    int16 dday, dyear, syear, sday, eyear, eday;
    int16 subsc_syear, subsc_sday, subsc_eyear, subsc_eday;
    int32 emsec, smsec, pix_sub, pix_start;
    int32 i, rec, cmsec, FSlines = 0;
    int32 *msec, cscan, *nflag, xfirst, xlast, xcenter;
    int32 start[3] = {0, 0, 0}, edge[3] = {0, 0, 0};
    div_t quot1, quot2, quot3;
    char string[MAXVAL], ptime[17], *prod_name;
    float32 maxlat = -999, minlat = 999, maxlon = -999, minlon = 999;
    float32 pct_flags[16];

    /*** read msec dataset and data time */
    if ((msec = (int32 *) calloc(nrec, sizeof (int32))) == NULL)
        return MALLOC_ERR;

    if ((rdslice(sdfid, MSEC, start, edge, (VOIDP) msec)) < 0)
        return FAIL;

    if ((s_flags = (uint8 *) calloc(nrec * 4, sizeof (uint8))) == NULL)
        return MALLOC_ERR;

    edge[0] = edge[1] = edge[2] = 0;
    if ((rdslice(sdfid, SFLAGS, start, edge, (VOIDP) s_flags)) < 0)
        return FAIL;
    /*
       get nflag information for deciding the start, center, and end scan line that
       have the good navigation  G.Fu 3/2/98
     */
    if ((nflag = (int32 *) calloc(nrec * 8, sizeof (int32))) == NULL)
        return MALLOC_ERR;

    edge[0] = edge[1] = edge[2] = 0;
    if ((rdslice(sdfid, NFLAG, start, edge, (VOIDP) nflag)) < 0)
        return FAIL;

    if ((rdattr(sdfid, SYEAR, &syear)) < 0)
        return FAIL;
    if ((rdattr(sdfid, SDAY, &sday)) < 0)
        return FAIL;
    if ((rdattr(sdfid, SMSEC, &smsec)) < 0)
        return FAIL;
    if ((rdattr(sdfid, EYEAR, &eyear)) < 0)
        return FAIL;
    if ((rdattr(sdfid, EDAY, &eday)) < 0)
        return FAIL;
    if ((rdattr(sdfid, EMSEC, &emsec)) < 0)
        return FAIL;
    if ((rdattr(sdfid, PIX_START, &pix_start)) < 0)
        return FAIL;
    if ((rdattr(sdfid, PIX_SUB, &pix_sub)) < 0)
        return FAIL;

    /* Removed, BAF, 12/2002
       smsec = msec[0];
       emsec = msec[nrec-1];
     */
    cscan = (nrec + 1) / 2;
    cmsec = (msec[cscan - 1]);


    /*** check if scene crosses day  */
    if (msec[0] < smsec) {
        subsc_syear = eyear;
        subsc_sday = eday;
        subsc_eyear = eyear;
        subsc_eday = eday;
    } else {
        subsc_syear = syear;
        subsc_sday = sday;
        if (msec[nrec - 1] < emsec) {
            subsc_eyear = eyear;
            subsc_eday = eday;
        } else {
            subsc_eyear = syear;
            subsc_eday = sday;
        }
    }

    /* add 2 lines, BAF, 12/2002 */
    smsec = msec[0];
    emsec = msec[nrec - 1];

    quot1 = div(smsec, MSECHOUR);
    quot2 = div(quot1.rem, MSECMIN);
    quot3 = div(quot2.rem, MSECSEC);
    sprintf(string, "%4.4d%3.3d%2.2d%2.2d%2.2d%3.3d", subsc_syear,
            subsc_sday, quot1.quot, quot2.quot, quot3.quot, quot3.rem);
    SDsetattr(sdfid, "Start Time", DFNT_CHAR, strlen(string) + 1,
            (VOIDP) string);
    quot1 = div(emsec, MSECHOUR);
    quot2 = div(quot1.rem, MSECMIN);
    quot3 = div(quot2.rem, MSECSEC);
    sprintf(string, "%4.4d%3.3d%2.2d%2.2d%2.2d%3.3d", subsc_eyear,
            subsc_eday, quot1.quot, quot2.quot, quot3.quot, quot3.rem);
    SDsetattr(sdfid, "End Time", DFNT_CHAR, strlen(string) + 1,
            (VOIDP) string);

    SDsetattr(sdfid, SYEAR, DFNT_INT16, 1, (VOIDP) & subsc_syear);
    SDsetattr(sdfid, SDAY, DFNT_INT16, 1, (VOIDP) & subsc_sday);
    SDsetattr(sdfid, SMSEC, DFNT_INT32, 1, (VOIDP) & msec[0]);
    SDsetattr(sdfid, EYEAR, DFNT_INT16, 1, (VOIDP) & subsc_eyear);
    SDsetattr(sdfid, EDAY, DFNT_INT16, 1, (VOIDP) & subsc_eday);
    SDsetattr(sdfid, EMSEC, DFNT_INT32, 1, (VOIDP) & msec[nrec - 1]);

    dyear = subsc_syear;
    dday = subsc_sday;
    if (subsc_sday != subsc_eday || subsc_syear != subsc_eyear) {
        if (cmsec < 43200000) {
            dday = subsc_eday;
            dyear = subsc_eyear;
        }
    }
    quot1 = div(cmsec, MSECHOUR);
    quot2 = div(quot1.rem, MSECMIN);
    quot3 = div(quot2.rem, MSECSEC);
    sprintf(string, "%4.4d%3.3d%2.2d%2.2d%2.2d%3.3d", dyear, dday,
            quot1.quot, quot2.quot, quot3.quot, quot3.rem);
    SDsetattr(sdfid, CTIME, DFNT_CHAR, strlen(string) + 1,
            (VOIDP) string);

    /*** write Data Quality global attributes */

    SDsetattr(sdfid, NSAMP, DFNT_INT32, 1, (VOIDP) & nsamp);

    SDsetattr(sdfid, NREC, DFNT_INT32, 1, (VOIDP) & nrec);

    pix_start = pix_start + (ssamp * pix_sub);
    SDsetattr(sdfid, PIX_START, DFNT_INT32, 1, (VOIDP) & pix_start);

    pix_sub = pix_sub * xsub;
    SDsetattr(sdfid, PIX_SUB, DFNT_INT32, 1, (VOIDP) & pix_sub);

    if (!l1aflag) {
        for (i = 0; i < 16; i++) {
            if (flags[i] > 0)
                pct_flags[i] = flags[i] / (nsamp * nrec * 1.0) * 100;
            else
                pct_flags[i] = 0;
        }
        SDsetattr(sdfid, PCT_FLAG, DFNT_FLOAT32, 16, (VOIDP) pct_flags);
    } else { /* update file metrics */
        SDsetattr(sdfid, L1SATG1, DFNT_INT32, 8, (VOIDP) fm_rec->satg1);

        SDsetattr(sdfid, L1SATG2, DFNT_INT32, 8, (VOIDP) fm_rec->satg2);

        SDsetattr(sdfid, L1NSATG1, DFNT_INT32, 8, (VOIDP) fm_rec->nsatg1);

        SDsetattr(sdfid, L1NSATG2, DFNT_INT32, 8, (VOIDP) fm_rec->nsatg2);

        SDsetattr(sdfid, L1ZEROS, DFNT_INT32, 8, (VOIDP) fm_rec->zeroes);

        SDsetattr(sdfid, L1MEANR1, DFNT_FLOAT32, 8, (VOIDP) fm_rec->meanr1);

        SDsetattr(sdfid, L1MEANR2, DFNT_FLOAT32, 8, (VOIDP) fm_rec->meanr2);
    }

    /*** write scene coordinates global attributes */
    /*
        the following global attributes will be decided based on scan lines with
        good navigation only (to be consistent with l1agen program  3/2/98 G. Fu
     */
    xfirst = -1;
    xlast = nrec - 1;
    xcenter = cscan - 1;
    for (i = 0; i < nrec; i++) {
        if (nflag[i * 8] == 0) {
            if (xfirst < 0)
                xfirst = i;
            xlast = i;
        }
    }

    if (xfirst >= 0) {
        xcenter = (xlast + xfirst) / 2;
        while ((xcenter > xfirst) && (nflag[xcenter * 8]) != 0)
            xcenter--;

        xcenter++; /* change to 1-rel */
    }
    else {
        xfirst = 0;
        xlast = nrec - 1;
        xcenter = (nrec + 1) / 2;
        printf("\nWarning - No good navigation record found for the data extracted\n");
    }
    SDsetattr(sdfid, NCREC, DFNT_INT32, 1, (VOIDP) & xcenter);
    SDsetattr(sdfid, SCLAT, DFNT_FLOAT32, 1, (VOIDP) & geo_rec->clat[xcenter - 1]);
    SDsetattr(sdfid, SCLON, DFNT_FLOAT32, 1, (VOIDP) & geo_rec->clon[xcenter - 1]);
    SDsetattr(sdfid, ULLAT, DFNT_FLOAT32, 1, (VOIDP) & geo_rec->slat[xfirst]);
    SDsetattr(sdfid, ULLON, DFNT_FLOAT32, 1, (VOIDP) & geo_rec->slon[xfirst]);
    SDsetattr(sdfid, URLAT, DFNT_FLOAT32, 1, (VOIDP) & geo_rec->elat[xfirst]);
    SDsetattr(sdfid, URLON, DFNT_FLOAT32, 1, (VOIDP) & geo_rec->elon[xfirst]);
    SDsetattr(sdfid, LLLAT, DFNT_FLOAT32, 1, (VOIDP) & geo_rec->slat[xlast]);
    SDsetattr(sdfid, LLLON, DFNT_FLOAT32, 1, (VOIDP) & geo_rec->slon[xlast]);
    SDsetattr(sdfid, LRLAT, DFNT_FLOAT32, 1, (VOIDP) & geo_rec->elat[xlast]);
    SDsetattr(sdfid, LRLON, DFNT_FLOAT32, 1, (VOIDP) & geo_rec->elon[xlast]);

    for (rec = 0; rec < nrec; rec++) {
        if (geo_rec->slat[rec] > maxlat)
            maxlat = geo_rec->slat[rec];
        if (geo_rec->elat[rec] < minlat)
            minlat = geo_rec->elat[rec];

        /*** also calculate filled scan lines */
        if (s_flags[rec * 4 + 1] > 0)
            FSlines++;
    }

    /*  check if the scene crosses dateline  */
    if ((geo_rec->slon[0] < 0) && (geo_rec->slon[nrec - 1] >= 0)) {
        for (rec = 0; rec < nrec; rec++) { /* find the least +ve no. */
            if ((geo_rec->slon[rec] < minlon) && (geo_rec->slon[rec] > 0))
                minlon = geo_rec->slon[rec];
            if ((geo_rec->elon[rec] > maxlon) && (geo_rec->elon[rec] < 0))
                maxlon = geo_rec->elon[rec];
        }
    } else {
        for (rec = 0; rec < nrec; rec++) {
            if (geo_rec->slon[rec] < minlon)
                minlon = geo_rec->slon[rec];
            if (geo_rec->elon[rec] > maxlon)
                maxlon = geo_rec->elon[rec];
        }
    }

    SDsetattr(sdfid, NFREC, DFNT_INT32, 1, (VOIDP) & FSlines);

    SDsetattr(sdfid, NORTHLAT, DFNT_FLOAT32, 1, (VOIDP) & maxlat);

    SDsetattr(sdfid, SOUTHLAT, DFNT_FLOAT32, 1, (VOIDP) & minlat);

    SDsetattr(sdfid, WESTLON, DFNT_FLOAT32, 1, (VOIDP) & minlon);

    SDsetattr(sdfid, EASTLON, DFNT_FLOAT32, 1, (VOIDP) & maxlon);

    SDsetattr(sdfid, BCLAT, DFNT_FLOAT32, 1, (VOIDP) & geo_rec->clat[0]);
    SDsetattr(sdfid, BCLON, DFNT_FLOAT32, 1, (VOIDP) & geo_rec->clon[0]);

    SDsetattr(sdfid, ECLAT, DFNT_FLOAT32, 1, (VOIDP) & geo_rec->clat[nrec - 1]);
    SDsetattr(sdfid, ECLON, DFNT_FLOAT32, 1, (VOIDP) & geo_rec->clon[nrec - 1]);

    get_time(ptime);

    SDsetattr(sdfid, PTIME, DFNT_CHAR, strlen(ptime) + 1, (VOIDP) ptime);

    if ((prod_name = strrchr(outfile, '/')) != NULL)
        prod_name++;
    else
        prod_name = outfile;
    SDsetattr(sdfid, PROD_NAME, DFNT_CHAR, strlen(prod_name) + 1, (VOIDP) prod_name);

    free(msec);
    free(s_flags);

    return SUCCEED;
}
