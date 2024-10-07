/* =========================================================== */
/* Module hdf_utils.c                                          */
/*                                                             */
/* HDF I/O utilities.                                          */
/*                                                             */
/* Written By:                                                 */
/*     Norman Kuring, NASA/GSFC                                */
/*                                                             */
/* Modification History:                                       */
/*     B. A. Franz, SAIC GSC, January 1999,                    */
/*     Moved general utilities from SWl01 to this module.      */
/*     Eliminated global usage of sds_id.                      */
/* =========================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <time.h>
#include <math.h>
#include <hdf4utils.h>

#include <hdf.h>
#include <mfhdf.h>

/************************************************************************
return the size of the HDF data type in bytes.  Return 0 if not found.
 ************************************************************************/
int32_t hdf_sizeof(int32_t dtype) {
    switch (dtype) {
    case DFNT_CHAR8:
    case DFNT_UCHAR8:
    case DFNT_INT8:
    case DFNT_UINT8:
    case DFNT_NINT8:
    case DFNT_NUINT8:
    case DFNT_LINT8:
    case DFNT_LUINT8:
        return 1;
    case DFNT_INT16:
    case DFNT_UINT16:
    case DFNT_NINT16:
    case DFNT_NUINT16:
    case DFNT_LINT16:
    case DFNT_LUINT16:
        return 2;
    case DFNT_INT32:
    case DFNT_UINT32:
    case DFNT_FLOAT32:
    case DFNT_NINT32:
    case DFNT_NUINT32:
    case DFNT_NFLOAT32:
    case DFNT_LINT32:
    case DFNT_LUINT32:
    case DFNT_LFLOAT32:
        return 4;
    case DFNT_INT64:
    case DFNT_UINT64:
    case DFNT_FLOAT64:
    case DFNT_NINT64:
    case DFNT_NUINT64:
    case DFNT_NFLOAT64:
    case DFNT_LINT64:
    case DFNT_LUINT64:
    case DFNT_LFLOAT64:
        return 8;
    default:
        return 0;
    }
    return 0;
}


/* -------------------------------------------------------- */
/* Create an SDS using the wrappers for the HDF routines     */

/* -------------------------------------------------------- */
int CreateSDS(
        int32_t sd_id,
        const char *sname, /* short name */
        const char *lname, /* long name */
        const char *standard_name, /* NetCDF standard name (not set if passed NULL or "") */
        const char *units, /* units (not set if passed NULL or "") */
        double low, /* low end of valid range */
        double high, /* high end of range (no range set if low >= high) */
        float slope, /* scale factor (not set if 1)  */
        float offset, /* scaling offset (not set if 0)  */
        int32_t nt, /* HDF number type */
        int32_t rank, /* number of dimensions (must be <= 3) */
        int32_t d0, /* size of 1st dimension */
        int32_t d1, /* size of 2nd dimension */
        int32_t d2, /* size of 3rd dimension (1 if rank < 2) */
        const char *dn0, /* name of 1st dimension */
        const char *dn1, /* name of 2nd dimension (NULL if rank < 2) */
        const char *dn2 /* name of 3rd dimension (NULL if rank < 3) */
        ) {

    int32_t sds_id;

    if (rank < 3) {
        d2 = 1;
        dn2 = NULL;
    }
    if (rank < 2) {
        d1 = 1;
        dn1 = NULL;
    }

    /* Create the SDS */
    PTB(sd_create(sd_id, sname, nt, rank, d0, d1, d2, &sds_id));

    /* Name its dimensions */
    if (d0 != SD_UNLIMITED) {
        PTB(sd_setdimnames(sds_id, dn0, dn1, dn2));
    } else if (SDsetblocksize(sds_id, 4194304) == FAIL) {
        printf("-E- %s line %d: Could not enlarge the HDF block size.\n", __FILE__, __LINE__);
        return (HDF_FUNCTION_ERROR);
    }

    /* Add a "long_name" attribute */
    PTB(sd_setattr(sds_id, "long_name", DFNT_CHAR, strlen(lname) + 1, lname));

    /* Add a "valid_range" attribute if one is specified */
    if (low < high) {
        switch (nt) { /* Use the appropriate number type */
        case DFNT_INT8:
        {
            int8_t vr[2];
            vr[0] = (int8_t) low;
            vr[1] = (int8_t) high;
            PTB(sd_setattr(sds_id, "valid_range", DFNT_INT8, 2, vr));
        }
            break;
        case DFNT_UINT8:
        {
            uint8_t vr[2];
            vr[0] = (uint8_t) low;
            vr[1] = (uint8_t) high;
            PTB(sd_setattr(sds_id, "valid_range", DFNT_UINT8, 2, vr));
        }
            break;
        case DFNT_INT16:
        {
            int16 vr[2];
            vr[0] = (int16) low;
            vr[1] = (int16) high;
            PTB(sd_setattr(sds_id, "valid_range", DFNT_INT16, 2, vr));
        }
            break;
        case DFNT_UINT16:
        {
            uint16_t vr[2];
            vr[0] = (uint16_t) low;
            vr[1] = (uint16_t) high;
            PTB(sd_setattr(sds_id, "valid_range", DFNT_UINT16, 2, vr));
        }
            break;
        case DFNT_INT32:
        {
            int32_t vr[2];
            vr[0] = (int32_t) low;
            vr[1] = (int32_t) high;
            PTB(sd_setattr(sds_id, "valid_range", DFNT_INT32, 2, vr));
        }
            break;
        case DFNT_UINT32:
        {
            uint32_t vr[2];
            vr[0] = (uint32_t) low;
            vr[1] = (uint32_t) high;
            PTB(sd_setattr(sds_id, "valid_range", DFNT_UINT32, 2, vr));
        }
            break;
        case DFNT_FLOAT32:
        {
            float32 vr[2];
            vr[0] = (float32) low;
            vr[1] = (float32) high;
            PTB(sd_setattr(sds_id, "valid_range", DFNT_FLOAT32, 2, vr));
        }
            break;
        default:
            fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
            fprintf(stderr, "Got unsupported number type (%d) ", nt);
            fprintf(stderr, "while trying to create SDS, \"%s\", ", sname);
            return (PROGRAMMER_BOOBOO);
        }
    }

    /* Add a "slope" attribute if one is specified, and also
    an intercept attribute */
    if (slope != 1.0 || offset != 0.0) {
        PTB(sd_setattr(sds_id, "slope", DFNT_FLOAT, 1, &slope));
        PTB(sd_setattr(sds_id, "intercept", DFNT_FLOAT, 1, &offset));
    }

    /* Add a "units" attribute if one is specified */
    if (units != NULL && *units != 0) {
        PTB(sd_setattr(sds_id, "units", DFNT_CHAR, strlen(units) + 1, units));
    }

    /* Add a "standard_name" attribute if one is specified */
    if (standard_name != NULL && *standard_name != 0) {
        PTB(sd_setattr(sds_id, "standard_name", DFNT_CHAR, strlen(standard_name) + 1, standard_name));
    }

    /* Release this SDS */
    PTB(sd_endaccess(sds_id));

    return (LIFE_IS_GOOD);
}

/************************************************************************
The following functions are just wrappers for the corresponding HDF
functions.
 ************************************************************************/
int sd_setattr(int32_t id, const char *nam, int32_t typ, int32_t cnt, const void* data) {
    if (SDsetattr(id, nam, typ, cnt, data)) {

        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "SDsetattr(%d,\"%s\",%d,%d,data) failed.  ",
                id, nam, typ, cnt);
        return (HDF_FUNCTION_ERROR);
    }
    return (LIFE_IS_GOOD);
}

int sd_create(
        int32_t id,
        const char *nam,
        int32_t typ,
        int32_t rank,
        int32_t d0,
        int32_t d1,
        int32_t d2,
        int32_t *sds_id
        ) {
    int32_t dimsizes[3];
    if (rank > 3) {
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "sd_create() expects to be passed a rank <= 3.  ");
        return (PROGRAMMER_BOOBOO);
    }
    dimsizes[0] = d0;
    dimsizes[1] = d1;
    dimsizes[2] = d2;
    if ((*sds_id = SDcreate(id, nam, typ, rank, dimsizes)) == FAIL) {
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "SDcreate(%d,\"%s\",%d,%d,[%d,%d,%d]) failed.  ",
                id, nam, typ, rank, d0, d1, d2);
        return (HDF_FUNCTION_ERROR);
    }
    return (LIFE_IS_GOOD);
}

int sd_endaccess(int32_t id) {
    if (SDendaccess(id)) {
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "SDendaccess(%d) failed.  ", id);
        return (HDF_FUNCTION_ERROR);
    }
    return (LIFE_IS_GOOD);
}

int sd_setdimnames(int32_t id, const char *d0, const char *d1, const char *d2) {
    PTB(sd_setdimname(id, 0, d0));
    if (d1 != NULL && *d1 != 0) {
        PTB(sd_setdimname(id, 1, d1));
    }
    if (d2 != NULL && *d2 != 0) {
        PTB(sd_setdimname(id, 2, d2));
    }
    return (LIFE_IS_GOOD);
}

int sd_setdimname(int32_t sds_id, int32_t dim_number, const char *name) {
    int32_t dim_id;
    dim_id = SDgetdimid(sds_id, dim_number);
    if (dim_id == FAIL) {
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "SDgetdimid(%d,%d) failed.\n",
                sds_id, dim_number);
        return (HDF_FUNCTION_ERROR);
    }
    if (SDsetdimname(dim_id, name)) {
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "SDsetdimname(%d,\"%s\") failed.\n",
                dim_id, name);
        return (HDF_FUNCTION_ERROR);
    }
    return (LIFE_IS_GOOD);
}

int sd_writedata(
        int32_t sd_id,
        const char *name,
        const void* data,
        int32_t s0,
        int32_t s1,
        int32_t s2,
        int32_t e0,
        int32_t e1,
        int32_t e2
        ) {
    int32_t sds_id, start[3], edge[3];

    PTB(sd_select(sd_id, name, &sds_id));
    start[0] = s0;
    edge[0] = e0;
    start[1] = s1;
    edge[1] = e1;
    start[2] = s2;
    edge[2] = e2;
    if (SDwritedata(sds_id, start, NULL, edge, (void*)data) == FAIL) {
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "SDwritedata(%d,[%d,%d,%d],NULL,[%d,%d,%d],0x%p) ",
                sds_id, s0, s1, s2, e0, e1, e2, data);
        fprintf(stderr, "failed.\n");
        return (HDF_FUNCTION_ERROR);
    }
    PTB(sd_endaccess(sds_id));
    return (LIFE_IS_GOOD);
}

int sd_readdata(
        int32_t sd_id,
        const char *name,
        VOIDP data,
        int32_t s0,
        int32_t s1,
        int32_t s2,
        int32_t e0,
        int32_t e1,
        int32_t e2
        ) {
    int32_t sds_id, start[3], edge[3];

    PTB(sd_select(sd_id, name, &sds_id));
    start[0] = s0;
    edge[0] = e0;
    start[1] = s1;
    edge[1] = e1;
    start[2] = s2;
    edge[2] = e2;
    if (SDreaddata(sds_id, start, NULL, edge, data) == FAIL) {
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "SDreaddata(%d,[%d,%d,%d],NULL,[%d,%d,%d],0x%p) ",
                sds_id, s0, s1, s2, e0, e1, e2, data);
        fprintf(stderr, "failed.\n");
        return (HDF_FUNCTION_ERROR);
    }
    PTB(sd_endaccess(sds_id));
    return (LIFE_IS_GOOD);
}

int sd_select(int32_t sd_id, const char *name, int32_t *sds_id) {

    int32_t index;

    index = SDnametoindex(sd_id, name);
    if (index == FAIL) {
        /*
            fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
            fprintf(stderr,"SDnametoindex(%d,\"%s\") failed.\n",
            sd_id, name);
         */
        return (HDF_FUNCTION_ERROR);
    }
    *sds_id = SDselect(sd_id, index);
    if (*sds_id == FAIL) {
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "SDselect(%d,%d) failed.\n",
                sd_id, index);
        return (HDF_FUNCTION_ERROR);
    }
    return (LIFE_IS_GOOD);
}

/****************************************************************************
Add the named SDS to the Vgroup specified by its Vgroup identifier.
This function uses the global variable sd_id which is set elsewhere
by a call to SDstart().
 *****************************************************************************/
int AddSdsToVgroup(int32_t sd_id, int32_t v_id, const char *name) {

    int32_t sds_ref, sds_id;

    PTB(sd_select(sd_id, name, &sds_id));
    sds_ref = SDidtoref(sds_id);
    if (sds_ref == FAIL) {
        fprintf(stderr,
                "-E- %s line %d: SDidtoref(%d) failed.\n",
                __FILE__, __LINE__, sds_id);
        return (HDF_FUNCTION_ERROR);
    }
    if (Vaddtagref(v_id, DFTAG_NDG, sds_ref) == FAIL) {
        fprintf(stderr,
                "-E- %s line %d: Vaddtagref(%d,%d,%d) failed.\n",
                __FILE__, __LINE__, v_id, DFTAG_SD, sds_ref);
        return (HDF_FUNCTION_ERROR);
    }
    return (LIFE_IS_GOOD);
}

int v_attach(int32_t h_id, int32_t *v_id) {
    *v_id = Vattach(h_id, -1, "w");
    if (*v_id == FAIL) {
        fprintf(stderr,
                "-E- %s line %d: Vattach(%d,-1,\"w\") failed.\n",
                __FILE__, __LINE__, h_id);
        return (HDF_FUNCTION_ERROR);
    }
    return (LIFE_IS_GOOD);
}


/* ------------------------------------------------------ */
/* rdSDS() - reads a SDS (scientific data set), and       */
/*           returns the data                             */
/*                                                        */

/* ------------------------------------------------------ */
int rdSDS(
        int32_t fileID,
        const char *sdsname,
        int32_t start1, /* 1st dimension of starting point */
        int32_t start2, /* 2nd dimension of starting point */
        int32_t edges1, /* 1st dim of length of subset( 0 if reading entire SDS) */
        int32_t edges2, /* 2nd dim of length of subset */
        VOIDP array_data
        ) {
    int32_t sds_id, numtype;
    int32_t sds_index, rank, dims[H4_MAX_VAR_DIMS], nattrs;
    int32_t start[2], edges[2];
    char tmp_sdsname[H4_MAX_NC_NAME];

    start[0] = 0;
    start[1] = 0;
    edges[0] = 0;
    edges[1] = 0;

    /* Get the SDS index */
    sds_index = SDnametoindex(fileID, sdsname);

    if (sds_index < 0) {
        printf("-E- %s: SDS \"%s\" not found.\n", "rdSDS", sdsname);
        SDend(fileID);
        exit(1);
    }


    /* Select the SDS */
    sds_id = SDselect(fileID, sds_index);

    /* Get the rank and number type of the array */
    SDgetinfo(sds_id, tmp_sdsname, &rank, dims, &numtype, &nattrs);

    /*
    if (strcmp(sdsname,"pxl") == 0) {
       ncol = dims[0];
    }
     */

    /* Define the location, pattern and size of the data to read */
    /* set 1st dimension */
    start[0] = start1;
    if (edges1 == 0) {
        edges[0] = dims[0];
    } else {
        edges[0] = edges1;
    }
    /* if rank > 1, set 2nd dimension */
    if (rank > 1) {
        start[1] = start2;
        if (edges2 == 0) {
            edges[1] = dims[1];
        } else {
            edges[1] = edges2;
        }
    }

    /* Based on number type, call the corresponding wrapper
    for the HDF SDreaddata function */
    SDreaddata(sds_id, start, NULL, edges, array_data);

    /* Terminate access to the array */
    SDendaccess(sds_id);

    return (0);
}

/* ------------------------------------------------------ */
/* getDims() - gets a SDS (scientific data set), and      */
/*           returns the dimensions of the data           */
/*                                                        */

/* ------------------------------------------------------ */
int getDims(
        int32_t fileID,
        const char sdsname[H4_MAX_NC_NAME],
        int32_t dims[H4_MAX_VAR_DIMS]
        ) {
    int32_t sds_id, numtype;
    int32_t sds_index, rank, nattrs;
    char tmp_sdsname[H4_MAX_NC_NAME];


    /* Get the SDS index */
    sds_index = SDnametoindex(fileID, sdsname);

    /* Check that the SDS exists */
    if (sds_index == -1) {
        printf("-E- %s:  Error seeking SDS\n", __FILE__);
        return (1);
    }

    /* Select the SDS */
    sds_id = SDselect(fileID, sds_index);

    /* Get the rank and number type of the array */
    SDgetinfo(sds_id, tmp_sdsname, &rank, dims, &numtype, &nattrs);

    /* Terminate access to the array */
    SDendaccess(sds_id);

    return (0);
}


/* ------------------------------------------------------ */
/* get_type() - gets a SDS (scientific data set), and     */
/*           returns the data type                        */
/*                                                        */

/* ------------------------------------------------------ */
int get_type(
        int32_t fileID,
        const char sdsname[H4_MAX_NC_NAME],
        int32_t *dtype
        ) {
    int32_t sds_id;
    int32_t sds_index, rank, nattrs;
    char tmp_sdsname[H4_MAX_NC_NAME];
    int32_t dims[H4_MAX_VAR_DIMS];

    /* Get the SDS index */
    sds_index = SDnametoindex(fileID, sdsname);

    /* Check that the SDS exists */
    if (sds_index == -1) {
        printf("-E- %s:  Error seeking SDS\n", __FILE__);
        return (1);
    }

    /* Select the SDS */
    sds_id = SDselect(fileID, sds_index);

    /* Get the rank and number type of the array */
    SDgetinfo(sds_id, tmp_sdsname, &rank, dims, dtype, &nattrs);

    /* Terminate access to the array */
    SDendaccess(sds_id);

    return (0);
}


/* ------------------------------------------------------ */
/* getHDFattr() - gets an HDF SDS attribute or            */
/*              file attribute (if sdsname is "")         */
/*                                                        */

/* ------------------------------------------------------ */
int getHDFattr(
        int32_t fileID,
        const char *attrname,
        const char *sdsname,
        VOIDP data
        ) {
    int32_t id, attr_index, sds_index;
    int32_t data_type, count;
    char tmp_attrname[H4_MAX_NC_NAME];

    /* get the SDS identifier for SDS attributes */
    if (strcmp(sdsname, "") != 0) {
        /* get the # of the SDS from the SDS name */
        sds_index = SDnametoindex(fileID, sdsname);

        /* Check that the SDS exists */
        if (sds_index == -1) {
            printf("-E- %s:  Error seeking SDS\n", __FILE__);
            return (1);
        }

        id = SDselect(fileID, sds_index);
    } else {
        /* identifier = fileID for file (global) attributes */
        id = fileID;
    }

    /* get the attribute index */
    attr_index = SDfindattr(id, attrname);
    if (attr_index == -1) {
        return (-1);
    }

    /* get the information about the file attribute */
    SDattrinfo(id, attr_index, tmp_attrname,
            &data_type, &count);

    /* read the attribute */
    SDreadattr(id, attr_index, data);

    /* Terminate access to the SDS */
    if (strcmp(sdsname, "") != 0) {
        SDendaccess(id);
    }
    return (0);
}

/*-----------------------------------------------------------------------------
    Function: rdvdata

    Returns: intn (status)

    Description:
        The function rdvdata reads the requested vdata into the given buffer
        and returns the status.

    Arguments: (in calling order)
      Type         Name      I/O     Description
      ----         ----      ---     -----------
      int32        vskey      I      ID number of the vdata
      char *       fields     I      field names of the data to read
      int32        start      I      start element position
      int32        nelt       I      number of elements to read
      uchar *      databuf    O      buffer to read the data

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/11/94    Original development

------------------------------------------------------------------------------*/
intn rdvdata(int32 vskey, const char *fields, int32 start, int32 nelt,
        unsigned char *databuf) {
    int32 ret;

    if ((ret = VSsetfields(vskey, fields)) < 0)
        return FAIL;

    if ((ret = VSseek(vskey, start)) < 0)
        return FAIL;

    if ((ret = VSread(vskey, databuf, nelt, FULL_INTERLACE)) < 0)
        return FAIL;

    return ret;

}

/*-----------------------------------------------------------------------------
    Function: attach_vdata

    Returns: intn (Status)

    Description:
        The function attach_vdata attaches to the requested vdata

    Arguments: (in calling order)
      Type       Name        I/O     Description
      ----       ----        ---     -----------
      int32      fid          I      HDF file ID
      char *     sname        I      vdata name

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/11/94    Original development


------------------------------------------------------------------------------*/
intn attach_vdata(int32 fid, const char *sname) {
    int32 vsid, vskey;

    vsid = VSfind(fid, sname);
    if ((vskey = VSattach(fid, vsid, "r")) < 0)
        return FAIL;
    return vskey;
}

/*-----------------------------------------------------------------------------
    Function: read_SDS

    Returns: intn (status)

    Description:
        The function read_SDS reads the requested SDS/NDG into the given 
        buffer and returns the status.

    Arguments: (in calling order)
      Type         Name      I/O     Description
      ----         ----      ---     -----------
      int32        sdfid      I      ID number 
      char *       sds_name   I      SDS name
      void *       buffer     O      data buffer  

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      06/07/94    Original development

------------------------------------------------------------------------------*/
int32 read_SDS(int32 sdfid, const char *sds_name, void *buffer) {

    int32 index, sdsid, rank, numbertype, nattrs;
    int32 dimsizes[2], start[2];
    char name[255];

    start[0] = start[1] = 0;

    if ((index = SDnametoindex(sdfid, sds_name)) < 0)
        return -2;

    if ((sdsid = SDselect(sdfid, index)) < 0)
        return -2;

    if ((SDgetinfo(sdsid, name, &rank, dimsizes, &numbertype, &nattrs)) < 0)
        return -2;

    if ((SDreaddata(sdsid, start, NULL, dimsizes, buffer)) < 0)
        return -2;

    if ((SDendaccess(sdsid)) < 0)
        return -2;

    return SUCCEED;
}


/* ------------------------------------------------------ */
/* GetFileDesc() - reads the HDF file descriptions into   */
/*                 a buffer.                              */
/*                                                        */

/* ------------------------------------------------------ */
char *GetFileDesc(const char *filename) {
    static char *desc_buffer = NULL;
    int32 file_id, desc_length, fds_len;

    /* Open the file */
    file_id = Hopen(filename, DFACC_READ, 0);

    /* Get the length of the file description */
    desc_length = DFANgetfdslen(file_id, 1);
    if (desc_length == 0)
        return (desc_buffer);

    /* Create a buffer for the file description */
    desc_buffer = HDmalloc(desc_length);
    if (desc_buffer == NULL)
        return (desc_buffer);
    /* Read the file description */
    fds_len = DFANgetfds(file_id, desc_buffer, desc_length, 1);
    if (fds_len == 0)
        return (desc_buffer);

    /* Add ';' if not last non-whitespace character */
    while (desc_buffer[fds_len] <= 32) fds_len--;
    if (desc_buffer[fds_len] != ';') desc_buffer[fds_len + 1] = ';';

    /* Close the file */
    Hclose(file_id);

    return (desc_buffer);
}


