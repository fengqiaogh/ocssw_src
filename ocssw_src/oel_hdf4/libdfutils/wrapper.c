/* =========================================================== */
/* Module wrapper.c                                            */
/*                                                             */
/* HDF/NCDF I/O wrappers                                       */
/*                                                             */
/* Written By:                                                 */
/*     Joel Gales, Futuretech                                  */
/*                                                             */
/* Modification History:                                       */
/*     Joel Gales, Futuretech, OBPG Project, Nov 2013.         */
/*           Add support for CF-compliant metadata             */
/* =========================================================== */

#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h> 
#include <time.h>
#include <math.h>
#include <productInfo.h>

#include <dfutils.h>
#include <genutils.h>
#include <nc4utils.h>

#include <hdf.h>
#include <mfhdf.h>


#define ROUND(x) ((((x) >= 0)?0.5:-0.5) + (x))

int s2u(const char *in, char *out) {
    int i;
    int l = strlen(in);
    for (i = 0; i < l; i++) {
        if (in[i] == ' ')
            out[i] = '_';
        else
            out[i] = in[i];
    }
    out[l] = 0;
    return 0;
}

int8_t findAttr(idDS ds_id, const char *nam) {
    int32_t attr_index = -1;
    if (ds_id.fftype == DS_HDF) {
        if (ds_id.fid > 0) {
            attr_index = SDfindattr(ds_id.fid, nam);
        } else {
            attr_index = SDfindattr(ds_id.sid, nam);
        }

    } else if (ds_id.fftype == DS_NCDF) {
        int status;
        if (ds_id.fid > 0) {
            status = nc_inq_attid(ds_id.fid, NC_GLOBAL, nam, (int32_t *) & attr_index);
        } else {
            status = nc_inq_attid(-ds_id.fid, ds_id.sid, nam, (int32_t *) & attr_index);
        }
        if (status != NC_NOERR) {
            attr_index = -1;
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return attr_index;
}

int readAttr(idDS ds_id, const char *nam, VOIDP data) {
    int status;
    if (ds_id.fftype == DS_HDF) {
        if (ds_id.fid > 0)
            status = getHDFattr(ds_id.fid, nam, "", data);
        else
            status = getHDFattr(ds_id.sid, nam, "", data);
    } else if (ds_id.fftype == DS_NCDF) {
        if (ds_id.fid > 0) {
            status = nc_get_att(ds_id.fid, NC_GLOBAL, nam, data);
        } else {
            status = nc_get_att(-ds_id.fid, ds_id.sid, nam, data);
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    return status;
}

/**
 * Read a string attribute.  Allocating memory for the value which must be freed.
 *
 * @param ds_id data set ID for the file or data set
 * @param name name of the attribute
 * @return pointer to value string or NULL if not found
 */
char* readAttrStr(idDS ds_id, const char *name) {
    char* data;
    int status = 0;
    char buf[128];

    if (ds_id.fftype == DS_HDF) {
        int32_t dtype;
        int32_t count;
        if (ds_id.fid > 0) {
            status = SDattrinfo(ds_id.fid, SDfindattr(ds_id.fid, name), buf, &dtype, &count);
            if (status)
                return NULL;
            if (dtype != DFNT_CHAR)
                return NULL;
            data = (char*) malloc(count + 1);
            if (!data)
                goto memoryError;
            status = getHDFattr(ds_id.fid, name, "", data);
            data[count] = 0; // NULL terminate the string
        } else {
            status = SDattrinfo(ds_id.sid, SDfindattr(ds_id.sid, name), buf, &dtype, &count);
            if (status)
                return NULL;
            if (dtype != DFNT_CHAR)
                return NULL;
            data = (char*) malloc(count + 1);
            if (!data)
                goto memoryError;
            status = getHDFattr(ds_id.sid, name, "", data);
            data[count] = 0; // NULL terminate the string
        }
    } else if (ds_id.fftype == DS_NCDF) {
        int32_t dtype;
        size_t count;
        if (ds_id.fid > 0) {
            status = nc_inq_att(ds_id.fid, NC_GLOBAL, name, &dtype, &count);
            if (status != NC_NOERR)
                return NULL;
            if(dtype == NC_CHAR) {
                data = (char*) malloc(count + 1);
                if (!data)
                    goto memoryError;
                status = nc_get_att(ds_id.fid, NC_GLOBAL, name, data);
                data[count] = 0; // NULL terminate the string
            } else if(dtype == NC_STRING) {
                if(count < 1)
                    return NULL;
                char *list[count];
                status = nc_get_att(ds_id.fid, NC_GLOBAL, name, list);
                data = strdup(list[0]);
                nc_free_string(count, list);
            } else
                return NULL;
        } else {
            status = nc_inq_att(-ds_id.fid, ds_id.sid, name, &dtype, &count);
            if (status != NC_NOERR)
                return NULL;
            if (dtype != NC_STRING && dtype != NC_CHAR)
                return NULL;
            data = (char*) malloc(count + 1);
            if (!data)
                goto memoryError;
            status = nc_get_att(-ds_id.fid, ds_id.sid, name, data);
            data[count] = 0; // NULL terminate the string
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    if (status) {
        free(data);
        return NULL;
    }
    return data;

memoryError:
    printf("-E- %s %d Could not allocate memory for reading attribute %s\n", __FILE__, __LINE__, name);
    exit(1);
}

/**
 * get information about an attribute.
 *
 * @param ds_id data set structure
 * @param nam name of the attribute to query
 * @param dtype datatype of the attribute
 * @param count number of items in the array
 * @return status of 0 if success
 */
int infoAttr(idDS ds_id, const char *nam, int32_t *dtype, int32_t *count) {
    int32_t attr_index;
    int status;
    char buf[128];
    if (ds_id.fftype == DS_HDF) {
        if (ds_id.fid > 0) {
            attr_index = SDfindattr(ds_id.fid, nam);
            status = SDattrinfo(ds_id.fid, attr_index, buf, dtype, count);
        } else {
            attr_index = SDfindattr(ds_id.sid, nam);
            status = SDattrinfo(ds_id.sid, attr_index, buf, dtype, count);
        }
    } else if (ds_id.fftype == DS_NCDF) {
        size_t cnt;
        if (ds_id.fid > 0) {
            status = nc_inq_att(ds_id.fid, NC_GLOBAL, nam, (int32_t *) dtype, &cnt);
            *count = cnt;
        } else {
            status = nc_inq_att(-ds_id.fid, ds_id.sid, nam, (int32_t *) dtype, &cnt);
            *count = cnt;
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return status;
}

int setAttr(idDS ds_id, const char *nam, int32_t typ, int32_t cnt, const void* data) {
    int status;
    if (ds_id.fftype == DS_HDF) {
        status = sd_setattr(ds_id.sid, (char*) nam, typ, cnt, (VOIDP) data);
    } else if (ds_id.fftype == DS_NCDF) {
        status = nc_put_att(ds_id.fid, ds_id.sid, nam, typ, cnt, data);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), nam);
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return status;
}

int SetChrGA(idDS ds_id, const char *name, const char *value) {
    int status;
    if (ds_id.fftype == DS_HDF) {
        status = sd_setattr(ds_id.fid, name, DFNT_CHAR, strlen(value) + 1, value);
    } else if (ds_id.fftype == DS_NCDF) {
        status = nc_put_att_text(ds_id.fid, NC_GLOBAL, name,
                strlen(value) + 1, value);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), name);
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return (status);
}

int SetF64GA(idDS ds_id, const char *name, float64 value) {
    int status;
    if (ds_id.fftype == DS_HDF) {
        status = sd_setattr(ds_id.fid, name, DFNT_FLOAT64, 1, (VOIDP) & value);
    } else if (ds_id.fftype == DS_NCDF) {
        status = nc_put_att_double(ds_id.fid, NC_GLOBAL, name,
                NC_DOUBLE, 1, (const double *) &value);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), name);
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return (status);
}

int SetF32GA(idDS ds_id, const char *name, float32 value) {
    int status;
    if (ds_id.fftype == DS_HDF) {
        status = sd_setattr(ds_id.fid, name, DFNT_FLOAT32, 1, (VOIDP) & value);
    } else if (ds_id.fftype == DS_NCDF) {
        status = nc_put_att_float(ds_id.fid, NC_GLOBAL, name,
                NC_FLOAT, 1, (const float *) &value);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), name);
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return (status);
}

int SetI16GA(idDS ds_id, const char *name, int16 value) {
    int status;
    if (ds_id.fftype == DS_HDF) {
        status = sd_setattr(ds_id.fid, name, DFNT_INT16, 1, (VOIDP) & value);
    } else if (ds_id.fftype == DS_NCDF) {
        status = nc_put_att_short(ds_id.fid, NC_GLOBAL, name,
                NC_SHORT, 1, (const short *) &value);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), name);
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return (status);
}

int SetI8GA(idDS ds_id, const char *name, uint8 value) {
    int status;
    if (ds_id.fftype == DS_HDF) {
        status = sd_setattr(ds_id.fid, name, DFNT_UINT8, 1, (VOIDP) & value);
    } else if (ds_id.fftype == DS_NCDF) {
        status = nc_put_att_uchar(ds_id.fid, NC_GLOBAL, name,
                NC_BYTE, 1, (const uint8 *) &value);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), name);
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return (status);
}

int SetI32GA(idDS ds_id, const char *name, int32_t value) {
    int status;
    if (ds_id.fftype == DS_HDF) {
        status = sd_setattr(ds_id.fid, name, DFNT_INT32, 1, (VOIDP) & value);
    } else if (ds_id.fftype == DS_NCDF) {
        status = nc_put_att_int(ds_id.fid, NC_GLOBAL, name,
                NC_INT, 1, (const int *) &value);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), name);
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return (status);
}

int createDS(idDS ds_id, int sensorId, const char *sname, int32_t dm[3],
        const char dm_name[3][80]) {
    static productInfo_t* p_info;

    if (p_info == NULL) {
        p_info = allocateProductInfo();
    }

    if (findProductInfo(sname, sensorId, p_info)) {
        return createDS2(ds_id, sname, p_info, dm, dm_name);
    } else {
        printf("%s not found in XML product table\n", sname);
        exit(1);
    }
}

int createDS2(idDS ds_id, const char *sname, productInfo_t* p_info, int32_t dm[3],
        const char dm_name[3][80]) {
    int i, status;
    int32_t nt; /* data type */

    if (ds_id.fftype == DS_HDF) {

        if (strcmp(p_info->dataType, "byte") == 0)
            nt = DFNT_INT8;
        else if (strcmp(p_info->dataType, "ubyte") == 0)
            nt = DFNT_UINT8;
        else if (strcmp(p_info->dataType, "short") == 0)
            nt = DFNT_INT16;
        else if (strcmp(p_info->dataType, "ushort") == 0)
            nt = DFNT_UINT16;
        else if (strcmp(p_info->dataType, "int") == 0)
            nt = DFNT_INT32;
        else if (strcmp(p_info->dataType, "uint") == 0)
            nt = DFNT_UINT32;
        else if (strcmp(p_info->dataType, "float") == 0)
            nt = DFNT_FLOAT32;
        else if (strcmp(p_info->dataType, "double") == 0)
            nt = DFNT_FLOAT64;
        else {
            printf("-E- %s %d: datatype %s is not valid\n", __FILE__, __LINE__, p_info->dataType);
            exit(1);
        }

        status = CreateSDS(ds_id.fid, sname, p_info->description, p_info->standardName,
                p_info->units, p_info->validMin, p_info->validMax, p_info->scaleFactor,
                p_info->addOffset, nt, p_info->rank, dm[0], dm[1], dm[2],
                dm_name[0], dm_name[1], dm_name[2]);

    } else if (ds_id.fftype == DS_NCDF) {
        int32_t dimids[3];
        char buf[512];
        for (i = 0; i < p_info->rank; i++) {
            s2u(dm_name[i], buf);
            status = nc_inq_dimid(ds_id.fid, buf, &dimids[i]);
            if (status != NC_NOERR) {
                printf("-E- %s %d: %s for %s\n",
                        __FILE__, __LINE__, nc_strerror(status), dm_name[i]);
                exit(1);
            }
        }

        if (strcmp(p_info->dataType, "byte") == 0)
            nt = NC_BYTE;
        else if (strcmp(p_info->dataType, "ubyte") == 0)
            nt = NC_UBYTE;
        else if (strcmp(p_info->dataType, "short") == 0)
            nt = NC_SHORT;
        else if (strcmp(p_info->dataType, "ushort") == 0)
            nt = NC_USHORT;
        else if (strcmp(p_info->dataType, "int") == 0)
            nt = NC_INT;
        else if (strcmp(p_info->dataType, "uint") == 0)
            nt = NC_UINT;
        else if (strcmp(p_info->dataType, "float") == 0)
            nt = NC_FLOAT;
        else if (strcmp(p_info->dataType, "double") == 0)
            nt = NC_DOUBLE;
        else {
            printf("-E- %s %d datatype %s not defined\n", __FILE__, __LINE__, p_info->dataType);
            exit(1);
        }

        status = CreateNCDF(ds_id, sname, p_info->description, p_info->standardName,
                p_info->reference, p_info->comment, p_info->units, p_info->validMin,
                p_info->validMax, p_info->scaleFactor, p_info->addOffset,
                p_info->fillValue, nt, p_info->rank, dimids);
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return status;
}

int32_t selectDS(idDS ds_id, const char *l2_prod_name) {
    int32_t outid;
    if (ds_id.fftype == DS_HDF) {
        outid = SDselect(ds_id.fid, SDnametoindex(ds_id.fid, l2_prod_name));
    } else if (ds_id.fftype == DS_NCDF) {
        int status;
        status = nc_inq_varid(ds_id.fid, l2_prod_name, (int32_t *) & outid);
        if (status != NC_NOERR) {
            outid = -1;
        }
    } else {
        printf("-E- %s %d fftype %d not defined\n", __FILE__, __LINE__, ds_id.fftype);
        exit(1);
    }
    return outid;
}

int32 checkDS(idDS ds_id, const char *l2_prod_name) {
    int32 outid;
    if (ds_id.fftype == DS_HDF) {
        outid = SDselect(ds_id.fid, SDnametoindex(ds_id.fid, l2_prod_name));
        if (outid != -1) {
            SDendaccess(outid);
        }
    } else if (ds_id.fftype == DS_NCDF) {
        int status;
        status = nc_inq_varid(ds_id.fid, l2_prod_name, (int32_t *) & outid);
        if (status != NC_NOERR) {
            outid = -1;
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return outid;
}

int writeDS(
        idDS ds_id,
        const char *name,
        const void* data,
        int32_t s0,
        int32_t s1,
        int32_t s2,
        int32_t e0,
        int32_t e1,
        int32_t e2) {

    idDS ds_id0 = ds_id;

    if (ds_id0.fftype == DS_HDF) {
        PTB(sd_writedata(ds_id0.fid, name, data, s0, s1, s2, e0, e1, e2));
    } else if (ds_id0.fftype == DS_NCDF) {
        int status = NC_NOERR;
        size_t startp[3] = {s0, s1, s2};
        size_t countp[3] = {e0, e1, e2};
        status = nc_inq_varid(ds_id0.fid, name, (int32_t *) & ds_id0.sid);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), name);
            exit(1);
        }
        status = nc_put_vara(ds_id0.fid, ds_id0.sid, startp, countp, data);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), name);
            exit(1);
        }

    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return 0;
}

int readDS(
        idDS ds_id,
        const char *name,
        int32_t *start,
        int32_t *stride,
        int32_t *count,
        VOIDP data) {

    idDS ds_id0 = ds_id;
    int32_t s0 = start[0], s1 = start[1], s2 = start[2];
    int32_t e0 = count[0], e1 = count[1], e2 = count[2];

    if (ds_id0.fftype == DS_HDF) {
        PTB(sd_readdata(ds_id0.fid, name, data, s0, s1, s2, e0, e1, e2));
    } else if (ds_id0.fftype == DS_NCDF) {
        int status = NC_NOERR;
        size_t startp[3] = {start[0], start[1], start[2]};
        size_t countp[3] = {count[0], count[1], count[2]};
        status = nc_inq_varid(ds_id0.fid, name, (int32_t *) & ds_id0.sid);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), name);
            return status;
        }
        status = nc_get_vara(ds_id0.fid, ds_id0.sid, startp, countp, data);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), name);
            return status;
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return 0;
}

/**
 * Open a data file.
 * Must close the file using endDS()
 *
 * @param filename path to file
 * @param format DS_HDF for HDF4 file, DS_NCDF for netCDF4 file
 * @param access DS_READ for reading, DS_WRITE for writing
 * @param deflate compression factor, 0 is no compression
 * @return data set structure, ds_id.fid = -1 if there was a problem.
 */
idDS startDS(const char *filename, ds_format_t format, ds_access_t access,
        int32_t deflate) {
    idDS ds_id;
    ds_id.fftype = format;
    ds_id.deflate = deflate;

    if (format == DS_HDF) {
        int32_t hdfAccess;
        if (access == DS_READ)
            hdfAccess = DFACC_RDONLY;
        else
            hdfAccess = DFACC_CREATE;
        ds_id.fid = SDstart(filename, hdfAccess);
        if (ds_id.fid == -1) {
            fprintf(stderr, "-E- %s line %d: SDstart failure, %s .\n",
                    __FILE__, __LINE__, filename);
        }
    } else if (format == DS_NCDF) {
        nc_init_chunk_cache();
        int status;
        if (access == DS_READ) {
            status = nc_open(filename, NC_NOWRITE, &ds_id.fid);
            if (status != NC_NOERR) {
                ds_id.fid = -1;
            }
        } else {
            status = nc_create(filename, NC_NETCDF4, &ds_id.fid);
            if (status != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: Could not create NCDF4 file, %s .\n",
                        __FILE__, __LINE__, filename);
                ds_id.fid = -1;
            }
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return ds_id;
}

/*
 * open a data file for reading.  This routine figures out if the
 * file is HDF4 or netCDF4.  Must use endDS() to close the file.
 *
 * @param filename path to file
 * @return a data set structure.  ds_id.fid = -1 if there was a problem.
 */
idDS openDS(const char *filename) {
    ds_format_t fileformat;

    if (Hishdf(filename)) {
        fileformat = DS_HDF;
    } else {
        fileformat = DS_NCDF;
    }
    return startDS(filename, fileformat, DS_READ, 0);
}

int endaccessDS(idDS ds_id) {
    if (ds_id.fftype == DS_HDF) {
        PTB(SDendaccess(ds_id.sid));
    }
    return 0;
}

int endDS(idDS ds_id) {
    if (ds_id.fftype == DS_HDF) {
        return SDend(ds_id.fid);
    } else if (ds_id.fftype == DS_NCDF) {
        return nc_close(ds_id.fid);
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return 0;
}

int getDimsDS(idDS ds_id, const char sdsname[H4_MAX_NC_NAME],
        int32_t dims[H4_MAX_VAR_DIMS]) {

    int status;
    if (ds_id.fftype == DS_HDF) {
        status = getDims(ds_id.fid, sdsname, dims);
    } else if (ds_id.fftype == DS_NCDF) {
        int32_t varid;
        status = nc_inq_varid(ds_id.fid, sdsname, &varid);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), sdsname);
        } else {
            int ndims;
            int dimids[H4_MAX_VAR_DIMS];
            status = nc_inq_var(ds_id.fid, varid, NULL, NULL, &ndims, dimids, NULL);
            if (status != NC_NOERR) {
                printf("-E- %s %d: %s for %s\n",
                        __FILE__, __LINE__, nc_strerror(status), sdsname);
            } else {
                int i;
                size_t dim;
                for (i = 0; i < ndims; i++) {
                    status = nc_inq_dimlen(ds_id.fid, dimids[i], &dim);
                    dims[i] = dim;
                }
            }
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return status;
}

int getTypeDS(idDS ds_id, const char sdsname[H4_MAX_NC_NAME], int32_t *dtype) {
    int status;
    if (ds_id.fftype == DS_HDF) {
        status = get_type(ds_id.fid, sdsname, dtype);
    } else if (ds_id.fftype == DS_NCDF) {
        int32_t varid;
        status = nc_inq_varid(ds_id.fid, sdsname, &varid);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), sdsname);
        } else {
            status = nc_inq_vartype(ds_id.fid, varid, dtype);
            if (status != NC_NOERR) {
                printf("-E- %s %d: %s for %s\n",
                        __FILE__, __LINE__, nc_strerror(status), sdsname);
            }
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return status;
}

int fileInfo(idDS ds_id, int32_t *n_datasets, int32_t *n_globalattr) {
    int status;
    if (ds_id.fftype == DS_HDF) {
        status = SDfileinfo(ds_id.fid, n_datasets, n_globalattr);
    } else if (ds_id.fftype == DS_NCDF) {
        status = nc_inq(ds_id.fid, NULL, n_datasets, n_globalattr, NULL);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s\n", __FILE__, __LINE__, nc_strerror(status));
        }
    } else {
        printf("-E- %s %d fftype not defined\n", __FILE__, __LINE__);
        exit(1);
    }
    return status;
}

int getProdlist(const char *fname, char **prodlist, int32_t *l2_flags_type) {

    size_t i;
    char buffer[2048 * 8];
    int32_t listlen = 0;
    char *cptr;

    if (filesize(fname) == -1) {
        printf("-E- File %s does not exist\n", fname);
        exit(1);
    }

    if (Hishdf(fname) == 1) {

        int32_t HDFfid;
        int32_t vg_ref;
        int32_t vgid;
        int32_t tag;
        int32_t ref;
        int32_t sd_id;
        int32_t sds_id;
        int32_t dims[8];
        int32_t rank;
        int32_t dtype;
        int32_t nattrs;

        HDFfid = Hopen(fname, DFACC_READ, 0);

        if (HDFfid == -1) {
            printf("\n%s not found or is corrupt.\n", fname);
            exit(1);
        }

        sd_id = SDstart(fname, DFACC_RDONLY);
        if (sd_id == -1) {
            printf("Error opening (SDstart) %s\n", fname);
            exit(-1);
        }

        Vstart(HDFfid);
        vg_ref = Vfind(HDFfid, "Geophysical Data");
        vgid = Vattach(HDFfid, vg_ref, "r");

        for (i = 0; i < Vntagrefs(vgid); i++) {
            Vgettagref(vgid, i, &tag, &ref);
            sds_id = SDselect(sd_id, SDreftoindex(sd_id, ref));
            if (sds_id == -1) {
                printf("Error accessing SDS (reference #: %d) in: %s .\n", ref,
                        fname);
                exit(-1);
            }
            SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);
            if (strcmp(buffer, "l2_flags") != 0) {
                listlen += strlen(buffer) + 1;
            }
            SDendaccess(sds_id);
        }

        *prodlist = (char *) calloc(listlen + 1, sizeof (char));
        cptr = *prodlist;

        for (i = 0; i < Vntagrefs(vgid); i++) {
            Vgettagref(vgid, i, &tag, &ref);
            sds_id = SDselect(sd_id, SDreftoindex(sd_id, ref));
            SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);
            if (strcmp(buffer, "l2_flags") != 0) {
                if (i == 0)
                    strcpy(*prodlist, buffer);
                else
                    strcat(*prodlist, buffer);
                strcat(*prodlist, ",");
            } else {
                *l2_flags_type = dtype;
            }
            SDendaccess(sds_id);
        }
        cptr[listlen - 1] = 0;

        SDend(sd_id);
        Vdetach(vgid);
        Vend(HDFfid);

        Hclose(HDFfid);

    } else {

        int ncid, grp_ncid, nvars;
        if (nc_open(fname, NC_NOWRITE, &ncid) == NC_NOERR) {
            nc_inq_ncid(ncid, "geophysical_data", &grp_ncid);
            nc_inq(grp_ncid, NULL, &nvars, NULL, NULL);
            for (i = 0; i < nvars; i++) {
                nc_inq_varname(grp_ncid, i, buffer);
                if (strcmp(buffer, "l2_flags") != 0) {
                    listlen += strlen(buffer) + 1;
                }
            }

            *prodlist = (char *) calloc(listlen + 1, sizeof (char));
            cptr = *prodlist;

            for (i = 0; i < nvars; i++) {
                nc_inq_varname(grp_ncid, i, buffer);
                if (strcmp(buffer, "l2_flags") != 0) {
                    if (i == 0)
                        strcpy(*prodlist, buffer);
                    else
                        strcat(*prodlist, buffer);
                    strcat(*prodlist, ",");
                } else {
                    nc_inq_vartype(grp_ncid, i, l2_flags_type);
                }
            }
        } else {
            printf("-E- Can not open NetCDF file %s\n", fname);
            exit(EXIT_FAILURE);
        }
        cptr[listlen - 1] = 0;
        nc_close(ncid);
    }

    return 0;
}

/* -------------------------------------------------------------------- */
/* Create an NCDF variable using the wrappers for the NCDF routines     */

/* -------------------------------------------------------------------- */
int CreateNCDF(
        idDS ds_id,
        const char *sname, /* short name */
        const char *lname, /* long name */
        const char *standard_name, /* NetCDF standard name (not set if passed NULL or "") */
        const char *reference,
        const char *comment,
        const char *units, /* units (not set if passed NULL or "") */
        double low, /* low end of valid range */
        double high, /* high end of range (no range set if low >= high) */
        float scale_factor, /* scale factor (not set if 1.0)  */
        float add_offset, /* scaling offset (not set if 0.0)  */
        int32_t fillValue, /* fill value */
        int32_t nt, /* NCDF number type */
        int32_t rank, /* number of dimensions (must be <= 3) */
        int32_t dimids[3]/* dimension ids */
        ) {

    int32_t nc_id = ds_id.fid;
    int32_t var_id;
    int i;
    int status;
    size_t chunksize[3] = {0, 0, 0};
    char *validnames[2] = {"valid_min", "valid_max"};

    /* Create the NCDF dataset */
    status = nc_def_var(nc_id, sname, nt, rank, dimids, &var_id);
    if (status != NC_NOERR) {
        printf("-E- %s %d: %s for %s\n",
                __FILE__, __LINE__, nc_strerror(status), sname);
        exit(1);
    }
    if (ds_id.deflate > 0)
        nc_init_compress(nc_id, var_id, dimids, rank, chunksize, (int) ds_id.deflate);

    /* Add a "long_name" attribute */
    status = nc_put_att_text(nc_id, var_id, "long_name", strlen(lname), lname);
    if (status != NC_NOERR) {
        printf("-E- %s %d: %s for %s\n",
                __FILE__, __LINE__, nc_strerror(status), "long_name");
        exit(1);
    }

    /* Add a "scale_factor" attribute and an "add_offset" attribute */
    if (nt != NC_FLOAT && nt != NC_DOUBLE) {
        if (scale_factor != 1.0 || add_offset != 0.0) {
            status = nc_put_att_float(nc_id, var_id, "scale_factor", NC_FLOAT, 1, &scale_factor);
            if (status != NC_NOERR) {
                printf("-E- %s %d: %s for %s\n",
                        __FILE__, __LINE__, nc_strerror(status), "scale_factor");
                exit(1);
            }
            status = nc_put_att_float(nc_id, var_id, "add_offset", NC_FLOAT, 1, &add_offset);
            if (status != NC_NOERR) {
                printf("-E- %s %d: %s for %s\n",
                        __FILE__, __LINE__, nc_strerror(status), "add_offset");
                exit(1);
            }
        }
    }

    /* Add a "units" attribute if one is specified */
    if (units != NULL
            && units[0] != 0
            && strcasecmp(units, "dimensionless") != 0
            && strcasecmp(units, "unitless") != 0) {
        status = nc_put_att_text(nc_id, var_id, "units", strlen(units), units);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), "units");
            exit(1);
        }
    }

    /* add coordinates if rank > 1 */
    if (rank > 1 && strcasecmp(sname, "longitude") != 0 && strcasecmp(sname, "latitude") != 0) {
        status = nc_put_att_text(nc_id, var_id, "coordinates", 18, "longitude latitude");
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), "coordinates");
            exit(1);
        }
    }

    /* Add a "standard_name" attribute */
    if (standard_name != NULL && strcmp(standard_name, "") != 0) {
        status = nc_put_att_text(nc_id, var_id, "standard_name",
                strlen(standard_name), standard_name);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), "standard_name");
            exit(1);
        }
    }

    if (strstr(sname, "flag_") == sname ||
            strstr(sname, "l2_flags") != NULL)
        goto skip_fill;

    /* Add a "_FillValue" attribute */
    switch (nt) { /* Use the appropriate number type */
    case NC_BYTE:
    {
        uint8_t fv_i8 = 255;
        status = nc_put_att(nc_id, var_id, "_FillValue", nt, 1, &fv_i8);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), "_FillValue");
            exit(1);
        }
        break;
    }
    case NC_UBYTE:
    {
        uint8_t fv_i8 = -128;
        status = nc_put_att(nc_id, var_id, "_FillValue", nt, 1, &fv_i8);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), "_FillValue");
            exit(1);
        }
        break;
    }
    case NC_USHORT:
    {
        uint16_t fv_i16 = (uint16_t) fillValue;
        status = nc_put_att(nc_id, var_id, "_FillValue", nt, 1, &fv_i16);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), "_FillValue");
            exit(1);
        }
        break;
    }
    case NC_SHORT:
    {
        int16_t fv_i16 = (int16_t) fillValue;
        status = nc_put_att(nc_id, var_id, "_FillValue", nt, 1, &fv_i16);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), "_FillValue");
            exit(1);
        }
        break;
    }
    case NC_INT:
    {
        status = nc_put_att(nc_id, var_id, "_FillValue", nt, 1, &fillValue);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), "_FillValue");
            exit(1);
        }
        break;
    }
    case NC_UINT:
    {
        status = nc_put_att(nc_id, var_id, "_FillValue", nt, 1, &fillValue);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), "_FillValue");
            exit(1);
        }
        break;
    }
    case NC_FLOAT:
    {
        float fv_f32 = (float) fillValue;
        status = nc_put_att(nc_id, var_id, "_FillValue", nt, 1, &fv_f32);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), "_FillValue");
            exit(1);
        }
        break;
    }
    case NC_DOUBLE:
    {
        double fv_f64 = (double) fillValue;
        status = nc_put_att(nc_id, var_id, "_FillValue", nt, 1, &fv_f64);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), "_FillValue");
            exit(1);
        }
        break;
    }
    default:
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "Got unsupported fill values number type (%d) ", nt);
        fprintf(stderr, "while trying to create NCDF variable, \"%s\", ", sname);
        return (PROGRAMMER_BOOBOO);
    }

skip_fill:
    /* Add "valid_min" & "valid_max" attributes if low is less than high
     * and they are not default values
     * (i.e. they were properly defined in the product.xml)
     */
    if ((low < high) && (low != PRODUCT_DEFAULT_validMin || high != PRODUCT_DEFAULT_validMax)) {
        switch (nt) { /* Use the appropriate number type */
        case NC_BYTE:
        {
            int8_t vr[2];
            vr[0] = (int8_t) ROUND((low - add_offset) / scale_factor);
            vr[1] = (int8_t) ROUND((high - add_offset) / scale_factor);
            for (i = 0; i < 2; i++) {
                status = nc_put_att_schar(nc_id, var_id, validnames[i], NC_BYTE,
                        1, (const signed char *) &vr[i]);
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), validnames[i]);
                    exit(1);
                }
            }
        }
            break;
        case NC_UBYTE:
        {
            uint8_t vr[2];
            vr[0] = (uint8_t) ROUND((low - add_offset) / scale_factor);
            vr[1] = (uint8_t) ROUND((high - add_offset) / scale_factor);
            for (i = 0; i < 2; i++) {
                status = nc_put_att_uchar(nc_id, var_id, validnames[i], NC_UBYTE,
                        1, (const unsigned char *) &vr[i]);
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), validnames[i]);
                    exit(1);
                }
            }
        }
            break;
        case NC_SHORT:
        {
            int16_t vr[2];
            vr[0] = (int16_t) ROUND((low - add_offset) / scale_factor);
            vr[1] = (int16_t) ROUND((high - add_offset) / scale_factor);
            for (i = 0; i < 2; i++) {
                status = nc_put_att_short(nc_id, var_id, validnames[i], NC_SHORT,
                        1, &vr[i]);
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), validnames[i]);
                    exit(1);
                }
            }
        }
            break;
        case NC_USHORT:
        {
            uint16_t vr[2];
            vr[0] = (uint16_t) ROUND((low - add_offset) / scale_factor);
            vr[1] = (uint16_t) ROUND((high - add_offset) / scale_factor);
            for (i = 0; i < 2; i++) {
                status = nc_put_att_ushort(nc_id, var_id, validnames[i], NC_USHORT,
                        1, &vr[i]);
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), validnames[i]);
                    exit(1);
                }
            }
        }
            break;
        case NC_INT:
        {
            int32_t vr[2];
            vr[0] = (int32_t) ROUND((low - add_offset) / scale_factor);
            vr[1] = (int32_t) ROUND((high - add_offset) / scale_factor);
            for (i = 0; i < 2; i++) {
                status = nc_put_att_int(nc_id, var_id, validnames[i], NC_INT,
                        1, &vr[i]);
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), validnames[i]);
                    exit(1);
                }
            }
        }
            break;
        case NC_UINT:
        {
            uint32_t vr[2];
            vr[0] = (uint32_t) ROUND((low - add_offset) / scale_factor);
            vr[1] = (uint32_t) ROUND((high - add_offset) / scale_factor);
            for (i = 0; i < 2; i++) {
                status = nc_put_att_uint(nc_id, var_id, validnames[i], NC_UINT,
                        1, &vr[i]);
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), validnames[i]);
                    exit(1);
                }
            }
        }
            break;
        case NC_FLOAT:
        {
            float vr[2];
            vr[0] = (float) low;
            vr[1] = (float) high;
            for (i = 0; i < 2; i++) {
                status = nc_put_att_float(nc_id, var_id, validnames[i], NC_FLOAT,
                        1, &vr[i]);
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), validnames[i]);
                    exit(1);
                }
            }
        }
            break;
        default:
            fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
            fprintf(stderr, "Got unsupported number type (%d) ", nt);
            fprintf(stderr, "while trying to create NCDF variable, \"%s\", ", sname);
            return (PROGRAMMER_BOOBOO);
        }
    }

    /* Add a "reference" attribute */
    if (reference != NULL && strcmp(reference, "") != 0) {
        status = nc_put_att_text(nc_id, var_id, "reference",
                strlen(reference), reference);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), "reference");
            exit(1);
        }
    }

    /* Add a "comment" attribute */
    if (comment != NULL && strcmp(comment, "") != 0) {
        status = nc_put_att_text(nc_id, var_id, "comment",
                strlen(comment), comment);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), "comment");
            exit(1);
        }
    }

    return (LIFE_IS_GOOD);
}
