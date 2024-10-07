#include "nc4utils.h"
#ifndef FAIL
#define FAIL 1
#endif
#ifndef SUCCESS
#define SUCCESS 0
#endif

/* ------------------------------------------------------------------ */
/* Utilities */

static const char* typename[] = {/* from netcdf.h */
    "notype", /*  0 = NC_NAT    : Not A Type */
    "byte", /*  1 = NC_BYTE   : signed 1-byte integer */
    "char", /*  2 = NC_CHAR   : ISO/ASCII character */
    "short", /*  3 = NC_SHORT  : signed 2-byte integer */
    "int", /*  4 = NC_INT    : signed 4-byte integer */
    "float", /*  5 = NC_FLOAT  : single precision floating point number */
    "double", /*  6 = NC_DOUBLE : double precision floating point number */
    "ubyte", /*  7 = NC_UBYTE  : unsigned 1-byte integer */
    "ushort", /*  8 = NC_USHORT : unsigned 2-byte integer */
    "uint", /*  9 = NC_UINT   : unsigned 4-byte integer */
    "int64", /* 10 = NC_INT64  : signed 8-byte integer */
    "uint64", /* 11 = NC_UINT64 : unsigned 8-byte integer */
    "string" /* 12 = NC_STRING : string */
};

void report_err(const int stat, const int line, const char *file) {
    if (stat != NC_NOERR) {
        (void) fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
        fflush(stderr);
    }
}

void check_err(const int stat, const int line, const char *file) {
    if (stat != NC_NOERR) {
        (void) fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
        fflush(stderr);
        exit(1);
    }
}

char* format_ncval(void *value, nc_type type, size_t nvals) {
    int i, nchar;
    char *result;
    char tmp[256];

    /*
      TBD: prevent buffer overflows.
      how much memory to allocate for tmp string?
      NC_CHAR:        nvals + "" + null = nvals+3
      how about numerical types?
     */

    if (type == NC_CHAR) {
        ((char *) value)[nvals] = '\0'; // null terminate
        nchar = nvals + 3; //add space for quotes
        TRYMEM(__FILE__, __LINE__,
                (result = malloc(nchar * sizeof (char))));
        result[0] = '\0';
        sprintf(result, "\"%s\"", (char *) value);

    } else {
        TRYMEM(__FILE__, __LINE__,
                (result = malloc(64 * nvals)));
        result[0] = '\0';
        for (i = 0; i < nvals; i++) {
            switch (type) {
            case NC_BYTE:
                sprintf(tmp, "\%db", ((signed char *) value)[i]);
                break;
            case NC_SHORT:
                sprintf(tmp, "\%ds", ((short *) value)[i]);
                break;
            case NC_INT:
                sprintf(tmp, "\%d", ((int *) value)[i]);
                break;
            case NC_FLOAT:
                sprintf(tmp, "\%gf", ((float *) value)[i]);
                break;
            case NC_DOUBLE:
                sprintf(tmp, "\%g", ((double *) value)[i]);
                break;
            case NC_UBYTE:
                sprintf(tmp, "\%uUB", ((unsigned char *) value)[i]);
                break;
            case NC_USHORT:
                sprintf(tmp, "\%huUS", ((unsigned short *) value)[i]);
                break;
            case NC_UINT:
                sprintf(tmp, "\%uU", ((unsigned int *) value)[i]);
                break;
            default:
                sprintf(tmp, "\tvalue of unknown type %d, nvals %d)",
                        (int) type, (int) nvals);
            } // TO DO: handle type NC_STRING
            strcat(result, tmp);
            if (i < nvals - 1) strcat(result, ", ");
        }
    }
    return result;
}

/* ------------------------------------------------------------------ */

/* Dimensions */


dim_str_nc *load_grpdims(int ncid) {
    int ndims = 0;
    dim_str_nc *dim = NULL;

    /* find number of dimensions visible to group */
    TRY_NC(__FILE__, __LINE__,
            nc_inq_dimids(ncid, &ndims, NULL, 1));

    /* load dimension info */
    if (ndims > 0) {
        int i; /* dimension index */
        int *dimids; /* dimension ID */

        TRYMEM(__FILE__, __LINE__,
                (dimids = malloc(ndims * sizeof (*dimids))));
        TRY_NC(__FILE__, __LINE__,
                nc_inq_dimids(ncid, NULL, dimids, 1)); // ids for each dimension

        TRYMEM(__FILE__, __LINE__,
                (dim = malloc(ndims * sizeof (*dim))));
        for (i = 0; i < ndims; i++)
            TRY_NC(__FILE__, __LINE__,
                nc_inq_dim(ncid, dimids[i], dim[i].name, &dim[i].len));

        free(dimids);
        dimids = NULL;
    }
    return dim;
}

void print_dims_nc(dim_str_nc *dim, int ndims) {
    int i; /* dimension index */
    for (i = 0; i < ndims; i++)
        printf("\t%s = %d ;\n",
            dim[i].name,
            (int) dim[i].len);
}

void free_dims_nc(dim_str_nc *dim) {
    if (dim) {
        free(dim);
        dim = NULL;
    }
}

/**
   Find ID of first dimension found in input list.
   @param[in] ncid NetCDF file or group ID,
   @param[in] dimnames Null-terminated list of possible dimension names.
   @param[out] dimid NetCDF dimension ID
   @return Error if no matching dimension found.
 */
int find_dimid_nc(int ncid, int *dimid, const char *dimnames[]) {
    int i;
    int status; /* error status */
    i = 0;
    while (dimnames[i] != NULL) {
        status = nc_inq_dimid(ncid, dimnames[i], dimid);
        if (status == NC_NOERR) break;
        i++;
    }
    return status;
}

/* ------------------------------------------------------------------ */

/* Attributes */

att_str_nc *load_atts_nc(int ncid, int varid) {
    int natts = 0;
    att_str_nc *att = NULL;

    /* find number of attributes for this group or variable */
    TRY_NC(__FILE__, __LINE__,
            nc_inq_varnatts(ncid, varid, &natts));

    if (natts > 0) {
        int i; /* attribute index */
        size_t typesize, nvals;

        /* initialize attribute structure */
        TRYMEM(__FILE__, __LINE__,
                (att = malloc(natts * sizeof (*att))));

        for (i = 0; i < natts; i++) {

            /* load info for this attribute */
            TRY_NC(__FILE__, __LINE__,
                    nc_inq_attname(ncid, varid, i, att[i].name));
            TRY_NC(__FILE__, __LINE__,
                    nc_inq_att(ncid, varid,
                    att[i].name,
                    &att[i].type,
                    &att[i].nvals));

            /* allocate space to read in all values of this attribute */
            TRY_NC(__FILE__, __LINE__,
                    nc_inq_type(ncid, att[i].type, NULL, &typesize));
            nvals = (att[i].type == NC_CHAR) ? att[i].nvals + 1 : att[i].nvals;
            TRYMEM(__FILE__, __LINE__,
                    (att[i].value = malloc(nvals * typesize)));

            /* read attribute values */
            TRY_NC(__FILE__, __LINE__,
                    nc_get_att(ncid, varid,
                    att[i].name,
                    att[i].value));
            if (att[i].type == NC_CHAR) // set terminating char
                ((char*) att[i].value)[att[i].nvals] = '\0';

        }
    }
    return att;
}

void print_atts_nc(att_str_nc *att, int natts, char *varname) {
    int i; /* attribute index */
    for (i = 0; i < natts; i++)
        if (att[i].nvals > 0)
            printf("\t\t%s:%s = %s ;\n",
                varname, att[i].name,
                format_ncval(att[i].value,
                att[i].type,
                att[i].nvals)
                );
}

void free_atts_nc(att_str_nc *att, int natts) {
    int i; /* attribute index */
    if (att) {
        for (i = 0; i < natts; i++) {
            if (att[i].value) {
                free(att[i].value);
                att[i].value = NULL;
            }
        }
        free(att);
        att = NULL;
    }
}

/* ------------------------------------------------------------------ */

/* Variables */

var_str_nc *load_vars_nc(int ncid) {
    int nvars = 0;
    var_str_nc *var = NULL;

    /* find number of variables inside group */
    TRY_NC(__FILE__, __LINE__,
            nc_inq_varids(ncid, &nvars, NULL)); // # vars inside group

    if (nvars > 0) {
        int i; /* variable index */
        int *varids; /* variable ID */

        /* initialize variable structure */
        TRYMEM(__FILE__, __LINE__,
                (varids = malloc(nvars * sizeof (*varids))));
        TRY_NC(__FILE__, __LINE__,
                nc_inq_varids(ncid, NULL, varids)); // ids for each variable
        TRYMEM(__FILE__, __LINE__,
                (var = malloc(nvars * sizeof (*var))));

        for (i = 0; i < nvars; i++) {
            var[i].grpid = ncid;
            var[i].id = varids[i];
            var[i].data = NULL;
            var[i].dim = NULL;

            /* load info for this variable */
            TRY_NC(__FILE__, __LINE__,
                    nc_inq_var(ncid, varids[i],
                    var[i].name,
                    &var[i].type,
                    &var[i].ndims,
                    var[i].dimids,
                    &var[i].natts));

            /* load dimensions */
            if (var[i].ndims > 0) {
                int idim; /* dimension index */
                dim_str_nc *dim = NULL;
                TRYMEM(__FILE__, __LINE__,
                        (dim = malloc(var[i].ndims * sizeof (*dim))));
                for (idim = 0; idim < var[i].ndims; idim++) {
                    TRY_NC(__FILE__, __LINE__,
                            nc_inq_dim(ncid,
                            var[i].dimids[idim],
                            dim[idim].name,
                            &dim[idim].len));
                }
                var[i].dim = dim;
            }
            /* load variable attributes */
            var[i].att = load_atts_nc(ncid, varids[i]);
        }
    }
    return var;
}

/* Returns pointer to specified variable in already-loaded group structure */

/* Works only for 1st-level groups */
var_str_nc* find_var_byname_nc(grp_str_nc nc,
        const char* varname, const char* grpname) {
    int32_t igrp, ivar;
    var_str_nc *varptr = NULL;

    for (igrp = 0; igrp < nc.ngrps; igrp++)
        if (strcmp(grpname, (char*) nc.grps[igrp].name) == 0)
            break;

    if (igrp < nc.ngrps) {
        for (ivar = 0; ivar < nc.grps[igrp].nvars; ivar++)
            if (strcmp(varname, (char*) nc.grps[igrp].var[ivar].name) == 0)
                break;
        if (ivar < nc.grps[igrp].nvars)
            varptr = &nc.grps[igrp].var[ivar];
    }

    return varptr;
}

int readall_var(var_str_nc *var) {
    int32_t status = SUCCESS;

    if (!var) {
        printf("-E- readall_var - Hey, you passed in a NULL pointer.\n");
        exit(EXIT_FAILURE);
    }

    if (var->ndims > 0) {
        int i; /* dimension index */
        size_t typesize;
        size_t nvals = 1;

        /* allocate space to read in all values of this variable */
        TRY_NC(__FILE__, __LINE__,
                nc_inq_type(var->id, var->type, NULL, &typesize));
        for (i = 0; i < var->ndims; i++)
            nvals *= var->dim[i].len;
        TRYMEM(__FILE__, __LINE__,
                (var->data = malloc(nvals * typesize)));

        /* read variable values */
        TRY_NC(__FILE__, __LINE__,
                nc_get_var(var->grpid, var->id, var->data));
    }
    return status;
}

void print_vars_nc(var_str_nc *var, int nvars) {
    int i; /* variable index */
    var_str_nc v;
    int idim;

    for (i = 0; i < nvars; i++) {
        v = var[i];
        printf("\t%s %s", typename[v.type], v.name);

        /* dimension array */
        if (v.ndims > 0)
            printf("(");
        for (idim = 0; idim < v.ndims; idim++) {
            printf("%s", v.dim[idim].name);
            printf("%s", idim < v.ndims - 1 ? ", " : ")");
        }

        printf(" ;\n");
        print_atts_nc(v.att, v.natts, v.name);
    }
}

void free_vars_nc(var_str_nc *var, int nvars) {
    int i; /* variable index */
    if (var) {
        for (i = 0; i < nvars; i++) {
            free_dims_nc(var[i].dim); // dimensions
            free_atts_nc(var[i].att, var[i].natts); // attributes
            if (var[i].data) { // data
                free(var[i].data);
                var[i].data = NULL;
            }
        }
        free(var);
        var = NULL;
    }
}

/**
   Find ID of first variable found in input list.
   @param[in] ncid NetCDF file or group ID,
   @param[in] varnames Null-terminated list of possible variable names.
   @param[out] varid NetCDF variable ID
   @return Error if no variable found.
 */
int find_varid_nc(int ncid, int *varid, const char *varnames[]) {
    int i;
    int status; /* error status */
    i = 0;
    while (varnames[i] != NULL) {
        status = nc_inq_varid(ncid, varnames[i], varid);
        if (status == NC_NOERR) break;
        i++;
    }
    return status;
}

/* ------------------------------------------------------------------ */

/* Group (including root level) */

int load_grp_nc(grp_str_nc *grp) {
    int status = SUCCESS; /* error status */
    int ncid = grp->id;

    /* load group info */
    TRY_NC(__FILE__, __LINE__,
            nc_inq_grpname(ncid, grp->name));
    TRY_NC(__FILE__, __LINE__,
            nc_inq(ncid,
            &grp->ndims,
            &grp->nvars,
            &grp->natts,
            NULL)); // don't bother with unlimited dimensions

    /* load dimensions, variables & group attributes */
    grp->dim = load_grpdims(ncid);
    grp->att = load_atts_nc(ncid, NC_GLOBAL);
    grp->var = load_vars_nc(ncid);

    /* load any subgroups */
    TRY_NC(__FILE__, __LINE__,
            nc_inq_grps(grp->id, &grp->ngrps, NULL)); // number of subgroups
    if (grp->ngrps > 0) {
        int j; /* subgroup index */
        int *grpids; /* subgroup ID */

        TRYMEM(__FILE__, __LINE__,
                (grpids = malloc(grp->ngrps * sizeof (*grpids))));
        TRY_NC(__FILE__, __LINE__,
                nc_inq_grps(ncid, NULL, grpids)); // ids for each subgroup

        TRYMEM(__FILE__, __LINE__,
                (grp->grps = malloc(grp->ngrps * sizeof (*grp->grps))));
        for (j = 0; j < grp->ngrps; j++) {
            grp->grps[j].id = grpids[j];
            load_grp_nc(&grp->grps[j]); // recursive
        }
        free(grpids);
    }

    return status;
}

void print_grp_nc(grp_str_nc grp) {
    int i;
    //static int level = 0;
    printf("\ngroup: %s {\n", grp.name);

    /* dimensions */
    if (grp.ndims > 0) printf("dimensions:\n");
    print_dims_nc(grp.dim, grp.ndims);

    /* variables */
    if (grp.nvars > 0) printf("variables:\n");
    print_vars_nc(grp.var, grp.nvars);

    /* group attributes */
    if (grp.natts > 0) printf("\n// attributes:\n");
    print_atts_nc(grp.att, grp.natts, "");

    /* subgroups */
    for (i = 0; i < grp.ngrps; i++) {
        //level++;
        print_grp_nc(grp.grps[i]); // recursive
        //level--;
    }

    printf("} // group %s\n", grp.name);
    printf("\n");
    fflush(stdout); // make sure everything has printed
}

void free_grp_nc(grp_str_nc *grp) {
    int i;

    if (grp) {

        /* dimensions */
        free_dims_nc(grp->dim);

        /* variables */
        free_vars_nc(grp->var, grp->nvars);

        /* group attributes */
        free_atts_nc(grp->att, grp->natts);

        /* subgroups */
        for (i = 0; i < grp->ngrps; i++)
            free_grp_nc(&grp->grps[i]); // recursive

        /* current group */
        if (grp->grps) {
            free(grp->grps);
            grp->grps = NULL;
        }
    }

}

/* ------------------------------------------------------------------ */
