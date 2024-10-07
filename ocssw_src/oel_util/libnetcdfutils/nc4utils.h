#ifndef _NC4UTILS_
#define _NC4UTILS_

#include <netcdf.h>

#ifdef __cplusplus
extern "C" {
#endif


#include <stddef.h>
#include <stdint.h>
//#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <math.h>
//#include <ctype.h>

/*----- Utility macros -----*/

#define TRYMEM(file,line,memstat) {   \
        if (memstat == NULL) {    \
            fprintf(stderr, "-E- %s:%d: "  \
                    "Memory allocation error.\n", \
                    file, line);   \
            fflush(stderr);    \
            exit(1); }     \
    }

#define TRY_NC(file,line,ncstat) {   \
        if (ncstat != NC_NOERR) {   \
            fprintf(stderr, "-E- %s:%d: "  \
                    "%s\n",    \
                    file, line, nc_strerror(ncstat)); \
            fflush(stderr);    \
            exit(1); }     \
    }

/*----- Datatype definitions -----*/

typedef struct dim_str_nc {
    char name[NC_MAX_NAME + 1]; /* dimension name */
    size_t len; /* dimension length */
} dim_str_nc;

typedef struct att_str_nc {
    char name[NC_MAX_NAME + 1]; /* attribute name */
    nc_type type; /* attribute type */
    size_t nvals; /* number of attribute values */
    void *value; /* pointer to untyped attribute value(s) */
} att_str_nc;

typedef struct var_str_nc {
    int grpid; /* parent group ID */
    int id; /* variable index */
    char name[NC_MAX_NAME + 1]; /* variable name */
    nc_type type; /* variable type */
    int ndims; /* number  of dimensions */
    dim_str_nc *dim; /* pointer to dimensions */
    int natts; /* number  of attributes */
    att_str_nc *att; /* pointer to attributes */
    int dimids[NC_MAX_VAR_DIMS]; /* array of dimension IDs */
    void *data; /* pointer to data */
} var_str_nc;

typedef struct grp_str_nc {
    int id; /* group ID */
    char name[NC_MAX_NAME + 1]; /* group name */
    int ndims; /* number  of dimensions */
    dim_str_nc *dim; /* pointer to dimensions */
    int natts; /* number  of attributes */
    att_str_nc *att; /* pointer to attributes */
    int nvars; /* number  of variables */
    var_str_nc *var; /* pointer to variables */
    int ngrps; /* number  of subgroups */
    struct grp_str_nc *grps; /* pointer to subgroups */
} grp_str_nc;

/*----- Function Prototypes----- */

dim_str_nc *load_dims_nc(int ncid);
void print_dims_nc(dim_str_nc *dim, int ndims);
void free_dims_nc(dim_str_nc *dim);
int find_dimid_nc(int ncid, int *dimid, const char *dimnames[]);

att_str_nc *load_atts_nc(int ncid, int varid);
void print_atts_nc(att_str_nc *att, int natts, char *varname);
void free_atts_nc(att_str_nc *att, int natts);

var_str_nc *load_vars_nc(int ncid);
var_str_nc* find_var_byname_nc(grp_str_nc nc,
        const char* varname, const char* grpname);
int readall_var(var_str_nc *var);
void print_vars_nc(var_str_nc *var, int nvars);
void free_vars_nc(var_str_nc *var, int nvars);
int find_varid_nc(int ncid, int *varid, const char *varnames[]);

int load_grp_nc(grp_str_nc *grp);
void print_grp_nc(grp_str_nc grp);
void free_grp_nc(grp_str_nc *grp);
void check_err(const int stat, const int line, const char *file);
void report_err(const int stat, const int line, const char *file);
int createNCDF(int ncid, const char *sname, const char *lname,
        const char *standard_name, const char *units,
        void *fill_value,
        const char *flag_values, const char *flag_meanings,
        double low, double high, int nt,
        int rank, int *dimids, size_t *chunksize);

/**
 * init the netCDF reading cache sizes to something reasonable
 */
void nc_init_chunk_cache();
void nc_init_compress(int32_t nc_id, int32_t var_id, int32_t *dimids, int32_t rank,
        size_t *chunksize, int deflate_level);
int check_chunksizes(size_t type_len, int32_t ndims, const size_t *chunksizes);
#ifdef __cplusplus
}
#endif

#endif /* _NC_UTILS_ */
