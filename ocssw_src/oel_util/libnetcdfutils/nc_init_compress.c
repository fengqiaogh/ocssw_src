#include <netcdf.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <nc4utils.h>

#define DEFAULT_CHUNK_SIZE 1000

#define NEW_CACHE_SIZE 16000000
#define NEW_CACHE_NELEMS 2003
#define NEW_CACHE_PREEMPTION .75 

void nc_init_chunk_cache() {
    /* Change chunk cache. */
    int status = nc_set_chunk_cache(NEW_CACHE_SIZE, NEW_CACHE_NELEMS, NEW_CACHE_PREEMPTION);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: Could not set NetCDF4 cache size.\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
}

/**
 * nc_init_compress
 *
 *	R. Healy 9/26/2016
 *
 * 	@param[in]nc_id         - int32_t  netcdf file (or group) id
 * 	@param[in]var_id        - int32_t  netcdf variable id
 * 	@param[in]dimids        - int32_t* netcdf dimension ids for variable with var_id
 * 	@param[in]rank          - int32_t  rank (number) of dimensions
 * 	@param[in]chunksize     - size_t*  chunk size for compression for each dimension with dimid
 * 	@param[in]deflate_level - int      deflation level 1=lowest(fastest), 9=highest(slowest)
 *
 * Note: chunksize needs to be positive non 0 for unlimited dimensions
 */
void nc_init_compress(int32_t nc_id, int32_t var_id, int32_t *dimids, int32_t rank,
        size_t *chunksize, int deflate_level) {
    int i, status;

    // if deflate level below 0 then don't deflate
    if (deflate_level < 1)
        return;

    size_t dimlength;
    if (rank > 3) {
        printf("Whoops!  Refusing to chunk/compress more than 3 dimensions");
        exit(EXIT_FAILURE);
    }
    size_t suggested_size[3] = {32,256,40};
    if (rank == 2) {
        suggested_size[0] = 256;
        suggested_size[1] = 2048;
    }
    // Set  chunksize for each dimension, if one has not been provided
    for (i = 0; i < rank; i++) {
        if (chunksize==NULL || chunksize[i]==0) {
            chunksize[i] = suggested_size[i];
        }
        status = nc_inq_dimlen(nc_id, dimids[i], &dimlength);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d: Could not read size of dimension.\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        if(chunksize[i] > dimlength)
            chunksize[i] = dimlength;
    }

    /* Set compression */
    /* First set chunking */
    status = nc_def_var_chunking(nc_id, var_id, NC_CHUNKED, chunksize);
    if (status != NC_NOERR) {
        printf("-E- %s %d: %s \n", __FILE__, __LINE__,
                nc_strerror(status));
        exit(EXIT_FAILURE);
    }
    /* Now we can set compression */
    status = nc_def_var_deflate(nc_id, var_id, NC_SHUFFLE, 9,
            deflate_level);
    if (status != NC_NOERR) {
        printf("-E- %s %d: %s \n", __FILE__, __LINE__,
                nc_strerror(status));
        exit(EXIT_FAILURE);
    }
}

/* Check a set of chunksizes to see if they specify a chunk that is too big. */
int check_chunksizes(size_t type_len, int32_t ndims, const size_t *chunksizes) {
    double dprod;
    int d;

    dprod = (double) type_len;
    for (d = 0; d < ndims; d++) {
        if (chunksizes[d] < 1)
            return NC_EINVAL;
        dprod *= (double) chunksizes[d];
    }

    if (dprod > (double) NC_MAX_UINT)
        return NC_EBADCHUNK;

    return NC_NOERR;
}
