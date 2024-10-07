#ifndef SEAPROTO_H
#define SEAPROTO_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*  Prototypes  */
int32_t get_beg_ext32(int32_t n_bins_write, int32_t *binnum_data,
        int32_t *basebin, int32_t nrows, int32_t *beg, int32_t *ext);

int64_t get_beg_ext(int32_t n_bins_write, int64_t *binnum_data,
        int64_t *basebin, int32_t nrows, int64_t *beg, int32_t *ext);

int32_t wr_vdata(char *outname, int32_t fileid_w, int32_t vgid, char *name,
        char *class1, int32_t n_flds, int32_t n_recs_to_write, char *fldname[],
        int32_t type[], int32_t noext, uint8_t *data, int32_t verbose);


#ifdef __cplusplus
}
#endif

#endif /* SEAPROTO_H */
