#ifndef L2EXTRACT_3D_WV_H
#define L2EXTRACT_3D_WV_H

#ifdef __cplusplus
extern "C" {
#endif
void getWave3Indexes(const char *inp_file, const char *wv_list, int **wv_vals_out, int **wv_indexes_out,
                     int *wv_num_to_pass, const char *prod_list);
#ifdef __cplusplus
}
#endif

#endif