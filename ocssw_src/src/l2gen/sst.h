#ifndef SST_H
#define SST_H
#include "l12_proto.h"
#ifdef __cplusplus
extern "C" {
#endif
void call_sst(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type, float **sst_out);
void call_sses_bias(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type,
                    float **sses_bias_out);
void call_sses_bias_mean(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type,
                         float **sses_bias_mean_out);
void call_sses_std(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type,
                   float **sses_std_out);
void call_sst_flags(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type, int16 **flags_out);
void call_qual_flags(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type, int8 **qual_out);
void call_sses_counts(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type,
                      int16_t **sses_counts);
void call_dust_correction(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type,
                          float **dust_correction);
#ifdef __cplusplus
}
#endif
#endif