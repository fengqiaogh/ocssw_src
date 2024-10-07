#ifndef TLM_PROTO_H_
#define TLM_PROTO_H_

extern char *get_tlm_name
PROTO((int Word_Offset, int Bit_Offset, int Bit_Length));

extern char *get_ana_tlm_name
PROTO((int index));

extern char *get_dis_tlm_name
PROTO((int index));

extern unsigned int get_ana_tlm
PROTO((byte *ptr, int index));

extern unsigned int get_dis_tlm
PROTO((byte *ptr, int index));

extern unsigned int get_tlm_mask
PROTO((int Bit_Length, int Bit_Offset));

extern unsigned int get_tlm
PROTO((byte *ptr, int Word_Offset, int Bit_Offset, int Bit_Length));

extern int set_tlm
PROTO((int val, byte *ptr, int Word_Offset, int Bit_Offset, int Bit_Length));

extern float get_tlm_eng_value
PROTO((byte *ptr, int Word_Offset, int Bit_Offset, int Bit_Length, int ecvtype, float slope, float intercept));

extern int get_all_ana_tlm
PROTO((byte *ptr, int *rawtlm, float *ecvtlm));

extern int get_all_dis_tlm
PROTO((byte *ptr, int *rawtlm, byte *ecvtlm));

extern int ANA_TLM_INDEX
PROTO((int woff, int boff, int blen));

extern int DIS_TLM_INDEX
PROTO((int woff, int boff, int blen));


#endif /* TLM_PROTO_H_ */
