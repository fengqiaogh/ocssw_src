#ifndef ENG_QUAL_PROTO_H_
#define ENG_QUAL_PROTO_H_

extern void put_eng_qual_by_index
PROTO((eng_qualType *eng_qual, int bflag, int nindex));

extern void set_eng_qual_by_index
PROTO((eng_qualType *eng_qual, int nindex));

extern void reset_eng_qual_by_index
PROTO((eng_qualType *eng_qual, int nindex));

extern int put_eng_qual
PROTO((eng_qualType *eng_qual, int bflag, int woff, int boff, int blen));

extern void set_eng_qual
PROTO((eng_qualType *eng_qual, int woff, int boff, int blen));

extern void reset_eng_qual
PROTO((eng_qualType *eng_qual, int woff, int boff, int blen));

extern int fill_eng_qual
PROTO((eng_qualType *eng_qual, float *inst_ana));

extern int get_eng_qual_bit_by_index
PROTO((eng_qualType eng_qual, int nindex));


#endif /* ENG_QUAL_PROTO_H_ */
