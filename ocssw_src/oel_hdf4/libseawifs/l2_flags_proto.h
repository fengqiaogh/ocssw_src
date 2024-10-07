#ifndef L2_FLAGS_PROTO_H_
#define L2_FLAGS_PROTO_H_

extern int maskbits
PROTO((char *masknames, char *flagnames, char *delimeter,
        unsigned short *bits));

extern int retrieve_flagnames
PROTO((int32 fid, char *sdsname, char *delimeter, char **flagnames));


#endif /* L2_FLAGS_PROTO_H_ */
