#ifndef STATPROTO_H
#define STATPROTO_H
#include "usrmac.h"
extern int32 read_cntldata
PROTO((char *cntl_file, cntl_str *grschk, cntl_str *statchk,
        flag_str *flgchk));

extern int32 l2file
PROTO((int32 sdfid, int32 *nsamp, int32 *nscans, char *dtype));

extern int32 chk_grsstat
PROTO((int32 sdfid, int32 nscans, int32 nsamp, cntl_str *grschk,
        cntl_str *statchk, char *extra_flags));

extern int16 set_mask
PROTO((int32 sdfid, char *extra_flags, uint16 *mask));

extern int32 chk_flg
PROTO((int32 sdfid, flag_str *flgchk));

extern void stat_exit
PROTO((int status));

extern int32 rdattr
PROTO((int32 sdfid, char *attr_name, void *buf));

extern int32 rdslice
PROTO((int32 sdfid, char *name, int32 *start, int32 *edge, void *buf));

extern int32 getattrsz
PROTO((int32 id, char *attr_name, int32 *nt, int32 *count));

#endif /* STATPROTO_H */
