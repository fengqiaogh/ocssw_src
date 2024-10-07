/*
 *  W. Robinson, SAIC, 10 Dec 2004  new for CZCS
 */
#ifndef L1_CZCSW_H
#define L1_CZCSW_H

int openl1_czcs(filehandle *file);
int readl1_czcs(filehandle *file, int32_t recnum, l1str *l1rec);
int closel1_czcs(filehandle *file);

#endif
