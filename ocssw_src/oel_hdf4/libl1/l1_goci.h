#ifndef L1_GOCI_H
#define L1_GOCI_H

int openl1_goci(filehandle *l1file);
int readl1_goci(filehandle *l1file, int recnum, l1str *l1rec, int lonlat);
int closel1_goci(filehandle *l1file);

#endif
