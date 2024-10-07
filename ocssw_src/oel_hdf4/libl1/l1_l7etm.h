#ifndef L1_L7ETM_H
#define L1_L7ETM_H

int openl1_l7etm( filehandle *l1file );
int readl1_l7etm( filehandle *l1file, int recnum, l1str *l1rec, int lonlat);
int closel1_l7etm( filehandle *l1file );

#endif
