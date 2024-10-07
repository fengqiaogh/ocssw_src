#ifndef L1_L5TM_H
#define L1_L5TM_H

int openl1_l5tm( filehandle *l1file );
int readl1_l5tm( filehandle *l1file, int recnum, l1str *l1rec, int lonlat);
int closel1_l5tm( filehandle *l1file );

#endif
