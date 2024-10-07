/*
 * time_utl.h  just the time conversion routine prototype defs
 */
#ifndef TIME_UTL_DEF
#define TIME_UTL_DEF
int jd_to_ydoy(int jd, int *y, int *doy);
int jd_to_yydoy(int jd, int *yy, int *doy);
int yymd_to_jd(int yy, int m, int d);
int yydoy_to_md70(int yy, int doy, int *m, int *d);
int yydoy_to_jd(int yy, int doy);
int ydoy_to_jd(int y, int doy);
int yydoy_to_md(int yy, int doy, int *m, int *d);
int leap(int);
#endif
/*
 * end
 */
