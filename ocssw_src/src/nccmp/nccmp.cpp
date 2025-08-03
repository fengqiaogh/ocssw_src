/*
   Copyright (C) 2004-2007,2009,2010,2012 Remik Ziemlinski @ noaa gov

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; see the file COPYING.
   If not, write to the Free Software Foundation,
   59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
 **
   Converted to C++ and generalized compare using templates
   January 2015 by R. Healy (richard.healy@nasa.gov)

   9/3/2015 - Added call to freestringlist in makecmpvarlist to
              make number of variables to compare calculated correctly.
              Previously, if there was more than one group in the file
              then the value to compare was always from the first group.
              Also added uint16_t(ushort) and ubyte to the types of netcdf variables

   1/6/2015 - In the case where type is not supported limited INFO message

   9/29/2016 - Added ability to compare compound types - rjh

 */

#include "nccmp.hpp"
#include "ncinfo.h"
#include "strlist.h"
#include <stdint.h>
#include <float.h>
#include "strlist.c"
#include "ncinfo.c"
#include "opt.c"

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#define NC_MAX_TYPES 64

using namespace std;

varstruct vars1[(int) NC_MAX_VARS], vars2[(int) NC_MAX_VARS];
dimstruct dims1[(int) NC_MAX_DIMS], dims2[(int) NC_MAX_DIMS];
static int notsupported[NC_MAX_TYPES];

size_t nrec1, nrec2;
int nvars1, nvars2, ndims1, ndims2, recid1, recid2;

vector<string> groupPath;

#define NCFORMATSTR(f)                                          \
    (f == NC_FORMAT_CLASSIC ? "NC_FORMAT_CLASSIC" :               \
      (f == NC_FORMAT_64BIT ? "NC_FORMAT_64BIT" :                 \
      (f == NC_FORMAT_NETCDF4 ? "NC_FORMAT_NETCDF4" :             \
     "NC_FORMAT_NETCDF4_CLASSIC")))                                 \

const char* getGroupPath() {
    static string path;

    if(groupPath.empty())
        return "";

    path = "/";
    for(vector<string>::iterator it = groupPath.begin(); it != groupPath.end(); ++it) {
        path += *it;
        path += "/";
    }
    return path.c_str();
}

template <typename T> int cmp_(T* in1, T* in2) {
    T const *p1, *p2;
    for (p1 = in1, p2 = in2; *p1 == *p2; ++p1, ++p2)
        continue;
    return p1 - in1;
}

template <typename T> int cmp_missing(T* in1, T* in2, T m1, T m2) {
    T const *p1, *p2;
    for (p1 = in1, p2 = in2; (*p1 == *p2) || ((*p1 == m1) && (*p2 == m2)); ++p1, ++p2)
        continue;
    return p1 - in1;
}

template <typename T> int cmp_nanequal(T* in1, T* in2) {
    T const *p1, *p2;
    for (p1 = in1, p2 = in2; (*p1 == *p2) || (std::isnan(*p1) && std::isnan(*p2)); ++p1, ++p2)
        continue;
    return p1 - in1;
}

template <typename T> int cmp_missing_nanequal(T* in1, T* in2, T m1, T m2) {
    T const *p1, *p2;
    for (p1 = in1, p2 = in2; (*p1 == *p2) || ((*p1 == m1) && (*p2 == m2)) || (std::isnan(*p1) && std::isnan(*p2)); ++p1, ++p2)
        continue;
    return p1 - in1;
}

template <typename T> void ToHex(T v, char* out) { \
  unsigned char *p = (unsigned char*) &v;
    int i;
    char tmp[3];

    strcpy(out, "0x");

    for (i = 0; i < (int) sizeof (T); ++i) {
        sprintf(tmp, "%02X", p[i]);
        strcat(out, tmp);
    }
}

template <typename T> void ToString(T v, char* out, char* formatprec) { \
  sprintf(out, formatprec, v);
}

template <typename T> int cmp_var(int ncid1, int ncid2, nccmpopts* opts, int rec, size_t *odomax,
        off_t nitems, size_t *count, size_t *start, varstruct *v1, varstruct *v2, T M1, T M2) {

    void getidxstr(varstruct* var, size_t* start, int curidx, char* out);
    void getidxstr_fortran(varstruct* var, size_t* start, int curidx, char* out);
    int status, cmplen, do_missing = 0, diffstatus = EXIT_SUCCESS;
    char idxstr[256];
    string message("DIFFER : VARIABLE : %s%s : POSITION : %s : VALUES : %s <> %s\n");
    char value1str[32], value2str[32];
    char printHex = strchr(opts->precision, 'X') || strchr(opts->precision, 'x');

    T *P1, *P2;
    T value1, value2;
    int diff;

    if (opts->missing && v1->hasmissing && v2->hasmissing) {
        do_missing = 1;
    }

    P1 = (T *) malloc(sizeof (T) * (nitems + 1));
    P2 = (T *) malloc(sizeof (T) * (nitems + 1));

    do {
        /* printf("start = %d %d %d %d, count = %d %d %d %d\n", (int)start[0], (int)start[1], (int)start[2], (int)start[3], (int)count[0], (int)count[1], (int)count[2], (int)count[3]); */ \
      status = nc_get_vara(ncid1, v1->varid, start, count, P1);
        status = nc_get_vara(ncid2, v2->varid, start, count, P2);
        handle_error(status);
        /* Sentinels. */
        P1[nitems] = 0;
        P2[nitems] = 1;

        cmplen = nitems;
        /* for(i=0; i<nitems; ++i) { \
            printf("nitems = %d, rec = %d, P1[%d] = %g, P2[%d] = %g\n", nitems, rec, i, P1[i], i, P2[i]); \
        } */
        if (do_missing)
            if (opts->nanequal)
                diff = cmp_missing_nanequal<T>(P1, P2, M1, M2);
            else
                diff = cmp_missing<T>(P1, P2, M1, M2);
        else
            if (opts->nanequal)
            diff = cmp_nanequal<T>(P1, P2);
        else
            diff = cmp_<T>(P1, P2);

        while (diff < cmplen && (opts->maxdiff == 0 || opts->diffcount < opts->maxdiff)) {
            //  printf("RJH: maxdiff=%ld cnt=%ld\n",opts->maxdiff,opts->diffcount);
            if (!opts->warn[NCCMP_W_ALL])
                diffstatus = EXIT_DIFFER;
            else
                diffstatus = EXIT_SUCCESS;

            if (opts->fortran)
                getidxstr_fortran(v1, start, nitems - cmplen + diff, idxstr);
            else
                getidxstr(v1, start, nitems - cmplen + diff, idxstr);

            value1 = P1[nitems - cmplen + diff];
            value2 = P2[nitems - cmplen + diff];
            //printf("1)Forcing...%d %d \n",P1[+nitems-cmplen+diff],P2[nitems-cmplen+diff]);
            if (printHex) {
                ToHex((T) value1, value1str);
                ToHex((T) value2, value2str);
            } else {
                ToString((T) value1, value1str, opts->tprintf);
                ToString((T) value2, value2str, opts->tprintf);
            }
            opts->diffcount++;
            //fprintf(stdout, message.c_str(), v1->name, idxstr, value1str, value2str);
            printf(message.c_str(), getGroupPath(), v1->name, idxstr, value1str, value2str);
            //printf("%d nitems=%d, cmplen=%d, diff=%d do_missing=%d\n", __LINE__, nitems, cmplen, diff, opts->missing);
            if (opts->force) {
                cmplen = cmplen - (diff + 1);
                diff = cmp_<T>((P1 + diff + 1), (P2 + diff + 1));
                //printf("Forcing...%d %d \n",P1[nitems-cmplen+diff],P2[nitems-cmplen+diff]);
            } else
                goto break_it;
        }
        /* Increment all but first (if record) and last dimensions. */
    } while (odometer(start, odomax, (int) (rec >= 0), (int) v1->ndims - 2));
break_it:
    free(P1);
    free(P2);
    return diffstatus;
}

template <typename T> int cmp_vartol(int ncid1, int ncid2, nccmpopts* opts, int rec, size_t *odomax,
        off_t nitems, size_t *count, size_t *start, varstruct *v1, varstruct *v2, T M1, T M2) {
    void getidxstr_fortran(varstruct* var, size_t* start, int curidx, char* out);
    void getidxstr(varstruct* var, size_t* start, int curidx, char* out);

    T *P1, *P2;

    int i, status, do_missing = 0, diffstatus = EXIT_SUCCESS;
    char idxstr[256];
    string message("DIFFER : VARIABLE : %s%s : POSITION : %s : VALUES : %s <> %s : PERCENT : %g\n");
    char value1str[32], value2str[32];
    char printHex = strchr(opts->precision, 'X') || strchr(opts->precision, 'x');
    double absdelta;

    T value1, value2;

    P1 = (T *) malloc(sizeof (T) * (nitems + 1));
    P2 = (T *) malloc(sizeof (T) * (nitems + 1));

    if (opts->missing && v1->hasmissing && v2->hasmissing) {
        do_missing = 1;
    }

    do {
        /* printf("start = %d %d %d %d, count = %d %d %d %d\n", (int)start[0], (int)start[1], (int)start[2], (int)start[3], (int)count[0], (int)count[1], (int)count[2], (int)count[3]); */ \
        if (opts->extent && opts->extentcount >= opts->extent) {
            break;
        }
      status = nc_get_vara(ncid1, v1->varid, start, count, P1);
        handle_error(status);
        status = nc_get_vara(ncid2, v2->varid, start, count, P2);
        handle_error(status);

        /* for(i=0; i<nitems; ++i) {
            printf("nitems = %d, rec = %d, P1[%d] = %g, P2[%d] = %g\n", nitems, rec, i, P1[i], i, P2[i]); \
        } */
        for (i = 0; i < nitems && 
                        (opts->maxdiff == 0 || opts->diffcount < opts->maxdiff) && 
                        (opts->extent == 0 || opts->extentcount < opts->extent); 
                        ++i) {
            if (do_missing) {
                if ((M1 == P1[i]) && (M2 == P2[i])) continue;
            }

            absdelta = fabs((double) (P1[i] - P2[i]));
            if (opts->abstolerance ? (absdelta > opts->tolerance) : (double) absdelta * 100. / (fabs((double) P1[i]) > fabs((double) P2[i]) ? fabs((double) P1[i]) : fabs((double) P2[i])) > opts->tolerance)  \
         {
                if (!opts->warn[NCCMP_W_ALL])
                    diffstatus = EXIT_DIFFER;
                else
                    diffstatus = EXIT_SUCCESS;

                if ((v1->ndims == 1) && (v1->hasrec)) {
                    if (opts->fortran)
                        getidxstr_fortran(v1, start, rec, idxstr);
                    else
                        getidxstr(v1, start, rec, idxstr);
                } else {
                    if (opts->fortran)
                        getidxstr_fortran(v1, start, i, idxstr);
                    else
                        getidxstr(v1, start, i, idxstr);
                }
                value1 = P1[i];
                value2 = P2[i];
                if (printHex) {
                    ToHex(value1, value1str);
                    ToHex(value2, value2str);
                } else {
                    ToString(value1, value1str, opts->tprintf);
                    ToString(value2, value2str, opts->tprintf);
                }
                fprintf(stdout, message.c_str(), getGroupPath(), v1->name, idxstr, value1str, value2str, (double) absdelta * 100. / (fabs((double) P1[i]) > fabs((double) P2[i]) ? fabs((double) P1[i]) : fabs((double) P2[i])));
                opts->diffcount++;
                if (!opts->force) {
                    goto break_;
                }
            }
            opts->extentcount++;
        }
    } while (odometer(start, odomax, (rec >= 0), v1->ndims - 2));
break_:
    free(P1);
    free(P2);

    return diffstatus;
}

template <typename T> int cmp_var_ut(int ncid1, int ncid2, nccmpopts* opts, int rec, size_t *odomax,
        off_t nitems, size_t *count, size_t *start, varstruct *v1, varstruct *v2, T M1, T M2) {

    void getidxstr(varstruct* var, size_t* start, int curidx, char* out);
    void getidxstr_fortran(varstruct* var, size_t* start, int curidx, char* out);
    int status, cmplen, do_missing = 0, diffstatus = EXIT_SUCCESS;
    char idxstr[256];
    string message("DIFFER : VARIABLE : %s%s : POSITION : %s : VALUES : %s <> %s\n");
    char value1str[32], value2str[32];
    char printHex = strchr(opts->precision, 'X') || strchr(opts->precision, 'x');

    T *P1, *P2;
    T value1, value2;
    int diff;

    if (opts->missing && v1->hasmissing && v2->hasmissing) {
        do_missing = 1;
    }

    P1 = (T *) malloc(sizeof (T) * (nitems + 1));
    P2 = (T *) malloc(sizeof (T) * (nitems + 1));

    do {
        /* printf("start = %d %d %d %d, count = %d %d %d %d\n", (int)start[0], (int)start[1], (int)start[2], (int)start[3], (int)count[0], (int)count[1], (int)count[2], (int)count[3]); */ \
      status = nc_get_vara(ncid1, v1->varid, start, count, P1);
        status = nc_get_vara(ncid2, v2->varid, start, count, P2);
        handle_error(status);
        /* Sentinels. */
        P1[nitems] = 0;
        P2[nitems] = 1;

        cmplen = nitems;
        /* for(i=0; i<nitems; ++i) { \
            printf("nitems = %d, rec = %d, P1[%d] = %g, P2[%d] = %g\n", nitems, rec, i, P1[i], i, P2[i]); \
        } */
        if (do_missing)
            if (opts->nanequal)
                diff = cmp_missing_nanequal<T>(P1, P2, M1, M2);
            else
                diff = cmp_missing<T>(P1, P2, M1, M2);
        else
            if (opts->nanequal)
            diff = cmp_nanequal<T>(P1, P2);
        else
            diff = cmp_<T>(P1, P2);

        while (diff < cmplen && (opts->maxdiff == 0 || opts->diffcount < opts->maxdiff)) {
            //  printf("RJH: maxdiff=%ld cnt=%ld\n",opts->maxdiff,opts->diffcount);
            if (!opts->warn[NCCMP_W_ALL])
                diffstatus = EXIT_DIFFER;
            else
                diffstatus = EXIT_SUCCESS;

            if (opts->fortran)
                getidxstr_fortran(v1, start, nitems - cmplen + diff, idxstr);
            else
                getidxstr(v1, start, nitems - cmplen + diff, idxstr);

            value1 = P1[nitems - cmplen + diff];
            value2 = P2[nitems - cmplen + diff];
            //printf("1)Forcing...%d %d \n",P1[+nitems-cmplen+diff],P2[nitems-cmplen+diff]);
            if (printHex) {
                ToHex((T) value1, value1str);
                ToHex((T) value2, value2str);
            } else {
                ToString((T) value1, value1str, opts->tprintf);
                ToString((T) value2, value2str, opts->tprintf);
            }
            opts->diffcount++;
            //fprintf(stdout, message.c_str(), v1->name, idxstr, value1str, value2str);
            printf(message.c_str(), getGroupPath(), v1->name, idxstr, value1str, value2str);
            //printf("%d nitems=%d, cmplen=%d, diff=%d do_missing=%d\n", __LINE__, nitems, cmplen, diff, opts->missing);
            if (opts->force) {
                cmplen = cmplen - (diff + 1);
                diff = cmp_<T>((P1 + diff + 1), (P2 + diff + 1));
                //printf("Forcing...%d %d \n",P1[nitems-cmplen+diff],P2[nitems-cmplen+diff]);
            } else
                goto break_it;
        }
        /* Increment all but first (if record) and last dimensions. */
    } while (odometer(start, odomax, (int) (rec >= 0), (int) v1->ndims - 2));
break_it:
    free(P1);
    free(P2);
    return diffstatus;
}

template <typename T> int cmp_vartol_ut(void *P1, void *P2, int offset1, int offset2, int size1, int size2, char* name, nccmpopts* opts, int rec, size_t *odomax,
        off_t nitems, size_t *count, size_t *start, varstruct *v1, varstruct *v2, T M1, T M2) {

    void getidxstr_fortran(varstruct* var, size_t* start, int curidx, char* out);
    void getidxstr(varstruct* var, size_t* start, int curidx, char* out);


    int i, do_missing = 0, diffstatus = EXIT_SUCCESS;
    char idxstr[256];
    string message("DIFFER : VARIABLE : %s%s(%s) : POSITION : %s : VALUES : %s <> %s : PERCENT : %g\n");
    char value1str[32], value2str[32];
    char printHex = strchr(opts->precision, 'X') || strchr(opts->precision, 'x');
    double absdelta;

    T value1, value2;

    if (opts->missing && v1->hasmissing && v2->hasmissing) {
        do_missing = 1;
    }

    do {

        /* for(i=0; i<nitems; ++i) {
            printf("nitems = %d, rec = %d, P1[%d] = %g, P2[%d] = %g\n", nitems, rec, i, P1[i], i, P2[i]); \
        } */
        for (i = 0; i < nitems && (opts->maxdiff == 0 || opts->diffcount < opts->maxdiff); ++i) {
            value1 = *((T *) ((unsigned char*) P1 + offset1) + i);
            value2 = *((T *) ((unsigned char*) P2 + offset2) + i);
            if (do_missing) {
                if ((M1 == value1) && (M2 == value2)) continue;
            }

            absdelta = fabs((double) (value1 - value2));
            if (opts->abstolerance ? (absdelta > opts->tolerance) : (double) absdelta * 100. / (fabs((double) value1) > fabs((double) value2) ? fabs((double) value1) : fabs((double) value2)) > opts->tolerance)  \
         {
                if (!opts->warn[NCCMP_W_ALL])
                    diffstatus = EXIT_DIFFER;
                else
                    diffstatus = EXIT_SUCCESS;

                if ((v1->ndims == 1) && (v1->hasrec)) {
                    if (opts->fortran)
                        getidxstr_fortran(v1, start, rec, idxstr);
                    else
                        getidxstr(v1, start, rec, idxstr);
                } else {
                    if (opts->fortran)
                        getidxstr_fortran(v1, start, i, idxstr);
                    else
                        getidxstr(v1, start, i, idxstr);
                }
                if (printHex) {
                    ToHex(value1, value1str);
                    ToHex(value2, value2str);
                } else {
                    ToString(value1, value1str, opts->tprintf);
                    ToString(value2, value2str, opts->tprintf);
                }
                fprintf(stdout, message.c_str(), getGroupPath(), v1->name, name, idxstr, value1str, value2str, (double) absdelta * 100. / (fabs((double) value1) > fabs((double) value2) ? fabs((double) value1) : fabs((double) value2)));
                opts->diffcount++;
                if (!opts->force) {
                    goto break_;
                }
            }
        }
    } while (odometer(start, odomax, (rec >= 0), v1->ndims - 2));
break_:

    return diffstatus;
}

/* *********************************************************** */

/* Returns formatted string of dimension indices.
  Intended to print out locations of value differences. 
  'out' should be pre-allocated large enough for output. */
void getidxstr(varstruct* var, size_t* start, int curidx, char* out) {
    int i;
    char tmp[8];
    memset(out, '\0', 32);
    for (i = 0; i < var->ndims - 1; ++i) {
        sprintf(tmp, "%d ", (int) start[i]);
        strcat(out, tmp);
    }
    sprintf(tmp, "%d", curidx);
    strcat(out, tmp);
}
/* *********************************************************** */

/* Same as above but with fortran style indices, which means they're
  1-based (i.e. first index is 1, not 0) and fast-varying dimension is
  printed first (reverse order compared with C-style indices) */
void getidxstr_fortran(varstruct* var, size_t* start, int curidx, char* out) {
    int i;
    char tmp[8];
    memset(out, '\0', 32);
    sprintf(tmp, "%d", curidx + 1);
    strcat(out, tmp);

    for (i = var->ndims - 2; i >= 0; --i) {
        sprintf(tmp, " %d", (int) start[i] + 1);
        strcat(out, tmp);
    }
}

/* *********************************************************** */
int
excludevars(int ncid1, int ncid2, char** finallist,
        int nfinal, char** excludelist, int nexclude) {
    int nvars, nmaxvars = NC_MAX_VARS;
    char** vars = NULL;
    int status = EXIT_SUCCESS;

    vars = NULL;

    if (nexclude == 0)
        return EXIT_SUCCESS;

    /* printf("%d: creating temporary var list array.\n", __LINE__); */
    /* get simple difference */
    if (newstringlist(&vars, &nvars, nmaxvars))
        status = EXIT_FATAL;

    /* printf("%d: getting all variables names from both input files.\n", __LINE__); */
    if (allvarnames(vars, nvars, ncid1, ncid2))
        status = EXIT_FATAL;

    /*printf("vars=");
     printstrlist(vars, nvars, stdout);
     */

    if (strlistsd(vars, excludelist, finallist,
            nvars, nexclude, nfinal))
        status = EXIT_FATAL;

    /*printf("excludelist=");
     printstrlist(excludelist, nexclude, stdout);
     */

    /*printf("finallist=");
     printstrlist(finallist, nfinal, stdout);
     */

    freestringlist(&vars, nvars);

    return status;
}

/* *********************************************************** */
void handle_error(int status) {
    if (status != NC_NOERR) {
        fprintf(stderr, "%s\n", nc_strerror(status));
        exit(-1);
    }
}

/* *********************************************************** */

/*
  Mimics incrementing odometer. 
  Returns 0 if rolled-over. 
  
  @param odo: the counters that are updated.
  @param limits: the maximum values allowed per counter.
  @param first: index of first counter to update.
  @param last: index of last counter to update.
 */
int odometer(size_t* odo, size_t* limits, int first, int last) {
    int i = last;
    while (i >= first) {
        odo[i] += 1;
        if (odo[i] > limits[i])
            odo[i] = 0;
        else
            break;

        --i;
    }

#ifdef __DEBUG__
    printf("DEBUG : %d : odo = ", __LINE__);
    for (i = first; i <= last; ++i) {
        printf("%d ", odo[i]);
    }
    printf("\n");
#endif

    /* Test for roll over. */
    for (i = first; i <= last; ++i) {
        if (odo[i] != 0)
            return 1;
    }

    /* If we get here then rolled over. */
    return 0;
}
/* *********************************************************** */

/* Pretty prints attribute values into a string. 
  Assumes 'str' is pre-allocated large enough to hold output string.
 */
void prettyprintatt(int ncid, char* varname, int varid, char* name, char* str) {
    int status, i;
    nc_type type;
    size_t len;
    char* pc;
    char charpr[3];
    int8_t* puc;
    uint8_t* upuc;
    short* ps;
    uint16_t* ups;
    int* pi;
    uint* upi;
    long* pl;
    uint64_t* upl;
    float* pf;
    double* pd;
    char tmpstr[32];


    strcpy(charpr, "%c");

    status = nc_inq_att(ncid, varid, name, &type, &len);
    if (status != NC_NOERR) {
        if (varid == NC_GLOBAL)
            fprintf(stderr, "ERROR : QUERYING GLOBAL ATTRIBUTE \"%s\"\n", name);
        else
            fprintf(stderr, "ERROR : QUERYING ATTRIBUTE \"%s\" FOR VARIABLE \"%s%s\"\n", name, getGroupPath(), varname);
        return;
    }

    str[0] = '\0';
    if (len < 1) {
        return;
    }

    switch (type) {
    case NC_BYTE:
        puc = XMALLOC(int8_t, len);

        status = nc_get_att(ncid, varid, name, puc);
        if (status != NC_NOERR) {
            XFREE(puc);
            return;
        }

        for (i = 0; i < (int) len; ++i) {
            ToString<unsigned char> (puc[i], str + 2 * i, charpr);
            str[2 * i + 1] = ',';
        }

        XFREE(puc);
        str[2 * (int) len - 1] = '\0';
        break;

    case NC_UBYTE:
        upuc = XMALLOC(uint8_t, len);

        status = nc_get_att_ubyte(ncid, varid, name, upuc);
        if (status != NC_NOERR) {
            XFREE(upuc);
            return;
        }

        for (i = 0; i < (int) len; ++i) {
            ToString<unsigned char> (upuc[i], str + 2 * i, charpr);
            str[2 * i + 1] = ',';
        }

        XFREE(upuc);
        str[2 * (int) len - 1] = '\0';
        break;


    case NC_CHAR:
        pc = XMALLOC(char, len);
        status = nc_get_att_text(ncid, varid, name, pc);
        if (status != NC_NOERR) {
            XFREE(pc);
            return;
        }

        for (i = 0; i < (int) len; ++i)
            ToString<char>(pc[i], str + i, charpr);

        XFREE(pc);
        str[(int) len] = '\0';
        break;

    case NC_SHORT:
        ps = XMALLOC(short, len);
        status = nc_get_att_short(ncid, varid, name, ps);
        if (status != NC_NOERR) {
            XFREE(ps);
            return;
        }

        for (i = 0; i < (int) len; ++i) {
            sprintf(tmpstr, "%d,", ps[i]);
            strcat(str, tmpstr);
        }

        XFREE(ps);
        str[strlen(str) - 1] = '\0'; // Remove last comma.
        break;

    case NC_USHORT:
        ups = XMALLOC(uint16_t, len);
        status = nc_get_att_ushort(ncid, varid, name, ups);
        if (status != NC_NOERR) {
            XFREE(ups);
            return;
        }

        for (i = 0; i < (int) len; ++i) {
            sprintf(tmpstr, "%d,", ups[i]);
            strcat(str, tmpstr);
        }

        XFREE(ups);
        str[strlen(str) - 1] = '\0'; // Remove last comma.
        break;

    case NC_INT:
        pi = XMALLOC(int, len);
        status = nc_get_att_int(ncid, varid, name, pi);
        if (status != NC_NOERR) {
            XFREE(pi);
            return;
        }

        for (i = 0; i < (int) len; ++i) {
            sprintf(tmpstr, "%d,", pi[i]);
            strcat(str, tmpstr);
        }

        XFREE(pi);
        str[strlen(str) - 1] = '\0'; // Remove last comma.
        break;

    case NC_UINT:
        upi = XMALLOC(uint, len);
        status = nc_get_att_uint(ncid, varid, name, upi);
        if (status != NC_NOERR) {
            XFREE(upi);
            return;
        }

        for (i = 0; i < (int) len; ++i) {
            sprintf(tmpstr, "%d,", upi[i]);
            strcat(str, tmpstr);
        }

        XFREE(upi);
        str[strlen(str) - 1] = '\0'; // Remove last comma.
        break;

    case NC_INT64:
        pl = XMALLOC(long, len);
        status = nc_get_att_long(ncid, varid, name, pl);
        if (status != NC_NOERR) {
            XFREE(pl);
            return;
        }

        for (i = 0; i < (int) len; ++i) {
            sprintf(tmpstr, "%ld,", pl[i]);
            strcat(str, tmpstr);
        }

        XFREE(pl);
        str[strlen(str) - 1] = '\0'; // Remove last comma.
        break;

    case NC_UINT64:
        upl = XMALLOC(uint64_t, len);
        status = nc_get_att_long(ncid, varid, name, (long*)upl);
        if (status != NC_NOERR) {
            XFREE(upl);
            return;
        }

        for (i = 0; i < (int) len; ++i) {
	  sprintf(tmpstr, "%llu,", (long long unsigned int)upl[i]);
            strcat(str, tmpstr);
        }

        XFREE(upl);
        str[strlen(str) - 1] = '\0'; // Remove last comma.
        break;

    case NC_FLOAT:
        pf = XMALLOC(float, len);
        status = nc_get_att_float(ncid, varid, name, pf);
        if (status != NC_NOERR) {
            XFREE(pf);
            return;
        }

        for (i = 0; i < (int) len; ++i) {
            sprintf(tmpstr, "%.9g,", pf[i]);
            strcat(str, tmpstr);
        }

        XFREE(pf);
        str[strlen(str) - 1] = '\0'; // Remove last comma.
        break;

    case NC_DOUBLE:
        pd = XMALLOC(double, len);
        status = nc_get_att_double(ncid, varid, name, pd);
        if (status != NC_NOERR) {
            XFREE(pd);
            return;
        }

        for (i = 0; i < (int) len; ++i) {
            sprintf(tmpstr, "%.17g,", pd[i]);
            strcat(str, tmpstr);
        }

        XFREE(pd);
        str[strlen(str) - 1] = '\0'; // Remove last comma.
        break;
    }
}

/* *********************************************************** */
int cmpatt(int ncid1, int ncid2, int varid1, int varid2,
        char* name, char* varname, nccmpopts* opts) {
    int ncstatus, status;
    nc_type type1, type2;
    size_t lenp1, lenp2;
    char typestr1[4096];
    char typestr2[4096];

    status = EXIT_SUCCESS;
    ncstatus = nc_inq_att(ncid1, varid1, name, &type1, &lenp1);
    if (ncstatus != NC_NOERR) {
        fprintf(stdout, "DIFFER : VARIABLE \"%s%s\" IS MISSING ATTRIBUTE WITH NAME \"%s\" IN FILE \"%s\"\n", getGroupPath(), varname, name, opts->file1);

        if (!opts->warn[NCCMP_W_ALL])
            status = EXIT_DIFFER;
        return status;
    }

    ncstatus = nc_inq_att(ncid2, varid2, name, &type2, &lenp2);
    if (ncstatus != NC_NOERR) {
        fprintf(stdout, "DIFFER : VARIABLE \"%s%s\" IS MISSING ATTRIBUTE WITH NAME \"%s\" IN FILE \"%s\"\n", getGroupPath(), varname, name, opts->file2);

        if (!opts->warn[NCCMP_W_ALL])
            status = EXIT_DIFFER;
        return status;
    }

    if (type1 != type2) {
        type2string(type1, typestr1);
        type2string(type2, typestr2);
        fprintf(stdout, "DIFFER : TYPES : ATTRIBUTE : %s : VARIABLE : %s%s : %s <> %s\n", name, getGroupPath(), varname, typestr1, typestr2);

        if (!opts->warn[NCCMP_W_ALL])
            status = EXIT_DIFFER;
        if (!opts->force) 
            return status;
    }

    if (lenp1 != lenp2) {
        prettyprintatt(ncid1, varname, varid1, name, typestr1);
        prettyprintatt(ncid2, varname, varid2, name, typestr2);

        fprintf(stdout, "DIFFER : LENGTHS : ATTRIBUTE : %s : VARIABLE : %s%s : %lu <> %lu : VALUES : ", name, getGroupPath(), varname, (unsigned long) lenp1, (unsigned long) lenp2);

        switch (type1) {
        case NC_CHAR:
            /* Quote strings. */
            fprintf(stdout, "\"%s\" : \"%s\"\n", typestr1, typestr2);
            if (strcmp(typestr1, typestr2) == 0) {
                /* Same text, but invisible trailing nulls because lengths differ. */
                if (opts->warn[NCCMP_W_EOS] || opts->warn[NCCMP_W_ALL]) {
                    /* Pass */
                } else {
                    status = EXIT_DIFFER;
                    if (!opts->force) 
                        return status;
                }
            }
            break;
        default:
            /* Unquoted. */
            fprintf(stdout, "%s : %s\n", typestr1, typestr2);
            if (!opts->warn[NCCMP_W_ALL]) {
                status = EXIT_DIFFER;
                if (!opts->force) 
                    return status;
            }
            break;
        }
    } else if (cmpattval(ncid1, ncid2, varid1, varid2, name, lenp1, type1) != NC_NOERR) {
        prettyprintatt(ncid1, varname, varid1, name, typestr1);
        prettyprintatt(ncid2, varname, varid2, name, typestr2);
        fprintf(stdout, "DIFFER : VARIABLE : %s%s : ATTRIBUTE : %s : VALUES : ", getGroupPath(), varname, name);

        switch (type1) {
        case NC_CHAR:
            /* Quote strings. */
            fprintf(stdout, "\"%s\" <> \"%s\"\n", typestr1, typestr2);
            break;
        default:
            /* Unquoted. */
            fprintf(stdout, "%s <> %s\n", typestr1, typestr2);
            break;
        }

        if (!opts->warn[NCCMP_W_ALL]) {
            status = EXIT_DIFFER;
            if (!opts->force) 
                return status;
        }
    }

    return EXIT_SUCCESS;
}


/* *********************************************************** */

/* Assumes that both attributes have same length.
 */
int cmpattval(int nc1, int nc2, int varid1, int varid2, char* name, int len, nc_type type) {
    char* c1;
    int8_t* uc1;
    uint8_t* uuc1;
    short* s1;
    uint16_t* us1;
    int* i1;
    uint* ui1;
    long* l1;
    uint64_t* ul1;
    float* f1;
    double* d1;
    char* c2;
    int8_t* uc2;
    uint8_t* uuc2;
    short* s2;
    uint16_t* us2;
    int* i2;
    uint* ui2;
    long* l2;
    uint64_t* ul2;
    float* f2;
    double* d2;
    int status;
    int i;

    if (name == NULL) return NC_EINVAL;

    switch (type) {
    case NC_BYTE:
        uc1 = XMALLOC(int8_t, len);
        uc2 = XMALLOC(int8_t, len);
        status = nc_get_att(nc1, varid1, name, uc1);
        if (status != NC_NOERR) {
            XFREE(uc1);
            XFREE(uc2);
            return NC_EINVAL;
        }
        status = nc_get_att(nc2, varid2, name, uc2);
        if (status != NC_NOERR) {
            XFREE(uc1);
            XFREE(uc2);
            return NC_EINVAL;
        }
        for (i = 0; i < len; ++i) {
            if (uc1[i] != uc2[i]) {
                XFREE(uc1);
                XFREE(uc2);
                return EXIT_DIFFER;
            }
        }
        XFREE(uc1);
        XFREE(uc2);
        break;
    case NC_UBYTE:
        uuc1 = XMALLOC(uint8_t, len);
        uuc2 = XMALLOC(uint8_t, len);
        status = nc_get_att_ubyte(nc1, varid1, name, uuc1);
        if (status != NC_NOERR) {
            XFREE(uuc1);
            XFREE(uuc2);
            return NC_EINVAL;
        }
        status = nc_get_att_ubyte(nc2, varid2, name, uuc2);
        if (status != NC_NOERR) {
            XFREE(uuc1);
            XFREE(uuc2);
            return NC_EINVAL;
        }
        for (i = 0; i < len; ++i) {
            if (uuc1[i] != uuc2[i]) {
                XFREE(uuc1);
                XFREE(uuc2);
                return EXIT_DIFFER;
            }
        }
        XFREE(uuc1);
        XFREE(uuc2);
        break;
    case NC_CHAR:
        c1 = XMALLOC(char, len);
        c2 = XMALLOC(char, len);
        status = nc_get_att_text(nc1, varid1, name, c1);
        if (status != NC_NOERR) {
            XFREE(c1);
            XFREE(c2);
            return NC_EINVAL;
        }
        status = nc_get_att_text(nc2, varid2, name, c2);
        if (status != NC_NOERR) {
            XFREE(c1);
            XFREE(c2);
            return NC_EINVAL;
        }
        if (strncmp(c1, c2, len) != 0) {
            XFREE(c1);
            XFREE(c2);
            return EXIT_DIFFER;
        }
        XFREE(c1);
        XFREE(c2);
        break;
    case NC_SHORT:
        s1 = XMALLOC(short, len);
        s2 = XMALLOC(short, len);
        status = nc_get_att_short(nc1, varid1, name, s1);
        if (status != NC_NOERR) {
            XFREE(s1);
            XFREE(s2);
            return NC_EINVAL;
        }
        status = nc_get_att_short(nc2, varid2, name, s2);
        if (status != NC_NOERR) {
            XFREE(s1);
            XFREE(s2);
            return NC_EINVAL;
        }
        for (i = 0; i < len; ++i) {
            if (s1[i] != s2[i]) {
                XFREE(s1);
                XFREE(s2);
                return EXIT_DIFFER;
            }
        }
        XFREE(s1);
        XFREE(s2);
        break;
    case NC_USHORT:
        us1 = XMALLOC(uint16_t, len);
        us2 = XMALLOC(uint16_t, len);
        status = nc_get_att_ushort(nc1, varid1, name, us1);
        if (status != NC_NOERR) {
            XFREE(us1);
            XFREE(us2);
            return NC_EINVAL;
        }
        status = nc_get_att_ushort(nc2, varid2, name, us2);
        if (status != NC_NOERR) {
            XFREE(us1);
            XFREE(us2);
            return NC_EINVAL;
        }
        for (i = 0; i < len; ++i) {
            if (us1[i] != us2[i]) {
                XFREE(us1);
                XFREE(us2);
                return EXIT_DIFFER;
            }
        }
        XFREE(us1);
        XFREE(us2);
        break;
    case NC_INT:
        i1 = XMALLOC(int, len);
        i2 = XMALLOC(int, len);
        status = nc_get_att_int(nc1, varid1, name, i1);
        if (status != NC_NOERR) {
            XFREE(i1);
            XFREE(i2);
            return NC_EINVAL;
        }
        status = nc_get_att_int(nc2, varid2, name, i2);
        if (status != NC_NOERR) {
            XFREE(i1);
            XFREE(i2);
            return NC_EINVAL;
        }
        for (i = 0; i < len; ++i) {
            if (i1[i] != i2[i]) {
                XFREE(i1);
                XFREE(i2);
                return EXIT_DIFFER;
            }
        }
        XFREE(i1);
        XFREE(i2);
        break;
    case NC_UINT:
        ui1 = XMALLOC(uint, len);
        ui2 = XMALLOC(uint, len);
        status = nc_get_att_uint(nc1, varid1, name, ui1);
        if (status != NC_NOERR) {
            XFREE(ui1);
            XFREE(ui2);
            return NC_EINVAL;
        }
        status = nc_get_att_uint(nc2, varid2, name, ui2);
        if (status != NC_NOERR) {
            XFREE(ui1);
            XFREE(ui2);
            return NC_EINVAL;
        }
        for (i = 0; i < len; ++i) {
            if (ui1[i] != ui2[i]) {
                XFREE(ui1);
                XFREE(ui2);
                return EXIT_DIFFER;
            }
        }
        XFREE(ui1);
        XFREE(ui2);
        break;
    case NC_INT64:
        l1 = XMALLOC(long, len);
        l2 = XMALLOC(long, len);
        status = nc_get_att_long(nc1, varid1, name, l1);
        if (status != NC_NOERR) {
            XFREE(l1);
            XFREE(l2);
            return NC_EINVAL;
        }
        status = nc_get_att_long(nc2, varid2, name, l2);
        if (status != NC_NOERR) {
            XFREE(l1);
            XFREE(l2);
            return NC_EINVAL;
        }
        for (i = 0; i < len; ++i) {
            if (l1[i] != l2[i]) {
                XFREE(l1);
                XFREE(l2);
                return EXIT_DIFFER;
            }
        }
        XFREE(l1);
        XFREE(l2);
        break;
    case NC_UINT64:
        ul1 = XMALLOC(uint64_t, len);
        ul2 = XMALLOC(uint64_t, len);
        status = nc_get_att_long(nc1, varid1, name, (long*)ul1);
        if (status != NC_NOERR) {
            XFREE(ul1);
            XFREE(ul2);
            return NC_EINVAL;
        }
        status = nc_get_att_long(nc2, varid2, name, (long*)ul2);
        if (status != NC_NOERR) {
            XFREE(ul1);
            XFREE(ul2);
            return NC_EINVAL;
        }
        for (i = 0; i < len; ++i) {
            if (ul1[i] != ul2[i]) {
                XFREE(ul1);
                XFREE(ul2);
                return EXIT_DIFFER;
            }
        }
        XFREE(ul1);
        XFREE(ul2);
        break;
    case NC_FLOAT:
        f1 = XMALLOC(float, len);
        f2 = XMALLOC(float, len);
        status = nc_get_att_float(nc1, varid1, name, f1);
        if (status != NC_NOERR) {
            XFREE(f1);
            XFREE(f2);
            return NC_EINVAL;
        }
        status = nc_get_att_float(nc2, varid2, name, f2);
        if (status != NC_NOERR) {
            XFREE(f1);
            XFREE(f2);
            return NC_EINVAL;
        }
        for (i = 0; i < len; ++i) {
            if (f1[i] != f2[i]) {
                XFREE(f1);
                XFREE(f2);
                return EXIT_DIFFER;
            }
        }
        XFREE(f1);
        XFREE(f2);
        break;
    case NC_DOUBLE:
        d1 = XMALLOC(double, len);
        d2 = XMALLOC(double, len);
        status = nc_get_att_double(nc1, varid1, name, d1);
        if (status != NC_NOERR) {
            XFREE(d1);
            XFREE(d2);
            return NC_EINVAL;
        }
        status = nc_get_att_double(nc2, varid2, name, d2);
        if (status != NC_NOERR) {
            XFREE(d1);
            XFREE(d2);
            return NC_EINVAL;
        }
        for (i = 0; i < len; ++i) {
            if (d1[i] != d2[i]) {
                XFREE(d1);
                XFREE(d2);
                return EXIT_DIFFER;
            }
        }
        XFREE(d1);
        XFREE(d2);
        break;
    }

    return EXIT_SUCCESS;
}

/* *********************************************************** */
void type2string(nc_type type, char* str) {
    switch (type) {
    case NC_BYTE:
        strcpy(str, "BYTE");
        break;
    case NC_UBYTE:
        strcpy(str, "UBYTE");
        break;
    case NC_CHAR:
        strcpy(str, "CHAR");
        break;
    case NC_SHORT:
        strcpy(str, "SHORT");
        break;
    case NC_USHORT:
        strcpy(str, "USHORT");
        break;
    case NC_INT:
        strcpy(str, "INT");
        break;
    case NC_UINT:
        strcpy(str, "UINT");
        break;
    case NC_INT64:
        strcpy(str, "LONG");
        break;
    case NC_UINT64:
        strcpy(str, "ULONG");
        break;
    case NC_FLOAT:
        strcpy(str, "FLOAT");
        break;
    case NC_DOUBLE:
        strcpy(str, "DOUBLE");
        break;
    default:
        strcpy(str, "");
        break;
    }
}

/* *********************************************************** */
int openfiles(nccmpopts* opts, int *ncid1, int *ncid2) {
    int status;

    status = nc_open(opts->file1, NC_NOWRITE, ncid1);
    handle_error(status);

    status = nc_open(opts->file2, NC_NOWRITE, ncid2);
    handle_error(status);

    return 0;
}
/* *********************************************************** */

/* Compares record names and lengths. */
int
nccmprecinfo(nccmpopts* opts, int ncid1, int ncid2) {
    char name1[256], name2[256];
    int status;

    status = EXIT_SUCCESS;

    if (opts->verbose)
        printf("INFO: Comparing record information.\n");

    status = nc_inq_unlimdim(ncid1, &recid1);
    handle_error(status);

    if (recid1 != -1) {
        status = nc_inq_dim(ncid1, recid1, name1, &nrec1);
        handle_error(status);
    } else {
        nrec1 = 0;
    }

    status = nc_inq_unlimdim(ncid2, &recid2);
    handle_error(status);

    if (recid2 != -1) {
        status = nc_inq_dim(ncid2, recid2, name2, &nrec2);
        handle_error(status);
    } else {
        nrec2 = 0;
    }

    if (instringlist(opts->excludelist, name1, opts->nexclude) ||
            instringlist(opts->excludelist, name2, opts->nexclude) ||
            !instringlist(opts->variablelist, name1, opts->nvariable) ||
            !instringlist(opts->variablelist, name2, opts->nvariable))
        return EXIT_SUCCESS;

    if (strcmp(name1, name2)) {
        fprintf(stdout, "DIFFER : NAMES OF RECORDS : %s <> %s\n", name1, name2);
        if (!opts->warn[NCCMP_W_ALL])
            status = EXIT_DIFFER;

        if (!opts->force) return status;
    }

    if (nrec1 != nrec2) {
        fprintf(stdout, "DIFFER : LENGTHS OF RECORDS : %s (%d) <> %s (%d)\n", name1, (int) nrec1, name2, (int) nrec2);
        if (!opts->warn[NCCMP_W_ALL])
            status = EXIT_DIFFER;

        if (!opts->force) return status;
    }

    return status;
}

void getgroupinfo(int ncid, vector<string> names, GROUP_NODE *groups) {
    int *gids;
    size_t i;
    size_t index;
    int numGroups;
    char name[NC_MAX_NAME];
    vector<string> groupNames;

    nc_inq_grps(ncid, &numGroups, NULL);
    gids = (int *) malloc(sizeof (int) * numGroups);
    nc_inq_grps(ncid, NULL, gids);

    for (i = 0; i < (size_t)numGroups; i++) {
        nc_inq_grpname(gids[i], name);
        groupNames.push_back(name);
    }

    for (index = 0; index < names.size(); index++) {
        for (i = 0; i < groupNames.size(); i++) {
            if(groupNames[i] == names[index]) {
                groups[index].groupID = gids[i];
                strcpy(groups[index].groupName, names[index].c_str());
                break;
            }
        }
    }
    free(gids);
}

/* *********************************************************** */

/* Get dim info for file. */
void getdiminfo(int ncid, dimstruct* dims, int* ndims) {
    int status, i;

    status = nc_inq_ndims(ncid, ndims);
    handle_error(status);

    // Query all dimids, which may not be from 0..N-1 in a HDF file.
    int dimids[(int) NC_MAX_DIMS];
    int include_parents = 1;
    status = nc_inq_dimids(ncid, ndims, dimids, include_parents);
    handle_error(status);

    for (i = 0; i < *ndims; ++i) {
        dims[i].dimid = dimids[i];
        status = nc_inq_dim(ncid, dimids[i], dims[i].name, &dims[i].len);
        handle_error(status);
    }
}
/* *********************************************************** */

/* Copy attribute type to same type as var, just in case of mismatch. */
void broadcast_missing(nc_type var_type, nc_type att_type, missing_struct *values) {
#define BROADCAST_MISSING(T) { \
    switch(att_type) { \
      case NC_CHAR: values->T = values->c; break; \
      case NC_BYTE: values->T = values->b; break; \
      case NC_UBYTE: values->T = values->ub; break; \
      case NC_SHORT: values->T = values->s; break; \
      case NC_USHORT: values->T = values->us; break; \
      case NC_INT: values->T = values->i; break; \
      case NC_UINT: values->T = values->ui; break; \
      case NC_INT64: values->T = values->l; break; \
      case NC_UINT64: values->T = values->ul; break; \
      case NC_FLOAT: values->T = values->f; break; \
      case NC_DOUBLE: values->T = values->d; break; \
    } \
  }

    switch (var_type) {
    case NC_CHAR: BROADCAST_MISSING(c);
        break;
    case NC_BYTE: BROADCAST_MISSING(b);
        break;
    case NC_SHORT: BROADCAST_MISSING(s);
        break;
    case NC_UBYTE: BROADCAST_MISSING(ub);
        break;
    case NC_USHORT: BROADCAST_MISSING(us);
        break;
    case NC_INT: BROADCAST_MISSING(i);
        break;
    case NC_UINT: BROADCAST_MISSING(ui);
        break;
    case NC_INT64: BROADCAST_MISSING(l);
        break;
    case NC_UINT64: BROADCAST_MISSING(ul);
        break;
    case NC_FLOAT: BROADCAST_MISSING(f);
        break;
    case NC_DOUBLE: BROADCAST_MISSING(d);
        break;
    }
}

/* *********************************************************** */
char get_missing(int ncid, varstruct * var, const char* attname) {
    nc_type att_type;
    int status;

    status = nc_inq_atttype(ncid, var->varid, attname, &att_type);
    if (status != NC_NOERR) return 0;

    var->missing.type = att_type;

    switch (att_type) {
    case NC_CHAR:
        status = nc_get_att_text(ncid, var->varid, attname, &var->missing.c);
        if (status != NC_NOERR) return 0;
        break;
    case NC_BYTE:
        status = nc_get_att(ncid, var->varid, attname, &var->missing.b);
        if (status != NC_NOERR) return 0;
        break;
    case NC_UBYTE:
        status = nc_get_att_ubyte(ncid, var->varid, attname, &var->missing.ub);
        if (status != NC_NOERR) return 0;
        break;
    case NC_SHORT:
        status = nc_get_att_short(ncid, var->varid, attname, &var->missing.s);
        if (status != NC_NOERR) return 0;
        break;
    case NC_USHORT:
        status = nc_get_att_ushort(ncid, var->varid, attname, &var->missing.us);
        if (status != NC_NOERR) return 0;
        break;
    case NC_INT:
        status = nc_get_att_int(ncid, var->varid, attname, &var->missing.i);
        if (status != NC_NOERR) return 0;
        break;
    case NC_UINT:
        status = nc_get_att_uint(ncid, var->varid, attname, &var->missing.ui);
        if (status != NC_NOERR) return 0;
        break;
    case NC_INT64:
        status = nc_get_att_longlong(ncid, var->varid, attname, (long long *)&var->missing.l);
        if (status != NC_NOERR) return 0;
        break;
    case NC_UINT64:
        status = nc_get_att_ulonglong(ncid, var->varid, attname, (unsigned long long *)&var->missing.ul);
        if (status != NC_NOERR) return 0;
        break;
    case NC_FLOAT:
        status = nc_get_att_float(ncid, var->varid, attname, &var->missing.f);
        if (status != NC_NOERR) return 0;
        break;
    case NC_DOUBLE:
        status = nc_get_att_double(ncid, var->varid, attname, &var->missing.d);
        if (status != NC_NOERR) return 0;
        break;
    default: return 0;
    }

    var->hasmissing = 1;
    broadcast_missing(var->type, att_type, &var->missing);

    return 1;
}
/* *****************************s****************************** */

/* Read all the file's metadata for variables. */
nccmp_user_type_t* getvarinfo(int ncid, varstruct* vars, int* nvars, int verbose, int *nuser_types) {
    int status, i, j, recid;
    nccmp_user_type_t* comp_types = NULL;

    status = nc_inq_nvars(ncid, nvars);
    handle_error(status);

    // see if it has any compound types
    comp_types = nccmp_load_group_usertype_array(ncid, nuser_types);
    
    status = nc_inq_unlimdim(ncid, &recid);
    handle_error(status);

    for (i = 0; i < *nvars; ++i) {
        vars[i].varid = i;
        status = nc_inq_var(ncid, i, vars[i].name, &vars[i].type,
                &vars[i].ndims, vars[i].dimids, &vars[i].natts);

        vars[i].user_type_idx = -1;

        if (*nuser_types > 0) {
            for (j = 0; j<*nuser_types; j++) {
                if (vars[i].type == comp_types[j].type_id)
                    vars[i].user_type_idx = j;
            }
        }
        handle_error(status);

        vars[i].len = 1;
        for (j = 0; j < vars[i].ndims; ++j) {
            status = nc_inq_dimlen(ncid, vars[i].dimids[j], &vars[i].dimlens[j]);
            handle_error(status);
#ifdef __DEBUG__
            printf("DEBUG : %d : %s dimid %d, len %d\n", __LINE__, vars[i].name, j, vars[i].dimlens[j]);
#endif
            vars[i].len *= vars[i].dimlens[j];
        }

        vars[i].hasrec = (vars[i].dimids[0] == recid);

        /* Get missing_value or _FillValue. */
        if (get_missing(ncid, &vars[i], "missing_value") == 0)
            get_missing(ncid, &vars[i], "_FillValue");

        if (verbose) {
            if (vars[i].hasmissing) {
                printf("INFO: \"%s\" missing value: ", vars[i].name);
                switch (vars[i].missing.type) {
                case NC_BYTE: printf("%hhd (byte)\n", vars[i].missing.b);
                    break;
                case NC_UBYTE: printf("%hhu (ubyte)\n", vars[i].missing.ub);
                    break;
                case NC_CHAR: printf("%c (char)\n", vars[i].missing.c);
                    break;
                case NC_SHORT: printf("%hd (short)\n", vars[i].missing.s);
                    break;
                case NC_USHORT: printf("%hu (ushort)\n", vars[i].missing.us);
                    break;
                case NC_INT: printf("%d (int)\n", vars[i].missing.i);
                    break;
                case NC_UINT: printf("%u (uint)\n", vars[i].missing.ui);
                    break;
                 case NC_INT64: printf("%lld (int64)\n", (long long)vars[i].missing.l);
                    break;
                case NC_UINT64: printf("%llu (uint64)\n", (unsigned long long)vars[i].missing.ul);
                    break;
               case NC_FLOAT: printf("%g (float)\n", vars[i].missing.f);
                    break;
                case NC_DOUBLE: printf("%g (double)\n", vars[i].missing.d);
                    break;
                }
            }
        }
    }
    return comp_types;
}
/* *********************************************************** */

/* free all the file's metadata for variables. */
void freevarinfo(int nuser_types, nccmp_user_type_t* comp_types) {
    int i;
    nccmp_user_type_t* ptr = comp_types;

    for(i=0; i<nuser_types; i++) {
        free(ptr->fields);
        ptr++;
    }

    free(comp_types);
}
/* *********************************************************** */

/* Returns index to varstruct in list, otherwise -1. */
int isinvarstructlist(char* name, varstruct* vars, int nvars) {
    int i;

    for (i = 0; i < nvars; ++i) {
        if (strcmp(name, vars[i].name) == 0)
            return i;
    }

    return -1;
}
/* *********************************************************** */

/* Get vars to use to do cmp based on input and exclusion lists. */
int makecmpvarlist(nccmpopts* opts, int ncid1, int ncid2) {
    int status, i;

    //  opts->ncmpvarlist = (int)NC_MAX_VARS;
    freestringlist(&opts->cmpvarlist, opts->ncmpvarlist);

    newstringlist(&opts->cmpvarlist, &opts->ncmpvarlist, (int) NC_MAX_VARS);
    if (opts->variable) {
        if (opts->verbose)
            printf("INFO: Using variables provided in list.\n");

        status = strlistu(opts->variablelist, opts->cmpvarlist, opts->cmpvarlist,
                opts->nvariable, opts->ncmpvarlist, opts->ncmpvarlist);
    } else if (opts->exclude) {
        if (opts->verbose)
            printf("INFO: Excluding variables in provided list.\n");

        status = excludevars(ncid1, ncid2, opts->cmpvarlist, opts->ncmpvarlist, opts->excludelist, opts->nexclude);
    } else {
        if (opts->verbose)
            printf("INFO: Using all variables.\n");

        status = allvarnames(opts->cmpvarlist, opts->ncmpvarlist, ncid1, ncid2);
    }

    opts->ncmpvarlist = getnumstrlist(opts->cmpvarlist, (int) NC_MAX_VARS);
    if (opts->verbose) {
        printf("INFO: Variables to compare (%d):\n", opts->ncmpvarlist);
        for (i = 0; i < opts->ncmpvarlist - 1; ++i)
            printf("%s, ", opts->cmpvarlist[i]);

        if (opts->ncmpvarlist)
            printf("%s\n", opts->cmpvarlist[i]);
    }

    return status;
}
/* *********************************************************** */

/* Gets list of all variable names in both input files. */
int
allvarnames(char** list, int nvars, int ncid1, int ncid2) {
    char** tmplist = NULL;
    int ntmp;

    newstringlist(&tmplist, &ntmp, NC_MAX_VARS);

    if (ncallvars(ncid1, tmplist, ntmp))
        return 1;

    /* printf("%d: ncallvars returned.\n", __LINE__); */

    if (strlistu(tmplist, list, list, ntmp, nvars, nvars))
        return 1;

    /* printf("%d: Variable names from file 1 collected.\n", __LINE__); */
    clearstringlist(tmplist, NC_MAX_VARS);

    if (ncallvars(ncid2, tmplist, ntmp))
        return 1;

    if (strlistu(tmplist, list, list, ntmp, nvars, nvars))
        return 1;

    /* printf("%d: Variable names from file 2 collected.\n", __LINE__); */
    freestringlist(&tmplist, ntmp);

    return EXIT_SUCCESS;
}

/* *********************************************************** */

int
nccmpformats(nccmpopts* opts, int ncid1, int ncid2) {
    int status, fmt1, fmt2;

    status = nc_inq_format(ncid1, &fmt1);
    handle_error(status);

    status = nc_inq_format(ncid2, &fmt2);
    handle_error(status);

    if (fmt1 != fmt2) {
        fprintf(stdout, "DIFFER : FILE FORMATS : %s <> %s\n", NCFORMATSTR(fmt1), NCFORMATSTR(fmt2));

        if (!opts->warn[NCCMP_W_ALL] &&
                !opts->warn[NCCMP_W_FORMAT])
            return EXIT_DIFFER;
    }

    return EXIT_SUCCESS;
}

/* *********************************************************** */

int
nccmpglobalatts(nccmpopts* opts, int ncid1, int ncid2) {
    int ngatts1, ngatts2, nattsex1, nattsex2, i, status, status2, attid1, attid2;
    nc_type type1, type2;
    size_t len1, len2;
    char name1[NC_MAX_NAME], name2[NC_MAX_NAME];
    char** processedatts = NULL;
    char typestr1[4096], typestr2[4096];
    const string hstry("history");

    status = EXIT_SUCCESS;
    status2 = status;
    if (!opts->global)
        return status;

    if (opts->history == 0)
        appendstringtolist(&opts->globalexclude, hstry.c_str(), &opts->nglobalexclude);

    /*  
    printf("globalexclude =");
    printstrlist(opts->globalexclude, opts->nglobalexclude, stdout);
     */

    /* Number of global atts to compare with exclusion taken into account. */
    nattsex1 = 0;
    nattsex2 = 0;

    status = nc_inq_natts(ncid1, &ngatts1);
    handle_error(status);

    status = nc_inq_natts(ncid2, &ngatts2);
    handle_error(status);

    for (i = 0; i < ngatts1; ++i) {
        attid1 = i;
        status = nc_inq_attname(ncid1, NC_GLOBAL, attid1, name1);
        handle_error(status);

        if (!instringlist(opts->globalexclude, name1, opts->nglobalexclude)) {
            ++nattsex1;
        }
    }

    for (i = 0; i < ngatts2; ++i) {
        attid2 = i;
        status = nc_inq_attname(ncid2, NC_GLOBAL, attid2, name2);
        handle_error(status);

        if (!instringlist(opts->globalexclude, name2, opts->nglobalexclude)) {
            ++nattsex2;
        }
    }

    if (nattsex1 != nattsex2) {
        fprintf(stdout, "DIFFER : NUMBER OF GLOBAL ATTRIBUTES : %d <> %d\n", nattsex1, nattsex2);

        if (!opts->warn[NCCMP_W_ALL])
            status2 = EXIT_DIFFER;

        if (!opts->force) return status2;
    }

    if (newstringlist(&processedatts, &i, NC_MAX_VARS)) {
        fprintf(stdout, "ERROR: Failed to allocate string list for comparing  global attributes.\n");
        return EXIT_FATAL;
    }

    for (i = 0; i < ngatts1; ++i) {
        attid1 = i;
        status = nc_inq_attname(ncid1, NC_GLOBAL, attid1, name1);
        handle_error(status);

        /* Log that this gatt was processed. */
        addstringtolist(processedatts, name1, NC_MAX_VARS);

        if (instringlist(opts->globalexclude, name1, opts->nglobalexclude))
            continue;

        status = nc_inq_att(ncid1, NC_GLOBAL, name1, &type1, &len1);
        if (status != NC_NOERR) {
            fprintf(stderr, "Query failed on global attribute in %s\n", opts->file1);
            if (!opts->warn[NCCMP_W_ALL])
                status2 = EXIT_DIFFER;

            if (opts->force) continue;
            else return status2;
        }

        status = nc_inq_att(ncid2, NC_GLOBAL, name1, &type2, &len2);
        if (status != NC_NOERR) {
            fprintf(stdout, "DIFFER : NAME OF GLOBAL ATTRIBUTE : %s : GLOBAL ATTRIBUTE DOESN'T EXIST IN \"%s\"\n", name1, opts->file2);
            if (!opts->warn[NCCMP_W_ALL])
                status2 = EXIT_DIFFER;

            if (opts->force) continue;
            else if (!opts->force) return status2;
        }

        if (type1 != type2) {
            type2string(type1, typestr1);
            type2string(type2, typestr2);
            fprintf(stdout, "DIFFER : GLOBAL ATTRIBUTE TYPES : %s : %s <> %s\n", name1, typestr1, typestr2);
            if (!opts->warn[NCCMP_W_ALL])
                status2 = EXIT_DIFFER;

            if (opts->force) continue;
            else if (!opts->force) return status2;
        }

        if (len1 != len2) {
            prettyprintatt(ncid1, NULL, NC_GLOBAL, name1, typestr1);
            prettyprintatt(ncid2, NULL, NC_GLOBAL, name1, typestr2);
            fprintf(stdout, "DIFFER : LENGTHS OF GLOBAL ATTRIBUTE : %s : %lu <> %lu : VALUES : %s <> %s\n", name1,
                    (unsigned long) len1, (unsigned long) len2, typestr1, typestr2);
            if (!opts->warn[NCCMP_W_ALL])
                status2 = EXIT_DIFFER;

            if (opts->force) continue;
            else if (!opts->force) return status2;
        }

        if (cmpattval(ncid1, ncid2, NC_GLOBAL, NC_GLOBAL, name1, len1, type1) != NC_NOERR) {
            /* Pretty print values. */
            prettyprintatt(ncid1, NULL, NC_GLOBAL, name1, typestr1);
            prettyprintatt(ncid2, NULL, NC_GLOBAL, name1, typestr2);
            fprintf(stdout, "DIFFER : VALUES OF GLOBAL ATTRIBUTE : %s : %s <> %s\n", name1, typestr1, typestr2);
            if (!opts->warn[NCCMP_W_ALL])
                status2 = EXIT_DIFFER;

            if (opts->force) continue;
            else if (!opts->force) return status2;
        }
    }

    for (i = 0; i < ngatts2; ++i) {
        attid2 = i;
        status = nc_inq_attname(ncid2, NC_GLOBAL, attid2, name2);
        if (status != NC_NOERR) {
            fprintf(stderr, "Query failed for global attribute name in %s\n", opts->file2);
            if (!opts->warn[NCCMP_W_ALL])
                status2 = EXIT_DIFFER;

            if (opts->force) continue;
            else if (!opts->force) return status2;
        }

        /* Skip if already processed (or excluded). */
        if (instringlist(processedatts, name2, NC_MAX_VARS))
            continue;

        /* Log that this att was processed. */
        addstringtolist(processedatts, name2, NC_MAX_VARS);

        status = nc_inq_att(ncid2, NC_GLOBAL, name2, &type2, &len2);
        if (status != NC_NOERR) {
            fprintf(stderr, "Query failed on global attribute in %s\n", opts->file2);
            if (!opts->warn[NCCMP_W_ALL])
                status2 = EXIT_DIFFER;

            if (opts->force) continue;
            else if (!opts->force) return status2;
        }

        status = nc_inq_att(ncid1, NC_GLOBAL, name2, &type1, &len1);
        if (status != NC_NOERR) {
            fprintf(stdout, "DIFFER : NAME OF GLOBAL ATTRIBUTE : %s : GLOBAL ATTRIBUTE DOESN'T EXIST IN %s\n", name2, opts->file1);
            if (!opts->warn[NCCMP_W_ALL])
                status2 = EXIT_DIFFER;

            if (opts->force) continue;
            else if (!opts->force) return status2;
        }

        if (type1 != type2) {
            type2string(type1, typestr1);
            type2string(type2, typestr2);
            fprintf(stdout, "DIFFER : GLOBAL ATTRIBUTE TYPE : %s : %s <> %s\n", name1, typestr1, typestr2);
            if (!opts->warn[NCCMP_W_ALL])
                status2 = EXIT_DIFFER;

            if (opts->force) continue;
            else if (!opts->force) return status2;
        }

        if (len1 != len2) {
            prettyprintatt(ncid1, NULL, NC_GLOBAL, name1, typestr1);
            prettyprintatt(ncid2, NULL, NC_GLOBAL, name1, typestr2);

            fprintf(stdout, "DIFFER : LENGTHS OF GLOBAL ATTRIBUTE : %s : %lu <> %lu : VALUES : ", name1, (unsigned long) len1, (unsigned long) len2);

            switch (type1) {
            case NC_CHAR:
                /* Quote strings. */
                fprintf(stdout, "\"%s\" : \"%s\"\n", typestr1, typestr2);
                if (strcmp(typestr1, typestr2) == 0) {
                    /* Same text, but invisible trailing nulls because lengths differ. */
                    if (opts->warn[NCCMP_W_EOS] || opts->warn[NCCMP_W_ALL]) {
                        /* Pass */
                    } else {
                        status2 = EXIT_DIFFER;
                        if (opts->force) continue;
                        else if (!opts->force) return status2;
                    }
                }
                break;
            default:
                /* No quotes. */
                fprintf(stdout, "%s : %s\n", typestr1, typestr2);
                if (!opts->warn[NCCMP_W_ALL]) {
                    status2 = EXIT_DIFFER;
                    if (opts->force) continue;
                    else if (!opts->force) return status2;
                }
                break;
            }
        }

        if (cmpattval(ncid1, ncid2, NC_GLOBAL, NC_GLOBAL, name1, len1, type1) != NC_NOERR) {
            /* Pretty print values. */
            prettyprintatt(ncid1, NULL, NC_GLOBAL, name1, typestr1);
            prettyprintatt(ncid2, NULL, NC_GLOBAL, name1, typestr2);
            fprintf(stdout, "DIFFER : VALUES OF GLOBAL ATTRIBUTE : %s : %s <> %s\n", name1, typestr1, typestr2);
            if (!opts->warn[NCCMP_W_ALL])
                status2 = EXIT_DIFFER;

            if (opts->force) continue;
            else if (!opts->force) return status2;
        }
    }

    /* Clear the list. */
    freestringlist(&processedatts, NC_MAX_VARS);
    processedatts = NULL;

    return status2;
}

/* *********************************************************** */
int
nccmpmetadata(nccmpopts *opts, int ncid1, int ncid2) {
    int i, j, j1, j2, status, ncstatus, dimid1, dimid2, tmp1, tmp2, attid1, attid2, natts1, natts2;
    size_t len1, len2;
    char name1[NC_MAX_NAME], name2[NC_MAX_NAME], recname1[NC_MAX_NAME], recname2[NC_MAX_NAME], typestr1[1024], typestr2[1024];
    char** processedatts = NULL;

    status = EXIT_SUCCESS;

    if (opts->verbose)
        printf("INFO: Comparing metadata.\n");

    if (opts->verbose)
        printf("INFO: Comparing number of dimensions.\n");

    if (ndims1 != ndims2) {
        fprintf(stdout, "DIFFER : NUMBER OF DIMENSIONS IN FILES : %d <> %d\n", ndims1, ndims2);
        if (!opts->warn[NCCMP_W_ALL])
            status = EXIT_DIFFER;

        if (!opts->force)
            return status;
    }

    if (opts->verbose)
        printf("INFO: Getting record dimension names, if they exist.\n");

    if (recid1 != -1) {
        ncstatus = nc_inq_dimname(ncid1, recid1, recname1);
        handle_error(ncstatus);
    } else
        strcpy(recname1, "");

    if (recid2 != -1) {
        ncstatus = nc_inq_dimname(ncid2, recid2, recname2);
        handle_error(ncstatus);
    } else
        strcpy(recname2, "");

    /* Dimensions */
    if (opts->verbose)
        printf("INFO: Comparing dimension lengths.\n");

    for (i = 0; i < ndims1; ++i) {
        dimid1 = dims1[i].dimid;
        ncstatus = nc_inq_dim(ncid1, dimid1, name1, &len1);
        if (ncstatus != NC_NOERR) {
            if (!opts->warn[NCCMP_W_ALL])
                status = EXIT_DIFFER;

            fprintf(stderr, "Failed to query dimension id %d in file %s.\n", dimid1, opts->file1);
            if (opts->force) continue;
            else return status;
        }

        ncstatus = nc_inq_dimid(ncid2, name1, &dimid2);
        if (ncstatus != NC_NOERR) {
            fprintf(stdout, "DIFFER : NAME : DIMENSION : %s : DIMENSION DOESN'T EXIST IN \"%s\"\n", name1, opts->file2);

            if (!opts->warn[NCCMP_W_ALL])
                status = EXIT_DIFFER;

            if (opts->force) continue;
            else return status;
        }

        ncstatus = nc_inq_dim(ncid2, dimid2, name2, &len2);
        if (ncstatus != NC_NOERR) {
            if (!opts->warn[NCCMP_W_ALL])
                status = EXIT_DIFFER;

            fprintf(stderr, "Failed to query dimension \"%s\" in file \"%s\".\n", name1, opts->file2);
            if (opts->force) continue;
            else return status;
        }

        if (len1 != len2) {
            fprintf(stdout, "DIFFER : LENGTHS : DIMENSION : %s : %lu <> %lu\n", name1,
                    (unsigned long) len1, (unsigned long) len2);

            if (!opts->warn[NCCMP_W_ALL])
                status = EXIT_DIFFER;

            if (opts->force) continue;
            else return status;
        }
    }

    /* Variables */
    if (opts->verbose)
        printf("INFO: Comparing variable datatypes and rank.\n");

    for (i = 0; i < opts->ncmpvarlist; ++i) {
        j1 = findvar(opts->cmpvarlist[i], vars1);
        if (j1 == -1) {
            fprintf(stdout, "DIFFER : NAME : VARIABLE : %s%s : VARIABLE DOESN'T EXIST IN \"%s\"\n", getGroupPath(), opts->cmpvarlist[i], opts->file1);

            if (!opts->warn[NCCMP_W_ALL])
                status = EXIT_DIFFER;

            if (opts->force) continue;
            else goto recover;
        }

        j2 = findvar(opts->cmpvarlist[i], vars2);
        if (j2 == -1) {
            fprintf(stdout, "DIFFER : NAME : VARIABLE : %s%s : VARIABLE DOESN'T EXIST IN \"%s\"\n", getGroupPath(), opts->cmpvarlist[i], opts->file2);

            if (!opts->warn[NCCMP_W_ALL])
                status = EXIT_DIFFER;

            if (opts->force) continue;
            else goto recover;
        }

        if (vars1[j1].type != vars2[j2].type) {
            type2string(vars1[j1].type, typestr1);
            type2string(vars2[j2].type, typestr2);
            fprintf(stdout, "DIFFER : TYPES : VARIABLE : %s%s : %s <> %s\n", getGroupPath(), opts->cmpvarlist[i], typestr1, typestr2);

            if (!opts->warn[NCCMP_W_ALL])
                status = EXIT_DIFFER;

            if (opts->force) continue;
            else goto recover;
        }

        if (vars1[j1].ndims != vars2[j2].ndims) {
            fprintf(stdout, "DIFFER : NUMBER : DIMENSIONS : VARIABLE : %s%s : %d <> %d\n", getGroupPath(), opts->cmpvarlist[i], ndims1, ndims2);

            if (!opts->warn[NCCMP_W_ALL])
                status = EXIT_DIFFER;

            if (opts->force) continue;
            else goto recover;
        }
    }

    /*printf("DEBUG : %d : \n", __LINE__);*/
    if (opts->verbose)
        printf("INFO: Comparing variables' dimension names.\n");

    for (i = 0; i < opts->ncmpvarlist; ++i) {
        j1 = findvar(opts->cmpvarlist[i], vars1);
        j2 = findvar(opts->cmpvarlist[i], vars2);

        if ((j1 == -1) || (j2 == -1))
            continue;

        /* dimensions */
        for (j = 0; j < vars1[j1].ndims; ++j) {
            dimid1 = vars1[j1].dimids[j];

            if (j < vars2[j2].ndims)
                dimid2 = vars2[j2].dimids[j];
            else
                break;

            /*printf("DEBUG : %d : %s,  %s, %s\n", __LINE__, opts->cmpvarlist[i], dims1[dimid1].name, dims2[dimid2].name);*/

            if (strcmp(dims1[dimid1].name, dims2[dimid2].name) != 0) {
                fprintf(stdout, "DIFFER : DIMENSION NAMES FOR VARIABLE %s%s : %s <> %s\n", getGroupPath(), opts->cmpvarlist[i], dims1[dimid1].name, dims2[dimid2].name);

                if (!opts->warn[NCCMP_W_ALL])
                    status = EXIT_DIFFER;

                if (opts->force) continue;
                else goto recover;
            }

            tmp1 = strcmp(dims1[dimid1].name, recname1);
            tmp2 = strcmp(dims2[dimid2].name, recname2);

            if ((tmp1 == 0) && (tmp2 != 0)) {
                fprintf(stdout, "DIFFER : VARIABLE : %s%s : DIMENSION %s IS RECORD IN FILE \"%s\" BUT NOT IN \"%s\"\n", getGroupPath(), vars1[j1].name, dims1[dimid1].name, opts->file1, opts->file2);

                if (!opts->warn[NCCMP_W_ALL])
                    status = EXIT_DIFFER;

                if (opts->force) continue;
                else goto recover;
            } else if ((tmp1 != 0) && (tmp2 == 0)) {
                fprintf(stdout, "DIFFER : VARIABLE : %s%s : DIMENSION %s IS RECORD IN FILE \"%s\" BUT NOT IN \"%s\"\n", getGroupPath(), vars1[j1].name, dims2[dimid2].name, opts->file2, opts->file1);

                if (!opts->warn[NCCMP_W_ALL])
                    status = EXIT_DIFFER;

                if (opts->force) continue;
                else goto recover;
            }
        }
    }

    if (opts->verbose)
        printf("INFO: Comparing variables' attributes.\n");

    /*printf("DEBUG : %d : \n", __LINE__);*/
    /* Pass in 'i' as dummy. */
    if (newstringlist(&processedatts, &i, NC_MAX_VARS)) {
        fprintf(stderr, "ERROR: Failed to allocate string list for comparing attributes.\n");
        return EXIT_FATAL;
    }
    /*printf("DEBUG : %d : \n", __LINE__); fflush(stdout);*/

    for (i = 0; i < opts->ncmpvarlist; ++i) {
        /*printf("DEBUG : %d : i = %d\n", __LINE__, i); fflush(stdout);*/
        j1 = findvar(opts->cmpvarlist[i], vars1);
        j2 = findvar(opts->cmpvarlist[i], vars2);

        if ((j1 == -1) || (j2 == -1))
            continue;

        natts1 = natts2 = 0;
        clearstringlist(processedatts, NC_MAX_VARS);

        /*printf("DEBUG : %d : var=%s\n", __LINE__, opts->cmpvarlist[i]); fflush(stdout);*/

        /* Attributes */
        for (attid1 = 0; attid1 < vars1[j1].natts; ++attid1) {
            ncstatus = nc_inq_attname(ncid1, vars1[j1].varid, attid1, name1);
            if (ncstatus != NC_NOERR) {
                if (!opts->warn[NCCMP_W_ALL])
                    status = EXIT_DIFFER;

                if (opts->force) continue;
                else goto recover;
            }

            if (instringlist(opts->excludeattlist, name1, opts->nexcludeatt) || instringlist(processedatts, name1, NC_MAX_VARS))
                continue;

            /* Log that this att was processed. */
            addstringtolist(processedatts, name1, NC_MAX_VARS);
            ++natts1;
            /*printf("natts1 %s, %d\n", name1, natts1);*/
            ncstatus = cmpatt(ncid1, ncid2, vars1[j1].varid, vars2[j2].varid, name1, vars1[j1].name, opts);
            if (ncstatus == EXIT_DIFFER) {
                status = ncstatus;
                if (opts->force) continue;
                else goto recover;
            }
        }

        /*printf("DEBUG : %d : \n", __LINE__);*/
        for (attid2 = 0; attid2 < vars2[j2].natts; ++attid2) {
            ncstatus = nc_inq_attname(ncid2, vars2[j2].varid, attid2, name2);
            if (ncstatus != NC_NOERR) {
                fprintf(stderr, "Failed to query variable %s%s attribute in file \"%s\"\n", getGroupPath(), vars2[j2].name, opts->file2);

                if (!opts->warn[NCCMP_W_ALL])
                    status = EXIT_DIFFER;

                if (opts->force) continue;
                else goto recover;
            }

            if (instringlist(opts->excludeattlist, name2, opts->nexcludeatt))
                continue;

            /* Count non-excluded attribute. */
            ++natts2;
            /*printf("natts2 %s, %d\n", name2, natts2);*/
            if (instringlist(processedatts, name2, NC_MAX_VARS))
                continue;

            /* Log that this att was processed. */
            addstringtolist(processedatts, name2, NC_MAX_VARS);

            /* Do comparison. */
            ncstatus = cmpatt(ncid1, ncid2, vars1[j1].varid, vars2[j2].varid, name2, vars2[j2].name, opts);
            if (ncstatus == EXIT_DIFFER) {
                status = ncstatus;
                if (opts->force) continue;
                else goto recover;
            }
        }

        if (natts1 != natts2) {
            fprintf(stdout, "DIFFER : NUMBER OF ATTRIBUTES : VARIABLE : %s%s : %d <> %d\n", getGroupPath(), opts->cmpvarlist[i], natts1, natts2);
            if (!opts->warn[NCCMP_W_ALL])
                status = EXIT_DIFFER;

            if (opts->force) continue;
            else goto recover;
        }
    }

recover:
    freestringlist(&processedatts, NC_MAX_VARS);
    processedatts = NULL;


    return status;
}
/* *********************************************************** */

/* Returns index into varstruct array if found using name otherwise -1. */
int findvar(char * name, varstruct *vars) {
    int i;
    for (i = 0; i < NC_MAX_VARS; ++i) {
        if (strcmp(name, vars[i].name) == 0)
            return i;
    }

    return -1;
}

/* *********************************************************** */

/* Do the comparision of variables.
   Record index (rec) is optional; -1 if not applicable.
   Returns comparison result success or failure.
 */
int cmpvar(char* name, int rec, nccmpopts* opts, int ncid1, int ncid2, nccmp_user_type_t *user_types1, nccmp_user_type_t *user_types2) {

    varstruct *v1, *v2;
    size_t idx1, idx2, status, i;
    int j;
    size_t start[NC_MAX_DIMS], count[NC_MAX_DIMS];
    off_t nitems;
    size_t odomax[NC_MAX_DIMS];
    int diffstatus = EXIT_SUCCESS;
    unsigned char *data1 = NULL, *data2 = NULL;

    idx1 = findvar(name, vars1);
    idx2 = findvar(name, vars2);

    v1 = &vars1[idx1];
    v2 = &vars2[idx2];

    /*printf("DEBUG : %s len : %d <> %d\n", name, v1->len, v2->len); */

    if (v1->len != v2->len) {
        fprintf(stdout, "DIFFER : SIZE OF VARIABLE \"%s%s\" : %d <> %d\n", getGroupPath(), name, (int) v1->len, (int) v2->len);

        if (!opts->warn[NCCMP_W_ALL]) \
				return EXIT_DIFFER;
        else
            return EXIT_SUCCESS;
    }

    // no need to compare if the length is 0, since netCDF will blow up
    if(v1->len == 0) {
        if (opts->verbose)
            printf("INFO: 0 length variables are equal.\n");
        return 0;
    }

    for (j = 0; j < v1->ndims; ++j) {
        start[j] = 0;
        odomax[j] = v1->dimlens[j] - 1;
    }

#ifdef __DEBUG__
    printf("DEBUG : %d : odomax = ", __LINE__);
    for (i = 0; i < v1->ndims; ++i) {
        printf("%d ", odomax[i]);
    }
    printf("\n");
#endif

    /* If has record dim. */
    if (v1->hasrec && (rec >= 0))
        start[0] = rec;

    /* Read in slab for last dimension at-a-time only.
We'll then loop over all the other outer dimensions. */
    for (j = 0; j < v1->ndims - 1; ++j) {
        count[j] = 1;
    }

    /* We'll always read in entire last dimension
   except if only dimension is record. */
    if ((v1->ndims == 1) && (v1->hasrec)) {
        nitems = 1;
    } else {
        nitems = v1->dimlens[v1->ndims - 1];
    }

    count[v1->ndims - 1] = nitems;

    /*printf("DEBUG : %d : nitems = %d\n", __LINE__, nitems);\*/
    int ut_nidx1 = v1->user_type_idx;
    int ut_nidx2 = v2->user_type_idx;
    int offset1, offset2, size1, size2;

    if (ut_nidx1 < 0) {

        switch (v1->type) {
        case NC_BYTE:
            strcpy(opts->tprintf, "%d");
            //cmp_var<int8_t>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.b, v2->missing.b);
            opts->tolerance = 0.0;
            diffstatus = cmp_vartol<int8_t>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.f, v2->missing.f);
            break;

        case NC_CHAR:
            strcpy(opts->tprintf, "%c");
            // cmp_var<char>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.c, v2->missing.c);
            opts->tolerance = 0.0;
            diffstatus = cmp_vartol<char>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.f, v2->missing.f);
            break;

        case NC_SHORT:
            strcpy(opts->tprintf, "%d");
            // cmp_var<short>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.s, v2->missing.s);
            opts->tolerance = 0.0;
            diffstatus = cmp_vartol<short>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.f, v2->missing.f);
            break;

        case NC_USHORT:
            strcpy(opts->tprintf, "%d");
            // cmp_var<short>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.s, v2->missing.s);
            opts->tolerance = 0.0;
            diffstatus = cmp_vartol<uint16_t>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.f, v2->missing.f);
            break;

        case NC_UBYTE:
            strcpy(opts->tprintf, "%d");
            // cmp_var<short>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.s, v2->missing.s);
            opts->tolerance = 0.0;
            diffstatus = cmp_vartol<uint8_t>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.f, v2->missing.f);
            break;

        case NC_INT:
            strcpy(opts->tprintf, "%d");
            opts->tolerance = 0.0;
            diffstatus = cmp_vartol<int>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.i, v2->missing.i);
            break;

        case NC_UINT:
            strcpy(opts->tprintf, "%d");
            opts->tolerance = 0.0;
            diffstatus = cmp_vartol<uint>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.ui, v2->missing.ui);
            break;

        case NC_INT64:
            strcpy(opts->tprintf, "%ld");
            opts->tolerance = 0.0;
            diffstatus = cmp_vartol<long>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.l, v2->missing.l);
            break;

        case NC_UINT64:
            strcpy(opts->tprintf, "%ld");
            opts->tolerance = 0.0;
            diffstatus = cmp_vartol<uint64_t>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.ul, v2->missing.ul);
            break;
            
        case NC_FLOAT:
            strcpy(opts->tprintf, opts->precision);
            if (opts->notolerance) {
                opts->tolerance = 0.0;
                diffstatus = cmp_vartol<float>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.f, v2->missing.f);
            } else {
                if (opts->verbose)
                    printf("INFO: Using maximum tolerance of 0.0001 for float comparisons.\n");
                opts->tolerance = 0.0001;
                diffstatus = cmp_vartol<float>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.f, v2->missing.f);
                opts->tolerance = 0;
            }
            break;

        case NC_DOUBLE:
            strcpy(opts->tprintf, opts->precision);
            if (opts->notolerance) {
                diffstatus = cmp_var<double>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.d, v2->missing.d);
            } else {
                if (opts->verbose)
                    printf("INFO: Using maximum tolerance of 0.0001 for double comparisons.\n");
                opts->tolerance = 0.0001;
                diffstatus = cmp_vartol<double>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.d, v2->missing.d);
                opts->tolerance = 0;
            }
            break;
        default:
            if (v1->type < NC_MAX_TYPES && !notsupported[v1->type]) {
                printf("INFO: This Type (%d) not supported .\n", v1->type);
                notsupported[v1->type] = 1;
            } else {
                if (v1->type >= NC_MAX_TYPES) {
                    printf("INFO: This Type (%d) not supported ... ", v1->type);
                    printf("INFO: Increase NC_MAX_TYPES > %d to avoid this message. \n", v1->type);
                }
            }
            diffstatus = -1;
            break;

        }
    } else {
        data1 = (unsigned char *) malloc(sizeof (unsigned char)*user_types1[ut_nidx1].size * (nitems + 1));
        data2 = (unsigned char *) malloc(sizeof (unsigned char)*user_types2[ut_nidx2].size * (nitems + 1));
        status = nc_get_vara(ncid1, v1->varid, start, count, data1);
        handle_error(status);
        status = nc_get_vara(ncid2, v2->varid, start, count, data2);
        handle_error(status);

        size1 = user_types1[ut_nidx1].root_size;
        size2 = user_types2[ut_nidx2].root_size;
        for (i = 0; i < user_types1[ut_nidx1].num_fields; ++i) {
            offset1 = user_types1[ut_nidx1].fields[i].offset;
            offset2 = user_types2[ut_nidx2].fields[i].offset;
            switch (user_types1[ut_nidx1].fields[i].type_id) {
            case NC_BYTE:
                strcpy(opts->tprintf, "%d");
                //cmp_var<int8_t>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.b, v2->missing.b);
                opts->tolerance = 0.0;
                diffstatus = cmp_vartol_ut<int8_t>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.b, v2->missing.b);
                break;

            case NC_CHAR:
                strcpy(opts->tprintf, "%c");
                // cmp_var<char>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.c, v2->missing.c);
                opts->tolerance = 0.0;
                diffstatus = cmp_vartol_ut<char>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.c, v2->missing.c);
                break;

            case NC_SHORT:
                strcpy(opts->tprintf, "%d");
                // cmp_var<short>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.s, v2->missing.s);
                opts->tolerance = 0.0;
                diffstatus = cmp_vartol<short>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.s, v2->missing.s);
                break;

            case NC_USHORT:
                strcpy(opts->tprintf, "%d");
                // cmp_var<short>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.s, v2->missing.s);
                opts->tolerance = 0.0;
                diffstatus = cmp_vartol<uint16_t>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.us, v2->missing.us);
                break;

            case NC_UBYTE:
                strcpy(opts->tprintf, "%d");
                // cmp_var<short>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.s, v2->missing.s);
                opts->tolerance = 0.0;
                diffstatus = cmp_vartol_ut<uint8_t>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.ub, v2->missing.ub);
                break;

            case NC_INT:
                strcpy(opts->tprintf, "%d");
                opts->tolerance = 0.0;
                diffstatus = cmp_vartol_ut<int>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.i, v2->missing.i);
                break;

            case NC_UINT:
                strcpy(opts->tprintf, "%d");
                opts->tolerance = 0.0;
                diffstatus = cmp_vartol_ut<uint>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.ui, v2->missing.ui);
                break;

             case NC_INT64:
                strcpy(opts->tprintf, "%ld");
                opts->tolerance = 0.0;
                diffstatus = cmp_vartol_ut<long>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.l, v2->missing.l);
                break;

            case NC_UINT64:
                strcpy(opts->tprintf, "%ld");
                opts->tolerance = 0.0;
                diffstatus = cmp_vartol_ut<uint64_t>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.ul, v2->missing.ul);
                break;

           case NC_FLOAT:
                strcpy(opts->tprintf, opts->precision);
                if (opts->notolerance) {
                    opts->tolerance = 0.0;
                    diffstatus = cmp_vartol_ut<float>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.f, v2->missing.f);
                } else {
                    if (opts->verbose)
                        printf("INFO: Using maximum tolerance of 0.0001 for float comparisons.\n");
                    opts->tolerance = 0.0001;
                    diffstatus = cmp_vartol_ut<float>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.f, v2->missing.f);
                    opts->tolerance = 0;
                }
                break;

            case NC_DOUBLE:
                strcpy(opts->tprintf, opts->precision);
                if (opts->notolerance) {
                    diffstatus = cmp_vartol_ut<double>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.d, v2->missing.d);
                } else {
                    if (opts->verbose)
                        printf("INFO: Using maximum tolerance of 0.0001 for double comparisons.\n");
                    opts->tolerance = 0.0001;
                    diffstatus = cmp_vartol_ut<double>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.d, v2->missing.d);
                    opts->tolerance = 0;
                }
                break;
            default:
                if (v1->type < NC_MAX_TYPES && !notsupported[user_types1[ut_nidx1].fields[i].type_id]) {
                    printf("INFO: This Type (%d) not supported .\n", user_types1[ut_nidx1].fields[i].type_id);
                    notsupported[v1->type] = 1;
                } else {
                    if (v1->type >= NC_MAX_TYPES) {
                        printf("INFO: This Type (%d) not supported ... ", user_types1[ut_nidx1].fields[i].type_id);
                        printf("INFO: Increase NC_MAX_TYPES > %d to avoid this message. \n", user_types1[ut_nidx1].fields[i].type_id);
                    }
                }
                diffstatus = -1;
                break;

            }
        }
    }
    return diffstatus;
}
/* *********************************************************** */

/* Do the comparision of variables.
   Record index (rec) is optional; -1 if not applicable.
   Returns comparison result success or failure.
 */
int cmpvartol(char* name, int rec, nccmpopts* opts, int ncid1, int ncid2, nccmp_user_type_t *user_types1, nccmp_user_type_t *user_types2) {

    varstruct *v1, *v2;
    size_t idx1, idx2, status, i;
    int j;
    size_t start[NC_MAX_DIMS], count[NC_MAX_DIMS];
    off_t nitems;
    size_t odomax[NC_MAX_DIMS];
    int diffstatus = EXIT_SUCCESS;
    unsigned char *data1 = NULL, *data2 = NULL;

    if (opts->verbose) {
        if (rec != -1)
            printf("INFO: Comparing data for variable \"%s%s\" at record %d.\n", getGroupPath(), name, (int) rec);
        else
            printf("INFO: Comparing non-record data for variable \"%s%s\".\n", getGroupPath(), name);
    }

    idx1 = findvar(name, vars1);
    idx2 = findvar(name, vars2);

    if (idx1 < 0) {
        if (!opts->metadata) /* This gets reported in cmpmeta. */
            fprintf(stdout, "DIFFER : Failed to find variable \"%s%s\" in file \"%s\".\n", getGroupPath(), name, opts->file1);

        if (!opts->warn[NCCMP_W_ALL])
            return EXIT_DIFFER;
        else
            return EXIT_SUCCESS;
    }

    if (idx2 < 0) {
        if (!opts->metadata) /* This gets reported in cmpmeta. */
            fprintf(stdout, "DIFFER : Failed to find variable \"%s%s\" in file \"%s\".\n", getGroupPath(), name, opts->file2);

        if (!opts->warn[NCCMP_W_ALL])
            return EXIT_DIFFER;
        else
            return EXIT_SUCCESS;
    }

    v1 = &vars1[idx1];
    v2 = &vars2[idx2];

    if (v1->len != v2->len) {
        fprintf(stdout, "DIFFER : SIZE OF VARIABLE \"%s%s\" : %d <> %d\n", getGroupPath(), name, (int) v1->len, (int) v2->len);

        if (!opts->warn[NCCMP_W_ALL]) \
				return EXIT_DIFFER;
        else
            return EXIT_SUCCESS;
    }

    int ut_nidx1 = v1->user_type_idx;
    int ut_nidx2 = v2->user_type_idx;
    int offset1, offset2, size1, size2;

    for (j = 0; j < v1->ndims; ++j) {
        start[j] = 0;
        odomax[j] = v1->dimlens[j] - 1;
    }

    /* If has record dim. */
    if (v1->hasrec && (rec >= 0))
        start[0] = rec;

    /* Read in slab for last dimension at-a-time only. 
      We'll then loop over all the other outer dimensions. */
    for (j = 0; j < v1->ndims - 1; ++j) {
        count[j] = 1;
    }

    /* We'll always read in entire last dimension
       except if only dimension is record. */
    if ((v1->ndims == 1) && (v1->hasrec)) {
        nitems = 1;
    } else
        nitems = v1->dimlens[v1->ndims - 1];

    count[v1->ndims - 1] = nitems;

    /* todo: make cmpvar and cmpvartol same function to re-use code immediately above; just use conditional to choose CMP_VAR or CMP_VARTOL macro. */

    if (ut_nidx1 < 0) {
        switch (v1->type) {
        case NC_BYTE:
            strcpy(opts->tprintf, "%d");
            diffstatus = cmp_vartol<int8_t>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.b, v2->missing.b);
            break;

        case NC_UBYTE:
            strcpy(opts->tprintf, "%d");
            diffstatus = cmp_vartol<uint8_t>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.ub, v2->missing.ub);
            break;

        case NC_CHAR:
            strcpy(opts->tprintf, "%c");
            diffstatus = cmp_vartol<char>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.c, v2->missing.c);
            break;

        case NC_SHORT:
            strcpy(opts->tprintf, "%d");
            diffstatus = cmp_vartol<short>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.s, v2->missing.s);
            break;

        case NC_USHORT:
            strcpy(opts->tprintf, "%d");
            diffstatus = cmp_vartol<uint16_t>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.us, v2->missing.us);
            break;

        case NC_INT:
            strcpy(opts->tprintf, "%d");
            diffstatus = cmp_vartol<int>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.i, v2->missing.i);
            break;

        case NC_UINT:
            strcpy(opts->tprintf, "%d");
            diffstatus = cmp_vartol<uint>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.ui, v2->missing.ui);
            break;

        case NC_INT64:
            strcpy(opts->tprintf, "%ld");
            diffstatus = cmp_vartol<long>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.l, v2->missing.l);
            break;

        case NC_UINT64:
            strcpy(opts->tprintf, "%ld");
            diffstatus = cmp_vartol<uint64_t>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.ul, v2->missing.ul);
            break;

        case NC_FLOAT:
            strcpy(opts->tprintf, "%g");
            diffstatus = cmp_vartol<float>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.f, v2->missing.f);
            break;

        case NC_DOUBLE:
            strcpy(opts->tprintf, "%g");
            diffstatus = cmp_vartol<double>(ncid1, ncid2, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.d, v2->missing.d);
            break;
        default:
            if (v1->type < NC_MAX_TYPES && !notsupported[v1->type]) {
                printf("INFO: This Type (%d) not supported .\n", v1->type);
                notsupported[v1->type] = 1;
            } else {
                if (v1->type >= NC_MAX_TYPES) {
                    printf("INFO: This Type (%d) not supported ... ", v1->type);
                    printf("INFO: Increase NC_MAX_TYPES > %d to avoid this message. \n", v1->type);
                }
            }
            diffstatus = -1;
            break;

        }
    } else {
        data1 = (unsigned char *) malloc(sizeof (unsigned char)*user_types1[ut_nidx1].size * (nitems + 1));
        data2 = (unsigned char *) malloc(sizeof (unsigned char)*user_types2[ut_nidx2].size * (nitems + 1));
        status = nc_get_vara(ncid1, v1->varid, start, count, data1);
        handle_error(status);
        status = nc_get_vara(ncid2, v2->varid, start, count, data2);
        handle_error(status);

        size1 = user_types1[ut_nidx1].root_size;
        size2 = user_types2[ut_nidx2].root_size;
        for (i = 0; i < user_types1[ut_nidx1].num_fields; ++i) {
            offset1 = user_types1[ut_nidx1].fields[i].offset;
            offset2 = user_types2[ut_nidx2].fields[i].offset;

            switch (user_types1[ut_nidx1].fields[i].type_id) {
            case NC_BYTE:
                strcpy(opts->tprintf, "%d");
                diffstatus = cmp_vartol_ut<int8_t>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.b, v2->missing.b);

                break;

            case NC_UBYTE:
                strcpy(opts->tprintf, "%d");
                diffstatus = cmp_vartol_ut<uint8_t>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.ub, v2->missing.ub);

                break;

            case NC_CHAR:
                strcpy(opts->tprintf, "%c");
                diffstatus = cmp_vartol_ut<char>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.c, v2->missing.c);
                break;

            case NC_SHORT:
                strcpy(opts->tprintf, "%d");
                diffstatus = cmp_vartol_ut<short>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.s, v2->missing.s);
                break;

            case NC_USHORT:
                strcpy(opts->tprintf, "%d");
                diffstatus = cmp_vartol_ut<uint16_t>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.us, v2->missing.us);
                break;

            case NC_INT:
                strcpy(opts->tprintf, "%d");
                diffstatus = cmp_vartol_ut<int>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.i, v2->missing.i);
                break;

            case NC_UINT:
                strcpy(opts->tprintf, "%d");
                diffstatus = cmp_vartol_ut<uint>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.ui, v2->missing.ui);
                break;

            case NC_INT64:
                strcpy(opts->tprintf, "%ld");
                diffstatus = cmp_vartol_ut<long>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.l, v2->missing.l);
                break;

            case NC_UINT64:
                strcpy(opts->tprintf, "%ld");
                diffstatus = cmp_vartol_ut<uint64_t>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.ul, v2->missing.ul);
                break;

            case NC_FLOAT:
                strcpy(opts->tprintf, "%g");
                diffstatus = cmp_vartol_ut<float>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.f, v2->missing.f);
                break;

            case NC_DOUBLE:
                strcpy(opts->tprintf, "%g");
                diffstatus = cmp_vartol_ut<double>(data1, data2, offset1, offset2, size1, size2, user_types1[ut_nidx1].fields[i].name, opts, rec, odomax, nitems, count, start, v1, v2, v1->missing.d, v2->missing.d);
                break;
            default:
                if (v1->type < NC_MAX_TYPES && !notsupported[user_types1[ut_nidx1].fields[i].type_id]) {
                    printf("INFO: This Type (%d) not supported .\n", user_types1[ut_nidx1].fields[i].type_id);
                    notsupported[v1->type] = 1;
                } else {
                    if (v1->type >= NC_MAX_TYPES) {
                        printf("INFO: This Type (%d) not supported ... ", user_types1[ut_nidx1].fields[i].type_id);
                        printf("INFO: Increase NC_MAX_TYPES > %d to avoid this message. \n", user_types1[ut_nidx1].fields[i].type_id);
                    }
                }
                diffstatus = -1;
                break;

            }
        }
    }
    return diffstatus;
}

/* *********************************************************** */
int
nccmpdatarecvartol(int ncid1, int ncid2, char* varname, nccmpopts* opts,
        size_t recstart, size_t recend, nccmp_user_type_t *user_types1, nccmp_user_type_t *user_types2) {
    int cmpstatus;
    int status = EXIT_SUCCESS;

    for (; recstart <= recend; ++recstart) {
        cmpstatus = cmpvartol(varname, recstart, opts, ncid1, ncid2, user_types1, user_types2) || status;
        if (cmpstatus != EXIT_SUCCESS) {
            status = cmpstatus;
            if (opts->force)
                continue;
            else
                break;
        }
    }

    return status;
}

/* *********************************************************** */
int
nccmpdatarecvar(int ncid1, int ncid2, char* varname, nccmpopts* opts,
        size_t recstart, size_t recend, nccmp_user_type_t *user_types1, nccmp_user_type_t *user_types2) {
    int cmpstatus;
    int status = EXIT_SUCCESS;

    for (; recstart <= recend; ++recstart) {
        cmpstatus = cmpvar(varname, recstart, opts, ncid1, ncid2, user_types1, user_types2) || status;
        if (cmpstatus != EXIT_SUCCESS) {
            status = cmpstatus;
            if (opts->force)
                continue;
            else
                break;
        }
    }

    return status;
}

/* *********************************************************** */
int nccmpdatatol(int ncid1, int ncid2, nccmpopts* opts, nccmp_user_type_t *user_types1, nccmp_user_type_t *user_types2) {
    int i, idx1;
    int status, cmpstatus, nprocessed;
    char** processed = NULL;

    status = EXIT_SUCCESS;

    if (opts->verbose)
        printf("INFO: Comparing data with tolerance.\n");

    if (newstringlist(&processed, &nprocessed, NC_MAX_VARS)) {
        fprintf(stderr, "ERROR: Failed to allocate string list for comparing data.\n");
        return EXIT_FATAL;
    }

    for (i = 0; i < opts->ncmpvarlist; ++i) {
        if (instringlist(processed, opts->cmpvarlist[i], nprocessed))
            /* Skip varnames already processed. */
            continue;

        if (opts->verbose)
            printf("INFO: Comparing data for variable \"%s%s\".\n", getGroupPath(), opts->cmpvarlist[i]);

        addstringtolist(processed, opts->cmpvarlist[i], nprocessed);

        /* Has rec? */
        idx1 = findvar(opts->cmpvarlist[i], vars1);
        if (vars1[idx1].hasrec) {

            /* Compare only if # recs are equal and not zero. */
            if ((nrec1 == nrec2) && (nrec1 + nrec2)) {
                cmpstatus = nccmpdatarecvartol(ncid1, ncid2, opts->cmpvarlist[i], opts, 0, nrec1 - 1, user_types1, user_types2);
                if (cmpstatus) {
                    status = cmpstatus;
                    if (opts->force)
                        continue;
                    else
                        break;
                }
            }
        } else {
            cmpstatus = cmpvartol(opts->cmpvarlist[i], -1, opts, ncid1, ncid2, user_types1, user_types2);
            if (cmpstatus) {
                status = cmpstatus;
                if (opts->force)
                    continue;
                else
                    break;
            }
        }
    }

    if (opts->verbose)
        printf("INFO: Finished comparing data.\n");

    return status;
}

/* *********************************************************** */
int nccmpdata(nccmpopts* opts, int ncid1, int ncid2, nccmp_user_type_t *user_types1, nccmp_user_type_t *user_types2) {
    int i, idx1, idx2;
    int status, cmpstatus, nprocessed;
    char** processed = NULL;
    char str1[32], str2[32];

    status = EXIT_SUCCESS;

    if (opts->verbose)
        printf("INFO: Comparing data.\n");

    if (opts->tolerance != 0)
        return nccmpdatatol(ncid1, ncid2, opts, user_types1, user_types2);

    if (newstringlist(&processed, &nprocessed, NC_MAX_VARS)) {
        fprintf(stderr, "ERROR: Failed to allocate string list for comparing data.\n");
        return EXIT_FATAL;
    }

    for (i = 0; i < opts->ncmpvarlist; ++i) {

        // reset diff counter for each variable
        opts->diffcount = 0;
        opts->extentcount = 0;
    
        if (instringlist(processed, opts->cmpvarlist[i], nprocessed))
            /* Skip varnames already processed. */
            continue;

        if (opts->verbose)
            printf("INFO: Comparing data for variable \"%s%s\".\n", getGroupPath(), opts->cmpvarlist[i]);

        addstringtolist(processed, opts->cmpvarlist[i], nprocessed);

        idx1 = findvar(opts->cmpvarlist[i], vars1);
        idx2 = findvar(opts->cmpvarlist[i], vars2);

        if (idx1 < 0) {
            if (!opts->metadata) /* This gets reported in cmpmeta. */
                fprintf(stdout, "DIFFER : Failed to find variable \"%s%s\" in file \"%s\".\n", getGroupPath(), opts->cmpvarlist[i], opts->file1);

            if (!opts->warn[NCCMP_W_ALL])
                status = EXIT_DIFFER;
            if (opts->force)
                continue;
            else
                break;
        }

        if (idx2 < 0) {
            if (!opts->metadata) /* This gets reported in cmpmeta. */
                fprintf(stdout, "DIFFER : Failed to find variable \"%s%s\" in file \"%s\".\n", getGroupPath(), opts->cmpvarlist[i], opts->file2);

            if (!opts->warn[NCCMP_W_ALL])
                status = EXIT_DIFFER;
            if (opts->force)
                continue;
            else
                break;
        }

        if (vars1[idx1].len != vars2[idx2].len) {
            fprintf(stdout, "DIFFER : SIZE OF VARIABLE \"%s%s\" : %d <> %d\n", getGroupPath(), opts->cmpvarlist[i], (int) vars1[idx1].len, (int) vars2[idx2].len);

            if (!opts->warn[NCCMP_W_ALL])
                status = EXIT_DIFFER;
            if (opts->force)
                continue;
            else
                break;
        }

        if (vars1[idx1].type != vars2[idx2].type) {
            type2string(vars1[idx1].type, str1);
            type2string(vars2[idx2].type, str2);
            fprintf(stdout, "DIFFER : TYPE OF VARIABLE \"%s%s\" : %s <> %s\n", getGroupPath(), opts->cmpvarlist[i], str1, str2);

            if (!opts->warn[NCCMP_W_ALL])
                status = EXIT_DIFFER;

            if (opts->force)
                continue;
            else
                break;
        }

        /* Has rec? */
        if (vars1[idx1].hasrec) {

            /* Compare only if # recs are equal and not zero. */
            if ((nrec1 == nrec2) && (nrec1 + nrec2)) {
                /* TODO: Check if ignorem missing. */
                cmpstatus = nccmpdatarecvar(ncid1, ncid2, opts->cmpvarlist[i], opts, 0, nrec1 - 1, user_types1, user_types2);
                if (cmpstatus) {
                    status = cmpstatus;
                    if (opts->force)
                        continue;
                    else
                        break;
                }
            }
        } else {
            /* TODO: Check if ignorem missing. */
            cmpstatus = cmpvar(opts->cmpvarlist[i], -1, opts, ncid1, ncid2, user_types1, user_types2);
            if (cmpstatus) {
                status = EXIT_FAILURE;
                if (opts->force)
                    continue;
                else
                    break;
            }
        }
    }

    if (opts->verbose)
        printf("INFO: Finished comparing data.\n");
    
    freestringlist(&processed, nprocessed);
    return status;
}

/* *********************************************************** */
int
nccmp(nccmpopts* opts) {
    int ncid1, ncid2;
    int status;
    int nuser_types1, nuser_types2;
    nccmp_user_type_t *user_types1, *user_types2;
    if (opts->verbose)
        printf("INFO: Opening input files.\n");

    status = openfiles(opts, &ncid1, &ncid2);
    if (status)
        return status;

    if (opts->verbose)
        printf("INFO: Creating variable comparison list.\n");

    if (opts->verbose)
        printf("INFO: Comparing file formats.\n");

    status += nccmpformats(opts, ncid1, ncid2);
    if (status && !opts->force)
        return status;

    if (opts->verbose)
        printf("INFO: Comparing global attributes.\n");
    status += nccmpglobalatts(opts, ncid1, ncid2);

    if (status && !opts->force)
        return status;

    if (opts->verbose)
        printf("INFO: Collecting dimension information for first file.\n");

    getdiminfo(ncid1, dims1, &ndims1);

    if (opts->verbose)
        printf("INFO: Collecting dimension information for second file.\n");

    getdiminfo(ncid2, dims2, &ndims2);
    status += compareGroup(opts, ncid1, ncid2);
    if (status && !opts->force)
        return status;

    // Check root group
    status += makecmpvarlist(opts, ncid1, ncid2);
    if (status && !opts->force)
        return status;

    if (opts->verbose)
        printf("INFO: Collecting variable information for first file.\n");

    user_types1 = getvarinfo(ncid1, vars1, &nvars1, opts->verbose, &nuser_types1);

    if (opts->verbose)
        printf("INFO: Collecting variable information for second file.\n");

    user_types2 = getvarinfo(ncid2, vars2, &nvars2, opts->verbose, &nuser_types2);

    if (nuser_types1 != nuser_types2)
        status++;

    if (status && !opts->force)
        return status;

    status += nccmprecinfo(opts, ncid1, ncid2);
    if (status && !opts->force)
        return status;

    if (opts->metadata) {
        status += nccmpmetadata(opts, ncid1, ncid2);
    }

    if (status && !opts->force)
        return status;

    if (opts->data) {
        status += nccmpdata(opts, ncid1, ncid2, user_types1, user_types2);
    }
    
    freevarinfo(nuser_types1, user_types1);
    freevarinfo(nuser_types2, user_types2);
    
    if (status && !opts->force)
        return status;

    if (opts->verbose)
        printf("INFO: Comparisons complete. Freeing memory.\n");

    return status;
}

int compareGroup(nccmpopts* opts, int ncid1, int ncid2) {

    int numgrps1, numgrps2;
    int i, status=0;
    GROUP_NODE *groups1 = NULL;
    GROUP_NODE *groups2 = NULL;
    int nuser_types1, nuser_types2;
    nccmp_user_type_t *user_types1, *user_types2;

    vector<string> groupNames1;
    vector<string> groupNames2;
    char name[NC_MAX_NAME];
    int *gids;
    
    // find the common list of group names
    nc_inq_grps(ncid1, &numgrps1, NULL);
    if (numgrps1 > 0) {
        gids = (int *) malloc(sizeof (int) * numgrps1);
        nc_inq_grps(ncid1, NULL, gids);

        for (i = 0; i < numgrps1; i++) {
            nc_inq_grpname(gids[i], name);
            groupNames1.push_back(name);
        }
        free(gids);
    }
    nc_inq_grps(ncid2, &numgrps2, NULL);
    if (numgrps2 > 0) {
        gids = (int *) malloc(sizeof (int) * numgrps2);
        nc_inq_grps(ncid2, NULL, gids);

        for (i = 0; i < numgrps2; i++) {
            nc_inq_grpname(gids[i], name);
            groupNames2.push_back(name);
        }
        free(gids);
    }

    // make sure names1 only has group names in it that are also in names2
    vector<string>::iterator it;
    for(i=0; i < (int)groupNames1.size(); i++) {
        it = find(groupNames2.begin(), groupNames2.end(), groupNames1[i]);
        if(it == groupNames2.end()) {
            printf("DIFFER : GROUP : %s%s : Does not exist in file 2\n", getGroupPath(), groupNames1[i].c_str());
            status++;
            if (status && !opts->force) {
                free(groups1);
                free(groups2);
                return status;
            }
            groupNames1.erase(groupNames1.begin() + i);
            i--;
        }
    }

    // check to see if that are groups in file2 that are not in file1
    for(i=0; i < (int)groupNames2.size(); i++) {
        it = find(groupNames1.begin(), groupNames1.end(), groupNames2[i]);
        if(it == groupNames1.end()) {
            printf("DIFFER : GROUP : %s%s : Does not exist in file 1\n", getGroupPath(), groupNames2[i].c_str());
            status++;
            if (status && !opts->force) {
                free(groups1);
                free(groups2);
                return status;
            }
        }
    }
    
    numgrps1 = groupNames1.size();
    
    if(numgrps1 > 0) {
        groups1 = (GROUP_NODE *) malloc(sizeof (GROUP_NODE) * numgrps1);
        getgroupinfo(ncid1, groupNames1, groups1);
        groups2 = (GROUP_NODE *) malloc(sizeof (GROUP_NODE) * numgrps1);
        getgroupinfo(ncid2, groupNames1, groups2);

        for (i = 0; i < numgrps1; i++) {
            opts->extentcount = 0;
            groupPath.push_back(groupNames1[i]);
            if (opts->verbose)
                printf("INFO: Comparing group %s:\n", getGroupPath());
            status += makecmpvarlist(opts, groups1[i].groupID, groups2[i].groupID);
            if (status && !opts->force) {
                free(groups1);
                free(groups2);
                return status;
            } else {
                // Recursively call groups within groups
                status += compareGroup(opts, groups1[i].groupID, groups2[i].groupID);
            }

            if (opts->verbose)
                printf("INFO: Comparing group attributes [%s]:\n", getGroupPath());
            status += nccmpglobalatts(opts, groups1[i].groupID, groups2[i].groupID);
            if (status && !opts->force) {
                free(groups1);
                free(groups2);
                return status;
            }
            
            if (opts->verbose)
                printf("INFO: Collecting variable information for first file.\n");

            user_types1 = getvarinfo(groups1[i].groupID, vars1, &nvars1, opts->verbose, &nuser_types1);

            if (opts->verbose)
                printf("INFO: Collecting variable information for second file.\n");

            user_types2 = getvarinfo(groups2[i].groupID, vars2, &nvars2, opts->verbose, &nuser_types2);
            status += (nuser_types1 != nuser_types2);
            if (status && !opts->force) {
                free(groups1);
                free(groups2);
                freevarinfo(nuser_types1, user_types1);
                freevarinfo(nuser_types2, user_types2);
                return (status);
            }


            status += nccmprecinfo(opts, groups1[i].groupID, groups2[i].groupID);
            if (status && !opts->force) {
                free(groups1);
                free(groups2);
                freevarinfo(nuser_types1, user_types1);
                freevarinfo(nuser_types2, user_types2);
                return status;
            }

            if (opts->metadata) {
                status += nccmpmetadata(opts, groups1[i].groupID, groups2[i].groupID);
            }

            if (status && !opts->force) {
                free(groups1);
                free(groups2);
                freevarinfo(nuser_types1, user_types1);
                freevarinfo(nuser_types2, user_types2);
                return status;
            }

            if (opts->data) {
                status += nccmpdata(opts, groups1[i].groupID, groups2[i].groupID, user_types1, user_types2);
            }
            
            freevarinfo(nuser_types1, user_types1);
            freevarinfo(nuser_types2, user_types2);
            
            if (status && !opts->force) {
                free(groups1);
                free(groups2);
                return status;
            }
            
            clearstringlist(opts->cmpvarlist, opts->ncmpvarlist);
            groupPath.pop_back();
        }
    }
    free(groups1);
    free(groups2);
    return status;
}

/* *********************************************************** */

int
main(int argc, char** argv) {
    /* Change chunk cache. */
    nc_set_chunk_cache(1024 * 1024 * 128, 2003, 0.75f); // 128MB
    
    int status;
    nccmpopts opts;

    status = EXIT_SUCCESS;

    initnccmpopts(&opts);

    /* parse command-line args. & options */
    status = getnccmpopts(argc, argv, &opts);

    if (status != EXIT_SUCCESS)
        goto end;

    if (opts.verbose)
        printf("INFO: Command-line options parsed.\n");

    status = nccmp(&opts);

end:
    if (opts.report_identical) {
        if (!(opts.help || opts.version) &&
                (status == EXIT_SUCCESS)) {
            printf("Files \"%s\" and \"%s\" are identical.\n",
                    opts.file1, opts.file2);
        }
    }

    freenccmpopts(&opts);

    exit(status);
}
