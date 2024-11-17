#include "nav.h"
#include "sun2000.h"
void ymdhms2jul_(int32_t *year, int32_t *month, int32_t *day, int32_t *hour, int32_t *minute, double *sec,
                 double *jul) {
    ymdhms2jul(*year, *month, *day, *hour, *minute, *sec, jul);
}

void jddate_(int *jd, int *i, int *j, int *k) {
    jddate(*jd, i, j, k);
}

void jdate_(int *jd, int *i, int *k) {
    jdate(*jd, i, k);
}

void julian_(double tin[2], double *tout) {
    *tout = julian(tin);
}

int32_t jd_(int32_t *i, int32_t *j, int32_t *k) {
    return jd(*i, *j, *k);
}

void gha2000_(int *iyr, double *day, double *gha) {
    gha2000(*iyr, *day, gha);
}
