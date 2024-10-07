#ifndef L1_SEABASS_H_
#define L1_SEABASS_H_

#include "filehandle.h"
#include "l1.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct seabass {
    FILE *fp;
    long int data_start;
    const char *delimiter;
    int *field_indexes;
    int lon_index;
    int lat_index;
    int year_index;
    int month_index;
    int day_index;
    int hour_index;
    int minute_index;
    int second_index;
    int current_row;
} seabass;

int open_seabass(filehandle *file);
int read_seabass(filehandle *file, l1str *l1rec);
int close_seabass(filehandle *file);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
