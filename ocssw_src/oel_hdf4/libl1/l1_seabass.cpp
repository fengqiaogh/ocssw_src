/* =========================================================== */
/* Module l1_seabass.c                                         */
/*                                                             */
/* Functions to open, close, read, and write  a level-1b file, */
/* with the format determined by the file handle content.      */
/*                                                             */
/* Written By:                                                 */
/*    J. Gales                                                 */
/*    Futuretech                                               */
/*    NASA/SIMBIOS Project                                     */
/*    02/17                                                    */
/*                                                             */
/* =========================================================== */

#include <ctype.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <regex.h>

#include <filetype.h>

#include "l1_seabass.h"
#include "l1.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

int open_seabass(filehandle *l1file){
    l1file->private_data = malloc(sizeof(seabass));
    seabass *priv = (seabass*)(l1file->private_data);
    priv->delimiter = ",";
    priv->lat_index = -1;
    priv->lon_index = -1;
    priv->year_index = -1;
    priv->month_index = -1;
    priv->day_index = -1;
    priv->hour_index = -1;
    priv->minute_index = -1;
    priv->second_index = -1;

    std::vector<int32_t>  sensor_lambdas;
    for(int32_t i=0; i<l1file->nbands + l1file->nbandsir; i++) {
        sensor_lambdas.push_back(l1file->iwave[i]);
    }

    int status = 0;

    l1file->npix = 1;
    priv->field_indexes = (int*)calloc(l1file->nbands + l1file->nbandsir, sizeof(int));

    char buffer[2048];
    buffer[2047] = '\0';  // make sure the buffer is string terminated
    priv->fp = fopen(l1file->name, "r");
    while (fgets(buffer, 2047, priv->fp)){
        if (strncmp(buffer, "/delimiter=", 11) == 0) {
            char *delim = &buffer[11];
            if (!strncmp(delim, "comma", 5)){
                priv->delimiter = ",";
            } else if (!strncmp(delim, "space", 5)){
                priv->delimiter = " ";
            } else if (!strncmp(delim, "tab", 3)){
                priv->delimiter = "\t";
            } else {
                status = 1;
            }
        } else if (strncmp(buffer, "/fields=", 8) == 0) {
            char *token = strtok(&buffer[8], ",");
            regex_t regx;
            regcomp(&regx, "^rrs_?([0-9][^_\n]*)\n?$", REG_EXTENDED | REG_ICASE);
            regmatch_t matches[2];

            std::vector<double> bands_found{};
            std::vector<int> band_indexes{};
            int field_i = 0;
            while (token != NULL) {
                if (!strcmp(token, "lat")){
                    priv->lat_index = field_i;
                } else if (!strcmp(token, "lon")){
                    priv->lon_index = field_i;
                } else if (!strcmp(token, "year")){
                    priv->year_index = field_i;
                } else if (!strcmp(token, "month")){
                    priv->month_index = field_i;
                } else if (!strcmp(token, "day")){
                    priv->day_index = field_i;
                } else if (!strcmp(token, "hour")){
                    priv->hour_index = field_i;
                } else if (!strcmp(token, "minute")){
                    priv->minute_index = field_i;
                } else if (!strcmp(token, "second")){
                    priv->second_index = field_i;
                } else if (regexec(&regx, token, 2, matches, 0) == 0) {
                    char value[matches[1].rm_eo - matches[1].rm_so + 1];
                    memcpy(value, &token[matches[1].rm_so], matches[1].rm_eo - matches[1].rm_so);
                    value[matches[1].rm_eo - matches[1].rm_so] = 0;
                    bands_found.push_back(std::stod(value));
                    band_indexes.push_back(field_i);
                }
                token = strtok(NULL, ",");
                ++field_i;
            } // field loop

            int nband_i = 0;
            for (auto l=sensor_lambdas.cbegin(); l != sensor_lambdas.cend(); l++, nband_i++){
                int closest_i = -1;
                double closest_distance = std::numeric_limits<double>::max();
                for (int field_i=0;field_i<(int)bands_found.size();field_i++){
                    double distance = std::abs(*l - bands_found[field_i]);
                    if (distance < closest_distance){
                        closest_i = field_i;
                        closest_distance = distance;
                    }
                }
                if (closest_distance <= 5){
                    priv->field_indexes[nband_i] = band_indexes[closest_i];
                } else {
                    priv->field_indexes[nband_i] = -1;
                }
            }
        } else if (strncmp(buffer, "/end_header", 11) == 0){
            break;
        }
    }

    priv->data_start = ftell(priv->fp);

    l1file->nscan = 0;
    while (fgets(buffer, 2047, priv->fp)) {
        if ((int32_t)strlen(buffer) > l1file->nbands){
            l1file->nscan++;
        }
    }

    fseek(priv->fp, priv->data_start, SEEK_SET);
    priv->current_row = 0;

    return status;
}

int read_seabass(filehandle *file, l1str *l1rec) {
    char buffer[2048];

    seabass *priv = (seabass*)(file->private_data);
    l1rec->is_l2 = true;
    l1rec->npix = 1;

    if (l1rec->iscan < priv->current_row){
        fseek(priv->fp, priv->data_start, SEEK_SET);
        priv->current_row = 0;
    }
    while (l1rec->iscan < priv->current_row){
        if (!fgets(buffer, 2047, priv->fp)){
            priv->current_row++;
            break;
        }
    }

    for (int i=0;i<file->nbands;i++) {
        l1rec->Lt[i] = BAD_FLT;
    }
    l1rec->lat[0] = BAD_FLT;
    l1rec->lon[0] = BAD_FLT;

    int year = 1970;
    int month = 1;
    int day= 1;
    int hour = 0;
    int minute = 0;
    double second = 0;

    if (fgets(buffer, 2047, priv->fp)) {
        priv->current_row++;

        int32_t field_i = 0;
        char value[128];
        char *token = strtok(buffer, priv->delimiter);

        while (token != NULL){
            strcpy(value, token);

            if (field_i == priv->lat_index) {
                l1rec->lat[0] = (float)atof(value);
            } else if (field_i == priv->lon_index) {
                l1rec->lon[0] = (float)atof(value);
            } else if (field_i == priv->year_index) {
                year = atoi(value);
            } else if (field_i == priv->month_index) {
                month = atoi(value);
            } else if (field_i == priv->day_index) {
                day = atoi(value);
            } else if (field_i == priv->hour_index) {
                hour = atoi(value);
            } else if (field_i == priv->minute_index) {
                minute = atoi(value);
            } else if (field_i == priv->second_index) {
                second = atoi(value);
            } else {
                for (int i=0;i<file->nbands;i++) {
                    if (priv->field_indexes[i] == field_i) {
                        l1rec->Lt[i] = (float) atof(value);
                        break;
                    }
                }
            }

            token = strtok(NULL, priv->delimiter);
            field_i++;
        }
    }

    // set the scan time
    double secs = second + 60 * (minute + 60 * hour);
    l1rec->scantime = ymds2unix(year, month, day, secs);

    return 0;
}

int close_seabass(filehandle *file){
    seabass *priv = (seabass*)(file->private_data);
    fclose(priv->fp);
    free(priv->field_indexes);
    free(file->private_data);

    return 0;
}
