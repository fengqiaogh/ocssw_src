#ifndef _SCENE_META_H
#define _SCENE_META_H

#include "l1.h"
#include <dfutils.h>

#define DAYSCENE       0
#define NIGHTSCENE     1
#define DAYANDNIGHT    2
#define UNKNOWNSCENE   3

#define ASCENDING      0
#define DSCENDING      1
#define UNKNOWNNODE    2

#ifdef __cplusplus
extern "C" {
#endif

typedef struct scene_meta_str {
    char start_node[32];
    char end_node[32];
    char daynight[32];
    double start_time;
    double center_time;
    double end_time;
    int start_year;
    int start_day;
    int start_msec;
    int end_year;
    int end_day;
    int end_msec;
    float earth_sun_dist_corr;
    float start_center_lon;
    float start_center_lat;
    float scene_center_lon;
    float scene_center_lat;
    float scene_center_solz;
    float end_center_lon;
    float end_center_lat;
    float northern_lat;
    float southern_lat;
    float eastern_lon;
    float western_lon;
    float upperleft_lon;
    float upperright_lon;
    float upperleft_lat;
    float upperright_lat;
    float lowerleft_lon;
    float lowerright_lon;
    float lowerleft_lat;
    float lowerright_lat;

} scnstr;

void scene_meta_put(l1str *l1rec);
scnstr *scene_meta_get(void);
void scene_meta_write(idDS ds_id);

#ifdef __cplusplus
}
#endif

#endif
