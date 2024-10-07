#ifndef _SWL0_PROTO_H
#define _SWL0_PROTO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include <unistd.h>

#include "swl0_parms.h"
#include "swl0_types.h"
#include "swl0_struc.h"
#include "swl0_utils.h"
#include "swl1_struc.h"
#include "swl1_hdf.h"
#include "input_s.h"
#include "navblk_s.h"
#include "navctl_s.h"
#include "navqc_s.h"
#include "orbit_s.h"
#include "tilt_s.h"
#include "genutils.h"


INT32 getl0indx(char *filename, swl0ctl *l0ctl, swl0indx *indx);
INT16 getl0scene(swl0indx *indx, swl0scene *scene);
INT32 getnavdata(swl0scene *scene, input_sType navinp[]);
INT32 getorbdata(swl0indx *indx, input_sType gps[]);
void printindx(swl0indx *indx);
void printnav(swl0scene *scene);
int getl0scene_nav(FLOAT32 xnodel, INT32 tnode,
        navblk_sType navblk[], tilt_states_sType *tiltblk, swl0scene *scene);
void printscene(int nscenes, swl0scene *scene);
INT32 getl1rec(INT16 iframe, swl0scene *scene, swl0ctl *l0ctl,
        input_sType navinp[], navblk_sType navblk[],
        tilt_states_sType *tiltblk, swl1rec l1rec[]);
void addL1Metrics(INT32 scanNum, swl1rec *l1rec);
l1met *getL1Metrics(void);
int mkmeta(const char *metafile, const char *l1afile, swl0scene *scene, swl0ctl *l0ctl);
INT32 getorbnum(FLOAT64 usec);
INT16 valid_instlm(INT16 mnftype, INT16 mnfnum, INT16 scanNum);
INT16 getEngQual(FLOAT32 ins_ana[], BYTE eng_qual[]);

INT32 locate_temporal_anomalies(swl0indx *indx, char *timefile);
timeTab *temporal_anomaly(INT16 mnftype, INT16 year, INT16 day, INT32 msec);


extern int conv_soh_(BYTE *soh, FLOAT32 *scana, BYTE *scdis);
extern int conv_ins_(BYTE *soh, FLOAT32 *insta, BYTE *instd);
extern int acs_block_(void);
extern int ins_block_(void);

extern int initnav_(input_sType input[],
        INT32 *nframes,
        navctl_sType *navctl,
        navqc_sType *navqc,
        orbit_sType *orbit,
        INT32 *status);

extern int swfnav_(navqc_sType *navqc,
        navctl_sType *navctl,
        input_sType input[],
        INT32 *nframes,
        orbit_sType *orbit,
        INT32 *nlines,
        navblk_sType navblk[],
        tilt_states_sType *tiltblk,
        FLOAT32 *xnodel,
        INT32 *tnode);

extern int geonav_(FLOAT32 pos[3],
        FLOAT32 rm[3][3],
        FLOAT32 coef[6],
        FLOAT32 sunref[3],
        INT32 *spix,
        INT32 *ipix,
        INT32 *npix,
        FLOAT32 lat[],
        FLOAT32 lon[],
        FLOAT32 solz[],
        FLOAT32 sola[],
        FLOAT32 senz[],
        FLOAT32 sena[]);

int get_swf_def_meta(char *seadas_vs, char *data_center);


int get_swf_def_meta1(char *seadas_vs,
        char *data_center,
        char *station_name,
        float *station_lat,
        float *station_lon,
        char *station_code);

int gen_soft_id(char *seadas_vs, char *prog_vs, char *soft_id);

#endif
