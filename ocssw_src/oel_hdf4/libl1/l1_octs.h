#ifndef  _L1_OCTS_H
#define  _L1_OCTS_H

#include <stdint.h>
#include "l1.h"

#define MAXOCLIN 6700      /* max # lines */
#define MAXOCPIX 2218      /* max # pixels */
#define MAXOCARR 10000
#define NOCBANDS 8

int openl1_octs(filehandle *l1file);
int readl1_octs(filehandle *l1file, int32_t recnum, l1str *l1rec);
int closel1_octs(filehandle *l1file);

int32_t get_octs_cal(char *file, int16_t year, int16_t day, int32_t msec[MAXOCLIN],
        int32_t recnum, int16_t npix, int32_t spix, int32_t tilt,
        int16_t gainset[MAXOCLIN], float inst_temp[MAXOCLIN],
        int16_t sample_table[3][8][2][400][2], int32_t scansPerScene,
        uint16_t l1acnts[NOCBANDS][MAXOCPIX], float l1brads[MAXOCPIX][NOCBANDS]);

int CalcViewAngle(float lon1, float lat1, float pos[3],float usun[3]);
int LeapCheck(int yr);
void reform_octs_time(char *time);
int navigation(int32_t fileID);

#endif

