#include <stdint.h>
#include "l12_parms.h"
#include "passthebuck.h"
#include "readL2scan.h"


#define  SHORTNAME   128

typedef char stname[SHORTNAME];
typedef char prname[PRODSTRLEN];

typedef struct mscal_struct {
    int32_t sensorID;
    int32_t nfiles;
    int32_t npixs;
    int32_t nbands;
    int32_t *Lambda;
    int32_t nprods;
    prname *l2prods;
    char *input_parms;
    unsigned char *data;
    stname *filenames;
    int16 *fileID;
    int16 *year;
    int16 *day;
    int32 *msec;
    int16 *iscan;
    uint8 *mside;
    uint8 *detnum;
    int16 *pixnum;
    float *lon;
    float *lat;
    float32 *ddata;
    char oformat[20];
} mscalstr;

typedef struct calinput_struct {
    char ifile [MAXNFILES][FILENAME_MAX];
    char ofile [FILENAME_MAX];
    char input_parms[20000];
    char flaguse [1024];
    int32_t spixl; /* starting pixel no. of the input (1-rel)  */
    int32_t epixl; /* ending pixel no. of the input (1-rel)    */
    int32_t dpixl; /* pixel subsampling increment              */
    int32_t sline; /* starting line no. of the input (1-rel)   */
    int32_t eline; /* ending line no. of the input (1-rel)     */
    int32_t dline; /* line subsampling increment               */
    char oformat[20];
} inputstr;





int crosscal_append(char *crosscalfile, mscalstr calstr);
int crosscal_read(char *crosscalfile, int32_t subsmpl, mscalstr *calstr);
int crosscal_npixs(char *xcalfile, int32_t subsmpl, int32_t *sensorID, int32_t *npixs, int32_t *ngranuls);
int32_t alloc_calstr(int32_t nfiles, int32_t npixs, mscalstr *calstr);
void free_calstr(mscalstr calstr, int all);
