#ifndef _SMI_CLIMATOLOGY_H
#define _SMI_CLIMATOLOGY_H

#define NPROD    3
#define ALPHA510 0
#define TAUA865  1
#define SST      2

int smi_climatology_init(char *file, int day, int prodID);
float smi_climatology(float lon, float lat, int prodID);


#endif
