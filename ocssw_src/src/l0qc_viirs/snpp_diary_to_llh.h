#define MAX_RECORDS 10000
#include <stdint.h>
int convert_diary(int npkts, uint8_t dstore[], double otime[], double orb[], double atime[], double quat[]);
int orb2lla(int nRecords, double orb[], double lon[], double lat[], double alt[]);