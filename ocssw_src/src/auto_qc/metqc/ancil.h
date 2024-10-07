/* ancil.h 
 * Header file for SeaWiFS ancillary related HDF programs
 */

#ifndef ANCIL_H_
#define ANCIL_H_

/* 
 * shared includes
 */

#include <stdint.h>

/*
 * Common flags
 */

#define DBG       1
#define ERROR     1
#define FAILURE   5
#define ISFIRST   1
#define NOTFIRST  0
#define LAST      999 
#define SUCCESS   0
#define OK        0

/*
 * HDF settings for ancillary data products
 * (see individual programs for other specific settings)
 */

#define LISTSIZE     1000
#define MAXLABLEN    80
#define MAXDESCLEN   1000
#define MAXNAMELNG   255
#define GEOMNAME     "Equal-Angle SDS"
#define GEOMCLASS    "Geometry"
#define GEOMSIZE     ((2*sizeof(int32_t))+(6*sizeof(float)))
#define AVGVGRPNAME  "Average"         /* inner HDF Average Vgroup name   */
#define STDVGRPNAME  "Std. Deviation"  /* inner HDF Std. Dev. Vgroup name */
#define OBSVGRPNAME  "No. of Obs."     /* inner HDF No. Obs. Vgroup name  */

/*
 * Flags for data bin origins (per Doug Ilg 9/21/93)
 */

#define UNKNOWN   0
#define NW        1
#define NORTH     2
#define NE        3
#define WEST      4
#define CENTER    5
#define EAST      6
#define SW        7
#define SOUTH     8
#define SE        9

/*
 * HDF metadata annotation struct 
 */

struct annotation {
    char label[MAXLABLEN];
    char descr[MAXDESCLEN];
    int32_t type;
};

extern struct annotation *annot;

#endif /*ANCIL_H_ */
