/*************************************************************/
/*               Sentinel 2A MSI L1C PROCESSING              */
/*                    Jiaying He 06/09/2016                  */
/*************************************************************/

#ifndef L1_MSI_H
#define L1_MSI_H

#include "l1.h"

/* -------------------------------------------------------------------------- */
/* l1c_msi Header */
/* Define */

#ifdef __cplusplus
extern "C" {
#endif
    /* Definition of all l1c_msi.c functions */

    int openl1_msi(filehandle *file);
    int readl1_msi(filehandle *file, int recnum, l1str *l1rec, int lonlat);
    int closel1_msi(filehandle *file);
    int readl1_msi_lonlat(filehandle *file, int recnum, l1str *l1rec);

#ifdef __cplusplus
}
#endif



#endif