/* -------------------------------------------------------------------- */
/* printscene - prints formatted report based on L0 scene structure     */
/*                                                                      */
/* ---------------------------------------------------------------------*/

#include <stdio.h>
#include <time.h>

#include "swl0_parms.h"
#include "swl0_types.h"
#include "swl0_struc.h"
#include "swl0_utils.h"

void printscene(int nscenes, swl0scene scene[]) {
    int is;
    time_t stime;

    printf("\nSeaWiFS Level-0 Scene Index Report\n");
    printf("Number of scenes found: %d\n", nscenes);
    for (is = 0; is < nscenes; is++) {

        stime = (time_t) scene[is].stime;

        printf("Scene: %3d, Type: %2d, Start: %4d, Frames: %4d, Time: %s",
                is, scene[is].mnftype,
                scene[is].srec, scene[is].nrec,
                asctime(gmtime(&stime)));
    }

    return;
}




