/* -------------------------------------------------------------------- */
/* printindx - prints formatted report based on L0 info structure       */
/*                                                                      */
/* ---------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

#include <timeutils.h>

#include "swl0_parms.h"
#include "swl0_types.h"
#include "swl0_struc.h"
#include "swl0_utils.h"

char *format_time(double usec) {
    static char tempstr[21];
    time_t temp = (time_t) usec;
    strftime(tempstr, 21, "%b %d %Y %T", gmtime((time_t *) & temp));
    return (tempstr);
}

void printindx(swl0indx *indx) {
    printf("\nSeaWiFS Level-0 Format Report\n");
    printf("File Name:               %s\n", indx->l0file);

    if (indx->type == HRPT)
        printf("Type:                    %s\n", "HRPT");
    else
        printf("Type:                    %s\n", "GAC/LAC");

    /*
    printf("AOS Time:                %24s",
        asctime(gmtime((time_t *)&(indx->timeAOS))));
    printf("LOS Time:                %24s",
        asctime(gmtime((time_t *)&(indx->timeLOS))));

    printf("Start_Time:              %s\n",ydhmsf((double) indx->timeFirst,'G'));
    printf("End_Time:                %s\n",ydhmsf((double) indx->timeLast ,'G'));
     */

    printf("AOS Time:                %s\n", format_time(indx->timeAOS));
    printf("LOS Time:                %s\n", format_time(indx->timeLOS));

    printf("Start_Time:              %s\n", format_time(indx->timeFirst));
    printf("End_Time:                %s\n", format_time(indx->timeLast));

    printf("Start_Zulu:              %s\n", ydhmsf((double) indx->timeFirst, 'G'));

    printf("Number of Frames:        %d\n", indx->nrecs);
    printf("GAC Frames:              %d\n", indx->ngac);
    printf("LAC Frames:              %d\n", indx->nlac);
    printf("Lunar Cal Frames:        %d\n", indx->nlun);
    printf("Solar Cal Frames:        %d\n", indx->nsol);
    printf("Intergain Cal Frames:    %d\n", indx->nigc);
    printf("TDI Cal Frames:          %d\n", indx->ntdi);
    printf("Frame Timing Errors:     %d\n", indx->timeError);
    printf("SCID Errors:             %d\n", indx->scidError);
    printf("SOH Header Errors:       %d\n", indx->sohError);
    printf("Reported Bit Error Rate: %e\n", indx->bitRateRep);
    printf("Computed Bit Error Rate: %e\n", indx->bitRateComp);

    return;
}
