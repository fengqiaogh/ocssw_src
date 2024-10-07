/* -------------------------------------------------------------------- */
/* printnav - prints formatted report based on scene navigation info    */
/*                                                                      */
/* ---------------------------------------------------------------------*/

#include <stdio.h>
#include <time.h>
#include "swl0_proto.h"

void printnav(swl0scene *scene) {
    int i;
    time_t utime;

    if (scene->nscan <= 0) {
        printf("Scene has no valid scan lines.\n");
        return;
    }

    printf("\nScene Navigation Summary Report\n");
    printf("Lines in scene:   %d\n", scene->nscan);
    printf("Data Type:        %d\n", scene->mnftype);
    utime = (time_t) scene->stime;
    printf("Start Time:       %s", asctime(gmtime(&utime)));
    utime = (time_t) scene->etime;
    printf("End   Time:       %s", asctime(gmtime(&utime)));
    printf("Start Node:       %s\n", scene->start_node);
    printf("End   Node:       %s\n", scene->end_node);
    printf("Node Longitude:   %0.2f\n", scene->node_lon);
    utime = (time_t) scene->node_time;
    printf("Node Time:        %s", asctime(gmtime(&utime)));
    printf("Orbit Number:     %d\n", scene->orbnum);

    printf("start center lon = %10.2f\n", scene->start_center_lon);
    printf("start center lat = %10.2f\n", scene->start_center_lat);
    printf("end center lon   = %10.2f\n", scene->end_center_lon);
    printf("end center lat   = %10.2f\n", scene->end_center_lat);
    printf("northern lat     = %10.2f\n", scene->northern_lat);
    printf("southern lat     = %10.2f\n", scene->southern_lat);
    printf("western lon      = %10.2f\n", scene->western_lon);
    printf("eastern lon      = %10.2f\n", scene->eastern_lon);

    printf("\nNumber of tilt states in scene: %5d\n", scene->ntilts);
    for (i = 0; i < scene->ntilts; i++) {
        printf("Tilt %2d -\n", i);
        printf("       Start Scan = %10d\n", scene->tilt_ranges[i][0]);
        printf("       End   Scan = %10d\n", scene->tilt_ranges[i][1]);
        printf("       Tilt  Flag = %10d\n", scene->tilt_flags[i]);
    }
    printf("\n\n");

    return;
}




