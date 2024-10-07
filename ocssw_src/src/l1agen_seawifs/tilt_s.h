#ifndef TILT_S_H_
#define TILT_S_H_

#include "swl0_parms.h"

typedef struct tilt_states_struct {
    float tilt[MAXFRAMES]; /* Tilts per line in scene   */
    int ntilts; /* Number of tilts in scene  */
    short int tilt_flags[20]; /* Tilt flag ??              */
    short int tilt_ranges[20][2]; /* Start and end line number */
} tilt_states_sType;

#endif /* TILT_S_H_ */
