/* 
 * File:   setupflags.h
 * Author: dshea
 *
 * Created on January 31, 2012, 8:47 AM
 */

#ifndef SETUPFLAGS_H
#define SETUPFLAGS_H

#include <stdint.h>


#ifdef __cplusplus
extern "C" {
#endif

void setupflags(char *flagsdef, char *flaguse, uint32_t *flagusemask,
        uint32_t *required, int *status, int32_t * BITS);


#ifdef __cplusplus
}
#endif

#endif /* SETUPFLAGS_H */
