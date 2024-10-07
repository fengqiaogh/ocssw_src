/*
 *  pml.c
 *
 *	Compute IOPs from rrs using PML IOP algorithm
 *
 *  Plymouth Marine Laboratory
 *  UK
 */

/*
 *  This code implements the IOP code of Smyth et al. "" Applied Optics
 */

#include <stdlib.h>
#include <math.h>
#include "pml.h"
#include "l12_proto.h"

static int iterate = 1;
static int idx410 = -1;
static int idx440 = -1;
static int idx490 = -1;
static int idx510 = -1;
static int idx555 = -1;
static int idx670 = -1;
static int initialized = 0;

static float *aw;
static float *bbw;

void pml_assert(int expr) {
    if (expr) return;
    else exit(1);
}

/**
 *  pml_is_initialized - 
 *  DESCRIPTION:
 *  Return 1 if pml_init has been previously called; 0, otherwise
 */

int pml_is_initialized(void) {
    return initialized;
}

/**
 *  pml_init - PML IOP algorithm
 *  @iter8: iterate (0-no,1-yes)
 *  @i410: 0-relative index in spectrum of 410 nm
 *  @i440: 0-relative index in spectrum of 440 nm
 *  @i490: 0-relative index in spectrum of 490 nm
 *  @i510: 0-relative index in spectrum of 510 nm
 *  @i555: 0-relative index in spectrum of 555 nm
 *  @i670: 0-relative index in spectrum of 670 nm
 *  DESCRIPTION:
 *  This routine should be called to initialize various parameters that may
 *  be adjusted in the PML algorithm or that must be known apriori.  
 */

int pml_init(int iter8, int i410, int i440, int i490, int i510, int i555, int i670, float *awptr, float *bbwptr) {
    iterate = iter8;
    idx410 = i410;
    idx440 = i440;
    idx490 = i490;
    idx510 = i510;
    idx555 = i555;
    idx670 = i670;
    aw = awptr;
    bbw = bbwptr;
    initialized = 1;
    return 0;
}

