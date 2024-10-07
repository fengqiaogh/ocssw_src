/*
 *  pml.h
 *
 *  Compute IOPs from RRS using PML IOP algorithm
 *
 *  Plymouth Marine Laboratory
 *  UK
 */

#ifndef _PML_H
#define _PML_H

int pml_is_initialized(void);
int pml_init(int iter8, int i410, int i440, int i490, int i510, int i555, int i640, float *aw, float *bbw);

#endif
