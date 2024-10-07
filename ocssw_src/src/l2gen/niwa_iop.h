/*
 * niwa_iop.h -- interface to NIWA / UoP / Moore IOP algorithm
 *               (Gerald Moore, Sam Lavender, Matt Pinkerton)
 *
 * Algorithm citation: Pinkerton et al., 2006. 'A method for estimating
 *      inherent optical properties of New Zealand continental shelf waters
 *      from satellite ocean colour measurements', New Zealand Journal of
 *      Marine and Freshwater Research 40, pp 227-247.
 */

#ifndef __NIWA_IOP_H__
#define __NIWA_IOP_H__

#include "l2_struc.h"


/* 
 * interface functions 
 */

/* 
 * niwa_iop()
 *  Calculate the IOP values for one scan line
 *
 *
 *  inputs:
 *      l2rec - level-2 structure containing one complete scan
 *              after atmospheric correction
 *  outputs:
 *      niwa_a    - total absorption coefficient - water absorption
 *      niwa_bb   - total backscatter coefficient - water backscatter
 *      niwa_iopf - per pixel error flags (0 if algorithm succeeded, refer flags_iop.h) 
 */
extern void niwa_iop(l2str *l2rec, float niwa_a[], float niwa_bb[], int16 niwa_iopf[]);


#endif  /* __NIWA_IOP_H__ */
