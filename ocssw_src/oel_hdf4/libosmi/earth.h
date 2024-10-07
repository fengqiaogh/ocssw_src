/*
 *----------------------------------------------------------------------
 * @(#) earth.h	1.0	30 Mar 93	<shc>
 *  Copyright (c) 1993, CSIRO Division of Oceanography.
 *----------------------------------------------------------------------
 *
 * Definitions of some basic earth-related constants.
 *
 * EMK		Gravitational constant (m**3 / s**2)
 * J2-J5	First four zonal gravitational harmonics
 * FLT		Flattening
 * FLTINV	1 / FLT
 * EQRAD	Equatorial radius (m)
 * POLRAD	Polar radius (m)
 * AEMEAN	Mean radius (m)
 * EAROT	Rotation rate in radians/sec
 * TEARTH	Minutes per revolution
 * ECCSQ	Eccentricity squared
 * NSPD		Number of seconds per (24 hour) day
 * NMPD		Number of minutes per (24 hour) day
 */

#ifndef __EARTH__
#define __EARTH__

#define EMK (3.986004415e14)
#define J2 (1.0862583e-3)
#define J3 (-2.5338975e-6)
#define J4 (-1.623821e-6)
#define J5 (-2.2963180e-7)
#define FLT (3.35281317789691443e-3)
#define FLTINV (298.257)
#define EQRAD (6378136.3)
#define POLRAD (EQRAD * (1 - FLT))
#define AEMEAN ((EQRAD + POLRAD)/2)
#define EAROT (7.2921159e-5)
#define TEARTH (1436.068176e0)
#define ECCSQ (FLT * (2 - FLT))
#define NSECPD (86400)
#define NMINPD (1440)

#endif /* __EARTH__ */
