/*
 * mph_flags.h
 *
 *  Created on: Sep 8, 2015
 *      Author: rhealy
 */

#ifndef SRC_L2GEN_MPH_FLAGS_H_
#define SRC_L2GEN_MPH_FLAGS_H_

#define NMPHFLAGS    4
#define MPH_FLOAT    0x0001       /* floating aquatic vegetation.  */
#define MPH_CYANO    0x0002       /* Cyanobacteria              .  */
#define MPH_ADJ      0x0004       /* Adjacency effect (stray light)*/
#define MPH_BADINPUT 0x0008       /* Bad Input */

#define NHABSFLAGS     6
#define HABS_CLOUD     0x0001       /* Cloud Flag */
#define HABS_NONWTR    0x0002       /* Not Water */
#define HABS_SNOWICE   0x0004       /* Snow/Ice */
#define HABS_ADJ       0x0008       /* Adjacency effect (stray light)*/
#define HABS_BADINPUT  0x0010       /* Bad Input */
#define HABS_NONCYANO  0x0020       /* Elevated non-cyano CI */

static const char * const mph_flag_lname[NMPHFLAGS] = {"FLOAT",
    "CYANO", "ADJACENCY", "BADINPUT"};
static const char * const habs_flag_lname[NHABSFLAGS] = {"CLOUD",
    "NONWATER", "SNOWICE", "ADJACENCY", "BADINPUT", "NONCYANO"};

#endif /* SRC_L2GEN_MPH_FLAGS_H_ */
