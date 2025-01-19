/**
 * @file l1bgen_oci.h
 * @brief Header file for OCI (Ocean Color Instrument) L1B data generation.
 *
 * This file contains definitions and structures used in the L1B data generation
 * process for the Ocean Color Instrument. It includes enums for device types
 * and gain dimensions, structures for calibration lookup tables, and utility
 * functions for environment variable expansion.
 *
 * @note This file uses C++11 features when available.
 *
 * @authors Jakob C Lindo (SSAI), Joel Gales (SAIC)
 * @date Aug 2024
 */

#ifndef _L1BGEN_OCI_H_
#define _L1BGEN_OCI_H_

#include <stdint.h>
#include <netcdf>
#include <genutils.h>
#include "geolocate_oci.h"
#include "corrections.hpp"
#include "calibrations.hpp"
#include "l1b_file.hpp"

#define CHUNK_CACHE_SIZE 256 * 1024 * 1024  // 256MiB of cache memory.
#define CHUNK_CACHE_NELEMS 1033
#define CHUNK_CACHE_PREEMPTION .75
#define VARCHUNK_CACHE_SIZE 4 * 1024 * 1024  // 4Mib of cache memory.
#define VARCHUNK_CACHE_NELEMS 1033
#define VARCHUNK_CACHE_PREEMPTION .75
#define CHUNKBANDS 40
#define CHUNKPIXELS 2000
#define CHUNKLINES 16

#endif  // _L1BGEN_OCI_H_
