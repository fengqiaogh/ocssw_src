/* =========================================================== */
/* Module l1_io.c                                              */
/*                                                             */
/* Functions to open, close, read, and write  a level-1b file, */
/* with the format determined by the file handle content.      */
/*                                                             */
/* Written By:                                                 */
/*                                                             */
/*     B. A. Franz                                             */
/*     SAIC General Sciences Corp.                             */
/*     NASA/SIMBIOS Project                                    */
/*     April 1998                                              */
/*                                                             */
/* Modifications By:                                           */
/*    J. Gales                                                 */
/*    Futuretech                                               */
/*    NASA/SIMBIOS Project                                     */
/*    10/00                                                    */
/*                                                             */
/*    Add support for OCTSL1A                                  */
/*    W. Robinson, SAIC, 10 Dec 2004  add CZCS, VIIRS support  */
/* =========================================================== */
#include "l1.h"

#include "scene_meta.h"

#include <stdio.h>

#include "l1_mos.h"
#include "l1_hdf_generic_read.h"
#include "l1_seawifs.h"
#include "l1_seawifs_netcdf.h"
#include "l1_octs.h"
#include "l1_osmi.h"
#include "l1_modis.h"
#include "l1_czcs.h"
#include "l1_xcal.h"
#include "l1_aci.h"
#include "l1_sgli.h"
#include "l1_ocm.h"
#include "l1_ocm2.h"
#include "l1_ocmdb.h"
#include "l1_meris_N1.h"
#include "l1_meris_CC.h"
#include "l1_safe.h"
#include "l1_viirs_h5.h"
#include "l1_viirs_l1b.h"
#include "l1_hico.h"
#include "l1_goci.h"
#include "l1_oli.h"
#include "l1_viirs_l1a.h"
#include "l1_ocia.h"
#include "l1_aviris.h"
#include "l1_prism.h"
#include "l1_nc_generic_read.h"
#include "l1_l5tm.h"
#include "l1_l7etm.h"
#include "l1_msi.h"
#include "l1_hawkeye.h"
#include "l1_misr.h"
#include "l1_oci.h"
#include "l1_ocis.h"
#include "l1_seabass.h"
#include "l1_l1c.h"
#include "l1_l1c_anc.h"
#include "l1_spexone.h"

/* ---------------------------------------------------------------- */
/* Close the level 1 file associated with the input file handle.    */

/* ---------------------------------------------------------------- */
void closel1(filehandle *l1file) {
    switch (l1file->format) {
    case FT_L1HDF:
        if (l1file->mode == READ)
            closel1_hdf_g(l1file);
        else
            closel1_generic(l1file);
        break;
    case FT_L1BNCDF:
        closel1_generic(l1file);
        break;
    case FT_MOSL1B:
        closel1_mos(l1file);
        break;
    case FT_SEAWIFSL1A:
        closel1_seawifs(l1file);
        break;
    case FT_SEAWIFSL1ANC:
        closel1_seawifs_netcdf(l1file);
        break;
    case FT_OCTSL1B:
        closel1_octs(l1file);
        break;
    case FT_OCTSL1A:
        closel1_octs(l1file);
        break;
    case FT_OSMIL1A:
        closel1_osmi(l1file);
        break;
    case FT_L1XCAL:
        closel1_xcal(l1file);
        break;
    case FT_CZCSL1A:
        closel1_czcs(l1file);
        break;
    case FT_MODISL1B:
        closel1_modis();
        break;
    case FT_CLASSAVHRR:
        closel1_aci(l1file);
        break;
    case FT_OCML1B:
        closel1_ocm(l1file);
        break;
    case FT_OCM2L1B:
        closel1_ocm2(l1file);
        break;
    case FT_OCML1BDB:
        closel1_ocmdb(l1file);
        break;
    case FT_MERISL1B:
        closel1_meris_N1(l1file);
        break;
    case FT_MERISCC:
        closel1_meris_CC(l1file);
        break;
    case FT_MERISL1BSAFE:
        closel1_safe(l1file);
        break;
    case FT_VIIRSL1B:
        closel1_viirs_h5(l1file);
        break;
    case FT_VIIRSL1BNC:
        closel1_viirs_l1b();
        break;
    case FT_VIIRSL1A:
        closel1_viirs_l1a(l1file);
        break;
    case FT_HICOL1B:
        closel1_hico(l1file);
        break;
    case FT_GOCIL1B:
        closel1_goci(l1file);
        break;
    case FT_OLIL1B:
        closel1_oli(l1file);
        break;
    case FT_OCIA:
        closel1_ocia(l1file);
        break;
    case FT_OCIL1B:
        closel1_oci(l1file);
        break;
    case FT_OCIS:
        closel1_ocis(l1file);
        break;
    case FT_AVIRIS: {
        aviris_t* data = (aviris_t*) l1file->private_data;
        if (data->isnetcdf)
            closel1_aviris_nc(l1file);
        else
            closel1_aviris(l1file);
        } break;
    case FT_PRISM:
        closel1_prism(l1file);
        break;
    case FT_OLCI:
        closel1_safe(l1file);
        break;
    case FT_SGLI:
        closel1_sgli(l1file);
        break;
    case FT_L5TML1B:
        closel1_l5tm(l1file);
        break;
    case FT_L7ETML1B:
        closel1_l7etm(l1file);
        break;
    case FT_MSIL1C:
        closel1_msi(l1file);
        break;
    case FT_HAWKEYEL1A:
        closel1_hawkeye(l1file);
        break;    
    case FT_MISR:
        closel1_misr(l1file);
        break;
    case FT_SEABASSRRS:
        close_seabass(l1file);
        break;
    case FT_L1C:
        closel1_l1c(l1file);
        break;
    case FT_L1CANC:
        closel1_l1c_anc(l1file);
        break;
    case FT_SPEXONE:
        closel1_spexone(l1file);
        break;
    default:
        fprintf(stderr,
                "-E- %s Line %d: l1close - Unknown L1 file format specifier: %d\n",
                __FILE__, __LINE__, l1file->format);
        break;
    };

    free(l1file->Fonom);
    l1file->Fonom = NULL;

    return;
}

/* ---------------------------------------------------------------- */
/* Open a level 1 file for reading, load file handle with metadata  */
/* from file header. At a minimum, the name and format fields of    */
/* file handle must be loaded before this function is called.       */

/* ---------------------------------------------------------------- */
int openl1(filehandle *l1file) {
    int status = 1;
    /*
    if (l1file->sensorID == SEABASS){
      return open_seabass(l1file);
    }
    */

    /* Get number of bands and band indexing from sensor table */
    l1file->nbands = rdsensorinfo(l1file->sensorID, l1_input->evalmask, NULL, NULL);
    if (l1file->nbands < 0) {
        printf("-E- %s line %d: Error reading sensor table\n",
                __FILE__, __LINE__);
        return (status);
    }
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "Bindx", (void **) &l1file->bindx);
    l1file->nbandsir = rdsensorinfo(l1file->sensorID, l1_input->evalmask, "NbandsIR", NULL);

    /* set wavelength index */
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "Lambda", (void **) &l1file->iwave);
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "fwave", (void **) &l1file->fwave);
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "Fobar", (void **) &l1file->Fobar);
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "Tau_r", (void **) &l1file->Tau_r);
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "k_oz", (void **) &l1file->k_oz);
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "k_no2", (void **) &l1file->k_no2);
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "aw", (void **) &l1file->aw);
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "bbw", (void **) &l1file->bbw);

    bindex_set(l1file->iwave, l1file->nbands + l1file->nbandsir, BANDW);

    if ((l1file->Fonom = (float*) calloc(l1file->nbands, sizeof (float))) == NULL) {
        printf("-E- (%s, %d) Cannot allocate space for l1file->Fonom\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    int i;
    for (i = 0; i < l1file->nbands; i++) {
        if (l1_input->outband_opt >= 2) {
            get_f0_thuillier_ext(l1file->iwave[i], BANDW, l1file->Fonom + i);
        } else {
            l1file->Fonom[i] = l1file->Fobar[i];
        }
    }

    /* Open L1 file for reading or writing, as requested */
    if (l1file->mode == READ) {

        switch (l1file->format) {
        case FT_MOSL1B:
            status = openl1_read_mos(l1file);
            break;
        case FT_SEAWIFSL1A:
            status = openl1_seawifs(l1file);
            break;
        case FT_SEAWIFSL1ANC:
            status = openl1_seawifs_netcdf(l1file);
            break;
        case FT_L1HDF:
            status = openl1_read_hdf_g(l1file);
            break;
        case FT_L1BNCDF:
            status = openl1_nc_generic(l1file);
            break;
        case FT_OCTSL1B:
            status = openl1_octs(l1file);
            break;
        case FT_OCTSL1A:
            status = openl1_octs(l1file);
            break;
        case FT_OSMIL1A:
            status = openl1_osmi(l1file);
            break;
        case FT_L1XCAL:
            status = openl1_xcal(l1file);
            break;
        case FT_CZCSL1A:
            status = openl1_czcs(l1file);
            break;
        case FT_MODISL1B:
            status = openl1_modis(l1file);
            break;
        case FT_CLASSAVHRR:
            status = openl1_aci(l1file);
            break;
        case FT_OCML1B:
            status = openl1_ocm(l1file);
            break;
        case FT_OCM2L1B:
            status = openl1_ocm2(l1file);
            break;
        case FT_OCML1BDB:
            status = openl1_ocmdb(l1file);
            break;
        case FT_MERISL1B:
            status = openl1_meris_N1(l1file);
            break;
        case FT_VIIRSL1B:
            status = openl1_viirs_h5(l1file);
            break;
        case FT_VIIRSL1BNC:
            status = openl1_viirs_l1b(l1file);
            break;
        case FT_VIIRSL1A:
            status = openl1_viirs_l1a(l1file);
            break;
        case FT_HICOL1B:
            status = openl1_hico(l1file);
            break;
        case FT_GOCIL1B:
            status = openl1_goci(l1file);
            break;
        case FT_MERISCC:
            status = openl1_meris_CC(l1file);
            break;
        case FT_MERISL1BSAFE:
            status = openl1_safe(l1file);
            break;
        case FT_OLIL1B:
            status = openl1_oli(l1file);
            break;
        case FT_OCIA:
            status = openl1_ocia(l1file);
            break;
        case FT_OCIL1B:
            status = openl1_oci(l1file);
            break;
        case FT_OCIS:
            status = openl1_ocis(l1file);
            break;
        case FT_AVIRIS: {
            aviris_t* data = (aviris_t*)l1file->private_data;
            if (data->isnetcdf)
                status = openl1_aviris_nc(l1file);
            else
                status = openl1_aviris(l1file);
            } break;
        case FT_PRISM:
            status = openl1_prism(l1file);
            break;
        case FT_OLCI:
            status = openl1_safe(l1file);
            break;
        case FT_SGLI:
            status = openl1_sgli(l1file);
            break;
        case FT_L5TML1B:
            status = openl1_l5tm(l1file);
            break;
        case FT_L7ETML1B:
            status = openl1_l7etm(l1file);
            break;
        case FT_MSIL1C:
            status = openl1_msi(l1file);
            break;
        case FT_HAWKEYEL1A:
            status = openl1_hawkeye(l1file);
            break;
        case FT_MISR:
          status = openl1_misr(l1file);
          break;
        case FT_SEABASSRRS:
            status = open_seabass(l1file);
            break;
        case FT_L1C:
            status = openl1_l1c(l1file);
            break;
        case FT_SPEXONE:
            status = openl1_spexone(l1file);
            break;
        case FT_L1CANC:
            status = openl1_l1c_anc(l1file);
            break;
        default:
        printf("openl1 - Unknown L1 input file format specifier: %d\n",
                l1file->format);
        break;
        };

    } else {

        switch (l1file->format) {
        case FT_L1HDF:
        case FT_L1BNCDF:
            status = openl1_write(l1file);
            break;
        default:
            printf("Unknown L1 output file format specifier: %d\n",
                    l1file->format);
            break;
        };
    }

    return (status);
}


/* ---------------------------------------------------------------- */
/* Read a specific level 1 record from the file pointed to by the   */
/* input file handle and load the data into the l1 record structure.*/

/* ---------------------------------------------------------------- */
int readl1(filehandle *l1file, int32_t recnum, l1str *l1rec) {
    int status=1;
    int32_t ip;
    aviris_t* data = (aviris_t*) l1file->private_data;

    /* Clear the L1 record */
    init_l1(l1rec);
    l1rec->tilt = 0.0;
    l1rec->mside = 0;
    l1rec->detnum = 0;
    /* Altitude of sensor should be set in l1_sensor.c
     * If not set it will be assumed sensor is above atmosphere
     * and no rayleigh correction will be done.
     */
    l1rec->alt = BAD_FLT;

    for (ip = 0; ip < l1file->npix; ip++) {
        l1rec->pixnum[ip] = ip;
        l1rec->slot[ip] = 0;
        l1rec->alpha[ip] = 0.0;
    }

    l1rec->iscan = recnum;

    l1rec->l1file = l1file;

    switch (l1file->format) {
    case FT_L1HDF:
        status = readl1_hdf_g(l1file, recnum, l1rec);
        break;
    case FT_L1BNCDF:
        status = readl1_nc_generic(l1file, recnum, l1rec);
        break;
    case FT_MOSL1B:
        status = readl1_mos(l1file, recnum, l1rec);
        break;
    case FT_SEAWIFSL1A:
        status = readl1_seawifs(l1file, recnum, l1rec);
        break;
    case FT_SEAWIFSL1ANC:
        status = readl1_seawifs_netcdf(l1file, recnum, l1rec);
        break;
    case FT_OCTSL1B:
        status = readl1_octs(l1file, recnum, l1rec);
        break;
    case FT_OCTSL1A:
        status = readl1_octs(l1file, recnum, l1rec);
        break;
    case FT_OSMIL1A:
        status = readl1_osmi(l1file, recnum, l1rec);
        break;
    case FT_L1XCAL:
        status = readl1_xcal(l1file, recnum, l1rec);
        break;
    case FT_CZCSL1A:
        status = readl1_czcs(l1file, recnum, l1rec);
        break; 
    case FT_MODISL1B:
        status = readl1_modis(l1file, recnum, l1rec, 0);
        break;
    case FT_CLASSAVHRR:
        status = readl1_aci(l1file, recnum, l1rec);
        break;
    case FT_OCML1B:
        status = readl1_ocm(l1file, recnum, l1rec);
        break;
    case FT_OCM2L1B:
        status = readl1_ocm2(l1file, recnum, l1rec);
        break;
    case FT_OCML1BDB:
        status = readl1_ocmdb(l1file, recnum, l1rec);
        break;
    case FT_MERISL1B:
        status = readl1_meris_N1(l1file, recnum, l1rec);
        break;
    case FT_MERISCC:
        status = readl1_meris_CC(l1file, recnum, l1rec);
        break;
    case FT_MERISL1BSAFE:
        status = readl1_safe(l1file, recnum, l1rec);
        break;
    case FT_VIIRSL1B:
        status = readl1_viirs_h5(l1file, recnum, l1rec, 0);
        break;
    case FT_VIIRSL1BNC:
        status = readl1_viirs_l1b(l1file, recnum, l1rec, 0);
        break;
    case FT_VIIRSL1A:
        status = readl1_viirs_l1a(l1file, recnum, l1rec);
        break;
    case FT_HICOL1B:
        status = readl1_hico(l1file, recnum, l1rec, 0);
        break;
    case FT_GOCIL1B:
        status = readl1_goci(l1file, recnum, l1rec, 0);
        break;
    case FT_OLIL1B:
        status = readl1_oli(l1file, recnum, l1rec, 0);
        break;
    case FT_OCIA:
        status = readl1_ocia(l1file, recnum, l1rec);
        break;
    case FT_OCIL1B:
        status = readl1_oci(l1file, recnum, l1rec, 0);
        break;
    case FT_OCIS:
        status = readl1_ocis(l1file, recnum, l1rec);
        break;
    case FT_AVIRIS:
        if (data->isnetcdf)
            status = readl1_aviris_nc(l1file, recnum, l1rec);
        else
            status = readl1_aviris(l1file, recnum, l1rec);
        //    status = 0;
        break;
    case FT_PRISM:
        status = readl1_prism(l1file, recnum, l1rec);
        //    status = 0;
        break;
    case FT_OLCI:
        status = readl1_safe(l1file, recnum, l1rec);
        //    status = 0;
        break;
    case FT_SGLI:
        status = readl1_sgli(l1file, recnum, l1rec);
        //    status = 0;
        break;
    case FT_L5TML1B:
        status = readl1_l5tm(l1file,recnum,l1rec, 0);
        break;
    case FT_L7ETML1B:
        status = readl1_l7etm(l1file,recnum,l1rec, 0);
        break;
    case FT_MSIL1C:
         status = readl1_msi(l1file,recnum,l1rec, 0);
         break;
    case FT_HAWKEYEL1A:
         status = readl1_hawkeye(l1file,recnum,l1rec);
         break;
    case FT_MISR:
         status = readl1_misr(l1file,l1rec);
         break;
    case FT_SEABASSRRS:
         status = read_seabass(l1file,l1rec);
         break;
    default:
        printf("readl1 - Unknown L1 input file format specifier: %d\n",
               l1file->format);
        break;
    };


    if (status != 0) {
        fprintf(stderr,
                "-E- %s Line %d: Error reading L1B.\n",
                __FILE__, __LINE__);
        return (1);
    }

    for (ip = 0; ip < l1file->npix; ip++) {
        if (l1rec->lon[ip] <= -180.0 || l1rec->lon[ip] >= 180.0 || isnan(l1rec->lon[ip]) ||
                l1rec->lat[ip] <= -90.0 || l1rec->lat[ip] >= 90.0 || isnan(l1rec->lat[ip])) {
            l1rec->navfail[ip] = 1;
        }
    }

    /* Reduce to sub-sample scan if requested */
    if (l1subpix(l1file, l1rec) != 0)
        return (1);

    setflagbits_l1(0, l1rec, -1);

    /* update scene meta */
    scene_meta_put(l1rec);

    return (0);
}




/* ---------------------------------------------------------------- */
/* Read lon and lat for a specific level 1 record from the file     */
/* pointed to by the input file handle and load the data into the   */
/* l1 record structure.                                             */
/*                                                                  */
/* Lonlat needs to read only the following variables                */
/*      Lon,                                                        */
/*      Lat,                                                        */
/*      Solar Zentih,                                               */
/*      Time                                                        */
/* ---------------------------------------------------------------- */
int readl1_lonlat(filehandle *l1file, int32_t recnum, l1str *l1rec) {
    int status;
    l1rec->iscan = recnum;

    for (int ip = 0; ip < l1file->npix; ip++) {
        l1rec->navfail[ip] = 0;
        l1rec->flags[ip] = 0;
    }

    switch (l1file->format) {
    case FT_SEAWIFSL1A:
        status = readl1_lonlat_seawifs(l1file, recnum, l1rec);
        break;
    case FT_SEAWIFSL1ANC:
        status = readl1_lonlat_seawifs_netcdf(l1file, recnum, l1rec);
        break;
    case FT_MODISL1B:
        status = readl1_modis(l1file, recnum, l1rec, 1);
        break;
    case FT_MERISL1B:
        status = readl1_lonlat_meris_N1(l1file, recnum, l1rec);
        break;
    case FT_VIIRSL1B:
        status = readl1_viirs_h5(l1file, recnum, l1rec, 1);
        break;
    case FT_VIIRSL1BNC:
        status = readl1_viirs_l1b(l1file, recnum, l1rec, 1);
        break;
    case FT_HICOL1B:
        status = readl1_hico(l1file, recnum, l1rec, 1);
        break;
    case FT_GOCIL1B:
        status = readl1_goci(l1file, recnum, l1rec, 1);
        break;
    case FT_OLIL1B:
        status = readl1_oli(l1file, recnum, l1rec, 1);
        break;
    case FT_L5TML1B:
        status = readl1_l5tm(l1file,recnum,l1rec, 1);
        break;
    case FT_L7ETML1B:
        status = readl1_l7etm(l1file,recnum,l1rec, 1);
        break;
    case FT_MSIL1C:
         status = readl1_msi(l1file, recnum, l1rec, 1);
         break;
    case FT_OCIL1B:
        status = readl1_oci(l1file, recnum, l1rec, 1);
        break;
    case FT_L1C:
         status = readl1_l1c(l1file, recnum, l1rec);
         break;
    case FT_L1CANC:
         status = readl1_l1c_anc(l1file, recnum, l1rec);
         break;
    case FT_SPEXONE:
        status = readl1_spexone(l1file, recnum, l1rec);
        break;
    default:
        status = readl1(l1file, recnum, l1rec);
        break;
    };

    if (status != 0) {
        fprintf(stderr,
                "-E- %s Line %d: Error reading L1B.\n",
                __FILE__, __LINE__);
        return (1);
    }

    for (int ip = 0; ip < l1file->npix; ip++) {
        if (l1rec->lon[ip] <= -180.0 || l1rec->lon[ip] >= 180.0 || isnan(l1rec->lon[ip]) ||
                l1rec->lat[ip] <= -90.0 || l1rec->lat[ip] >= 90.0 || isnan(l1rec->lat[ip])) {
            l1rec->navfail[ip] = 1;
            l1rec->flags[ip] |= NAVFAIL;
        }
    }

    return (0);
}

