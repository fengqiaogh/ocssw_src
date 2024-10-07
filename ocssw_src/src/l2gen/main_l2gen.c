/* =====================================================================*/
/*                                                                      */
/* Program: MSl12 - multi-sensor level-1b to level-2 conversion         */
/*                                                                      */
/* Usage:                                                               */
/*     See msl12_usage().                                               */
/*                                                                      */
/* Written By:                                                          */
/*                                                                      */
/*     Bryan A. Franz                                                   */
/*     SAIC General Sciences Corp.                                      */
/*     NASA/SIMBIOS Project                                             */
/*     March 1998                                                       */
/*                                                                      */
/* Modified By:                                                         */
/*                                                                      */
/*     Joel M. Gales                                                    */
/*     Futuretech                                                       */
/*     NASA/SIMBIOS Project                                             */
/*     September 1999                                                   */
/*     Add support for SeaWifs specific L2 metadata                     */
/*                                                                      */
/* =====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <libgen.h>
#include "l12_proto.h"
#include "lonlat2pixline.h"
#include "version.h"
#include "geo_region.h"
/* -------------------------------------------------------------------- */
/*                            main                                      */

/* -------------------------------------------------------------------- */
int main(int argc, char* argv[]) {
    int32_t iscan = 0; /* input scan number                  */
    int32_t oscan = 0; /* output scan number                 */
    int32_t npix = 0; /* input number pixels per scan       */
    int32_t spix = 0; /* start pixel for subscene process   */
    int32_t epix = -1; /* end pixel for subscene process     */
    int32_t dpix = 1; /* pixel increment for sub-sampling   */
    int32_t sscan = 0; /* start scan for subscene process    */
    int32_t escan = -1; /* end scan for subscene process      */
    int32_t dscan = 1; /* scan subsampling increment         */

    l1str *l1rec; /* generic level-1b scan structure    */
    l2str *l2rec; /* generic level-2  scan structure    */
    tgstr *tgrec; /* structure to store target values   */
    aestr *aerec; /* structure to store aerosol values  */

    filehandle l1file; /* input l1 file handle               */
    filehandle tgfile; /* input target file handle           */
    filehandle aefile; /* input aerosol file handle          */
    filehandle ofile[MAX_OFILES]; /* output file handles                */

    double start_time;
    int num_ofiles = 0;
    int32_t i;

    if (argc == 1) {
        l2gen_usage("l2gen");
        return 0;
    }

    for (i = 0; i < argc; i++) {
        // see if help on command line
        if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-help") == 0)) {
            l2gen_usage("l2gen");
            return 0;
        }

        // see if prodxmlfile is on command line
        if (strncmp(argv[i], "prodxmlfile=", 12) == 0) {
            char fname[FILENAME_MAX];

            parse_file_name(argv[i] + 12, fname);
            init_l2prod();
            printf("Writing product information to XML file %s\n", fname);
            write_product_XML_file(fname);
            return 0;
        }
    }

    setvbuf(stdout, NULL, _IOLBF, 0);
    setvbuf(stderr, NULL, _IOLBF, 0);

    // allocate structures
    l1rec = (l1str*) malloc(sizeof (l1str));
    l2rec = (l2str*) malloc(sizeof (l2str));
    tgrec = (tgstr*) malloc(sizeof (tgstr));
    aerec = (aestr*) malloc(sizeof (aestr));
    msl12_input_init();

    if (!l1rec || !l2rec || !tgrec || !aerec || !input) {
        printf("-E- %s %d: Error allocating data structures.\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    /* Initialize file handles */
    cdata_();
    filehandle_init(&l1file);
    filehandle_init(&tgfile);
    filehandle_init(&aefile);
    for (i = 0; i < MAX_OFILES; i++)
        filehandle_init(&ofile[i]);

    /* Parse input parameters */
    if (msl12_input(argc, argv, "msl12", &l1file) != 0) {
        printf("-E- %s: Error parsing input parameters.\n", argv[0]);
        exit(FATAL_ERROR);
    }

    if (access(input->ifile[0], F_OK) || access(input->ifile[0], R_OK)) {
        printf("-E- %s: Input file '%s' does not exist or cannot open.\n",
                argv[0], input->ifile[0]);
        exit(FATAL_ERROR);
    }

    /* if north, south, east, and west used, convert to pix and line */
    if (input->north != -999 || input->south != -999
            || input->east != -999 || input->west != -999
            || input->xbox != -1 || input->ybox != -1) {

        int result;
        int save_verbose = want_verbose;

        // turn off the extra output
        want_verbose = 0;
        if (input->north != -999 && input->south != -999
                && input->east != -999 && input->west != -999
                && input->xbox == -1 && input->ybox == -1) {

            // north, south, east and west are set
            result = lonlat2pixline1(input->ifile[0], input->geofile,
                    l1_input->resolution, input->west, input->south, input->east, input->north,
                    &l1_input->spixl, &l1_input->epixl, &l1_input->sline, &l1_input->eline);

        } else if (input->north == -999 && input->south != -999
                && input->east == -999 && input->west != -999
                && input->xbox != -1 && input->ybox != -1) {

            // south, west, xbox, ybox are set
            result = lonlat2pixline2(input->ifile[0], input->geofile,
                    l1_input->resolution, input->west, input->south, input->xbox, input->ybox,
                    &l1_input->spixl, &l1_input->epixl, &l1_input->sline, &l1_input->eline);

        } else {
            printf("-E- %s: set either \n", argv[0]);
            printf("    1) north, south, east and west or\n");
            printf("    2) south, west, xbox, ybox\n");
            exit(FATAL_ERROR);
        }
        want_verbose = save_verbose;

        if (result == 120) { // requested box includes the whole file
            l1_input->spixl = 1;
            l1_input->epixl = -1;
            l1_input->sline = 1;
            l1_input->eline = -1;
        } else {
            if (result != 0 && result != 110) {
                printf("-E- %s: north, south, east, west box not in this file.\n", argv[0]);
                exit(LONLAT_ERROR);
            }
        }
    }

    /*                                                                  */
    /* Determine number of output files.                                */
    /*                                                                  */
    for (i = 0, num_ofiles = 0; i < MAX_OFILES; i++, num_ofiles++)
        if (input->ofile[i][0] == '\0')
            break;

    if (num_ofiles == 0) {
        printf("-E- %s: No output file name given\n", argv[0]);
        exit(FATAL_ERROR);
    }

    /**
     * check if the georegion has been supplied and check if it exists and can be opened.
     * 
     */
    if(input->georegionfile[0]) {
        if (access(input->georegionfile, F_OK) || access(input->georegionfile, R_OK)) {
            printf("-E- %s: Input file '%s' does not exist or cannot open.\n",
                    argv[0], input->georegionfile);
            exit(FATAL_ERROR);
        }
    }

    /*                                                                  */
    /* Open input file and get sensor and scan information from handle. */
    /*                                                                  */
    if (openl1(&l1file) != 0) {
        printf("-E- %s: Error opening %s for reading.\n",
                argv[0], l1file.name);
        exit(FATAL_ERROR);
    }

    npix = l1file.npix;

    // Open aerosol file if provided

    if (input->aerfile[0] != '\0') {
        aerec->mode = ON;
        strcpy(aefile.name, input->aerfile);
        if (open_aer(&aefile) != 0) {
            printf("-E- %s: Error opening %s for reading.\n",
                    argv[0], aefile.name);
            exit(FATAL_ERROR);
        }
        if (aefile.npix != npix) {
            printf("-E- %s: Incompatible scan length between %s and %s: %d, %d.\n",
                    argv[0], l1file.name, aefile.name, aefile.npix, npix);
            exit(FATAL_ERROR);
        }
    } else {
        aerec->mode = OFF;
    }

    /*                                                                  */
    /* Allocate memory for L1 and L2 scan data and opional input recs   */
    /*                                                                  */
    if (alloc_l1(&l1file, l1rec) == 0) {
        printf("-E- %s: Unable to allocate L1 record.\n", argv[0]);
        exit(FATAL_ERROR);
    }
    if (alloc_l2(l1rec, l2rec) == 0) {
        printf("-E- %s: Unable to allocate L2 record.\n", argv[0]);
        exit(FATAL_ERROR);
    }    // Open calibration target file if provided

    if (input->tgtfile[0] != '\0') {
        tgrec->mode = ON;
        strcpy(tgfile.name, input->tgtfile);
        file_format format = getFormat(tgfile.name);
        tgfile.format = format.type;
        tgfile.sensorID = format.sensor_id;
        tgfile.subsensorID = format.subsensor_id;
        tgrec->sensorID = format.sensor_id;
        tgrec->subsensorID = format.subsensor_id;

        if (tgfile.format != FT_L3BIN) {
            if (open_target(&tgfile) != 0) {
                printf("-E- %s: Error opening %s for reading.\n",
                        argv[0], tgfile.name);
                exit(FATAL_ERROR);
            }
            if (tgfile.npix != npix) {
                printf("-E- %s: Incompatible scan length between %s and %s.\n",
                        argv[0], l1file.name, tgfile.name);
                exit(FATAL_ERROR);
            }
        }
        l2rec->tgrec = tgrec;
    } else {
        tgrec->mode = OFF;
        l2rec->tgrec = NULL;
    }


    if (aerec->mode == ON) {
        if (alloc_aer(npix, l1file.nbands, aerec) == 0) {
            printf("-E- %s: Unable to allocate aerfile record.\n", argv[0]);
            exit(FATAL_ERROR);
        }
    }
    if (tgrec->mode == ON) {
        if (alloc_target(npix, l1file.nbands, tgrec) == 0) {
            printf("-E- %s: Unable to allocate target record.\n", argv[0]);
            exit(FATAL_ERROR);
        }
    }

    /* Set the end pixel if it was not set by command argument          */
    if (l1_input->epixl == -1 || l1_input->epixl > l1file.npix)
        l1_input->epixl = l1file.npix;
    if (l1_input->eline == -1 || l1_input->eline > l1file.nscan)
        l1_input->eline = l1file.nscan;
    if (l1_input->spixl < 1)
        l1_input->spixl = 1;
    if (l1_input->sline < 1)
        l1_input->sline = 1;

    spix = MAX(l1_input->spixl - 1, 0);
    epix = MIN(l1_input->epixl - 1, l1file.npix - 1);
    dpix = MAX(l1_input->dpixl, 1);
    sscan = MAX(l1_input->sline - 1, 0);
    escan = MIN(l1_input->eline - 1, l1file.nscan - 1);
    dscan = MAX(l1_input->dline, 1);

    if (sscan > escan || spix > epix) {
        printf("-E- %s: scan and pixel limits make no sense.\n", argv[0]);
        printf(" start scan  = %d\n", sscan + 1);
        printf(" end   scan  = %d\n", escan + 1);
        printf(" start pixel = %d\n", spix + 1);
        printf(" end   pixel = %d\n", epix + 1);
        exit(FATAL_ERROR);
    }

    /* Note: for the L1 file, npix is still the native scan pixel count */
    l1file.spix = spix; /* start pixel rel to L1 scan    */
    l1file.epix = epix; /* end   pixel rel to L1 scan    */

    /*                                                                  */
    /* Open output file(s)                                              */
    /*                                                                  */
    if (input->mode != FORWARD) {

        /*                                                              */
        /* Transfer sensor and scan info to recal filehandle and open   */
        /*                                                              */
        strcpy(ofile[0].name, input->ofile[0]);
        ofile[0].format = FT_L1HDF;
        ofile[0].mode = WRITE;
        ofile[0].sensorID = l1file.sensorID;
        ofile[0].subsensorID = l1file.subsensorID;
        ofile[0].nbands = l1file.nbands;
        ofile[0].nbandsir = l1file.nbandsir;
        ofile[0].nlvl = l1file.nlvl;
        ofile[0].bindx = l1file.bindx;
        ofile[0].ndets = l1file.ndets;
        ofile[0].spix = spix;
        ofile[0].epix = epix;
        ofile[0].npix = (epix - spix) / dpix + 1;
        ofile[0].length = l1rec->length;
        ofile[0].nscan = (escan - sscan) / dscan + 1;

        printf("Opening L1B output file: %s\n", ofile[0].name);

        if (openl1(&ofile[0]) != 0) {
            printf("-E- %s: Error opening %s for writing.\n",
                    argv[0], ofile[0].name);
            exit(FATAL_ERROR);
        }

    } else {

        /*                                                              */
        /* Transfer sensor and scan info to output filehandles and open */
        /*                                                              */
        printf("\n");
        for (i = 0; i < num_ofiles; i++) {

            strcpy(ofile[i].name, input->ofile[i]);
            if (strcmp(input->oformat, "netCDF4") == 0)
                ofile[i].format = FT_L2NCDF;
            else
                ofile[i].format = FT_L2HDF;
            ofile[i].mode = WRITE;
            ofile[i].sensorID = l1file.sensorID;
            ofile[i].nbands = l1file.nbands;
            ofile[i].nbandsir = l1file.nbandsir;
            ofile[i].nlvl = l1file.nlvl;
            ofile[i].bindx = l1file.bindx;
            ofile[i].ndets = l1file.ndets;
            ofile[i].subsensorID = l1file.subsensorID;
            strcpy(ofile[i].spatialResolution, l1file.spatialResolution);
            ofile[i].spix = spix;
            ofile[i].epix = epix;
            ofile[i].npix = (epix - spix) / dpix + 1;
            ofile[i].length = l2rec->length;
            ofile[i].nscan = (escan - sscan) / dscan + 1;
            strcpy(ofile[i].l2prod, input->l2prod[i]);
            strcpy(ofile[i].def_l2prod, input->def_l2prod[i]);

            ofile[i].orbit_node_lon = l1file.orbit_node_lon;
            ofile[i].orbit_number = l1file.orbit_number;
            ofile[i].node_crossing_time = l1file.node_crossing_time;
            ofile[i].private_data = l1file.private_data;
 
            printf("Opening: %s\n", ofile[i].name);
            if (openl2(&ofile[i]) != 0) {
                printf("-E- %s: Error opening %s for writing.\n",
                        argv[0], ofile[i].name);
                exit(FATAL_ERROR);
            }
        }
        printf("\n");
    }

    /*                                                                  */
    /* Transfer any additional header info to the record headers        */
    /*                                                                  */
    l2rec->bindx = l1file.bindx;

    printf("\n\nBegin %s Version %d.%d.%d-%s Processing\n", PROGRAM, VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH, GITSHA);
    printf("Sensor is %s\n", sensorId2SensorName(l1file.sensorID));
    printf("Sensor ID is %d\n", l1file.sensorID);
    printf("Sensor has %d reflective bands\n", l1file.nbands);
    printf("Sensor has %d emissive bands\n", l1file.nbandsir);
    printf("Number of along-track detectors per band is %d\n", l1file.ndets);
    printf("Number of input pixels per scan is %d\n", l1file.npix);
    printf("Processing pixels %d to %d by %d\n", spix + 1, epix + 1, dpix);
    printf("Processing scans %d to %d by %d\n", sscan + 1, escan + 1, dscan);

    if (input->proc_ocean != 0)
        printf("Ocean processing enabled\n");
    else
        printf("Ocean processing disabled\n");

    if (input->proc_land != 0)
        printf("Land processing enabled\n");
    else
        printf("Land processing disabled\n");

    if (input->atmocor != 0)
        printf("Atmospheric correction enabled\n");
    else
        printf("Atmospheric correction disabled\n");

    if (aerec->mode == ON) {
        printf("Aerosol parameters will be extracted from %s\n",
                aefile.name);
    }

    start_time = now();
    printf("\nBegin MSl12 processing at %s\n\n", ydhmsf(start_time, 'L'));


    /*                                                                  */
    /*  Read file scan by scan, convert to L2, and write.               */
    /*                                                                  */
    for (iscan = sscan, oscan = 0; iscan <= escan; iscan += dscan, oscan++) {

        /*                                                              */
        /* This call returns the specified record, but it internally    */
        /* buffers enough records to facilitate L1B filtering.          */
        /*                                                              */
        if (getl1rec(iscan, dscan, l1rec) != 0) {
            exit(FATAL_ERROR);
        }

        if ((oscan % 50) == 0)
            printf("Processing scan #%6d (%d of %d) after %6.0f seconds\n",
                iscan, iscan - sscan + 1, escan - sscan + 1,
                now() - start_time);

        if (aerec->mode == ON) {
            if (read_aer(&aefile, iscan, aerec) != 0) {
                printf("-E- %s: Error reading %s at scan %d.\n",
                        argv[0], aefile.name, iscan);
                exit(FATAL_ERROR);
            }
        }

        if (tgrec->mode == ON) {
            if (tgfile.format == FT_L3BIN) {
                read_target_l3(&tgfile, l1rec, l1file.nbands, tgrec);
            } else {
                if (read_target(&tgfile, iscan, tgrec) != 0) {
                    printf("-E- %s: Error reading %s at scan %d.\n",
                            argv[0], tgfile.name, iscan);
                    exit(FATAL_ERROR);
                }
            }
        }

        /*                                                              */
        /* Convert the L1B radiances to L2                              */
        /*                                                              */
        convl12(l1rec, l2rec, 0, l1rec->npix - 1, aerec);

        if (input->mode != FORWARD) {

            /*                                                          */
            /* Recalibration mode. Read target nLw's for this scan and  */
            /* copy into L2 record.  Then reconstruct L1 radiances      */
            /* using the target nLw's and the components of the         */
            /* previous atmospheric correction.                         */
            /*                                                          */
            convl21(l2rec, tgrec, 0, l1rec->npix - 1, l1rec->Lt, NULL);

            /*                                                          */
            /* Write the new L1B record to output file.                 */
            /*                                                          */
            if (writel1(&ofile[0], oscan, l1rec) != 0) {
                printf("-E- %s: error writing to %s\n",
                        argv[0], ofile[0].name);
                exit(FATAL_ERROR);
            }

        } else {
/* Alex S. : writel2 is not just writing the output files but calls prodgen */
            /*                                                          */
            /* Forward mode. Write output record to file(s).            */
            /*                                                          */
            for (i = 0; i < num_ofiles; i++)
                //                if (writel2_hdf( &ofile[i], oscan, l2rec) != 0) {
                if (writel2(&ofile[i], oscan, l2rec, i) != 0) {
                    printf("-E- %s: error writing to %s\n",
                            argv[0], ofile[i].name);
                    exit(FATAL_ERROR);
                }


        }

    }

    printf("\nEnd MSl12 processing at %s\n", ydhmsf(now(), 'L'));
    printf("Processing Rate = %f scans/sec\n\n", ofile[0].nscan / (now() - start_time));

    /*                                                                  */
    /* Write SeaWifs-specific data if appropriate                       */
    /*                                                                  */
    if (l1file.sensorID == SEAWIFS &&
            (l1file.format == FT_SEAWIFSL1A || l1file.format == FT_L1HDF) &&
            (ofile[0].format == FT_L2HDF || ofile[0].format == FT_L1HDF)) {
        printf("Writing SeaWiFS-specific meta-data\n");

        unsigned char genBuf[8192];

        bzero(genBuf, sizeof (genBuf));
        PTB(getHDFattr(l1file.sd_id, "Mission", "", (VOIDP) & genBuf));
        for (i = 0; i < num_ofiles; i++)
            PTB(sd_setattr(ofile[i].sd_id, "Mission", DFNT_CHAR, strlen((const char *) genBuf) + 1, (VOIDP) genBuf));

        bzero(genBuf, sizeof (genBuf));
        PTB(getHDFattr(l1file.sd_id, "Mission Characteristics", "", (VOIDP) & genBuf));
        for (i = 0; i < num_ofiles; i++)
            PTB(sd_setattr(ofile[i].sd_id, "Mission Characteristics", DFNT_CHAR, strlen((const char *) genBuf) + 1, (VOIDP) genBuf));

        bzero(genBuf, sizeof (genBuf));
        PTB(getHDFattr(l1file.sd_id, "Sensor", "", (VOIDP) & genBuf));
        for (i = 0; i < num_ofiles; i++)
            PTB(sd_setattr(ofile[i].sd_id, "Sensor", DFNT_CHAR, strlen((const char *) genBuf) + 1, (VOIDP) genBuf));

        bzero(genBuf, sizeof (genBuf));
        PTB(getHDFattr(l1file.sd_id, "Sensor Characteristics", "", (VOIDP) & genBuf));
        for (i = 0; i < num_ofiles; i++)
            PTB(sd_setattr(ofile[i].sd_id, "Sensor Characteristics", DFNT_CHAR, strlen((const char *) genBuf) + 1, (VOIDP) genBuf));

        bzero(genBuf, sizeof (genBuf));
        PTB(getHDFattr(l1file.sd_id, "Data Type", "", (VOIDP) & genBuf));
        for (i = 0; i < num_ofiles; i++)
            PTB(sd_setattr(ofile[i].sd_id, "Data Type", DFNT_CHAR, strlen((const char *) genBuf) + 1, (VOIDP) genBuf));

        int32 l1a_spixl, l1a_dpixl;
        PTB(getHDFattr(l1file.sd_id, "LAC Pixel Start Number", "", (VOIDP) & l1a_spixl));
        PTB(getHDFattr(l1file.sd_id, "LAC Pixel Subsampling", "", (VOIDP) & l1a_dpixl));

        int32 spixl = (l1a_dpixl * spix) + l1a_spixl;
        int32 dpixl = l1a_dpixl * dpix;

        for (i = 0; i < num_ofiles; i++) {
            PTB(sd_setattr(ofile[i].sd_id, "LAC Pixel Start Number", DFNT_INT32, 1, (VOIDP) & spixl));
            PTB(sd_setattr(ofile[i].sd_id, "LAC Pixel Subsampling", DFNT_INT32, 1, (VOIDP) & dpixl));
        }

    }

    /*                                                                  */
    /* Close all files                                                  */
    /*                                                                  */
    // closing the geo_mask file
    if(input->georegionfile[0]){
        close_georegion_file();
    }
    closel1(&l1file);

    if (input->mode != FORWARD) {
        closel1(&ofile[0]);
        /*
        if (input->mode != INVERSE_ZERO)
            close_target();
         */
    } else
        for (i = 0; i < num_ofiles; i++)
            closel2(&ofile[i]);

    free_l1(l1rec);
    free_l2(l2rec);
    free_l1q();
    free_deminfo();

    // free structures
    free(l1rec);
    free(l2rec);
    free(tgrec);
    free(aerec);

    // free input internals
    for(i=0;i<input->fctl.nfilt; i++)
        free(input->fctl.f[i].kernel);
    free(input->gsm_aphs);
    free(input->gsm_aphw);
    free(input->giop_wave);
    free(input->giop_rrs_unc);
    free(input->taua);
    free(input->vcal_nLw);
    free(input->vcal_Lw);
    free(input);

    printf("\nProcessing Completed\n");

    return (SUCCESS);
}
