#include "l1_nc_generic_read.h"

#include "l1.h"
#include <genutils.h>
#include <timeutils.h>
#include <libnav.h>

#include <netcdf.h>
#include <stdlib.h>

static int16_t nline, npix;
static int32_t spix = 0;

static void scaleFloat(float *buf, int num, float scale, float offset) {
    int i;
    if (scale != 1.0 || offset != 0.0) {
        for (i = 0; i < num; i++) {
            if (buf[i] != -32767.0) {
                buf[i] = buf[i] * scale + offset;
            }
        }
    }
}

int openl1_nc_generic(filehandle *file) {
    size_t source_w, source_h;
    int32_t nscan;
    int fileID, xid, yid, retval;

    // Open the netcdf4 input file
    retval = nc_open(file->name, NC_NOWRITE, &fileID);
    if (retval != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }
    // Get pixel and scan dimensions
    retval = nc_inq_dimid(fileID, "number_of_lines", &yid);
    retval = nc_inq_dimid(fileID, "pixels_per_line", &xid);
    retval = nc_inq_dimlen(fileID, xid, &source_w);
    retval = nc_inq_dimlen(fileID, yid, &source_h);

    if (want_verbose) {
        printf("L1B Npix  :%d Nscans:%d\n", (int) source_w,
                (int) source_h);
    } // want_verbose
    npix = (int32_t) source_w;
    nscan = (int32_t) source_h;
    nline = nscan;

    file->sd_id = fileID;
    file->npix = npix;
    file->nscan = nscan;
    file->terrain_corrected = 1; //assumed
    return (LIFE_IS_GOOD);
}

int readl1_nc_generic(filehandle *file, int32_t scan, l1str *l1rec) {
    static int firstCall = 1;

    static int scanGrp = -1;
    static int dataGrp = -1;
    static int navigationGrp = -1;
    static int ancillaryGrp = -1;

    static int scantimeVar = -1;
    static int yearVar = -1;
    static int dayVar = -1;
    static int msecVar = -1;
    static int lonVar = -1;
    static int latVar = -1;
    static int senaVar = -1;
    static int senzVar = -1;
    static int solaVar = -1;
    static int solzVar = -1;
    static int wsVar = -1;
    static int wdVar = -1;
    static int prVar = -1;
    static int ozVar = -1;
    static int rhVar = -1;
    static int wvVar = -1;
    static int *ltVar;

    static float senaScale = 1.0;
    static float senaOffset = 0.0;
    static float senzScale = 1.0;
    static float senzOffset = 0.0;
    static float solaScale = 1.0;
    static float solaOffset = 0.0;
    static float solzScale = 1.0;
    static float solzOffset = 0.0;

    static int *pixnum;
    static float *rad_data;

    int retval;

    int32_t ip, ib, ipb;
    int32_t nbands = l1rec->l1file->nbands;
    size_t start[3], count[3];
    int i;
    int msec;
    int scan_year, scan_day;

    if (firstCall) {
        if (want_verbose)
            printf("file->nbands = %d\n", (int) file->nbands);
        firstCall = 0;

        int sensorGrp;
        int pixnumVar;
        retval = nc_inq_ncid(file->sd_id, "sensor_band_parameters", &sensorGrp);
        if (retval == NC_NOERR) {
            retval = nc_inq_varid(sensorGrp, "pixnum", &pixnumVar);
            if (retval == NC_NOERR) {
                pixnum = (int*) malloc(npix * sizeof (int));
                start[0] = 0;
                count[0] = npix;
                retval = nc_get_vara_int(sensorGrp, pixnumVar, start, count, pixnum);
                if (retval != NC_NOERR) {
                    fprintf(stderr,
                            "-E- %s line %d: nc_get_vara_int failed for file, %s  group, %s var, %s.\n",
                            __FILE__, __LINE__, file->name, "sensor_band_parameters", "pixnum");
                    exit(EXIT_FAILURE);
                }
            }
        }

        // find all of the group IDs
        retval = nc_inq_ncid(file->sd_id, "scan_line_attributes", &scanGrp);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_inq_ncid failed for file, %s  group, %s.\n",
                    __FILE__, __LINE__, file->name, "scan_line_attributes");
            exit(EXIT_FAILURE);
        }

        retval = nc_inq_ncid(file->sd_id, "geophysical_data", &dataGrp);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_inq_ncid failed for file, %s  group, %s.\n",
                    __FILE__, __LINE__, file->name, "geophysical_data");
            exit(EXIT_FAILURE);
        }

        retval = nc_inq_ncid(file->sd_id, "navigation_data", &navigationGrp);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_inq_ncid failed for file, %s  group, %s.\n",
                    __FILE__, __LINE__, file->name, "navigation_data");
            exit(EXIT_FAILURE);
        }

        retval = nc_inq_ncid(file->sd_id, "ancillary_data", &ancillaryGrp);
        if (retval != NC_NOERR)
            ancillaryGrp = -1;


        // setup all of the variable IDs
        retval = nc_inq_varid(scanGrp, "scantime", &scantimeVar);
        if (retval != NC_NOERR) {
            scantimeVar = -1;
            retval = nc_inq_varid(scanGrp, "year", &yearVar);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: variable \"scantime\" or \"year,day,msec\" required in group \"%s\"\n",
                        __FILE__, __LINE__, "scan_line_attributes");
                exit(EXIT_FAILURE);
            }
            retval = nc_inq_varid(scanGrp, "day", &dayVar);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: variable \"scantime\" or \"year,day,msec\" required in group \"%s\"\n",
                        __FILE__, __LINE__, "scan_line_attributes");
                exit(EXIT_FAILURE);
            }
            retval = nc_inq_varid(scanGrp, "msec", &msecVar);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: variable \"scantime\" or \"year,day,msec\" required in group \"%s\"\n",
                        __FILE__, __LINE__, "scan_line_attributes");
                exit(EXIT_FAILURE);
            }
        } // scantime not found

        retval = nc_inq_varid(navigationGrp, "longitude", &lonVar);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_inq_varid failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, file->name, "longitude");
            exit(EXIT_FAILURE);
        }
        retval = nc_inq_varid(navigationGrp, "latitude", &latVar);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_inq_varid failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, file->name, "latitude");
            exit(EXIT_FAILURE);
        }
        retval = nc_inq_varid(navigationGrp, "sena", &senaVar);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_inq_varid failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, file->name, "sena");
            exit(EXIT_FAILURE);
        }
        retval = nc_inq_varid(navigationGrp, "senz", &senzVar);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_inq_varid failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, file->name, "senz");
            exit(EXIT_FAILURE);
        }
        retval = nc_inq_varid(navigationGrp, "sola", &solaVar);
        if (retval != NC_NOERR)
            solaVar = -1;
        retval = nc_inq_varid(navigationGrp, "solz", &solzVar);
        if (retval != NC_NOERR) {
            if (solaVar == -1) {
                solzVar = -1;
            } else {
                fprintf(stderr,
                        "-E- %s line %d: nc_inq_varid failed for file, %s  field, %s.\n",
                        __FILE__, __LINE__, file->name, "solz");
                fprintf(stderr, "    solz must exist if sola is defined.\n");
                exit(EXIT_FAILURE);
            }
        } else {
            if (solaVar == -1) {
                fprintf(stderr,
                        "-E- %s line %d: nc_inq_varid failed for file, %s  field, %s.\n",
                        __FILE__, __LINE__, file->name, "sola");
                fprintf(stderr, "    sola must exist if solz is defined.\n");
                exit(EXIT_FAILURE);
            }
        }

        if (ancillaryGrp != -1) {
            retval = nc_inq_varid(ancillaryGrp, "windspeed", &wsVar);
            if (retval != NC_NOERR) {
                wsVar = -1;
            }
            retval = nc_inq_varid(ancillaryGrp, "windangle", &wdVar);
            if (retval != NC_NOERR) {
                wdVar = -1;
            }
            retval = nc_inq_varid(ancillaryGrp, "pressure", &prVar);
            if (retval != NC_NOERR) {
                prVar = -1;
            }
            retval = nc_inq_varid(ancillaryGrp, "ozone", &ozVar);
            if (retval != NC_NOERR) {
                ozVar = -1;
            }
            retval = nc_inq_varid(ancillaryGrp, "relhumid", &rhVar);
            if (retval != NC_NOERR) {
                rhVar = -1;
            }
            retval = nc_inq_varid(ancillaryGrp, "watervapor", &wvVar);
            if (retval != NC_NOERR) {
                wvVar = -1;
            }
        } // ancillary group exists

        char bandName[512];
        ltVar = (int*) malloc((nbands + file->nbandsir) * sizeof (int));
        rad_data = (float *) malloc(npix * sizeof (float));
        for (i = 0; i < nbands; i++) {
            sprintf(bandName, "Lt_%d", file->iwave[i]);
            retval = nc_inq_varid(dataGrp, bandName, ltVar + i);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: could not find \"%s\" variable in file, %s\n",
                        __FILE__, __LINE__, bandName, file->name);
                exit(EXIT_FAILURE);
            }
        }
        for (i = nbands; i < (nbands + file->nbandsir); i++) {
            sprintf(bandName, "Lt_%d", file->iwave[i]);
            retval = nc_inq_varid(dataGrp, bandName, ltVar + i);
            if (retval != NC_NOERR) {
                ltVar[i] = -1;
            }
        }

        // read scale factors since sol and sen angles can be scaled shorts
        nc_get_att_float(navigationGrp, senaVar, "scale_factor", &senaScale);
        nc_get_att_float(navigationGrp, senaVar, "add_offset", &senaOffset);

        nc_get_att_float(navigationGrp, senzVar, "scale_factor", &senzScale);
        nc_get_att_float(navigationGrp, senzVar, "add_offset", &senzOffset);

        if (solaVar != -1) {
            nc_get_att_float(navigationGrp, solaVar, "scale_factor", &solaScale);
            nc_get_att_float(navigationGrp, solaVar, "add_offset", &solaOffset);
        }

        if (solzVar != -1) {
            nc_get_att_float(navigationGrp, solzVar, "scale_factor", &solzScale);
            nc_get_att_float(navigationGrp, solzVar, "add_offset", &solzOffset);
        }

    } // first call


    if (pixnum) {
        for (ip = 0; ip < npix; ip++)
            l1rec->pixnum[ip] = pixnum[ip];
    }

    start[0] = scan;
    count[0] = 1;

    if (scantimeVar != -1) {
        nc_get_vara_double(scanGrp, scantimeVar, start, count, &l1rec->scantime);
    } else {
        retval = nc_get_vara_int(scanGrp, yearVar, start, count, &scan_year);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_int failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, file->name, "year");
            exit(EXIT_FAILURE);
        }
        retval = nc_get_vara_int(scanGrp, dayVar, start, count, &scan_day);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_int failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, file->name, "day");
            exit(EXIT_FAILURE);
        }
        retval = nc_get_vara_int(scanGrp, msecVar, start, count, &msec);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_int failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, file->name, "msec");
            exit(EXIT_FAILURE);
        }

        l1rec->scantime = yds2unix((int16_t) scan_year, (int16_t) scan_day, (double) (msec / 1.e3));
    }

    start[0] = scan;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = npix;
    count[2] = 1;

    retval = nc_get_vara_float(navigationGrp, lonVar, start, count, l1rec->lon);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "longitude");
        exit(EXIT_FAILURE);
    }
    retval = nc_get_vara_float(navigationGrp, latVar, start, count, l1rec->lat);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "latitude");
        exit(EXIT_FAILURE);
    }
    retval = nc_get_vara_float(navigationGrp, senaVar, start, count, l1rec->sena);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "sena");
        exit(EXIT_FAILURE);
    }
    scaleFloat(l1rec->sena, npix, senaScale, senaOffset);
    retval = nc_get_vara_float(navigationGrp, senzVar, start, count, l1rec->senz);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "senz");
        exit(EXIT_FAILURE);
    }
    scaleFloat(l1rec->senz, npix, senzScale, senzOffset);
    if (solaVar != -1) {
        retval = nc_get_vara_float(navigationGrp, solaVar, start, count, l1rec->sola);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, file->name, "sola");
            exit(EXIT_FAILURE);
        }
        scaleFloat(l1rec->sola, npix, solaScale, solaOffset);
        retval = nc_get_vara_float(navigationGrp, solzVar, start, count, l1rec->solz);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, file->name, "solz");
            exit(EXIT_FAILURE);
        }
        scaleFloat(l1rec->solz, npix, solzScale, solzOffset);
    } else {
        int iyear, iday;
        float gmt;
        int16_t syear, sday;
        double secs;
        int i;

        unix2yds(l1rec->scantime, &syear, &sday, &secs);
        iyear = syear;
        iday = sday;
        gmt = secs / 3600.0;
        for (i = 0; i < npix; i++) {
            sunangs_(&iyear, &iday, &gmt, l1rec->lon + i, l1rec->lat + i, l1rec->solz + i, l1rec->sola + i);
        }
    }

    // read in ancillary data
    if (ancillaryGrp) {

        if (wsVar != -1) {
            retval = nc_get_vara_float(ancillaryGrp, wsVar, start, count, l1rec->ws);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                        __FILE__, __LINE__, file->name, "windspeed");
                exit(EXIT_FAILURE);
            }
        }

        if (wdVar != -1) {
            retval = nc_get_vara_float(ancillaryGrp, wdVar, start, count, l1rec->wd);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                        __FILE__, __LINE__, file->name, "windangle");
                exit(EXIT_FAILURE);
            }
        }

        if (prVar != -1) {
            retval = nc_get_vara_float(ancillaryGrp, prVar, start, count, l1rec->pr);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                        __FILE__, __LINE__, file->name, "pressure");
                exit(EXIT_FAILURE);
            }
        }

        if (ozVar != -1) {
            retval = nc_get_vara_float(ancillaryGrp, ozVar, start, count, l1rec->oz);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                        __FILE__, __LINE__, file->name, "ozone");
                exit(EXIT_FAILURE);
            }
        }

        if (rhVar != -1) {
            retval = nc_get_vara_float(ancillaryGrp, rhVar, start, count, l1rec->rh);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                        __FILE__, __LINE__, file->name, "relhumid");
                exit(EXIT_FAILURE);
            }
        }

        if (wvVar != -1) {
            retval = nc_get_vara_float(ancillaryGrp, wvVar, start, count, l1rec->wv);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                        __FILE__, __LINE__, file->name, "watervapor");
                exit(EXIT_FAILURE);
            }
        }
    } // found ancillary group


    // read in radiance data
    for (ib = 0; ib < nbands; ib++) {

        retval = nc_get_vara_float(dataGrp, ltVar[ib], start, count, rad_data);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  Lt %d nm.\n",
                    __FILE__, __LINE__, file->name, file->iwave[ib]);
            exit(EXIT_FAILURE);
        }

        // copy to Lt record.
        for (ip = spix; ip < npix; ip++) {
            ipb = ip * nbands + ib;
            l1rec->Lt[ipb] = rad_data[ip] / 10.0;

            if (l1rec->Lt[ipb] <= 0.0)
                l1rec->Lt[ipb] = 0.0001;
        }
    } // for ib

    // now for IR bands
    for (ib = nbands; ib < (nbands + file->nbandsir); ib++) {
        if (ltVar[ib] != -1) {
            retval = nc_get_vara_float(dataGrp, ltVar[ib], start, count, rad_data);
            if (retval == NC_NOERR) {
                // copy to Lt record.
                for (ip = spix; ip < npix; ip++) {
                    ipb = ip * file->nbandsir + ib - nbands;
                    l1rec->Ltir[ipb] = rad_data[ip] / 10.0;

                    if (l1rec->Ltir[ipb] <= 0.0)
                        l1rec->Ltir[ipb] = 0.0001;
                }
            }
        }

    } // for ib

    l1rec->npix = file->npix;

    return (LIFE_IS_GOOD);
}

int closel1_nc_generic(filehandle *file) {
    int retval;

    retval = nc_close(file->sd_id);
    if (retval != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_close failed for file, %s.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }

    return (LIFE_IS_GOOD);
}

