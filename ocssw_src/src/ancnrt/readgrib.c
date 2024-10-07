#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <timeutils.h>

int readgrib(char *file, int npix, int nlin, float *data, int *year,
        int *month, int *day, int *hour)
/*******************************************************************

   read_grib

   purpose: read in the data from a binary grib file with header
      and pds

   Returns type: void

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      binary grib file name
      int               npix             I      # pixels in grid
      int               nlin             I      # lines in grid
      float *           data             O      returned data
      int *             year             O      year
      int *             month            O      month
      int *             day              O      day
      int *             hour             O      hour

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       17-Sep-2003     Original development - replacement
                                        for fortran version and upgrade
                                        for use with wgrib program
      W. Robinson, SAIC 24 Jan 2008     Adapt to call rd_grib_grid to 
                                        get the data size and data grid
      W. Robinson, SAIC 24 Jun 2014     have variable # pixels, lines

 *******************************************************************/ {
    FILE *stream;
    int retcode;
    int rd_grib_pds(FILE *, int *, int *, int *, int *);
    int rd_grib_grid(FILE *, int, int, float *);
    /*
     *  First, open the file
     */
    retcode = 0;
    if ((stream = fopen(file, "r")) == NULL) {
        printf("Unable to open GRIB file: %s, Exiting\n", file);
        retcode = 1;
    } else {
        /*
         *  read the PDS / GDS and get the date and time
         */
        if (rd_grib_pds(stream, year, month, day, hour) != 0) {
            printf(
                    "Error on PDS portion of grib file read for file: %s, Exiting\n", file);
            retcode = 1;
        } else {
            /*
             * check size and read grid
             */
            if (rd_grib_grid(stream, npix, nlin, data) != 0) {
                printf("%s: rd_grib_grid failed for file: %s\n", __FILE__, file);
                retcode = 1;
            }
        }
    }
    fclose(stream);

    /*
     *  close the data file and return hopefully with the data
     */
    return retcode;
}

int readgrib2(char *file, int npix, int nlin, int rgmode,
        float *data)
/*******************************************************************

   readgrib2

   purpose: read in the data from a binary grib 2 file generated 
      using wgrib2 with the -bin option

   Returns type: int 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      binary grib file name
      int               npix             I      # pixels in grid
      int               nlin             I      # lines in grid
      int               rgmode           I      read mode: 0 do not invert in 
                                                latitude (aquarius), 1 invert 
                                                (std ancnrt)
      float *           data             O      returned data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 24-Jan-2008     Original development 
      W. Robinson, SAIC 18 Sep 2009     add a rgmode: 0 do not invert in 
                                      latitude (aquarius), 1 invert (std ancnrt)
      W. Robinson, SAIC 10-Dec-2009     use variable grid size and extract
                                        time computation
      W. Robinson, SAIC 18-Mar-2015   remove the grib string interpretation

 *******************************************************************/ {
    FILE *stream;
    int retcode, ilin;
    float *xfr_arr;
    int rd_grib_grid(FILE *, int, int, float *);
    /*
     *  First, open the file
     */
    retcode = 0;
    if ((stream = fopen(file, "r")) == NULL) {
        printf("Unable to open GRIB file: %s, Exiting\n", file);
        retcode = 1;
    } else {
        /*
         * check size and read grid
         */
        if (rd_grib_grid(stream, npix, nlin, data) != 0) {
            printf("%s: rd_grib_grid failed for file: %s\n", __FILE__, file);
            retcode = 1;
        }
        /*
         *  close the data file
         */
        fclose(stream);
    }
    /*
     *  for ancnrt use (rgmode = 1) and if no earlier problems, the grid 
     *  must be reversed in latitude to match grib 1 grid
     */
    if ((rgmode == 1) && (retcode == 0)) {
        if ((xfr_arr = (float *) malloc(npix * sizeof ( float))) == NULL) {
            printf("%s: Transfer space allocation failed\n", __FILE__);
            retcode = 1;
        } else {
            for (ilin = 0; ilin < nlin / 2; ilin++) {
                memcpy((void *) xfr_arr, (void *) (data + ilin * npix),
                        npix * sizeof ( float));
                memcpy((void *) (data + ilin * npix),
                        (void *) (data + (nlin - 1 - ilin) * npix),
                        npix * sizeof ( float));
                memcpy((void *) (data + (nlin - 1 - ilin) * npix),
                        (void *) xfr_arr, npix * sizeof ( float));
            }
            free(xfr_arr);
        }
    }
    return retcode;
}

int grib2_t(char *grib2_t_str, int *year, int *doy, int *hour, int *npix,
        int *nlin, int *h_fcst)
/*******************************************************************

   grib2_t

   purpose: get the data time, and maybe size from information for a 
      GRIB2 file

   Returns type: int 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            grib2_t_str      I      string containing grib file
                                                date and time (from wgrib2 -t)
      int *             year             O      year
      int *             doy              O      day of year
      int *             hour             O      hour
      int *             npix             O      # pixels, if not gotten 
                                                  here, return -1
      int *             nlin             O      # lines, if not gotten
                                                  here, return -1

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 10-Dec-2009     break out of readgrib2...
      W. Robinson, SAIC 18-Mar-2015     enhance to just read the time data
                       gotten with  wgrib2 -t <file>
                       or time and grid size, forecast time gotten with
                       wgrib2 -t -nxny -ftime
                       the time is adjusted by the forecast time.

 *******************************************************************/ {
    int retcode, nchr, month, day;
    char *equal_ptr, *nex_ptr, *nex_ptr2, tmp_str[200];
    double unix_t, sec;
    retcode = 0;
    /*
     * The form of the time string can be: <chars>=YYYYMMDDHH
     * and optionally followed by:  :(<# pix> x <# lin>):<fcst>
     * with <# pix>, <# lin> an integer grid size in oix, lin and
     *  <fcst> either 'anl' for 0 h forecast or an integer # hours forecast
     */
    /*
     *  The time id of form YYYYMMDDHH after the last '=' in the string
     */
    if ((equal_ptr = strrchr(grib2_t_str, '=')) == NULL) {
        printf("%s: Unable to find start of time in GRIB 2 time string of: %s\n",
                __FILE__, grib2_t_str);
        retcode = 1;
    } else {
        if (sscanf(equal_ptr, "=%4d%2d%2d%2d", year, &month, &day, hour) != 4) {
            printf("%s: Unable to decode time in GRIB 2 time string of: %s\n",
                    __FILE__, grib2_t_str);
            retcode = 1;
        }
    }
    /*
     *  OK, for the optional portion
     */
    if ((nex_ptr = strchr(equal_ptr, '(')) == NULL) {
        *h_fcst = 0; /*  set to no forecast and fill for size */
        *npix = -1;
        *nlin = -1;
    } else {
        /*  find, read # pixels, lines  */
        if ((nex_ptr2 = strchr(nex_ptr, ')')) == NULL) {
            printf("%s, %d: wgrib2 info string has no ')' in grid size portion\n",
                    __FILE__, __LINE__);
            return 1;
        }
        nchr = nex_ptr2 - nex_ptr + 1;
        strncpy(tmp_str, nex_ptr, nchr);
        if (sscanf(tmp_str, "(%d x %d)", npix, nlin) != 2) {
            printf("%s, %d: Unable to convert grid size info\n",
                    __FILE__, __LINE__);
            return 1;
        }
        /*  get forecast time */
        if ((nex_ptr = strchr(nex_ptr2, ':')) == NULL) {
            printf("%s, %d: wgrib2 info string has no forecast time info\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if (strcmp(nex_ptr, ":anl") == 0)
            *h_fcst = 0;
        else {
            if (sscanf(nex_ptr, ":%d", h_fcst) != 1) {
                printf("%s, %d: Unable to convert forecast time info\n",
                        __FILE__, __LINE__);
                return 1;
            }
        }
    }
    /*
     *  The only thing left is to account for the forecast time or
     *  at least make day of year time
     */
    sec = (double) *hour * 3600;
    unix_t = ymds2unix(*year, month, day, sec);
    unix_t += (double) *h_fcst * 3600;
    unix2yds(unix_t, (int16_t *) year, (int16_t *) doy, &sec);
    *hour = sec / 3600;
    /*  and end */
    return retcode;
}

int readgrib2_3d(char *file, char *grib2_t_str, int npix, int nlin,
        float *data, int nprfl, int *year, int *month, int *day, int *hour)
/*******************************************************************

   readgrib2_3d

   purpose: read in the data from a 3D binary grib 2 file generated 
      using wgrib2 with the -bin option

   Returns type: int 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      binary grib file name
      char *            grib2_t_str      I      string containing grib file
                                                date and time (from wgrib2
                                                -t)
      int               npix             I      # pixels in grid
      int               nlin             I      # lines in grid
      float *           data             O      returned data
      int               nprfl            I      # levels in grid
      int *             year             O      year
      int *             month            O      month
      int *             day              O      day
      int *             hour             O      hour

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 24-Jan-2008     Original development 
      J. Gales,    F.T. 15-May-2008     Add support for profile (3d) files
      W. Robinson, SAIC 10-Dec-2009     use variable grid size and extract
                                        time computation
      W. Robinson, SAIC 19 Mar 2015     adapt for new grib2_t call seq

 *******************************************************************/ {
    FILE *stream;
    int retcode, iprfl, offset, doy, npix_n, nlin_n, h_fcst;
    int rd_grib_grid(FILE *, int, int, float *);
    int grib2_t(char *, int *, int *, int *, int *, int *, int *);
    /*
     *  First, open the file
     */
    retcode = 0;
    if ((stream = fopen(file, "r")) == NULL) {
        printf("Unable to open GRIB file: %s, Exiting\n", file);
        retcode = 1;
    } else {
        /*
         *  and get the date and time from the grib 2 time string
         */
        if (grib2_t(grib2_t_str, year, &doy, hour, &npix_n, &nlin_n, &h_fcst) != 0)
            retcode = 1;
        else {
            /*
             * check size and read grid
             */
            yd2md(*year, doy, (int16_t *) month, (int16_t *) day);
            offset = 0;
            for (iprfl = 0; iprfl < nprfl; iprfl++) {
                if (rd_grib_grid(stream, npix, nlin, &data[offset]) != 0) {
                    printf("%s: rd_grib_grid failed for file: %s\n", __FILE__, file);
                    retcode = 1;
                }
                //	    offset += (nlin * npix) + 2;
                offset += (nlin * npix);
            }
        }
    }
    /*
     *  close the data file
     */
    fclose(stream);

    return retcode;
}

int rd_grib_grid(FILE *stream, int npix, int nlin, float *data)
/*******************************************************************

   rd_grib_grid

   purpose: with an open binary file with grid, read / check the size in 
      bytes and read the data grid

   Returns type: int - 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      FILE *            stream           I      ID of opened file
      int               npix             I      # pixels in grid
      int               nlin             I      # lines in grid
      float *           data             O      data grid read

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       24 Jan 2008     Original development 
 *******************************************************************/ {
    size_t rd_count;
    int expect_ct, retcode, in_nbyt, nbyt;
    int32_t npt_grid;
    void rpt_err_typ(int, int, FILE *);
    /*
     *  read the dimensions
     */
    retcode = 0;
    nbyt = npix * nlin * sizeof ( float);
    expect_ct = 1; /* note that wgrib2 puts out a 4-byte integer 
                       count of the # bytes of data in the binary data -
                       that is checked first */
    if ((rd_count = fread(&in_nbyt, sizeof ( int), expect_ct, stream))
            != expect_ct) {
        printf("Error on header read of GRIB file, Exiting\n");
        rpt_err_typ(rd_count, expect_ct, stream);
        retcode = 1;
    } else {
        /*
         *  We'll assume size 144, 73 now and check that
         */
        if (in_nbyt != nbyt) {
            printf("current grid size of %d bytes is not the expected size of 144 X 73 = %d bytes, Exiting\n", in_nbyt, nbyt);
            retcode = 1;
        } else {
            /*
             *  read the grid
             */
            npt_grid = npix * nlin;
            if ((rd_count = fread(data, sizeof ( float), npt_grid, stream)) != npt_grid) {
                printf("Error reading grid from GRIB file, Exiting\n");
                rpt_err_typ(rd_count, npt_grid, stream);
                retcode = 1;
            }
        }

        if ((rd_count = fread(&in_nbyt, sizeof ( int), expect_ct, stream))
                != expect_ct) {
            printf("Error on footer read of GRIB file, Exiting\n");
            rpt_err_typ(rd_count, expect_ct, stream);
            retcode = 1;
        }
    }
    return retcode;
}

int rd_grib_pds(FILE *stream, int* year, int* month, int* day, int* hour)
/*******************************************************************

   read_grib

   purpose: 1 - read the date and time from the PDS in the grib data and, 
            2 - position the file so the rest of the data can be read

   Returns type: int 0 if OK read, 1 if problem

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      FILE *            stream           I      file stream id
      int *             year             O      year
      int *             month            O      month
      int *             day              O      day
      int *             hour             O      hour

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       23-Sep-2003     Original development 
      W. Robinson, SAIC 16-Mar-2010     # pixel, line possibility

 *******************************************************************/ {
    int retcode, expect_ct, rd_count, lenarr[2];
    unsigned char pds_gds[100];
    void rpt_err_typ(int, int, FILE *);
    /*
     * read the length of the PDS + 4 (bytes) and the char str 'PDS '
     * into 2 4 byte integers
     */
    retcode = 0;
    expect_ct = 2;
    if ((rd_count = fread(lenarr, sizeof ( int), expect_ct, stream)) !=
            expect_ct) {
        printf("Error on 1st PDS read\n");
        rpt_err_typ(rd_count, expect_ct, stream);
        retcode = 1;
    } else {
        /*
         * read the pds as bytes and the duplicate length
         * the char str of 'PDS ' from prev read could be checked but...
         */
        expect_ct = lenarr[0];
        if ((rd_count =
                fread(pds_gds, sizeof ( unsigned char), expect_ct, stream))
                != expect_ct) {
            printf("Error on 2nd PDS read\n");
            rpt_err_typ(rd_count, expect_ct, stream);
            retcode = 1;
        } else {
            /*
             *  The time info is in the 4th set of 8 bytes
             */
            *year = pds_gds[12];
            *month = pds_gds[13];
            *day = pds_gds[14];
            *hour = pds_gds[15];
            /*
             *  lastly, read the GDS info and the GDS bytes
             */
            expect_ct = 2;
            if ((rd_count = fread(lenarr, sizeof ( int), expect_ct, stream)) !=
                    expect_ct) {
                printf("Error on 3rd PDS read\n");
                rpt_err_typ(rd_count, expect_ct, stream);
                retcode = 1;
            } else {
                expect_ct = lenarr[0];
                if ((rd_count =
                        fread(pds_gds, sizeof ( unsigned char), expect_ct, stream))
                        != expect_ct) {
                    printf("Error on 4th PDS read\n");
                    rpt_err_typ(rd_count, expect_ct, stream);
                    retcode = 1;
                }
                /*
                 *  # pixels, lines available in GDS (for reg grid) and this will
                 *  show the getting of them if needed in future
                 */
                //npix = 256 * (int)pds_gds[6] + (int)pds_gds[7];
                //nlin = 256 * (int)pds_gds[8] + (int)pds_gds[9];
                /*
                 printf( "%s, %d, TEST - from GDS, # pixels: %d, # lines: %d\n",
                   __FILE__, __LINE__, npix, nlin );
                 */
            }
        }
    }
    return retcode;
}

/*
 *  Just to report the errors in fread
 *  int rd_count - count of what was read
 *  int expect_ct - count of what was expected
 *  FILE *stream - file id from open
 */
void rpt_err_typ(int rd_count, int expect_ct, FILE *stream) {
    int errcod;
    if ((errcod = ferror(stream)) != 0) {
        printf("Error code of %d encountered\n", errcod);
    } else if (feof(stream) != 0) {
        printf("EOF encountered\n");
    } else {
        printf("Only %d values were read instead of the expected %d values\n",
                rd_count, expect_ct);
    }
    return;
}
