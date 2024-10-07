#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include "hdfio.h"

int rd_smmr_orb(char *sfile, int *nrec, int *syear, int *sday,
        double **orbvec, double **time, float **pos_err)
/*******************************************************************

   rd_smmr_orb

   purpose: read in all important smmr orbit file information

   Returns type: int - 0 if no problems
                       if no file exists, return nrec = 0 and don't allocate 
                       orbvec, time

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            sfile            I      name of smmr hdf file
      int *             nrec             O      # of orbit records
      int *             syear            O      year YYYY
      int *             sday             O      day of year DDD
      double **         orbvec           O      pointer to 
                                                size [ nrec, 6 ] array of
                                                position (0-2) and 
                                                velocity (3-5)
      double **         time             O      pointer to
                                                time tag, seconds of the day
      float **          pos_err          O      extimated position error in m

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 15-Sep-2005     Original development

 *******************************************************************/ {
    int i, iret, yrday;
    hdfio_struc hdf_info;
    float *pos, *vel;
    short *min;
    int32 ntype, count;
    struct stat buf;
    /*
     *  open the file
     */
    ntype = DFNT_INT32;
    count = 1;
    iret = 0;

    if (stat(sfile, &buf) != 0) {
        *nrec = 0;
        return iret;
    }
    if (hdfio_open(sfile, &hdf_info) != 0) {
        iret = 1;
        return iret;
    }
    /*
     *  get the attributes first  
     */
    if (hdfio_rd_gattr(hdf_info, "Number of values", ntype, count,
            (void *) nrec) != 0) {
        printf("Unable to read global attribute 'Number_of_values'\n");
        printf("From file %s\n", sfile);
        iret = 1;
        return iret;
    }
    if (hdfio_rd_gattr(hdf_info, "Year Day", ntype, count,
            (void *) &yrday) != 0) {
        printf("Unable to read global attribute 'Year_Day'\n");
        printf("From file %s\n", sfile);
        iret = 1;
        return iret;
    }
    *syear = yrday / 1000 + 1900;
    *sday = yrday % 1000;
    /*
     *  set up the position and velocity arrays and the final orbvec and time arrs
     */
    if ((pos = malloc(3 * *nrec * sizeof ( float))) == NULL) {
        printf("Unable to allocate pos array\n");
        iret = 1;
        return iret;
    }
    if ((vel = malloc(3 * *nrec * sizeof ( float))) == NULL) {
        printf("Unable to allocate pos array\n");
        iret = 1;
        return iret;
    }
    if ((min = malloc(*nrec * sizeof ( short))) == NULL) {
        printf("Unable to allocate min array\n");
        iret = 1;
        return iret;
    }
    /*  */
    if ((*orbvec = malloc(6 * *nrec * sizeof ( double))) == NULL) {
        printf("Unable to allocate orbvec array\n");
        iret = 1;
        return iret;
    }
    if ((*time = malloc(*nrec * sizeof ( double))) == NULL) {
        printf("Unable to allocate time array\n");
        iret = 1;
        return iret;
    }
    if ((*pos_err = malloc(*nrec * sizeof ( float))) == NULL) {
        printf("Unable to allocate pos_err array\n");
        iret = 1;
        return iret;
    }
    /*
     *  read in the time and set it up for output
   NOTE that read_sd dosen't read anything but float, so be ready for enhancements!!
     */
    if (hdfio_rd_sd(hdf_info, "time", (void *) min) != 0) {
        printf(" Unable to read sd 'time'\n");
        iret = 1;
        return iret;
    }

    if (hdfio_rd_sd(hdf_info, "pos", (void *) pos) != 0) {
        printf(" Unable to read sd 'pos'\n");
        iret = 1;
        return iret;
    }
    if (hdfio_rd_sd(hdf_info, "vel", (void *) vel) != 0) {
        printf(" Unable to read sd 'pos'\n");
        iret = 1;
        return iret;
    }
    if (hdfio_rd_sd(hdf_info, "position error", (void *) (*pos_err)) != 0) {
        printf(" Unable to read sd 'position_error'\n");
        iret = 1;
        return iret;
    }
    for (i = 0; i < *nrec; i++) {
        *(*time + i) = ((double) *(min + i));
        *(*orbvec + i * 6) = *(pos + i * 3);
        *(*orbvec + i * 6 + 1) = *(pos + i * 3 + 1);
        *(*orbvec + i * 6 + 2) = *(pos + i * 3 + 2);
        *(*orbvec + i * 6 + 3) = *(vel + i * 3);
        *(*orbvec + i * 6 + 4) = *(vel + i * 3 + 1);
        *(*orbvec + i * 6 + 5) = *(vel + i * 3 + 2);
    }
    /*
     *  and end
     */
    hdfio_close(hdf_info);
    return iret;
}
