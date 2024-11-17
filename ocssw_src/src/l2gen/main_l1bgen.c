/* ========================================================================
 * MSl1bgen - multi-sensor L1B generator
 * 
 * Synopsis:
 *
 *   MSl1bgen par=parfile
 *
 * Description:
 * 
 * Modification history:
 *
 *     Programmer     Organization      Date       Description of change
 *   --------------   ------------    --------     ---------------------
 *   Bryan A. Franz   GSC             21 July 1998 Original development
 *   Joel M. Gales    Futuretech      20 Sept 1999 Generate standard L1B
 *                                                 SeaWIFS output
 *   Bryan A. Franz   GSC             28 Jan  2000 switch to parfile and
 *                                                 added stray-light control
 *
 * ======================================================================== */

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h>

#include "l12_proto.h"

#define INT32   int32_t 
#define FLOAT32 float
#define BYTE    unsigned char

int msl1bgen_usage(char *prog);

/* -------------------------------------------------------------------- *
 *                              main                                    *
 * -------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
    int32_t iscan = 0; /* input scan number                  */
    int32_t oscan = 0; /* input scan number                  */
    int32_t spix = 0; /* start pixel for subscene process   */
    int32_t epix = -1; /* end pixel for subscene process     */
    int32_t dpix = 1; /* pixel subsampling increment        */
    int32_t sscan = 0; /* start scan for subscene process    */
    int32_t escan = -1; /* end scan for subscene process      */
    int32_t dscan = 1; /* scan subsampling increment         */
    int32_t npix = 0; /* Number of output pixels per scan   */

    l1str l1rec; /* generic level-1b scan structure    */
    filehandle ifile; /* input file handle                  */
    filehandle ofile; /* output file handle                 */

    int32_t ip, ib, ipb, i;

    if (argc == 1) {
        l2gen_usage("l1bgen_generic");
        return 1;
    }

    // see if help on command line
    for (i = 0; i < argc; i++) {
        if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-help") == 0)) {
            l2gen_usage("l1bgen_generic");
        }
    }

    // cdata_();
    filehandle_init(&ifile);
    filehandle_init(&ofile);
    msl12_input_init();

    /* Parse input parameters */
    if (msl12_input(argc, argv, "l1bgen_generic", &ifile) != 0) {
        printf("-E- %s: Error parsing input parameters.\n", argv[0]);
        exit(FATAL_ERROR);
    }

    if (access(input->ifile[0], F_OK) || access(input->ifile[0], R_OK)) {
        printf("-E- %s: Input file '%s' does not exist or cannot open.\n",
                argv[0], input->ifile[0]);
        exit(FATAL_ERROR);
    }

    /*									*/
    /* Open input file and get sensor and scan information from handle  */
    /*                                                                  */
    if (openl1(&ifile) != 0) {
        printf("-E- %s: Error opening %s for reading.\n", argv[0], ifile.name);
        exit(1);
    }

    /*									*/
    /* Allocate memory for L1 scan data			  	        */
    /*									*/
    if (alloc_l1(&ifile, &l1rec) == 0) {
        printf("-E- %s: Unable to allocate L1 record.\n", argv[0]);
        exit(FATAL_ERROR);
    }


    /* Set the end pixel if it was not set by command argument	        */
    if (l1_input->epixl == -1)
        l1_input->epixl = ifile.npix;
    if (l1_input->eline == -1)
        l1_input->eline = ifile.nscan;

    spix = MAX(l1_input->spixl - 1, 0);
    epix = MIN(l1_input->epixl - 1, ifile.npix - 1);
    dpix = MAX(l1_input->dpixl, 1);
    sscan = MAX(l1_input->sline - 1, 0);
    escan = MIN(l1_input->eline - 1, ifile.nscan - 1);
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
    ifile.spix = spix; /* start pixel rel to L1 scan    */
    ifile.epix = epix; /* end   pixel rel to L1 scan    */

    npix = (epix - spix) / dpix + 1;

    /*			                                           */
    /* Transfer sensor and scan info to output filehandle and open */
    /*                                                             */
    strcpy(ofile.name, input->ofile[0]);
    if (strcmp(input->oformat, "netCDF4") == 0)
        ofile.format = FT_L1BNCDF;
    else
        ofile.format = FT_L1HDF;
    ofile.mode = WRITE;
    ofile.sensorID = ifile.sensorID;
    ofile.nbands = ifile.nbands;
    ofile.bindx = ifile.bindx;
    ofile.spix = spix;
    ofile.epix = epix;
    ofile.npix = (epix - spix) / dpix + 1;
    ofile.nscan = (escan - sscan) / dscan + 1;

    strcpy(ofile.spatialResolution, ifile.spatialResolution);


    printf("Opening L1B output file: %s\n", ofile.name);

    if (openl1(&ofile) != 0) {
        printf("-E- %s: Error opening %s for writing.\n", argv[0], ofile.name);
        exit(1);
    }

    /*					 			*/
    /* Read file scan by scan, scale radiances,  and write.     */
    /*								*/
    for (iscan = sscan; iscan <= escan; iscan += dscan, oscan++) {

        if ((iscan % 100) == 0) printf("Processing scan %d\n", iscan);

        readl1(&ifile, iscan, &l1rec);

        /* apply vicarious calibration */
        for (ip = 0; ip < npix; ip++) {
            for (ib = 0; ib < ifile.nbands; ib++) {
                i = ifile.bindx[ib];
                ipb = ip * ifile.nbands + i;
                l1rec.Lt[ipb] = l1rec.Lt[ipb] * l1_input->gain[ib] + l1_input->offset[ib];
            }
        }

        /* write this record */
        if (writel1(&ofile, oscan, &l1rec) != 0) {
            printf("-E- %s: error writing to %s.\n",
                    argv[0], ofile.name);
            exit(FATAL_ERROR);
        }

    }


    //    if (ifile.sensorID == SEAWIFS) {
    //      DPTB(l1b_seawifs(&ifile, &ofile,
    //			 sscan, escan, dscan,
    //			 spix,  epix,  dpix));
    //    }

    closel1(&ifile);
    closel1(&ofile);

    printf("\nProcessing Completed\n");

    exit(0);
}

