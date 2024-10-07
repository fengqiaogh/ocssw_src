#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <sstream>
#include <iomanip> 
#include <stdint.h>
#include <iostream>
#include <fstream>

#include "cdl_utils.h"
#include "netcdf.h"
#include "nc4utils.h"
#include "nc_gridutils.h"
#include "prism.h"
#include "rdprism.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort_double.h>
#include <timeutils.h>

#define SENSOR "ocip"
#define CDLFILE "OCIP_Level-1B_Data_Structure.cdl"
#define NBANDS 285 // number of OCI bands

extern "C" void handle_error(int status);

using namespace std;
//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     FutureTech     12/18/14 0.01  Original development
//  Joel Gales     FutureTech     01/12/15 0.10  Add support for writing
//                                               interpolated ORCA Lt values
//                                               and navigation
//  Rick Healy     SAIC           06/06/16 1.00  Integrated L2gen L1 reader for Aviris
//  Rick Healy     SAIC           07/15/16 1.00  Integrated L2gen L1 reader for Prism


#define VERSION "1.00"

size_t interpLt(size_t startJ, size_t npixels, size_t k, size_t nbMax,
        float ociaWavelength, size_t iwave, prism4ocia_t* data,
        gsl_interp_accel* acc, float* Lt) {
    size_t j, n;
    size_t i0 = startJ, i1 = startJ + 5;
    float minWave, maxWave;

    double* x = (double*) (calloc(nbMax, sizeof (double)));
    double* y = (double*) (calloc(nbMax, sizeof (double)));
    size_t iter = 0;
    while (i0 >= 0 && i1 < nbMax) {
        n = 0;
        minWave = 99999;
        maxWave = -1;
        for (j = i0; j < i1; j++) {
            if (data->Lt[j * npixels + k] > 0) {
                x[n] = data->wave[j];
                y[n] = data->Lt[j * npixels + k];
                n++;
                if (minWave > data->wave[j]) minWave = data->wave[j];
                if (maxWave < data->wave[j]) maxWave = data->wave[j];
            }
        }

        if (n > 2 && minWave <= ociaWavelength && maxWave >= ociaWavelength) break;
        if (minWave > ociaWavelength || n < 3)
            i0--;
        if (maxWave < ociaWavelength || n < 3)
            i1++;
        iter++;
    }
    if (n > 2) {
        // Allocate space for spline
        gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, n);
        gsl_sort2(x, 1, y, 1, n);
        if (ociaWavelength > x[0]) {
            // Sort wavelengths and corresponding Lt values
            gsl_sort2(x, 1, y, 1, n);
            // Initiate spline
            gsl_spline_init(spline, x, y, n);
            // Generate interpolated values
            double yi = gsl_spline_eval(spline, ociaWavelength, acc);
            Lt[iwave * npixels + k] = yi;
        } else {
            Lt[iwave * npixels + k] = -9999;
        }
        iwave++;
        gsl_spline_free(spline);
    }
    free(x);
    free(y);
    return iwave;
}

int main(int argc, char* argv[]) {
    int status;
    size_t sline = 0, eline = 0;
    char filename[FILENAME_MAX];
    char *filedir;
    double scantime_first, scantime_last;
    static int first_good_scantime = 0;
    char isodatetime_first[30], isodatetime_last[30];

    cout << "prismbil2oci " << VERSION << " ("
            << __DATE__ << " " << __TIME__ << ")" << endl;

    if (argc < 3) {
        cout << endl <<
                "prismbil2oci input_PRISM_file output_ORCA_file <sline (default=0)> <eline (default=last row)>"
                << endl;
        return 0;
    }

    if (argc >= 4) {
        sscanf(argv[3], "%ld", &sline);
        if (sline < 0) sline = 0;
    }

    if (argc >= 5) {
        sscanf(argv[4], "%ld", &eline);
        if (eline < 0) eline = 0;
    }
    //  prism4orca_t *&data;

    if ((filedir = getenv("OCDATAROOT")) == NULL) {
        printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
        return (-1);
    }
    strcpy(filename, filedir);

    strcat(filename, "/");
    strcat(filename, SENSOR);
    strcat(filename, "/");

    strcat(filename, CDLFILE);

    prism4ocia_t **data, *temp;
    float fwhm[NBANDS];
    data = (prism4ocia_t **) malloc(sizeof (prism4ocia_t *));
    temp = open_prism(argv[1], data);
    size_t nscans, npixels = temp->npix;

    if (eline <= 0) eline = temp->nscan;

    nscans = eline - sline;

    int ncid, grpID;
    //char** dim_names, uint32_t* dim_size, size_t n_dims
    const char *dimNames[] = {"number_of_lines", "pixels_per_line", "number_of_bands"};
    size_t dimSize[] = {nscans, npixels, NBANDS};
    size_t n_dims = 3;
    char time_coverage_start[64] = "yyyy-mm-ddThh:mm:ssZ";

    // Create ORCA output file
    static ncdfFile ociafile;
    ociafile.cdlCreateDim(argv[2], filename, dimNames, dimSize, n_dims, nscans);

    //  nc_def_dim(orcafile.ncid, "number_of_lines", nscans, &dimid[0]);
    //  nc_def_dim(orcafile.ncid, "pixels_per_line", nscans, &dimid[1]);
    //
    // Copy prism nav fields to ORCA
    //
    ncid = ociafile.getNcid();

    sprintf(time_coverage_start, "%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%2.2dZ", temp->year, temp->month, temp->day, temp->hour, temp->min, (int) temp->sec);
    status = nc_put_att_text(ncid, NC_GLOBAL, "time_coverage_start", strlen(time_coverage_start) + 1, time_coverage_start);
    //  status = nc_put_att_text(ncid, NC_GLOBAL, "time_coverage_end", strlen(time_coverage_start)+1, time_coverage_start);
    if (status != NC_NOERR) {
        printf("-E- %s %d: %s for %s\n",
                __FILE__, __LINE__, nc_strerror(status), "time_coverage_start");
        exit(1);
    }

    sprintf(time_coverage_start, "%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%2.2dZ", temp->year, temp->month, temp->day, temp->hour, temp->min, (int) temp->sec);
    status = nc_put_att_text(ncid, NC_GLOBAL, "time_coverage_end", strlen(time_coverage_start) + 1, time_coverage_start);



    size_t nbands_av = temp->numBands;
    size_t nbands_ocia = NBANDS;

    size_t n = 0, iwave = 0;
    istringstream istr;
    string str;


    // Allocate space for prism Lt (short) values

    gsl_interp_accel *acc = gsl_interp_accel_alloc();

    float *Lt = (float *) calloc(npixels*nbands_ocia, sizeof (float));
    float *ociaWavelength = (float *) calloc(nbands_ocia, sizeof (float));
    double* utc = (double*) (calloc(nscans, sizeof (double)));
    size_t j;

    const char *navfieldsOCIP[] = {"sena", "senz", "sola", "solz", "lat", "lon"};

    int varID;

    size_t startOCIP[2] = {0, 0};
    size_t countOCIP[2] = {1, npixels};
    size_t i;


    for (size_t iscan = sline; iscan < eline; iscan++) {

        grpID = ociafile.getGid("observation_data");

        // Read prism Lt's for this scan
        if (read_prism(temp, iscan)) {

            cout << "Can't read file " << argv[1] << endl;
            exit(1);

        }

        if ((iscan % 100) == 0) {
            cout << "Processing scan: " << iscan << " out of " << nscans << endl;
            printf("year=%d doy=%d msec=%d month=%d day=%d hour=%d min=%d sec=%f\n ", temp->year, temp->doy, temp->msec, temp->month, temp->day, temp->hour, temp->min, temp->sec);
        }
        // For each pixel in scan ...
        for (size_t k = 0; k < npixels; k++) {

            // Get number of good prism wavelengths for each pixel
            n = 0;
            j = 0;
            for (j = 0; j < nbands_av; j++) {
                if (temp->Lt[j * npixels + k] > 0) n++;
            }

            // Bail if no good prism Lt values
            if (n == 0) continue;


            // Generate x & y arrays for spline
            // x - prism wavelengths
            // y - prism Lt values
            iwave = 0;
            for (j = 0; j < nbands_av; j++) {
                ociaWavelength[iwave] = temp->wave[j];
                Lt[iwave * npixels + k] = temp->Lt[k * nbands_av + j];
                //              Lt[iwave*npixels + k] = temp->Lt[j*npixels+k];
                fwhm[iwave] = temp->fwhm[j];
                iwave++;
            }
            //      ociaWavelength[iwave] = 936.436901;
            //      fwhm[iwave] = 32.4702967;
            //      Lt[iwave*npixels + k] = (11.9269624*temp->Lt[222*npixels+k] +
            //                               10.78153  *temp->Lt[223*npixels+k] +
            //                               9.76180425*temp->Lt[224*npixels+k])/fwhm[iwave];
            //      iwave++;
            //
            //      ociaWavelength[iwave] = 1071.049452;
            //      fwhm[iwave] = 20.3179033;
            //      Lt[iwave*npixels + k] = (10.1448143*temp->Lt[271*npixels+k] +
            //                               10.173089 *temp->Lt[272*npixels+k])/fwhm[iwave];
            //      iwave++;
            //
            //      Lt[iwave*npixels + k] = temp->Lt[273*npixels+k];
            //      fwhm[iwave] = temp->fwhm[273];
            //      ociaWavelength[iwave] = temp->wave[273];
            //      iwave++;
            //
            //      ociaWavelength[iwave] = 1085.272670;
            //      fwhm[iwave] = 53.6097915;
            //      Lt[iwave*npixels + k] = (10.7227049*temp->Lt[274*npixels+k] +
            //                               10.7221941*temp->Lt[275*npixels+k] +
            //                               10.7218208*temp->Lt[276*npixels+k] +
            //                               10.721585 *temp->Lt[277*npixels+k] +
            //                               10.7214867*temp->Lt[278*npixels+k])/fwhm[iwave];
            //      iwave++;
            //      //iwave = interpLt(133, npixels, k, nbands_av, orcaWavelength[iwave], iwave, temp, acc,Lt);
            //

        } // pixel loop

        if (NBANDS != iwave) {
            printf("Oops, something's wrong iwave (%d) != nbands (%d) \n", (int) iwave, NBANDS);
            exit(-1);
        }
        // Write Lt's
        for (size_t iwave = 0; iwave < nbands_ocia; iwave++) {
            int varID;

            //      printf("Lambda(%d) = %d\n",iwave+1,(int)(orcaWavelength[iwave]+.5));
            status = nc_inq_varid(grpID, "Lt", &varID);
            check_err(status, __LINE__, __FILE__);

            size_t start[3] = {iwave, iscan - sline, 0};
            size_t count[3] = {1, 1, npixels};
            status = nc_put_vara_float(grpID, varID, start, count,
                    &Lt[iwave * npixels]);
            check_err(status, __LINE__, __FILE__);
        }

        utc[iscan] = temp->scantime;
        //if (utc[iscan]>0) printf("utc=%f\n",utc[iscan]);

        grpID = ociafile.getGid("navigation_data");

        startOCIP[0] = iscan - sline;

        for (i = 0; i < 6; i++) {
            //          cout << "Copying " << navfieldsPRISM[i] << " to " <<
            //                  navfieldsOCIP[i] << endl;

            status = nc_inq_varid(grpID, navfieldsOCIP[i], &varID);
            check_err(status, __LINE__, __FILE__);

            switch (i) {
            case 0:
                status = nc_put_vara_float(grpID, varID, startOCIP, countOCIP, &temp->sena[0]);
                check_err(status, __LINE__, __FILE__);
                break;
            case 1:
                status = nc_put_vara_float(grpID, varID, startOCIP, countOCIP, &temp->senz[0]);
                check_err(status, __LINE__, __FILE__);
                break;
            case 2:
                status = nc_put_vara_float(grpID, varID, startOCIP, countOCIP, &temp->sola[0]);
                check_err(status, __LINE__, __FILE__);
                break;
            case 3:
                status = nc_put_vara_float(grpID, varID, startOCIP, countOCIP, &temp->solz[0]);
                check_err(status, __LINE__, __FILE__);
                break;
            case 4:
                status = nc_put_vara_double(grpID, varID, startOCIP, countOCIP, &temp->lat[0]);
                check_err(status, __LINE__, __FILE__);
                break;
            case 5:
                status = nc_put_vara_double(grpID, varID, startOCIP, countOCIP, &temp->lon[0]);
                check_err(status, __LINE__, __FILE__);
                break;
            default:
                break;
            }
        }
        if (temp->scantime > 0) {
            scantime_last = temp->scantime;
            if (!first_good_scantime) {
                first_good_scantime = 1;
                scantime_first = temp->scantime;
            }
        }

    } // scan loop

    if (scantime_first < scantime_last) {
        strncpy(isodatetime_first, unix2isodate(scantime_first, 'G'), 30);
        strncpy(isodatetime_last, unix2isodate(scantime_last, 'G'), 30);
    } else {
        strncpy(isodatetime_last, unix2isodate(scantime_first, 'G'), 30);
        strncpy(isodatetime_first, unix2isodate(scantime_last, 'G'), 30);
    }
    status = nc_put_att_text(ncid, NC_GLOBAL, "time_coverage_start", strlen(isodatetime_first) + 1, isodatetime_first);
    if (status != NC_NOERR) {
        printf("-E- %s %d: %s for %s\n",
                __FILE__, __LINE__, nc_strerror(status), "time_coverage_start");
        exit(1);
    }
    status = nc_put_att_text(ncid, NC_GLOBAL, "time_coverage_end", strlen(isodatetime_last) + 1, isodatetime_last);
    if (status != NC_NOERR) {
        printf("-E- %s %d: %s for %s\n",
                __FILE__, __LINE__, nc_strerror(status), "time_coverage_end");
        exit(1);
    }

    gsl_interp_accel_free(acc);

    size_t startw[2] = {0, 0};
    size_t countw[2] = {nbands_ocia, 0};

    grpID = ociafile.getGid("sensor_band_parameters");

    status = nc_inq_varid(grpID, "wavelength", &varID);

    status = nc_put_vara_float(grpID, varID, startw, countw, ociaWavelength);
    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(grpID, "fwhm", &varID);

    status = nc_put_vara_float(grpID, varID, startw, countw, fwhm);
    check_err(status, __LINE__, __FILE__);

    size_t starts[2] = {0, 0};
    size_t counts[2] = {nscans, 0};

    const char *nscanfieldsPRISM[] = {"utc", "alt"};
    const char *nscanfieldsOCIP[] = {"scan_start_time", "altitude"};
    short* alt = (short*) (calloc(nscans, sizeof (short)));

    for (i = 0; i < nscans; i++) {
        //utc[i] = temp->scantime;
        alt[i] = temp->alt;
    }
    grpID = ociafile.getGid("scan_line_attributes");
    countOCIP[1] = 0;
    for (i = 0; i < 2; i++) {
        cout << "Copying " << nscanfieldsPRISM[i] << " to " <<
                nscanfieldsOCIP[i] << endl;

        status = nc_inq_varid(grpID, nscanfieldsOCIP[i], &varID);
        check_err(status, __LINE__, __FILE__);

        switch (i) {
        case 0:
            status = nc_put_vara_double(grpID, varID, starts, counts, utc);
            check_err(status, __LINE__, __FILE__);
            break;
        case 1:
            status = nc_put_vara_short(grpID, varID, starts, counts, alt);
            check_err(status, __LINE__, __FILE__);
            break;
        default:
            break;
        }
    }

    //  outfile.close();

    status = close_prism(temp);

    return 0;
}
