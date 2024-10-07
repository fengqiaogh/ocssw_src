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
#include "aviris.h"
#include <timeutils.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort_double.h>

#define SENSOR "ocia"
#define CDLFILE "OCI_Level-1B_Data_Structure.cdl"
#define NBANDS 61 // number of OCI bands
#define FILL -32767

extern "C" void handle_error(int status);

using namespace std;


#define VERSION "1.00"

size_t interpLt(size_t startJ, size_t npixels, size_t k, size_t nbMax,
        float ociaWavelength, size_t iwave, aviris4ocia_t* data,
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
            Lt[iwave * npixels + k] = FILL;
        }
        iwave++;
        gsl_spline_free(spline);
    }
    free(x);
    free(y);
    return iwave;
}

//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     FutureTech     12/18/14 0.01  Original development
//  Joel Gales     FutureTech     01/12/15 0.10  Add support for writing
//                                               interpolated ORCA Lt values
//                                               and navigation
//  Rick Healy     SAIC           06/06/16 1.00  Integrated L2gen L1 reader for Aviris

int main(int argc, char* argv[]) {
    int status;
    size_t sline = 0, eline = 0;
    char filename[FILENAME_MAX], hdrfile[FILENAME_MAX], imgfile[FILENAME_MAX],
            navfile[FILENAME_MAX], gainfile[FILENAME_MAX];
    ;
    char *filedir;
    double scantime_first, scantime_last;
    static int first_good_scantime = 0;
    cout << "avirisbil2oci " << VERSION << " ("
            << __DATE__ << " " << __TIME__ << ")" << endl;

    if (argc < 3) {
        cout << endl <<
                "avirisbil2oci input_AVIRIS_file output_OCIA_file <sline (default=0)> <eline (default=last row) <img file(optional)> <nav file(optional)> <gain file(optional)>>"
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
    if (argc >= 6) {
        sscanf(argv[5], "%s", imgfile);
        if (eline < 0) eline = 0;
    }
    if (argc >= 7) {
        sscanf(argv[6], "%s", navfile);
        if (eline < 0) eline = 0;
    }
    if (argc >= 8) {
        sscanf(argv[7], "%s", gainfile);
        if (eline < 0) eline = 0;
    }
    //  aviris4orca_t *&data;

    if ((filedir = getenv("OCDATAROOT")) == NULL) {
        printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
        return (-1);
    }
    strcpy(filename, filedir);

    strcat(filename, "/");
    strcat(filename, SENSOR);
    strcat(filename, "/");

    strcat(filename, CDLFILE);

    if (!checkAvProcessFile(argv[1], hdrfile, imgfile, navfile, gainfile, FILENAME_MAX))
        strncpy(hdrfile, argv[1], FILENAME_MAX);

    aviris4ocia_t **data, *temp;
    float fwhm[NBANDS];
    data = (aviris4ocia_t **) malloc(sizeof (aviris4ocia_t *));
    temp = open_aviris(hdrfile, imgfile, navfile, gainfile, data);
    size_t nscans, ip, npixels = temp->npix;

    if (eline <= 0) eline = temp->nscans;

    if (eline == sline && sline > 0) sline--;
    nscans = eline - sline;

    int ncid, grpID;
    //char** dim_names, uint32_t* dim_size, size_t n_dims
    const char *dimNames[] = {"number_of_scans", "number_of_pixels", "number_of_bands"};
    size_t dimSize[] = {nscans, npixels, NBANDS};
    size_t n_dims = 3;

    // Create ORCA output file
    static ncdfFile ociafile;
    ociafile.cdlCreateDim(argv[2], filename, dimNames, dimSize, n_dims, nscans);

    //  nc_def_dim(orcafile.ncid, "number_of_lines", nscans, &dimid[0]);
    //  nc_def_dim(orcafile.ncid, "pixels_per_line", nscans, &dimid[1]);
    //
    // Copy AVIRIS nav fields to ORCA
    //
    ncid = ociafile.getNcid();

    size_t nbands_av = temp->numBands;
    size_t nbands_ocia = NBANDS;

    size_t n = 0, iwave = 0;
    istringstream istr;
    string str;
    char isodatetime_first[30], isodatetime_last[30];

    // Allocate space for AVIRIS Lt (short) values
    //  size_t start[3] = {0,0,0};
    //  size_t count[3] = {1,npixels,0};

    gsl_interp_accel *acc = gsl_interp_accel_alloc();

    float *Lt = (float *) calloc(npixels*nbands_ocia, sizeof (float));
    float *ociaWavelength = (float *) calloc(nbands_ocia, sizeof (float));
    double* utc = (double*) (calloc(nscans, sizeof (double)));
    float* alt = (float*) (calloc(nscans, sizeof (float)));
    size_t j;

    grpID = ociafile.getGid("observation_data");

    for (size_t iscan = sline; iscan < eline; iscan++) {


        // Read AVIRIS Lt's for this scan
        if (read_aviris(temp, iscan)) {

            cout << "Can't read file " << argv[1] << endl;
            exit(1);

        }
        for (ip = 0; ip < npixels; ip++) for (iwave = 0; iwave < nbands_ocia; iwave++) Lt[iwave * npixels + ip] = FILL;

        if ((iscan % 100) == 0) {
            cout << "Processing scan: " << iscan << " out of " << nscans << endl;
            //printf("year=%d doy=%d msec=%d month=%d day=%d hour=%d min=%d sec=%f\n ",temp->year, temp->doy, temp->msec,temp->month,temp->day,temp->hour,temp->min,temp->sec);
        }
        // For each pixel in scan ...
        for (size_t k = 0; k < npixels; k++) {

            // Get number of good AVIRIS wavelengths for each pixel
            n = 0;
            j = 0;
            for (j = 0; j < nbands_av; j++) {
                if (temp->Lt[j * npixels + k] > 0) n++;
            }

            // Bail if no good AVIRIS Lt values
            if (n == 0) continue;


            // Generate x & y arrays for spline
            // x - AVIRIS wavelengths
            // y - AVIRIS Lt values
            iwave = 0;
            for (j = 0; j < 31; j++) {
                if (temp->Lt[j * npixels + k] > 0) {
                    ociaWavelength[iwave] = temp->wave[j];
                    Lt[iwave * npixels + k] = temp->Lt[j * npixels + k];
                    fwhm[iwave] = temp->fwhm[j];
                }
                iwave++;
            }
            for (j = 33; j < 57; j++) {
                if (temp->Lt[j * npixels + k] > 0) {
                    ociaWavelength[iwave] = temp->wave[j];
                    Lt[iwave * npixels + k] = temp->Lt[j * npixels + k];
                    fwhm[iwave] = temp->fwhm[j];
                }
                iwave++;
            }

            ociaWavelength[iwave] = 936.646;
            fwhm[iwave] = 32.4702967;
            Lt[iwave * npixels + k] = (11.9269624 * temp->Lt[60 * npixels + k] +
                    10.78153 * temp->Lt[61 * npixels + k] +
                    9.76180425 * temp->Lt[62 * npixels + k]) / fwhm[iwave];
            iwave++;

            ociaWavelength[iwave] = 1238.750;
            fwhm[iwave] = 20.3179033;
            Lt[iwave * npixels + k] = (10.1448143 * temp->Lt[92 * npixels + k] +
                    10.173089 * temp->Lt[93 * npixels + k]) / fwhm[iwave];
            iwave++;

            Lt[iwave * npixels + k] = temp->Lt[109 * npixels + k];
            fwhm[iwave] = temp->fwhm[109];
            ociaWavelength[iwave] = temp->wave[109];
            iwave++;

            ociaWavelength[iwave] = 1641.530;
            fwhm[iwave] = 53.6097915;
            Lt[iwave * npixels + k] = (10.7227049 * temp->Lt[133 * npixels + k] +
                    10.7221941 * temp->Lt[134 * npixels + k] +
                    10.7218208 * temp->Lt[135 * npixels + k] +
                    10.721585 * temp->Lt[136 * npixels + k] +
                    10.7214867 * temp->Lt[137 * npixels + k]) / fwhm[iwave];
            iwave++;
            //iwave = interpLt(133, npixels, k, nbands_av, orcaWavelength[iwave], iwave, temp, acc,Lt);
            ociaWavelength[iwave] = 2126.795;
            fwhm[iwave] = 54.1850639;
            Lt[iwave * npixels + k] = (10.8702638 * temp->Lt[184 * npixels + k] +
                    10.8539551 * temp->Lt[185 * npixels + k] +
                    10.8373296 * temp->Lt[186 * npixels + k] +
                    10.8203873 * temp->Lt[187 * npixels + k] +
                    10.8031281 * temp->Lt[188 * npixels + k]) / fwhm[iwave];
            iwave++;
            //iwave = interpLt(184, npixels, k, nbands_av, orcaWavelength[iwave], iwave, temp, acc,Lt);
            ociaWavelength[iwave] = 2246.658;
            fwhm[iwave] = 53.0639591;
            Lt[iwave * npixels + k] = (10.6536474 * temp->Lt[196 * npixels + k] +
                    10.6335365 * temp->Lt[197 * npixels + k] +
                    10.6131087 * temp->Lt[198 * npixels + k] +
                    10.592364 * temp->Lt[199 * npixels + k] +
                    10.5713025 * temp->Lt[200 * npixels + k]) / fwhm[iwave];
            iwave++;
            //iwave = interpLt(196, npixels, k, nbands_av, orcaWavelength[iwave], iwave, temp, acc,Lt);


        } // pixel loop

        if (NBANDS != iwave) {
            printf("Oops, something's wrong iwave (%d) != nbands (%d) \n", (int) iwave, NBANDS);
            exit(-1);
        }
        // Write Lt's
        //    for (size_t iwave=0; iwave<nbands_ocia; iwave++ ) {
        int varID;

        //      printf("Lambda(%d) = %d\n",iwave+1,(int)(orcaWavelength[iwave]+.5));
        status = nc_inq_varid(grpID, "Lt", &varID);
        check_err(status, __LINE__, __FILE__);

        size_t start[3] = {0, iscan, 0};
        size_t count[3] = {nbands_ocia, 1, npixels};
        status = nc_put_vara_float(grpID, varID, start, count,
                &Lt[0]);
        //      &Lt[iwave*npixels]);
        check_err(status, __LINE__, __FILE__);
        //    }
        alt[iscan] = temp->alt[iscan];
        utc[iscan] = temp->scantime;
        if (temp->scantime > 0) {
            scantime_last = temp->scantime;
            if (!first_good_scantime) {
                first_good_scantime = 1;
                scantime_first = temp->scantime;
                //    		printf("utc=%f\n",utc[iscan]);
                //        	strncpy(isodatetime_first,unix2isodate(scantime_first, 'G'),30);
                //        	printf("%s\n",isodatetime_first);
            }
            //    	strncpy(isodatetime_last,unix2isodate(scantime_last, 'G'),30);
            //    	printf("Last=%s\n",isodatetime_last);
        }


    } // scan loop
    int varID;

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

    const char *navfieldsAVIRIS[] = {"to-sensor_azimuth", "to-sensor_zenith", "to-sun_azimuth", "to-sun_zenith", "lat", "lon"};
    const char *navfieldsOCIA[] = {"sensor_azimuth", "sensor_zenith", "solar_azimuth", "solar_zenith", "latitude", "longitude"};

    grpID = ociafile.getGid("geolocation_data");

    size_t startOCIA[2] = {0, 0};
    size_t countOCIA[2] = {nscans, npixels};
    size_t i;

    for (i = 0; i < 6; i++) {
        cout << "Copying " << navfieldsAVIRIS[i] << " to " <<
                navfieldsOCIA[i] << endl;

        status = nc_inq_varid(grpID, navfieldsOCIA[i], &varID);
        check_err(status, __LINE__, __FILE__);

        switch (i) {
        case 0:
            status = nc_put_vara_float(grpID, varID, startOCIA, countOCIA, &temp->sena[sline * npixels]);
            check_err(status, __LINE__, __FILE__);
            break;
        case 1:
            status = nc_put_vara_float(grpID, varID, startOCIA, countOCIA, &temp->senz[sline * npixels]);
            check_err(status, __LINE__, __FILE__);
            break;
        case 2:
            status = nc_put_vara_float(grpID, varID, startOCIA, countOCIA, &temp->sola[sline * npixels]);
            check_err(status, __LINE__, __FILE__);
            break;
        case 3:
            status = nc_put_vara_float(grpID, varID, startOCIA, countOCIA, &temp->solz[sline * npixels]);
            check_err(status, __LINE__, __FILE__);
            break;
        case 4:
            status = nc_put_vara_double(grpID, varID, startOCIA, countOCIA, &temp->lat[sline * npixels]);
            check_err(status, __LINE__, __FILE__);
            break;
        case 5:
            status = nc_put_vara_double(grpID, varID, startOCIA, countOCIA, &temp->lon[sline * npixels]);
            check_err(status, __LINE__, __FILE__);
            break;
        default:
            break;
        }
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

    const char *nscanfieldsAVIRIS[] = {"utc", "alt"};
    const char *nscanfieldsOCIA[] = {"scan_start_time","altitude"};

    for (i = 0; i < nscans; i++) {
        //utc[i] = temp->scantime;
        //alt[i] = temp->alt[sline + i]*1000;
    }
    grpID = ociafile.getGid("scan_line_attributes");
    countOCIA[1] = 0;
    for (i = 0; i < 2; i++) {
        cout << "Copying " << nscanfieldsAVIRIS[i] << " to " <<
                nscanfieldsOCIA[i] << endl;

        status = nc_inq_varid(grpID, nscanfieldsOCIA[i], &varID);
        check_err(status, __LINE__, __FILE__);

        switch (i) {
        case 0:
            status = nc_put_vara_double(grpID, varID, starts, counts, utc);
            check_err(status, __LINE__, __FILE__);
            break;
        case 1:
            status = nc_put_vara_float(grpID, varID, starts, counts, alt);
            check_err(status, __LINE__, __FILE__);
            break;
        default:
            break;
        }
    }

    //  outfile.close();

    status = close_aviris(temp);

    return 0;
}
