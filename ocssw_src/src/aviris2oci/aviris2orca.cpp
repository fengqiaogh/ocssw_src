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
#include "dfutils.h"
#include "nc_gridutils.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort_double.h>


extern "C" void handle_error(int status);

using namespace std;


#define VERSION "0.10"

//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     FutureTech     12/18/14 0.01  Original development
//  Joel Gales     FutureTech     01/12/15 0.10  Add support for writing
//                                               interpolated ORCA Lt values
//                                               and navigation

int main(int argc, char* argv[]) {
    int status;

    cout << "aviris2orca " << VERSION << " ("
            << __DATE__ << " " << __TIME__ << ")" << endl;

    if (argc == 1) {
        cout << endl <<
                "aviris2orca input_AVIRIS_file output_ORCA_file ORCA_CDL_structure_file"
                << endl;
        return 0;
    }

    // Open AVIRIS file
    idDS aviris_id = openDS(argv[1]);

    size_t nscans, npixels;
    int dimid;
    int grpID;

    // Get number of pixels (per scan) in AVIRIS file
    status = nc_inq_dimid(aviris_id.fid, "x", &dimid);
    if (status != NC_NOERR) {
        cout << "Can't find 'x' dimension in: " << argv[1] << endl;
        exit(1);
    }
    status = nc_inq_dimlen(aviris_id.fid, dimid, &npixels);

    // Get number of scans in AVIRIS file
    status = nc_inq_dimid(aviris_id.fid, "y", &dimid);
    if (status != NC_NOERR) {
        cout << "Can't find 'y' dimension in: " << argv[1] << endl;
        exit(1);
    }
    status = nc_inq_dimlen(aviris_id.fid, dimid, &nscans);

    // Create ORCA output file
    static ncdfFile orcafile;
    orcafile.cdlCreate(argv[2], argv[3], nscans);

    //
    // Copy AVIRIS nav fields to ORCA
    //
    char *navfieldsAVIRIS[] = {"to-sensor_azimuth", "to-sensor_zenith",
        "to-sun_azimuth", "to-sun_zenith",
        "lat", "lon"};
    char *navfieldsORCA[] = {"sena", "senz", "sola", "solz", "lat", "lon"};

    float *navbuffer = (float *) calloc(nscans*npixels, sizeof (float));
    grpID = orcafile.getGid("navigation_data");
    int varID;

    int startAVIRIS[2] = {0, 0};
    int countAVIRIS[2] = {nscans, npixels};
    size_t startORCA[2] = {0, 0};
    size_t countORCA[2] = {nscans, npixels};

    for (size_t i = 0; i < 6; i++) {
        cout << "Copying " << navfieldsAVIRIS[i] << " to " <<
                navfieldsORCA[i] << endl;

        PTB(readDS(aviris_id, navfieldsAVIRIS[i],
                startAVIRIS, NULL, countAVIRIS, navbuffer));

        status = nc_inq_varid(grpID, navfieldsORCA[i], &varID);
        check_err(status, __LINE__, __FILE__);

        status = nc_put_vara_float(grpID, varID, startORCA, countORCA, navbuffer);
        check_err(status, __LINE__, __FILE__);
    }

    // Compute attitude from path_length and cos(sensor zenith angle)
    float *plbuffer = (float *) calloc(nscans*npixels, sizeof (float));
    PTB(readDS(aviris_id, "path_length",
            startAVIRIS, NULL, countAVIRIS, plbuffer));
    PTB(readDS(aviris_id, "to-sensor_zenith",
            startAVIRIS, NULL, countAVIRIS, navbuffer));
    float deg2rad = 3.1415927 / 180;
    for (size_t j = 0; j < nscans * npixels; j++) {
        if (navbuffer[j] == -9999)
            plbuffer[j] = -9999;
        else
            plbuffer[j] *= cos(deg2rad * navbuffer[j]);
    }
    status = nc_inq_varid(grpID, "altitude", &varID);
    check_err(status, __LINE__, __FILE__);

    status = nc_put_vara_float(grpID, varID, startORCA, countORCA, plbuffer);
    check_err(status, __LINE__, __FILE__);


    // Get number of datasets in AVIRIS file
    int n_datasets;
    status = nc_inq(aviris_id.fid, NULL, &n_datasets, NULL, NULL);

    // Allocate space for AVIRIS Lt dataset IDs, wavelengths, and scale factors
    int *varLtid = (int *) calloc(n_datasets, sizeof (int));
    float *avirisWavelength = (float *) calloc(n_datasets, sizeof (float));
    double *scale_factor = (double *) calloc(n_datasets, sizeof (double));

    size_t n = 0;
    istringstream istr;
    string str;
    int att_type;
    size_t typeSize;
    char nambuf[NC_MAX_NAME];

    for (size_t i = 0; i < n_datasets; i++) {
        nc_inq_varname(aviris_id.fid, i, nambuf);

        // If Lt dataset ...
        if (strncmp(nambuf, "Lt_", 3) == 0) {

            // Store dataset ID
            varLtid[n] = i;

            printf("Nambuf=%s %d\n", nambuf, varLtid[n]);
            // Store AVIRIS wavelength
            str.assign(nambuf);
            str.replace(str.find_first_of('_', 3), 1, ".");
            istr.clear();
            istr.str(str.substr(3, string::npos));
            istr >> avirisWavelength[n];

            // status = nc_inq_atttype( aviris_id.fid, i, "scale_factor", &att_type);
            //status = nc_inq_type( aviris_id.fid, att_type, nambuf, &typeSize);

            // Store scale factors
            status = nc_get_att_double(aviris_id.fid, i, "scale_factor",
                    &scale_factor[n]);
            n++;
        }
    }
    int n_Lt_datasets = n;

    // Allocate space for AVIRIS Lt (short) values
    int start[2] = {0, 0};
    int count[2] = {1, npixels};
    int16_t *sval = (int16_t *) calloc(npixels*n_Lt_datasets, sizeof (int16_t));

    gsl_interp_accel *acc = gsl_interp_accel_alloc();

    size_t n_refl_bands = 91; // Wavelength <= 800nm
    size_t n_ir_bands = 6; // Wavelength >= 940nm

    int orcaWavelength[] ={350, 355, 360, 365, 370, 375, 380, 385, 390, 395, 400, 405, 410, 415,
        420, 425, 430, 435, 440, 445, 450, 455, 460, 465, 470, 475, 480, 485,
        490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555,
        560, 565, 570, 575, 580, 585, 590, 595, 600, 605, 610, 615, 620, 625,
        630, 635, 640, 645, 650, 655, 660, 665, 670, 675, 680, 685, 690, 695,
        700, 705, 710, 715, 720, 725, 730, 735, 740, 745, 750, 755, 760, 765,
        770, 775, 780, 785, 790, 795, 800, 940, 1240, 1378, 1640, 2130, 2250, 0};

    float *visnir = (float *) calloc(npixels*n_refl_bands, sizeof (float));
    float *ir = (float *) calloc(npixels*n_ir_bands, sizeof (float));

    grpID = orcafile.getGid("earth_view_data");

    for (size_t iscan = 0; iscan < nscans; iscan++) {

        if ((iscan % 100) == 0)
            cout << "Processing scan: " << iscan << " out of " << nscans << endl;

        for (size_t j = 0; j < npixels * n_refl_bands; j++) visnir[j] = -999.0;

        // Read AVIRIS Lt datasets for this scan
        start[0] = iscan;
        for (size_t j = 0; j < n_Lt_datasets; j++) {
            printf("max=%d j=%d varid=%d\n", n_Lt_datasets, j, varLtid[j]);
            nc_inq_varname(aviris_id.fid, varLtid[j], nambuf);
            PTB(readDS(aviris_id, nambuf, start, NULL, count, &sval[j * npixels]));
        }

        // For each pixel in scan ...
        for (size_t k = 0; k < npixels; k++) {

            // Get number of good AVIRIS wavelengths for each pixel
            n = 0;
            for (size_t j = 0; j < n_Lt_datasets; j++) {
                if (sval[j * npixels + k] > 0) n++;
            }

            // Bail if no good AVIRIS Lt values
            if (n == 0) continue;

            // Allocate space for spline
            gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, n);
            double *x = (double *) calloc(n, sizeof (double));
            double *y = (double *) calloc(n, sizeof (double));

            // Generate x & y arrays for spline
            // x - AVIRIS wavelengths
            // y - AVIRIS Lt values
            n = 0;
            for (size_t j = 0; j < n_Lt_datasets; j++) {
                if (sval[j * npixels + k] > 0) {
                    x[n] = avirisWavelength[j];
                    y[n] = sval[j * npixels + k] * scale_factor[j];
                    n++;
                }
            }

            // Sort wavelenghts and corresponding Lt values
            gsl_sort2(x, 1, y, 1, n);

            // Initiate spline
            gsl_spline_init(spline, x, y, n);

            // Generate interpolated values
            size_t iwave = 0;
            while (orcaWavelength[iwave] != 0) {
                if (orcaWavelength[iwave] > x[0] &&
                        orcaWavelength[iwave] < x[n - 1]) {
                    double yi = gsl_spline_eval(spline, orcaWavelength[iwave], acc);
                    if (orcaWavelength[iwave] <= 800)
                        visnir[iwave * npixels + k] = yi;
                    else
                        ir[(iwave - n_refl_bands) * npixels + k] = yi;
                    //          cout << setw(6);
                    //cout << orcaWavelength[iwave] << " ";
                    //cout << setw(10) << setprecision(5) << fixed << yi << endl;
                }
                iwave++;
            }

            free(x);
            free(y);
            gsl_spline_free(spline);

        } // pixel loop

        // Write vis & nir bands to Lt_visnir
        for (size_t iwave = 0; iwave < n_refl_bands; iwave++) {
            int varID;
            status = nc_inq_varid(grpID, "Lt_visnir", &varID);
            check_err(status, __LINE__, __FILE__);

            size_t start[3] = {iwave, iscan, 0};
            size_t count[3] = {1, 1, npixels};
            status = nc_put_vara_float(grpID, varID, start, count,
                    &visnir[iwave * npixels]);
            check_err(status, __LINE__, __FILE__);
        }

        // Write irbands to Lt_xxxx (xxxx=940,1240,1378,1640,2130,2250)
        char *irLt[] = {"Lt_940", "Lt_1240", "Lt_1378",
            "Lt_1640", "Lt_2130", "Lt_2250"};

        for (size_t iwave = 0; iwave < n_ir_bands; iwave++) {
            int varID;
            status = nc_inq_varid(grpID, irLt[iwave], &varID);
            check_err(status, __LINE__, __FILE__);

            size_t start[2] = {iscan, 0};
            size_t count[2] = {1, npixels};
            status = nc_put_vara_float(grpID, varID, start, count,
                    &ir[iwave * npixels]);
            check_err(status, __LINE__, __FILE__);
        }

    } // scan loop

    gsl_interp_accel_free(acc);

    //  outfile.close();

    status = endDS(aviris_id);

    free(visnir);
    free(ir);

    free(varLtid);
    free(avirisWavelength);
    free(scale_factor);
    free(sval);

    return 0;
}
