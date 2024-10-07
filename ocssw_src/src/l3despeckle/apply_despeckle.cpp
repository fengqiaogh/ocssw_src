#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <list>
#include "hdf4_bin.h"
#include <timeutils.h>

using namespace std;

int main(int argc, char* argv[]) {
    char ptime[17];
    char buf[1024];

    get_time(ptime);

    // Open input binfile
    static Hdf::hdf4_bin input_binfile;
    input_binfile.open(argv[1]); //, H5F_ACC_RDONLY);


    int32 n_input_bins = input_binfile.read(argv[3]);

    // Allocate BinList structure buffer
    Hdf::binListStruct *inBinList =
            (Hdf::binListStruct *) calloc(n_input_bins, sizeof (Hdf::binListStruct));

    // Allocate data buffer
    float32 *inData = (float32 *) calloc(n_input_bins, sizeof (float32));

    // Read BinList & Data
    n_input_bins = input_binfile.read(inData, inBinList);

    // Allocate bins to copy
    int32 *binsToCopy = (int32 *) calloc(n_input_bins, sizeof (float32));

    // Open speckled bins datafile
    FILE *fp;
    int32 n_speckle_bins, *speckle_bins;
    float32 *stdev, *mean;
    fp = fopen(argv[2], "rb");
    fread(&n_speckle_bins, sizeof (int32), 1, fp);
    speckle_bins = (int32 *) calloc(n_speckle_bins, sizeof (int32));
    mean = (float32 *) calloc(n_speckle_bins, sizeof (float32));
    stdev = (float32 *) calloc(n_speckle_bins, sizeof (float32));
    for (size_t i = 0; i < n_speckle_bins; i++) {
        fread(&speckle_bins[i], sizeof (int32), 1, fp);
        fread(&mean[i], sizeof (float32), 1, fp);
        fread(&stdev[i], sizeof (float32), 1, fp);
    }



    // Determine bins to copy
    int32 lasti = 0;
    int32 n_bins_removed = 0;
    for (size_t j = 0; j < n_speckle_bins; j++) {
        //    if ( (j % 1000) == 0) printf("%d\n", j);
        for (size_t i = lasti; i < n_input_bins; i++) {
            if (inBinList[i].bin_num == speckle_bins[j]) {
                if ((fabs(inData[i] - mean[j]) / stdev[j]) > 3.0) {
                    inBinList[i].bin_num = -1;
                    n_bins_removed++;
                }
                lasti = i + 1;
                break;
            }
            if (inBinList[i].bin_num > speckle_bins[j]) {
                break;
            }
        }
    }

    int32 n_output_bins = 0;
    for (size_t i = 0; i < n_input_bins; i++) {
        if (inBinList[i].bin_num != -1) {
            binsToCopy[n_output_bins++] = i;
        }
    }

    static Hdf::hdf4_bin output_binfile;

    // Copy metadata
    output_binfile.meta_l3b = input_binfile.meta_l3b;
    strcpy(output_binfile.meta_l3b.ptime, ptime);
    strcat(output_binfile.meta_l3b.product_name, "_DS");
    strcat(output_binfile.meta_l3b.soft_name, "|apply_despeckle");
    strcat(output_binfile.meta_l3b.soft_ver, "|0.1");
    strcat(output_binfile.meta_l3b.proc_con, "|apply_despeckle ");
    strcat(output_binfile.meta_l3b.proc_con, argv[1]);
    strcat(output_binfile.meta_l3b.proc_con, " ");
    strcat(output_binfile.meta_l3b.proc_con, argv[2]);
    strcat(output_binfile.meta_l3b.proc_con, " ");
    strcat(output_binfile.meta_l3b.proc_con, argv[3]);
    strcat(output_binfile.meta_l3b.proc_con, " ");
    strcat(output_binfile.meta_l3b.proc_con, argv[4]);

    // Create output bin file
    strcpy(buf, argv[4]);
    output_binfile.create(buf, input_binfile.nrows);

    // Determine size of input file product list
    char *fullprodlist = (char *) malloc(input_binfile.query());

    // Read input file product list
    input_binfile.query(fullprodlist);

    // Copy specified bins
    output_binfile.copy(fullprodlist, n_output_bins, binsToCopy,
            inBinList, &input_binfile);

    // Close bin files
    output_binfile.close();
    input_binfile.close();

    free(fullprodlist);
    free(inBinList);
    free(binsToCopy);
    free(speckle_bins);
    free(mean);
    free(stdev);

    printf("--Number of input bins:   %8d\n", n_input_bins);
    printf("--Number of bins removed: %8d\n", n_bins_removed);
    printf("--Percent removed: %7.2f\n", n_bins_removed * 100. / n_input_bins);

    cout << "Normal Completion" << endl << endl;

    return 0;
}


