#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <list>
#include "hdf4_bin.h"

/*
 typedef struct binListStruct {
    int32 bin_num;
    int16 nobs;
    int16 nscenes;
    int16 time_rec;
    float32 weights;
    uchar8 sel_cat;
    int32 flags_set;
    float32 lat;
    float32 lon;
  } binListStructure;

 */

using namespace std;

int main(int argc, char* argv[]) {
    char ptime[17];
    char buf[1024];

    get_time(ptime);

    // Open input binfile
    static Hdf::hdf4_bin input_binfile;
    input_binfile.open(argv[1]); //, H5F_ACC_RDONLY);

    // Get sigma value
    float32 sigma = atof(argv[2]);

    // Alternative product array
    char **prod_array;
    //  int nprod = input_binfile.query( &prod_array);
    // for (size_t i=0; i<nprod; i++)
    //cout << i << " " << prod_array[i] << endl;

    // Get number of input file bins (data_records)
    int32 n_input_bins = input_binfile.n_data_records;

    // Get total number of data elements for two nLw fields & allocate buffer
    int n_data_elem = input_binfile.read(argv[3]);
    float32 *inData = (float32 *) calloc(n_data_elem, sizeof (float32));
    float32 *inVar = (float32 *) calloc(n_data_elem, sizeof (float32));

    // Allocate BinList structure buffer
    Hdf::binListStruct *inBinList =
            (Hdf::binListStruct *) calloc(n_input_bins, sizeof (Hdf::binListStruct));

    // Read BinList & Data
    int32 nread = input_binfile.read(inData, inVar, inBinList);

    // Allocate bins to copy, means, std dev
    int32 *binsSpeckled = (int32 *) calloc(n_input_bins, sizeof (int32));
    float32 *meansSpeckled = (float32 *) calloc(n_input_bins, sizeof (float32));
    float32 *sdevsSpeckled = (float32 *) calloc(n_input_bins, sizeof (float32));

    // Open 9km near bin datafile
    FILE *fp;
    int32 n_index, n_near_bins, *index, *near_bins;

    char *tmp_str;
    if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
        printf("OCDATAROOT environment variable is not defined.\n");
        return (1);
    }
    strcpy(buf, tmp_str);
    strcat(buf, "/common/9km_near.dat");

    fp = fopen(buf, "rb");
    if (fp == NULL) {
        cout << buf << " not found." << endl;
        exit(1);
    }

    fread(&n_index, sizeof (int32), 1, fp);
    fread(&n_near_bins, sizeof (int32), 1, fp);
    index = (int32 *) calloc(n_index, sizeof (int32));
    near_bins = (int32 *) calloc(n_near_bins, sizeof (int32));

    fread(index, sizeof (int32), n_index, fp);
    fread(near_bins, sizeof (int32), n_near_bins, fp);
    fclose(fp);

    // Compute bin_arr (bin # to entry in inBinList[].bin_num)
    int32 totbins = input_binfile.totbins;
    int32 *bin_arr = (int32 *) calloc(totbins, sizeof (int32));
    for (size_t i = 0; i < totbins; i++) bin_arr[i] = -1;
    for (size_t i = 0; i < n_input_bins; i++) {
        int32 bin_zerobased = inBinList[i].bin_num - 1;
        bin_arr[bin_zerobased] = i;
    }

    // Main Loop
    float32 cen_val, cen_var;
    float32 neigh_val[32];
    float32 wgt_val[32];
    int32 n_speckle_bins = 0;
    for (size_t i = 0; i < n_input_bins; i++) {
        int32 bin_zerobased = inBinList[i].bin_num - 1;
        int32 n_close_bins = index[bin_zerobased + 1] - index[bin_zerobased];
        int32 idx = index[bin_zerobased];

        int32 k = 0;
        for (size_t j = 0; j < n_close_bins; j++) {
            int32 ibin = near_bins[idx + j]; // 1-based
            if (bin_arr[ibin - 1] != -1) {
                if (ibin == inBinList[i].bin_num) {
                    cen_val = inData[i];
                    cen_var = inVar[i];
                } else {
                    neigh_val[k] = inData[bin_arr[ibin - 1]];
                    wgt_val[k] = inBinList[bin_arr[ibin - 1]].weights;

                    k++;
                }
            } // if close bin found
        } // j-loop (n_close_bins)

        bool keep = true;
        float32 mean, stdev;
        if (k >= 2) {
            float32 sum = 0;
            float32 sum2 = 0;
            float32 sumw = 0;
            float32 sumw2 = 0;
            for (size_t j = 0; j < k; j++) {
                sum += neigh_val[j] * wgt_val[j];
                sum2 += neigh_val[j] * neigh_val[j] * wgt_val[j];
                sumw += wgt_val[j];
                sumw2 += wgt_val[j] * wgt_val[j];
            }
            mean = sum / sumw;
            stdev = sqrt((sum2 * sumw - sum * sum) / (sumw * sumw - sumw2));
            if (fabs(cen_val - mean) / stdev > sigma) {
                keep = false;
                //	printf("Speckle Bin: %7d %3d %3d %8.5f %8.5f %8.5f %8.5f %7.2f\n", 
                //      inBinList[i].bin_num, inBinList[i].nobs, inBinList[i].nscenes, 
                //      cen_val, sqrt(cen_var), mean, stdev, 
                //      fabs(cen_val - mean) / stdev);

                // This is the entry number within the BinList (0-based),
                // NOT the bin number itself which is 1-based.
                binsSpeckled[n_speckle_bins] = i;
                meansSpeckled[n_speckle_bins] = mean;
                sdevsSpeckled[n_speckle_bins++] = stdev;
            }
        }

    } // i-loop (n_input_bins)

    cout << n_input_bins << " " << n_speckle_bins << endl;


    // Open speckled bins output file (binary)
    FILE *fp_out;
    fp_out = fopen(argv[4], "wb");
    fwrite(&n_speckle_bins, sizeof (int32), 1, fp_out);
    for (size_t i = 0; i < n_speckle_bins; i++) {
        int32 binnum = inBinList[binsSpeckled[i]].bin_num;
        fwrite(&binnum, sizeof (int32), 1, fp_out);
        fwrite(&meansSpeckled[i], sizeof (float32), 1, fp_out);
        fwrite(&sdevsSpeckled[i], sizeof (float32), 1, fp_out);
    }

    free(inData);
    free(inVar);
    free(inBinList);
    free(binsSpeckled);
    free(meansSpeckled);
    free(sdevsSpeckled);
    cout << "Normal Completion" << endl << endl;

    return 0;
}

