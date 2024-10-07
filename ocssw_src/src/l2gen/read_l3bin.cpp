// =====================================================================
// read_l3bin - level-3 bin reader for l2gen
// B. Franz, NASA/OBPG, March 2013            
// =====================================================================

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <list>
#include "hdf_bin.h"
#include <strings.h>

#include "read_l3bin.h"
#include "l12_proto.h"

using namespace std;

/**
 * read_l3bin reads "target" information for vicarious calibration
 * from a L3 bin file (nLw_vvv or Rrs_vvv, chl, tau)
 *
 * @param file - input target L3 bin file
 * @param l3bin - structure to populate
 * @return
 */
extern "C" int32_t read_l3bin(char *file, l3binstr *l3bin, int32_t nbands) {
    int32_t i, j, k, iprod, band;

    // Open input binfile and to get bin dimension
    cout << endl << "Opening " << file << " for reading" << endl;
    Hdf::hdf_bin *binfile;
    binfile = Hdf::openBinObject(file);
    if (binfile == NULL) {
        fprintf(stderr, "-E- Unable to open %s for reading.\n", file);
        exit(EXIT_FAILURE);
    }
    int32_t nbins = binfile->n_data_records;
    int32_t nrows = binfile->nrows;
    l3bin->sensorID = sensorName2SensorId(binfile->meta_l3b.sensor_name); //binfile->meta_l3b.sensorID;
    char **prodlist;
    int32 nprods = binfile->query(&prodlist);

    // Allocate space for output and work space for i/o        
    l3bin->data = (float **) malloc(nbins * sizeof (float *));
    for (i = 0; i < nbins; i++) {
        if ((l3bin->data[i] = (float *) calloc(nbands, sizeof (float))) == NULL) {
            printf("-E- Cannot allocate memory to hold L3 data\n");
            exit(EXIT_FAILURE);
        }
    }
    l3bin->nobs = (int32_t *) calloc(nbins, sizeof (int32_t));
    l3bin->nscenes = (int32_t *) calloc(nbins, sizeof (int32_t));
    l3bin->chl = (float *) calloc(nbins, sizeof (float));
    l3bin->tau = (float *) calloc(nbins, sizeof (float));
    l3bin->bins = (int32_t *) calloc(nbins, sizeof (int32_t));
    l3bin->basebin = (int32_t *) calloc(nrows, sizeof (int32_t));
    l3bin->numbin = (int32_t *) calloc(nrows, sizeof (int32_t));
    l3bin->wavelengths = (float *) calloc(nbands, sizeof(float));

    float **inData;
    if ((l3bin->prods = (char **) malloc(nprods * sizeof (char *))) == NULL) {
        printf("-E- Error allocating memory to the L3 product list .\n");
        exit(EXIT_FAILURE);
    }
    if ((inData = (float **) malloc(nprods * sizeof (float *))) == NULL) {
        printf("-E- Error allocating memory to the L3 input data array.\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < nprods; i++) {
        if ((inData[i] = (float *) calloc(2 * MAXPIX, sizeof (float))) == NULL) {
            printf("-E- Error allocating memory to the L3 input data array.\n");
            exit(EXIT_FAILURE);
        }
        if ((l3bin->prods[i] = (char *) calloc(32, sizeof (char))) == NULL) {
            printf("-E- Error allocating memory to the L3 product list.\n");
            exit(EXIT_FAILURE);
        }
    }

    char buf[FILENAME_MAX];
    char prodstr[L1_PRODSTRLEN] = "";
    char *wavelengthStr;
    j = 0;
    l3bin->hasRrs = 0;
    l3bin->nwave = 0;
    for (i = 0; i < nprods; i++) {
        if ((strcmp(prodlist[i], "chlor_a") == 0) ||
                (strncmp(prodlist[i], "aot", 3) == 0) ||
                (strncmp(prodlist[i], "tau", 3) == 0) ||
                (strncmp(prodlist[i], "nLw", 3) == 0) ||
                (strncmp(prodlist[i], "Rrs", 3) == 0)) {
            strcpy(l3bin->prods[j++], prodlist[i]);
            l3bin->nprods++;
            if ((strncmp(prodlist[i], "nLw", 3) == 0) ||
                    (strncmp(prodlist[i], "Rrs", 3) == 0)) {
                strcpy(buf,prodlist[i]);
                wavelengthStr = strtok(buf, "_");
                wavelengthStr = strtok(NULL,"\0");
                l3bin->wavelengths[l3bin->nwave] = (float) atof(wavelengthStr);
                l3bin->nwave++;
                if (strncmp(prodlist[i], "Rrs", 3) == 0)
                    l3bin->hasRrs = 1;
            }
            strcat(prodstr, prodlist[i]);
            if (i < nprods - 1)
                strcat(prodstr, ",");
        }
    }
    binfile->read(prodstr);

    // Save info to return struct
    l3bin->file = strdup(file);
    l3bin->nbins = nbins;
    l3bin->nrows = nrows;
    int bindx = 0;

    for (int krow = 0; krow < nrows; krow++) {
        int basebin = binfile->get_basebin(krow);
        int numbin = binfile->get_numbin(krow);
        l3bin->basebin[krow] = basebin;
        l3bin->numbin[krow] = numbin;
        binfile->readBinIndex(krow);

        int ext = binfile->get_ext();
        int beg = binfile->get_beg();
        // if the row has no filled bins, skip it
        if (beg == 0)
            continue;

        int status = binfile->readBinList(ext);
        // if something went horribly wrong...
        if (status == -1) {
            printf("Unable to read bin numbers...: %d\n", ext);
            exit(EXIT_FAILURE);
        }

        // fill the input data array with the necessary input data
        iprod = 0;

        for (int j = 0; j < binfile->nprod(); j++) {

            if (binfile->active_data_prod[j] == true) {
                binfile->get_prodname(j, buf);
                if (strcmp(prodlist[j], buf) == 0) {
                    binfile->readSums(inData[iprod], ext, j);
                    iprod++;
                }
            }
        }

        if (binfile->isHDF5 || binfile->isCDF4)
            binfile->setDataPtr(ext);

        // Apply the weights to get the mean from the sum
        for (k = 0; k < ext; k++) {
            float weight = binfile->get_weights(k);
            band = 0;
            for (iprod = 0; iprod < l3bin->nprods; iprod++) {
                if (strcmp(l3bin->prods[iprod], "chlor_a") == 0) {
                    l3bin->chl[bindx] = inData[iprod][2 * k] / weight;
                } else if ((strncmp(l3bin->prods[iprod], "tau", 3) == 0) ||
                        (strncmp(l3bin->prods[iprod], "aot", 3) == 0)) {
                    l3bin->tau[bindx] = inData[iprod][2 * k] / weight;
                } else if ((strncmp(l3bin->prods[iprod], "nLw", 3) == 0) ||
                        (strncmp(l3bin->prods[iprod], "Rrs", 3) == 0)) {
                    l3bin->data[bindx][band++] = inData[iprod][2 * k] / weight;
                }
            }
            l3bin->bins[bindx] = binfile->get_bin_num(k);
            l3bin->nobs[bindx] = binfile->get_nobs(k);
            l3bin->nscenes[bindx] = binfile->get_nscenes(k);
            bindx++;
        }

    }
    // Clean-up
    binfile->close();
    for (i = 0; i < nprods; i++)
        free(inData[i]);
    free(inData);

    return (nbins);
}

