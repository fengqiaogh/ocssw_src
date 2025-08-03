/**
 * @file lib3merge.cpp
 * @brief This tool merges multiple level-3 binned files into a single level-3 binned file.
 * @authors J. Gales, S. Bailey, R. Healy
 **/
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <libgen.h>
#include <sys/types.h>

#include "netcdf.h"

#include "l3bin_input.h"
#include <timeutils.h>
#include <genutils.h>

#include "hdf_bin.h"

#include <hdf.h>
#include <mfhdf.h>

#define MAXNFILES 128

#define BYTE unsigned char

#define BINCHECK -1

#define VERSION "2.1"

static instr input;

/*

 Revision 2.1  12/5/2016
 Updated with the CLO library
 R. Healy

 Revision 2.0  08/29/2016
 Converted to 64 bit integer version
 R. Healy

 Revision 1.03 06/19/2014
 Added support for netCDF4 files
 S. Bailey

 Revision 1.02 08/16/12
 Add "-u" switch (Union of input bins)
 J. Gales

 Revision 1.01 03/21/12
 Replace last_prod with prod_offset[]
 J. Gales

 Revision 1.00 03/18/12
 Initial Version
 J. Gales
 */

int main(int argc, char **argv) {
    intn i;

    int32_t irow;
    int32_t kbin;
    int32_t iprod;

    int32_t ifile, nfiles;

    int32_t nread;
    int32_t offset;
    int32_t offset_out;
    int32_t offmin;
    int32_t offmax;
    int64_t bin_num;
    int64_t bin_num_out;
    int32_t row_write;
    int32_t n_write_total = 0;
    int32_t n_write = 0;
    int32_t nprod[MAXNFILES];
    int32_t ncols;
    int32_t ncols_out;

    int32_t nrows;

    float wgt;
    float *sort_array[MAXNVDATA];

    char buf[FILENAME_MAX];

    /**
     * @brief input buffer for sum and sum-squared data
     **/
    float *in_sum_buf[MAXNVDATA - 2];
    /**
     * @brief output buffer for sum and sum-squared data
     **/
    float *out_sum_buf[MAXNVDATA - 2];

    int status;

    char *ptime = ydhmsf(now(), 'G');
    /**
     * @brief string containing the command and command line options
     **/
    char proc_con[2048];

    char prodname[MAXNVDATA - 2][FILENAME_MAX];

    /**
     * @brief catch all float-32 variable, temporarily holds the sum data 
     **/
    float f32;

    float minlon;
    float maxlon;
    float minlat;
    float maxlat;

    float lat;

    time_t tnow;
    struct tm *tmnow;

    FILE *fp;

    uint8 *nfiles_contribute;

    setlinebuf(stdout);

    printf("%s %s (%s %s)\n", "L3BINMERGE", VERSION, __DATE__, __TIME__);

    if (l3bin_input(argc, argv, &input, "l3binmerge", VERSION) != 0) {
        exit(EXIT_FAILURE);
    }

    //    get_time(ptime);

    input.unit_wgt = 1;
    /*
     * @brief True if output file contains the union of input bins
     **/
    bool case_u = input.union_bins;

    strcpy(proc_con, argv[0]);
    for (i = 1; i < argc; i++) {
        strcat(proc_con, " ");
        strcat(proc_con, argv[i]);
    }

    if (input.loneast <= input.lonwest) {
        printf("loneast: %f must be greater than lonwest: %f.\n", input.loneast, input.lonwest);
        exit(EXIT_FAILURE);
    }

    if (input.latnorth <= input.latsouth) {
        printf("latnorth: %f must be greater than latsouth: %f.\n", input.latnorth, input.latsouth);
        exit(EXIT_FAILURE);
    }

    /* Get lon/lat limits */
    minlon = input.lonwest;
    maxlon = input.loneast;
    minlat = input.latsouth;
    maxlat = input.latnorth;

    /* Determine number of input files */
    /* ------------------------------- */
    nfiles = 0;

    bool isHDF4 = false;
    bool isHDF5 = false;
    bool isCDF4 = false;

    Hdf::hdf_bin *input_binfile[MAXNFILES];

    fp = fopen(input.infile, "r");
    if (fp == NULL) {
        printf("Input listing file: \"%s\" not found.\n", input.infile);
        return -1;
    }
    while (fgets(buf, 256, fp) != NULL)
        nfiles++;
    fclose(fp);
    printf("%d input files\n", nfiles);

    /* Open L3 input files */
    /* ------------------- */
    uint32_t tot_nprod = 0;
    uint32_t prod_offset[256];
    fp = fopen(input.infile, "r");
    for (ifile = 0; ifile < nfiles; ifile++) {
        fgets(buf, 256, fp);
        buf[strlen(buf) - 1] = 0;

        if (ifile == 0) {
            if (Hishdf(buf) == TRUE)
                isHDF4 = true;
            if (H5Fis_hdf5(buf) == TRUE) {
                int ncid;
                char nam_buf[256];
                status = nc_open(buf, NC_NOWRITE, &ncid);
                if (status != NC_NOERR) {
                    isHDF5 = true;
                } else {
                    status = nc_get_att(ncid, NC_GLOBAL, "Mission", nam_buf);
                    if (strcmp(nam_buf, "SAC-D Aquarius") == 0) {
                        isHDF5 = true;
                    } else {
                        isCDF4 = true;
                    }
                    nc_close(ncid);
                }
            }
        }

        if (isHDF4)
            input_binfile[ifile] = new Hdf::hdf4_bin;
        if (isHDF5)
            input_binfile[ifile] = new Hdf::hdf5_bin;
        if (isCDF4)
            input_binfile[ifile] = new Hdf::cdf4_bin;

        if (ifile == 0)
            prod_offset[ifile] = 0;
        else
            prod_offset[ifile] = prod_offset[ifile - 1] + nprod[ifile - 1];

        printf("%d %s\n", ifile, buf);
        input_binfile[ifile]->open(buf);
        nprod[ifile] = input_binfile[ifile]->nprod();
        for (iprod = 0; iprod < nprod[ifile]; iprod++) {
            input_binfile[ifile]->get_prodname(iprod, prodname[tot_nprod]);
            input_binfile[ifile]->active_data_prod[iprod] = true;
            ;
            tot_nprod++;
        }

    } /* ifile loop */

    fclose(fp);

    nrows = input_binfile[0]->nrows;

    /* Create output file */
    /* ------------------ */
    Hdf::hdf_bin *output_binfile;

    if (getFileFormatName(input.oformat) == NULL) {
        if (isHDF4)
            strcpy(input.oformat, "netCDF4");
        if (isHDF5)
            strcpy(input.oformat, "HDF5");
        if (isCDF4)
            strcpy(input.oformat, "netCDF4");
    }

    strcpy(buf, input.ofile);

    if (strcmp(input.oformat, "HDF5") == 0)
        output_binfile = new Hdf::hdf5_bin;
    if (strcmp(input.oformat, "netCDF4") == 0)
        output_binfile = new Hdf::cdf4_bin;

    char *tmpProdNames[tot_nprod];
    for (uint32_t i = 0; i < tot_nprod; i++) {
        tmpProdNames[i] = prodname[i];
    }
    output_binfile->setProductList(tot_nprod, tmpProdNames);
    output_binfile->create(buf, nrows);
    if (isCDF4)
        output_binfile->deflate = input.deflate;

    if (input.noext == 0)
        strcat(buf, ".main");

    /* Allocate I/O buffers */
    /* -------------------- */
    ncols = 2 * nrows;
    ncols_out = 2 * nrows;

    for (iprod = 0; iprod < (int)tot_nprod; iprod++) {
        in_sum_buf[iprod] = (float *)calloc(ncols, 2 * sizeof(float));
        out_sum_buf[iprod] = (float *)calloc(ncols_out, 2 * sizeof(float));
    } /* iprod loop */

    nfiles_contribute = (uint8 *)calloc(ncols_out, sizeof(uint8));

    /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
    /* For each scan ... (Main Loop) */
    /* ----------------------------- */
    for (irow = 0; irow < nrows; irow++) {
        if ((irow % 500) == 0) {
            time(&tnow);
            tmnow = localtime(&tnow);
            printf("irow:%6d of %8d %s", irow, nrows, asctime(tmnow));
        }

        // Get basebin and numbin for this input row
        int64_t basebin = input_binfile[0]->get_basebin(irow);
        int64_t numbin = input_binfile[0]->get_numbin(irow);

        /* Clear output binlist and sum buffers */
        /* ------------------------------------ */
        row_write = 0;

        output_binfile->clear_binlist();
        for (iprod = 0; iprod < (int32_t)tot_nprod; iprod++) {
            memset(&out_sum_buf[iprod][0], 0, ncols_out * 2 * sizeof(float));
        } /* iprod loop */

        if (case_u) {
            for (iprod = 0; iprod < (int32_t)tot_nprod; iprod++) {
                for (i = 0; i < ncols_out * 2; i++) {
                    out_sum_buf[iprod][i] = -32767.0;
                }
            }
        }

        memset(nfiles_contribute, 0, ncols_out * sizeof(uint8));

        /* Get bin info */
        /* ------------ */
        for (ifile = 0; ifile < nfiles; ifile++) {
            input_binfile[ifile]->readBinIndex(irow);

            int ext = input_binfile[ifile]->get_ext();
            if (ext == 0)
                continue;

            /* Read BinList */
            /* ------------ */
            nread = input_binfile[ifile]->readBinList(ext);

            if (nread == -1) {
                printf("Unable to read bin numbers...: %d\n", ext);
            }
        } /* ifile loop */

        /* For each file ... */
        /* ----------------- */
        for (ifile = 0; ifile < nfiles; ifile++) {
            int64_t beg = input_binfile[ifile]->get_beg();
            int ext = input_binfile[ifile]->get_ext();

            /* ========== If row has data ... ========== */
            /* ----------------------------------------- */
            if (beg != 0) {
                /*	printf("row has data: %d\n", irow);*/

                /* Determine lon kbin limits */
                /* ------------------------- */
                offmin = (int32_t)((minlon + 180) * (numbin / 360.0) + 0.5);
                offmax = (int32_t)((maxlon + 180) * (numbin / 360.0) + 0.5);

                /* Get data values (sum, sum_sq) for each filled bin in row */
                /* -------------------------------------------------------- */
                int nbins_to_read = ext;
                for (iprod = 0; iprod < nprod[ifile]; iprod++) {
                    input_binfile[ifile]->readSums(&in_sum_buf[iprod][0], nbins_to_read, iprod);
                } /* iprod loop */

                if (isHDF5 || isCDF4)
                    input_binfile[ifile]->setDataPtr(nbins_to_read);

                /* Skip row if not between minlat & maxlat */
                lat = ((irow + 0.5) / nrows) * 180.0 - 90.0;
                if (lat < minlat || lat > maxlat) {
                    // row_write = 0;
                    continue;
                }

                /* Fill output buffers with input bin data */
                /* --------------------------------------- */
                for (kbin = 0; kbin < ext; kbin++) {
                    /* Store bin number */
                    /* ---------------- */
                    bin_num = input_binfile[ifile]->get_bin_num(kbin);
                    offset = bin_num - basebin;

                    /* If bin outside lon range then skip */
                    /* ---------------------------------- */
                    if (offset < offmin || offset > offmax)
                        continue;

                    float weights = input_binfile[ifile]->get_weights(kbin);

                    /* Skip if not good enough */
                    /* ----------------------- */
                    // Assign output offset & bin number
                    offset_out = offset;
                    bin_num_out = bin_num;

                    if (offset_out >= ncols_out) {
                        printf("Bad write to BINLIST: %d %d %d %d\n", ifile, irow, ncols_out, offset_out);
                        exit(EXIT_FAILURE);
                    }

                    nfiles_contribute[offset_out]++;

                    output_binfile->set_bin_num(offset_out, bin_num_out);
                    row_write = 1;

                    /* Sum & store number of observations,nscenes */
                    /* ------------------------------------------ */
                    int nobs = input_binfile[ifile]->get_nobs(kbin);
                    output_binfile->inc_nobs(offset_out, nobs);
                    int nscenes = input_binfile[ifile]->get_nscenes(kbin);
                    output_binfile->inc_nscenes(offset_out, nscenes);

                    /* Sum & store weights */
                    /* ------------------- */
                    if (input.unit_wgt || input.median) {
                        output_binfile->set_weights(offset_out, 1);
                    } else {
                        output_binfile->inc_weights(offset_out, weights);
                    }

                    /* Product loop */
                    /* ------------ */
                    for (iprod = 0; iprod < nprod[ifile]; iprod++) {
                        if (input.unit_wgt) {
                            wgt = weights;
                            f32 = in_sum_buf[iprod][2 * kbin];
                            if (out_sum_buf[prod_offset[ifile] + iprod][2 * offset_out] == -32767.0) {
                                out_sum_buf[prod_offset[ifile] + iprod][2 * offset_out] = f32 / wgt;
                                out_sum_buf[prod_offset[ifile] + iprod][2 * offset_out + 1] =
                                    (f32 / wgt) * (f32 / wgt);
                            } else {
                                out_sum_buf[prod_offset[ifile] + iprod][2 * offset_out] += f32 / wgt;
                                out_sum_buf[prod_offset[ifile] + iprod][2 * offset_out + 1] +=
                                    (f32 / wgt) * (f32 / wgt);
                            }
                        } else {
                            /* Add new sum to accumulated sum & sum2 */
                            /* ------------------------------------- */
                            /* This else block is never excecuted, input.unit_wgt is hard-coded to 1 */
                            f32 = in_sum_buf[iprod][2 * kbin];
                            if (out_sum_buf[prod_offset[ifile] + iprod][2 * offset_out] == -32767.0) {
                                out_sum_buf[prod_offset[ifile] + iprod][2 * offset_out] = f32;
                            } else {
                                out_sum_buf[prod_offset[ifile] + iprod][2 * offset_out] += f32;
                            }

                            f32 = in_sum_buf[iprod][2 * kbin + 1];
                            if (out_sum_buf[prod_offset[ifile] + iprod][2 * offset_out + 1] == -32767.0) {
                                out_sum_buf[prod_offset[ifile] + iprod][2 * offset_out + 1] = f32;
                            } else {
                                out_sum_buf[prod_offset[ifile] + iprod][2 * offset_out + 1] += f32;
                            }
                        } /* input.unit_wgt */
                    } /* product loop */

                } /* kbin loop */

            } /* ========== if row has data ========== */

        } /* ifile loop */

        /* Write output vdatas */
        /* ------------------- */
        if (row_write) {
            n_write = 0;
            for (kbin = 0; kbin < output_binfile->get_numbin(irow); kbin++) {
                bin_num = output_binfile->get_bin_num(kbin);

                if (bin_num != 0) {
                    // Skip rows where not all files contribute
                    if ((nfiles_contribute[bin_num - basebin] != nfiles) and !case_u)
                        continue;

                    /* Loop over data products */
                    /* ----------------------- */
                    for (iprod = 0; iprod < (int32_t)tot_nprod; iprod++) {
                        /* Remove "blank" bin records */
                        /* -------------------------- */
                        if (n_write != kbin)
                            memcpy(&out_sum_buf[iprod][2 * n_write], &out_sum_buf[iprod][2 * kbin], 8);
                    } /* iprod loop */

                    /* Remove "blank" bin records */
                    /* -------------------------- */
                    if (n_write != kbin)
                        output_binfile->copy_binlist(kbin, n_write);

                    n_write++;
                    n_write_total++;

                } /* bin_num != 0 */
            } /* kbin loop */

            /* Write BinList & Data Products */
            /* ----------------------------- */
            if (n_write > 0) {
                output_binfile->writeBinList(n_write);

                for (iprod = 0; iprod < (int32_t)tot_nprod; iprod++) {
                    strcpy(buf, prodname[iprod]);
                    output_binfile->writeSums(&out_sum_buf[iprod][0], n_write, buf);
                }
                if (isHDF5 || isCDF4)
                    output_binfile->incNumRec(n_write);
            }
        }
        // if median free storage arrays
        if (input.median) {
            for (iprod = 0; iprod < (int32_t)tot_nprod; iprod++) {
                free(sort_array[iprod]);
            }
        }

    } /* irow loop */
    /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

    for (iprod = 0; iprod < (int32_t)tot_nprod; iprod++) {
        free(in_sum_buf[iprod]);
        free(out_sum_buf[iprod]);
    } /* iprod loop */

    free(nfiles_contribute);

    // Copy metadata from input to output binfile
    strcpy(buf, input.ofile);
    strcpy(output_binfile->meta_l3b.product_name, buf);

    input.tflag = toupper(input.tflag);
    if (input.tflag == 'D') {
        strcpy(output_binfile->meta_l3b.prod_type, "day");
    } else if (input.tflag == 'W') {
        strcpy(output_binfile->meta_l3b.prod_type, "8-day");
    } else if (input.tflag == 'M') {
        strcpy(output_binfile->meta_l3b.prod_type, "month");
    } else if (input.tflag == 'Y') {
        strcpy(output_binfile->meta_l3b.prod_type, "year");
    } else if (input.tflag == 'O') {
        strcpy(output_binfile->meta_l3b.prod_type, "other");
    }
    strcpy(output_binfile->meta_l3b.pversion, input.pversion);
    strcpy(output_binfile->meta_l3b.soft_name, "L3BINMERGE");
    strcpy(output_binfile->meta_l3b.soft_ver, VERSION);
    strcpy(output_binfile->meta_l3b.proc_con, proc_con);
    strcpy(output_binfile->meta_l3b.input_parms, input.parms);

    char buf2[5000];
    if (Hishdf(input.infile) == TRUE || H5Fis_hdf5(input.infile) == TRUE) {
        strcpy(buf2, input.infile);
    } else {
        fp = fopen(input.infile, "r");
        buf2[0] = 0;
        for (ifile = 0; ifile < nfiles; ifile++) {
            fgets(buf, 256, fp);
            buf[strlen(buf) - 1] = 0;

            strcat(buf2, buf);
            strcat(buf2, ",");
        } /* ifile loop */
        fclose(fp);
        buf2[strlen(buf2) - 1] = 0;
    }
    strcpy(output_binfile->meta_l3b.infiles, buf2);

    output_binfile->copymeta(nfiles, &input_binfile[0]);

    if (strcmp(input.oformat, "netCDF4") == 0) {
        ptime = unix2isodate(now(), 'G');
        strcpy(output_binfile->meta_l3b.ptime, ptime);
    } else {
        strcpy(output_binfile->meta_l3b.ptime, ptime);
    }

    for (ifile = 0; ifile < nfiles; ifile++)
        input_binfile[ifile]->close();
    output_binfile->close();

    return 0;
}
