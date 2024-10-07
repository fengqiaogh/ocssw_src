#include <stdio.h>
#include <math.h>
#include <time.h>
#include <libgen.h>
#include <sys/types.h>

#include "netcdf.h"

#include "l3bin_input.h"
#include <timeutils.h>
#include <genutils.h>
#include "hdf_bin.h"
#include "sensorInfo.h"

#include <hdf.h>
#include <mfhdf.h>

using namespace std;

#define MAXNFILES 256

#define BYTE    unsigned char

#define BINCHECK -1
#define L3BIN_CACHE_SIZE 8 * 1024  // 8 kb of cache memory. 
#define L3BIN_CACHE_NELEMS 512
#define L3BIN_CACHE_PREEMPTION .75
#define VERSION "5.14"

static instr input;

#define EXIT_STATUS(func,status,...) {int status = func; if(status!=NC_NOERR) {printf("--Error--: %s returned non-zero exit code. \n",#func); printf(__VA_ARGS__); exit(EXIT_FAILURE);} }

int main(int argc, char **argv) {
    intn i;
    EXIT_STATUS(nc_set_chunk_cache(L3BIN_CACHE_SIZE, L3BIN_CACHE_NELEMS, L3BIN_CACHE_PREEMPTION),status,"Setting cache size failed. See %s, line %d\n. The program ran on %s at %s. Exiting.",__FILE__,__LINE__, __DATE__, __TIME__);
    int retval = 0;
    int status;
    int32_t irow;
    int32_t kbin;
    int32_t iprod;
    int activeProdId;

    int32_t ifile, nfiles;

    int32_t nread;
    int32_t offset;
    int64_t offset_out;
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

    int16_t composite_prod_index[MAXNFILES];
    int16_t composite_outprod_index;

    int32_t nrows;

    int32_t reduce_fac;

    float wgt;
    float *sort_array[MAXNVDATA];

    char buf[FILENAME_MAX];
    char filename_L3[FILENAME_MAX];

    float *in_sum_buf[MAXNVDATA - 2];
    float *out_sum_buf[MAXNVDATA - 2];
    uint8 * in_qual_buf[MAXNFILES + 1], *uint8_buf;
    char ptime[17];
    char proc_con[2048];

    float f32;

    float minlon;
    float maxlon;
    float minlat;
    float maxlat;

    float lat;

    time_t tnow;
    struct tm *tmnow;

    FILE *fp, *fp_L3;

    void insertion_sort(float a[], int length);

    setlinebuf(stdout);

    printf("%s %s (%s %s)\n", "L3BIN", VERSION, __DATE__, __TIME__);

    if (l3bin_input(argc, argv, &input, "l3bin", VERSION) != 0) {
        exit(1);
    }

    get_time(ptime);


    strcpy(proc_con, argv[0]);
    for (i = 1; i < argc; i++) {
        strcat(proc_con, " ");
        strcat(proc_con, argv[i]);
    }

    if (input.loneast <= input.lonwest) {
        printf("loneast: %f must be greater than lonwest: %f.\n",
                input.loneast, input.lonwest);
        exit(-1);
    }

    if (input.latnorth <= input.latsouth) {
        printf("latnorth: %f must be greater than latsouth: %f.\n",
                input.latnorth, input.latsouth);
        exit(-1);
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

    Hdf::hdf_bin * input_binfile[MAXNFILES];

    /* Single HDF input */
    /* ---------------- */
    if (Hishdf(input.infile) == TRUE || H5Fis_hdf5(input.infile) == TRUE) {
        printf("Single input file\n");
        nfiles = 1;

        if (Hishdf(input.infile) == TRUE) {
            isHDF4 = true;
            input_binfile[0] = new Hdf::hdf4_bin;
        }

        if (H5Fis_hdf5(input.infile) == TRUE) {
            int ncid;
            char nam_buf[256];
            bzero(nam_buf, 256);
            status = nc_open(input.infile, NC_NOWRITE, &ncid);
            if (status != NC_NOERR) {
                isHDF5 = true;
                input_binfile[0] = new Hdf::hdf5_bin;
            } else {
                status = nc_get_att(ncid, NC_GLOBAL, "Mission", nam_buf);
                if (status != NC_NOERR)
                    status = nc_get_att(ncid, NC_GLOBAL, "mission", nam_buf);

                if ((status == NC_NOERR) &&
                        ((strcmp(nam_buf, "SAC-D Aquarius") == 0) ||
                        (strcmp(nam_buf, "SMAP") == 0))) {
                    nc_close(ncid);
                    isHDF5 = true;
                    input_binfile[0] = new Hdf::hdf5_bin;
                } else {
                    isCDF4 = true;
                    nc_get_att(ncid, NC_GLOBAL, "title", nam_buf);
                            
                    char *output = NULL;
                    output = strstr(nam_buf,"Level-3 Binned Data");
                    if(!output) {
                        printf("Input file must be a Level-3 file. Please verify and retry\n");
                        exit(EXIT_FAILURE);
                    }

                    input_binfile[0] = new Hdf::cdf4_bin;
                }
            }
        }

        input_binfile[0]->open(input.infile);
        nprod[0] = input_binfile[0]->nprod();
    } else {
        // Work with a list of files
        fp = fopen(input.infile, "r");
        if (fp == NULL) {
            printf("Input listing file: \"%s\" not found.\n", input.infile);
            return -1;
        }
        while (fgets(buf, 256, fp) != NULL) {
            // Verify that the file (in this list of files) does indeed exist
            parse_file_name(buf, filename_L3);
            fp_L3 = fopen(filename_L3, "r");
            if (fp_L3 == NULL) {
                printf("Input file: \"%s\" not found.\n", filename_L3);
                return -1;
            }
            fclose(fp_L3);
            // Increment count of files (in this list of files)
            nfiles++;
        }
        fclose(fp);
        printf("%d input files\n", nfiles);


        /* Open L3 input files */
        /* ------------------- */
        fp = fopen(input.infile, "r");
        for (ifile = 0; ifile < nfiles; ifile++) {
            fgets(buf, 256, fp);
            buf[strcspn(buf, "\n")] = 0;
            buf[strlen(buf)] = 0;
            // buf[strlen(buf) - 1] = 0;

            if (ifile == 0) {
                if (Hishdf(buf) == TRUE) {
                    isHDF4 = true;
                } else if (H5Fis_hdf5(buf) == TRUE) {
                    int ncid;
                    char nam_buf[256];
                    status = nc_open(buf, NC_NOWRITE, &ncid);
                    if (status != NC_NOERR) {
                        isHDF5 = true;
                    } else {
                        status = nc_get_att(ncid, NC_GLOBAL, "Mission", nam_buf);
                        if (status != NC_NOERR)
                            status = nc_get_att(ncid, NC_GLOBAL, "mission", nam_buf);
                        if (status == NC_NOERR) {
                            if ((strcmp(nam_buf, "SAC-D Aquarius") == 0) ||
                                    (strcmp(nam_buf, "SMAP") == 0)) {
                                isHDF5 = true;
                            } else {
                                isCDF4 = true;
                                nc_get_att(ncid, NC_GLOBAL, "title", nam_buf);
                                char *output = NULL;
                                output = strstr(nam_buf,"Level-3 Binned Data");
                                if(!output) {
                                    printf("Input files must be Level-3 files. Please verify and retry\n");
                                    exit(EXIT_FAILURE);
                                }
                            }
                        } else {
                            isCDF4 = true;
                            nc_get_att(ncid, NC_GLOBAL, "title", nam_buf);
                            char *output = NULL;
                            output = strstr(nam_buf,"Level-3 Binned Data");
                            if(!output) {
                                printf("Input files must be Level-3 files. Please verify and retry\n");
                                exit(EXIT_FAILURE);
                            }
                        }
                        nc_close(ncid);
                    }
                }
            }
            if (isHDF4) 
                input_binfile[ifile] = new Hdf::hdf4_bin;
            else if (isHDF5)
                input_binfile[ifile] = new Hdf::hdf5_bin;
            else if (isCDF4)
                input_binfile[ifile] = new Hdf::cdf4_bin;

            printf("%d %s\n", ifile, buf);
            input_binfile[ifile]->open(buf);
            nprod[ifile] = input_binfile[ifile]->nprod();

            //printf("open status: %d\n", status);

        } /* ifile loop */

        fclose(fp);
    }

    nrows = input_binfile[0]->nrows;

    if(input.resolve[0] != '\0') {
        int32_t out_nrows = resolve2binRows(input.resolve);
        if(out_nrows <= 0) {
            printf("-E- unknown resolve param = %s\n", input.resolve);
            exit(EXIT_FAILURE);
        }
        reduce_fac = nrows / out_nrows;

    } else {
        reduce_fac = input.reduce_fac;
    }

    if (reduce_fac != 1 &&
            reduce_fac != 2 &&
            reduce_fac != 4 &&
            reduce_fac != 8 &&
            reduce_fac != 16) {
        printf("Reduction factor must be power of 2 less than 16\n");
        exit(EXIT_FAILURE);
    }

    /* Generate output product list from 1st input file if DEFAULT */
    /* ----------------------------------------------------------- */
    if (strcmp(input.out_parm, ":DEFAULT:") == 0) {

        // Read input file product list
        input_binfile[0]->query(input.out_parm);

        //    printf("out_parm: %s\n", &input.out_parm[1]);
    } else {
        strcpy(buf, &input.out_parm[1]);
        buf[strlen(buf) - 1] = 0;
        for (i = 0; i < (intn) strlen(buf); i++)
            if (buf[i] == ':') buf[i] = ',';
        strcpy(input.out_parm, buf);
    }

    /* Determine active products */
    /* ------------------------- */
    for (ifile = 0; ifile < nfiles; ifile++) {
        int status = input_binfile[ifile]->read(input.out_parm);
        if (status == -1) {
            printf("Not all output products found in input file %d\n",
                    ifile);
            exit(-1);
        }
    }

    /* Check whether Composite product exists in L3 inputfiles */
    /* ------------------------------------------------------- */
    if (input.composite_prod[0] != 0) {
        for (ifile = 0; ifile < nfiles; ifile++) {
            int tmpNumProducts = input_binfile[ifile]->n_active_prod;
            int i;
            for (i = 0; i < tmpNumProducts; i++)
                if (strcmp(input.composite_prod,
                        input_binfile[0]->getActiveProdName(i)) == 0) break;

            composite_prod_index[ifile] = i;

            if (i == tmpNumProducts) {
                printf("Composite product: \"%s\" not found in L3 dataset %d\n",
                        input.composite_prod, ifile);
                exit(-1);
            }
        }
    }

#if 0
    /* Check whether NDVI product exists in L3 inputfiles */
    /* -------------------------------------------------- */
    if (input.composite_prod[0] != 0) {
        for (ifile = 0; ifile < nfiles; ifile++) {
            int tmpNumProducts = input_binfile[ifile]->n_active_prod;
            int i;
            for (i = 0; i < tmpNumProducts; i++)
                if (strcmp("ndvi",
                        input_binfile[0]->getActiveProdName(i)) == 0) break;

            if (i == tmpNumProducts) {
                printf("NDVI product not found in L3 dataset %d\n", ifile);
                exit(-1);
            }
        }
    }
#endif

    /* Create output file */
    /* ------------------ */
    Hdf::hdf_bin *output_binfile;

    if (getFileFormatName(input.oformat) == NULL) {
        if (isHDF4) strcpy(input.oformat, "HDF4");
        if (isHDF5) strcpy(input.oformat, "HDF5");
        if (isCDF4) strcpy(input.oformat, "netCDF4");
    } else {
        // normalize the oformat name
        strcpy(input.oformat, getFileFormatName(input.oformat));
    }

    strcpy(buf, input.ofile);

    if (strcmp(input.oformat, "HDF4") == 0) {
        output_binfile = new Hdf::hdf4_bin;
        if (input.noext == 0) {
            output_binfile->hasNoext = false;
            strcat(buf, ".main");
        } else {
            output_binfile->hasNoext = true;
        }

    } else if (strcmp(input.oformat, "HDF5") == 0) {
        output_binfile = new Hdf::hdf5_bin;
    } else if (strcmp(input.oformat, "netCDF4") == 0) {
        output_binfile = new Hdf::cdf4_bin;
        output_binfile->deflate = input.deflate;
    }

    int tmpNumProducts = input_binfile[0]->n_active_prod;
    char* tmpProductNames[tmpNumProducts];
    for (int i = 0; i < tmpNumProducts; i++) {
        tmpProductNames[i] = (char*) input_binfile[0]->getActiveProdName(i);
    }


    output_binfile->setProductList(tmpNumProducts, tmpProductNames);
    if (input_binfile[0]->has_qual()) {
        output_binfile->hasQual = true;
    }
    output_binfile->create(buf, nrows / reduce_fac);


    /* Allocate I/O buffers */
    /* -------------------- */
    ncols = 2 * nrows;
    ncols_out = 2 * nrows / reduce_fac;

    activeProdId = 0;
    for (iprod = 0; iprod < nprod[0]; iprod++) {
        if (input_binfile[0]->active_data_prod[iprod] == true) {
            in_sum_buf[activeProdId] = (float *) calloc(ncols, 2 * sizeof (float));
            out_sum_buf[activeProdId] = (float *) calloc(ncols_out, 2 * sizeof (float));

            if (input.composite_prod[0] != 0) {
                if (strcmp(input.composite_prod,
                        input_binfile[0]->getProdName(iprod)) == 0)
                    composite_outprod_index = activeProdId;
            }
            activeProdId++;
        }
    } /* iprod loop */


    /* Allocate quality buffer */
    /* ----------------------- */
    for (ifile = 0; ifile < nfiles; ifile++) {
        in_qual_buf[ifile] = (uint8 *) calloc(ncols, sizeof (uint8));
    }
    in_qual_buf[nfiles] = (uint8 *) calloc(ncols, sizeof (uint8));

    uint8_buf = (uint8 *) calloc(ncols, sizeof (uint8));


    /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
    /* For each scan ... (Main Loop) */
    /* ----------------------------- */
    for (irow = 0; irow < nrows; irow++) {

        if ((irow % 500) == 0) {
            time(&tnow);
            tmnow = localtime(&tnow);
            printf("irow:%6d of %8d %s", irow, nrows, asctime(tmnow));
        }

        int32_t max_out_kbin;
        if ((irow % reduce_fac) == 0)
            max_out_kbin = -1;

        // Get basebin and numbin for this input row
        int64_t basebin = input_binfile[0]->get_basebin(irow);
        int64_t numbin = input_binfile[0]->get_numbin(irow);

        double ratio = 1.0;
        if (reduce_fac > 1)
            ratio = 1.0 * input_binfile[0]->get_numbin(irow) /
            output_binfile->get_numbin(irow / reduce_fac);


        // If median allocate and initialize storage
        if (input.median) {
            for (activeProdId = 0; activeProdId < input_binfile[0]->n_active_prod; activeProdId++) {
                if (input_binfile[0]->active_data_prod[activeProdId] == true) {
                    sort_array[activeProdId] =
                            (float *) calloc(nfiles*numbin, sizeof (float));
                    for (int64_t i = 0; i < nfiles * numbin; i++) sort_array[activeProdId][i] = -999;
                }
            }
        }


        /* Clear output binlist, sum and quality buffers */
        /* --------------------------------------------- */
        if ((irow % reduce_fac) == 0) {
            row_write = 0;

            output_binfile->clear_binlist();
            for (activeProdId = 0; activeProdId < input_binfile[0]->n_active_prod; activeProdId++) {
                memset(&out_sum_buf[activeProdId][0], 0, ncols_out * 2 * sizeof (float));
            } /* iprod loop */

            memset(in_qual_buf[nfiles], 255, ncols);
            for (ifile = 0; ifile < nfiles; ifile++) {
                memset(in_qual_buf[ifile], 255, ncols);
            }
        }

        /* Get bin & qual info */
        /* ------------------- */
        for (ifile = 0; ifile < nfiles; ifile++) {
            input_binfile[ifile]->readBinIndex(irow);

            int ext = input_binfile[ifile]->get_ext();

            /* Read BinList */
            /* ------------ */
            // Get current binlist pointer
            int32_t list_ptr = input_binfile[ifile]->get_list_ptr();

            if (ext > 0) {
                nread = input_binfile[ifile]->readBinList(ext);

                if (nread == -1) {
                    printf("Unable to read bin numbers...: %d\n", ext);
                }
            }

            /* Read quality if vdata exists */
            /* Note: in_qual_buf is "uncompressed" along row */
            if (input_binfile[ifile]->has_qual()) {
                if ((irow % reduce_fac) == 0) {
                    int32_t ext_qual = ext;
                    for (int32_t i = 0; i < reduce_fac; i++) {
                        switch (i) {
                        case 0:
                            input_binfile[ifile]->readQual(uint8_buf, ext_qual, list_ptr);
                            break;
                        default:
                            input_binfile[ifile]->readBinIndex(irow + i);
                            ext_qual = input_binfile[ifile]->get_ext();
                            if (ext_qual != 0) {
                                nread = input_binfile[ifile]->readBinList(ext_qual);
                                input_binfile[ifile]->readQual(uint8_buf, ext_qual);
                            }
                        } // switch

                        // Update in_qual_buf for (irow+i)
                        if (ext_qual != 0) {
                            int64_t basebin = input_binfile[0]->get_basebin(irow + i);

                            double ratio = 1.0 * input_binfile[0]->get_numbin(irow + i) /
                                    output_binfile->get_numbin(irow / reduce_fac);

                            for (kbin = 0; kbin < ext_qual; kbin++) {
                                bin_num = input_binfile[ifile]->get_bin_num(kbin);
                                offset = bin_num - basebin;
                                if (offset < 0) {
                                    cout << "bin_num - basebin is negative for ifile: " << ifile;
                                    cout << " irow: " << irow << " kbin: " << kbin << endl;
                                    cout << "bin_num: " << bin_num << "  basebin: " << basebin
                                            << endl;
                                    exit(1);
                                }

                                int32_t j = reduce_fac * (int32_t) ((offset / ratio) + 0.0);
                                if ((j < ncols) && (in_qual_buf[ifile][j] > uint8_buf[kbin])) {
                                    in_qual_buf[ifile][j] = uint8_buf[kbin];
                                }
                            }
                        } // ext_qual != 0
                    } // i (reduce_fac) loop

                    // Reset to current row if reduce_fac != 1
                    if (reduce_fac != 1) {
                        input_binfile[ifile]->readBinIndex(irow);
                        nread = input_binfile[ifile]->readBinList(ext, list_ptr);
                    }

                } // if ( (irow % reduce_fac) == 0)
            } // if ( input_binfile[ifile]->has_qual())
        } /* ifile loop */


        /* Find best quality */
        /* ----------------- */
        for (kbin = 0; kbin < ncols; kbin++) {
            for (ifile = 0; ifile < nfiles; ifile++) {
                int32_t j = reduce_fac * (int32_t) ((kbin / ratio) + 0.0);
                if (j < ncols) {
                    if (in_qual_buf[ifile][j] < in_qual_buf[nfiles][j])
                        in_qual_buf[nfiles][j] = in_qual_buf[ifile][j];
                }
            }
        }

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
                offmin =
                        (int32_t) ((minlon + 180) * (numbin / 360.0) + 0.5);
                offmax =
                        (int32_t) ((maxlon + 180) * (numbin / 360.0) + 0.5);


                /* Get data values (sum, sum_sq) for each filled bin in row */
                /* -------------------------------------------------------- */
                int nbins_to_read = ext;
                activeProdId = 0;
                for (iprod = 0; iprod < nprod[ifile]; iprod++) {
                    if (input_binfile[ifile]->active_data_prod[iprod] == true) {

                        input_binfile[ifile]->readSums(&in_sum_buf[activeProdId][0],
                                nbins_to_read, iprod);
                        activeProdId++;
                    }
                } /* iprod loop */
                if (isHDF5 || isCDF4)
                    input_binfile[ifile]->setDataPtr(nbins_to_read);
                //	row_write = 1;


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
                    if (offset < 0) {
                        cout << "bin_num - basebin is negative for ifile: " << ifile;
                        cout << " irow: " << irow << " kbin: " << kbin << endl;
                        cout << "bin_num: " << bin_num << "  basebin: " << basebin << endl;
                        exit(1);
                    }

                    /* If bin outside lon range then skip */
                    /* ---------------------------------- */
                    if (offset < offmin || offset > offmax)
                        continue;

                    float weights = input_binfile[ifile]->get_weights(kbin);
                    float time_rec = input_binfile[ifile]->get_time_rec(kbin);

                    /* Skip if not good enough */
                    /* ----------------------- */
                    int32_t j = reduce_fac * (int32_t) ((offset / ratio) + 0.0);
                    if (j < ncols) {
                        if (in_qual_buf[ifile][j] > in_qual_buf[nfiles][j])
                            continue;
                    }

                    // Assign output offset & bin number
                    if (reduce_fac > 1) {
                        offset_out = (int64_t) ((offset / ratio) + 0.0);
                        //	    if ( offset_out == ncols_out)
                        // offset_out = ncols_out - 1;
                        bin_num_out =
                                offset_out + output_binfile->get_basebin(irow / reduce_fac);
                    } else {
                        offset_out = offset;
                        bin_num_out = bin_num;
                    }

                    if (offset_out >= ncols_out) {
                        printf("Bad write to BINLIST: %d %d %d %ld\n",
                                ifile, irow, ncols_out, (long int) offset_out);
                        exit(1);
                    }

                    if (offset_out > max_out_kbin)
                        max_out_kbin = offset_out;

                    output_binfile->set_bin_num(offset_out, bin_num_out);
                    row_write = 1;

                    /* Sum & store number of observations,nscenes */
                    /* ------------------------------------------ */
                    int nobs, nscenes;
                    if (input.composite_prod[0] == 0) {
                        nobs = input_binfile[ifile]->get_nobs(kbin);
                        output_binfile->inc_nobs(offset_out, nobs);
                        nscenes = input_binfile[ifile]->get_nscenes(kbin);
                        output_binfile->inc_nscenes(offset_out, nscenes);

                        /* Sum & store weights */
                        /* ------------------- */
                        if (input.unit_wgt || input.median) {
                            output_binfile->set_weights(offset_out, 1);
                        } else {
                            output_binfile->inc_weights(offset_out, weights);
                        }

                        output_binfile->inc_time_rec(offset_out, time_rec);
                    } else {
                        if (output_binfile->get_nobs(offset_out) != 0) {
                            f32 = in_sum_buf[composite_prod_index[ifile]][2 * kbin] / weights;
                            if (strcmp(input.composite_scheme, "max") == 0) {
                                if (f32 < (out_sum_buf[composite_outprod_index][2 * offset_out] / output_binfile->get_weights(offset_out))) continue;
                            } else {
                                if (f32 > (out_sum_buf[composite_outprod_index][2 * offset_out] / output_binfile->get_weights(offset_out))) continue;
                            }
                        }
                        nobs = input_binfile[ifile]->get_nobs(kbin);
                        output_binfile->set_nobs(offset_out, nobs);
                        nscenes = input_binfile[ifile]->get_nscenes(kbin);
                        output_binfile->set_nscenes(offset_out, nscenes);
                        output_binfile->set_weights(offset_out, weights);
                    }

                    /* Product loop */
                    /* ------------ */
                    for (activeProdId = 0; activeProdId < input_binfile[0]->n_active_prod; activeProdId++) {
                        if (input.unit_wgt) {
                            wgt = weights;
                            f32 = in_sum_buf[activeProdId][2 * kbin];
                            out_sum_buf[activeProdId][2 * offset_out] += f32 / wgt;
                            out_sum_buf[activeProdId][2 * offset_out + 1] += (f32 / wgt)*(f32 / wgt);
                        } else if (input.median) {
                            wgt = weights;
                            f32 = in_sum_buf[activeProdId][2 * kbin];
                            sort_array[activeProdId][nfiles * offset + ifile] = f32 / wgt;
                        } else if (input.composite_prod[0] != 0) {
                            f32 = in_sum_buf[activeProdId][2 * kbin];
                            out_sum_buf[activeProdId][2 * offset_out] = f32;
                            f32 = in_sum_buf[activeProdId][2 * kbin + 1];
                            out_sum_buf[activeProdId][2 * offset_out + 1] = f32;
                        } else {
                            /* Add new sum to accumulated sum & sum2 */
                            /* ------------------------------------- */
                            f32 = in_sum_buf[activeProdId][2 * kbin];
                            out_sum_buf[activeProdId][2 * offset_out] += f32;
                            f32 = in_sum_buf[activeProdId][2 * kbin + 1];
                            out_sum_buf[activeProdId][2 * offset_out + 1] += f32;
                        } /* input.unit_wgt */

                    } /* product loop */

                } /* kbin loop */

            } /* ========== if row has data ========== */

        } /* ifile loop */


        // Compute median
        if (row_write && input.median) {
            for (kbin = 0; kbin < numbin; kbin++) {
                for (activeProdId = 0; activeProdId < input_binfile[0]->n_active_prod; activeProdId++) {
                    float *sort_buf = (float *) calloc(nfiles, sizeof (float));
                    int nsort = 0;
                    for (ifile = 0; ifile < nfiles; ifile++) {
                        f32 = sort_array[activeProdId][nfiles * kbin + ifile];
                        if (f32 != -999)
                            sort_buf[nsort++] = sort_array[activeProdId][nfiles * kbin + ifile];
                    }
                    // Call insertion sort
                    if (nsort > 0) {
                        insertion_sort(sort_buf, nsort);
                        out_sum_buf[activeProdId][2 * (kbin / reduce_fac)] = sort_buf[nsort / 2];
                        out_sum_buf[activeProdId][2 * (kbin / reduce_fac) + 1] = 1.0;
                        free(sort_buf);
                    }
                }
            }
        }

        /* Write output vdatas */
        /* ------------------- */
        if (row_write && (irow % reduce_fac) == reduce_fac - 1) {

            n_write = 0;
            for (kbin = 0; kbin <= max_out_kbin; kbin++) {

                int64_t bin_num = output_binfile->get_bin_num(kbin);
                if (bin_num != 0) {

                    /* Loop over data products */
                    /* ----------------------- */
                    for (activeProdId = 0; activeProdId < input_binfile[0]->n_active_prod; activeProdId++) {

                        /* Remove "blank" bin records */
                        /* -------------------------- */
                        if (n_write != kbin)
                            memcpy(&out_sum_buf[activeProdId][2 * n_write],
                                &out_sum_buf[activeProdId][2 * kbin], 8);
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

                activeProdId = 0;
                for (iprod = 0; iprod < nprod[0]; iprod++) {
                    if (input_binfile[0]->active_data_prod[iprod] == true) {

                        input_binfile[0]->get_prodname(iprod, buf);
                        output_binfile->writeSums(&out_sum_buf[activeProdId][0], n_write, buf);
                        activeProdId++;
                    }
                }
                if (strcmp(input.oformat, "HDF5") == 0 ||
                        strcmp(input.oformat, "netCDF4") == 0)
                    output_binfile->incNumRec(n_write);
                //	if ( isHDF5 || isCDF4) output_binfile->incNumRec( n_write);
            }

            /* Write to Quality Vdata */
            /* ---------------------- */
            i = 0;
            int64_t nbin;
            if (irow < nrows / 2)
                nbin = input_binfile[0]->get_numbin(irow);
            else
                nbin = input_binfile[0]->get_numbin(irow - (reduce_fac / 2));

            offmin =
                    (int32_t) ((minlon + 180) * (nbin / 360.0) + 0.5);
            offmax =
                    (int32_t) ((maxlon + 180) * (nbin / 360.0) + 0.5);

            for (kbin = 0; kbin < ncols; kbin++) {
                if (kbin < offmin || kbin > offmax) continue;
                if (in_qual_buf[nfiles][kbin] != 255) {
                    in_qual_buf[nfiles][i++] = in_qual_buf[nfiles][kbin];
                }
            }

            if (input_binfile[0]->has_qual()) {
                if ((i - n_write) > 2) {
                    cout << "Problem with Quality write: irow: " << irow << " i: " <<
                            i << " n_write: " << n_write << " " << max_out_kbin << endl;
                    exit(1);
                }

                output_binfile->writeQual(in_qual_buf[nfiles], n_write);
            } // has_qual() == true
        } /* row_write = 1 */


        // if median free storage arrays
        if (input.median) {
            for (activeProdId = 0; activeProdId < input_binfile[0]->n_active_prod; activeProdId++) {
                free(sort_array[activeProdId]);
            }
        }

    } /* irow loop (Main loop) */
    /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

    for (ifile = 0; ifile <= nfiles; ifile++) {
        free(in_qual_buf[ifile]);
    }
    free(uint8_buf);

    activeProdId = 0;
    for (activeProdId = 0; activeProdId < input_binfile[0]->n_active_prod; activeProdId++) {
        free(in_sum_buf[activeProdId]);
        free(out_sum_buf[activeProdId]);
    } /* iprod loop */


    if(output_binfile->n_data_records > 0) {
    
        // Copy metadata from input to output binfile
        output_binfile->copymeta(nfiles, &input_binfile[0]);
        int sensorID = sensorName2SensorId(output_binfile->meta_l3b.sensor_name);
        output_binfile->meta_l3b.sensorID = sensorID;

        // put in a fix to the missing Mission
        if (strlen(output_binfile->meta_l3b.mission) == 0) {
            if (sensorID != -1) {
                strcpy(output_binfile->meta_l3b.mission, sensorId2PlatformName(sensorID));
            }
        }

        strcpy(buf, input.ofile);
        strcpy(output_binfile->meta_l3b.product_name, buf);
        strcpy(output_binfile->meta_l3b.pversion, input.pversion);
        strcpy(output_binfile->meta_l3b.soft_name, "L3BIN");
        strcpy(output_binfile->meta_l3b.soft_ver, VERSION);
        if (isCDF4 == 1) {
            // 1994-11-05T13:15:30Z
            strcpy(output_binfile->meta_l3b.ptime, unix2isodate(tnow, 'G'));
        } else {
            // yyyydddhhmmssmmm
            strcpy(output_binfile->meta_l3b.ptime, ydhmsf(tnow, 'L'));
        }
        strcpy(output_binfile->meta_l3b.ptime, ptime);
        strcpy(output_binfile->meta_l3b.proc_con, proc_con);
        strcpy(output_binfile->meta_l3b.input_parms, input.parms);

        if(strlen(input.doi)) {
            strcpy(output_binfile->meta_l3b.doi, input.doi);
        } else {
            strcpy(output_binfile->meta_l3b.doi, "");
        }

        char buf2[LG_ATTRSZ];
        if (Hishdf(input.infile) == TRUE || H5Fis_hdf5(input.infile) == TRUE) {
            strcpy(output_binfile->meta_l3b.infiles, basename(input.infile));

        } else {

            fp = fopen(input.infile, "r");
            buf2[0] = 0;
            for (ifile = 0; ifile < nfiles; ifile++) {
                fgets(buf, 256, fp);
                buf[strlen(buf) - 1] = 0;

                strcat(buf2, basename(buf));
                strcat(buf2, ",");
            } /* ifile loop */
            fclose(fp);
            buf2[strlen(buf2) - 1] = 0;
            strncpy(output_binfile->meta_l3b.infiles, buf2, LG_ATTRSZ - 1);
        }
    } else {
        retval = 110;
        printf("No valid data was binned\n");
    }
    
    
    for (ifile = 0; ifile < nfiles; ifile++) {
        input_binfile[ifile]->close();
        delete input_binfile[ifile];
    }
    
    output_binfile->close();
    delete output_binfile;

    if(retval == 110) {
        string cmd = "rm -f ";
        cmd += input.ofile;
        system(cmd.c_str());
    }
    
    printf("Done\n");
    return retval;
}


/* Copyright (c) 2009 the authors listed at the following URL, and/or
the authors of referenced articles or incorporated external code:
http://en.literateprograms.org/Insertion_sort_(C)?action=history&offset=20081205204844

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Retrieved from: http://en.literateprograms.org/Insertion_sort_(C)?oldid=15530
 */

/* Sort an array of integers */
void insertion_sort(float a[], int length) {
    int i;
    for (i = 0; i < length; i++) {
        /* Insert a[i] into the sorted sublist */
        int j;
        float v = a[i];

        for (j = i - 1; j >= 0; j--) {
            if (a[j] <= v) break;
            a[j + 1] = a[j];
        }
        a[j + 1] = v;

    }
}

