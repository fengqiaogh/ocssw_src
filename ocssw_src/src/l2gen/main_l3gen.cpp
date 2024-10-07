// =====================================================================
// l3gen - level-3 to level-3 ocean algorithm processor                 
// B. Franz, NASA/OBPG, August 2008            
//                         
// Modification History
// --------------------
// Use hdf_bin class rather than hdf4_bin class
// J. Gales     Futuretech     07/18/11
// =====================================================================

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <list>
#include "netcdf.h"
#include "hdf_bin.h"
#include "version.h"

#define NBINREAD MAXPIX
#define FLAGMASK PRODFAIL

#include "l12_proto.h"


using namespace std;

typedef Hdf::binListStruct blstr;

void l3gen_usage(char *prog) {
    l2gen_usage("l3gen");
}

void load_l3file_handle(int sensorID, int filenum, filehandle *file) {
    strcpy(file->name, (const char *) input->ofile[filenum]);
    strcpy(file->l2prod, (const char *) input->l2prod[filenum]);
    strcpy(file->def_l2prod, input->def_l2prod[filenum]);

    file->format = FT_L3BIN;
    file->mode = WRITE;
    file->sensorID = sensorID;

    file->tot_prod = prodlist(file->sensorID, l1_input->evalmask, file->l2prod,
            file->def_l2prod, file->l2_prod_names);
}

void load_input_filehandle(filehandle *l1file) {
    int32 sensorID = l1file->sensorID;
    int32 evalmask = l1_input->evalmask;

    l1file->nbands = rdsensorinfo(sensorID, evalmask, NULL, NULL);
    l1file->nbandsir = rdsensorinfo(sensorID, evalmask, "NbandsIR", NULL);
    l1file->ndets = 1;

    rdsensorinfo(sensorID, evalmask, "Bindx", (void **) &l1file->bindx);
    rdsensorinfo(sensorID, evalmask, "Lambda", (void **) &l1file->iwave);
    rdsensorinfo(sensorID, evalmask, "fwave", (void **) &l1file->fwave);
    rdsensorinfo(sensorID, evalmask, "Fobar", (void **) &l1file->Fobar);
    rdsensorinfo(sensorID, evalmask, "Tau_r", (void **) &l1file->Tau_r);
    rdsensorinfo(sensorID, evalmask, "k_oz", (void **) &l1file->k_oz);
    rdsensorinfo(sensorID, evalmask, "k_no2", (void **) &l1file->k_no2);
    rdsensorinfo(sensorID, evalmask, "aw", (void **) &l1file->aw);
    rdsensorinfo(sensorID, evalmask, "bbw", (void **) &l1file->bbw);

    bindex_set(l1file->iwave, l1file->nbands + l1file->nbandsir, BANDW);

    if ((l1file->Fonom = (float*) calloc(l1file->nbands, sizeof (float))) == NULL) {
        printf("-E- (%s, %d) Cannot allocate space for l1file->Fonom\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    int i;
    for (i = 0; i < l1file->nbands; i++) {
        if (l1_input->outband_opt >= 2) {
            get_f0_thuillier_ext(l1file->iwave[i], BANDW, l1file->Fonom + i);
        } else {
            l1file->Fonom[i] = l1file->Fobar[i];
        }
    }

}

void load_l12(l2str *l2rec) {
    int32 iw;

    l1str* l1rec = l2rec->l1rec;
    filehandle* l1file = l1rec->l1file;
    int32 nbands = l1file->nbands;

    //memset(l1rec->data, 0, l1rec->length);
    //memset(l2rec->data, 0, l2rec->length);

    for (iw = 0; iw < nbands; iw++) {
        l1rec->Fo[iw] = l1file->Fobar[iw]; // for the moment
    }

    l1rec->tilt = 0.0;
    l1rec->mside = 0;
    l1rec->detnum = 0;

    int ip;
    for (ip = 0; ip < l1rec->npix; ip++) {
        l1rec->pixnum[ip] = 0;
        l1rec->alpha[ip] = 0.0;
        l1rec->flags[ip] = 0;
        l1rec->mask[ip] = 0;
    }

}

// -------------------------------------------------------------------- 
//                            main                                      
// -------------------------------------------------------------------- 

int main(int argc, char* argv[]) {
    static l1str *l1rec; // generic level-1b scan structure
    static l2str *l2rec; // generic level-2  scan structure
    static filehandle ifile; // input file handle 
    static filehandle ofile; // output file handle 
    int iprod;

    char *ptime = ydhmsf(now(), 'G');
    double start_time = now();

    int16 year, day;
    int32_t syear;
    int32_t sday;
    double dsec;
    double mtime;

    char soft_id[200];
    float gmt;
    float solz;
    float sola;
    char buf[FILENAME_MAX];

    int i;
    int ip;
    int offset;
    int64_t bin_num;
    int nobs;
    int nscenes;
    int nwrite;
    static bool atLeastOne = false;
    int mainReturnCode = EXIT_SUCCESS;

    char units_string[MD_ATTRSZ];

    if (argc == 1) {
        l3gen_usage(argv[0]);
        return EXIT_SUCCESS;
    }

    // see if help on command line
    for (i = 0; i < argc; i++) {
        if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-help") == 0)) {
            l3gen_usage(argv[0]);
            return EXIT_SUCCESS;
        }
    }

    setvbuf(stdout, NULL, _IOLBF, 0);
    setvbuf(stderr, NULL, _IOLBF, 0);

    // Initialize file handles
    msl12_input_init();

    cdata_();
    filehandle_init(&ifile);
    filehandle_init(&ofile);

    // Parse input parameters
    if (msl12_input(argc, argv, "l3gen", &ifile) != 0) {
        printf("-E- %s: Error parsing input parameters.\n", argv[0]);
        exit(1);
    }

    load_input_filehandle(&ifile);

    char proc_con[2048];
    strcpy(proc_con, basename(argv[0]));
    for (i = 1; i < argc; i++) {
        strcat(proc_con, " ");
        strcat(proc_con, argv[i]);
    }

    // Transfer input info to file handles
    load_l3file_handle(ifile.sensorID, 0, &ofile);

    // Build output product list

    cout << endl << "The following products will be included in " << ofile.name
            << endl;
    int32 nprods_out = ofile.tot_prod;
    char prodlist_out[L1_PRODSTRLEN] = "";
    for (iprod = 0; iprod < ofile.tot_prod; iprod++) {
        cout << iprod + 1 << " " << ofile.l2_prod_names[iprod] << endl;
        strcat(prodlist_out, ofile.l2_prod_names[iprod]);
        if (iprod < ofile.tot_prod - 1)
            strcat(prodlist_out, ",");
    }

    // Open input binfile and get full input product list

    static Hdf::hdf_bin *input_binfile;
    const char *inputfile;

    inputfile = input->ifile[0];
    input_binfile = Hdf::openBinObject(inputfile);

    int len = input_binfile->query();
    char *fullprodlist = (char *) malloc(len);
    input_binfile->query(fullprodlist);
    int nrows = input_binfile->nrows;

    // Get input bin dimension
    int32 npix = NBINREAD;

    // Get mean time
    mtime = (input_binfile->meta_l3b.startTime + input_binfile->meta_l3b.endTime) / 2.0;

    unix2yds(mtime, &year, &day, &dsec);
    //yes, silly, but easier than fixing the various routines that need either a 16 or 32bit int for year/day
    syear = year;
    sday = day;

    // Allocate memory for L1 and L2 scan data (the L2 record shares
    // space with the L1 record, hence we need to allocate both)
    l1rec = (l1str*) malloc(sizeof (l1str));
    l2rec = (l2str*) malloc(sizeof (l2str));

    ifile.npix = npix;
    ifile.nscan = 1;
    if (alloc_l1(&ifile, l1rec) == 0) {
        printf("-E- %s: Unable to allocate L1 record.\n", argv[0]);
        exit(1);
    }
    if (alloc_l2(l1rec, l2rec) == 0) {
        printf("-E- %s: Unable to allocate L2 record.\n", argv[0]);
        exit(1);
    }

    // Add meta-data to L1 and L2 structures - use a single pixel just to get things
    // started - load_l12 will be called for each bin row with an appropriate # of pixels
    load_l12(l2rec);

    // Build input product request list

    int32 n_vvv = MIN(windex(700., ifile.fwave, ifile.nbands) + 1, ifile.nbands);
    char prodstr[L1_PRODSTRLEN] = "";

    //   what we have

    char **prodlist_in;
    int32 nprods_in = input_binfile->query(&prodlist_in);

    //   what we want

    char prodlist_wt[L1_MAXPROD][32];
    int32 nprods_wt = prodlist(ifile.sensorID, l1_input->evalmask,
            "Rrs_vvv nLw_vvv Rrs_vc_vvv solz sola senz sena relaz", "",
            prodlist_wt);

    enum geo_order {
        SOLZ, SOLA, SENZ, SENA, RELAZ
    };

    //   what we have that we want

    int32 *rrs_bindex;
    int32 *rrs_windex;
    int32 n_rrs = 0;
    int32 *nlw_bindex;
    int32 *nlw_windex;
    int32 n_nlw = 0;
    int32 geo_bindex[5];
    int32 geo_windex[5];
    int32 n_geo = 0;
    int32 have_solz = 0;
    int32 have_sola = 0;
    int32 have_sena = 0;
    int32 have_relaz = 0;

    if ((rrs_bindex = (int32 *) malloc(ifile.nbands * sizeof (int32))) == NULL) {
        printf("-E- %s: Error allocating memory to the input data array.\n",
                argv[0]);
        exit(FATAL_ERROR);
    }
    if ((rrs_windex = (int32 *) malloc(ifile.nbands * sizeof (int32))) == NULL) {
        printf("-E- %s: Error allocating memory to the input data array.\n",
                argv[0]);
        exit(FATAL_ERROR);
    }
    if ((nlw_bindex = (int32 *) malloc(ifile.nbands * sizeof (int32))) == NULL) {
        printf("-E- %s: Error allocating memory to the input data array.\n",
                argv[0]);
        exit(FATAL_ERROR);
    }
    if ((nlw_windex = (int32 *) malloc(ifile.nbands * sizeof (int32))) == NULL) {
        printf("-E- %s: Error allocating memory to the input data array.\n",
                argv[0]);
        exit(FATAL_ERROR);
    }

    cout << endl << "Checking input for desired products: " << endl;
    iprod = 0;
    for (int32 iwant = 0; iwant < nprods_wt; iwant++) {
        for (int32 ihave = 0; ihave < nprods_in; ihave++) {
            if (strcmp(prodlist_in[ihave], prodlist_wt[iwant]) == 0) {
                strcat(prodstr, prodlist_in[ihave]);
                strcat(prodstr, ",");
                cout << "found " << prodlist_wt[iwant] << endl;
                if (iwant < n_vvv) {
                    rrs_bindex[n_rrs] = iprod; // posiiton in bin read
                    rrs_windex[n_rrs] = iwant; // position in sensor wavelengths
                    n_rrs++;
                    iprod++;
                } else if (iwant < n_vvv * 2) {
                    nlw_bindex[n_nlw] = iprod; // positon in bin read
                    nlw_windex[n_nlw] = iwant - n_vvv; // position in sensor wavelengths
                    n_nlw++;
                    iprod++;
                } else if (iwant < n_vvv * 3) { // Rrs_vc products overwriting Rrs
                    rrs_bindex[n_rrs] = iprod; // posiiton in bin read
                    rrs_windex[n_rrs] = iwant - 2 * n_vvv; // position in sensor wavelengths
                    n_rrs++;
                    iprod++;
                } else {
                    geo_bindex[n_geo] = iprod;
                    geo_windex[n_geo] = iwant - 3 * n_vvv;
                    switch (geo_windex[n_geo]) {
                    case SOLZ:
                        have_solz = 1;
                        break;
                    case SOLA:
                        have_sola = 1;
                        break;
                    case SENA:
                        have_sena = 1;
                        break;
                    case RELAZ:
                        have_relaz = 1;
                        break;
                    }
                    n_geo++;
                    iprod++;
                }
                break;
            }
        }
    }
    nprods_in = iprod;

    // remove trailing comma
    len = strlen(prodstr);
    prodstr[len - 1] = '\0';

    // Allocate data arrays and bin list and masking flag
    float **outData;
    float **tmpData;
    float **inData;
    if ((inData = (float **) malloc(nprods_in * sizeof (float *))) == NULL) {
        printf("-E- %s: Error allocating memory to the input data array.\n",
                argv[0]);
        exit(FATAL_ERROR);
    }
    if ((outData = (float **) malloc(nprods_out * sizeof (float *))) == NULL) {
        printf("-E- %s: Error allocating memory to the output data array.\n",
                argv[0]);
        exit(FATAL_ERROR);
    }
    if ((tmpData = (float **) malloc(nprods_out * sizeof (float *))) == NULL) {
        printf("-E- %s: Error allocating memory to the temp data array.\n",
                argv[0]);
        exit(FATAL_ERROR);
    }
    for (i = 0; i < nprods_in; i++) {
        if ((inData[i] = (float *) calloc(2 * npix, sizeof (float))) == NULL) {
            printf("-E- %s: Error allocating memory to the output data array.\n",
                    argv[0]);
            exit(FATAL_ERROR);
        }
    }
    for (i = 0; i < nprods_out; i++) {
        if ((outData[i] = (float *) calloc(2 * npix, sizeof (float))) == NULL) {
            printf(
                    "-E- %s: Error allocating memory to the output data array.\n",
                    argv[0]);
            exit(FATAL_ERROR);
        }
        if ((tmpData[i] = (float *) calloc(npix, sizeof (float))) == NULL) {
            printf(
                    "-E- %s: Error allocating memory to the output data array.\n",
                    argv[0]);
            exit(FATAL_ERROR);
        }
    }
    char *outMask = (char *) calloc(npix, sizeof (char));

    // set up the bin file reading for the products we need
    input_binfile->read(prodstr);

    /* Create output file */
    /* ------------------ */
    Hdf::hdf_bin *output_binfile;

    /*
     * If the input structure does not set the output format,
     * make the output format match the input L3 format
     */
    if (getFileFormatName(input->oformat) == NULL) {
        if (input_binfile->isHDF5) {
            strcpy(input->oformat, "HDF5");
        } else if (input_binfile->isCDF4) {
            strcpy(input->oformat, "netCDF4");
        } else {
            strcpy(input->oformat, "HDF4");
        }
    }

    if (strcmp(input->oformat, "HDF4") == 0) {
        output_binfile = new Hdf::hdf4_bin;
        output_binfile->hasNoext = true;
    }
    if (strcmp(input->oformat, "HDF5") == 0)
        output_binfile = new Hdf::hdf5_bin;
    if (strcmp(input->oformat, "netCDF4") == 0)
        output_binfile = new Hdf::cdf4_bin;

    output_binfile->deflate = input->deflate;

    strcpy(output_binfile->meta_l3b.product_name, ofile.name);
    strncpy(output_binfile->meta_l3b.ptime, ptime, 16);
    strcpy(output_binfile->meta_l3b.proc_con, proc_con);
    strcpy(output_binfile->meta_l3b.input_parms, l1_input->input_parms);

    char* tmpProdNames[ofile.tot_prod];
    for (int i = 0; i < ofile.tot_prod; i++) {
        tmpProdNames[i] = ofile.l2_prod_names[i];
    }
    output_binfile->setProductList(ofile.tot_prod, tmpProdNames);
    output_binfile->create(ofile.name, input_binfile->nrows);

    // Begin processing
    for (int32 iscan = 0; iscan < nrows; iscan++) {
        if (iscan % 50 == 0) {
            int secs = (int) now() - start_time;
            cout << "Reading row " << iscan << " of " << nrows << " after " << secs << " seconds" << endl;
        }
        // Read the inputs
        // Get basebin and numbin for this input row
        int64_t basebin = input_binfile->get_basebin(iscan);
        input_binfile->readBinIndex(iscan);
        int ext = input_binfile->get_ext();
        int64_t beg = input_binfile->get_beg();
        // if the row has no filled bins, skip it
        if (beg == 0)
            continue;

        input_binfile->readBinList(ext);

        // fill the input data array with the necessary input data
        i = 0;
        for (int j = 0; j < input_binfile->nprod(); j++) {

            if (input_binfile->active_data_prod[j] == true) {
                input_binfile->get_prodname(j, buf);
                if (strcmp(prodlist_in[j], buf) == 0) {
                    input_binfile->readSums(inData[i], ext, j);
                    i++;
                }
            }
        }
        if (input_binfile->isHDF5 || input_binfile->isCDF4)
            input_binfile->setDataPtr(ext);


        // Stuff the L2 record
        // only process as many pixels as the row as filled bins
        l1rec->npix = ext;
        ifile.npix = ext;
        ifile.epix = ext - 1;
        ifile.terrain_corrected = 1;

        memset(outMask, '\0', npix);
        init_l1(l1rec);
        init_l2(l2rec, ifile.nbands);

        load_l12(l2rec);

        l1rec->iscan = iscan;
        l1rec->scantime = mtime;

        // Fill in the l1/2 structures with data from input bin file
        // nLw, Rrs, geometries - if available.
        for (ip = 0; ip < ext; ip++) {
            float weight = input_binfile->get_weights(ip);
            int32 ipw;
            float lat;
            float lon;
            float flat, flon;

            bin_num = input_binfile->get_bin_num(ip);
            input_binfile->bin2latlon(bin_num, lat, lon);
            l1rec->lon[ip] = lon;
            l1rec->lat[ip] = lat;
            l1rec->nobs[ip] = input_binfile->get_nobs(ip);
            for (int32 iw = 0; iw < n_nlw; iw++) {
                ipw = ip * ifile.nbands + nlw_windex[iw];
                l2rec->nLw[ipw] = inData[nlw_bindex[iw]][2 * ip] / weight;
                l2rec->Rrs[ipw] = MAX(l2rec->nLw[ipw] / ifile.Fonom[iw], BAD_FLT);
            }
            for (int32 iw = 0; iw < n_rrs; iw++) {
                ipw = ip * ifile.nbands + rrs_windex[iw];
                l2rec->Rrs[ipw] = inData[rrs_bindex[iw]][2 * ip] / weight;
                l2rec->nLw[ipw] = MAX(l2rec->Rrs[ipw] * ifile.Fonom[iw], BAD_FLT);
            }
            l1rec->solz[ip] = 0.0;
            l1rec->sola[ip] = 0.0;
            l1rec->senz[ip] = 0.0;
            l1rec->sena[ip] = 90.0;
            l1rec->delphi[ip] = 90.0;

            for (int32 iw = 0; iw < n_geo; iw++) {
                switch (geo_windex[iw]) {
                case SOLZ:
                    l1rec->solz[ip] = inData[geo_bindex[iw]][2 * ip] / weight;
                    break;
                case SOLA:
                    l1rec->sola[ip] = inData[geo_bindex[iw]][2 * ip] / weight;
                    break;
                case SENZ:
                    l1rec->senz[ip] = inData[geo_bindex[iw]][2 * ip] / weight;
                    break;
                case SENA:
                    l1rec->sena[ip] = inData[geo_bindex[iw]][2 * ip] / weight;
                    break;
                case RELAZ:
                    l1rec->delphi[ip] = inData[geo_bindex[iw]][2 * ip] / weight;
                    break;
                }
            }
            if (!have_solz) {
                // assume noon orbit
                flon = l1rec->lon[ip];
                flat = l1rec->lat[ip];
                gmt = 12 - lon / 15;
                sunangs_(&syear, &sday, &gmt, &flon, &flat, &solz, &sola);
                l1rec->solz[ip] = solz;
                l1rec->sola[ip] = sola;
            }
            if (have_sola && have_sena && !have_relaz) {
                l1rec->delphi[ip] = l1rec->sena[ip] - 180.0 - l1rec->sola[ip];
                if (l1rec->delphi[ip] < -180)
                    l1rec->delphi[ip] += 360.0;
            } else if (!have_sola && !have_sena && have_relaz) {
                l1rec->sena[ip] = 90.0;
                l1rec->sola[ip] = l1rec->sena[ip] - 180.0
                        - l1rec->delphi[ip];
                if (l1rec->sola[ip] < -180)
                    l1rec->sola[ip] += 360.0;
            }
        }

        // Add ancillary, Rayleigh, etc.
        loadl1(&ifile, l1rec);

        // clear out the masks
        for (ip = 0; ip < ext; ip++)
            l1rec->mask[ip] = 0;

        // Add a default chl

        for (ip = 0; ip < ext; ip++)
            l2rec->chl[ip] = get_default_chl(l2rec, &l2rec->Rrs[ip * ifile.nbands]);

        // Add default inherent optical properties

        if (input->iop_opt > 0 && (input->proc_ocean != 0))
            get_iops(l2rec, input->iop_opt);

        // Compute output products and store in temporary buffer

        for (iprod = 0; iprod < nprods_out; iprod++) {

            // get the product index record
            for (i = 0; i < ext; i++)
                tmpData[iprod][i] = BAD_FLT;

            l2prodstr *p;

            if ((p = get_l2prod_index(ofile.l2_prod_names[iprod],
                    ifile.sensorID, ifile.nbands + ifile.nbandsir, ext,
                    ifile.nscan, ifile.iwave)) == NULL) {
                printf("-E- %s line %d: product index failure.\n", __FILE__,
                        __LINE__);
                exit(1);
            };

            // Compute or extract the product & copy to output buffer

            VOIDP pbuf = prodgen(p, l2rec);

            // build the units attribute string

            productInfo_t *p_info;
            p_info = allocateProductInfo();

            if (!findProductInfo(ofile.l2_prod_names[iprod], ifile.sensorID, p_info)) {
                printf("-E- product %s not found in XML product table\n",
                        ofile.l2_prod_names[iprod]);
                exit(EXIT_FAILURE);
            }

            if (iprod == 0){
                strcpy(units_string,getProductNameFull(p_info));
                strcat(units_string,":");
                strcat(units_string,p->units);
            }
            else {
                strcat(units_string,getProductNameFull(p_info));
                strcat(units_string,":");
                strcat(units_string,p->units);
            }
            if (iprod < nprods_out -1)
                strcat(units_string,",");

            freeProductInfo(p_info);

            memcpy(tmpData[iprod], (float *) pbuf,
                    ext * sizeof (float));

            // Check flags and set masking

            for (ip = 0; ip < ext; ip++) {
                if ((tmpData[iprod][ip] == BAD_FLT)
                        || (l1rec->flags[ip] & FLAGMASK) != 0) {
                    outMask[ip] = 1;
                }
            }

        }

        /* Write output */
        /* ------------ */
        output_binfile->clear_binlist();
        nwrite = 0;

        for (ip = 0; ip < ext; ip++) {
            if (!outMask[ip]) {
                atLeastOne = true;
                bin_num = input_binfile->get_bin_num(ip);
                offset = bin_num - basebin;
                output_binfile->set_bin_num(offset, bin_num);
                nobs = input_binfile->get_nobs(ip);
                output_binfile->inc_nobs(offset, nobs);
                nscenes = input_binfile->get_nscenes(ip);
                output_binfile->inc_nscenes(offset, nscenes);
                output_binfile->set_weights(offset, 1);

                /* Loop over data products */
                /* ----------------------- */
                for (iprod = 0; iprod < nprods_out; iprod++) {
                    outData[iprod][2 * nwrite] = tmpData[iprod][ip];

                } /* iprod loop */
                if (nwrite != offset)
                    output_binfile->copy_binlist(offset, nwrite);

                nwrite++;
            }
        } /* ip loop */


        /* Write BinList & Data Products */
        /* ----------------------------- */
        if (nwrite > 0) {
            output_binfile->writeBinList(nwrite);
            for (iprod = 0; iprod < nprods_out; iprod++) {
                strcpy(buf, ofile.l2_prod_names[iprod]);
                output_binfile->writeSums(outData[iprod], nwrite, buf);
            }

            if (strcmp(input->oformat, "HDF5") == 0
                    || strcmp(input->oformat, "netCDF4") == 0)
                output_binfile->incNumRec(nwrite);
        }
    }

    output_binfile->copymeta(1, &input_binfile);
    if (strcmp(input->oformat, "netCDF4")) {
        strcpy(output_binfile->meta_l3b.sensor, sensorId2InstrumentName(ifile.sensorID));
    } else {
        strcpy(output_binfile->meta_l3b.sensor_name, sensorId2SensorName(ifile.sensorID));
    }
    strcpy(output_binfile->meta_l3b.mission, sensorId2PlatformName(ifile.sensorID));
    strcpy(output_binfile->meta_l3b.infiles, basename((char *) input->ifile[0]));
    strcpy(output_binfile->meta_l3b.soft_name, "l3gen");
    sprintf(soft_id, "%d.%d.%d-%s", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH, GITSHA);
    strcpy(output_binfile->meta_l3b.soft_ver, soft_id);
    if (strcmp(input->oformat, "netCDF4") == 0) {
        ptime = unix2isodate(now(), 'G');
        strcpy(output_binfile->meta_l3b.ptime, ptime);
    } else {
        strcpy(output_binfile->meta_l3b.ptime, ptime);

    }
    strcpy(output_binfile->meta_l3b.proc_con, proc_con);
    //update units attribute to match the output products
    strcpy(output_binfile->meta_l3b.units,units_string);
    // Some missions have suites with DOIs and keywords...

    strcpy(output_binfile->meta_l3b.doi, input->doi);
    //if suite is somehow NOT set, clear it out as it's better than lying
    // since we copy from the source
    strcpy(output_binfile->meta_l3b.keywords, "");
    if (input->suite[0]) {
        const char* keywordStr = getGCMDKeywords(input->suite);
        if (keywordStr) {
             strcpy(output_binfile->meta_l3b.keywords, keywordStr);
        }
    }

    // Close files and free memory
    if (atLeastOne == false) {
        cout << "No valid bins to output!" << endl << "...seems a selected product may have resulted in 100% PRODFAIL..." << endl;
        mainReturnCode = 110;  // no bins filled
    }
    output_binfile->close();
    delete(output_binfile);
    input_binfile->close();
    delete(input_binfile);

    free(fullprodlist);

    for (i = 0; i < nprods_in; i++)
        free(inData[i]);
    free(inData);

    free(outMask);

    for (i = 0; i < nprods_out; i++) {
        free(outData[i]);
        free(tmpData[i]);
    }
    free(outData);
    free(tmpData);

    free_l1(l1rec);
    free_l2(l2rec);
    free(l1rec);
    free(l2rec);
    // free input internals
    for(i=0;i<input->fctl.nfilt; i++)
        free(input->fctl.f[i].kernel);
    free(input->gsm_aphs);
    free(input->gsm_aphw);
    free(input->giop_wave);
    free(input->giop_rrs_unc);
    free(input->taua);
    free(input->vcal_nLw);
    free(input->vcal_Lw);
    free(input);

    free(rrs_bindex);
    free(rrs_windex);
    free(nlw_bindex);
    free(nlw_windex);

    cout << "Processing Complete at " << ptime << endl << endl;

    return (mainReturnCode);
}

