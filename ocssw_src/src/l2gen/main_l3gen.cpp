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
#include "geo_region.h"
#include "get_dataday.h"
#include <unordered_map>
#include <unordered_set>
#include <regex>
#define NBINREAD MAXPIX
#define FLAGMASK (PRODFAIL | NAVFAIL)

#include "L3File.h"
#include "l12_proto.h"


using namespace std;

typedef Hdf::binListStruct blstr;

void l3gen_usage(char *prog) {
    clo_optionList_t* list;
    list = clo_createList();
    l2gen_init_options(list, prog);
    clo_addOption(list, "averaging_scheme", CLO_TYPE_STRING, NULL, "Averaging scheme to apply when averaging data\n        Usage - 'averaging_scheme=prod1:method1,prod2:method2,...'\n"
                  "        Averaging method specification:\n                arithmetic: arithmetic mean\n        "
                  "        geometric: geometric mean\n                harmonic: harmonic mean\n        "
                  "        quadratic: quadratic mean\n");
    // add averaging scheme
    clo_printUsage(list);
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

    int status = set_solar_irradiance(l1_input->f0file);
    int i;
    for (i = 0; i < l1file->nbands; i++) {
        if (l1_input->outband_opt >= 2) {
            if (status) {
                printf("-E- %s:%d : error reading solar_irradiance LUT \"%s\".\n", __FILE__, __LINE__, l1_input->f0file);
                exit(EXIT_FAILURE);
            }
            status = get_f0(l1file->iwave[i],BANDW,l1file->Fonom + i);
            if (status) {
                printf("-E- %s:%d - Can't calculate solar irradiance.\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
        } else {
            l1file->Fonom[i] = l1file->Fobar[i];
        }
    }

}

static bool matches_prefix_int_suffix_regex(const std::string& prefix, const std::string& s) {
  const std::regex re("^" + prefix + "_[1-9][0-9]*$");    
  return std::regex_match(s, re);
  
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

/**
 * Determine the averaging scheme to be applied for each output product.  If "ALL" is specified, then the same scheme is applied to all products.  Otherwise, the scheme must be specified for each product.  If a product is not specified, then the default of arithmetic mean is applied.
 * @param ofile Output file handle containing the list of products to be included in the output
 * @param output_averaging_schemes Output vector of averaging schemes to be applied for each product
 * @param output_averaging_schemes_names Output vector of averaging scheme names to be applied for each product
 * @param averaging_schemes Map of product name to averaging scheme and averaging scheme name specified on the command line
 * @return true on success, false if a product specified in the averaging_schemes map is not found in the output file's product list
 */
static bool get_output_averaging_schemes(
    filehandle *ofile, std::vector<AveragingScheme> &output_averaging_schemes, std::vector<std::string> &output_averaging_schemes_names,
    const std::unordered_map<std::string, std::pair<AveragingScheme, std::string>> &averaging_schemes) {
    if (averaging_schemes.count("ALL")) {
        std::fill(output_averaging_schemes.begin(), output_averaging_schemes.end(),
                  averaging_schemes.at("ALL").first);
        std::fill(output_averaging_schemes_names.begin(), output_averaging_schemes_names.end(), averaging_schemes.at("ALL").second);
        return true;
    }
    // check that all specified products are found
    std::unordered_set<std::string> requested_products;
    for (int iprod = 0; iprod < ofile->tot_prod; iprod++) {
        requested_products.insert(ofile->l2_prod_names[iprod]);
    }
    for (const auto &product : averaging_schemes) {
        if (requested_products.count(product.first) == 0) {
            return false;
        }
    }
    for (int iprod = 0; iprod < ofile->tot_prod; iprod++) {
        if (averaging_schemes.count(ofile->l2_prod_names[iprod])) {
            output_averaging_schemes.at(iprod) = averaging_schemes.at(ofile->l2_prod_names[iprod]).first;
            output_averaging_schemes_names.at(iprod) = averaging_schemes.at(ofile->l2_prod_names[iprod]).second;
        }
        else {
            output_averaging_schemes.at(iprod) = ARITHMETIC_MEAN;
            output_averaging_schemes_names.at(iprod) = "arithmetic";
        }
    }
    return true;
}
// -------------------------------------------------------------------- 
//                            main                                      
// -------------------------------------------------------------------- 
class ProductMatching {
   public:
    productInfo_t* productInfo = nullptr;
    int32_t sensorID{-1};
    size_t nprods_total{0};
    bool found_solz{false};
    bool found_sola{false};
    bool found_sena{false};
    bool found_senz{false};
    bool found_relaz{false};
    size_t index_solz{}, index_senz{}, index_sena{}, index_sola{}, index_relaz{};

    class ThreeDimsVariable
    {
        public:
        std::vector<int32_t> match_l3b_index{};
        std::vector<int32_t> bindex_found{};
        std::vector<int32_t> wave_found{};
        std::unordered_map<std::string, int32_t> wavelengths{};
        std::unordered_map<std::string, int32_t> bindex{};

    };

    void set_requested_3d_products(ThreeDimsVariable& threeDimsVariable, const std::string& prodname_requested) {
        int wave = get_wavelenth(prodname_requested);
        if (wave != -1) {
            int index = bindex_get(wave);
            if (index < 0) {
                fprintf(stderr, "-E-: %s:%d wavelength %d not found for %s or bindex has not been set\n",
                        __FILE__, __LINE__, wave, prodname_requested.c_str());
                exit(EXIT_FAILURE);
            }
            threeDimsVariable.wavelengths[prodname_requested] = wave;
            threeDimsVariable.bindex[prodname_requested] = index;
        }
    }

    void match_products_3d(ThreeDimsVariable& threeDimsVariable, const std::string& prodname_requested,
                           size_t iprod_provided) {
        if (threeDimsVariable.wavelengths.count(prodname_requested)) {
            int32_t wave = threeDimsVariable.wavelengths.at(prodname_requested);
            int32_t bindex = threeDimsVariable.bindex.at(prodname_requested);
            threeDimsVariable.bindex_found.push_back(bindex);
            threeDimsVariable.match_l3b_index.push_back(iprod_provided);
            threeDimsVariable.wave_found.push_back(wave);
        }
    }

    ThreeDimsVariable rrs_var{};
    ThreeDimsVariable rrs_vc_var{};
    ThreeDimsVariable rrs_unc_var{};
    ThreeDimsVariable rhos_var{};
    ThreeDimsVariable nlw_var{};
    
    // active products 
    std::string prodstring{};
    std::vector<size_t> active_products_index{};
    std::unordered_map<std::string, size_t> additional_prodname_to_index{};
    explicit ProductMatching(char prodlist_requested[L1_MAXPROD][32], char** prodlist_provided,
                             size_t number_of_product_requested, size_t number_of_products_provided,
                             int32_t sensorID) {
        this->sensorID = sensorID;
        productInfo = allocateProductInfo();
        // first, determine all the sets of wavelengths requested by the user
         std::unordered_set<std::string> requested_products_set;
        for (size_t iprod = 0; iprod < number_of_product_requested; ++iprod) {
            std::string prodname_requested = prodlist_requested[iprod];
            if(requested_products_set.count(prodname_requested))
                continue;
            requested_products_set.insert(prodname_requested);
            if (matches_prefix_int_suffix_regex("Rrs_vc", prodname_requested)) {
                set_requested_3d_products(rrs_vc_var, prodname_requested);
            }
            if (matches_prefix_int_suffix_regex("Rrs", prodname_requested)) {
                set_requested_3d_products(rrs_var, prodname_requested);
            }
            if (matches_prefix_int_suffix_regex("Rrs_unc", prodname_requested)) {
                set_requested_3d_products(rrs_unc_var, prodname_requested);
            }
            if (matches_prefix_int_suffix_regex("nLw", prodname_requested)) {
                set_requested_3d_products(nlw_var, prodname_requested);
            }
            if (matches_prefix_int_suffix_regex("rhos", prodname_requested)) {
                set_requested_3d_products(rhos_var, prodname_requested);
            }
        }

        // parse through all products
        for (size_t iprod_provided = 0; iprod_provided < number_of_products_provided; ++iprod_provided) {
            std::string prodname_provided = prodlist_provided[iprod_provided];
            if (requested_products_set.count(prodname_provided)) {
                match_products_3d(rrs_vc_var, prodname_provided, nprods_total);
                match_products_3d(rrs_var, prodname_provided, nprods_total);
                match_products_3d(rrs_unc_var, prodname_provided, nprods_total);
                match_products_3d(rhos_var, prodname_provided, nprods_total);
                match_products_3d(nlw_var, prodname_provided, nprods_total);
                if (prodname_provided == "solz") {
                    found_solz = 1;
                    index_solz = nprods_total;
                } else if (prodname_provided == "senz") {
                    found_senz = 1;
                    index_senz = nprods_total;
                } else if (prodname_provided == "sena") {
                    found_sena = 1;
                    index_sena = nprods_total;
                } else if (prodname_provided == "sola") {
                    found_sola = 1;
                    index_sola = nprods_total;
                } else if (prodname_provided == "relaz") {
                    found_relaz = 1;
                    index_relaz = nprods_total;
                } else {
                    additional_prodname_to_index[prodname_provided] = nprods_total;
                }
                active_products_index.push_back(iprod_provided);
                if (prodstring.empty()) {
                    prodstring = prodname_provided;
                } else {
                    prodstring += "," + prodname_provided;
                }
                nprods_total++;
            }
        }
    }
    int32_t get_wavelenth(const std::string& product_name) {
        int32_t res = findProductInfo(product_name.c_str(), sensorID, productInfo);
        if (res != 1) {
            return -1;
        }
        if (std::string(productInfo->paramDesignator) != "none") {
            int32_t wave = productInfo->prod_ix;
            if (wave > 0)
                return wave;
            else {
                fprintf(stderr, "-E-: %s:%d wavelength not found for %s\n", __FILE__, __LINE__,
                        product_name.c_str());
                exit(EXIT_FAILURE);
            }
        }
        return -1;
    }
    bool validate_matched_products(const ThreeDimsVariable& threeDimsVariable,
                                   const filehandle* l1file) const {
        for (size_t iw = 0; iw < threeDimsVariable.bindex_found.size(); iw++) {
            size_t band_index = threeDimsVariable.bindex_found[iw];
            int32_t wave = threeDimsVariable.wave_found[iw];
            if (l1file->iwave[band_index] != wave) {
                return false;
            }
        }
        return true;
    }

    void validate_matched_products(const filehandle* l1file, Hdf::hdf_bin* input_binfile) const {
        bool validation_status = true;
        bool status = validate_matched_products(rrs_var, l1file);
        validation_status &= status;
        if (!status) {
            fprintf(stderr, "-E-: %s:%d Wrong wavelength for Rrs \n", __FILE__, __LINE__);
        }
        status = validate_matched_products(rrs_unc_var, l1file);
        validation_status &= status;
        if (!status) {
            fprintf(stderr, "-E-: %s:%d Wrong wavelength for Rrs_unc \n", __FILE__, __LINE__);
        }
        status = validate_matched_products(rrs_vc_var, l1file);
        validation_status &= status;
        if (!status) {
            fprintf(stderr, "-E-: %s:%d Wrong wavelength for Rrs_vc \n", __FILE__, __LINE__);
        }
        status = validate_matched_products(rhos_var, l1file);
        validation_status &= status;
        if (!status) {
            fprintf(stderr, "-E-: %s:%d Wrong wavelength for rhos \n", __FILE__, __LINE__);
        }
        status = validate_matched_products(nlw_var, l1file);
        validation_status &= status;
        if (!status) {
            fprintf(stderr, "-E-: %s:%d Wrong wavelength for nLw \n", __FILE__, __LINE__);
        }
        if (input_binfile->n_active_prod != (int32_t)nprods_total) {
            fprintf(stderr,
                    "-E-: %s:%d number of active products in input bin file (%d) does not match number of "
                    "products found in input bin file (%ld)\n",
                    __FILE__, __LINE__, input_binfile->n_active_prod, nprods_total);
            validation_status = false;
        }
        size_t prod_counter = 0;
        for (size_t j = 0; j < (size_t)input_binfile->nprod(); j++) {
            if (input_binfile->active_data_prod[j]) {
                if (active_products_index.at(prod_counter) != j) {
                    fprintf(stderr, "-E-: %s:%d active product index mismatch: expected %ld but found %ld\n",
                            __FILE__, __LINE__, active_products_index.at(prod_counter), j);
                    validation_status = false;
                }
                prod_counter++;
            }
        }
        if (!validation_status) {
            exit(EXIT_FAILURE);
        }
    }

    ~ProductMatching() {
        if(productInfo) {
            freeProductInfo(productInfo);
        }
    }
    // we need an array for requested Rrr
    // int - order, -> get int
};

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

    // cdata_();
    filehandle_init(&ifile);
    filehandle_init(&ofile);

    // Parse input parameters
    if (msl12_input(argc, argv, "l3gen", &ifile) != 0) {
        printf("-E- %s: Error parsing input parameters.\n", argv[0]);
        exit(1);
    }

    load_input_filehandle(&ifile);

    if(input->georegionfile[0]) {
        if (access(input->georegionfile, F_OK) || access(input->georegionfile, R_OK)) {
            printf("-E-  georegionfile '%s' does not exist or cannot be opened.\n",
                    input->georegionfile);
            exit(FATAL_ERROR);
        }
        set_georegion_filename(input->georegionfile);
    }

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

    char prodstr[L1_PRODSTRLEN] = "";

    //   what we have

    char **prodlist_in;
    int32 nprods_in = input_binfile->query(&prodlist_in);

    //   what we want

    char prodlist_wt[L1_MAXPROD][32];
    int32 nprods_wt = prodlist(ifile.sensorID, l1_input->evalmask,
            "Rrs_vvv rhos_vvv Rrs_unc_vvv nLw_vvv Rrs_vc_vvv solz sola senz sena relaz", "",
            prodlist_wt);

    //   what we have that we want

    int32 n_rrs = 0;
    int32 n_nlw = 0;
    size_t n_rrs_vc = 0;
    size_t n_rrs_unc = 0;
    size_t n_rhos = 0;
    int32 have_solz = 0;
    int32 have_sola = 0;
    int32 have_sena = 0;
    int32 have_relaz = 0;
    int32 have_senz = 0;
    cout << endl << "Checking input for desired products: " << endl;
    float equatorialCrossingTime;
    int32_t plusDay;
    getEquatorCrossingTime(ifile.sensorID, false, mtime, &equatorialCrossingTime, &plusDay);
    const ProductMatching productMatching(prodlist_wt, prodlist_in, nprods_wt, nprods_in, ifile.sensorID);
    nprods_in = productMatching.nprods_total;
    n_nlw = productMatching.nlw_var.bindex_found.size();
    n_rrs = productMatching.rrs_var.bindex_found.size();
    n_rrs_unc = productMatching.rrs_unc_var.bindex_found.size();
    n_rrs_vc = productMatching.rrs_vc_var.bindex_found.size();
    n_rhos = productMatching.rhos_var.bindex_found.size();
    have_solz = productMatching.found_solz;
    have_sola = productMatching.found_sola;
    have_sena = productMatching.found_sena;
    have_relaz = productMatching.found_relaz;
    have_senz = productMatching.found_senz;
    strncpy(prodstr, productMatching.prodstring.c_str(), (sizeof prodstr) - 1 );
    // if n_rrs_unc > 0, then allocate Rrs_unc array in L2 record
    if (n_rrs_unc > 0 && l2rec->Rrs_unc == NULL) {
        l2rec->Rrs_unc = (float*)malloc(ifile.nbands * npix * sizeof(float));
        if (l2rec->Rrs_unc == NULL) {
            printf("-E- %s: Unable to allocate Rrs_unc array in L2 record.\n", argv[0]);
            exit(1);
        }
        // set to BAD_FLT to indicate not set
        for (int32_t ip = 0; ip < ifile.nbands * npix; ip++) {
            l2rec->Rrs_unc[ip] = BAD_FLT;
        }
    }
    // if n_rhos > 0, then allocate rhos array in L1 record
    if (n_rhos > 0 && l2rec->l1rec->rhos == NULL) {
        l2rec->l1rec->rhos = (float*)malloc(ifile.nbands * npix * sizeof(float));
        if (l2rec->l1rec->rhos == NULL) {
            printf("-E- %s: Unable to allocate rhos array in l1 record.\n", argv[0]);
            exit(1);
        }
        // set to BAD_FLT to indicate not set
        for (int32_t ip = 0; ip < ifile.nbands * npix; ip++) {
            l2rec->l1rec->rhos[ip] = BAD_FLT;
        }
    }

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
    std::unordered_map<std::string, std::pair<AveragingScheme, std::string>> averaging_schemes =  get_averaging_scheme_per_product(input->averaging_scheme);
    std::vector<AveragingScheme> output_averaging_schemes(ofile.tot_prod);
    std::vector<std::string> output_averaging_schemes_names(ofile.tot_prod);
    bool status = get_output_averaging_schemes(&ofile,output_averaging_schemes,output_averaging_schemes_names,averaging_schemes);
    if (!status) {
        fprintf(stderr,"-E- %s:%d Error setting output averaging schemes\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    char* tmpProdNames[ofile.tot_prod];
    for (int i = 0; i < ofile.tot_prod; i++) {
        tmpProdNames[i] = ofile.l2_prod_names[i];
    }
    // we need to add attributes here
    // determine output averaging_scheme
    char(*tmpAttVals[ofile.tot_prod])[MAXNUMATTR][MAXATTRLEN];
    char(*tmpAttNames[ofile.tot_prod])[MAXNUMATTR][MAXATTRNAMELEN];
    int(*tmpAttTypes[ofile.tot_prod])[MAXNUMATTR];
    size_t(*tmpAttSize[ofile.tot_prod])[MAXNUMATTR];
    int tempAttNums[ofile.tot_prod];
    // set just one attribute
    for (int i_prod = 0; i_prod < ofile.tot_prod; i_prod++) {
        tempAttNums[i_prod] = 1;
        tmpAttNames[i_prod] = ( char(*)[MAXNUMATTR][MAXATTRNAMELEN]) malloc(sizeof(char[MAXNUMATTR][MAXATTRNAMELEN]));  // sizeof char[MAXNUMATTR][MAXATTRNAMELEN]
        tmpAttVals[i_prod] = ( char(*)[MAXNUMATTR][MAXATTRLEN]) malloc(sizeof(char[MAXNUMATTR][MAXATTRLEN]));  // sizeof char[MAXNUMATTR][MAXATTRLEN]
        tmpAttTypes[i_prod] = ( int(*)[MAXNUMATTR]) malloc(sizeof(int[MAXNUMATTR]));  // sizeof int[MAXNUMATTR]
        tmpAttSize[i_prod] = ( size_t(*)[MAXNUMATTR]) malloc(sizeof(size_t[MAXNUMATTR]));  // sizeof size_t[MAXNUMATTR]
        strncpy((*tmpAttVals[i_prod])[0], output_averaging_schemes_names[i_prod].c_str(),MAXATTRLEN);
        strncpy((*tmpAttNames[i_prod])[0], "averaging_scheme",MAXATTRNAMELEN);
        (*tmpAttSize[i_prod])[0] = strlen(output_averaging_schemes_names[i_prod].c_str());
        (*tmpAttTypes[i_prod])[0] = NC_CHAR;
    }

    output_binfile->setProductList(ofile.tot_prod, tmpProdNames, tmpAttVals,tmpAttNames,tmpAttSize,tempAttNums,tmpAttTypes);

    //cleanup
    for (int i_prod = 0; i_prod < ofile.tot_prod; i_prod++) {
        free(tmpAttNames[i_prod]);
        free(tmpAttVals[i_prod]);
        free(tmpAttTypes[i_prod]);
        free(tmpAttSize[i_prod]);
    }

    output_binfile->create(ofile.name, input_binfile->nrows);
    // read product L3B attributes
    std::vector<ProductL3Attributes> prod_attrs(input_binfile->nprod());
    for (int j = 0; j < input_binfile->nprod(); j++) {
        input_binfile->get_prodname(j, buf);
        prod_attrs[j] = readProductL3Attributes(input_binfile,buf);
    }
    
    // validate matched products
    productMatching.validate_matched_products(l1rec->l1file,input_binfile);

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

            if (input_binfile->active_data_prod[j]) {
                input_binfile->get_prodname(j, buf);
                if (strcmp(prodlist_in[j], buf) == 0) {
                    input_binfile->readSums(inData[i], ext, j);
                    for (int size = 0; size < 2 * ext; size++) {
                        apply_reverse_averaging_scheme(inData[i][size], prod_attrs[j].averaging_scheme);
                    }
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

            if(input->georegionfile[0]){
                if(!get_georegion(l1rec->lat[ip], l1rec->lon[ip])) {
                    l2rec->l1rec->flags[ip] |= NAVFAIL | ATMFAIL;
                    l2rec->l1rec->navfail[ip] = 1;
                    continue;
                }
            }

            l1rec->nobs[ip] = input_binfile->get_nobs(ip);
            for (int32 iw = 0; iw < n_nlw; iw++) {
                size_t band_index = productMatching.nlw_var.bindex_found[iw];
                size_t bin_index = productMatching.nlw_var.match_l3b_index[iw];
                ipw = ip * ifile.nbands + band_index;
                l2rec->nLw[ipw] = inData[bin_index][2 * ip] / weight;
                l2rec->Rrs[ipw] = MAX(l2rec->nLw[ipw] / ifile.Fonom[band_index], BAD_FLT);
            }
            for (int32 iw = 0; iw < n_rrs; iw++) {
                size_t band_index = productMatching.rrs_var.bindex_found[iw];
                size_t bin_index = productMatching.rrs_var.match_l3b_index[iw];
                ipw = ip * ifile.nbands + band_index;
                l2rec->Rrs[ipw] = inData[bin_index][2 * ip] / weight;
                if (l2rec->nLw[ipw] ==
                    BAD_FLT)  // if nLw not set through an nLw product, calculate it from Rrs and F0
                    l2rec->nLw[ipw] = MAX(l2rec->Rrs[ipw] * ifile.Fonom[band_index], BAD_FLT);
            }
            for (size_t iw = 0; iw < n_rrs_vc; iw++) {
                size_t band_index = productMatching.rrs_vc_var.bindex_found[iw];
                size_t bin_index = productMatching.rrs_vc_var.match_l3b_index[iw];
                ipw = ip * ifile.nbands + band_index;
                l2rec->Rrs[ipw] = inData[bin_index][2 * ip] / weight;
                if (l2rec->nLw[ipw] ==
                    BAD_FLT)  // if nLw not set through an nLw product, calculate it from Rrs and F0
                    l2rec->nLw[ipw] = MAX(l2rec->Rrs[ipw] * ifile.Fonom[band_index], BAD_FLT);
            }
            for (size_t iw = 0; iw < n_rrs_unc; iw++) {
                size_t band_index = productMatching.rrs_unc_var.bindex_found[iw];
                size_t bin_index = productMatching.rrs_unc_var.match_l3b_index[iw];
                ipw = ip * ifile.nbands + band_index;
                l2rec->Rrs_unc[ipw] = inData[bin_index][2 * ip] / weight;
            }
            for (size_t iw = 0; iw < n_rhos; iw++) {
                size_t band_index = productMatching.rhos_var.bindex_found[iw];
                size_t bin_index = productMatching.rhos_var.match_l3b_index[iw];
                ipw = ip * ifile.nbands + band_index;
                l2rec->l1rec->rhos[ipw] = inData[bin_index][2 * ip] / weight;
            }
            l1rec->solz[ip] = 0.0;
            l1rec->sola[ip] = 0.0;
            l1rec->senz[ip] = 0.0;
            l1rec->sena[ip] = 90.0;
            l1rec->delphi[ip] = 90.0;
            if (have_sola) {
                l1rec->sola[ip] = inData[productMatching.index_sola][2 * ip] / weight;
            }
            if (have_sena) {
                l1rec->sena[ip] = inData[productMatching.index_sena][2 * ip] / weight;
            }
            if (have_solz) {
                l1rec->solz[ip] = inData[productMatching.index_solz][2 * ip] / weight;
            }
            if (have_relaz) {
                l1rec->delphi[ip] = inData[productMatching.index_relaz][2 * ip] / weight;
            }
            if (have_senz) {
                l1rec->senz[ip] = inData[productMatching.index_senz][2 * ip] / weight;
            }
            if (!have_solz) {
                // assume noon orbit
                flon = l1rec->lon[ip];
                flat = l1rec->lat[ip];
                gmt = equatorialCrossingTime - lon / 15;
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
        } //ip loop close

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
                tmpData[iprod][ip] = apply_averaging_scheme(tmpData[iprod][ip],output_averaging_schemes.at(iprod));
                outMask[ip] |= (tmpData[iprod][ip] == BAD_FLT);
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
                    outData[iprod][2 * nwrite + 1] = outData[iprod][2 * nwrite] * outData[iprod][2 * nwrite];
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
    snprintf(soft_id,sizeof soft_id,  "%d.%d.%d-%s", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH, GITSHA);
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
    if(input->georegionfile[0]) {
        close_georegion_file();
    }
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
    cout << "Processing Complete at " << ptime << endl << endl;

    return (mainReturnCode);
}
