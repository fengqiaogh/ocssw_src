#ifndef _L1AGEN_OCI_H_
#define _L1AGEN_OCI_H_

#include <stdint.h>
#include <fstream>
#include <timeutils.h>
#include <netcdf>

#include "common.h"

#define PBUFFER_SIZE 32768
#define ANCSIZE 104
#define TLMSIZE 3200
// Chunking the arrays
#define CHUNK_CACHE_SIZE 96 * 1024 * 1024  // 32MiB of cache memory.
#define CHUNK_CACHE_NELEMS 251
#define CHUNK_CACHE_PREEMPTION 1
#define CHUNKBANDS 36
#define CHUNKPIXELS 256
#define CHUNKLINES 128
// Do something like this
//   vector<size_t> chunkVec{CHUNKBANDS,CHUNKPIXELS,CHUNKLINES};
// for dimensions less than these defaults, fix it to the dimension, e.g.:
//      if (dimBands < CHUNKBANDS)
//        chunkVec[0] = dimBands;

const int16_t SWIR_LOFF_DEFAULT[9] = {-16, 80, 64, 0, -16, 64, 0, 80, 96};
const int16_t SWIR_LOFF_ETU[9] = {80, 88, 72, 64, 80, 72, 0, 0, 0};

// From: https://stackoverflow.com/questions/2782725/
//       converting-float-values-from-big-endian-to-little-endian
float ReverseFloat(const float inFloat) {
    float retVal;
    char *floatToConvert = (char *)&inFloat;
    char *returnFloat = (char *)&retVal;

    // swap the bytes into a temporary buffer
    returnFloat[0] = floatToConvert[3];
    returnFloat[1] = floatToConvert[2];
    returnFloat[2] = floatToConvert[1];
    returnFloat[3] = floatToConvert[0];

    return retVal;
}

typedef struct {
    int32_t iyear;
    int32_t iday;
    double sec;
} time_struct;

class l1aFile {
    netCDF::NcFile *l1afile;

    std::string fileName;

    int ngrps;
    int ndims;

    netCDF::NcDim ncDims[1000];
    netCDF::NcGroup ncGrps[10];

   public:
    l1aFile();
    ~l1aFile();

    netCDF::NcFile *ncfile() {
        return l1afile;
    }

    std::string platform;
    int apktsize;
    int bpktsize;
    int EV_APIDs;

    int createl1(char *l1_filename, uint16_t maxsc, uint16_t ncps, uint16_t nbbs, uint16_t nrbs,
                 uint16_t nsps, uint16_t ndcs);

    int parseDims(std::string dimString, std::vector<netCDF::NcDim> &varDims);

    int write_oci_science_data(uint32_t isc, uint16_t nbbs, uint16_t nrbs, uint16_t nswb, uint16_t ncps,
                               uint16_t nsps, uint16_t **bsci, uint16_t **rsci, uint32_t **ssci,
                               int8_t *sfrms);

    int write_processing_control(std::string hktList, std::string l0List, std::string time_start, std::string maxgap, std::string nametag,
                                    std::string swir_loff_set, std::string outlist, std::string outfile, std::string doi, 
                                    std::string pversion, std::string isSPW);

    int write_oci_cal_data(uint32_t isc, uint16_t nbbs, uint16_t nrbs, uint16_t nswb, uint16_t ndcs,
                           uint16_t ndss, uint16_t *dark_b, uint16_t *dark_r, uint32_t *dark_s,
                           int8_t *sdfrms);

    int write_oci_scan_metadata(uint32_t isc, uint8_t *ancdata, uint8_t *seqerr, int8_t *linerr,
                                int32_t *spinID, time_struct &starttime);

    int write_oci_global_metadata(time_struct &starttime, time_struct &endtime, std::string l1a_name,
                                  std::string sdir, std::string edir, uint32_t isc, short dtype,
                                  uint16_t smode, uint16_t cdsmode, std::ofstream &fout);

    int write_oci_ancil_data(uint32_t isc, uint8_t *ancdata);

    int write_oci_tlm_data(itab *itable, uint32_t ntlm, uint8_t (*tlmdata)[TLMSIZE], int32_t *spinID,
                           uint16_t &cdsmode, uint32_t isc, const time_struct &starttime);
    int write_navigation(std::string hktlist, time_struct &starttime, time_struct &endtime);
    int close();
};

int make_oci_line_index(itab *itable, int16_t *cindex, int16_t *sindex, int16_t *cdindex, int16_t *sdindex,
                        int16_t *swir_loff);

int unpack_oci_sci(uint32_t npkts, int32_t spin, uint16_t ncps, uint16_t nsps, uint16_t msps,
                   uint16_t &nbands, uint16_t btaps[16], uint16_t rtaps[16], uint8_t (*pbuffer)[PKTSIZE],
                   uint16_t **bbands, uint16_t **rbands, uint32_t **sbands, int16_t *blines, int16_t *rlines,
                   int16_t *slines, uint16_t &btype, uint16_t bagg[16], uint16_t &rtype, uint16_t ragg[16],
                   int8_t *sfrm, int &iret);

int unpack_ccd_packet(uint8_t *packet, uint16_t btaps[16], uint16_t rtaps[16], uint16_t &ccdid,
                      uint32_t &line, uint16_t &dtype, uint16_t &iagg, uint16_t jagg[16], uint16_t &nbands,
                      uint16_t **ccddata, uint16_t ossdata[16]);

int check_load_oci_data(short dtype, uint16_t ncps, uint16_t nsps, uint16_t ndcs, uint16_t ndss,
                        uint16_t nbbs, uint16_t nrbs, uint16_t nswb, int16_t *cindex, int16_t *sindex,
                        int16_t *cdindex, int16_t *sdindex, uint16_t **bbands, uint16_t **rbands,
                        uint32_t **sbands, int16_t *blines, int16_t *rlines, int16_t *slines, uint16_t **bsci,
                        uint16_t **rsci, uint32_t **ssci, uint16_t **bdark, uint16_t **rdark,
                        uint32_t **sdark, int8_t &linerr, int &icheck);

int unpack_swir_packet(uint8_t *packet, int16_t *slines, uint8_t *swirfrm, uint32_t *swirdata);

uint8_t check_sum(int32_t nc, uint8_t *dat, uint8_t *chk);

int eight20(uint8_t *inbytes, uint32_t *outsamples);

int createNCDF(netCDF::NcGroup &ncGrp, const char *sname, const char *lname, const char *standard_name,
               const char *units, void *fill_value, const char *flag_values, const char *flag_meanings,
               const char *reference, double low, double high, int nt, std::vector<netCDF::NcDim> &varVec);

int find_nav_index(char *hkt_t_start, double usec_l1a_start, double usec_l1a_end, double *nav_t, size_t n_nav,
                   int &ind_start, int &ind_end);

inline int expandEnvVar(std::string *sValue) {
    if ((*sValue).find_first_of("$") == std::string::npos)
        return 0;
    std::string::size_type posEndIdx = (*sValue).find_first_of("/");
    if (posEndIdx == std::string::npos)
        return 0;
    const std::string envVar = sValue->substr(1, posEndIdx - 1);
    char *envVar_str = getenv(envVar.c_str());
    if (envVar_str == 0x0) {
        printf("Environment variable: %s not defined.\n", sValue->c_str());
        exit(1);
    }
    *sValue = envVar_str + (*sValue).substr(posEndIdx);

    return 0;
}

l1aFile::l1aFile() {
}

l1aFile::~l1aFile() {
    delete(l1afile);
}

#endif  // _L1AGEN_OCI_H_
