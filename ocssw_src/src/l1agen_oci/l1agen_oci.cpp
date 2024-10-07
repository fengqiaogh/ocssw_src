#include <stdio.h>
// #include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <regex>

#include <sstream>
#include <fstream>
#include <iomanip>
#include <getopt.h>
#include <libgen.h>
#include <boost/algorithm/string.hpp>

#include "nc4utils.h"
#include "global_attrs.h"
#include "l1agen_oci.h"
#include "genutils.h"
#include "l0stream.hpp"

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

// FIXME
// we need to use the timeutils function ccsds_to_yds() instead
int ccsds_sec_to_yds(uint8_t *cctime, int32_t *iyear, int32_t *iday, double *sec) {
    uint32_t ui32;
    uint32_t ccsec;

    memcpy(&ui32, cctime, 4);
    ccsec = SWAP_4(ui32);

    double dsec = (double)ccsec;
    int leap = leapseconds_since_1993(dsec);
    ccsec -= (leap + 27);

    *iday = ccsec / 86400;
    int32_t jday = *iday + 2436205;  // Jan. 1, 1958 is Julian day 2436205
    jdate(jday, iyear, iday);

    // Get milliseconds
    int32_t msec = cctime[4] * 256 + cctime[5];
    int32_t isec = ccsec % 86400;
    *sec = isec + msec / 65536.0;

    return 0;
}

#define VERSION "1.18.00_2024-06-12"

//    Modification history:
//  Programmer     Organization   Date     Ver     Description of change
//  ----------     ------------   ----     ---     ---------------------
//  Joel Gales     FutureTech     09/20/18 0.10    Original development
//                                                 based on IDL routines
//                                                 developed by F. Patt
//  Joel Gales     SAIC           10/26/18 0.20    Complete initial alpha version
//  Joel Gales     SAIC           12/20/18 0.30    Complete initial version
//                                                 of SWIR bands
//  Joel Gales     SAIC           04/26/19 0.40    Implement code changes from
//                                                 F. Patt
//  Joel Gales     SAIC           05/20/19 0.50    Implement code changes from
//                                                 F. Patt
//  Joel Gales     SAIC           06/04/19 0.60    Implement code changes from
//                                                 F. Patt
//  Joel Gales     SAIC           07/01/19 0.70    Add support for outlist
//  Joel Gales     SAIC           07/23/19 0.75    Implement code changes from
//                                                 F. Patt
//  Joel Gales     SAIC           08/02/19 0.80    Flag and ignore packets
//                                                 greater than 1200 bytes.
//                                                 Remove common routines and
//                                                 Place in common.cpp
//  Joel Gales     SAIC           08/09/19 0.81    Fix memory overwrite bug in
//                                                 unpack_ccd_packet() in the
//                                                 ossdata array.  Initialize
//                                                 "lines" and "bands" arrays.
//                                                 Add support for granules with
//                                                 missing blue/red bands.
//  Joel Gales     SAIC           10/25/19 0.82    Exit if EOF before finding
//                                                 good packet
//  Joel Gales     SAIC           11/20/19 0.85    Implement code changes from
//                                                 F. Patt
//  Joel Gales     SAIC           11/25/19 0.86    Change 60 to maxsc fpr dspn
//                                                 comparison
//  Joel Gales     SAIC           11/27/19 0.87    Change L1A name to standard
//                                                 Trap granules with no ancillary
//                                                 granules
//  Joel Gales     SAIC           12/05/19 0.88    Add nametage command line
//                                                 option
//  Joel Gales     SAIC           01/03/20 0.90    Update and add additional
//                                                 telemetry fields
//  Joel Gales     SAIC           01/22/20 0.91    Add scomp/maxgap
//  Joel Gales     SAIC           01/23/20 0.92    Check for EOF when checking
//                                                 for zero science pixels
//                                                 Zero out sci/dark fields
//  Joel Gales     SAIC           01/28/20 0.93    Bug fixes and updates
//  Joel Gales     SAIC           02/07/20 0.94    Implement changes from F.Patt
//                                                 (020720)
//  Joel Gales     SAIC           02/20/20 0.95    Fixed bugs with bagerr, ragerr,
//                                                 digerr and dautemp
//  Joel Gales     SAIC           02/25/20 0.96    Fixed bug in ccsds..
//                                                 Writing 6 bytes into uint32_t
//  Joel Gales     SAIC           03/04/20 0.97    Implement changes from F.Patt
//                                                 (022820)
//  Joel Gales     SAIC           03/24/20 0.98    Fix HAM_side issue
//  Joel Gales     SAIC           04/02/20 0.99    Implement SWIR mode changes
//  Liang Hong     SAIC           04/28/20 0.9901  Fixed last scan_line_attributes
//                                                 records in last output granule
//  Liang Hong     SAIC           05/20/20 0.9902  Handle dark zone science data
//  Liang Hong     SAIC           08/11/20 0.9903  Handle 0 pix 0 band sci/cal data in l1a
//                                                 "no data" and order variation in input
//  Liang Hong     SAIC           08/24/20 0.9904  fixed a bug to close file with 0-scan
//                                                 ; demand non-negative indices
//  Liang Hong     SAIC           09/02/20 0.9905  HKT and science data overlap check; order
//                                                 SWIR bands in ascending wavelength
//  Liang Hong     SAIC           10/28/20 0.9906  handle fill values in SWIR
//  Liang Hong     SAIC           10/29/20 0.9907  APID for the MCE HK packet changed to 713
//  Liang Hong     SAIC           11/23/20 0.9908  fixed rare start/end time error; SWIR band
//                                                 -specific pixel shifts as input option;
//                                                 run with optional granule start time;
//                                                 fixed science packet sequence error flag
//                                                 fixed HAM_side value after index nmce
//  Liang Hong     SAIC           12/01/20 0.99.00 fixed number_of_filled_scans in metadata
//  Liang Hong     SAIC           04/22/21 0.99.10 fixed no ancillary data exit conditions
//  Liang Hong     SAIC           06/17/21 0.99.20 generate 1 telemetry entry when no data
//                                                 bug fix in duplicated reading of last tlm
//  Liang Hong     SAIC           07/27/21 0.99.21 return 110 for No ancillary packets
//  Liang Hong     SAIC           01/07/22 1.00.00 temperature fields update in HKT packets
//                                                 SWAP_4 in common.h updated
//  Liang Hong     SAIC           01/11/22 1.00.01 OCI SWIR Raw mode reading correction
//  Liang Hong     SAIC           01/25/22 1.00.11 blue and red spectral mode order correction
//  Liang Hong     SAIC           03/11/22 1.01.00 added telemetry for the solar calibrator;
//                                                 references in metadata; write ancillary_tlm
//  Liang Hong     SAIC           03/30/22 1.02.00 SPCA and lunar stare data types added
//  Liang Hong     SAIC           04/14/22 1.03.00 update error checking and handling
//  Liang Hong     SAIC           04/21/22 1.03.01 updated packet # threshold; mode table check
//  Liang Hong     SAIC           04/24/22 1.03.02 exit read packets when endfile
//  Liang Hong     SAIC           05/11/22 1.04.00 update of APID list and check all for spin #
//                                                 added navigation data from hkt input
//  Liang Hong     SAIC           05/27/22 1.04.01 limit packet reading to maxpkts
//  Liang Hong     SAIC           06/24/22 1.05.00 added CCD masks; fixed bugs in nav data
//  Liang Hong     SAIC           11/18/22 1.06.00 noSPW option set for OCI test data only
//  Liang Hong     SAIC           11/22/22 1.07.00 granule starts with specified start time
//  Liang Hong     SAIC           11/28/22 1.08.01 updated start time, outfile and --noSPW options
//  Liang Hong     SAIC           12/05/22 1.08.03 clear outlist if no L1A generated, return 111
//  Liang Hong     SAIC           01/10/23 1.09.04 added data type to last column of outlist
//  Gwyn Fireman   SAIC           02/02/23 1.09.05 CF-compliant output; new CDL file; history attribute
//  Liang Hong     SAIC           02/28/23 1.10.05 update OCI solar cal datatypes
//  Liang Hong 	   SAIC           05/01/23 1.11.00 navigation data read and fill update
//  Liang Hong     SAIC           05/15/23 1.11.01 bug fix in PACE navigation data read
//  Liang Hong     SAIC           05/30/23 1.12.01 added usage of tilt flag
//  Liang Hong     SAIC           06/23/23 1.13.00 Fill value update; fixed cross date data search issue
//  Liang Hong     SAIC           06/26/23 1.14.00 Added maxgap as an option allowed missing scans/file
//  Gwyn Fireman   SAIC           07/10/23 1.14.01  Read global metadata from json file
//  Liang Hong     SAIC           09/12/23 1.15.00 Added SCA diffuser; 1 granule given start, granule_len
//  Wei Jiang      SAIC           03/15/24 1.16.00 maxtime (mtime) does not update if there's a data type
//  change Wei Jiang      SAIC           03/18/24 1.16.00 handles single .oci file that have 64 byte header
//  intact

void print_usage() {
    cout << endl
         << "l1agen_oci OCI_packet_file granule_len\n"
         << "    [-g | --maxgap maximum missing scans allowed in an output file]\n"
         <<  // allowed number of missing scans in one output file, 0 means no limit
        "    [-k | --hktlist SC_HKT_input_list]\n"
         <<  // file containing list of Spacecraft housekeeping telemetry file names
        "    [-s | --swir_loff_set SWIR_LOFF_config, list of 9 comma separated integers]\n"
         << "    [-t | --start_time YYYYmmddTHHMMSS or YYYY-mm-ddTHH:MM:SS]\n"
         << "    [-o | --outlist output_list_file]\n"
         << "    [-f | --outfile output_file_name]\n"
         << "    [-p | --nametag first part of filename (default=PACE_OCI)]\n"
         << "    [-d | --doi doi_string]\n"
         << "    [-v | --pversion processing_version]\n"
         << "    [-n | --noSPW no space wire header]\n\n"
         << "Return   1 Fatal error\n"
         << "       110 No ancillary packets found in file\n"
         << "       120 No L1A file was generated\n";
}

ofstream tempOut;

void verify_packet_len(L0Stream *tfileStream, uint32_t &len, uint32_t apid, int32_t endfile, bool isSPW) {
    read_packet(tfileStream, NULL, len, apid, endfile, isSPW);
    if (len > PKTSIZE) {
        cout << "Packet too big (" << len << ") for buffer (" << PKTSIZE << ")" << endl;
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    cout << "l1agen_oci " << VERSION << " (" << __DATE__ << " " << __TIME__ << ")" << endl;

    if (argc == 1) {  // TODO: Usage function
        print_usage();
        return 0;
    }
    string history = call_sequence(argc, argv);

    int maxgap = 10;    // by default, the max allowed missing scan  = 10
    bool isSPW = true;  // by default, OCI data is from DSB and has space wire header
    int return_status = 0;
    int c;
    string hktlist = "";
    string swir_loff_set = "";
    string time_start = "";
    string outlist = "";
    string outfname = "";
    string nametag = "";
    string doi = "";
    string pversion = "Unspecified";

    while (1) {  // TODO: Make loop not infinite
        static struct option long_options[] = {{"maxgap", required_argument, 0, 'g'},
                                               {"hktlist", required_argument, 0, 'k'},
                                               {"swir_loff_set", optional_argument, 0, 's'},
                                               {"start_time", required_argument, 0, 't'},
                                               {"outlist", required_argument, 0, 'o'},
                                               {"outfile", required_argument, 0, 'f'},
                                               {"nametag", required_argument, 0, 'p'},
                                               {"doi", required_argument, 0, 'd'},
                                               {"pversion", required_argument, 0, 'v'},
                                               {"noSPW", no_argument, 0, 'n'},
                                               {0, 0, 0, 0}};

        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "g:k:s:t:o:f:p:d:v:n", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
            case 'g':
                maxgap = atoi(optarg);
                if (maxgap == 0)
                    maxgap = 65535;
                break;

            case 'k':
                hktlist.assign(optarg);
                break;

            case 's':
                swir_loff_set.assign(optarg);
                break;

            case 't':
                time_start.assign(optarg);
                break;

            case 'o':
                //        printf ("option -o with value `%s'\n", optarg);
                outlist.assign(optarg);
                break;

            case 'f':
                outfname.assign(optarg);
                break;

            case 'p':
                nametag.assign(optarg);
                break;

            case 'd':
                doi.assign(optarg);
                break;

            case 'v':
                pversion.assign(optarg);
                break;

            case 'n':
                isSPW = false;  // noSPW indicates data is from OCI instrument testing
                break;

            default:
                abort();
        }
    }

    string l0Path = argv[optind];
    string granuleLengthStr = argv[optind + 1];
    string time_start_copy = time_start;

    nc_set_chunk_cache(CHUNK_CACHE_SIZE, CHUNK_CACHE_NELEMS, CHUNK_CACHE_PREEMPTION);

    // Error if output is not just one granule and an output file name is specified
    if ((strcmp(granuleLengthStr.c_str(), "0") != 0) && (time_start.empty()) && (outfname.compare("") != 0)) {
        cout << "Set granule length parameter to 0 or provide start time when specifying output file name."
             << endl;
        exit(1);
    }

    struct tm tm_start;
    if (!time_start.empty()) {
        std::regex time_pattern0("([0-9]{4})-([0-9]{2})-([0-9]{2})T([0-9]{2}):([0-9]{2}):([0-9]{2})");
        std::regex time_pattern1("([0-9]{4})([0-9]{2})([0-9]{2})T([0-9]{2})([0-9]{2})([0-9]{2})");

        if (regex_match(time_start, time_pattern0)) {
            strptime(time_start.c_str(), "%Y-%m-%dT%H:%M:%S", &tm_start);
        } else if (regex_match(time_start, time_pattern1)) {
            strptime(time_start.c_str(), "%Y%m%dT%H%M%S", &tm_start);
        } else {
            cout << "Start time format is wrong. Use  YYYY-mm-ddTHH:MM:SS  or YYYYmmddTHHMMSS." << endl;
            exit(1);
        }
    }

    ofstream fout;
    if (outlist.compare("") != 0)
        fout.open(outlist.c_str());

    vector<string> fileNames;
    // Read S/C telemetry data
    // read_pace_telemetry(hktlist,iyrs,idays,otime,pos,vel,atime,quat,arate,ttime,tilt);

    fileNames = readFileList(l0Path);
    L0Stream tfileStream(fileNames);

    // if the input is empty
    if (fileNames.size() == 0) {
        cout << "Error - No L0 files found" << endl;
        return 0;
    }
    if (tfileStream.fail())
        cout << "Failed to open tfileStream at line " << __LINE__ << endl;

    uint32_t apid = 0;
    uint32_t len = 0;
    int32_t endfile = 0;
    uint8_t fpacket[PKTSIZE];
    uint8_t apacket[ANCSIZE];
    uint8_t apacket0[ANCSIZE];
    uint32_t maxpkts = 30000;  // LH, 4/21/2022, v1.03.01
    int acomp;
    // int maxgap = 10;

    uint8_t **pbuffer0 = new uint8_t *[maxpkts];
    pbuffer0[0] = new uint8_t[PKTSIZE * maxpkts];
    for (size_t i = 1; i < maxpkts; i++)
        pbuffer0[i] = pbuffer0[i - 1] + PKTSIZE;

    uint8_t **pbuffer1 = new uint8_t *[maxpkts];
    pbuffer1[0] = new uint8_t[PKTSIZE * maxpkts];
    for (size_t i = 1; i < maxpkts; i++)
        pbuffer1[i] = pbuffer1[i - 1] + PKTSIZE;

    uint32_t npkts0;
    int32_t ancind0;
    //  uint32_t spn0, spn;
    int32_t spn0, spn, spnp;
    vector<int32_t> tlmind0;

    // Get first science or ancillary packet
    verify_packet_len(&tfileStream, len, apid, endfile, isSPW);
    read_packet(&tfileStream, fpacket, len, apid, endfile, isSPW);

    // Grab science or ancillary packets
    apid = (fpacket[0] % 8) * 256 + fpacket[1];
    int nsk = 0;
    while (apid != 636 && apid != 700 && apid != 720 && !endfile) {
        nsk++;
        verify_packet_len(&tfileStream, len, apid, endfile, isSPW);
        read_packet(&tfileStream, fpacket, len, apid, endfile, isSPW);
        apid = (fpacket[0] % 8) * 256 + fpacket[1];
    }
    if (endfile) {
        cout << "No science packets found in file" << endl;
        exit(1);
    }
    if (nsk > 0)
        cout << nsk << " packets skipped" << endl;

    // Read first scan and check for ancillary packet
    ancind0 = -1;
    tlmind0.clear();

    int32_t iyear, iday;
    double stime;

    // Get granule period in minutes
    int32_t mper;
    istringstream(granuleLengthStr) >> mper;

    itab itable[10];
    string dtypes[] = {"",      "",      "_DARK", "_SOL", "_SPCA",   "_LIN",    "_LUN",
                       "_DIAG", "_STAT", "_SPEC", "",     "_SNAP-X", "_SNAP-I", "_NBSB"};
    string smodes[] = {"", "_SDIAG", "_SRAW", "_STEST"};

    uint8_t **pbuffer = new uint8_t *[maxpkts];
    pbuffer[0] = new uint8_t[PKTSIZE * maxpkts];
    for (size_t i = 1; i < maxpkts; i++)
        pbuffer[i] = pbuffer[i - 1] + PKTSIZE;

    uint16_t ncps, nbbs, nrbs, nsps, ndcs, ndss, btaps[16], rtaps[16], msps;

    int8_t *linerr = new int8_t[maxpkts];
    uint8_t *seqerr = new uint8_t[maxpkts];
    for (size_t i = 0; i < maxpkts; i++) {
        linerr[i] = 0;
        seqerr[i] = 0;
    }
    uint8_t noseq = 255;
    uint16_t smode = 0;

    while ((ancind0 == -1) && !endfile) {
        read_oci_scan_packets(&tfileStream, fpacket, (uint8_t(*)[PKTSIZE]) & pbuffer0[0][0], npkts0, spn0,
                              ancind0, tlmind0, noseq, endfile, isSPW);
    }

    if (endfile) {
        cout << "No ancillary packets found in file" << endl;
        exit(110);  // ver 0.99.21
    }


    // time variables -- change when refactoring!
    int32_t ltime = -1;     // time used to name the file
    int32_t mtime = -1;     // max time
    uint16_t maxsc;         // max scans
    double scanp = 0.1755;  // changed from 1.0 / 5.737 to 0.1755 in v1.04.00


    ////////////////////// Main Loop ///////////////////
    while (!endfile) {
        // ncps  - number of CCD band pixels (ccd_pixels)
        // nsps  - number of SWIR band pixels
        // ndcs  - number of dark collect pixels
        // ndss  - number of dark SWIR pixels
        // nbbs  - number of blue bands
        // nrbs  - number of red bands
        // btaps - Blue CCD tap enable flags (1 = enabled)
        // rtaps - Red  CCD tap enable flags (1 = enabled)

        // Check for zero science pixels
        if (ancind0 != -1) {
            memcpy(apacket0, &pbuffer0[ancind0][0], ANCSIZE);
            get_band_dims(apacket0, ncps, nbbs, nrbs, nsps, ndcs, ndss, btaps, rtaps, itable);
        }

        while ((ncps == 1 || ancind0 == -1) && !endfile) {
            if (ancind0 != -1)
                cout << "Ancillary packet has zero science pixels at spin " << spn0 << endl;
            read_oci_scan_packets(&tfileStream, fpacket, (uint8_t(*)[PKTSIZE]) & pbuffer0[0][0], npkts0, spn0,
                                  ancind0, tlmind0, noseq, endfile, isSPW);
            if (ancind0 != -1) {
                memcpy(apacket0, &pbuffer0[ancind0][0], ANCSIZE);
                get_band_dims(apacket0, ncps, nbbs, nrbs, nsps, ndcs, ndss, btaps, rtaps, itable);
            }
        }

        if (ancind0 == -1) {
            cout << "No ancillary packets for last granule" << endl;
            break;
        }

        cout << endl;
        cout << "ncps - number of CCD band pixels (ccd_pixels): " << ncps << endl;
        cout << "nsps - number of SWIR band pixels:             " << nsps << endl;
        cout << "ndcs - number of dark collect pixels:          " << ndcs << endl;
        cout << "ndss - number of dark SWIR pixels:             " << ndss << endl;
        cout << "nbbs - number of blue bands:                   " << nbbs << endl;
        cout << "nrbs - number of red bands:                    " << nrbs << endl;

        // get the first scan time of the current packet and save it to stime
        get_anc_packet_time(apacket0, iyear, iday, stime);
        double stimp = stime - scanp;

        time_struct starttime, endtime;
        starttime.iyear = iyear;
        starttime.iday = iday;
        starttime.sec = stime;

        /*
            --start_time was given, set the file time name (ltime) to be this time.
            time_start string will be cleared after so that if there's a non-science
            data inbetween the science data, the non-science data won't use --start_time
            and the next science data file will also not use it
            
        */ 
        if (!time_start.empty()) {  

            ltime = (int32_t)(tm_start.tm_hour * 3600 + tm_start.tm_min * 60 + tm_start.tm_sec);
            // adjust for different dates
            ltime += (jday(tm_start.tm_year + 1900, 1, tm_start.tm_yday + 1) - jday(iyear, 1, iday)) * 86400;

            // set max time for the L0 files and max scans depending on granuel length.
            // mper == granuel length in minutes
            // if minutes is 0, then max time is +600 by default to read the entire file.
            // else, add (granuel length * 60) to start_time to get max time for this file. 
            mtime = mper > 0 ? ltime + mper * 60 : ltime + 600;

            // update max scans based on granduel length
            maxsc = mper > 0 
                        ? (uint16_t)((mper * 60 / scanp) + 2) 
                        : 3600;
            
            // **clear time_start so this does not run again **
            time_start.clear();

            // moves packet pointer forward if the current scantime is earlier than when we want to start
            // specified by ltime 
            while (((stime < ltime) || (ncps == 1)) && !endfile) {
                read_oci_scan_packets(&tfileStream, fpacket, (uint8_t(*)[PKTSIZE]) & pbuffer0[0][0], npkts0,
                                      spn0, ancind0, tlmind0, noseq, endfile, isSPW);
                if (ancind0 != -1) {
                    memcpy(apacket0, &pbuffer0[ancind0][0], ANCSIZE);
                    get_anc_packet_time(apacket0, iyear, iday, stime);

                    get_band_dims(apacket0, ncps, nbbs, nrbs, nsps, ndcs, ndss, btaps, rtaps, itable);
                    if (ncps == 1)
                        cout << "Ancillary packet has zero science pixels at spin " << spn0 << endl;
                }
            }

            if (endfile) {
                cout << "No science data in time range" << endl;
                return 110;
            }

            starttime.sec = stime;
        } 

        /*  if --start_time was not set, then use the first scantime of the packet.
        
            This covers instances where the science data gets split because non-science data
            is in-between. The first science L1A will use the --start_time and the second will use the
            first scantime of the packet when science data resumes. 
        */
        else {
            // when naming, use stime
            ltime = (int32_t)stime; 

            // update mtime only if it has not been set before. if set, then dont touch it
            if (mtime == -1) {
                mtime = mper > 0 ? ltime + mper * 60 : ltime + 600;
            }
            // update max scans based on granduel length
            maxsc = mper > 0 
                        ? (uint16_t)((mper * 60 / scanp) + 2) 
                        : 3600;
        }

        acomp = 1;

        // number of SWIR bands
        nsps = ((nsps + 7) / 8) *
               8;  // Round up SWIR number of pixels to a multiple of 8,  LH, 4/15/2022, v1.03.00
        unsigned short nswb = 9;

        int16_t cindex[32768];
        int16_t sindex[32768 * nswb];
        int16_t cdindex[32768];
        int16_t sdindex[32768];
        int16_t swir_loff[9];

        for (size_t i = 0; i < 32768; i++) {
            cindex[i] = -1;
            for (size_t j = 0; j < nswb; j++)
                sindex[i * nswb + j] = -1;
            cdindex[i] = -1;
            sdindex[i] = -1;
        }

        if (!swir_loff_set.empty()) {
            if (swir_loff_set.compare("ETU") == 0) {
                memcpy(swir_loff, SWIR_LOFF_ETU, 9 * sizeof(int16_t));
            } else {
                vector<string> parts;
                boost::split(parts, swir_loff_set, boost::is_any_of(","));
                if (parts.size() != 9) {
                    cout << "Processing with single file option" << endl;
                    exit(EXIT_FAILURE);
                }
                for (size_t i = 0; i < 9; i++) {
                    swir_loff[i] = atoi(parts[i].c_str());
                }
            }
        } else {
            memcpy(swir_loff, SWIR_LOFF_DEFAULT, 9 * sizeof(int16_t));
        }

        // int32_t ldark = 32768;

        // if (stime < 15299.458) {
        //  cout << "make_oci_line_index" << endl;
        //  make_oci_line_index( itable, cindex, sindex, ldark);
        // }
        //  make_oci_line_index( itable, cindex, sindex, ldark);
        make_oci_line_index(itable, cindex, sindex, cdindex, sdindex, swir_loff);

        // Get SWIR band data mode
        // uint16_t smode = 0;      // LH, 09/10/2020
        get_swir_mode((uint8_t(*)[PKTSIZE]) & pbuffer0[0][0], npkts0, smode);
        uint16_t smodep = smode;
        int scomp = 1;

        uint16_t cdsmode;

        // Determine start and end time of granule
        uint32_t jd0 = jday(iyear, 1, iday);
        int16_t yr16 = (int16_t)iyear;
        int16_t doy = (int16_t)iday;
        int16_t month, dom;
        yd2md(yr16, doy, &month, &dom);

        int32_t ih = 0;
        int32_t mn = 0;
        int32_t isec = 0;

        // ---- SETTING TIME NAMING ----
        ih = (int32_t)(ltime / 3600);
        mn = (int32_t)((ltime - ih * 3600) / 60);
        isec = (int32_t)(ltime - ih * 3600 - mn * 60);
       

        stringstream timestr, datestr;
        timestr << setfill('0') << setw(2) << ih << setfill('0') << setw(2) << mn << setfill('0') << setw(2)
                << isec;

        datestr << setfill('0') << setw(4) << iyear << setfill('0') << setw(2) << month << setfill('0')
                << setw(2) << dom;

        static l1aFile outfile;
        string basenme;
        char* tmpFileStr = strdup(l0Path.c_str());  // need this since basename may modify the char* passed in
        basenme.assign(basename(tmpFileStr));
        free(tmpFileStr);

        if (nametag == "")
            nametag.assign("PACE_OCI");
        short dtype = itable[1].dtype;  // ;  Get data type for file name and metadata
        short maxdtype = 2;
        for (size_t i = 0; i < 10; i++) {
            if (itable[i].dtype > maxdtype)
                maxdtype = itable[i].dtype;
        }
        // if (maxdtype > 2) dtype = maxdtype;
        if ((maxdtype != 2) && (maxdtype != 10)) {
            dtype = maxdtype;  // LH, 8/24/2020
            dtype = maxdtype;
            cout << "\nWARNING: Non-Science Data is now being processed. Type: " << dtypes[dtype] << "\n"
                 << endl;
        }

        // data type mod for ETU before June 2020
        if ((jd0 < 2459000) && dtype == 11)
            dtype = 9;

        string l1a_name = outfname;
        if (outfname.compare("") == 0) {
            l1a_name = nametag + dtypes[dtype] + smodes[smode];
            // keep input filename substrings for OCI instrument test data
            if (!isSPW)
                l1a_name += "_" + basenme.substr(0, 4) + basenme.substr(5, 3) + basenme.substr(9, 3);
            l1a_name += "." + datestr.str() + "T" + timestr.str() + ".L1A.nc";
        }

        // Initialize data arrays

        // blue band dark collect data for granule
        uint16_t **dark_b = new uint16_t *[maxsc];
        dark_b[0] = new uint16_t[ndcs * (nbbs > 0 ? nbbs : 1) * maxsc];
        for (size_t i = 1; i < maxsc; i++)
            dark_b[i] = dark_b[i - 1] + ndcs * (nbbs > 0 ? nbbs : 1);

        // red band dark collect data for granule
        uint16_t **dark_r = new uint16_t *[maxsc];
        dark_r[0] = new uint16_t[ndcs * (nrbs > 0 ? nrbs : 1) * maxsc];
        for (size_t i = 1; i < maxsc; i++)
            dark_r[i] = dark_r[i - 1] + ndcs * (nrbs > 0 ? nrbs : 1);

        // swir band dark collect data for granule
        uint32_t **dark_s = new uint32_t *[maxsc];
        dark_s[0] = new uint32_t[ndss * (nswb > 0 ? nswb : 1) * maxsc];
        for (size_t i = 1; i < maxsc; i++)
            dark_s[i] = dark_s[i - 1] + ndss * (nswb > 0 ? nswb : 1);

        uint16_t **bdark = new uint16_t *[nbbs];
        uint16_t **rdark = new uint16_t *[nrbs];
        uint32_t **sdark = new uint32_t *[nswb];

        uint16_t **bsci = new uint16_t *[nbbs];
        bsci[0] = new uint16_t[ncps * (nbbs > 0 ? nbbs : 1)];
        for (size_t i = 1; i < nbbs; i++)
            bsci[i] = bsci[i - 1] + ncps;

        uint16_t **rsci = new uint16_t *[nrbs];
        rsci[0] = new uint16_t[ncps * (nrbs > 0 ? nrbs : 1)];
        for (size_t i = 1; i < nrbs; i++)
            rsci[i] = rsci[i - 1] + ncps;

        uint32_t **ssci = new uint32_t *[nswb];
        ssci[0] = new uint32_t[nsps * (nswb > 0 ? nswb : 1)];
        for (size_t i = 1; i < nswb; i++)
            ssci[i] = ssci[i - 1] + nsps;

        uint8_t **ancdata = new uint8_t *[maxsc + 1];
        ancdata[0] = new uint8_t[ANCSIZE * (maxsc + 1)];
        for (size_t i = 1; i < (size_t)(maxsc + 1); i++)
            ancdata[i] = ancdata[i - 1] + ANCSIZE;

        uint32_t maxtlm = maxsc * 10;
        uint8_t **tlmdata = new uint8_t *[maxtlm];
        tlmdata[0] = new uint8_t[TLMSIZE * (maxtlm)];
        for (size_t i = 1; i < (size_t)(maxtlm); i++)
            tlmdata[i] = tlmdata[i - 1] + TLMSIZE;

        // uint16_t ncpt = ncps + ndcs;
        // uint16_t nspt0 = nsps + ndss;
        uint16_t ncpt = ncps;   // Dark pixels now included in science pixels
        uint16_t nspt0 = nsps;  // Dark pixels now included in science pixels

        // Round up SWIR number of pixels to a multiple of 8
        uint16_t nspt = ((nspt0 + 7) / 8) * 8;

        // Note: "bands" arrays here are the reverse order as IDL versions
        uint16_t **bbands = new uint16_t *[ncpt];
        uint16_t **rbands = new uint16_t *[ncpt];
        for (size_t i = 0; i < ncpt; i++) {
            bbands[i] = NULL;
            rbands[i] = NULL;
        }

        int16_t *blines = new int16_t[ncpt];
        int16_t *rlines = new int16_t[ncpt];

        for (size_t i = 0; i < ncpt; i++)
            blines[i] = -1;
        for (size_t i = 0; i < ncpt; i++)
            rlines[i] = -1;

        msps = 4096;  // May be full-scan SWIR pixels in ETU data
        uint32_t **sbands = new uint32_t *[msps];
        for (size_t i = 0; i < msps; i++) {
            sbands[i] = new uint32_t[nswb];
        }

        int16_t *slines = new int16_t[msps];
        //    for (size_t i=0; i<msps; i++) slines[i] = 0;

        // swir frame type for dark collect
        int8_t **sdfrms = new int8_t *[maxsc];
        sdfrms[0] = new int8_t[ndss * maxsc];
        for (size_t i = 1; i < maxsc; i++)
            sdfrms[i] = sdfrms[i - 1] + ndss;
        memset(sdfrms[0], -1, sizeof(int8_t) * ndss * maxsc);

        int8_t *sfrm = new int8_t[msps];
        int8_t *sfrms = new int8_t[msps];

        // Initialize output file and get object IDs for EV data
        cout << "Creating: " << l1a_name.c_str() << endl;
        cout << "Starting at spin number " << spn0 << endl;
        outfile.createl1((char *)l1a_name.c_str(), maxsc, ncps, nbbs, nrbs, nsps, ndcs);

        std::string maxgap_string = std::to_string(maxgap);
        std::string noSPW_string = std::to_string(isSPW);
        // write process control group
        outfile.write_processing_control(hktlist, l0Path, time_start_copy, maxgap_string, nametag, swir_loff_set, outlist, l1a_name, doi, pversion, noSPW_string);
        // Write output filename to outlist if needed.
        long outlistpos;
        if (outlist.compare("") != 0) {
            outlistpos = fout.tellp();
            fout << l1a_name.c_str() << " ";
        }

        // Read and process OCI scans
        int iret = 0;
        int icheck = 0;
        uint32_t isc = 0;
        uint32_t npkts, npkt1;
        int32_t ancind;
        vector<int32_t> tlmind;
        spnp = spn0;
        int32_t enddata = 0;
        int dspn = 1;
        uint32_t ntlm = 0;

        uint16_t btype, rtype;
        uint16_t bagg[16], ragg[16];
        uint16_t nbands;

        //    tempOut.open ("bsci.bin", ios::out | ios::trunc | ios::binary);

        // ---------------------------------------------------------------------
        // ---------------------------------------------------------------------
        // ---------------------------------------------------------------------

        while (stime < mtime && acomp && !enddata && scomp && (dspn <= maxgap) && isc < maxsc) {
            // sometimes, data contain more scans than possible, set isc<maxsc to avoid overflow, LH 8/26/2020
            if ((isc % 100) == 0) {
                cout << "Processing scan " << isc << endl;
            }
            // cout << "Processing scan " << isc << endl;
            // cout << tfileStream.tellg() << endl;

            memcpy(&pbuffer1[0][0], &pbuffer0[0][0], PKTSIZE * maxpkts);
            ancind = ancind0;
            tlmind = tlmind0;
            npkt1 = npkts0;
            spn = spn0;
            enddata = endfile;

            // Read next scan

            read_oci_scan_packets(&tfileStream, fpacket, (uint8_t(*)[PKTSIZE]) & pbuffer0[0][0], npkts0, spn0,
                                  ancind0, tlmind0, seqerr[isc], endfile, isSPW);

            // Disabled maxsc check for I&T
            if (dspn >= 1 && npkt1 > 1) {  // ver 1.03.00

                // Save ancillary packet in array
                if (ancind != -1) {
                    memcpy(ancdata[isc], pbuffer1[ancind], ANCSIZE);
                } else {
                    // insert spin # if missing anc packet
                    uint32_t ui32;
                    ui32 = SWAP_4(spn);
                    memcpy(&ancdata[isc][24], &ui32, 4);

                    // set apid to 0 where there is missing anc packet, LH 4/28/2020
                    ui32 = SWAP_4(0);
                    memcpy(&ancdata[isc][0], &ui32, 4);
                }

                // Save telemetry packet in array
                int ntind = tlmind.size();
                if (ntind > 0 && ntlm < maxtlm) {
                    int itt = 0;
                    while (itt < ntind && ntlm < maxtlm) {
                        memcpy(tlmdata[ntlm], pbuffer1[tlmind[itt]], TLMSIZE);
                        ntlm++;
                        itt++;
                    }
                    if (ntlm >= maxtlm)
                        cout << "Maximum number of telemetry packets exceeded at spin: " << spn << endl;
                }

                // Zero out slines, sfrm, and sbands arrays
                for (size_t i = 0; i < msps; i++) {
                    slines[i] = -1;
                    sfrm[i] = 0;
                    for (size_t j = 0; j < nswb; j++)
                        sbands[i][j] = 0;
                }

                // Unpack science data from packets
                npkts = npkt1 + npkts0;
                if (npkts > maxpkts) {  // LH, 5/27/2022, ver 1.04.01
                    cout << "Packets exceed max [" << maxpkts << "]." << endl;
                    npkts = maxpkts;
                    npkts0 = npkts - npkt1;
                }
                memcpy(&pbuffer[0][0], &pbuffer1[0][0], PKTSIZE * npkt1);
                if (npkts0 > 0)
                    memcpy(&pbuffer[npkt1][0], &pbuffer0[0][0],
                           PKTSIZE * npkts0);  // pbuffer[*,npkt1:npkts-1] = pbuffer2[*,0:npkt2-1]

                unpack_oci_sci(npkts, spn, ncpt, nspt, msps, nbands, btaps, rtaps,
                               (uint8_t(*)[PKTSIZE]) & pbuffer[0][0], bbands, rbands, sbands, blines, rlines,
                               slines, btype, bagg, rtype, ragg, sfrm, icheck);
                if (icheck != 0)
                    cout << "Science data unpacking error in spin: " << spn << endl;
                iret = iret | icheck;

                bdark[0] = dark_b[isc];
                for (size_t i = 1; i < nbbs; i++)
                    bdark[i] = bdark[i - 1] + ndcs;

                rdark[0] = dark_r[isc];
                for (size_t i = 1; i < nrbs; i++)
                    rdark[i] = rdark[i - 1] + ndcs;

                sdark[0] = dark_s[isc];
                for (size_t i = 1; i < nswb; i++)
                    sdark[i] = sdark[i - 1] + ndss;

                // Initialize science/dark array with fill value
                std::fill(bsci[0], bsci[0] + ncps * ((nbbs > 0) ? nbbs : 1), 65535);
                std::fill(rsci[0], rsci[0] + ncps * ((nrbs > 0) ? nrbs : 1), 65535);
                std::fill(ssci[0], ssci[0] + nsps * ((nswb > 0) ? nswb : 1), 1048575);
                memset(sfrms, -1, sizeof(int8_t) * msps);

                // Initialize with fill value, LH 8/6/2020
                std::fill(bdark[0], bdark[0] + ndcs * ((nbbs > 0) ? nbbs : 1), 65535);
                std::fill(rdark[0], rdark[0] + ndcs * ((nrbs > 0) ? nrbs : 1), 65535);
                std::fill(sdark[0], sdark[0] + ndss * ((nswb > 0) ? nswb : 1), 1048575);

                // If science data in scan, check data for gaps or inconsistencies and determine data types
                if ((blines[0] != -1) || (rlines[0] != -1) || (slines[0] != -1)) {
                    check_load_oci_data(
                        dtype, ncpt, nsps, ndcs, ndss, nbbs, nrbs, nswb,  // msps -> nsps, LH, ver 1.03.00
                        cindex, sindex, cdindex, sdindex, bbands, rbands, sbands, blines, rlines, slines,
                        bsci, rsci, ssci, bdark, rdark, sdark, linerr[isc], icheck);
                    if (icheck >= 2)
                        cout << "Science data checking error in spin: " << spn << endl;
                    iret = iret | icheck;

                    // Store SWIR band frame types
                    for (size_t i = 0; i < msps; i++) {
                        if (sindex[slines[i] * nswb] > -1) {  // LH, 8/24/2020; use SWIR band0, LH 11/20/2020
                            sfrms[sindex[slines[i] * nswb]] = sfrm[i];
                        }
                    }

                    for (size_t i = 0; i < msps; i++) {
                        if (sdindex[slines[i]] > -1) {  // LH, 8/24/2020
                            sdfrms[isc][sdindex[slines[i]]] = sfrm[i];
                        }
                    }

                    // Write science data to file
                    outfile.write_oci_science_data(isc, nbbs, nrbs, nswb, ncps, nsps, bsci, rsci, ssci,
                                                   sfrms);

                    //  Check for instrument HKT packets in scan (placeholder for now)

                    isc++;
                    stimp = stime;
                    endtime.iyear = iyear;
                    endtime.iday = iday;
                    endtime.sec = stime;
                    while (endtime.sec > 86400)
                        endtime.sec -= 86400;
                    spnp = spn;
                }
            } else {
                ih = (int32_t)(stime / 3600);
                mn = (int32_t)((stime - ih * 3600) / 60);
                isec = (int32_t)(stime - ih * 3600 - mn * 60);

                if (stime <= stimp && ancind != -1)
                    cout << "Scan " << spn << " out of order at" << ih << " " << mn << " " << isec << endl;
            }  // if (dspn >= 1 && npkts >

            // Get scan time and band dimensions for next scan
            if (!enddata && ancind0 != -1) {
                memcpy(apacket, &pbuffer0[ancind0][0], ANCSIZE);
                get_anc_packet_time(apacket, iyear, iday, stime);
                uint32_t jd = jday(iyear, 1, iday);
                stime += (jd - jd0) * 86400;

                get_swir_mode((uint8_t(*)[PKTSIZE]) & pbuffer0[0][0], npkts0, smode);

                scomp = (smode == smodep);

                if (!scomp) {
                    ih = (int32_t)(stime / 3600);
                    mn = (int32_t)((stime - ih * 3600) / 60);
                    isec = (int32_t)(stime - ih * 3600 - mn * 60);
                    cout << "SWIR data mode change at: " << ih << " " << mn << " " << isec << endl;
                }
                acomp = anc_compare(apacket0, apacket);
                //         cout << "acomp: " << acomp << endl;
            }

            dspn = spn0 - spnp;
            if (dspn > maxgap)
                cout << "Spin number gap at spin  " << spnp << endl;

        }  // while ( stime < mtime && acomp && !enddata)
        // ---------------------------------------------------------------------
        // ---------------------------------------------------------------------
        // ---------------------------------------------------------------------

        if (mper > 0 && isc < (size_t)(maxsc - 10) && acomp) {
            iret = iret | 1;
        }

        cout << "Scans in file: " << isc << endl;
        cout << "Complete flag: " << iret << endl;

        if (isc > 0) {
            // Need to include ancillary packet from next scan
            if (ancind0 != -1)
                memcpy(ancdata[isc], apacket, ANCSIZE);

            // Write calibration data to file
            outfile.write_oci_cal_data(isc, nbbs, nrbs, nswb, ndcs, ndss, dark_b[0], dark_r[0], dark_s[0],
                                       sdfrms[0]);

            // Get scan metadata and write to file
            // Save spid id computed in this function for later use
            int32_t *spinID = new int32_t[isc + 100];
            for (size_t i = 0; i < (isc + 100); i++)
                spinID[i] = BAD_INT;
            outfile.write_oci_scan_metadata(isc, ancdata[0], seqerr, linerr, spinID, starttime);

            // Write ancillary data to file
            outfile.write_oci_ancil_data(isc, ancdata[0]);

            // Unpack OCI telemetry data and write to file
            cdsmode = 0;
            int ntind = tlmind0.size();
            if (tlmind0.size() != 0 && dspn <= maxgap && !endfile) {  // 0.99.20
                int itt = 0;
                while (itt < ntind && ntlm < maxtlm) {
                    memcpy(tlmdata[ntlm], pbuffer0[tlmind0[itt]], TLMSIZE);
                    ntlm++;
                    itt++;
                }
            }

            if (ntlm > 0) {
                cout << ntlm << " HKT packets" << endl;
                outfile.write_oci_tlm_data(itable, ntlm, (uint8_t(*)[TLMSIZE]) & tlmdata[0][0], spinID,
                                           cdsmode, isc, starttime);
            }
            delete[] spinID;

            // Locate navigation data and write to file
            if (hktlist.compare("") != 0)
                outfile.write_navigation(hktlist, starttime, endtime);

            //      if ((endfile && mper <= 0) || !acomp)
            //  mtime = (uint32_t) floor(stimp);

            // Generate granule metadata and write to file
            string sdir = "Ascending";  // Hard-code until we have spacecraft data
            string edir = "Ascending";  // Hard-code until we have spacecraft data
            outfile.write_oci_global_metadata(starttime, endtime, l1a_name, sdir, edir, isc,
                                              dtype,  // itable[1].dtype changed to dtype
                                              smode, cdsmode, fout);

            // Write complete flag to outlist if needed
            // write data type to the last column
            if (outlist.compare("") != 0) {
                fout << iret;

                if (smode > 0) {
                    fout << " DIAG";
                } else if (dtype == 0) {
                    fout << " _";
                } else if (dtypes[dtype].substr(0, 1) == "_") {
                    fout << " " << dtypes[dtype].substr(1);
                }

                fout << "\n";
            }

            // Write common global metadata
            set_global_attrs(outfile.ncfile(), history, doi, pversion);

            // Close file
            outfile.close();

        } else {
            // Remove 0-scan file
            outfile.close();  // LH, 08/19/2020
            int status = remove(l1a_name.c_str());
            if (status == 0) {
                cout << "Removing 0-scan file: " << l1a_name.c_str() << endl;
                if (outlist.compare("") != 0)
                    fout.seekp(outlistpos);
            } else {
                cout << "Error removing " << l1a_name.c_str() << endl;
                exit(1);
            }
        }

        delete[] dark_b[0];
        delete[] dark_b;
        delete[] dark_r[0];
        delete[] dark_r;
        delete[] dark_s[0];
        delete[] dark_s;

        delete[] bdark;
        delete[] rdark;
        delete[] sdark;

        delete[] bsci[0];
        delete[] bsci;
        delete[] rsci[0];
        delete[] rsci;
        delete[] ssci[0];
        delete[] ssci;

        for (size_t i = 0; i < ncpt; i++)
            if (bbands[i] != NULL)
                delete[] bbands[i];
        delete[] bbands;
        bbands = NULL;
        for (size_t i = 0; i < ncpt; i++)
            if (rbands[i] != NULL)
                delete[] rbands[i];
        delete[] rbands;
        rbands = NULL;
        for (size_t i = 0; i < msps; i++)
            if (sbands[i] != NULL)
                delete[] sbands[i];
        delete[] sbands;
        sbands = NULL;

        delete[] blines;
        delete[] rlines;
        delete[] slines;

        delete[] sdfrms[0];
        delete[] sdfrms;

        delete[] sfrm;
        delete[] sfrms;

        delete[] ancdata[0];
        delete[] ancdata;

        delete[] tlmdata[0];
        delete[] tlmdata;

        if (stime > mtime)
            break;
    }  // while (!endfile)

    //  tempOut.close();

    /////////////////// End Main Loop ///////////////////

    if (outlist.compare("") != 0) {
        if (fout.tellp() == 0) {
            fout.close();
            fout.open(outlist.c_str());
            fout << "No L1A file was generated.";  // clear outlist if no L1A is generated.
            return_status = 120;
        }

        fout.close();
    }

    // tfileStream.close();

    // Deallocate
    delete[] pbuffer0[0];
    delete[] pbuffer0;
    delete[] pbuffer1[0];
    delete[] pbuffer1;
    delete[] pbuffer[0];
    delete[] pbuffer;

    delete[] linerr;
    delete[] seqerr;

    return return_status;
}

int make_oci_line_index(itab *itable, int16_t *cindex, int16_t *sindex,
                        // int32_t &ldark) {
                        // cdindex: CCD dark line index array
                        // sdindex: SWIR dark line index array
                        int16_t *cdindex, int16_t *sdindex, int16_t *swir_loff) {
    const uint16_t maxlines = 32768;
    const uint16_t nagg[4] = {1, 2, 4, 8};
    const uint16_t nswb = 9;
    const size_t n_table = 10;
    const size_t n_aggr_size = 8;

    // Loop through data zones
    uint16_t loff = 0;
    uint16_t cpix = 0;
    uint16_t spix = 0;
    uint16_t cdpix = 0;
    uint16_t sdpix = 0;
    int16_t sloff = 0;

    for (size_t i = 0; i < n_table; i++) {
        // If not "no data" type, add indices to array
        if ((itable[i].dtype != 0) && (itable[i].dtype != n_table)) {
            if (itable[i].lines > 0) {
                // CCD pixel index
                if (itable[i].iagg > 3 || itable[i].iagg < 0) {
                    printf(
                        "--Error-- : the value of itable at i = %d is an out of range index %d. See %s at "
                        "%d\n\n",
                        int(i), itable[i].iagg, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
                uint16_t iagg = nagg[itable[i].iagg];
                uint16_t lines = itable[i].lines;
                // Check for total lines within limit
                if ((loff + lines) > maxlines) {
                    cout << "Mode table entry " << i << " exceeds max lines" << endl;
                    lines = maxlines - loff - iagg;
                }

                uint16_t cp = lines / iagg;
                for (size_t j = 0; j < cp; j++) {
                    uint16_t cind = loff + j * iagg;
                    cindex[cind] = cpix + j;
                }
                cpix += cp;

                // SWIR pixel index (fixed aggregation of n_aggr_size)
                uint16_t sp = lines / n_aggr_size;
                for (size_t j = 0; j < sp; j++) {
                    uint16_t sind = (loff / n_aggr_size + j) * n_aggr_size;
                    // sindex[sind] = spix + j;
                    for (size_t k = 0; k < nswb; k++) {
                        sloff = 0;
                        // If not dark view, use line offsets
                        if (itable[i].dtype != 2) {
                            sloff = swir_loff[k];
                        }
                        if ((sind - sloff) < 0) {
                            printf("--Error-- : %s:%d : SWIR line offset caused a negative index\n", __FILE__,
                                   __LINE__);
                            exit(EXIT_FAILURE);
                        }
                        sindex[(sind - sloff) * nswb + k] = spix + j;
                    }
                }
                spix += sp;

                // Dark collect ( if dtype equals 2)
                if (itable[i].dtype == 2) {
                    for (size_t j = 0; j < cp; j++) {
                        uint16_t cind = loff + j * iagg;
                        cdindex[cind] = cdpix + j;
                    }
                    cdpix += cp;

                    for (size_t j = 0; j < sp; j++) {
                        uint16_t sind = loff + j * n_aggr_size;
                        sdindex[sind] = sdpix + j;
                    }
                    sdpix += sp;
                }
            } else {
                cout << "Data zone " << i << " type " << itable[i].dtype << " has zero lines" << endl;
            }
        }
        loff += itable[i].lines;
        if (loff >= maxlines)
            break;
    }

    return 0;
}

int unpack_oci_sci(uint32_t npkts, int32_t spin, uint16_t ncps, uint16_t nsps, uint16_t msps,
                   uint16_t &nbands, uint16_t btaps[16], uint16_t rtaps[16], uint8_t (*pbuffer)[PKTSIZE],
                   uint16_t **bbands, uint16_t **rbands, uint32_t **sbands, int16_t *blines, int16_t *rlines,
                   int16_t *slines, uint16_t &btype, uint16_t bagg[16], uint16_t &rtype, uint16_t ragg[16],
                   int8_t *sfrm, int &iret) {
    // Unpack all of the science data for an OCI scan.
    // The SWIR bands are re-ordered in ascending wavelength order.
    // Order of dual-gain bands is (SG, HG).
    for (size_t i = 0; i < ncps; i++) {
        blines[i] = -1;
        rlines[i] = -1;
    }
    for (size_t i = 0; i < msps; i++)
        slines[i] = -1;

    std::vector<int> iswav = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    int ibpix = 0;
    int irpix = 0;
    int ispix = 0;
    iret = 0;

    uint16_t ccdid;
    uint32_t line;
    uint16_t iagg;
    uint16_t jagg[16];

    uint16_t **ccddata = new uint16_t *;
    uint16_t ossdata[16];

    uint16_t nbbnds;
    uint16_t nrbnds;
    uint16_t nb;

    uint32_t apid;
    uint32_t ui32;
    int32_t spn;

    // tempOut.open ("ccddata.bin", ios::out | ios::trunc | ios::binary);

    // SWIR band sorting indices for science mode
    uint16_t smode = 0;
    get_swir_mode((uint8_t(*)[PKTSIZE]) & pbuffer[0][0], npkts, smode);

    if (smode == 0) {
        iswav = {3, 0, 1, 2, 8, 6, 7, 5, 4};
    }

    for (size_t ipkt = 0; ipkt < npkts; ipkt++) {
        apid = (pbuffer[ipkt][0] % 8) * 256 + pbuffer[ipkt][1];
        uint16_t dtype = (pbuffer[ipkt][12] % 64) / 4;
        memcpy(&ui32, &pbuffer[ipkt][6], 4);
        spn = SWAP_4(ui32);

        // If CCD science APID
        if (apid == 700 && dtype > 0 && spn == spin) {
            unpack_ccd_packet(pbuffer[ipkt], btaps, rtaps, ccdid, line, dtype, iagg, jagg, nbands, ccddata,
                              ossdata);

            //  tempOut.write((char *) &nbands, sizeof(uint16_t));
            // tempOut.write((char *) &(*ccddata)[0], nbands*sizeof(uint16_t));

            // If blue
            if (ccdid) {
                if (ibpix == 0) {
                    nbbnds = nbands;
                    btype = dtype;
                    memcpy(bagg, jagg, 16 * sizeof(uint16_t));
                    for (size_t i = 0; i < ncps; i++) {
                        if (bbands[i])
                            delete[] bbands[i];
                        bbands[i] = new uint16_t[nbands];
                        for (size_t j = 0; j < nbands; j++)
                            bbands[i][j] = 65535;  // LH, 11/23/2020
                    }

                } else {
                    // Check for inconsistent non-dark data type or spectral aggregation
                    uint8_t agg = 0;
                    for (size_t i = 0; i < 16; i++)
                        if (jagg[i] != bagg[i])
                            agg++;

                    if ((dtype != btype && dtype != 2 && dtype != 5) ||
                        (agg > 0)) {  // Ignore dark and linearity
                        cout << "Data type or spectral aggregation error, CCDID: " << ccdid
                             << "  Line: " << line << endl;
                        iret = 4;
                    }
                }
                if (nbands <= nbbnds)
                    nb = nbands;
                else
                    nb = nbbnds;
                if (ibpix < ncps) {
                    memcpy(bbands[ibpix], &(*ccddata)[0], nb * sizeof(uint16_t));
                    blines[ibpix] = (line / iagg) * iagg;
                } else {
                    // memcpy( &ui32, &pbuffer[ipkt][6], 4);
                    // int32_t spn = SWAP_4( ui32);
                    cout << "Number of blue pixels exceeded in spin: " << spin
                         << " in packet (0-based): " << ipkt << endl;
                }
                ibpix++;

            } else {
                if (irpix == 0) {
                    nrbnds = nbands;
                    rtype = dtype;
                    memcpy(ragg, jagg, 16 * sizeof(uint16_t));
                    for (size_t i = 0; i < ncps; i++) {
                        if (rbands[i])
                            delete[] rbands[i];
                        rbands[i] = new uint16_t[nbands];
                        for (size_t j = 0; j < nbands; j++)
                            rbands[i][j] = 65535;  // LH, 11/23/2020
                    }
                } else {
                    // Check for inconsistent non-dark data type or spectral aggregation
                    uint8_t agg = 0;
                    for (size_t i = 0; i < 16; i++)
                        if (jagg[i] != ragg[i])
                            agg++;

                    if ((dtype != rtype && dtype != 2 && dtype != 5) ||
                        (agg > 0)) {  // Ignore dark and linearity
                        cout << "Data type or spectral aggregation error, CCDID: " << ccdid
                             << "  Line: " << line << endl;
                        iret = 4;
                    }
                }

                if (nbands <= nrbnds)
                    nb = nbands;
                else
                    nb = nrbnds;
                if (irpix < ncps) {
                    memcpy(rbands[irpix], &(*ccddata)[0], nb * sizeof(uint16_t));
                    rlines[irpix] = (line / iagg) * iagg;
                } else {
                    // memcpy( &ui32, &pbuffer[ipkt][6], 4);
                    // int32_t spn = SWAP_4( ui32);
                    cout << "Number of red pixels exceeded in spin: " << spin
                         << " in packet (0-based): " << ipkt << endl;
                    // cout << "irpix="<<irpix<<", ncps="<<ncps<<endl;
                    // exit(1);  // ver 1.03.02
                }
                irpix++;
            }  // if (ccdid)

            // cout << "delete ccddata" << endl;
            delete[] *ccddata;

        }  // if (apid == 700 && dtype > 0)

        // If SWIR science APID
        if (apid == 720 && (ispix + 7) < msps && spn == spin) {
            int16_t lines[8];
            uint8_t swirfrm[8];
            uint32_t swirdata[8 * 9];

            // swirdata: 8 rows of 9 columns
            unpack_swir_packet(pbuffer[ipkt], lines, swirfrm, swirdata);

            for (size_t i = 0; i < 8; i++) {
                slines[ispix + i] = (lines[i] / 8) * 8;
                for (size_t j = 0; j < 9; j++)
                    sbands[ispix + i][j] = swirdata[iswav[j] + 9 * i];  // bands in ascending order
            }
            memcpy(&sfrm[ispix], swirfrm, 8 * sizeof(uint8_t));

            ispix += 8;
        }  // if (apid == 720 && (ispix+7) < msps)

        if ((ibpix <= 0) || (ibpix > ncps))
            blines[0] = -1;
        if ((irpix <= 0) || (irpix > ncps))
            rlines[0] = -1;
        if ((ispix <= 0) || (ispix > msps))
            slines[0] = -1;

    }  // ipkt loop

    delete ccddata;
    //  tempOut.close();

    return 0;
}

int unpack_ccd_packet(uint8_t *packet, uint16_t btaps[16], uint16_t rtaps[16], uint16_t &ccdid,
                      uint32_t &line, uint16_t &dtype, uint16_t &iagg, uint16_t jagg[16], uint16_t &nbands,
                      uint16_t **ccddata, uint16_t ossdata[16]) {
    // Get CCD ID, line number, data type and spatial aggregation
    ccdid = (packet[12] & 64) / 64;
    line = packet[10] * 256 + packet[11];
    dtype = (packet[12] % 64) / 4;
    uint16_t oss = packet[17] % 16;
    iagg = packet[12] % 4;
    uint16_t agg[4] = {1, 2, 4, 8};
    iagg = agg[iagg];

    // Get tap enable flags for focal plane
    uint16_t *ftaps;
    if (ccdid)
        ftaps = btaps;
    else
        ftaps = rtaps;

    // Get spectral aggregation factors for all taps
    uint16_t its[4] = {0, 4, 8, 12};
    uint32_t taps[16];
    for (size_t i = 0; i < 4; i++) {
        jagg[its[i]] = agg[(packet[13 + i] & 192) / 64];
        jagg[its[i] + 1] = agg[(packet[13 + i] & 48) / 16];
        jagg[its[i] + 2] = agg[(packet[13 + i] & 12) / 4];
        jagg[its[i] + 3] = agg[(packet[13 + i] & 3)];
    }
    for (size_t i = 0; i < 16; i++)
        taps[i] = 32 * ftaps[i] / jagg[i];

    // Allocate output buffer
    nbands = 0;
    for (size_t i = 0; i < 16; i++)
        nbands += taps[i];

    // cout << "new ccddata" << endl;
    *ccddata = new uint16_t[nbands];

    int ioff = 18;
    uint16_t ibnd = nbands - 1;

    // For each tap (1-16):
    for (size_t j = 0; j < 16; j++) {
        // Copy band data from packet to output buffer
        // Packet data are in reverse spectral order, so swap order here
        if (ftaps[j]) {
            for (size_t i = 0; i < taps[j]; i++) {
                uint16_t ui16;
                memcpy(&ui16, &packet[ioff + 2 * i], 2);
                ui16 = SWAP_2(ui16);
                memcpy(&(*ccddata)[ibnd - i], &ui16, 2);
            }
            ibnd -= taps[j];
            ioff += 2 * taps[j];
        }
    }

    if (oss > 0) {
        for (size_t j = 0; j < 16; j++) {
            uint16_t ui16;
            memcpy(&ui16, &packet[ioff + 2 * j], 2);
            ui16 = SWAP_2(ui16);
            memcpy(&ossdata[j], &ui16, 2);
        }
    }

    return 0;
}

int unpack_swir_packet(uint8_t *packet, int16_t *slines, uint8_t *swirfrm, uint32_t *swirdata) {
    // Program to unpack one OCI SWIR data packet
    // Each packet contains data for eight science pixels
    // Reference: OCI-ELEC-SPEC-0028

    // Allocate output buffers (9 bands for 8 pixels)

    int ioff = 10;
    uint8_t smeta;
    ;

    for (size_t i = 0; i < 8; i++) {
        // Get line number and metadata
        if ((unsigned(packet[ioff]) == 255) && (unsigned(packet[ioff + 1]) == 255)) {
            // handle fill values in SWIR, LH 10/28/2020
            slines[i] = 0;
        } else {
            slines[i] = packet[ioff] * 256 + packet[ioff + 1];
        }

        smeta = packet[ioff + 2];
        swirfrm[i] = smeta % 8;

        // Extract SWIR band data
        eight20(&packet[ioff + 3], &swirdata[9 * i]);

        ioff += 26;
    }

    return 0;
}

int check_load_oci_data(short dtype, uint16_t ncpt, uint16_t nspt0, uint16_t ndcs, uint16_t ndss,
                        uint16_t nbbs, uint16_t nrbs, uint16_t nswb, int16_t *cindex, int16_t *sindex,
                        int16_t *cdindex, int16_t *sdindex, uint16_t **bbands, uint16_t **rbands,
                        uint32_t **sbands, int16_t *blines, int16_t *rlines, int16_t *slines, uint16_t **bsci,
                        uint16_t **rsci, uint32_t **ssci, uint16_t **bdark, uint16_t **rdark,
                        uint32_t **sdark, int8_t &linerr, int &icheck) {
    icheck = 0;

    // CCD line check except in linearity mode
    bool linchk = (dtype != 5);

    // Check for presence of blue and red focal plane data
    // dbb = size(bbands)
    bool bb = ((bbands[0] != NULL) && (blines[0] != -1));
    // drb = size(rbands)
    bool rb = ((rbands[0] != NULL) && (rlines[0] != -1));
    // drs = size(sbands)
    bool sb = ((sbands[0] != NULL) && (slines[0] != -1));

    //////////////////// Blue Bands ////////////////////
    uint16_t nbb;
    if (bbands[0] == NULL)
        nbb = 0;
    else
        nbb = nbbs;
    int mb = -1;

    // Check for invalid line numbers
    if (bb && linchk) {
        for (size_t i = 0; i < ncpt; i++) {
            if (blines[i] == -1)
                continue;
            if (cindex[blines[i]] == -1)
                mb++;
        }
        for (size_t i = 0; i < ncpt; i++)
            if (cindex[blines[i]] != (cindex[blines[0]] + (int)i)) {
                cout << "Blue CCD line sequence error" << endl;
                linerr = 1;
                break;
            }
    }

    // Identify pixels and load into output arrays
    if (bb) {
        uint16_t nbs = 0;
        for (size_t i = 0; i < ncpt; i++) {
            if (blines[i] == -1)
                continue;
            // if (blines[i] < ldark && cindex[blines[i]] != -1) nbs++;
            if (cindex[blines[i]] > -1)
                nbs++;  // LH, 8/24/2020, replaced != with >
        }
        if ((nbs < ncpt - ndcs) && linchk) {
            cout << "Missing blue band science pixels" << endl;
            icheck = 2;
        }

        for (size_t i = 0; i < ncpt; i++) {
            if (blines[i] == -1)
                continue;
            // if (blines[i] < ldark && cindex[blines[i]] != -1) {
            if (cindex[blines[i]] > -1) {  // LH, 8/24/2020
                for (size_t j = 0; j < nbb; j++) {
                    bsci[j][cindex[blines[i]]] = bbands[i][j];
                }
            }
        }

        // Identify dark count pixels and load into output arrays
        uint16_t nbd = 0;
        for (size_t i = 0; i < ncpt; i++) {
            if (blines[i] == -1)
                continue;
            // if (blines[i] >= ldark && cindex[blines[i]] != -1) nbd++;
            if (cdindex[blines[i]] > -1)
                nbd++;  // LH, 8/24/2020
        }
        if (nbd < ndcs) {
            cout << "Missing blue band dark pixels" << endl;
            icheck = 2;
        }

        for (size_t i = 0; i < ncpt; i++) {
            if (blines[i] == -1)
                continue;
            // if (blines[i] >= ldark && cindex[blines[i]] != -1) {
            if (cdindex[blines[i]] > -1) {  // LH, 8/24/2020
                // int k = cindex[blines[i]]-cindex[ldark];
                int k = cdindex[blines[i]];
                if (k >= ndcs || k < 0) {
                    // cout << "bdark indice out of bounds: " << k << endl;
                    continue;
                    //          exit(1);
                }
                for (size_t j = 0; j < nbb; j++) {
                    bdark[j][k] = bbands[i][j];
                }
            }
        }
    }

    //////////////////// Red Bands ////////////////////
    uint16_t nrb;
    if (rbands[0] == NULL)
        nrb = 0;
    else
        nrb = nrbs;
    int mr = -1;

    // Check for invalid line numbers
    if (rb && linchk) {
        for (size_t i = 0; i < ncpt; i++) {
            if (rlines[i] == -1)
                continue;
            if (cindex[rlines[i]] == -1)
                mr++;
        }
        for (size_t i = 0; i < ncpt; i++) {
            if (rlines[i] == -1)
                continue;
            if (cindex[rlines[i]] != (cindex[rlines[0]] + (int)i)) {
                cout << "Red CCD line sequence error" << endl;
                linerr = 1;
                break;
            }
        }
    }

    // Identify pixels and load into output arrays
    if (rb) {
        uint16_t nrs = 0;
        for (size_t i = 0; i < ncpt; i++) {
            if (rlines[i] == -1)
                continue;
            // if (rlines[i] < ldark && cindex[rlines[i]] != -1) nrs++;
            if (cindex[rlines[i]] > -1)
                nrs++;  // LH, 8/24/2020
        }

        if ((nrs < ncpt - ndcs) && linchk) {
            cout << "Missing red band science pixels" << endl;
            icheck = 2;
        }

        for (size_t i = 0; i < ncpt; i++) {
            // if (rlines[i] < ldark && cindex[rlines[i]] != -1) {
            if (rlines[i] == -1)
                continue;
            if (rlines[i] > -1 && cindex[rlines[i]] > -1) {  // LH, 8/24/2020
                for (size_t j = 0; j < nrb; j++) {
                    rsci[j][cindex[rlines[i]]] = rbands[i][j];
                }
            }
        }

        // Identify dark count pixels and load into output arrays
        uint16_t nrd = 0;
        for (size_t i = 0; i < ncpt; i++) {
            if (rlines[i] == -1)
                continue;
            // if (rlines[i] >= ldark && cindex[rlines[i]] != -1) nrd++;
            if (cdindex[rlines[i]] > -1)
                nrd++;  // LH, 8/24/2020
        }

        if (nrd < ndcs) {
            cout << "Missing red band dark pixels" << endl;
            icheck = 2;
        }

        for (size_t i = 0; i < ncpt; i++) {
            if (rlines[i] == -1)
                continue;
            // if (rlines[i] >= ldark && cindex[rlines[i]] != -1) {
            if (cdindex[rlines[i]] > -1) {  // LH, 8/24/2020
                // int k = cindex[rlines[i]]-cindex[ldark];
                int k = cdindex[rlines[i]];
                if (k >= ndcs || k < 0) {
                    // cout << "rdark indice out of bounds: " << k << endl;
                    continue;
                    //          exit(1);
                }
                for (size_t j = 0; j < nrb; j++) {
                    rdark[j][k] = rbands[i][j];
                }
            }
        }
    }

    if (mb != -1 || mr != -1) {
        // Disable SWIR line check for ETU
        cout << "Invalid line numbers" << endl;
        icheck = 4;
    }

    //////////////////// SWIR Bands ////////////////////

    // Identify pixels and load into output arrays
    // SWIR science data line index has band-by-band offsets, LH, 11/20/2020
    if (sb) {
        uint16_t nss = 0;
        for (size_t i = 0; i < nspt0; i++) {
            if (slines[i] > -1) {
                nss++;
                for (size_t j = 0; j < nswb; j++) {
                    if (sindex[slines[i] * nswb + j] >
                        -1) {  // LH, 8/24/2020; sindex expanded by * nswb bands, LH 11/20/2020
                        ssci[j][sindex[slines[i] * nswb + j]] = sbands[i][j];
                    }
                }
            }
        }

        if (nss < nspt0 - ndss) {
            cout << "Missing SWIR band science pixels" << endl;
            icheck = 2;
        }

        // Identify dark count pixels and load into output arrays
        uint16_t nsd = 0;
        for (size_t i = 0; i < nspt0; i++) {
            if (slines[i] == -1)
                continue;
            // if (slines[i] >= ldark && sindex[slines[i]] != -1) nsd++;
            if (sdindex[slines[i]] > -1)
                nsd++;  // LH, 8/24/2020
        }
        // if (nsd < ndss) {
        //   cout << "Missing SWIR band dark pixels" << endl;
        //   icheck = 1;
        // }

        for (size_t i = 0; i < nspt0; i++) {
            // if (slines[i] >= ldark && sindex[slines[i]] != -1) {
            if (slines[i] == -1)
                continue;
            if (sdindex[slines[i]] > -1) {  // LH, 8/24/2020
                // int k = sindex[slines[i]]-sindex[ldark];
                int k = sdindex[slines[i]];
                if (k >= ndss || k < 0) {
                    cout << "sdark indice out of bounds: " << k << endl;
                    exit(1);
                }

                for (size_t j = 0; j < nswb; j++) {
                    sdark[j][k] = sbands[i][j];
                }
            }
        }
    }
    return 0;
}

uint8_t check_sum(int32_t nc, uint8_t *dat, uint8_t *chk) {
    // Function to check data against checksum
    // Checksum is 4 bytes computed by XOR

    uint8_t chks[4], tmp[4];
    memcpy(chks, &dat[0], 4);

    for (int i = 1; i < nc; i++) {
        memcpy(&tmp, &dat[4 * i], 4);
        for (int j = 0; j < 4; j++)
            chks[j] = chks[j] ^ tmp[j];
    }

    uint8_t check[4];
    for (int i = 0; i < 4; i++)
        check[i] = (chk[i] == chks[i]);

    uint8_t check_sum = (check[0] & check[1] & check[2] & check[3]);

    return check_sum;
}

/*----------------------------------------------------------------- */
/* Create an Generic NETCDF4 level1 file                            */
/* ---------------------------------------------------------------- */
int l1aFile::createl1(char *l1_filename, uint16_t maxsc, uint16_t ncps, uint16_t nbbs, uint16_t nrbs,
                      uint16_t nsps, uint16_t ndcs) {
    try {
        l1afile = new NcFile(l1_filename, NcFile::replace);
    } catch (NcException &e) {
        cerr << e.what() << "\nFailure creating OCI L1A file: " << l1_filename << endl;
        exit(1);
    }

    fileName.assign(l1_filename);

    ifstream oci_l1a_data_structure;
    string line;
    string dataStructureFile;
    dataStructureFile.assign("$OCDATAROOT/oci/OCI_Level-1A_Data_Structure.cdl");
    expandEnvVar(&dataStructureFile);

    oci_l1a_data_structure.open(dataStructureFile.c_str(), ifstream::in);
    if (oci_l1a_data_structure.fail() == true) {
        cout << "\"" << dataStructureFile.c_str() << "\" not found" << endl;
        exit(1);
    }

    // Find "dimensions" section of CDL file
    while (1) {
        getline(oci_l1a_data_structure, line);
        boost::trim(line);
        size_t pos = line.find("dimensions:");
        if (pos == 0)
            break;
    }

    // Define dimensions from "dimensions" section of CDL file
    ndims = 0;

    while (1) {
        getline(oci_l1a_data_structure, line);
        boost::trim(line);
        if (line.substr(0, 2) == "//")
            continue;

        size_t pos = line.find(" = ");
        size_t semi = line.find(" ;");
        if (pos == string::npos)
            break;

        uint32_t dimSize;
        istringstream iss;
        //     istringstream iss(line.substr(pos+2, string::npos));
        string dimString = line.substr(pos + 2, semi - (pos + 2));
        if (dimString.find("UNLIMITED") == string::npos) {
            iss.str(dimString);
            iss >> dimSize;
        } else {
            dimSize = NC_UNLIMITED;
        }

        iss.clear();
        iss.str(line);
        iss >> skipws >> line;

        //     cout << "Dimension Name: " << line.c_str() << " Dimension Size: "
        //	  << dimSize << endl;

        if (line.compare("ccd_pixels") == 0) {
            dimSize = ncps;
        }

        if (line.compare("SWIR_pixels") == 0) {
            dimSize = nsps;
        }

        if (line.compare("DC_pixels") == 0) {
            dimSize = ndcs;
        }

        if (line.compare("blue_bands") == 0) {
            dimSize = nbbs;
        }

        if (line.compare("red_bands") == 0) {
            dimSize = nrbs;
        }

        try {
            ncDims[ndims++] = l1afile->addDim(line, dimSize);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure creating dimension: " << line.c_str() << endl;
            exit(1);
        }

    }  // while loop

    // Read global attributes (string attributes only)
    while (1) {
        getline(oci_l1a_data_structure, line);
        boost::trim(line);
        size_t pos = line.find("// global attributes");
        if (pos == 0)
            break;
    }

    while (1) {
        getline(oci_l1a_data_structure, line);
        size_t pos = line.find(" = ");
        if (pos == string::npos)
            break;

        string attValue = line.substr(pos + 3);

        // Remove any leading and trailing quotes
        attValue.erase(attValue.length() - 2);  // skip final " ;"
        size_t begQuote = attValue.find('"');
        size_t endQuote = attValue.find_last_of('"');
        if (begQuote == string::npos)
            continue;  // Skip non-string global attributes
        attValue = attValue.substr(begQuote + 1, endQuote - begQuote - 1);

        istringstream iss(line.substr(pos + 2));
        iss.clear();
        iss.str(line);
        iss >> skipws >> line;

        // Skip commented out attributes
        if (line.substr(0, 2) == "//")
            continue;

        string attName;
        attName.assign(line.substr(1).c_str());

        try {
            l1afile->putAtt(attName, attValue);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure creating attribute: " + attName << endl;
            exit(1);
        }

    }  // while(1)

    ngrps = 0;
    // Loop through groups
    while (1) {
        getline(oci_l1a_data_structure, line);

        // Check if end of CDL file
        // If so then close CDL file and return
        if (line.substr(0, 1).compare("}") == 0) {
            oci_l1a_data_structure.close();
            return 0;
        }

        // Check for beginning of new group
        size_t pos = line.find("group:");

        // If found then create new group and variables
        if (pos == 0) {
            // Parse group name
            istringstream iss(line.substr(6, string::npos));
            iss >> skipws >> line;

            ncGrps[ngrps++] = l1afile->addGroup(line);

            int numDims = 0;
            string sname;
            string lname;
            string standard_name;
            string units;
            string flag_values;
            string flag_meanings;
            string reference;
            double valid_min = 0.0;
            double valid_max = 0.0;
            double fill_value = 0.0;

            vector<NcDim> varVec;

            int ntype = 0;
            NcType ncType;

            // Loop through datasets in group
            getline(oci_l1a_data_structure, line);  // skip "variables:"
            while (1) {
                getline(oci_l1a_data_structure, line);
                boost::trim(line);

                if (line.substr(0, 2) == "//")
                    continue;
                if (line.length() == 0)
                    continue;
                if (line.substr(0, 1).compare("\r") == 0)
                    continue;
                if (line.substr(0, 1).compare("\n") == 0)
                    continue;

                size_t pos = line.find(":");

                // No ":" found, new dataset or empty line or end-of-group
                if (pos == string::npos) {
                    if (numDims > 0) {
                        // Create previous dataset

                        createNCDF(ncGrps[ngrps - 1], sname.c_str(), lname.c_str(), standard_name.c_str(),
                                   units.c_str(), (void *)&fill_value, flag_values.c_str(),
                                   flag_meanings.c_str(), reference.c_str(), valid_min, valid_max, ntype,
                                   varVec);

                        flag_values.assign("");
                        flag_meanings.assign("");
                        reference.assign("");
                        units.assign("");
                        varVec.clear();
                    }

                    valid_min = 0.0;
                    valid_max = 0.0;
                    fill_value = 0.0;

                    if (line.substr(0, 10).compare("} // group") == 0)
                        break;

                    // Parse variable type
                    string varType;
                    istringstream iss(line);
                    iss >> skipws >> varType;

                    // Get corresponding NC variable type
                    if (varType.compare("char") == 0)
                        ntype = NC_CHAR;
                    else if (varType.compare("byte") == 0)
                        ntype = NC_BYTE;
                    else if (varType.compare("short") == 0)
                        ntype = NC_SHORT;
                    else if (varType.compare("int") == 0)
                        ntype = NC_INT;
                    else if (varType.compare("long") == 0)
                        ntype = NC_INT;
                    else if (varType.compare("float") == 0)
                        ntype = NC_FLOAT;
                    else if (varType.compare("real") == 0)
                        ntype = NC_FLOAT;
                    else if (varType.compare("double") == 0)
                        ntype = NC_DOUBLE;
                    else if (varType.compare("ubyte") == 0)
                        ntype = NC_UBYTE;
                    else if (varType.compare("ushort") == 0)
                        ntype = NC_USHORT;
                    else if (varType.compare("uint") == 0)
                        ntype = NC_UINT;
                    else if (varType.compare("ulong") == 0)
                        ntype = NC_UINT;
                    else if (varType.compare("int64") == 0)
                        ntype = NC_INT64;
                    else if (varType.compare("uint64") == 0)
                        ntype = NC_UINT64;

                    // Parse short name (sname)
                    pos = line.find("(");
                    size_t posSname = line.substr(0, pos).rfind(" ");
                    sname.assign(line.substr(posSname + 1, pos - posSname - 1));
                    //           cout << "sname: " << sname.c_str() << endl;

                    // Parse variable dimension info
                    this->parseDims(line.substr(pos + 1, string::npos), varVec);
                    numDims = varVec.size();

                } else {
                    // Parse variable attributes
                    size_t posEql = line.find("=");
                    size_t pos1qte = line.find("\"");
                    size_t pos2qte = line.substr(pos1qte + 1, string::npos).find("\"");
                    // cout << line.substr(pos+1, posEql-pos-2).c_str() << endl;

                    string attrName = line.substr(pos + 1, posEql - pos - 2);

                    // Get long_name
                    if (attrName.compare("long_name") == 0) {
                        lname.assign(line.substr(pos1qte + 1, pos2qte));
                        //             cout << "lname: " << lname.c_str() << endl;
                    }

                    // Get units
                    else if (attrName.compare("units") == 0) {
                        units.assign(line.substr(pos1qte + 1, pos2qte));
                        //             cout << "units: " << units.c_str() << endl;
                    }

                    // Get _FillValue
                    else if (attrName.compare("_FillValue") == 0) {
                        iss.clear();
                        iss.str(line.substr(posEql + 1, string::npos));
                        iss >> fill_value;
                    }

                    // Get flag_values
                    else if (attrName.compare("flag_values") == 0) {
                        flag_values.assign(line.substr(pos1qte + 1, pos2qte));
                    } else if (attrName.compare("flag_masks") == 0) {
                        flag_values.assign(line.substr(pos1qte + 1, pos2qte));
                    }

                    // Get flag_meanings
                    else if (attrName.compare("flag_meanings") == 0) {
                        flag_meanings.assign(line.substr(pos1qte + 1, pos2qte));
                    }

                    // Get reference
                    else if (attrName.compare("reference") == 0) {
                        reference.assign(line.substr(pos1qte + 1, pos2qte));
                    }

                    // Get valid_min
                    else if (attrName.compare("valid_min") == 0) {
                        iss.clear();
                        iss.str(line.substr(posEql + 1, string::npos));
                        iss >> valid_min;
                        //             cout << "valid_min: " << valid_min << endl;
                    }

                    // Get valid_max
                    else if (attrName.compare("valid_max") == 0) {
                        iss.clear();
                        iss.str(line.substr(posEql + 1, string::npos));
                        iss >> valid_max;
                        //             cout << "valid_max: " << valid_max << endl;
                    }

                }  // if ( pos == string::npos)
            }      // datasets in group loop
        }          // New Group loop
    }              // Main Group loop

    return 0;
}

int l1aFile::parseDims(string dimString, vector<NcDim> &varDims) {
    size_t curPos = 0;
    //  char dimName[NC_MAX_NAME+1];
    string dimName;

    while (1) {
        size_t pos = dimString.find(",", curPos);
        if (pos == string::npos)
            pos = dimString.find(")");

        string varDimName;
        istringstream iss(dimString.substr(curPos, pos - curPos));
        iss >> skipws >> varDimName;

        for (int i = 0; i < ndims; i++) {
            try {
                dimName = ncDims[i].getName();
            } catch (NcException &e) {
                e.what();
                cerr << "Failure accessing dimension: " + dimName << endl;
                exit(1);
            }

            if (varDimName.compare(dimName) == 0) {
                //    cout << "     " << dimName << " " << ncDims[i].getSize() << endl;
                varDims.push_back(ncDims[i]);
                break;
            }
        }
        if (dimString.substr(pos, 1).compare(")") == 0)
            break;

        curPos = pos + 1;
    }

    return 0;
}

int l1aFile::write_oci_science_data(uint32_t isc, uint16_t nbbs, uint16_t nrbs, uint16_t nswb, uint16_t ncps,
                                    uint16_t nsps, uint16_t **bsci, uint16_t **rsci, uint32_t **ssci,
                                    int8_t *sfrms) {
    // Writes one scan at a time

    NcVar blu_bands;
    NcVar red_bands;
    NcVar swir_bands;
    NcVar swir_frms;
    vector<size_t> start;
    vector<size_t> count_b;
    vector<size_t> count_r;
    vector<size_t> count_s;
    vector<size_t> count_0;  // for 1-D fill value variables

    vector<size_t> start_frms;
    vector<size_t> count_frms;

    NcGroup gid = l1afile->getGroup("science_data");
    blu_bands = gid.getVar("sci_blue");
    red_bands = gid.getVar("sci_red");
    swir_bands = gid.getVar("sci_SWIR");
    swir_frms = gid.getVar("frm_type_SWIR");

    start.push_back(0);
    start.push_back(0);
    start.push_back(0);

    count_b.push_back(1);
    count_b.push_back(nbbs);
    count_b.push_back(ncps);

    count_r.push_back(1);
    count_r.push_back(nrbs);
    count_r.push_back(ncps);

    count_s.push_back(1);
    count_s.push_back(nswb);
    count_s.push_back(nsps);

    count_0.push_back(1);
    count_0.push_back(1);
    count_0.push_back(1);

    // NcDim dim = blu_bands.getDim(0);
    // cout << dim.getName() << endl;
    // cout << dim.getSize() << endl;
    // cout << dim.isUnlimited() << endl;

    start[0] = isc;

    // create 1-D sci_blue/sci_red/sci_SWIR with fill value if number of blue bands is 0
    (count_b[1] > 0) ? blu_bands.putVar(start, count_b, &bsci[0][0])
                     : blu_bands.putVar(start, count_0, &bsci[0][0]);
    (count_r[1] > 0) ? red_bands.putVar(start, count_r, &rsci[0][0])
                     : red_bands.putVar(start, count_0, &rsci[0][0]);
    (count_s[1] > 0) ? swir_bands.putVar(start, count_s, &ssci[0][0])
                     : swir_bands.putVar(start, count_0, &ssci[0][0]);

    start_frms.push_back(isc);
    start_frms.push_back(0);

    count_frms.push_back(1);
    count_frms.push_back(nsps);

    if (count_frms[1] > 0)
        swir_frms.putVar(start_frms, count_frms, sfrms);

    return 0;
}

int l1aFile::write_oci_cal_data(uint32_t isc, uint16_t nbbs, uint16_t nrbs, uint16_t nswb, uint16_t ndcs,
                                uint16_t ndss, uint16_t *dark_b, uint16_t *dark_r, uint32_t *dark_s,
                                int8_t *sdfrms) {
    NcVar blu_dark;
    NcVar red_dark;
    NcVar swr_dark;
    vector<size_t> start;
    vector<size_t> count_b;
    vector<size_t> count_r;
    vector<size_t> count_s;
    vector<size_t> count_0;  // for 1-D fill value variables

    start.push_back(0);
    start.push_back(0);
    start.push_back(0);

    count_b.push_back(isc);
    count_b.push_back(nbbs);
    count_b.push_back(ndcs);

    count_r.push_back(isc);
    count_r.push_back(nrbs);
    count_r.push_back(ndcs);

    count_s.push_back(isc);
    count_s.push_back(nswb);
    count_s.push_back(ndss);

    count_0.push_back(1);
    count_0.push_back(1);
    count_0.push_back(1);

    NcGroup gid = l1afile->getGroup("onboard_calibration_data");
    blu_dark = gid.getVar("DC_blue");
    red_dark = gid.getVar("DC_red");
    swr_dark = gid.getVar("DC_SWIR");

    // create 1-D DC_blue/DC_red/DC_SWIR with fill value if number of blue bands is 0
    (count_b[1] > 0) ? blu_dark.putVar(start, count_b, dark_b) : blu_dark.putVar(start, count_0, dark_b);
    (count_r[1] > 0) ? red_dark.putVar(start, count_r, dark_r) : red_dark.putVar(start, count_0, dark_r);
    (count_s[1] > 0) ? swr_dark.putVar(start, count_s, dark_s) : swr_dark.putVar(start, count_0, dark_s);

    NcVar sdfrm_dark;
    vector<size_t> start_sdfrm;
    vector<size_t> count_sdfrm;

    start_sdfrm.push_back(0);
    start_sdfrm.push_back(0);
    count_sdfrm.push_back(isc);
    count_sdfrm.push_back(ndss);

    sdfrm_dark = gid.getVar("frm_type_DC_SWIR");
    // if ( count_sdfrm[1] > 0 && sdfrm_dark.isNull())
    if (count_sdfrm[1] > 0)
        sdfrm_dark.putVar(start_sdfrm, count_sdfrm, sdfrms);

    return 0;
}

int l1aFile::write_oci_scan_metadata(uint32_t isc, uint8_t *ancdata, uint8_t *seqerr, int8_t *linerr,
                                     int32_t *spinID, time_struct &starttime) {
    vector<size_t> start;
    vector<size_t> count;
    start.push_back(0);
    count.push_back(isc);

    uint32_t ancap = 636;

    int32_t iyear, iday;
    double stime;

    // Extract and convert times to seconds of day
    // Extract CCSDS scan start times
    double *stimes = new double[isc];
    short int toff = 24;
    uint32_t *scss = new uint32_t[isc];
    int32_t *scsu = new int32_t[isc];

    double firstGoodTime = BAD_FLT;
    for (size_t i = 0; i < isc; i++) {
        uint32_t apid = (ancdata[i * ANCSIZE] % 8) * 256 + ancdata[i * ANCSIZE + 1];
        if (apid == ancap) {
            get_anc_packet_time(&ancdata[i * ANCSIZE], iyear, iday, stime);
            if (firstGoodTime == BAD_FLT)
                firstGoodTime = stime;
            if (stime < firstGoodTime)
                stimes[i] = stime + SECONDS_IN_DAY;  // Ensure that stimes[i] stays relative to the same day
                                                     // that the granule started
            else
                stimes[i] = stime;

            uint32_t ui32;
            memcpy(&ui32, &ancdata[i * ANCSIZE + toff + 4], 4);
            scss[i] = SWAP_4(ui32);
            memcpy(&ui32, &ancdata[i * ANCSIZE + toff + 8], 4);
            scsu[i] = SWAP_4(ui32) / 4096;
        } else {
            stimes[i] = BAD_FLT;
            scss[i] = 0;
            scsu[i] = BAD_INT;
        }
    }

    NcGroup gid = l1afile->getGroup("scan_line_attributes");

    NcVar var;

    // Scan start time (seconds of day)
    var = gid.getVar("scan_start_time");
    var.putVar(start, count, stimes);
    double tmpTime = yds2unix(starttime.iyear, starttime.iday, starttime.sec);
    string timeStr = unix2isodate(tmpTime, 'G');
    timeStr = timeStr.substr(0, 10);
    timeStr.insert(0, "seconds since ");
    var.putAtt("units", timeStr);

    // Scan start time (CCSDS)

    // Seconds since 1958
    var = gid.getVar("scan_start_CCSDS_sec");
    var.putVar(start, count, scss);

    // Microseconds
    var = gid.getVar("scan_start_CCSDS_usec");
    var.putVar(start, count, scsu);

    // Extract and write HAM side
    uint8_t *hamSide = new uint8_t[isc];

    for (size_t i = 0; i < isc; i++) {
        // This part is modified to fix bogus values for the last two records in last chunk of L1A
        // Liang Hong, 4/24/2020
        // int ham =
        //  (ancdata[(i+1)*ANCSIZE+toff+86] % 128) * 256 +
        //  ancdata[(i+1)*ANCSIZE+toff+87];
        int ham = 255;
        // ham /= 16384; // Need to revisit this when better known
        hamSide[i] = (uint8_t)ham;
    }
    var = gid.getVar("HAM_side");
    var.putVar(start, count, hamSide);

    // Extract and write instrument spin ID
    //  int32_t *spinID = new int32_t[isc];

    for (size_t i = 0; i < isc; i++) {
        uint32_t ui32;
        memcpy(&ui32, &ancdata[i * ANCSIZE + toff], 4);
        spinID[i] = SWAP_4(ui32);
    }
    var = gid.getVar("spin_ID");
    var.putVar(start, count, spinID);

    // Packet and line number sequence error flags
    var = gid.getVar("pseq_flag");
    var.putVar(start, count, seqerr);
    var = gid.getVar("line_flag");
    var.putVar(start, count, linerr);

    delete[] stimes;
    delete[] scss;
    delete[] scsu;
    delete[] hamSide;
    // delete[] spinID;

    return 0;
}

int l1aFile::write_oci_ancil_data(uint32_t isc, uint8_t *ancdata) {
    vector<size_t> start;
    vector<size_t> count;
    start.push_back(0);

    uint16_t ui16;

    // Extract ancillary telemetry status flags and error counts
    short ioff = 12;
    int16_t *anctlm = new int16_t[6 * isc];
    for (size_t i = 0; i < 6 * isc; i++)
        anctlm[i] = BAD_INT;
    for (size_t i = 0; i < isc; i++) {
        for (size_t j = 0; j < 6; j++) {
            memcpy(&ui16, &ancdata[i * ANCSIZE + ioff + j * 2], 2);
            anctlm[i * 6 + j] = SWAP_2(ui16);
        }
    }

    // Extract spatial aggregation / data type table
    short dtype[10];
    short iagg[10];
    short lines[10];
    short nagg[4] = {1, 2, 4, 8};
    ioff = 36;

    int32_t iyear, iday;
    double stime;
    get_anc_packet_time(&ancdata[0], iyear, iday, stime);
    uint32_t jd = jday(iyear, 1, iday);

    for (size_t i = 0; i < 10; i++) {
        dtype[i] = ancdata[ioff + 3] % 16;
        if ((jd < 2459000) && (dtype[i] == 11))
            dtype[i] = 9;  // data type mod for ETU before June 2020
        iagg[i] = ancdata[ioff + 2] % 4;
        if (dtype[i] != 0)
            iagg[i] = nagg[iagg[i]];
        lines[i] = ancdata[ioff] * 256 + ancdata[ioff + 1];
        ioff += 4;
    }
    ioff += 4;

    // Extract spectral aggregation and compute numbers of bands

    // Tap enable flags
    uint16_t btap = ancdata[ioff + 2] * 256 + ancdata[ioff + 3];
    uint16_t rtap = ancdata[ioff + 0] * 256 + ancdata[ioff + 1];

    // Tap aggregation factors
    uint32_t ui32;
    int32_t i32;

    memcpy(&ui32, &ancdata[ioff + 8], 4);
    uint32_t bagg = SWAP_4(ui32);

    memcpy(&ui32, &ancdata[ioff + 4], 4);
    uint32_t ragg = SWAP_4(ui32);

    int16_t baggs[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int16_t raggs[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    // Compute number of bands for enabled taps
    uint16_t ken = 1;
    uint32_t kag = 3;
    uint32_t lag = 1;

    for (size_t i = 0; i < 16; i++) {  // Put tap information in ascending spectral order
        uint16_t btaps = (btap & ken) / ken;
        if (btaps)
            baggs[15 - i] = nagg[(bagg & kag) / lag];
        uint16_t rtaps = (rtap & ken) / ken;
        if (rtaps)
            raggs[15 - i] = nagg[(ragg & kag) / lag];

        ken *= 2;
        kag *= 4;
        lag *= 4;
    }

    NcGroup gid;
    NcVar var;

    // Write to file

    gid = l1afile->getGroup("spatial_spectral_modes");

    count.push_back(10);

    // Data type
    var = gid.getVar("spatial_zone_data_type");
    var.putVar(start, count, dtype);

    // Spatial aggregation
    var = gid.getVar("spatial_aggregation");
    var.putVar(start, count, iagg);

    // Number of lines
    var = gid.getVar("spatial_zone_lines");
    var.putVar(start, count, lines);

    count.pop_back();
    count.push_back(16);

    // Blue spectral aggregation
    var = gid.getVar("blue_spectral_mode");
    var.putVar(start, count, baggs);

    // Red spectral aggregation
    var = gid.getVar("red_spectral_mode");
    var.putVar(start, count, raggs);

    // Extract MCE data and write to file
    // Use ancillary packets AFTER science data

    gid = l1afile->getGroup("engineering_data");

    count.pop_back();
    count.push_back(isc);

    // write ancillary telemetry status flags and error counts
    start.push_back(0);
    count.push_back(6);
    var = gid.getVar("ancillary_tlm");
    var.putVar(start, count, anctlm);
    delete[] anctlm;
    start.pop_back();
    count.pop_back();

    // Aggregation control
    int32_t *aggcon = new int32_t[isc];
    for (size_t i = 0; i < isc; i++) {
        memcpy(&i32, &ancdata[i * ANCSIZE + ioff - 4], 4);
        aggcon[i] = SWAP_4(i32);
    }
    var = gid.getVar("agg_control");
    var.putVar(start, count, aggcon);
    delete[] aggcon;

    // Aggregation errors
    uint16_t *bagerr = new uint16_t[isc];
    uint16_t *ragerr = new uint16_t[isc];
    for (size_t i = 0; i < isc; i++) {
        memcpy(&ui16, &ancdata[(i + 1) * ANCSIZE + ioff + 12], 2);
        bagerr[i] = SWAP_2(ui16);

        memcpy(&ui16, &ancdata[(i + 1) * ANCSIZE + ioff + 14], 2);
        ragerr[i] = SWAP_2(ui16);
    }
    var = gid.getVar("blue_agg_error");
    var.putVar(start, count, bagerr);
    var = gid.getVar("red_agg_error");
    var.putVar(start, count, ragerr);
    delete[] bagerr;
    delete[] ragerr;

    // Digital card error status
    int32_t *digerr = new int32_t[isc];
    for (size_t i = 0; i < isc; i++) {
        memcpy(&i32, &ancdata[(i + 1) * ANCSIZE + ioff + 20], 4);
        digerr[i] = SWAP_4(i32);
    }
    var = gid.getVar("dig_card_error");
    var.putVar(start, count, digerr);
    delete[] digerr;

    start.push_back(0);
    count.push_back(9);

    return 0;
}

int l1aFile::write_oci_tlm_data(itab *itable, uint32_t ntlm, uint8_t (*tlmdata)[TLMSIZE], int32_t *spinID,
                                uint16_t &cdsmode, uint32_t isc, const time_struct &starttime) {
    // ntlm:  Number of telemetry packets
    // tlmdata: OCI telemetry data packets
    // cdsmode: CDS mode (0 = enabled, 1 = reset, 2 = video)
    // tlmzs(2): Telmetry zone start times (msec)
    // tlmzd(2): Telmetry zone durations (msec)
    // tditime: CCD data line TDI time (clock cycles)

    // Set up output array and extract/convert data from DAU packets
    const uint16_t dauapid = 723;
    const uint16_t ddcapid = 701;
    const uint16_t mceapid = 713;  // 711->713, changed by Liang Hong, 10/29/2020
    const uint16_t encapid = 712;
    const uint16_t scaapid = 717;
    const uint16_t senapid = 716;
    const uint16_t tcapid = 656;
    const uint16_t ddctapid = 703;
    const uint16_t dauctapid = 744;
    const uint16_t icdumcetapid = 745;

    // uint16_t nfptmps = 14;
    // uint16_t nsdtmps = 16;
    // uint16_t nabtmps = 9;
    // uint16_t ndtmps = 14;
    uint16_t nicthrm = 74;
    uint16_t ndctmps = 69;
    uint16_t nimtmps = 16;
    uint16_t nadclat = 4;
    uint16_t ndau = 0;
    uint16_t nddc = 0;
    uint16_t nmce = 0;
    uint16_t nenc = 0;
    uint16_t nsca = 0;
    uint16_t nsen = 0;
    uint16_t ntc = 0;
    uint16_t ndct = 0;
    uint16_t nddt = 0;
    uint16_t nimt = 0;
    int16_t *tlmzs = new int16_t[2];
    int16_t *tlmzd = new int16_t[2];
    uint16_t tditime = 0;

    double l1afile_epoch = yds2unix(starttime.iyear, starttime.iday, 0.0);
    string l1afile_epoch_str = unix2isodate(l1afile_epoch, 'G');
    l1afile_epoch_str = l1afile_epoch_str.substr(0, 10);
    l1afile_epoch_str.insert(0, "seconds since ");

    double *dausec = new double[ntlm];
    for (size_t i = 0; i < ntlm; i++)
        dausec[i] = BAD_FLT;
    int32_t *dauspin = new int32_t[ntlm];
    for (size_t i = 0; i < ntlm; i++)
        dauspin[i] = BAD_INT;
    uint8_t *dautlm = new uint8_t[620 * ntlm];
    for (size_t i = 0; i < 620 * ntlm; i++)
        dautlm[i] = 0;
    double *ddcsec = new double[ntlm];
    for (size_t i = 0; i < ntlm; i++)
        ddcsec[i] = BAD_FLT;
    uint8_t *ddctlm = new uint8_t[524 * ntlm];
    for (size_t i = 0; i < 524 * ntlm; i++)
        ddctlm[i] = 0;
    int32_t *mcespin = new int32_t[ntlm];
    for (size_t i = 0; i < ntlm; i++)
        mcespin[i] = BAD_INT;
    uint8_t *mcetlm = new uint8_t[480 * ntlm];
    for (size_t i = 0; i < 480 * ntlm; i++)
        mcetlm[i] = 0;
    int32_t *encspin = new int32_t[ntlm];
    for (size_t i = 0; i < ntlm; i++)
        encspin[i] = BAD_INT;
    int16_t *encoder = new int16_t[4 * 200 * ntlm];
    for (size_t i = 0; i < 4 * 200 * ntlm; i++)
        encoder[i] = BAD_INT;
    int16_t *auxparms = new int16_t[33];
    for (size_t i = 0; i < 33; i++)
        auxparms[i] = BAD_INT;
    int32_t *scaspin = new int32_t[ntlm];
    for (size_t i = 0; i < ntlm; i++)
        scaspin[i] = BAD_INT;
    uint8_t *scatlm = new uint8_t[480 * ntlm];
    for (size_t i = 0; i < 480 * ntlm; i++)
        scatlm[i] = 0;
    double *scasec = new double[ntlm];
    for (size_t i = 0; i < ntlm; i++)
        scasec[i] = BAD_FLT;
    float *scapos = new float[ntlm];
    for (size_t i = 0; i < ntlm; i++)
        scapos[i] = BAD_FLT;
    int32_t *senspin = new int32_t[ntlm];
    for (size_t i = 0; i < ntlm; i++)
        senspin[i] = BAD_INT;
    int16_t *sencoder = new int16_t[4 * 200 * ntlm];
    for (size_t i = 0; i < 4 * 200 * ntlm; i++)
        sencoder[i] = BAD_INT;
    double *tcsec = new double[ntlm];
    for (size_t i = 0; i < ntlm; i++)
        tcsec[i] = BAD_FLT;
    uint8_t *tctlm = new uint8_t[1216 * ntlm];
    for (size_t i = 0; i < 1216 * ntlm; i++)
        tctlm[i] = 0;
    double *dauctsec = new double[ntlm];
    for (size_t i = 0; i < ntlm; i++)
        dauctsec[i] = BAD_FLT;
    double *icdumcetsec = new double[ntlm];
    for (size_t i = 0; i < ntlm; i++)
        icdumcetsec[i] = BAD_FLT;
    uint8_t *icdumcettlm = new uint8_t[76 * ntlm];
    for (size_t i = 0; i < 76 * ntlm; i++)
        icdumcettlm[i] = 0;
    // uint16_t *redtemp = new uint16_t[nfptmps*ntlm];
    // for (size_t i=0;i<nfptmps*ntlm;i++) redtemp[i]=0;
    // uint16_t *blutemp = new uint16_t[nfptmps*ntlm];
    // for (size_t i=0;i<nfptmps*ntlm;i++) blutemp[i]=0;
    // uint16_t *sdstemp = new uint16_t[nsdtmps*ntlm];
    // for (size_t i=0;i<nsdtmps*ntlm;i++) sdstemp[i]=0;
    // uint16_t *aobtemp = new uint16_t[nabtmps*ntlm];
    // for (size_t i=0;i<nabtmps*ntlm;i++) aobtemp[i]=0;
    // uint16_t *dautemp = new uint16_t[ndtmps*ntlm];
    // for (size_t i=0;i<ndtmps*ntlm;i++) dautemp[i]=0;
    float *icdutherm = new float[nicthrm * ntlm];
    for (size_t i = 0; i < nicthrm * ntlm; i++)
        icdutherm[i] = BAD_FLT;
    float *dauctemp = new float[ndctmps * ntlm];
    for (size_t i = 0; i < ndctmps * ntlm; i++)
        dauctemp[i] = BAD_FLT;
    float *icdumcetemp = new float[nimtmps * ntlm];
    for (size_t i = 0; i < nimtmps * ntlm; i++)
        icdumcetemp[i] = BAD_FLT;
    uint8_t *adclat = new uint8_t[nadclat * ntlm];
    for (size_t i = 0; i < nadclat * ntlm; i++)
        adclat[i] = 0;
    uint8_t *cdsdis = new uint8_t[ntlm];
    for (size_t i = 0; i < ntlm; i++)
        cdsdis[i] = 0;
    uint8_t *hamside = new uint8_t[ntlm];
    uint32_t redmask[16];
    uint32_t bluemask[16];
    double sca_enc_scal = 360.0 / pow(2, 32);  // SCA diffuser position scale factor

    for (size_t i = 0; i < ntlm; i++) {
        uint32_t apid = ((uint8_t)tlmdata[i][0] % 8) * 256 + (uint8_t)tlmdata[i][1];
        int32_t iy, idy;
        double sc;
        uint8_t cctime[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        uint32_t ui32;
        // uint16_t ui16;
        int16_t i16;
        float f32;
        cctime[6] = 0;
        cctime[7] = 0;

        switch (apid) {
            case dauapid:
                memcpy(cctime, &tlmdata[i][6], 6);
                ccsds_sec_to_yds(cctime, &iy, &idy, &sc);
                dausec[ndau] = yds2unix(iy, idy, sc) - l1afile_epoch;
                memcpy(&ui32, &tlmdata[i][12], 4);
                dauspin[ndau] = SWAP_4(ui32);
                memcpy(&dautlm[ndau * 620], &tlmdata[i][16], 620);

                // Extract and convert temperatures
                /*
                for (size_t j=0; j<4; j++) {
                  memcpy( &ui16, &tlmdata[i][312+2*j], 2);
                  redtemp[ndau*nfptmps+j] = SWAP_2(ui16);

                  memcpy( &ui16, &tlmdata[i][376+2*j], 2);
                  blutemp[ndau*nfptmps+j] = SWAP_2(ui16);
                }
                for (size_t j=0; j<3; j++) {
                  memcpy( &ui16, &tlmdata[i][434+2*j], 2);
                  redtemp[ndau*nfptmps+4+j] = SWAP_2(ui16);

                  memcpy( &ui16, &tlmdata[i][466+2*j], 2);
                  blutemp[ndau*nfptmps+4+j] = SWAP_2(ui16);
                }
                for (size_t j=0; j<7; j++) {
                  memcpy( &ui16, &tlmdata[i][482+2*j], 2);
                  redtemp[ndau*nfptmps+7+j] = SWAP_2(ui16);

                  memcpy( &ui16, &tlmdata[i][546+2*j], 2);
                  blutemp[ndau*nfptmps+7+j] = SWAP_2(ui16);
                }

                for (size_t j=0; j<nsdtmps; j++) {
                  memcpy( &ui16, &tlmdata[i][162+2*j], 2);
                  sdstemp[ndau*nsdtmps+j] = SWAP_2(ui16);
                }

                for (size_t j=0; j<nabtmps; j++) {
                  memcpy( &ui16, &tlmdata[i][194+2*j], 2);
                  aobtemp[ndau*nabtmps+j] = SWAP_2(ui16);
                }

                for (size_t j=0; j<4; j++) {
                  memcpy( &ui16, &tlmdata[i][154+2*j], 2);
                  dautemp[ndau*ndtmps+j] = SWAP_2(ui16);
                }
                for (size_t j=0; j<4; j++) {
                  memcpy( &ui16, &tlmdata[i][212+2*j], 2);
                  dautemp[ndau*ndtmps+4+j] = SWAP_2(ui16);
                }
                for (size_t j=0; j<3; j++) {
                  memcpy( &ui16, &tlmdata[i][402+2*j], 2);
                  dautemp[ndau*ndtmps+8+j] = SWAP_2(ui16);
                }
                for (size_t j=0; j<3; j++) {
                  memcpy( &ui16, &tlmdata[i][608+2*j], 2);
                  dautemp[ndau*ndtmps+11+j] = SWAP_2(ui16);
                }
                */

                ndau++;
                break;

            case ddcapid:
                // DDC telemetry
                memcpy(cctime, &tlmdata[i][6], 6);
                ccsds_sec_to_yds(cctime, &iy, &idy, &sc);
                ddcsec[nddc] = yds2unix(iy, idy, sc) - l1afile_epoch;
                memcpy(&ddctlm[nddc * 524], &tlmdata[i][12], 524);
                cdsdis[nddc] = tlmdata[i][29];
                memcpy(&adclat[nddc * nadclat], &tlmdata[i][176], nadclat);
                if (nddc == 0) {
                    for (size_t j = 0; j < 16; j++) {
                        memcpy(&ui32, &tlmdata[i][196 + j * 4], 4);
                        redmask[j] = (uint32_t)SWAP_4(ui32);

                        memcpy(&ui32, &tlmdata[i][260 + j * 4], 4);
                        bluemask[j] = (uint32_t)SWAP_4(ui32);
                    }
                }
                nddc++;
                break;

            case mceapid:
                // RTA/HAM MCE telemetry
                memcpy(&ui32, &tlmdata[i][12], 4);
                mcespin[nmce] = SWAP_4(ui32);
                memcpy(&mcetlm[nmce * 480], &tlmdata[i][16], 480);
                hamside[nmce] = (tlmdata[i][49] & 8) / 8;
                nmce++;
                break;

            case encapid:
                // RTA/HAM MCE encoder data
                memcpy(&ui32, &tlmdata[i][12], 4);
                encspin[nenc] = SWAP_4(ui32);

                for (size_t j = 0; j < 4 * 200; j++) {
                    memcpy(&i16, &tlmdata[i][16 + 2 * j], 2);
                    encoder[nenc * 4 * 200 + j] = SWAP_2(i16);
                }
                nenc++;
                break;

            case scaapid:
                // SCA MCE telemetry
                memcpy(cctime, &tlmdata[i][6], 6);
                ccsds_sec_to_yds(cctime, &iy, &idy, &sc);
                scasec[nsca] = yds2unix(iy, idy, sc) - l1afile_epoch;
                memcpy(&ui32, &tlmdata[i][12], 4);
                scaspin[nsca] = SWAP_4(ui32);
                memcpy(&scatlm[nsca * 480], &tlmdata[i][16], 480);
                memcpy(&ui32, &tlmdata[i][348], 4);
                float sca_abs_enc_count;
                sca_abs_enc_count =
                    (float)SWAP_4(ui32);  // Convert byte array to unsigned long and little-endian
                scapos[nsca] = 330.0 - sca_abs_enc_count * sca_enc_scal;  // Perform conversion from Capon
                nsca++;
                break;

            case senapid:
                // SCA MCE encoder data
                memcpy(&ui32, &tlmdata[i][12], 4);
                senspin[nsen] = SWAP_4(ui32);
                for (size_t j = 0; j < 4 * 200; j++) {
                    memcpy(&i16, &tlmdata[i][16 + 2 * j], 2);
                    sencoder[nsen * 4 * 200 + j] = SWAP_2(i16);
                }
                nsen++;
                break;

            case tcapid:
                // ICDU TC telemetry and thermistors
                memcpy(cctime, &tlmdata[i][6], 6);
                ccsds_sec_to_yds(cctime, &iy, &idy, &sc);
                tcsec[ntc] = yds2unix(iy, idy, sc) - l1afile_epoch;
                memcpy(&tctlm[ntc * 1216], &tlmdata[i][12], 1216);
                for (size_t j = 0; j < nicthrm; j++) {
                    memcpy(&f32, &tlmdata[i][184 + 4 * j], 4);
                    icdutherm[ntc * nicthrm + j] = ReverseFloat(f32);
                }
                ntc++;
                break;

            case ddctapid:
                // DDC table telemetry
                // Only need one set of these
                if (nddt == 0) {
                    for (size_t j = 0; j < 33; j++) {
                        memcpy(&ui32, &tlmdata[i][64 + j * 4], 4);
                        auxparms[j] = (int16_t)SWAP_4(ui32);
                    }
                    nddt = 1;
                }
                break;

            case dauctapid:
                // FSW DAUC temperatures
                memcpy(cctime, &tlmdata[i][6], 6);
                ccsds_sec_to_yds(cctime, &iy, &idy, &sc);
                dauctsec[ndct] = yds2unix(iy, idy, sc) - l1afile_epoch;
                for (size_t j = 0; j < ndctmps; j++) {
                    memcpy(&f32, &tlmdata[i][12 + 4 * j], 4);
                    dauctemp[ndct * ndctmps + j] = ReverseFloat(f32);
                }
                ndct++;
                break;

            case icdumcetapid:
                // FSW ICCU/MCE telemetry and temperatures
                memcpy(cctime, &tlmdata[i][6], 6);
                ccsds_sec_to_yds(cctime, &iy, &idy, &sc);
                icdumcetsec[nimt] = yds2unix(iy, idy, sc) - l1afile_epoch;
                memcpy(&icdumcettlm[nimt * 76], &tlmdata[i][12], 76);
                for (size_t j = 0; j < nimtmps; j++) {
                    memcpy(&f32, &tlmdata[i][24 + 4 * j], 4);
                    icdumcetemp[nimt * nimtmps + j] = ReverseFloat(f32);
                }
                nimt++;
                break;
        }
    }  // tlm loop

    NcGroup gid;
    NcVar var;

    vector<size_t> start;
    vector<size_t> count;

    // Write to file

    gid = l1afile->getGroup("engineering_data");

    if (ndau == 0)
        ndau = 1;  // 0.99.20

    // if (ndau > 0) {
    start.push_back(0);
    count.push_back(ndau);

    var = gid.getVar("DAU_tlm_time");
    var.putVar(start, count, dausec);
    var.putAtt("units", l1afile_epoch_str);

    var = gid.getVar("DAU_spin_ID");
    var.putVar(start, count, dauspin);

    start.push_back(0);
    count.push_back(620);

    var = gid.getVar("DAU_telemetry");
    var.putVar(start, count, dautlm);

    // count.pop_back();
    // count.push_back(nfptmps);

    // var = gid.getVar("blue_FPA_temperatures");
    // var.putVar(start, count, blutemp);
    // var = gid.getVar("red_FPA_temperatures");
    // var.putVar(start, count, redtemp);

    // count.pop_back();
    // count.push_back(nsdtmps);

    // var = gid.getVar("SDS_det_temperatures");
    // var.putVar(start, count, sdstemp);

    // count.pop_back();
    // count.push_back(nabtmps);

    // var = gid.getVar("AOB_temperatures");
    // var.putVar(start, count, aobtemp);

    // count.pop_back();
    // count.push_back(ndtmps);

    // var = gid.getVar("DAU_temperatures");
    // var.putVar(start, count, dautemp);

    // Get telemetry zone start times and durations
    tlmzs[0] = dautlm[123];
    tlmzs[1] = dautlm[125];
    tlmzd[0] = dautlm[122];
    tlmzd[1] = dautlm[124];
    //}

    if (nddc == 0) {
        for (size_t i = 0; i < 16; i++) {
            redmask[i] = 4294967295;
            bluemask[i] = 4294967295;
        }
        nddc = 1;  // 0.99.20
    }
    // if (nddc > 0) {
    start.clear();
    count.clear();
    start.push_back(0);
    count.push_back(nddc);

    // DDC tlm time
    var = gid.getVar("DDC_tlm_time");
    var.putVar(start, count, ddcsec);
    var.putAtt("units", l1afile_epoch_str);

    // CDS disable flag
    var = gid.getVar("CDS_disable");
    var.putVar(start, count, cdsdis);

    start.push_back(0);
    count.push_back(524);

    // DDC telemetry
    var = gid.getVar("DDC_telemetry");
    var.putVar(start, count, ddctlm);

    count.pop_back();
    count.push_back(nadclat);

    // ADC latency
    var = gid.getVar("ADC_latency");
    var.putVar(start, count, adclat);

    // Red and blue channel masks
    start.clear();
    count.clear();
    start.push_back(0);
    count.push_back(16);
    var = gid.getVar("blue_channel_mask");
    var.putVar(start, count, bluemask);
    var = gid.getVar("red_channel_mask");
    var.putVar(start, count, redmask);

    // Get CDS mode for metadata
    cdsmode = cdsdis[0] * (adclat[0] - 14);

    //  Get TDI time
    uint16_t ui16;
    memcpy(&ui16, &ddctlm[346], 2);
    tditime = SWAP_2(ui16);
    //}

    if (nmce == 0)
        nmce = 1;  // 0.99.20
                   // if (nmce > 0) {
    start.clear();
    count.clear();
    start.push_back(0);
    count.push_back(nmce);

    // RTA/HAM MCE spin ID
    var = gid.getVar("MCE_spin_ID");
    var.putVar(start, count, mcespin);

    // HAM side
    // if (nmce > 1) {
    uint8_t *mside = new uint8_t[isc];
    for (size_t i = 0; i < isc; i++)
        mside[i] = 255;
    for (size_t i = 0; i < isc; i++) {
        for (size_t j = 0; j < ntlm; j++) {
            if ((int)mcespin[j] == spinID[i]) {
                mside[i] = hamside[j];
                break;
            } else {
                mside[i] = 1 - mside [i-1];
            }
        }
    }
    count.pop_back();
    count.push_back(isc);  // LH, 11/20/2020
    NcGroup scgid = l1afile->getGroup("scan_line_attributes");

    var = scgid.getVar("HAM_side");
    //      var.putVar(start, count, &hamside[1]);
    var.putVar(start, count, mside);
    delete[] mside;
    //}

    start.push_back(0);
    count.clear();          // LH, 11/20/2020
    count.push_back(nmce);  // 1.00.00
    count.push_back(480);

    // RTA/HAM MCE telemetry
    var = gid.getVar("MCE_telemetry");
    var.putVar(start, count, mcetlm);
    //}

    if (nenc == 0)
        nenc = 1;  // 0.99.20
                   // if (nenc > 0) {
    start.clear();
    count.clear();
    start.push_back(0);
    count.push_back(nenc);

    // RTA/HAM encoder spin ID
    var = gid.getVar("encoder_spin_ID");
    var.putVar(start, count, encspin);

    start.push_back(0);
    count.push_back(200);
    start.push_back(0);
    count.push_back(4);

    // RTA/HAM encoder data
    var = gid.getVar("MCE_encoder_data");
    var.putVar(start, count, encoder);
    //}

    if (nsca == 0)
        nsca = 1;  // 1.01.00
    start.clear();
    count.clear();
    start.push_back(0);
    count.push_back(nsca);
    // SCA MCE spin ID
    var = gid.getVar("SCA_spin_ID");
    var.putVar(start, count, scaspin);

    // SCA MCE telemetry
    start.clear();
    start.push_back(0);
    count.clear();
    count.push_back(nsca);
    count.push_back(480);

    var = gid.getVar("SCA_telemetry");
    var.putVar(start, count, scatlm);

    // SCA telemetry time
    start.clear();
    count.clear();
    start.push_back(0);
    count.push_back(nsca);

    var = gid.getVar("SCA_tlm_time");
    var.putVar(start, count, scasec);
    var.putAtt("units", l1afile_epoch_str);

    // SCA diffuser position
    var = gid.getVar("SCA_diffuser_position");
    var.putVar(start, count, scapos);

    if (nsen == 0)
        nsen = 1;
    start.clear();
    count.clear();
    start.push_back(0);
    count.push_back(nsen);

    // RTA/HAM encoder spin ID
    var = gid.getVar("SCA_encoder_spin_ID");
    var.putVar(start, count, senspin);

    // SCA encoder data
    start.push_back(0);
    count.push_back(200);
    start.push_back(0);
    count.push_back(4);

    // RTA/HAM encoder data
    var = gid.getVar("SCA_encoder_data");
    var.putVar(start, count, sencoder);

    if (ntc == 0)
        ntc = 1;  // 0.99.20
                  // if (ntc > 0) {
    start.clear();
    count.clear();
    start.push_back(0);
    count.push_back(ntc);

    // TC tlm time
    var = gid.getVar("TC_tlm_time");
    var.putVar(start, count, tcsec);
    var.putAtt("units", l1afile_epoch_str);

    start.push_back(0);
    count.push_back(1216);

    // TC telemetry
    var = gid.getVar("TC_telemetry");
    var.putVar(start, count, tctlm);

    count.pop_back();
    count.push_back(nicthrm);

    // ICDU thermisters
    var = gid.getVar("ICDU_thermisters");
    var.putVar(start, count, icdutherm);
    //}

    // Auxilary parameter table
    gid = l1afile->getGroup("spatial_spectral_modes");

    // if (nddt > 0) {
    start.clear();
    count.clear();
    start.push_back(0);
    count.push_back(33);

    var = gid.getVar("aux_param_table");
    var.putVar(start, count, auxparms);
    //}

    if (ndct == 0)
        ndct = 1;
    start.clear();
    count.clear();
    start.push_back(0);
    count.push_back(ndct);

    gid = l1afile->getGroup("engineering_data");

    // FSW DAUC temperature time
    var = gid.getVar("DAUC_temp_time");
    var.putVar(start, count, dauctsec);
    var.putAtt("units", l1afile_epoch_str);

    start.push_back(0);
    count.push_back(ndctmps);

    // FSW DAUC temperatures
    var = gid.getVar("DAUC_temperatures");
    var.putVar(start, count, dauctemp);

    if (nimt == 0)
        nimt = 1;
    start.clear();
    count.clear();
    start.push_back(0);
    count.push_back(nimt);

    // FSW ICDU/MCE temperature time
    var = gid.getVar("ICDU_MCE_temp_time");
    var.putVar(start, count, icdumcetsec);
    var.putAtt("units", l1afile_epoch_str);

    start.push_back(0);
    count.push_back(76);

    // FSW ICDU/MCE temperature telemetry
    var = gid.getVar("ICDU_MCE_temp_tlm");
    var.putVar(start, count, icdumcettlm);

    count.pop_back();
    count.push_back(nimtmps);

    // FSW ICDU/MCE temperatures
    var = gid.getVar("ICDU_MCE_temperatures");
    var.putVar(start, count, icdumcetemp);

    // check the science data and telemetry collection zones for conflicts and write the results
    // as an attribute to the engineering data group.  There are two telemetry collection zones
    string conflict = "No";
    gid = l1afile->getGroup("engineering_data");

    if ((tlmzd[0] == 0) && (tlmzd[1] == 0)) {
        cout << "No telemetry zone fields in file" << endl;
        gid.putAtt("science_telemetry_zone_conflict", conflict);
    } else {
        float clock = 136000;                    // Master clock frequency in msec
        float secpline = (tditime + 1) / clock;  // msec per line

        // Find no-data zones
        float *znd = new float[2 * 10];
        float *zndcp = new float[2 * 10];
        int ndz = 0;
        uint16_t line = 0;
        for (size_t i = 0; i < 10; i++) {
            if (((itable[i].dtype == 0) || (itable[i].dtype == 10)) && (itable[i].lines > 0)) {
                znd[ndz] = line * secpline;
                znd[1 * 10 + ndz] = (line + itable[i].lines) * secpline;
                // Check for consecutive no-data zones
                if ((ndz > 0) && (znd[ndz] == znd[1 * 10 + ndz - 1])) {
                    znd[1 * 10 + ndz - 1] = znd[1 * 10 + ndz];
                } else {
                    ndz++;
                }
            }
            line += itable[i].lines;
        }
        memcpy(zndcp, znd, 20);
        // znd new dimension: 2*ndz
        for (int i = 0; i < ndz; i++)
            znd[ndz + i] = znd[10 + i];  // znd = znd(*,0:ndz-1)

        //  If first no-data zone is at start of spin, extend last zone by that amount.
        if ((itable[0].dtype == 0) || (itable[0].dtype == 10))
            znd[1 * ndz + ndz - 1] += znd[1 * ndz + 0];

        if (ndz > 0) {
            // Loop through telemetry zones
            for (size_t i = 0; i < 2; i++) {
                // Find overlapping no-data zone
                int16_t tlmze = tlmzs[i] + tlmzd[i] + 3;
                int k;
                for (k = ndz - 1; k >= 0; k--) {
                    // k = max(where(tlmze gt znd(0,*)))
                    if (tlmze > znd[k])
                        break;
                }
                if ((k == -1) || (tlmzs[i] < znd[k]) || (tlmze > znd[1 * ndz + k])) {
                    conflict = "Yes";
                    cout << "Telemetry zone " << tlmzs[i] << ", " << tlmze << endl;
                    cout << "No-data zone " << znd[k] << ", " << znd[1 * ndz + k] << endl;
                }
            }
        } else
            conflict = "Yes";

        //  Write attribute to engineering_data_group
        //   Open group
        NcGroup gid = l1afile->getGroup("engineering_data");
        gid.putAtt("science_telemetry_zone_conflict", conflict);

        delete[] znd;
        delete[] zndcp;
    }  // if (tlmzd[0] == 0) && (tlmzd[1] == 0) else

    delete[] dausec;
    delete[] dauspin;
    delete[] dautlm;
    delete[] ddcsec;
    delete[] ddctlm;
    delete[] mcespin;
    delete[] mcetlm;
    delete[] encspin;
    delete[] encoder;
    delete[] auxparms;
    delete[] tcsec;
    delete[] tctlm;
    delete[] dauctsec;
    delete[] icdumcetsec;
    delete[] icdumcettlm;
    delete[] icdumcetemp;
    delete[] dauctemp;
    delete[] scatlm;
    delete[] sencoder;
    // delete[] redtemp;
    // delete[] blutemp;
    // delete[] sdstemp;
    // delete[] aobtemp;
    // delete[] dautemp;
    delete[] icdutherm;
    delete[] adclat;
    delete[] cdsdis;
    delete[] hamside;
    delete[] tlmzs;
    delete[] tlmzd;
    delete[] scaspin;
    delete[] scapos;
    delete[] senspin;
    delete[] scasec;

    return 0;
}

int find_nav_index(char *hkt_t_start, double usec_l1a_start, double usec_l1a_end, double *nav_t, size_t n_nav,
                   int &ind_start, int &ind_end) {
    double usec_hkt;
    int hkt_yr, hkt_mon, hkt_day;
    sscanf(hkt_t_start, "%4d-%2d-%2d", &hkt_yr, &hkt_mon, &hkt_day);

    for (size_t i = 0; i < n_nav; i++) {
        usec_hkt = ymds2unix(hkt_yr, hkt_mon, hkt_day, nav_t[i]);
        if ((usec_hkt > usec_l1a_start - 10) && (ind_start == 1e6)) {
            ind_start = i;
            break;
        }
    }

    for (size_t i = n_nav - 1; i >= 0; i--) {
        usec_hkt = ymds2unix(hkt_yr, hkt_mon, hkt_day, nav_t[i]);
        if ((usec_hkt < usec_l1a_end + 10) && (ind_end == -1)) {
            ind_end = i;
            break;
        }
    }
    return 0;
}

/*
 * Specify all the L0 and hkt files used in the processing of this L1A file
 * L0 is given as a vector because
 * @return 1 on success and 0 on failure
 */
int l1aFile::write_processing_control(std::string hktList, std::string l0List, std::string time_start, std::string maxgap,  std::string nametag,
                                    std::string swir_loff_set, std::string outlist, std::string outfile, std::string doi, 
                                    std::string pversion, std::string isSPW) {
    try {
        // add the processing control group
        netCDF::NcGroup processCtrlGrp = l1afile->addGroup("processing_control");

        // within process control group, add input parameter group as its child
        netCDF::NcGroup inputParaGrp = processCtrlGrp.addGroup("input_parameter");

        // add current input parameter file names to the input parameter group.
        // this only contains the important command line args for l1agen_oci
        inputParaGrp.putAtt("OCI_packet_file", l0List);
        inputParaGrp.putAtt("maxgap", maxgap);
        inputParaGrp.putAtt("hktlist_iFile", hktList);
        inputParaGrp.putAtt("swir_loff_set", swir_loff_set);
        inputParaGrp.putAtt("start_time", time_start);
        inputParaGrp.putAtt("outlist", outlist);
        inputParaGrp.putAtt("outfile", outfile);
        inputParaGrp.putAtt("nametag", nametag);
        inputParaGrp.putAtt("doi", doi);
        inputParaGrp.putAtt("pversion", pversion);
        inputParaGrp.putAtt("isSPW", isSPW);
        
        // add current software details
        processCtrlGrp.putAtt("software_name", "l1agen_oci");
        processCtrlGrp.putAtt("software_version", VERSION);

        // given hkt list path, convert it into a vector and add
        // all the hkt files used separated by commas (,)
        vector<std::string> hktFiles = readFileList(hktList);
        string hktString = "";
        for (size_t i = 0; i < hktFiles.size(); i++) {
            //
            if (i < hktFiles.size() - 1) {
                hktString.append(hktFiles[i]);
                hktString.append(", ");
                continue;
            }
            // dont end with a  comma if last file
            hktString.append(hktFiles[i]);
        }
        processCtrlGrp.putAtt("hkt_list", hktString);

        // do the same for l0 file list
        vector<std::string> l0Files = readFileList(l0List);
        string l0String = "";
        for (size_t i = 0; i < l0Files.size(); i++) {
            if (i < l0Files.size() - 1) {
                l0String.append(l0Files[i]);
                l0String.append(", ");
                continue;
            }
            l0String.append(l0Files[i]);
        }
        processCtrlGrp.putAtt("l0_list", l0String);

    } catch (std::exception &e) {
        std::cerr << "--Error-- writing processing control group. what: " << e.what() << endl;
        return 0;
    }

    return 1;
}

int l1aFile::write_navigation(std::string hktlist, time_struct &starttime, time_struct &endtime) {
    uint16_t nnavmax = 10000;
    size_t iatt = 0, iorb = 0, itilt = 0;

    double usec_l1a_start, usec_l1a_end;
    usec_l1a_start = yds2unix(starttime.iyear, starttime.iday, starttime.sec);
    usec_l1a_end = yds2unix(endtime.iyear, endtime.iday, endtime.sec);
    double l1afile_epoch = yds2unix(starttime.iyear, starttime.iday, 0.0);

    string l1afile_epoch_str = unix2isodate(l1afile_epoch, 'G');
    l1afile_epoch_str = l1afile_epoch_str.substr(0, 10);
    l1afile_epoch_str.insert(0, "seconds since ");

    // read HKT file list, loop through HKT files
    ifstream file(hktlist);
    string strHKTfile;

    double *atime = new double[nnavmax];
    for (size_t i = 0; i < nnavmax; i++)
        atime[i] = BAD_FLT;
    double *otime = new double[nnavmax];
    for (size_t i = 0; i < nnavmax; i++)
        otime[i] = BAD_FLT;
    double *ttime = new double[nnavmax];
    for (size_t i = 0; i < nnavmax; i++)
        ttime[i] = BAD_FLT;

    float *arate = new float[3 * nnavmax];
    for (size_t i = 0; i < 3 * nnavmax; i++)
        arate[i] = BAD_FLT;
    float *quat = new float[4 * nnavmax];
    for (size_t i = 0; i < 4 * nnavmax; i++)
        quat[i] = BAD_FLT;
    float *pos = new float[3 * nnavmax];
    for (size_t i = 0; i < 3 * nnavmax; i++)
        pos[i] = -9999999;
    float *vel = new float[3 * nnavmax];
    for (size_t i = 0; i < nnavmax; i++)
        vel[i] = -9999999;
    for (size_t i = nnavmax; i < 3 * nnavmax; i++)
        vel[i] = BAD_FLT;
    float *tlt = new float[nnavmax];
    for (size_t i = 0; i < nnavmax; i++)
        tlt[i] = BAD_FLT;
    uint8_t *tltflg = new uint8_t[nnavmax];
    for (size_t i = 0; i < nnavmax; i++)
        tltflg[i] = 0;

    NcVar var;

    char *hkt_t_start = new char[100];
    char *hkt_t_end = new char[100];

    while (getline(file, strHKTfile)) {
        try {
            // read HKT netcdf file
            int status, ncid, gid;
            int atttid, attqid, attrid;
            int orbpid, orbtid, orbvid;
            int tid, ttid, tfid;
            size_t natt = 0, norb = 0, ntilt = 0;

            nc_open(strHKTfile.c_str(), NC_NOWRITE, &ncid);

            // find HKT data year, month, date from global attributes
            nc_get_att_text(ncid, NC_GLOBAL, "time_coverage_start",
                            hkt_t_start);  // e.g. 2022-04-20T17:26:00.000000Z
            size_t length;
            nc_inq_attlen(ncid, NC_GLOBAL, "time_coverage_start", &length);
            hkt_t_start[length] = '\0';

            nc_get_att_text(ncid, NC_GLOBAL, "time_coverage_end", hkt_t_end);
            nc_inq_attlen(ncid, NC_GLOBAL, "time_coverage_end", &length);
            hkt_t_end[length] = '\0';

            // skip if PACE HKT data doesn't overlap OCI data
            int hkt_yr, hkt_mon, hkt_day, hkt_hr, hkt_min, hk_sec;
            double usec_hkt;
            sscanf(hkt_t_start, "%4d-%2d-%2dT%2d:%2d:%2d", &hkt_yr, &hkt_mon, &hkt_day, &hkt_hr, &hkt_min,
                   &hk_sec);
            usec_hkt = ymds2unix(hkt_yr, hkt_mon, hkt_day, hkt_hr * 3600 + hkt_min * 60 + hk_sec);
            double hktfile_epoch = ymds2unix(hkt_yr, hkt_mon, hkt_day, 0.0);
            if (usec_hkt > usec_l1a_end + 10)
                continue;
            sscanf(hkt_t_end, "%4d-%2d-%2dT%2d:%2d:%2d", &hkt_yr, &hkt_mon, &hkt_day, &hkt_hr, &hkt_min,
                   &hk_sec);
            usec_hkt = ymds2unix(hkt_yr, hkt_mon, hkt_day, hkt_hr * 3600 + hkt_min * 60 + hk_sec);
            if (usec_hkt < usec_l1a_start - 10)
                continue;

            // cout<<"Reading PACE HKT file "<<strHKTfile<<"..."<<endl;

            nc_inq_grp_ncid(ncid, "navigation_data", &gid);

            nc_type attt_type, orbt_type, tiltt_type;
            int attt_ndims, orbt_ndims, tiltt_ndims;
            int attt_dimids[NC_MAX_VAR_DIMS], orbt_dimids[NC_MAX_VAR_DIMS], tiltt_dimids[NC_MAX_VAR_DIMS];
            int attt_natts, orbt_natts, tiltt_natts;

            status = nc_inq_varid(gid, "att_time", &atttid);
            if (status == NC_NOERR) {
                nc_inq_var(gid, atttid, 0, &attt_type, &attt_ndims, attt_dimids, &attt_natts);
                nc_inq_dimlen(gid, attt_dimids[0], &natt);
            }

            status = nc_inq_varid(gid, "orb_time", &orbtid);
            if (status == NC_NOERR) {
                nc_inq_var(gid, orbtid, 0, &orbt_type, &orbt_ndims, orbt_dimids, &orbt_natts);
                nc_inq_dimlen(gid, orbt_dimids[0], &norb);
            }

            status = nc_inq_varid(gid, "tilt_time", &ttid);
            if (status == NC_NOERR) {
                nc_inq_var(gid, ttid, 0, &tiltt_type, &tiltt_ndims, tiltt_dimids, &tiltt_natts);
                nc_inq_dimlen(gid, tiltt_dimids[0], &ntilt);
            }

            // read data from HKT file
            // find values within L1A file time range, i.e. [starttime-10 sec, endtime+10 sec]
            int ind_att_start = 1e6, ind_att_end = -1;
            int ind_orb_start = 1e6, ind_orb_end = -1;
            int ind_tilt_start = 1e6, ind_tilt_end = -1;

            if (natt > 0) {
                double *at = new double[natt];
                nc_get_var(gid, atttid, at);
                float *r = new float[3 * natt];
                float *q = new float[4 * natt];
                find_nav_index(hkt_t_start, usec_l1a_start, usec_l1a_end, at, natt, ind_att_start,
                               ind_att_end);

                if ((ind_att_end - ind_att_start) > 0) {
                    nc_inq_varid(gid, "att_quat", &attqid);
                    nc_get_var(gid, attqid, q);
                    nc_inq_varid(gid, "att_rate", &attrid);
                    nc_get_var(gid, attrid, r);

                    // fix the time if l1a and hkt files have different time epoch
                    if (l1afile_epoch != hktfile_epoch) {
                        double epoch_diff = hktfile_epoch - l1afile_epoch;
                        for (size_t i = 0; i < natt; i++) {
                            at[i] += epoch_diff;
                        }
                    }

                    // concatenate arrays
                    size_t nnav = ind_att_end - ind_att_start + 1;
                    memcpy(atime + iatt, at + ind_att_start, nnav * sizeof(double));
                    memcpy(quat + iatt * 4, q + ind_att_start * 4, 4 * nnav * sizeof(float));
                    memcpy(arate + iatt * 3, r + ind_att_start * 3, 3 * nnav * sizeof(float));

                    iatt += nnav;
                }
                delete[] at;
                delete[] q;
                delete[] r;
            }

            if (norb > 0) {
                double *ot = new double[norb];
                nc_get_var(gid, orbtid, ot);
                float *p = new float[3 * norb];
                float *v = new float[3 * norb];
                find_nav_index(hkt_t_start, usec_l1a_start, usec_l1a_end, ot, norb, ind_orb_start,
                               ind_orb_end);

                if ((ind_orb_end - ind_orb_start) > 0) {
                    nc_inq_varid(gid, "orb_pos", &orbpid);
                    nc_get_var(gid, orbpid, p);
                    nc_inq_varid(gid, "orb_vel", &orbvid);
                    nc_get_var(gid, orbvid, v);

                    // fix the time if l1a and hkt files have different time epoch
                    if (l1afile_epoch != hktfile_epoch) {
                        double epoch_diff = hktfile_epoch - l1afile_epoch;
                        for (size_t i = 0; i < norb; i++) {
                            ot[i] += epoch_diff;
                        }
                    }

                    // concatenate arrays
                    size_t nnav = ind_orb_end - ind_orb_start + 1;
                    memcpy(otime + iorb, ot + ind_orb_start, nnav * sizeof(double));
                    memcpy(pos + iorb * 3, p + ind_orb_start * 3, 3 * nnav * sizeof(float));
                    memcpy(vel + iorb * 3, v + ind_orb_start * 3, 3 * nnav * sizeof(float));

                    iorb += nnav;
                }

                delete[] ot;
                delete[] p;
                delete[] v;
            }

            if (ntilt > 0) {
                double *tt = new double[ntilt];
                nc_get_var(gid, ttid, tt);
                float *t = new float[ntilt];
                uint8_t *tf = new uint8_t[ntilt];
                find_nav_index(hkt_t_start, usec_l1a_start, usec_l1a_end, tt, ntilt, ind_tilt_start,
                               ind_tilt_end);

                if ((ind_tilt_end - ind_tilt_start) > 0) {
                    nc_inq_varid(gid, "tilt", &tid);
                    nc_get_var(gid, tid, t);

                    // fix the time if l1a and hkt files have different time epoch
                    if (l1afile_epoch != hktfile_epoch) {
                        double epoch_diff = hktfile_epoch - l1afile_epoch;
                        for (size_t i = 0; i < ntilt; i++) {
                            tt[i] += epoch_diff;
                        }
                    }

                    // concatenate arrays
                    size_t nnav = ind_tilt_end - ind_tilt_start + 1;
                    memcpy(ttime + itilt, tt + ind_tilt_start, nnav * sizeof(double));
                    memcpy(tlt + itilt, t + ind_tilt_start, nnav * sizeof(float));

                    // if tilt flag is included in the HKT data
                    status = nc_inq_varid(gid, "tilt_flag", &tfid);
                    if (status == NC_NOERR) {
                        nc_get_var(gid, tfid, tf);
                        memcpy(tltflg + itilt, tf + ind_tilt_start, nnav);

                        // set tilt time and tile angle to fill value where the tilt flag is set
                        for (size_t itf = itilt; itf < itilt + nnav; itf++) {
                            if (tltflg[itf] > 0) {
                                ttime[itf] = BAD_FLT;
                                tlt[itf] = BAD_FLT;
                            }
                        }
                    }

                    itilt += nnav;
                }

                delete[] tt;
                delete[] t;
                delete[] tf;
            }

            nc_close(ncid);

        } catch (NcException &e) {
            // e.what();
            cout << "Error reading " << strHKTfile << "!" << endl;
            file.close();
            return NC2_ERR;
        }

    }  // while (getline(file, strHKTfile))
    file.close();

    // write to L1A file
    vector<size_t> start;
    vector<size_t> count;

    NcGroup l1a_gid = l1afile->getGroup("navigation_data");

    if (iatt == 0)
        iatt = 1;
    start.push_back(0);
    count.push_back(iatt);
    var = l1a_gid.getVar("att_time");
    var.putVar(start, count, atime);
    var.putAtt("units", l1afile_epoch_str);

    start.push_back(0);
    count.push_back(4);
    var = l1a_gid.getVar("att_quat");
    var.putVar(start, count, quat);

    count.pop_back();
    count.push_back(3);
    var = l1a_gid.getVar("att_rate");
    var.putVar(start, count, arate);
    start.clear();
    count.clear();
    start.push_back(0);
    count.push_back(iorb);
    var = l1a_gid.getVar("orb_time");
    var.putVar(start, count, otime);
    var.putAtt("units", l1afile_epoch_str);

    start.push_back(0);
    count.push_back(3);
    var = l1a_gid.getVar("orb_pos");
    var.putVar(start, count, pos);

    var = l1a_gid.getVar("orb_vel");
    var.putVar(start, count, vel);

    start.clear();
    count.clear();
    start.push_back(0);
    count.push_back(itilt);
    var = l1a_gid.getVar("tilt_time");
    var.putVar(start, count, ttime);
    var.putAtt("units", l1afile_epoch_str);

    var = l1a_gid.getVar("tilt");
    var.putVar(start, count, tlt);

    delete[] hkt_t_start;
    delete[] hkt_t_end;

    delete[] atime;
    delete[] otime;
    delete[] ttime;
    delete[] arate;
    delete[] quat;
    delete[] pos;
    delete[] vel;
    delete[] tlt;
    delete[] tltflg;

    return 0;
}

int l1aFile::write_oci_global_metadata(time_struct &starttime, time_struct &endtime, std::string l1a_name,
                                       std::string sdir, std::string edir, uint32_t isc, short dtype,
                                       uint16_t smode, uint16_t cdsmode, std::ofstream &fout) {
    string dtypes[] = {"",
                       "Earth Collect",
                       "Dark Collect",
                       "Solar Cal",
                       "SPCA Cal",
                       "Response Curve",
                       "Lunar Cal",
                       "Diagnostic",
                       "Static",
                       "Earth Spectral",
                       "",
                       "External Snapshot Trigger",
                       "Internal Snapshot Trigger",
                       "Earth Collect with Non-Baseline Spectral Bands"};

    string smodes[] = {"Science", "Diagnostic", "Single-image raw", "Test pattern"};

    string cdsmodes[] = {"CDS", "Reset", "Video"};

    // Write start, end, create time attributes
    int16_t mon, idm;
    int32_t ih, mn;
    stringstream ss;

    yd2md((int16_t)starttime.iyear, (int16_t)starttime.iday, &mon, &idm);

    starttime.sec = floor(starttime.sec * 1000) / 1000;  // LH, 11/18/2020
    ih = (int)(starttime.sec / 3600);
    starttime.sec -= ih * 3600;
    mn = (int)(starttime.sec / 60);
    starttime.sec -= mn * 60;

    // yyyy-mn-dyThr:mn:ss.sss
    ss = stringstream();
    ss << setw(4) << to_string(starttime.iyear) << "-";
    ss << setw(2) << setfill('0') << mon << "-";
    ss << setw(2) << setfill('0') << idm << "T";

    ss << setw(2) << setfill('0') << ih << ":";
    ss << setw(2) << setfill('0') << mn << ":";
    ss << fixed << setw(6) << setprecision(3) << setfill('0') << starttime.sec;

    // cout << ss.str() << endl;
    l1afile->putAtt("time_coverage_start", ss.str() + "Z");

    // Write to outlist if needed
    if (fout.is_open())
        fout << ss.str().c_str() << " ";

    yd2md((int16_t)endtime.iyear, (int16_t)endtime.iday, &mon, &idm);

    endtime.sec =
        floor(endtime.sec * 1000) / 1000;  // LH, 11/18/2020, to handle 2020-11-05T20:06:59.999980000
    ih = (int)(endtime.sec / 3600);
    endtime.sec -= ih * 3600;
    mn = (int)(endtime.sec / 60);
    endtime.sec -= mn * 60;

    // yyyy-mn-dyThr:mn:ss.sss
    ss = stringstream();
    ss << setw(4) << to_string(endtime.iyear) << "-";
    ss << setw(2) << setfill('0') << mon << "-";
    ss << setw(2) << setfill('0') << idm << "T";

    ss << setw(2) << setfill('0') << ih << ":";
    ss << setw(2) << setfill('0') << mn << ":";
    ss << fixed << setw(6) << setprecision(3) << setfill('0') << endtime.sec;

    // cout << ss.str() << endl;
    l1afile->putAtt("time_coverage_end", ss.str() + "Z");

    // Write product file name
    l1afile->putAtt("product_name", l1a_name.c_str());

    // Write orbit direction
    l1afile->putAtt("startDirection", sdir.c_str());
    l1afile->putAtt("endDirection", edir.c_str());

    // Write data collect type and SWIR mode
    l1afile->putAtt("data_collect_mode", dtypes[dtype].c_str());
    l1afile->putAtt("SWIR_data_mode", smodes[smode].c_str());
    l1afile->putAtt("CDS_mode", cdsmodes[cdsmode].c_str());

    // Write to outlist if needed
    if (fout.is_open())
        fout << ss.str().c_str() << " ";

    return 0;
}

int l1aFile::close() {
    try {
        l1afile->close();
    } catch (NcException &e) {
        cout << e.what() << endl;
        cerr << "Failure closing: " + fileName << endl;
        exit(1);
    }

    return 0;
}

int createNCDF(NcGroup &ncGrp, const char *sname, const char *lname, const char *standard_name,
               const char *units, void *fill_value, const char *flag_values, const char *flag_meanings,
               const char *reference, double low, double high, int nt, vector<NcDim> &varVec) {
    /* Create the NCDF dataset */
    NcVar ncVar;
    try {
        ncVar = ncGrp.addVar(sname, nt, varVec);
    } catch (NcException &e) {
        cout << e.what() << endl;
        cerr << "Failure creating variable: " << sname << endl;
        exit(1);
    }

    // Set fill value
    double fill_value_dbl;
    memcpy(&fill_value_dbl, fill_value, sizeof(double));

    int8_t i8;
    uint8_t ui8;
    int16_t i16;
    uint16_t ui16;
    int32_t i32;
    uint32_t ui32;
    float f32;

    if (low != fill_value_dbl) {
        if (nt == NC_BYTE) {
            i8 = fill_value_dbl;
            ncVar.setFill(true, (void *)&i8);
        } else if (nt == NC_UBYTE) {
            ui8 = fill_value_dbl;
            ncVar.setFill(true, (void *)&ui8);
        } else if (nt == NC_SHORT) {
            i16 = fill_value_dbl;
            ncVar.setFill(true, (void *)&i16);
        } else if (nt == NC_USHORT) {
            ui16 = fill_value_dbl;
            ncVar.setFill(true, (void *)&ui16);
        } else if (nt == NC_INT) {
            i32 = fill_value_dbl;
            ncVar.setFill(true, (void *)&i32);
        } else if (nt == NC_UINT) {
            ui32 = fill_value_dbl;
            ncVar.setFill(true, (void *)&ui32);
        } else if (nt == NC_FLOAT) {
            f32 = fill_value_dbl;
            ncVar.setFill(true, (void *)&f32);
        } else {
            ncVar.setFill(true, (void *)&fill_value_dbl);
        }
    }

    /* vary chunck size based on dimensions */
    int do_deflate = 0;
    vector<size_t> chunkVec{CHUNKLINES, CHUNKBANDS, CHUNKPIXELS};
    if (varVec.size() == 3 && (strncmp(sname, "sci_", 4) == 0)) {
        // size_t dimLines = varVec[0].getSize();
        size_t dimBands = varVec[1].getSize();
        size_t dimPixels = varVec[2].getSize();

        // if (dimLines < CHUNKLINES)
        // chunkVec[0] = dimLines;
        if (dimBands < CHUNKBANDS)
            chunkVec[1] = dimBands;
        if (dimPixels < CHUNKPIXELS)
            chunkVec[2] = dimPixels;

        do_deflate = 1;
    }

    /* Set compression */
    if (do_deflate) {
        /* First set chunking */
        try {
            ncVar.setChunking(ncVar.nc_CHUNKED, chunkVec);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure setting chunking: " << sname << endl;
            exit(1);
        }

        try {
            ncVar.setCompression(true, true, 5);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure setting compression: " << sname << endl;
            exit(1);
        }
    }

    /* Add a "long_name" attribute */
    try {
        ncVar.putAtt("long_name", lname);
    } catch (NcException &e) {
        e.what();
        cerr << "Failure creating 'long_name' attribute: " << lname << endl;
        exit(1);
    }

    if (strcmp(flag_values, "") != 0) {
        size_t curPos = 0;

        string fv;
        fv.assign(flag_values);
        size_t pos = fv.find("=", curPos);
        fv = fv.substr(pos + 1);

        size_t semicln = fv.find(";");
        pos = 0;

        int8_t vec[1024];
        int n = 0;
        while (pos != semicln) {
            pos = fv.find(",", curPos);
            if (pos == string::npos)
                pos = semicln;

            string flag_value;
            istringstream iss(fv.substr(curPos, pos - curPos));
            iss >> skipws >> flag_value;
            vec[n++] = atoi(flag_value.c_str());
            curPos = pos + 1;
        }

        try {
            ncVar.putAtt("flag_values", nt, n, vec);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure creating 'flag_values' attribute: " << lname << endl;
            exit(1);
        }
    }

    /* Add a "flag_meanings" attribute if specified*/
    if (strcmp(flag_meanings, "") != 0) {
        try {
            ncVar.putAtt("flag_meanings", flag_meanings);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure creating 'flag_meanings' attribute: " << flag_meanings << endl;
            exit(1);
        }
    }

    /* Add a "reference" attribute if specified*/
    if (strcmp(reference, "") != 0) {
        try {
            ncVar.putAtt("reference", reference);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure creating 'reference' attribute: " << reference << endl;
            exit(1);
        }
    }

    /* Add "valid_min/max" attributes if specified */
    if (low < high) {
        switch (nt) { /* Use the appropriate number type */
            case NC_BYTE: {
                uint8_t vr[2];
                vr[0] = (uint8_t)low;
                vr[1] = (uint8_t)high;

                try {
                    ncVar.putAtt("valid_min", NC_BYTE, 1, &vr[0]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
                    exit(1);
                }

                try {
                    ncVar.putAtt("valid_max", NC_BYTE, 1, &vr[1]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
                    exit(1);
                }
            } break;
            case NC_UBYTE: {
                uint8_t vr[2];
                vr[0] = (uint8_t)low;
                vr[1] = (uint8_t)high;

                try {
                    ncVar.putAtt("valid_min", NC_UBYTE, 1, &vr[0]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
                    exit(1);
                }

                try {
                    ncVar.putAtt("valid_max", NC_UBYTE, 1, &vr[1]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
                    exit(1);
                }
            } break;
            case NC_SHORT: {
                int16_t vr[2];
                vr[0] = (int16_t)low;
                vr[1] = (int16_t)high;

                try {
                    ncVar.putAtt("valid_min", NC_SHORT, 1, &vr[0]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
                    exit(1);
                }

                try {
                    ncVar.putAtt("valid_max", NC_SHORT, 1, &vr[1]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
                    exit(1);
                }
            } break;
            case NC_USHORT: {
                uint16_t vr[2];
                vr[0] = (uint16_t)low;
                vr[1] = (uint16_t)high;

                try {
                    ncVar.putAtt("valid_min", NC_USHORT, 1, &vr[0]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
                    exit(1);
                }

                try {
                    ncVar.putAtt("valid_max", NC_USHORT, 1, &vr[1]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
                    exit(1);
                }
            } break;
            case NC_INT: {
                int32_t vr[2];
                vr[0] = (int32_t)low;
                vr[1] = (int32_t)high;

                try {
                    ncVar.putAtt("valid_min", NC_INT, 1, &vr[0]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
                    exit(1);
                }

                try {
                    ncVar.putAtt("valid_max", NC_INT, 1, &vr[1]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
                    exit(1);
                }

            } break;
            case NC_UINT: {
                uint32_t vr[2];
                vr[0] = (uint32_t)low;
                vr[1] = (uint32_t)high;

                try {
                    ncVar.putAtt("valid_min", NC_UINT, 1, &vr[0]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
                    exit(1);
                }

                try {
                    ncVar.putAtt("valid_max", NC_UINT, 1, &vr[1]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
                    exit(1);
                }

            } break;
            case NC_FLOAT: {
                float vr[2];
                vr[0] = (float)low;
                vr[1] = (float)high;

                try {
                    ncVar.putAtt("valid_min", NC_FLOAT, 1, &vr[0]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
                    exit(1);
                }

                try {
                    ncVar.putAtt("valid_max", NC_FLOAT, 1, &vr[1]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
                    exit(1);
                }
            } break;
            case NC_DOUBLE: {
                double vr[2];
                vr[0] = low;
                vr[1] = high;

                try {
                    ncVar.putAtt("valid_min", NC_DOUBLE, 1, &vr[0]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
                    exit(1);
                }

                try {
                    ncVar.putAtt("valid_max", NC_DOUBLE, 1, &vr[1]);
                } catch (NcException &e) {
                    e.what();
                    cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
                    exit(1);
                }
            } break;
            default:
                fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
                fprintf(stderr, "Got unsupported number type (%d) ", nt);
                fprintf(stderr, "while trying to create NCDF variable, \"%s\", ", sname);
                return (1);
        }
    }

    /* Add a "units" attribute if one is specified */
    if (units != NULL && *units != 0) {
        try {
            ncVar.putAtt("units", units);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure creating 'units' attribute: " << units << endl;
            exit(1);
        }
    }

    /* Add a "standard_name" attribute if one is specified */
    if (standard_name != NULL && *standard_name != 0) {
        try {
            ncVar.putAtt("standard_name", standard_name);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure creating 'standard_name' attribute: " << standard_name << endl;
            exit(1);
        }
    }

    return 0;
}

int eight20(uint8_t *inbytes, uint32_t *outsamples) {
    // This routine takes groups of 9 20-bit samples that are packed into
    // 23 bytes and unpacks them into a 4-byte integer array

    for (size_t i = 0; i < 5; i++) {
        outsamples[i * 2] = 4096 * inbytes[i * 5] + 16 * inbytes[i * 5 + 1] + inbytes[i * 5 + 2] / 16;
    }

    for (size_t i = 0; i < 4; i++) {
        outsamples[i * 2 + 1] =
            65536 * (inbytes[i * 5 + 2] % 16) + 256 * inbytes[i * 5 + 3] + inbytes[i * 5 + 4];
    }

    return 0;
}
