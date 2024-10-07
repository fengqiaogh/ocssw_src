#include <cstdio>
#include <sstream>
#include <string.h>
#include <unistd.h>

#include <iomanip>

#include "l0info_oci.h"

using namespace std;

#define VERSION "0.68"

//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     SAIC           04/26/19 0.10  Original development
//                                               based on l1agen_oci
//  Joel Gales     SAIC           06/07/19 0.20  Incorporate code changes from
//                                               F.Patt
//  Joel Gales     SAIC           06/19/19 0.30  Fix single file processing
//  Joel Gales     SAIC           07/24/19 0.40  Incorporate code changes from
//                                               F.Patt
//  Joel Gales     SAIC           08/02/19 0.50  Flag and ignore packets
//                                               greater than 1200 bytes.
//                                               Remove common routines and
//                                               Place in common.cpp
//  Joel Gales     SAIC           10/25/19 0.51  Exit if EOF before finding
//                                               good packet
//  Joel Gales     SAIC           11/20/19 0.60  Implement code changes from
//                                               F. Patt
//  Joel Gales     SAIC           11/25/19 0.61  Change 60 to maxsc fpr dspn
//                                               comparison
//  Joel Gales     SAIC           11/27/19 0.62  Change L1A name to standard
//                                               Trap granules with no ancillary
//                                               granules
//  Joel Gales     SAIC           01/03/20 0.63  Update call to
//                                               read_oci_scan_packets
//  Joel Gales     SAIC           01/22/20 0.65  Don't break for gaps
//  Joel Gales     SAIC           01/23/20 0.66  Check for EOF when checking
//                                               for zero science pixels
//  Joel Gales     SAIC           01/28/20 0.67  Bug fixes and updates
//  Joel Gales     SAIC           04/02/20 0.68  Fix scomp statement

int main(int argc, char *argv[]) {
    cout << "l0info_oci " << VERSION << " (" << __DATE__ << " " << __TIME__ << ")" << endl;

    if (argc == 1) {
        cout << endl << "l0info_oci OCI_packet_file granule_len" << endl;
        return 0;
    }

    fstream tfileStream;

    // OCI packet file
    tfileStream.open(argv[optind + 0], fstream::in | fstream::binary);
    if (tfileStream.fail()) {
        cout << argv[optind + 0] << " not found" << endl;
        exit(1);
    }

    uint32_t apid = 0;
    uint32_t len = 0;
    int32_t endfile = 0;
    uint8_t fpacket[PKTSIZE];
    uint8_t apacket[ANCSIZE];
    uint8_t apacket0[ANCSIZE];
    uint32_t maxpkts = 15000;
    int acomp = 1;
    //  int maxgap = 2000;

    uint8_t **pbuffer0 = new uint8_t *[maxpkts];
    pbuffer0[0] = new uint8_t[PKTSIZE * maxpkts];
    for (size_t i = 1; i < maxpkts; i++)
        pbuffer0[i] = pbuffer0[i - 1] + PKTSIZE;

    uint32_t npkts0;
    int32_t ancind0;
    int32_t spn0, spn;
    vector<int32_t> tlmind0;

    // Get first science or ancillary packet
    read_packet(&tfileStream, NULL, len, apid, endfile);
    if (len > PKTSIZE) {
        cout << "Packet too big (" << len << ") for buffer (" << PKTSIZE << ")" << endl;
        exit(1);
    }
    read_packet(&tfileStream, fpacket, len, apid, endfile);
    apid = (fpacket[0] % 8) * 256 + fpacket[1];
    int nsk = 0;
    while (apid != 636 && apid != 700 && apid != 720 && !endfile) {
        nsk++;
        read_packet(&tfileStream, NULL, len, apid, endfile);
        if (len > PKTSIZE) {
            cout << "Packet too big (" << len << ") for buffer (" << PKTSIZE << ")" << endl;
            exit(1);
        }
        read_packet(&tfileStream, fpacket, len, apid, endfile);
        apid = (fpacket[0] % 8) * 256 + fpacket[1];
        // cout << "apid: " << apid << " " << nsk << endl;
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

    if (argc == 3) {
        string str = argv[optind + 1];
        istringstream(str) >> mper;
    } else {
        mper = 0;
    }

    itab itable[10];
    string dtypes[] = {"",     "",      "_DARK",   "_SOL-D",  "_SOL-M", "_LIN",
                       "_LUN", "_DIAG", "_SNAP-T", "_SNAP-S", "_STAT",  "_SPEC"};
    string smodes[] = {"", "_SDIAG", "_SRAW", "_STEST"};

    uint8_t **pbuffer = new uint8_t *[maxpkts];
    pbuffer[0] = new uint8_t[PKTSIZE * maxpkts];
    for (size_t i = 1; i < maxpkts; i++)
        pbuffer[i] = pbuffer[i - 1] + PKTSIZE;

    uint16_t ncps, nbbs, nrbs, nsps, ndcs, ndss, btaps[16], rtaps[16];

    uint8_t *seqerr = new uint8_t[maxpkts];
    for (size_t i = 1; i < maxpkts; i++) {
        seqerr[i] = 0;
    }
    uint8_t noseq = 255;

    ////////////////////// Main Loop ///////////////////
    while (!endfile) {
        while (ancind0 == -1) {
            read_oci_scan_packets(&tfileStream, fpacket, (uint8_t(*)[PKTSIZE]) & pbuffer0[0][0], npkts0, spn0,
                                  ancind0, tlmind0, noseq, endfile);
            if (endfile) {
                cout << "No ancillary packets found in file" << endl;
                exit(1);
            }
        }

        // Check for zero science pixels
        memcpy(apacket0, &pbuffer0[ancind0][0], ANCSIZE);
        get_band_dims(apacket0, ncps, nbbs, nrbs, nsps, ndcs, ndss, btaps, rtaps, itable);

        while (((ncps == 1 && ndcs == 1) || ancind0 == -1) && !endfile) {
            //    while ((ncps == 0 || ancind0 == -1) && !endfile) {
            read_oci_scan_packets(&tfileStream, fpacket, (uint8_t(*)[PKTSIZE]) & pbuffer0[0][0], npkts0, spn0,
                                  ancind0, tlmind0, noseq, endfile);
            // cout << ancind0 << " " << endfile << endl;
            if (ancind0 != -1) {
                memcpy(apacket0, &pbuffer0[ancind0][0], ANCSIZE);
                get_band_dims(apacket0, ncps, nbbs, nrbs, nsps, ndcs, ndss, btaps, rtaps, itable);
            }
        }
        if (ancind0 == -1) {
            cout << "No ancillary packets for last granule" << endl;
            break;
        }

        get_anc_packet_time(apacket0, iyear, iday, stime);

        double scanp = 1.0 / 5.737;
        double stimp = stime - scanp;

        time_struct starttime, endtime;
        starttime.iyear = iyear;
        starttime.iday = iday;
        starttime.sec = stime;

        uint32_t ltime;
        // uint32_t mtime;
        int32_t maxsc;

        if (mper > 0) {
            ltime = (((int32_t)(stime)) / 60 / mper) * (mper * 60);
            // mtime = ltime + mper*60;
            maxsc = mper * 400;  // increased for I&T
            if (!acomp) {
                ltime = (int32_t)stime;
            }
        } else {
            cout << "Processing with single file option" << endl;
            ltime = (int32_t)stime;
            // mtime = ltime + 600*25;
            maxsc = 3600 * 25;
        }
        acomp = 1;

        // Get SWIR band data mode
        uint16_t smode;
        get_swir_mode((uint8_t(*)[PKTSIZE]) & pbuffer0[0][0], npkts0, smode);
        uint16_t smodep = smode;
        int scomp = 1;

        // Determine start and end time of granule
        uint32_t jd0 = jday(iyear, 1, iday);

        int16_t yr16 = (int16_t)iyear;
        int16_t doy = (int16_t)iday;
        int16_t month, dom;
        yd2md(yr16, doy, &month, &dom);

        int32_t ih = (int32_t)(ltime / 3600);
        int32_t mn = (int32_t)((ltime - ih * 3600) / 60);
        int32_t isec = (int32_t)(ltime - ih * 3600 - mn * 60);

        stringstream timestr, datestr;
        timestr << setfill('0') << setw(2) << ih << setfill('0') << setw(2) << mn << setfill('0') << setw(2)
                << isec;

        //    datestr << setfill('0') << setw(4) << iyear
        //      << setfill('0') << setw(3) << iday;
        datestr << setfill('0') << setw(4) << iyear << setfill('0') << setw(2) << month << setfill('0')
                << setw(2) << dom;

        string l1a_name = string("PACE_OCI") + dtypes[itable[1].dtype] + smodes[smode] + "." + datestr.str() +
                          "T" + timestr.str() + ".L1A.nc";

        //    string l1a_name = string("OCI") + datestr.str() + timestr.str() +
        // ".L1A_PACE" + dtypes[itable[1].dtype] + ".nc";
        cout << endl << l1a_name.c_str() << endl;

        uint8_t **ancdata = new uint8_t *[maxsc + 1];
        ancdata[0] = new uint8_t[ANCSIZE * (maxsc + 1)];
        for (size_t i = 1; i < (size_t)(maxsc + 1); i++)
            ancdata[i] = ancdata[i - 1] + ANCSIZE;

        // Read and process OCI scans
        uint32_t isc = 0;
        uint32_t npkts;
        int32_t ancind;
        int32_t enddata = 0;
        int dspn = 1;

        // while ( stime < mtime && acomp && !enddata && scomp && (dspn <= maxgap)) {
        while (!enddata && scomp) {
            // while ( !enddata) {
            //         if ((isc % 100) == 0) cout << "Processing scan " << isc << endl;

            memcpy(&pbuffer[0][0], &pbuffer0[0][0], PKTSIZE * maxpkts);
            ancind = ancind0;
            npkts = npkts0;
            spn = spn0;
            enddata = endfile;

            // Read next scan
            read_oci_scan_packets(&tfileStream, fpacket, (uint8_t(*)[PKTSIZE]) & pbuffer0[0][0], npkts0, spn0,
                                  ancind0, tlmind0, seqerr[isc], endfile);

            // Check for backward time jump
            //       if (stime > stimp && npkts > 1 && isc < maxsc) {
            if (stime > stimp && npkts > 1) {
                // Save ancillary packet in array
                if (ancind != -1)
                    memcpy(ancdata[isc], pbuffer[ancind], ANCSIZE);

                isc++;
                stimp = stime;
                endtime.iyear = iyear;
                endtime.iday = iday;
                endtime.sec = stime;
            }  // if (stime > stimp ...

            //       cout << "stime: " << stime << endl;

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
            }

            dspn = spn0 - spn;
            if (dspn > maxsc)
                cout << "Spin number gap: " << spn << " " << spn0 << endl;
        }  // while ( stime < mtime && acomp && !enddata)

        cout << "Scans in file: " << isc << endl;

        if (isc > 0) {
            // Need to include ancillary packet from next scan
            if (ancind0 != -1)
                memcpy(ancdata[isc], apacket, ANCSIZE);

            //      if ((endfile && mper <= 0) || !acomp)
            //  mtime = (uint32_t) floor(stime);
        }

        string startstring, endstring;
        gen_time_string(starttime, endtime, startstring, endstring);
        cout << "Start Time=" << startstring.c_str() << endl;
        cout << "End Time  =" << endstring.c_str() << endl;

        delete[] ancdata[0];
        delete[] ancdata;

    }  // while (!endfile)

    /////////////////// End Main Loop ///////////////////

    tfileStream.close();

    return 0;
}

int gen_time_string(time_struct &starttime, time_struct &endtime, string &startstring, string &endstring) {
    int16_t mon, idm;
    int32_t ih, mn;
    stringstream ss;

    yd2md((int16_t)starttime.iyear, (int16_t)starttime.iday, &mon, &idm);

    ih = (int)(starttime.sec / 3600);
    starttime.sec -= ih * 3600;
    mn = (int)(starttime.sec / 60);
    starttime.sec -= mn * 60;

    // yyyy-mn-dyThr:mn:ss
    ss = stringstream();
    ss << setw(4) << to_string(starttime.iyear) << "-";
    ss << setw(2) << setfill('0') << mon << "-";
    ss << setw(2) << setfill('0') << idm << "T";

    ss << setw(2) << setfill('0') << ih << ":";
    ss << setw(2) << setfill('0') << mn << ":";
    ss << fixed << setw(6) << setprecision(3) << setfill('0') << starttime.sec;

    startstring = ss.str();

    yd2md((int16_t)endtime.iyear, (int16_t)endtime.iday, &mon, &idm);

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

    endstring = ss.str();

    return 0;
}
