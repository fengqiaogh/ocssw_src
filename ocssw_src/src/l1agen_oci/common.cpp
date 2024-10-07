#include <algorithm>
#include <iterator>

#include "common.h"

using namespace std;

#define BAGGBASE 2863311530
#define RAGGBASE 2774100650

//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     FutureTech     08/02/19 0.10  Original development
//  Liang Hong     SAIC           04/14/22 0.20  l1agen_oci v1.03.00
//  Liang Hong     SAIC           04/24/22 0.21  l1agen_oci v1.03.02
//  Liang Hong	   SAIC           05/05/22 0.30  l1agen_oci v1.04.00
//  Liang Hong     SAIC           11/18/22 0.40  l1agen_oci v1.06.00
//  Jakob Lindo    SSAI           03/16/24 0.50  Added aggregation indication

int get_band_dims(uint8_t *apacket, uint16_t &ncp, uint16_t &nbb, uint16_t &nrb, uint16_t &nsp, uint16_t &ndc,
                  uint16_t &nds, uint16_t *btaps, uint16_t *rtaps, itab *itable) {
    // Extract spatial aggregation table and compute numbers of pixels
    short ioff = 36;
    short nagg[4] = {1, 2, 4, 8};
    ncp = 0;
    nsp = 0;
    ndc = 0;
    nds = 0;

    for (size_t i = 0; i < 10; i++) {
        itable[i].dtype = apacket[ioff + 3] % 16;
        itable[i].iagg = apacket[ioff + 2] % 4;
        if (itable[i].dtype == 5)
            itable[i].iagg = 0;
        itable[i].lines = apacket[ioff] * 256 + apacket[ioff + 1];
        ioff += 4;

        if (itable[i].dtype > 0 && itable[i].dtype <= 14 &&
            itable[i].dtype != 10) {  // Changed dtype<=12 to 14, LH, 3/30/2022
            if (itable[i].dtype == 2) {
                ndc += itable[i].lines / nagg[itable[i].iagg];
                nds += itable[i].lines / 8;
            }  // else {    Changed to include dark pixels, Liang Hong, 5/12/2020
            ncp += itable[i].lines / nagg[itable[i].iagg];
            nsp += itable[i].lines / 8;
            // }
        }
    }

    // to ensure that the science and dark count arrays will be created with at least one pixel
    if (ncp == 0) {
        ncp = 1;
        nsp = 1;
    }

    if (ndc == 0) {
        ndc = 1;
        nds = 1;
    }

    ioff += 4;

    // Extract spectral aggregation and compute numbers of bands
    // Tap enable flags
    short btap = apacket[ioff + 2] * 256 + apacket[ioff + 3];
    short rtap = apacket[ioff] * 256 + apacket[ioff + 1];

    // Tap aggregation factors
    uint32_t bagg;
    uint32_t ragg;
    memcpy(&bagg, &apacket[ioff + 8], sizeof(uint32_t));
    memcpy(&ragg, &apacket[ioff + 4], sizeof(uint32_t));
    bagg = SWAP_4(bagg);
    ragg = SWAP_4(ragg);

    if (itable[1].dtype == 1 && (bagg != BAGGBASE || ragg != RAGGBASE))
        itable[1].dtype = 13;

    // Compute number of bands for enabled taps
    // Bands are in reverse spectral order
    nbb = 0;
    nrb = 0;
    uint16_t ken = 1;
    uint32_t kag = 3;
    uint32_t lag = 1;

    //  for (size_t i=0; i<16; i++) {
    for (int i = 15; i >= 0; i--) {
        btaps[i] = (btap & ken) / ken;
        if (btaps[i])
            nbb += 32 / nagg[(bagg & kag) / lag];
        rtaps[i] = (rtap & ken) / ken;
        if (rtaps[i])
            nrb += 32 / nagg[(ragg & kag) / lag];
        ken *= 2;
        kag *= 4;
        lag *= 4;
    }

    return 0;
}

int get_anc_packet_time(uint8_t *apacket, int32_t &iyear, int32_t &iday, double &stime) {
    // Unpack and convert the CCSDS segmented time code
    // from the OCI ancillary packet

    // Get day count since Jan. 1, 1958 (Julian day 2436205)
    // sec58 = seconds since 0 UT Jan. 1, 1958
    short int toff = 28;
    uint32_t sec58;
    memcpy(&sec58, &apacket[toff], sizeof(uint32_t));
    sec58 = SWAP_4(sec58);

    double dbl58 = (double)sec58;
    int leap = leapseconds_since_1993(dbl58) + 27;
    sec58 -= leap;

    // Convert to year and day
    iday = sec58 / 86400;
    int32_t jd = iday + 2436205;  // Jan. 1, 1958 is Julian day 2436205
    jdate(jd, &iyear, &iday);

    // Get microseconds
    uint32_t usec;
    memcpy(&usec, &apacket[toff + 4], sizeof(uint32_t));
    usec = SWAP_4(usec);
    usec = usec / 4096;  // 20 MSBs used, last 12 are spares

    uint32_t isec = sec58 % 86400;
    stime = isec + ((double)usec) * 1e-6;

    return 0;
}

int read_oci_scan_packets(L0Stream *tfileStream, uint8_t *apacket, uint8_t (*pbuffer)[PKTSIZE],
                          uint32_t &npkts, int32_t &spnum, int32_t &ancind, vector<int32_t> &tlmind,
                          uint8_t &seqerr, int32_t &endfile, bool isSPW) {
    // int pos = tfileStream->tellg();
    if (endfile)
        return 0;  // LH, 4/24/2022, v1.03.02

    uint32_t apidmin = 636;  // ver
    uint32_t apidmax = 745;  // ver 1.00.00, LH, 1/7/2022
    uint32_t ui32;
    uint32_t len;
    uint32_t maxpkts = 30000;  // ver 0.20, LH, 4/14/2022

    // Ancillary packet APID
    static uint32_t apida = 636;

    // APIDs with time fields and spin numbers
    static uint32_t apidt[8] = {
        711, 712, 713, 715,
        716, 717, 721, 723};  //  telemetry APIDs with time fields and spin numbers, ver 0.30, LH, 5/5/2022

    // APIDs with spin numbers but no time fields
    static uint32_t apidn[2] = {700, 720};

    // Get OCI spin number from first packet of next scan
    uint32_t apid = (apacket[0] % 8) * 256 + apacket[1];

    npkts = 0;
    tlmind.clear();
    ancind = -1;
    // int spne = 0;

    if (apid == apida) {
        // Ancillary packet
        memcpy(&ui32, &apacket[24], 4);
        spnum = SWAP_4(ui32);
    } else if (apid == apidn[0] || apid == apidn[1]) {
        // Science Packet without time field, ignore SWIR packets for now
        memcpy(&ui32, &apacket[6], 4);
        spnum = SWAP_4(ui32);
        // } else if (apid == apidt[0] || apid == apidt[1] || apid == apidt[2] ||apid == apidt[3]) {
    } else if (std::find(std::begin(apidt), std::end(apidt), apid) != std::end(apidt)) {
        // Packet with time field
        memcpy(&ui32, &apacket[12], 4);
        spnum = SWAP_4(ui32);
    }

    //  uint32_t spn = spnum;
    // uint32_t spnt = spnum;
    int32_t spn = spnum;
    int32_t spnt = spnum;
    uint8_t *packet = apacket;

    len = apacket[4] * 256 + apacket[5] + 7;

    int prev_pseq = 0;
    seqerr = 0;

    // Read all of the packets with this spin number
    // and store in the packet buffer
    while (spn <= spnum && spnt <= (spnum + 1) && (npkts < maxpkts)) {
        // Check for packet out of order
        // Load packet into buffer
        memcpy(&pbuffer[npkts][0], packet, len);
        apid = (packet[0] % 8) * 256 + packet[1];

        // Check science packet sequence numbers, LH 11/20/2020
        if (apid == apidn[0]) {
            int pseq = (pbuffer[npkts][2] % 64) * 256 + pbuffer[npkts][3];
            if ((prev_pseq > 0) && ((pseq - prev_pseq) != 1 && (pseq - prev_pseq) != -16383)) {
                // cout<<"npkts="<<npkts<<"pseq,pseq-prev_pseq="<<pseq<<","<<pseq-prev_pseq<<endl;
                seqerr = 1;
            }
            prev_pseq = pseq;
        }

        if (apid == apida)
            ancind = npkts;
        //      if (pos >=  527223898) {
        //  pos = tfileStream->tellg();
        // cout << pos << " " << apid << " " << apida << " " << ancind
        //     << spnum << " " << spn << " " << spnt << endl;
        //}

        if (apid >= apidmin & apid <= apidmax && apid != apidn[0] && apid != apidn[1]) {
            tlmind.push_back(npkts);
        }
        npkts++;
        if (npkts == maxpkts) {
            cout << "Maximum number of packets: " << maxpkts << " exceeded." << endl;
            // exit (1);  // Liang Hong, 08/19/2020
        }
        // cout << apid << " " << npkts << endl;

        // Read next packet
        read_packet(tfileStream, NULL, len, apid, endfile, isSPW);
        if (endfile)
            return 0;

        if (len > PKTSIZE) {
            cout << "Packet size > " << PKTSIZE << " (" << len << ")" << endl;
            uint8_t *dummy_buf = new uint8_t[len];
            read_packet(tfileStream, dummy_buf, len, apid, endfile, isSPW);
            delete[] dummy_buf;
            //  exit(1);
        } else {
            read_packet(tfileStream, packet, len, apid, endfile, isSPW);
            apid = (packet[0] % 8) * 256 + packet[1];
        }

        while (apid < apidmin || apid > apidmax) {
            read_packet(tfileStream, NULL, len, apid, endfile, isSPW);
            if (endfile)
                return 0;

            if (len > PKTSIZE) {
                cout << "Packet size > " << PKTSIZE << " (" << len << ")" << endl;
                uint8_t *dummy_buf = new uint8_t[len];
                read_packet(tfileStream, dummy_buf, len, apid, endfile, isSPW);
                delete[] dummy_buf;
                //        exit(1);
            } else {
                read_packet(tfileStream, packet, len, apid, endfile, isSPW);
                apid = (packet[0] % 8) * 256 + packet[1];
            }
        }

        // Get spin number
        if (apid == apida) {
            memcpy(&ui32, &packet[24], 4);
            spn = SWAP_4(ui32);
        } else if (apid == apidn[0] || apid == apidn[1]) {  // ver 1.00.01
            memcpy(&ui32, &packet[6], 4);
            spn = SWAP_4(ui32);
            //    } else if (apid == apidt[0] || apid == apidt[1] || apid == apidt[2] || apid == apidt[3]) {
        } else if (std::find(std::begin(apidt), std::end(apidt), apid) != std::end(apidt)) {
            memcpy(&ui32, &packet[12], 4);
            spnt = SWAP_4(ui32);
        }  // if (apid == apidn[0] || apid == apidn[1]) {

    }  // while (spn <= spnum && spnt <= spnum)

    // Report spin out of order
    //  if (spn < spnum)
    //  cout << "Next packet spin number: " << spn
    //       << " out of order for spin: " << spnum << endl;
    //  cout << "npkts: " << npkts << endl;

    return 0;
}

int anc_compare(uint8_t *apacket0, uint8_t *apacket) {
    int anc_compare = 1;

    double stime;
    int32_t iyear, iday;
    uint32_t ui32;

    if (apacket[15] > 0)
        return anc_compare;  // ancillary packet without valid mode table, v0.20

    // Compare spatial data collection fields
    int ioff = 36;
    size_t ilen = 40;

    uint16_t aggi = 0;
    for (size_t i = 0; i < ilen; i++)
        if (apacket[ioff + i] != apacket0[ioff + i])
            aggi++;

    if (aggi > 0) {
        get_anc_packet_time(apacket, iyear, iday, stime);
        memcpy(&ui32, &apacket[24], 4);
        int32_t spn = SWAP_4(ui32);
        uint16_t ih = (uint16_t)floor(stime / 3600);
        uint16_t mn = (uint16_t)floor((stime - ih * 3600) / 60);
        uint16_t sec = (uint16_t)floor(stime - ih * 3600 - mn * 60);
        cout << "Spatial table change at: spin=" << spn << ", " << ih << ":" << mn << ":" << sec << endl;
        anc_compare = 0;
    }

    // Compare spectral data collection fields
    //  int joff = 76;
    // size_t jlen = 16;
    int joff = 80;
    size_t jlen = 12;

    uint16_t aggj = 0;
    for (size_t i = 0; i < jlen; i++)
        if (apacket[joff + i] != apacket0[joff + i])
            aggj++;

    if (aggj > 0) {
        get_anc_packet_time(apacket, iyear, iday, stime);
        memcpy(&ui32, &apacket[24], 4);
        int32_t spn = SWAP_4(ui32);
        uint16_t ih = (uint16_t)floor(stime / 3600);
        uint16_t mn = (uint16_t)floor((stime - ih * 3600) / 60);
        uint16_t sec = (uint16_t)floor(stime - ih * 3600 - mn * 60);
        cout << "Spectral table change at: spin=" << spn << ", " << ih << ":" << mn << ":" << sec << endl;
        anc_compare = 0;
    }

    return anc_compare;
}

int read_packet(L0Stream *tfileStream, uint8_t *packet, uint32_t &len, uint32_t &apid, int32_t &endfile,
                bool isSPW) {
    if (tfileStream->eof()) {
        endfile = 1;
        len = 0;
        cout << "End of packet file" << endl;
        return 0;
    }

    // Read spacewire header if needed
    if (isSPW) {
        uint8_t spwhead[2];
        tfileStream->read((char *)&spwhead, 2);
    }

    // Read packet header
    uint8_t phead[6];
    if (packet == NULL) {
        tfileStream->read((char *)&phead, 6);

        // Get length of packet body and APID
        len = phead[4] * 256 + phead[5] + 1 + 6;
        apid = (phead[0] % 8) * 256 + phead[1];

        if (tfileStream->tellg() == -1)
            endfile = 1;
        tfileStream->seekg(-6, ios_base::cur);
        if (isSPW)
            tfileStream->seekg(-2, ios_base::cur);
        return 0;
    }
    tfileStream->read((char *)packet, len);
    // cout << tfileStream->tellg() << endl;

    return 0;
}

int get_swir_mode(uint8_t (*pbuffer)[PKTSIZE], uint32_t npkts, uint16_t &smode) {
    // Look for a SWIR band packet

    // smode = 0;
    for (size_t ipkt = 0; ipkt < npkts; ipkt++) {
        uint32_t apid = (pbuffer[ipkt][0] % 8) * 256 + pbuffer[ipkt][1];
        //    if (npkts == 5) cout << "apid: " << apid << endl;
        if (apid == 720) {
            //  Get mode from packet
            uint8_t smeta = pbuffer[ipkt][12];
            smode = (smeta % 64) / 16;
            return 0;
        }
    }

    return 0;
    //  cout << "SWIR mode not found" << endl;
    // exit(1);
}
