#include <stdint.h>
#include <fstream>
#include <timeutils.h>
#include <netcdf>
#include <vector>
#include "l0stream.hpp"

#define PKTSIZE 2048  // changed from 3200 to 2048, 3/10/2022, LH

typedef struct {
    short dtype;
    short iagg;
    short lines;
} itab;

#define SWAP_2(x) ((((x) & 0xff) << 8) | ((unsigned short)(x) >> 8))

#define SWAP_4(x)                                                                         \
    ((((x) << 24) & 0xFF000000) | (((x) << 8) & 0x00FF0000) | (((x) >> 8) & 0x0000FF00) | \
     (((x) >> 24) & 0x000000FF))

int get_band_dims(uint8_t *apacket, uint16_t &ncp, uint16_t &nbb, uint16_t &nrb, uint16_t &nsp, uint16_t &ndc,
                  uint16_t &nds, uint16_t *btaps, uint16_t *rtaps, itab *itable);

int get_anc_packet_time(uint8_t *apacket, int32_t &iyear, int32_t &iday, double &stime);

int read_oci_scan_packets(L0Stream *tfileStream, uint8_t *apacket, uint8_t (*pbuffer)[PKTSIZE],
                          uint32_t &npkts, int32_t &spnum, int32_t &ancind, std::vector<int32_t> &tlmind,
                          uint8_t &seqerr, int32_t &endfile, bool isSPW);

int anc_compare(uint8_t *apacket0, uint8_t *apacket);

int read_packet(L0Stream *tfileStream, uint8_t *packet, uint32_t &len, uint32_t &apid, int32_t &endfile,
                bool isSPW);

int get_swir_mode(uint8_t (*pbuffer)[PKTSIZE], uint32_t npkts, uint16_t &smode);
