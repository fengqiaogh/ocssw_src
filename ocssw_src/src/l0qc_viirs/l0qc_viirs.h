#include <stdint.h>
#include <stdio.h>
int init_packet_read(FILE *infile, uint8_t epacket[], int *len_packet, long int *endfile);
int read_packet(FILE *infile, uint8_t packet[], int *len, long int *endfile);
int get_viirs_packet_time(uint8_t p1[], int32_t *pyear, int32_t *eyear, int32_t *syear, int32_t *pday, int32_t *eday, int32_t *sday, double *stime, double *ptime, double *etime);
int read_viirs_scan_packets(FILE *infile, uint8_t epacket[], int *len_packet, uint8_t pbuffer[], int *npkts, long int *endfile);
