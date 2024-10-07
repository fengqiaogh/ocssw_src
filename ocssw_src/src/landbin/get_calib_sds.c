#include "proto.h"
#include "mfhdf.h"

int get_calib_sds(int32 sd_id, int16 *dark_rest, int16 *gain, int16 *tdi, int16 *scan_temp, int16 *side, int32 *msec) {

    /* Read dark restore pixel */
    if (read_sds(sd_id, "dark_rest", dark_rest) == -1) return -1;

    /* Read gain */
    if (read_sds(sd_id, "gain", gain) == -1) return -1;

    /* Read time delay and integration settings */
    if (read_sds(sd_id, "tdi", tdi) == -1) return -1;

    /* Read scan temperature */
    if (read_sds(sd_id, "scan_temp", scan_temp) == -1) return -1;

    /* Read mirror side */
    if (read_sds(sd_id, "side", side) == -1) return -1;

    /* Read scan time */
    if (read_sds(sd_id, "msec", msec) == -1) return -1;

    return 0;
}
