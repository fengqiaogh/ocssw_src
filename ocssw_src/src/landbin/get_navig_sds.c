#include "proto.h"
#include "mfhdf.h"

int get_navig_sds_line(int32 sd_id, int iline, float32 *sen_mat, float32 *orb_vec, float32 *sun_ref, float32 *scan_ell, char newfile) {
    static int32 sds_id_sen_mat, sds_id_orb_vec, sds_id_sun_ref, sds_id_scan_ell;

    if (newfile) {
        if ((sds_id_sen_mat = open_sds(sd_id, "sen_mat")) == -1 ||
                (sds_id_orb_vec = open_sds(sd_id, "orb_vec")) == -1 ||
                (sds_id_sun_ref = open_sds(sd_id, "sun_ref")) == -1 ||
                (sds_id_scan_ell = open_sds(sd_id, "scan_ell")) == -1)
            return -1;
    }

    /* Read sensor matrix */
    if (read_sds_block(sds_id_sen_mat, iline, 1, sen_mat) == -1) return -1;

    /* Read orbit position vector  */
    if (read_sds_block(sds_id_orb_vec, iline, 1, orb_vec) == -1) return -1;

    /* Read reference sun vector */
    if (read_sds_block(sds_id_sun_ref, iline, 1, sun_ref) == -1) return -1;

    /* Read scan-track ellipse coefs */
    if (read_sds_block(sds_id_scan_ell, iline, 1, scan_ell) == -1) return -1;

    return 0;

}

int get_navig_sds(int32 sd_id, float32 *sen_mat, float32 *orb_vec, float32 *sun_ref, float32 *scan_ell, int32 *msec) {

    /* Read sensor matrix */
    if (read_sds(sd_id, "sen_mat", sen_mat) == -1) return -1;

    /* Read orbit position vector  */
    if (read_sds(sd_id, "orb_vec", orb_vec) == -1) return -1;

    /* Read reference sun vector */
    if (read_sds(sd_id, "sun_ref", sun_ref) == -1) return -1;

    /* Read scan-track ellipse coefs */
    if (read_sds(sd_id, "scan_ell", scan_ell) == -1) return -1;

    /* Read scan line time */
    if (read_sds(sd_id, "msec", msec) == -1) return -1;

    return 0;
}
