#ifndef NAVBLK_S_H_
#define NAVBLK_S_H_

typedef struct navblk_struct {
    float orb_vec[3];
    float l_vert[3];
    float sun_ref[3];
    float att_ang[3];
    float sen_mat[3][3];
    float scan_ell[6];
    int nflag[8];
} navblk_sType;

#endif /* NAVBLK_S_H_ */
