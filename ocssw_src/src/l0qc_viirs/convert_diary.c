// procedure to convert spacecraft diary data into time-tagged
//  orbit and attitude data

//       Arguments
//
//       Name     Type    I/O    Description
//       ----     ----    ---    -----------
//       npkts    int      I     Number of SC diary packets
//       dstore   byte     I     71 x npkts array of SC diary packets
//       otime    double   O     Time tags of orbit (seconds of day)
//       orb      double    O     6 x npkts array of orbit vectors
//       atime    double   O     Time tags of attitude (seconds of day) 
//       quat     double    O     4 x npkts array of quaternions

// Liang Hong, Jan 14, 2016, Ported from IDL
// V0.1, Feb 1, 2016

#include <timeutils.h>

int convert_diary(int npkts, uint8_t dstore[], double otime[], double orb[], double atime[], double quat[]) {
    // Set up output arrays
    int i, j, m, ioff;
    int32_t jd0, ccsds_iy, idy, jd;
    int16_t iy;
    double secd;
    uint8_t cctime[8];

    // initialize output arrays
    for (i = 0; i < npkts; i++) {
        otime[i] = 0;
        atime[i] = 0;
        for (j = 0; j < 6; j++) orb[j + i * 6] = 0;
        for (j = 0; j < 4; j++) quat[j + i * 4] = 0;
    }

    // Get start year and day
    //jd0 = swap_endian(fix(dstore(6:7,0),0)) + 2436205 ;Days since 1/1/1958
    jd0 = (dstore[0 * 71 + 6] << 8) + dstore[0 * 71 + 7] + 2436205;

    // Loop through packets
    for (i = 0; i < npkts; i++) {
        // Get orbit time
        //ccsds_to_yds,dstore(15:22,i),iy,idy,sec
        for (m = 0; m < 8; m++) cctime[m] = dstore[i * 71 + m + 15];
        ccsds_to_yds(cctime, &ccsds_iy, &idy, &secd);
        iy = (int16_t) ccsds_iy;
        jd = jday(iy, 1, idy);
        otime[i] = secd + (jd - jd0)*86400.0;

        // Convert orbit vectors to doubles
        ioff = 23;
        for (j = 0; j < 6; j++) {
            float tmpa;
            *((uint8_t*) (&tmpa) + 0) = dstore[i * 71 + ioff + 3];
            *((uint8_t*) (&tmpa) + 1) = dstore[i * 71 + ioff + 2];
            *((uint8_t*) (&tmpa) + 2) = dstore[i * 71 + ioff + 1];
            *((uint8_t*) (&tmpa) + 3) = dstore[i * 71 + ioff];
            orb[j + i * 6] = (double) tmpa;
            ioff += 4;
        }

        // Get attitude time
        //ccsds_to_yds,dstore(47:54,i),iy,idy,sec
        for (m = 0; m < 8; m++) cctime[m] = dstore[i * 71 + m + 47];
        ccsds_to_yds(cctime, &ccsds_iy, &idy, &secd);
        iy = (int16_t) ccsds_iy;
        jd = jday(iy, 1, idy);
        atime[i] = secd + (jd - jd0)*86400.0;

        // Convert attitude quaternion to doubles
        ioff = 55;
        for (j = 0; j < 4; j++) {
            float tmpa;
            *((uint8_t*) (&tmpa) + 0) = dstore[i * 71 + ioff + 3];
            *((uint8_t*) (&tmpa) + 1) = dstore[i * 71 + ioff + 2];
            *((uint8_t*) (&tmpa) + 2) = dstore[i * 71 + ioff + 1];
            *((uint8_t*) (&tmpa) + 3) = dstore[i * 71 + ioff];
            quat[j + i * 4] = (double) tmpa;
            ioff += 4;
        }
    }
    return 0;
}
