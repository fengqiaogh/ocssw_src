/*

 $Header$
 $Log$



From yeh@calval  Tue Mar 19 16:04:26 1996
Received: from calval.gsfc.nasa.gov by manua.gsfc.nasa.gov via ESMTP (950911.SGI.8.6.12.PATCH825/931108.SGI.AUTO.ANONFTP)
        for <frank@manua.gsfc.nasa.gov> id QAA24118; Tue, 19 Mar 1996 16:04:25 -0500
Received: by calval.gsfc.nasa.gov (950911.SGI.8.6.12.PATCH825/931108.SGI.AUTO.ANONFTP)
        for frank@manua id QAA03350; Tue, 19 Mar 1996 16:04:24 -0500
From: yeh@calval (Eueng-nan Yeh)
Message-Id: <199603192104.QAA03350@calval.gsfc.nasa.gov>
Subject: out-of-band subroutine
To: frank
Date: Tue, 19 Mar 1996 16:04:18 -0500 (EST)
Reply-To: yeh@calval (Eueng-nan Yeh)
Organization: NASA/Goddard/Laboratory for Hydrospheric Processes/SeaWiFS
X-Mailer: ELM [version 2.4 PL0]
MIME-Version: 1.0
Content-Type: text/plain; charset=US-ASCII
Content-Length: 2794      
Status: RO

 */


#include <stdio.h>
#define NBAND 8
#define NSAMPLE 1285


#ifdef __STDC__
int out_band_corr(float r[8][1285], float oxygen, int xsample)
#else

int out_band_corr(r, oxygen, xsample)
/************************************************************2*******
/To perform out-of-band correction
/-------------------------------------------------------------
/Inputs:
/  r  float[8][1285]  level-1 radiances (band 1 to 8) per line.
/  oxygen  float      Oxygen correction factor; =1 if r(765nm) already 
/                     corrected; =1.12 if not.
/  xsample int        Number of data points per line (max 1285).
/Outputs:
/  r  float[8][1285]  Out-of-band corrected radiances per line.
/  Return value=-1 for open data file failure; otherwise O.K.
/-------------------------------------------------------------
/
/  Created by Eueng-nan Yeh      GSC/SAIC          10/10/95
/
/%***********************************************************2******/
float r[NBAND][NSAMPLE];
float oxygen;
int xsample;
#endif /* __STDC__ */
{
    int i, m, id = 1;

    float wl[9] = {412., 443., 490., 510., 555., 670., 765., 865., 965.};
    /*Spectral response of band 8 divided into six components. */
    float s[6] = {0.583, 0.829, 6.251, 13.533, 1013.592, 9.374};
    float s5[3] = {5.261, 488.215, 7.475};
    float s6[3] = {2.592, 491.802, 2.112};
    /*Out of band ratio==within-band response/total response*/
    float br[8] = {0.9943, 0.9949, 0.9930, 0.9934, 0.9734, 0.9856, 0.9844, 0.9417};
    float x, ra, r5, r6, r8, r9, br5, br6, br8;
    if (id == 1) { /*simple method */
        ra = (wl[8] - wl[7]) / (wl[7] - wl[6]);
        for (m = 0; m < xsample; m++) { /*process one line*/
            r5 = r[4][m];
            r6 = r[5][m];
            r8 = r[7][m];
            for (i = 0; i < 8; i++) r[i][m] *= br[i]; /*adjust 0-7 first */
            r9 = r[7][m] + ra * (r[7][m] - oxygen * r[6][m]); /*r[6] for wl=765nm */
            x = r[7][m] * s[4];
            br8 = x / (r[0][m] * s[0] + r[3][m] * s[1] + r[4][m] * s[2]
                    + oxygen * r[6][m] * s[3] + x + r9 * s[5]);
            r[7][m] = r8* br8;
            x = r[5][m] * s6[1]; /* band 6*/
            br6 = x / (r[2][m] * s6[0] + x + r[7][m] * s6[2]);
            r[5][m] = r6*br6;
            x = r[4][m] * s5[1]; /* band 5*/
            br5 = x / (r[1][m] * s5[0] + x + r[5][m] * s5[2]);
            r[4][m] = r5*br5;
        }
    }
    return (0);
}

#ifdef TEST_OOB

void main() {
    int i, j;
    /* Typical SeaWiFS radiances TM 1 Table 1 or TM 22 Table 26.  */
    float rad[NBAND][NSAMPLE] = {
        {9.10, 9.10, 9.10},
        {8.41, 8.41, 8.41}
        ,
        {6.56, 6.56, 6.56},
        {5.64, 5.64, 5.64}
        ,
        {4.57, 4.57, 4.57},
        {2.46, 2.46, 2.46}
        ,
        {1.61, 1.61, 1.61},
        {1.0, 1.0, 1.09}};
    float r[NBAND][NSAMPLE], oxygen = 1.12; /*oxygen absorption correc. factor*/
    float ra[NBAND];

    for (i = 0; i < 8; i++) printf("%7.5f  ", rad[i][2]);
    printf("\n");
    for (i = 0; i < 8; i++) for (j = 0; j < 3; j++) r[i][j] = rad[i][j];
    out_band_corr(r, oxygen, 3);
    for (i = 0; i < 8; i++) printf("%7.5f  ", r[i][2]);
    printf("\n");
    for (i = 0; i < 8; i++) ra[i] = r[i][2] / rad[i][2];
    for (i = 0; i < 8; i++) printf("%7.5f  ", ra[i]);
    printf("\n\n");

    exit(0);
}
#endif /* TEST_OOB */

/* Test Result

From yeh@calval  Tue Mar 19 16:09:54 1996
Received: from calval.gsfc.nasa.gov by manua.gsfc.nasa.gov via ESMTP (950911.SGI.8.6.12.PATCH825/931108.SGI.AUTO.ANONFTP)
        for <frank@manua.gsfc.nasa.gov> id QAA24798; Tue, 19 Mar 1996 16:09:53 -0500
Received: by calval.gsfc.nasa.gov (950911.SGI.8.6.12.PATCH825/931108.SGI.AUTO.ANONFTP)
        for frank@manua id QAA04098; Tue, 19 Mar 1996 16:09:52 -0500
From: yeh@calval (Eueng-nan Yeh)
Message-Id: <199603192109.QAA04098@calval.gsfc.nasa.gov>
Subject: out-of-band result
To: frank
Date: Tue, 19 Mar 1996 16:09:47 -0500 (EST)
Reply-To: yeh@calval (Eueng-nan Yeh)
Organization: NASA/Goddard/Laboratory for Hydrospheric Processes/SeaWiFS
X-Mailer: ELM [version 2.4 PL0]
MIME-Version: 1.0
Content-Type: text/plain; charset=US-ASCII
Content-Length: 220       
Status: RO

9.10000  8.41000  6.56000  5.64000  4.57000  2.46000  1.61000  1.09000  
9.04813  8.36711  6.51408  5.60278  4.44292  2.42131  1.58488  1.02651  
0.99430  0.99490  0.99300  0.99340  0.97219  0.98427  0.98440  0.94175  

 */
