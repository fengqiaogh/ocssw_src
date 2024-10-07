/*-----------------------------------------------------------------------------
    Program:   viirs_sim_sdr

    Description:  Create a VIIRS SDR file set using simulated data

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       argc          I    count of command line args - 3
        char *[]  argv          I    command line arguments:
                                     [1] input geo data file name
                                     [2] id to use in the output file name
                                     and options

    Modification history:

    W. Robinson, SAIC  08 Oct 2008  Original development
    W. Robinson, SAIC  21 Sep 2010  switch input to clo library

----------------------------------------------------------------------------*/
#include "viirs_sim_sdr.h"
#include "l12_parms.h"
#include <stdio.h>

int main(int argc, char *argv[]) {
    int iscn, ipx, ibnd, idet, n_bad, msec;
    ctl_struc ctl;
    sdr_info_struc sdr_info;
    in_rec_struc in_rec;
    out_rec_struc out_rec;
    FILE *meta_id;
    static float f0bar[MAX_BND];
    double esdist;
    double esdist_(int *, int *, int *);

    /*
     *  shorten output buffering
     */
    setvbuf(stdout, NULL, _IOLBF, 0);
    setvbuf(stderr, NULL, _IOLBF, 0);

    if (viirs_sim_input(argc, argv, &ctl) != 0) exit(1);
    /*
     *  set the origin (distributor) and domain 
     */
    strcpy(sdr_info.origin, ctl.id_origin);
    strcpy(sdr_info.domain, ctl.id_domain);
    /*
     *  read in the VIIRS band center wavelengths, average F0
     */
    if (bnd_ix_2_sen_info("Lambda", (void *) in_rec.lam_band) < 0) {
        printf("%s, %d: failure to read sensor information\n",
                __FILE__, __LINE__);
    }
    out_rec.lam_band = in_rec.lam_band;
    if (bnd_ix_2_sen_info("Fobar", (void *) f0bar) < 0) {
        printf("%s, %d: failure to read sensor information\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  set up the input file (get attributes, ready for dataset reading)
     */
    if (rd_sim_init(&ctl, &sdr_info, &in_rec) != 0) return 1;
    /*
     *  get an initial fo for use in init_sdr
     */
    in_rec.year = sdr_info.year;
    in_rec.yday = sdr_info.day;
    msec = *sdr_info.scan_time * 1000.;
    esdist = esdist_(&(in_rec.year), &(in_rec.yday), &msec);
    for (ibnd = 0; ibnd < in_rec.nbnd; ibnd++)
        in_rec.f0[ibnd] = RAD_CGS_2_MKS * pow(1.0 / esdist, 2) * f0bar[ibnd];
    out_rec.f0 = &(in_rec.f0[0]);
    printf("%s, %d: info: earth sun distance correction is %f\n",
            __FILE__, __LINE__, esdist);
    /*
     *  proceed to create the SDR file with waiting data arrays
     */
    if (init_sdr(&ctl, &sdr_info, &in_rec, &out_rec) != 0) return 1;
    /*
     * loop through the lines and create the data for the SDR
     */
    for (iscn = 0; iscn < in_rec.nscan; iscn++) {
        /*
         *  the SDR scan read will also get and combine land and cloud info
         */
        if (rd_sdr_scan(iscn, &ctl, &sdr_info, &in_rec) != 0) return 1;
        /*
         *  set some times in scan and create the f0 for the scan time
         */
        msec = *(sdr_info.scan_time + iscn) * 1000.;
        in_rec.msec = msec;
        esdist = esdist_(&(in_rec.year), &(in_rec.yday), &msec);
        /*  note F0 conversion below from loadl1.c  */
        for (ibnd = 0; ibnd < in_rec.nbnd; ibnd++)
            in_rec.f0[ibnd] = RAD_CGS_2_MKS * pow(1.0 / esdist, 2) * f0bar[ibnd];

        out_rec.year = in_rec.year;
        out_rec.yday = in_rec.yday;
        out_rec.msec = in_rec.msec;
        out_rec.f0 = &(in_rec.f0[0]);
        /*
         *  note any bad values in the scan - only for VNIR
         */
        n_bad = 0;
        for (idet = 0; idet < in_rec.ndet_scan; idet++)
            for (ipx = 0; ipx < in_rec.npix; ipx++)
                for (ibnd = 0; ibnd < N_VNIR_BND; ibnd++) {
                    if (*(in_rec.bnd_q[ibnd] + ipx + in_rec.npix * idet) != 0)
                        n_bad++;
                }
        printf("%s %d: scan: %d, # bad samples: %d\n", __FILE__, __LINE__,
                iscn, n_bad);
        /*
         *  This is where any artifact addition is done and possible conversion
         *  to aggregated
         */
        if ((ctl.any_artifact == 1) && (in_rec.scn_fmt != 0))
            if (mod_artifact(&ctl, &in_rec) != 0) return 1;
        /*
         *  convert scan format
         */
        if (scan_cvt(&in_rec, &out_rec) != 0)
            return 1;
        /*
         *  should we want to remove a vicarious gain, to be added in l2gen,
         *  do it here
         */
        if (ctl.vic_cal_chg == 1) {
            for (ibnd = 0; ibnd < in_rec.nbnd; ibnd++)
                for (idet = 0; idet < out_rec.ndet_scan; idet++)
                    for (ipx = 0; ipx < out_rec.npix; ipx++)
                        if (*(out_rec.bnd_q[ ibnd ] + ipx + idet * out_rec.ndet_scan)
                                == 0) {
                            *(out_rec.bnd_lt[ ibnd ] + ipx + idet * out_rec.npix) -=
                                    ctl.offset[ibnd];
                            *(out_rec.bnd_lt[ ibnd ] + ipx + idet * out_rec.npix) /=
                                    ctl.gain[ibnd];
                        }
        }
        /*
         *  if requested, add the bow tie artifact
         *  Only do so if scan format indicates aggregated data
         */
        if ((ctl.bowtie_opt == 1) && (out_rec.scn_fmt == 0)) {
            for (idet = 0; idet < NDET; idet++) {
                switch (idet) {
                case 0:
                case 15:
                    for (ibnd = 0; ibnd < in_rec.nbnd; ibnd++) {
                        for (ipx = 0; ipx < 1008; ipx++)
                            *(out_rec.bnd_q[ibnd] + idet * 3200 + ipx) = 1;
                        for (ipx = 2192; ipx < 3200; ipx++)
                            *(out_rec.bnd_q[ibnd] + idet * 3200 + ipx) = 1;
                    }
                    break;
                case 1:
                case 14:
                    for (ibnd = 0; ibnd < in_rec.nbnd; ibnd++) {
                        for (ipx = 0; ipx < 640; ipx++)
                            *(out_rec.bnd_q[ibnd] + idet * 3200 + ipx) = 1;
                        for (ipx = 2560; ipx < 3200; ipx++)
                            *(out_rec.bnd_q[ibnd] + idet * 3200 + ipx) = 1;
                    }
                    break;
                }
            }
        }
        /*
         *  and output the data
         */
        if (wr_sdr_scan(iscn, &out_rec) != 0) return 1;
    }
    /*
     *  finish, close files, write file info and exit
     */
    if (fin_sdr(&ctl, &in_rec, &out_rec) != 0) return 1;

    for (iscn = 0; iscn < out_rec.nbnd + 1; iscn++)
        printf("Created file # %d: %s\n", iscn, sdr_info.sdr_files[iscn]);
    printf("Path: %s\n", ctl.out_loc);
    printf("Start date: %s\n", sdr_info.st_date);
    printf("Start time: %s\n", sdr_info.st_time);
    printf("End date: %s\n", sdr_info.en_date);
    printf("End time: %s\n", sdr_info.en_time);
    if (ctl.meta_use != 0) {
        if ((meta_id = fopen(ctl.meta_file, "w")) == NULL) {
            printf("%s, line %d: Failure to create the metadata file\n",
                    __FILE__, __LINE__);
            return 1;
        }
        for (iscn = 0; iscn < out_rec.nbnd + 1; iscn++)
            fprintf(meta_id, "%s\n", sdr_info.sdr_files[iscn]);
        fprintf(meta_id, "%s\n", ctl.out_loc);
        fprintf(meta_id, "%s\n", sdr_info.st_date);
        fprintf(meta_id, "%s\n", sdr_info.st_time);
        fprintf(meta_id, "%s\n", sdr_info.en_date);
        fprintf(meta_id, "%s\n", sdr_info.en_time);
        fclose(meta_id);
    }

    printf("Successfully completed simulated SDR generation\n");
    return 0;
}
