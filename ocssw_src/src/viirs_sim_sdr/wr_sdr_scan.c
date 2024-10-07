#include "viirs_sim_sdr.h"

int wr_sdr_scan(int iscn, out_rec_struc *out_rec)
/*-----------------------------------------------------------------------------
    Program:   wr_sdr_scan.c

    Description:  write all line oriented data to the VIIRS SDRs

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       iscn          I    scan of lines to write
        out_rec_struc * out_rec   I/O  output record controls

    Modification history:

    W. Robinson, SAIC  21 Oct 2008  Original development
    W. Robinson, SAIC  15 Mar 2010  recast to scan oriented version

----------------------------------------------------------------------------*/ {
    /*
     *  write the geoloc SDR
     */
    if (wr_geo_scan(iscn, out_rec) != 0) return 1;
    /*
     *  Proceed to writing out the band files
     */
    if (wr_bnd_scan(iscn, out_rec) != 0) return 1;

    return 0;
}
