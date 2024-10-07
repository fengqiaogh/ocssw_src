#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <timeutils.h>

#include "hdf_bin.h"
#include <L3ShapeIsine.h>

#define VERSION "2.0"

using namespace std;

/*
  Revision 2.0  08/31/16
  Converted to 64 bit binner
  R. Healy

  Revision 1.03 04/05/12
  Use (int32_t *) cast for meta_l3b.smsec and meta_l3b.emsec
  J. Gales

  Revision 1.02 03/22/12
  Populate Start & End Time metadata
  J. Gales

  Revision 1.01 03/20/12
  Add support for products other that CHL (!)
  Clean up metadata
  J. Gales

  Revision 1.00 03/19/12
  Initial Version
  J. Gales
 */

void usage(const char *progname) {
    printf("%s %s (%s %s)\n", progname, VERSION, __DATE__, __TIME__);

    printf("\nUsage: %s input_meris_binfile output_hdf4_binfile meris_prodname hdf4_prodname\n",
            progname);

    exit(0);
}

int main(int argc, char **argv) {
    int ncid;
    int row_varid, col_varid, binid;
    int rowid, lon_varid, lat_varid, lon_step_varid;
    size_t nbins, nrows;

    short int *row, *col;
    float *lon, *lat, *lon_step;

    char ncbin_filename[256];

    /* Loop index. */
    int i;

    if (argc == 1) usage("l3rebin_meris");

    strcpy(ncbin_filename, argv[1]);
    printf("Reading product:  %s in %s\n", argv[3], ncbin_filename);

    /* Open the file. */
    nc_open(ncbin_filename, NC_NOWRITE, &ncid);

    nc_inq_dimid(ncid, "bin", &binid);
    nc_inq_dimlen(ncid, binid, &nbins);
    nc_inq_dimid(ncid, "row", &rowid);
    nc_inq_dimlen(ncid, rowid, &nrows);

    row = (short int *) calloc(nbins, sizeof (short int));
    col = (short int *) calloc(nbins, sizeof (short int));

    lon = (float *) calloc(nrows, sizeof (float));
    lat = (float *) calloc(nrows, sizeof (float));
    lon_step = (float *) calloc(nrows, sizeof (float));

    l3::L3ShapeIsine binShape(nrows);
    
    /* Read the row/col data. */
    nc_inq_varid(ncid, "row", &row_varid);
    nc_inq_varid(ncid, "col", &col_varid);
    nc_get_var_short(ncid, row_varid, row);
    nc_get_var_short(ncid, col_varid, col);

    nc_inq_varid(ncid, "center_lon", &lon_varid);
    nc_inq_varid(ncid, "center_lat", &lat_varid);
    nc_inq_varid(ncid, "lon_step", &lon_step_varid);
    nc_get_var_float(ncid, lon_varid, lon);
    nc_get_var_float(ncid, lat_varid, lat);
    nc_get_var_float(ncid, lon_step_varid, lon_step);

    int prod_varid, nobs_varid;
    float *prod = (float *) calloc(nbins, sizeof (float));
    short int *nobs = (short int *) calloc(nbins, sizeof (short int));

    int nvars;
    nc_inq_nvars(ncid, &nvars);

    string varname;
    string cntname;
    varname = argv[3];
    varname.append("_mean");
    cntname = argv[3];
    cntname.append("_count");

    nc_inq_varid(ncid, varname.c_str(), &prod_varid);
    nc_inq_varid(ncid, cntname.c_str(), &nobs_varid);

    nc_get_var_float(ncid, prod_varid, prod);
    nc_get_var_short(ncid, nobs_varid, nobs);

    Hdf::hdf_bin *output_binfile;
    output_binfile = new Hdf::hdf4_bin;
    output_binfile->hasNoext = true;
    printf("Creating product: %s in %s\n\n", argv[4], argv[2]);
    output_binfile->create(argv[2], 4320);
    output_binfile->clear_binlist();

    int64_t bin_num;
    int32_t irow, icol;
    float out_sum_buf[2 * 2 * 4320];
    int32_t last_row;
    int32_t tot_write = 0;
    int32_t row_off = row[0];
    for (i = 0; i < (int) nbins; i++) {
        if ((i % 1000000) == 0) cout << i << " out of " << nbins << endl;

        //p lat[row[0]]
        //p lon[row[0]] + lon_step[row[0]]*col[0]
        int32_t r = row[i] - row_off;
        binShape.latlon2rowcol(lat[r], lon[r] + lon_step[r] * col[i], irow, icol);
        bin_num = binShape.rowcol2bin(irow, icol);
        irow--;
        icol--;
        if (i == 0) last_row = irow;

        if (irow != last_row) {
            int32_t n_write = 0;
            int32_t kbin_max = output_binfile->get_numbin(last_row);
            for (int32_t kbin = 0; kbin <= kbin_max; kbin++) {

                if (output_binfile->get_bin_num(kbin) != 0) {

                    /* Remove "blank" bin records */
                    /* -------------------------- */
                    if (n_write != kbin)
                        memcpy(&out_sum_buf[2 * n_write], &out_sum_buf[2 * kbin], 8);

                    /* Remove "blank" bin records */
                    /* -------------------------- */
                    if (n_write != kbin)
                        output_binfile->copy_binlist(kbin, n_write);

                    n_write++;
                    tot_write++;
                } // bin_num != 0
            } /* kbin loop */

            /* Write BinList & Data Products */
            /* ----------------------------- */
            output_binfile->writeBinList(n_write);
            output_binfile->writeSums(&out_sum_buf[0], n_write, argv[4]);
            last_row = irow;

            output_binfile->clear_binlist();
        }

        int64_t basebin = output_binfile->get_basebin(irow);
        int32_t offset = bin_num - basebin;
        if (offset < 0) {
            cout << "Negative offset: " << offset <<
                    " for bin_num: " << bin_num << endl;
            exit(1);
        }
        //     cout << i << " " << bin_num << " " << offset << " " << irow << 
        //" " << lat[row[i]] << " " << lon[row[i]] + lon_step[row[i]]*col[i] << 
        //" " << row[i] << " " << col[i] << endl;

        output_binfile->set_bin_num(offset, bin_num);
        output_binfile->inc_nobs(offset, nobs[i]);
        output_binfile->inc_weights(offset, (float) nobs[i]);
        out_sum_buf[2 * icol] = prod[i] * nobs[i];
        out_sum_buf[2 * icol + 1] = 0.0;
    }


    // Fill metadata
    strcpy(output_binfile->meta_l3b.product_name, argv[2]);
    strcpy(output_binfile->meta_l3b.sensor_name, "MERIS");
    strcpy(output_binfile->meta_l3b.mission, "");
    strcpy(output_binfile->meta_l3b.mission_char, "");
    strcpy(output_binfile->meta_l3b.sensor, "");
    strcpy(output_binfile->meta_l3b.sensor_char, "");
    strcpy(output_binfile->meta_l3b.prod_type, "");
    strcpy(output_binfile->meta_l3b.pversion, "");
    strcpy(output_binfile->meta_l3b.soft_name, "l3rebin_meris");
    strcpy(output_binfile->meta_l3b.soft_ver, VERSION);
    output_binfile->meta_l3b.start_orb = 0;
    output_binfile->meta_l3b.end_orb = 0;
    strcpy(output_binfile->meta_l3b.ptime, "");
    strcpy(output_binfile->meta_l3b.proc_con, "");
    strcpy(output_binfile->meta_l3b.input_parms, "");
    strcpy(output_binfile->meta_l3b.infiles, "");
    strcpy(output_binfile->meta_l3b.flag_names, "");
    //  strcpy( output_binfile->meta_l3b.stime, "");
    //  strcpy( output_binfile->meta_l3b.etime, "");
    //  output_binfile->meta_l3b.bin_syear = 0;
    //  output_binfile->meta_l3b.bin_sday = 0;
    //  output_binfile->meta_l3b.bin_eyear = 0;
    //  output_binfile->meta_l3b.bin_eday = 0;

    char buf[256];
    string str;
    istringstream istr;
    int32_t yy, mm, dd, hh, mn, sc;
    int32_t year, day;
    int32_t msec;
    nc_get_att_text(ncid, NC_GLOBAL, "start_time", buf);
    str = buf;
    istr.clear();
    istr.str(str.substr(0, 4));
    istr >> yy;
    istr.clear();
    istr.str(str.substr(4, 2));
    istr >> mm;
    istr.clear();
    istr.str(str.substr(6, 2));
    istr >> dd;
    istr.clear();
    istr.str(str.substr(9, 2));
    istr >> hh;
    istr.clear();
    istr.str(str.substr(11, 2));
    istr >> mn;
    istr.clear();
    istr.str(str.substr(13, 2));
    istr >> sc;

    ymdhms2ydmsec(yy, mm, dd, hh, mn, sc, &year, &day, &msec);

    output_binfile->meta_l3b.startTime = yds2unix((int16_t) year, (int16_t) day, (double) msec);

    nc_get_att_text(ncid, NC_GLOBAL, "end_time", buf);
    str = buf;
    istr.clear();
    istr.str(str.substr(0, 4));
    istr >> yy;
    istr.clear();
    istr.str(str.substr(4, 2));
    istr >> mm;
    istr.clear();
    istr.str(str.substr(6, 2));
    istr >> dd;
    istr.clear();
    istr.str(str.substr(9, 2));
    istr >> hh;
    istr.clear();
    istr.str(str.substr(11, 2));
    istr >> mn;
    istr.clear();
    istr.str(str.substr(13, 2));
    istr >> sc;

    ymdhms2ydmsec(yy, mm, dd, hh, mn, sc, &year, &day, &msec);
    output_binfile->meta_l3b.startTime = yds2unix((int16_t) year, (int16_t) day, (double) msec);


    strcpy(output_binfile->meta_l3b.units, "");

    /* Close the file. */
    nc_close(ncid);
    output_binfile->close();

    return 0;
}




