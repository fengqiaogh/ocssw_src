/*
  Modification history:
  Programmer       Organization      Date      Description of change
  --------------   ------------    --------    ---------------------
  Joel Gales       Futuretech      08/10/03    Original Development

 */

/*
  Revision 0.76 2018-04-30
  Return status 110 when cannot extract specified area.
  Truncate area when requested beyond granule, rather than give up.
  G. Fireman

  Revision 0.75 2017-03-03
  Check for "empty" fields
  This avoids writing day-only EV fields
  G. Fireman

  Revision 0.70 2011-02-08
  Make sure all engineering data required by L1B has been collected
  Complete overhaul for maintainability
  G. Fireman

  2009-10-14
  Add values of any pre-existing pixel & scan offsets for correct geolocation
  G. Fireman

  Revision 0.62 01/25/07
  Remove clear page command in usage
  J. Gales

  Revision 0.61 05/30/06
  Emphasis Line # (not scan #) in usage.
  J. Gales
 */


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "hdf.h"
#include "mfhdf.h"

#define VERSION "0.76"
#define BOUNDS_ERROR 110

int main(int argc, char *argv[]) {
    int32 i, j;
    int32 status;

    /* pixel & scan  */
    int32 spixl, epixl, npixl;
    int32 sscan, escan, nscan, iscan;
    char *scan_type;
    int32 nscan_day = 0, nscan_night = 0;
    int32 spixl_old = 0, sscan_old = 0, nscan_old = 0;

    /* specific to input file */
    char *infile;
    int32 HDFfid_r;
    int32 sd_id_r;
    int32 sds_id_r;
    int32 vd_id_r;
    int32 dims_r[H4_MAX_VAR_DIMS];
    int32 start_r[H4_MAX_VAR_DIMS];
    int32 n_datasets_r;
    intn emptySDS;

    /* specific to output file */
    char *outfile;
    int32 HDFfid_w;
    int32 sd_id_w;
    int32 sds_id_w;
    int32 vd_id_w;
    int32 dims_w[H4_MAX_VAR_DIMS];
    int32 start_w[H4_MAX_VAR_DIMS];
    int32 n_datasets_w;

    /* attribute-specific */
    int32 nattrs;
    int32 nglobal_attr;
    static char attr_buf[2048];

    /* SDS-specific */
    static char sds_name[H4_MAX_NC_NAME];
    int32 sds_index;
    int32 ndims; /* same as rank */
    int32 dtype;
    int32 count;
    unsigned char *data;
    int32 nelem;

    /* VD-specific */
    static char vdata_name[VSNAMELENMAX];
    int32 vd_index;
    int32 n_vdata;
    int32 *ref_array;
    char fieldlist[VSFIELDMAX * FIELDNAMELENMAX];
    int32 vdata_size, n_records, interlace_mode;
    int32 n_fields;
    int32 nrec_r;
    int32 nrec_w;

    /* dimension info for expected SDSs */
    typedef struct {
        char name[H4_MAX_NC_NAME];
        int32 scandim;
        int32 scanfactor;
        int32 pixldim;
        int32 pixlfactor;
    } sds_str;

    sds_str sds_list[] = {
        {"Scan number", 0, 1, -1, -1},
        {"Frame count array", 0, 1, -1, -1},
        {"Scan Type", 0, 1, -1, -1},
        {"SD start time", 0, 1, -1, -1},
        {"SRCA start time", 0, 1, -1, -1},
        {"BB start time", 0, 1, -1, -1},
        {"SV start time", 0, 1, -1, -1},
        {"EV start time", 0, 1, -1, -1},
        {"SRCA calibration mode", 0, 1, -1, -1},
        {"Packet scan count", 0, 1, -1, -1},
        {"CCSDS Application Identifiers", 0, 1, -1, -1},
        {"Packet expedited data flag", 0, 1, -1, -1},
        {"Mirror side", 0, 1, -1, -1},
        {"Scan quality array", 0, 1, -1, -1},
        {"SD sector Pixel quality", 0, 1, -1, -1},
        {"SRCA sector Pixel quality", 0, 1, -1, -1},
        {"BB sector Pixel quality", 0, 1, -1, -1},
        {"SV sector Pixel quality", 0, 1, -1, -1},
        {"Earth sector Pixel quality", 0, 1, 1, 1},
        {"SD_250m", 0, 40, -1, -1},
        {"SD_500m", 0, 20, -1, -1},
        {"SD_1km_day", 0, 10, -1, -1},
        {"SD_1km_night", 0, 10, -1, -1},
        {"SRCA_250m", 0, 40, -1, -1},
        {"SRCA_500m", 0, 20, -1, -1},
        {"SRCA_1km_day", 0, 10, -1, -1},
        {"SRCA_1km_night", 0, 10, -1, -1},
        {"BB_250m", 0, 40, -1, -1},
        {"BB_500m", 0, 20, -1, -1},
        {"BB_1km_day", 0, 10, -1, -1},
        {"BB_1km_night", 0, 10, -1, -1},
        {"SV_250m", 0, 40, -1, -1},
        {"SV_500m", 0, 20, -1, -1},
        {"SV_1km_day", 0, 10, -1, -1},
        {"SV_1km_night", 0, 10, -1, -1},
        {"EV_250m", 0, 40, 2, 4},
        {"EV_500m", 0, 20, 2, 2},
        {"EV_1km_day", 0, 10, 2, 1},
        {"EV_1km_night", 0, 10, 2, 1},
        {"fpa_aem_config", 0, 1, -1, -1},
        {"science_state", 0, 1, -1, -1},
        {"science_abnormal", 0, 1, -1, -1},
        {"fpa_dcr_offset", 0, 1, -1, -1},
        {"raw_mir_enc", 0, 1, -1, -1},
        {"raw_vs_def", 0, 1, -1, -1},
        {"raw_vs_act", 0, 1, -1, -1},
        {"raw_sci_eng", 0, 1, -1, -1},
        {"raw_hk_telem", 0, 1, -1, -1},
        {"raw_sc_ancil", 0, 1, -1, -1},
        {"raw_param", 0, 1, -1, -1},
        {"raw_pv_gains", 0, 1, -1, -1},
    };
    n_datasets_w = sizeof (sds_list) / sizeof (sds_str);


    /* ------------------------------------------------------------------------ */

    printf("This is version %s of %s (compiled on %s %s)\n",
            VERSION, argv[0], __DATE__, __TIME__);

    /*** check usage ***/
    if (argc != 7) {
        printf("\nUsage: %s ", argv[0]);
        printf("infile spixl epixl sline eline outfile"
                "\n   where:"
                "\n\tinfile   - input MODIS L1A datafile"
                "\n\tspixl    - start pixel number (1-based)"
                "\n\tepixl    - end pixel number (1-based)"
                "\n\tsline    - start line (1-based)"
                "\n\teline    - end line (1-based)"
                "\n\toutfile  - output file name");
        printf("\n\nNote: Enter line number NOT scan number!\n");
        exit(EXIT_FAILURE);
    }

    /* load input parameters into local variables */
    infile = argv[1];
    spixl = atoi(argv[2]);
    epixl = atoi(argv[3]);
    sscan = atoi(argv[4]);
    escan = atoi(argv[5]);
    outfile = argv[6];

    printf("Input file:  %s\n", infile);
    printf("Output file: %s\n", outfile);

    /* Open HDF input file */
    HDFfid_r = Hopen(infile, DFACC_READ, 0);
    sd_id_r = SDstart(infile, DFACC_RDONLY);

    /* Convert from line to scan (line = 10*scan) */
    sscan = (sscan / 10) + (1 * ((sscan % 10) != 0));
    escan = (escan / 10) + (1 * ((escan % 10) != 0));

    if ((spixl > epixl) || (sscan > escan)) {
        printf("\nInvalid range requested:\n");
	npixl = epixl - spixl + 1;
	nscan = escan - sscan + 1;
	printf("spixl: %4d  epixl: %4d  npixl: %4d\n", spixl, epixl, npixl);
	printf("sscan: %4d  escan: %4d  nscan: %4d\n", sscan, escan, nscan);
        exit(BOUNDS_ERROR);
    }

    /* retrieve any previous pixel & scan offsets */
    if (SDfindattr(sd_id_r, "Extract Pixel Offset") != FAIL)
        SDreadattr(sd_id_r, SDfindattr(sd_id_r, "Extract Pixel Offset"), &spixl_old);
    if (SDfindattr(sd_id_r, "Extract Line Offset") != FAIL)
        SDreadattr(sd_id_r, SDfindattr(sd_id_r, "Extract Line Offset"), &sscan_old);

    /* retrieve dimension info */
    sds_id_r = SDselect(sd_id_r, SDnametoindex(sd_id_r, "EV_1km_day"));
    status = SDgetinfo(sds_id_r, sds_name, &ndims, dims_r, &dtype, &nattrs);
    status = SDreadattr(sd_id_r, SDfindattr(sd_id_r, "Number of Scans"), &nscan_old);

    /* To make sure that all engineering data required by L1B has been collected,
       scan number 7 (or later) of the original L1A file must be included. */
    if (nscan_old < 7) {
            printf("Only %d scans in granule: not enough engineering data to continue processing.", nscan_old);
        exit(BOUNDS_ERROR);
    }
    if ((sscan_old + escan) < 7) {
        printf("Adjusting end scan to ensure inclusion of required engineering data.\n");
        escan = 7 - sscan_old;
    }

    /* Check requested number of pixels & scan lines  */
    if (spixl > dims_r[2]) {
        printf("spixl: %d greater than # of pixels per scan: %d\n", spixl, dims_r[2]);
        exit(BOUNDS_ERROR);
    }
    if (epixl > dims_r[2]) {
        printf("epixl: %d greater than # of pixels per scan: %d; adjusting.\n",
               epixl, dims_r[2]);
        epixl = dims_r[2];
    }

    if (sscan > nscan_old) {
        printf("sscan: %d greater than # of scan lines: %d\n", sscan, nscan_old);
        exit(BOUNDS_ERROR);
    }
    if (escan > nscan_old) {
        printf("escan: %d greater than # of scan lines: %d; adjusting.\n",
               escan, nscan_old);
        escan = nscan_old;
    }

    npixl = epixl - spixl + 1;
    nscan = escan - sscan + 1;
    printf("spixl: %4d  epixl: %4d  npixl: %4d\n", spixl, epixl, npixl);
    printf("sscan: %4d  escan: %4d  nscan: %4d\n", sscan, escan, nscan);
    printf("Actual line start (0-based): %d\n", (sscan - 1) * 10);

    /* Create output HDF file & open SDS interface */
    HDFfid_w = Hopen(outfile, DFACC_CREATE, 0);
    sd_id_w = SDstart(outfile, DFACC_RDWR);

    /* ------------------------------------------------------------------------ */
    /* Write SDS subsets to output file */

    printf("\nWriting Science Data Sets...\n");

    /* Step through input datasets */
    SDfileinfo(sd_id_r, &n_datasets_r, &nglobal_attr);
    for (sds_index = 0; sds_index < n_datasets_r; sds_index++) {

        /* Get info for this SDS */
        sds_id_r = SDselect(sd_id_r, sds_index);
        status = SDgetinfo(sds_id_r, sds_name, &ndims, dims_r, &dtype, &nattrs);
        for (j = 0; j < n_datasets_w; j++) {
            if (strcmp(sds_list[j].name, sds_name) == 0) break;
        }
        if (j == n_datasets_w) {
            printf("Unexpected input SDS: \"%s\"\n", sds_name);
            exit(EXIT_FAILURE);
        }

        /* Find input & output dimension parameters */
        nelem = 1;
        for (i = 0; i < ndims; i++) {
            start_r[i] = 0;
            dims_w[i] = dims_r[i];
            if (i == sds_list[j].scandim) {
                start_r[i] = sds_list[j].scanfactor * (sscan - 1);
                dims_w[i] = sds_list[j].scanfactor * nscan;
            }
            if (i == sds_list[j].pixldim) {
                start_r[i] = sds_list[j].pixlfactor * (spixl - 1);
                dims_w[i] = sds_list[j].pixlfactor * npixl;
            }
            nelem *= dims_w[i];
        }

        /* Print info */
        printf("%s\n    rank %d, type %d, dims", sds_name, ndims, dtype);
        for (i = 0; i < ndims; i++) printf("[%d]", dims_r[i]);
        printf("; extracting ");
        for (i = 0; i < ndims; i++) printf("[%d-%d]", start_r[i], start_r[i] + dims_w[i] - 1);
        printf("\n");

        /* Create SDS */
        sds_id_w = SDcreate(sd_id_w, sds_name, dtype, ndims, dims_w);
        if (sds_id_w == -1) {
            printf("Field: %s cannot be created\n", sds_name);
            exit(EXIT_FAILURE);
        }

        /* Read and set SDS attributes */
        for (i = 0; i < nattrs; i++) {
            status = SDreadattr(sds_id_r, i, attr_buf);
            status = SDattrinfo(sds_id_r, i, sds_name, &dtype, &count);
            status = SDsetattr(sds_id_w, sds_name, dtype, count, attr_buf);
        }

        /* Read and write SDS values */
        data = calloc(nelem, DFKNTsize(dtype));
        status = SDreaddata(sds_id_r, start_r, NULL, dims_w, data);
        status = SDcheckempty(sds_id_r, &emptySDS);
        if (emptySDS == 0)
            status = SDwritedata(sds_id_w, start_w, NULL, dims_w, data);
        if (status == -1) {
            printf("write status: %d\n\n", status);
            exit(EXIT_FAILURE);
        }

        /* Count Day and Night Scans (for global attributes) */
        if (strcmp(sds_name, "Scan Type") == 0) {
            for (iscan = 0; iscan < nscan; iscan++) {
                scan_type = (char *) (data + 10 * iscan);
                /* printf("scan_type[%ld] = %s\n", iscan, scan_type); */
                if (strcmp(scan_type, "Day") == 0) nscan_day++;
                if (strcmp(scan_type, "Night") == 0) nscan_night++;
            }
        }

        free(data);
        SDendaccess(sds_id_w);
        SDendaccess(sds_id_r);
    }

    /* ------------------------------------------------------------------------ */
    /* Write Vdata to output file */

    printf("\nWriting VData...\n");
    status = Vstart(HDFfid_r);
    status = Vstart(HDFfid_w);

    /* Step through each lone vdata */
    n_vdata = VSlone(HDFfid_r, NULL, 0);
    ref_array = calloc(n_vdata, sizeof (int32));
    n_vdata = VSlone(HDFfid_r, ref_array, n_vdata);
    for (vd_index = 0; vd_index < n_vdata; vd_index++) {

        /* Get info for input VData */
        vd_id_r = VSattach(HDFfid_r, ref_array[vd_index], "r");
        status = VSinquire(vd_id_r, &n_records, &interlace_mode,
                fieldlist, &vdata_size, vdata_name);

        /* Set info for output VData */
        vd_id_w = VSattach(HDFfid_w, -1, "w");
        VSsetname(vd_id_w, vdata_name);
        VSsetinterlace(vd_id_w, interlace_mode);

        /* Define output fields */
        n_fields = VFnfields(vd_id_r);
        printf("%s\n    n_records=%d, n_fields=%d\n", vdata_name, n_records, n_fields);
        for (i = 0; i < n_fields; i++) {
            VSfdefine(vd_id_w,
                    VFfieldname(vd_id_r, i),
                    VFfieldtype(vd_id_r, i),
                    VFfieldorder(vd_id_r, i));
        }
        status = VSsetfields(vd_id_w, fieldlist);

        /* Set pointer to read only desired scans */
        /* (except for temperatures, which don't update often enough) */
        status = VSsetfields(vd_id_r, fieldlist);
        if ((n_records == nscan_old) &&
                strcmp(vdata_name, "Telemetry Major Cycle 0 of 63") != 0) {
            VSseek(vd_id_r, sscan - 1);
            n_records = nscan;
        }

        /* Read and write VData values */
        data = calloc(n_records * VSsizeof(vd_id_r, fieldlist), 1);
        nrec_r = VSread(vd_id_r, (unsigned char *) data, n_records, interlace_mode);
        nrec_w = VSwrite(vd_id_w, (unsigned char *) data, n_records, interlace_mode);
        if (nrec_r != nrec_w) {
            printf("Write Error: %d %d %d\n", i, nrec_r, nrec_w);
        }

        free(data);
        VSdetach(vd_id_w);
        VSdetach(vd_id_r);
    }

    free(ref_array);
    Vend(HDFfid_w);
    Vend(HDFfid_r);


    /* ------------------------------------------------------------------------ */
    /* Write Global Attributes to output file */

    printf("\nWriting Global Attributes...\n");
    for (i = 0; i < nglobal_attr; i++) {
        status = SDattrinfo(sd_id_r, i, attr_buf, &dtype, &count);
        data = calloc(count, DFKNTsize(dtype));
        status = SDreadattr(sd_id_r, i, data);

        if (strcmp(attr_buf, "Number of Scans") == 0) memcpy(data, &nscan, 4);
        if (strcmp(attr_buf, "Number of Day mode scans") == 0) memcpy(data, &nscan_day, 4);
        if (strcmp(attr_buf, "Number of Night mode scans") == 0) memcpy(data, &nscan_night, 4);

        status = SDsetattr(sd_id_w, attr_buf, dtype, count, data);
        free(data);
    }

    /* update pixel offset & count */
    spixl = spixl_old + spixl - 1;
    status = SDsetattr(sd_id_w, "Extract Pixel Offset", DFNT_INT32, 1, &spixl);
    status = SDsetattr(sd_id_w, "Extract Pixel Count", DFNT_INT32, 1, &npixl);

    /* update scan offset & count */
    sscan = sscan_old + sscan - 1;
    status = SDsetattr(sd_id_w, "Extract Line Offset", DFNT_INT32, 1, &sscan);
    status = SDsetattr(sd_id_w, "Extract Line Count", DFNT_INT32, 1, &nscan);


    /* ------------------------------------------------------------------------ */
    /* Close input and output files */

    SDend(sd_id_w);
    SDend(sd_id_r);
    Hclose(HDFfid_w);
    Hclose(HDFfid_r);

    return EXIT_SUCCESS;
}
