#include "l1czcs.h"
#include "hdf.h"
#include "mfhdf.h"

#define csds_fail( prog, sds ) \
  {   \
  fprintf(stderr, "%s: create_sds failed on %s\n", prog, sds ); \
  return FAIL;  \
  }
#define setdim_fail( prog, sds ) \
  {   \
  fprintf(stderr, "%s: set_dim_names failed on %s\n", prog, sds ); \
  return FAIL;  \
  }
#define setatt_fail( prog, sds, atr_name )  \
  {  \
  fprintf(stderr, "%s: SDsetattr failed on SDS %s, attrib name %s\n", prog, sds, atr_name );  \
  return FAIL;  \
  }

int32 set_czcs_ctl_data(int32 sdfid, int32 vid, gattr_struc gattr,
        l1_data_struc l1_data)
/*-----------------------------------------------------------------------------
    Function:  set_czcs_ctl_data 

    Returns:   int32 (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        Write the nav and control points data to CZCS HDF file. 

    Parameters: (in calling order)
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int32     sdfid         I   SD interface ID of the output file
        int32     vid           I   group ID the SD belongs to
        gattr_struc  gattr      I   attribute structure
        l1_data_struc  l1_data  I   L1 data structure with sds data
    
    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  --------------------- 
    Gary Fu        SAIC GSC       05/13/98  Original development
    W. Robinson, SAIC             26May04   add other data for output and
                                            modify call structure
----------------------------------------------------------------------------*/ {
    int32 create_sds(int32, char *, int32, int32, int32 *, int32, VOIDP *);
    int32 set_dim_names(int32, char *, char *, char *);
    int32 dims[3] = {0, 0, 0};
    int sdsid;
    int rank, iret;
    char *pname = "set_czcs_ctl_data";
    char sds_name[128], sds_attr[128];
    char *num_scn_ln = "Number of Scan Lines";
    char *num_pix_ctl = "Number of Pixel Control Points";
    char *attr_nm_long = "long_name";
    char *attr_nm_valid = "valid range";
    char *attr_nm_units = "units";

    rank = 2;
    dims[0] = gattr.scan_lines;
    dims[1] = gattr.n_ctl_pt;
    strcpy(sds_name, "latitude");
    if ((sdsid = create_sds(sdfid, sds_name, DFNT_FLOAT32, rank, dims, vid,
            (VOIDP) l1_data.ctl_pt_lat)) < 0)
        csds_fail(pname, sds_name);

    if ((iret = set_dim_names(sdsid, num_scn_ln, num_pix_ctl, NULL)) < 0)
        setdim_fail(pname, sds_name);

    strcpy(sds_attr, "Latitude");
    if ((iret = SDsetattr(sdsid, attr_nm_long, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_long);

    strcpy(sds_attr, "(-90.,90.)");
    if ((iret = SDsetattr(sdsid, attr_nm_valid, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_valid);

    strcpy(sds_attr, "degrees");
    if ((iret = SDsetattr(sdsid, attr_nm_units, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_units);

    /*****************************************/
    strcpy(sds_name, "longitude");
    if ((sdsid = create_sds(sdfid, sds_name, DFNT_FLOAT32, rank, dims, vid,
            (VOIDP) l1_data.ctl_pt_lon)) < 0)
        csds_fail(pname, sds_name);

    if ((iret = set_dim_names(sdsid, num_scn_ln, num_pix_ctl, NULL)) < 0)
        setdim_fail(pname, sds_name);

    strcpy(sds_attr, "Longitude");
    if ((iret = SDsetattr(sdsid, attr_nm_long, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_long);

    strcpy(sds_attr, "(-180.,180.)");
    if ((iret = SDsetattr(sdsid, attr_nm_valid, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_valid);

    strcpy(sds_attr, "degrees");
    if ((iret = SDsetattr(sdsid, attr_nm_units, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_units);

    /*****************************************/
    rank = 1;
    dims[0] = gattr.scan_lines;

    strcpy(sds_name, "cntl_pt_rows");
    if ((sdsid = create_sds(sdfid, sds_name, DFNT_INT32, rank, dims, vid,
            (VOIDP) l1_data.ctl_pt_rows)) < 0)
        csds_fail(pname, sds_name);

    if ((iret = set_dim_names(sdsid, num_scn_ln, NULL, NULL)) < 0)
        setdim_fail(pname, sds_name);

    strcpy(sds_attr, "Scan Control Points");
    if ((iret = SDsetattr(sdsid, attr_nm_long, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_long);

    strcpy(sds_attr, "none");
    if ((iret = SDsetattr(sdsid, attr_nm_units, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_units);

    /*****************************************/
    dims[0] = gattr.n_ctl_pt;

    strcpy(sds_name, "cntl_pt_cols");
    if ((sdsid = create_sds(sdfid, sds_name, DFNT_INT32, rank, dims, vid,
            (VOIDP) l1_data.ctl_pt_cols)) < 0)
        csds_fail(pname, sds_name);

    if ((iret = set_dim_names(sdsid, num_pix_ctl, NULL, NULL)) < 0)
        setdim_fail(pname, sds_name);

    strcpy(sds_attr, "Control Points for Columns");
    if ((iret = SDsetattr(sdsid, attr_nm_long, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_long);

    strcpy(sds_attr, "none");
    if ((iret = SDsetattr(sdsid, attr_nm_units, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_units);
    /*
     *  orb_vec
     */
    rank = 2;
    dims[0] = gattr.scan_lines;
    dims[1] = 3;

    strcpy(sds_name, "orb_vec");
    if ((sdsid = create_sds(sdfid, sds_name, DFNT_FLOAT32, rank, dims, vid,
            (VOIDP) l1_data.orb_vec)) < 0)
        csds_fail(pname, sds_name);

    if ((iret = set_dim_names(sdsid, num_scn_ln, "3", NULL)) < 0)
        setdim_fail(pname, sds_name);

    strcpy(sds_attr, "Orbit position vector at scan line time");
    if ((iret = SDsetattr(sdsid, attr_nm_long, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_long);

    strcpy(sds_attr, "(-7200.,7200.)");
    if ((iret = SDsetattr(sdsid, attr_nm_valid, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_valid);

    strcpy(sds_attr, "kilometers");
    if ((iret = SDsetattr(sdsid, attr_nm_units, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_units);
    /*
     *  gain
     */
    rank = 1;
    dims[0] = gattr.scan_lines;

    strcpy(sds_name, "gain");
    if ((sdsid = create_sds(sdfid, sds_name, DFNT_INT16, rank, dims, vid,
            (VOIDP) l1_data.gain)) < 0)
        csds_fail(pname, sds_name);

    if ((iret = set_dim_names(sdsid, num_scn_ln, NULL, NULL)) < 0)
        setdim_fail(pname, sds_name);

    strcpy(sds_attr, "Gain setting at scan line time");
    if ((iret = SDsetattr(sdsid, attr_nm_long, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_long);

    strcpy(sds_attr, "( 1, 4 )");
    if ((iret = SDsetattr(sdsid, attr_nm_valid, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_valid);

    strcpy(sds_attr, "none");
    if ((iret = SDsetattr(sdsid, attr_nm_units, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_units);
    /*
     * slope and intercept
     */
    rank = 2;
    dims[0] = gattr.scan_lines;
    dims[1] = 6;

    strcpy(sds_name, "slope");
    if ((sdsid = create_sds(sdfid, sds_name, DFNT_FLOAT32, rank, dims, vid,
            (VOIDP) l1_data.slope)) < 0)
        csds_fail(pname, sds_name);

    if ((iret = set_dim_names(sdsid, num_scn_ln, "6", NULL)) < 0)
        setdim_fail(pname, sds_name);

    strcpy(sds_attr, "Calibration slope at scan line time");
    if ((iret = SDsetattr(sdsid, attr_nm_long, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_long);

    strcpy(sds_attr, "( -20., 20. )");
    if ((iret = SDsetattr(sdsid, attr_nm_valid, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_valid);

    strcpy(sds_attr, "mW cm^-2 um^-1 sr^-1 count^-1");
    if ((iret = SDsetattr(sdsid, attr_nm_units, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_units);
    /*   */
    strcpy(sds_name, "intercept");
    if ((sdsid = create_sds(sdfid, sds_name, DFNT_FLOAT32, rank, dims, vid,
            (VOIDP) l1_data.intercept)) < 0)
        csds_fail(pname, sds_name);

    if ((iret = set_dim_names(sdsid, num_scn_ln, "6", NULL)) < 0)
        setdim_fail(pname, sds_name);

    strcpy(sds_attr, "Calibration intercept at scan line time");
    if ((iret = SDsetattr(sdsid, attr_nm_long, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_long);

    strcpy(sds_attr, "( -20., 20. )");
    if ((iret = SDsetattr(sdsid, attr_nm_valid, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_valid);

    strcpy(sds_attr, "mW cm^-2 um^-1 sr^-1");
    if ((iret = SDsetattr(sdsid, attr_nm_units, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_units);
    /*
     * att_ang - attitude (yaw, roll, pitch )
     */
    dims[1] = 3;

    strcpy(sds_name, "att_ang");
    if ((sdsid = create_sds(sdfid, sds_name, DFNT_FLOAT32, rank, dims, vid,
            (VOIDP) l1_data.att_ang)) < 0)
        csds_fail(pname, sds_name);

    if ((iret = set_dim_names(sdsid, num_scn_ln, "3", NULL)) < 0)
        setdim_fail(pname, sds_name);

    strcpy(sds_attr, "Computed yaw, roll, pitch at scan line time");
    if ((iret = SDsetattr(sdsid, attr_nm_long, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_long);

    strcpy(sds_attr, "( -180., 180. )");
    if ((iret = SDsetattr(sdsid, attr_nm_valid, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_valid);

    strcpy(sds_attr, "degrees");
    if ((iret = SDsetattr(sdsid, attr_nm_units, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_units);

    /*
     * pos_err - position error in orb_vec orbit, due to interpolation error
     */
    rank = 1;

    strcpy(sds_name, "pos_err");
    if ((sdsid = create_sds(sdfid, sds_name, DFNT_FLOAT32, rank, dims, vid,
            (VOIDP) l1_data.pos_err)) < 0)
        csds_fail(pname, sds_name);

    if ((iret = set_dim_names(sdsid, num_scn_ln, NULL, NULL)) < 0)
        setdim_fail(pname, sds_name);

    strcpy(sds_attr, "Orbit position error");
    if ((iret = SDsetattr(sdsid, attr_nm_long, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_long);

    strcpy(sds_attr, "( -10., 7200.)");
    if ((iret = SDsetattr(sdsid, attr_nm_valid, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_valid);

    strcpy(sds_attr, "meters");
    if ((iret = SDsetattr(sdsid, attr_nm_units, DFNT_CHAR, strlen(sds_attr) + 1,
            (VOIDP) sds_attr)) < 0)
        setatt_fail(pname, sds_name, attr_nm_units);

    return SUCCEED;
}

