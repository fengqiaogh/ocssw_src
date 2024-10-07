/*******************************************************************

  ncio.c - a set of convenience netcdf functions

  (not all netcdf operations require a convenience function, but in 
   some cases, these are nice)

  ncio_dim_siz - for a dimension name, return a dimension size 
  ncio_grab_stdscl_ds - get an entire dataset and scale it using netcdf 
                        standards
 *******************************************************************/
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <netcdf.h>
#include <stdlib.h>

int ncio_dim_siz(int ncid, char *dim_nm)
/*******************************************************************

   ncio_dim_siz

   purpose: find the size of a named dimension

   Returns type: int - size of dimension or -1 if problem

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               ncid             I      netcdf id of file
      char *            dim_nm           I      name of dimension

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       1 Aug 2013     original development

 *******************************************************************/ {
    int dim_id, status, ret;
    size_t dim_len;

    status = nc_inq_dimid(ncid, dim_nm, &dim_id);
    if (status != NC_NOERR) {
        printf("%s, %d: nc_inq_dimid failure\n", __FILE__, __LINE__);
        ret = -1;
    } else {
        status = nc_inq_dimlen(ncid, dim_id, &dim_len);
        if (status != NC_NOERR) {
            printf("%s, %d: nc_inq_dimid failure\n", __FILE__, __LINE__);
            ret = -1;
        } else
            ret = (int) dim_len;
    }
    return ret;
}

int ncio_grab_f_ds(int ncid, char *ds_name, float *data)
/*******************************************************************

   ncio_grab_f_ds

   purpose: grab dataset and return it in float format

   Returns type: int -  0 if OK, -1 if problem

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               ncid             I      netcdf id of file
      char *            ds_name          I      name of dataset to read
      float *           data             O      pointer to a pre-allocated 
                                                array to receive the data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       1 Aug 2013     original development

 *******************************************************************/ {
    int status;
    int var_id;
    /*
     *  get ID of the dataset
     */
    if ((status = nc_inq_varid(ncid, ds_name, &var_id)) != NC_NOERR) {
        printf("%s, %d: nc_inq_varid returned error %d\n", __FILE__, __LINE__,
                status);
        return -1;
    }
    /*
     *  read the dataset as a float 
     */
    if ((status = nc_get_var_float(ncid, var_id, data)) != NC_NOERR) {
        printf("%s, %d: nc_get_var_float returned error %d\n", __FILE__, __LINE__,
                status);
        return -1;
    }
    return 0;
}

int ncio_grab_stdsclf_ds(int ncid, char *ds_name, float fill, float *data)
/*******************************************************************

   ncio_grab_stdsclf_ds

   purpose: grab dataset with standard scaling, floating point values
     read an entire dataset that is scaled and return as unscaled
     float

   This uses the netcdf standard dataset attributes scale_factor and add_offset
   to convert the scaled values in the dataset into real values and uses 
   attributes _FillValue and missing_value to identify non-values and replace
   them with a value you choose

   Returns type: int -  0 if OK, -1 if problem

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               ncid             I      netcdf id of file
      char *            ds_name          I      name of dataset to read
      float             fill             I      a value to put in the output 
                                                dataset for any fill or 
                                                missing values - preferably
                                                something not in normal data
                                                range
      float *           data             O      pointer to a pre-allocated 
                                                array to receive the data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       1 Aug 2013     original development

 *******************************************************************/ {
    int ndim, dim_ids[100], status, ntot, i;
    int var_id;
    size_t dim_len;
    float scale, offset, fillv, missv;
    char dim_nam[NC_MAX_NAME + 1];
    /*
     *  get ID of the dataset
     */
    if ((status = nc_inq_varid(ncid, ds_name, &var_id)) != NC_NOERR) {
        printf("%s, %d: nc_inq_varid returned error %d\n", __FILE__, __LINE__,
                status);
        return -1;
    }
    /*
     *  get all the needed attribute for the dataset
     */
    if ((status = nc_get_att_float(ncid, var_id, "scale_factor", &scale))
            != NC_NOERR) {
        printf("%s, %d: nc_get_att_float returned error %d\n", __FILE__, __LINE__,
                status);
        return -1;
    }
    if ((status = nc_get_att_float(ncid, var_id, "add_offset", &offset))
            != NC_NOERR) {
        printf("%s, %d: nc_get_att_float returned error %d\n", __FILE__, __LINE__,
                status);
        return -1;
    }
    if ((status = nc_get_att_float(ncid, var_id, "_FillValue", &fillv))
            != NC_NOERR) {
        printf("%s, %d: nc_get_att_float returned error %d\n", __FILE__, __LINE__,
                status);
        return -1;
    }
    status = nc_get_att_float(ncid, var_id, "missing_value", &missv);
    if ((status = nc_get_att_float(ncid, var_id, "missing_value", &missv))
            != NC_NOERR) {
        printf("%s, %d: nc_get_att_float returned error %d\n", __FILE__, __LINE__,
                status);
        return -1;
    }
    /*
     *  get the size for de-scaling use
     */
    if ((status = nc_inq_var(ncid, var_id, NULL, NULL, &ndim, dim_ids, NULL))
            != NC_NOERR) {
        printf("%s, %d: nc_inq_var returned error %d\n", __FILE__, __LINE__,
                status);
        return -1;
    }
    ntot = 1;
    for (i = 0; i < ndim; i++) {
        if ((status = nc_inq_dim(ncid, *(dim_ids + i), dim_nam, &dim_len))
                != NC_NOERR) {
            printf("%s, %d: nc_inq_dim returned error %d\n", __FILE__, __LINE__,
                    status);
            return -1;
        }
        ntot *= dim_len;
    }
    /*
     *  read the dataset as a float and convert the data values to unscaled values
     *  and fill, missing values to the desired fill
     */
    if ((status = nc_get_var_float(ncid, var_id, data)) != NC_NOERR) {
        printf("%s, %d: nc_get_var_float returned error %d\n", __FILE__, __LINE__,
                status);
        return -1;
    }
    for (i = 0; i < ntot; i++)
        *(data + i) = (*(data + i) == fillv) || (*(data + i) == missv) ?
        fill : *(data + i) * scale + offset;
    return 0;
}
