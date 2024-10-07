#include "h5io.h"
#include <stdlib.h>

int h5io_openr(char *file, int opt, h5io_str *id)
/*******************************************************************

   h5io_openr

   purpose: open a general hdf 5 file for reading (option to write in 
     existing file too)

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      name of file to open
      int32             opt              I      open option: 0 for standard
                                                open with no write, 1 to 
                                                open to read and write 
      h5io_str *        id               O      id for the opened file

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 23 Sep 2008     Original development

 *******************************************************************/ {
    unsigned int open_opt, opt_sel[2] = {H5F_ACC_RDONLY, H5F_ACC_RDWR};
    htri_t h5ret;
    H5E_auto_t old_func;
    void *old_client_data;
    /*
     *  for the read set proper open value
     */
    if ((opt != 0) && (opt != 1)) {
        printf("%s: Inproper open option selected\n", __FILE__);
        return 1;
    }
    /*
     *  check if file is HDF 5.  The H5E calls are to disable / ensble ugly error
     *  reporting if file is not HDF 5 (I know it should not for this, but it does
     */
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    h5ret = H5Fis_hdf5(file);

    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);
    if (h5ret <= 0) {
        return 1;
    }
    open_opt = opt_sel[ opt ];
    if ((id->file_id = H5Fopen(file, open_opt, H5P_DEFAULT))
            < 0) return 1;
    /*
     *  Also, find the root group associated with the file
     */
    if ((id->grp_id = H5Gopen1(id->file_id, "/")) < 0) return 1;
    /*
     *  file is open
     */
    id->type = H5IO_TYP_FILE_ID;
    return 0;
}

int h5io_openw(char *file, int opt, h5io_str *id)
/*******************************************************************

   h5io_openw

   purpose: open a general hdf 5 file for write

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      name of file to open
      int32             opt              I      open option: 0 for standard
                                                open, 1 to overwrite any old
                                                file by that name (dangerous)
      h5io_str *        id               O      id for the opened file

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 23 Sep 2008     Original development

 *******************************************************************/ {
    unsigned int open_opt, opt_sel[2] = {H5F_ACC_EXCL, H5F_ACC_TRUNC};
    /*
     *  for the read set proper open value
     */
    if ((opt != 0) && (opt != 1)) {
        printf("%s: Improper open option selected\n", __FILE__);
        return 1;
    }
    open_opt = opt_sel[ opt ];
    if ((id->file_id = H5Fcreate(file, open_opt, H5P_DEFAULT, H5P_DEFAULT))
            < 0) return 1;
    /*
     *  Also, find the root group associated with the file
     */
    if ((id->grp_id = H5Gopen1(id->file_id, "/")) < 0) return 1;
    /*
     *  file is open
     */
    id->type = H5IO_TYP_FILE_ID;
    return 0;
}

int h5io_close(h5io_str *id)
/*******************************************************************

   h5io_close

   purpose: close a h5io id for accessing the file, group, or dataset
     in the file, this frees the storage associated with the id 

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        id               I      id to close

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 23 Sep 2008     Original development

 *******************************************************************/ {
    /*
     *  based on the type and the read mode, clear storage of the 
     *  HDF 5 ID values
     */
    switch (id->type) {
    case H5IO_TYP_FILE_ID:
        if (H5Fclose(id->file_id) < 0) return 1;
        if (H5Gclose(id->grp_id) < 0) return 1;
        break;
    case H5IO_TYP_GRP_ID:
        if (H5Gclose(id->grp_id) < 0) return 1;
        break;
    case H5IO_TYP_DAT_ID:
        if (H5Dclose(id->dat_id) < 0) return 1;
        /*
         *  if the slice read / write mode was used, also close the data 
         *  space ids and indicate dataset is not open for slice I/O
         */
        if (id->dat_rw_mode == H5IO_DAT_RWSLICE) {
            if (H5Sclose(id->fil_space_id) < 0) return 1;
            if (H5Sclose(id->mem_space_id) < 0) return 1;
            if (H5Tclose(id->native_typ) < 0) return 1;
            id->dat_rw_mode = H5IO_DAT_RWALL; /* just in case, say not ready */
            /* for slices */
        }
        break;
    default:
        printf("%s: Unknown or un-initialized id\n", __FILE__);
        return 1;
    }
    /*
     *  and end
     */
    return 0;
}

int h5io_info(h5io_str *id, char *attr_name, H5T_class_t *class,
        hid_t *native_typ, int* ndim, int *dim_siz, int *sto_len)
/*******************************************************************

   h5io_info

   purpose: Get information about a hdf 5 attribute or dataset from the 
     file, group (for attribute), or dataset (for attribute or dataset)

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        id               I      id for the opened file,
                                                group, or dataset
      char *            attr_name        I      name of attribute to get
                                                If NULL, id must be a dataset
                                                and dataset info will be found
      H5T_class_t *     class            O      class found, see below
      hid_t *           native_typ       O      the full type of the variable
      int *             ndim             O      # dimensions of the array
      int *             dim_siz          O      length of each dimension
      int *             sto_len          O      storage length in bytes (use
                                                to get string length needed)

   Notes
     HDF classes ( class )
       These are the classes or general types of data that can be handled
       by this package (self descriptive)
     H5T_INTEGER, H5T_FLOAT, H5T_STRING 

     HDF native type
       The native type is more descriptive of the data including not only
       the class, but the size (ie int with 1, 2, 4 bytes)
       Also, to comlpetely free the storage associated with this, you should 
       do a H5Tclose( native_typ )

     The common native types are:
     Code                 Corresponding C Type
     -------------------  --------------------
     H5T_NATIVE_CHAR      char
     H5T_NATIVE_SCHAR     signed char
     H5T_NATIVE_UCHAR     unsigned char
     H5T_NATIVE_SHORT     short
     H5T_NATIVE_USHORT    unsigned short
     H5T_NATIVE_INT       int
     H5T_NATIVE_UINT      unsigned
     H5T_NATIVE_LONG      long
     H5T_NATIVE_ULONG     unsigned long
     H5T_NATIVE_LLONG     long long
     H5T_NATIVE_ULLONG    unsigned long long
     H5T_NATIVE_FLOAT     float
     H5T_NATIVE_DOUBLE    double
     H5T_NATIVE_LDOUBLE   long double

     all information on the types and classes are in the  HDF5 Datatypes 
     web page (version 1.6)
     http://www.hdfgroup.org/HDF5/doc1.6/UG/UG_frame.html.

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 30 Sep 2008     Original development
      W. Robinson, SAIC 10 Feb 2009     Change test for null attr name to
                                        if " == NULL

 *******************************************************************/ {
    hid_t *base_grp, spaceid, attr_datyp, attr_id;
    hsize_t dims[H5IO_MAXDIM], maxdims[H5IO_MAXDIM];
    size_t lsto_len;
    int indx;
    /*
     * set the proper ID and make sure attribute name is non null for a group
     */
    base_grp = (id->type == H5IO_TYP_DAT_ID) ?
            &(id->dat_id) : &(id->grp_id);
    if ((id->type != H5IO_TYP_DAT_ID) && (attr_name == NULL)) {
        printf("%s: no attribute entered for group information query\n",
                __FILE__);
        return 1;
    }
    if ((id->type == H5IO_TYP_DAT_ID) && (attr_name == NULL)) {
        /*
         *  dataset and no attribute name - get dataset specific information
         */
        if ((spaceid = H5Dget_space(*base_grp)) < 0) return 1;
        if ((attr_datyp = H5Dget_type(*base_grp)) < 0) return 1;
    } else {
        /*
         *  We have an attribute to open and get specific information
         */
        if ((attr_id = H5Aopen_name(*base_grp, attr_name)) < 0) return 1;
        if ((spaceid = H5Aget_space(attr_id)) < 0) return 1;
        if ((attr_datyp = H5Aget_type(attr_id)) < 0) return 1;
    }
    /*
     *  get the data size, class (basic type), and naitve data type
     */
    if ((*ndim = H5Sget_simple_extent_dims(spaceid, dims, maxdims)) < 0)
        return 1;
    if (*ndim > H5IO_MAXDIM) {
        printf("%s: # dimensions of attribute: %s is > software maximum of %d\n",
                __FILE__, attr_name, H5IO_MAXDIM);
        return 1;
    }
    for (indx = 0; indx < *ndim; indx++)
        dim_siz[indx] = dims[indx];
    if ((*native_typ = H5Tget_native_type(attr_datyp, H5T_DIR_ASCEND)) < 0)
        return 1;
    if ((*class = H5Tget_class(attr_datyp)) == -1) return 1;
    if ((lsto_len = H5Tget_size(attr_datyp)) <= 0) return 1;
    *sto_len = lsto_len;
    /*  that's all there is to know, close the stuff  */

    /*
     *  based on the class, read the data in
     */
    switch (*class) {
    case H5T_STRING:
    case H5T_INTEGER:
    case H5T_FLOAT:
        break;
    case H5T_TIME:
    case H5T_BITFIELD:
    case H5T_OPAQUE:
    case H5T_COMPOUND:
    case H5T_REFERENCE:
    case H5T_ENUM:
    case H5T_VLEN:
    case H5T_ARRAY:
        printf("%s: Currently, classes beyond STRING, INTEGER, and FLOAT \n",
                "may not be fully characterized here\n");
        break;
    default:
        printf("%s: Unable to handle the class for attribute: %s\n",
                __FILE__, attr_name);
        return 1;
        break;
    }
    /*
     *  close the allocated spaces before leaving
     */
    if (H5Sclose(spaceid) < 0) return 1;
    if (H5Tclose(attr_datyp) < 0) return 1;
    if (attr_name != NULL)
        if (H5Aclose(attr_id) < 0) return 1;
    return 0;
}

int h5io_set_ds(h5io_str *id, char *path_name, h5io_str *ds_id)
/*******************************************************************

   h5io_rd_attr

   purpose: get the access to a particular dataset

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        id               I      id for the opened file,
                                                or group
      char *            path_name        I      path to set from current 
                                                location with the dataset name
                                                as the last part of the path
      h5io_str *        ds_id            O      id of the dataset found

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 25 Sep 2008     Original development

 *******************************************************************/ {
    /*
     *  make sure the incoming id is OK
     */
    if (id->type == H5IO_TYP_DAT_ID) {
        printf("%s: Cannot set to a dataset under a dataset\n", __FILE__);
        return 1;
    }
    /*
     *  just find the id of the dataset
     */
    if ((ds_id->dat_id = H5Dopen1(id->grp_id, path_name)) < 0) return 1;
    /*
     *  and set type for a group id and default to entire data array read
     */
    ds_id->type = H5IO_TYP_DAT_ID;
    ds_id->dat_rw_mode = H5IO_DAT_RWALL;
    return 0;
}

int h5io_set_grp(h5io_str *id, char *path_name, h5io_str *grp_id)
/*******************************************************************

   h5io_rd_attr

   purpose: get the access to a particular group

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        id               I      id for the opened file,
                                                or group
      char *            path_name        I      path to set from current 
                                                location
      h5io_str *        grp_id           O      id of the group found

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 25 Sep 2008     Original development

 *******************************************************************/ {
    /*
     *  make sure the incoming id is OK
     */
    if (id->type == H5IO_TYP_DAT_ID) {
        printf("%s: Cannot set group under a dataset\n", __FILE__);
        return 1;
    }
    /*
     *  just find the id of the group
     */
    if ((grp_id->grp_id = H5Gopen1(id->grp_id, path_name)) < 0) return 1;
    /*
     *  and set type for a group id
     */
    grp_id->type = H5IO_TYP_GRP_ID;
    return 0;
}

int h5io_rd_attr(h5io_str *id, char *attr_name, void *data)
/*******************************************************************

   h5io_rd_attr

   purpose: read a hdf 5 attribute from the file, group, or dataset

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        id               I      id for the opened file,
                                                group, or dataset
      char *            attr_name        I      name of attribute to get
      void *            data             O      pointer to space made for 
                                                the data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 23 Sep 2008     Original development

 *******************************************************************/ {
    hid_t *base_grp, spaceid, attr_datyp, attr_id, native_typ;
    hsize_t dims[H5IO_MAXDIM], maxdims[H5IO_MAXDIM];
    H5T_class_t d_class;
    int ndim, nval, idim;
    /*
     *  open the attribute
     */
    base_grp = (id->type == H5IO_TYP_DAT_ID) ?
            &(id->dat_id) : &(id->grp_id);
    if ((attr_id = H5Aopen_name(*base_grp, attr_name)) < 0) return 1;
    /*
     *  get the data size, class (basic type), and naitve data type
     */
    if ((spaceid = H5Aget_space(attr_id)) < 0) return 1;
    if ((ndim = H5Sget_simple_extent_dims(spaceid, dims, maxdims)) < 0)
        return 1;
    if (ndim > H5IO_MAXDIM) {
        printf("%s: # dimensions of attribute: %s is > software maximum of %d\n",
                __FILE__, attr_name, H5IO_MAXDIM);
        return 1;
    }
    if ((attr_datyp = H5Aget_type(attr_id)) < 0) return 1;
    if ((d_class = H5Tget_class(attr_datyp)) == -1) return 1;
    if ((native_typ = H5Tget_native_type(attr_datyp, H5T_DIR_ASCEND)) < 0)
        return 1;

    nval = 1;
    idim = 0;
    if (ndim > 0) {
        while (idim < ndim)
            nval *= dims[idim++];
    }
    /*
     *  based on the class, read the data in
     */
    switch (d_class) {
    case H5T_STRING:
        if (H5Aread(attr_id, native_typ, data) < 0) return 1;
        hsize_t length = H5Aget_storage_size(attr_id);
        ((char *) data)[length] = 0;
        break;
    case H5T_INTEGER:
    case H5T_FLOAT:
        if (H5Aread(attr_id, native_typ, data) < 0) return 1;
        break;
    case H5T_TIME:
    case H5T_BITFIELD:
    case H5T_OPAQUE:
    case H5T_COMPOUND:
    case H5T_REFERENCE:
    case H5T_ENUM:
    case H5T_VLEN:
    case H5T_ARRAY:
    default:
        printf("Unable to handle the class for attribute: %s\n", attr_name);
        return 1;
        break;
    }
    /*
     *  close the allocated spaces before leaving
     */
    if (H5Sclose(spaceid) < 0) return 1;
    if (H5Aclose(attr_id) < 0) return 1;
    if (H5Tclose(attr_datyp) < 0) return 1;
    if (H5Tclose(native_typ) < 0) return 1;
    return 0;
}

int h5io_rd_ds(h5io_str *ds_id, void *data)
/*******************************************************************

   h5io_rd_ds

   purpose: read an entire dataset into an array

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        ds_id            I      id for the dataset in the file
      void *            data             O      allocated space for the data 
                                                to be read into

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 25 Sep 2008     Original development

 *******************************************************************/ {
    hid_t spaceid, ds_datyp, native_typ;
    hsize_t dims[H5IO_MAXDIM], maxdims[H5IO_MAXDIM];
    H5T_class_t d_class;
    int ndim, nval, idim;
    /*
     *  make sure the incoming id is OK - it must be s dataset
     */
    if (ds_id->type != H5IO_TYP_DAT_ID) {
        printf("%s: A dataset id must be passed in for this routine\n", __FILE__);
        return 1;
    }
    /*
     *  The code here is very similar to the h5io_rd_attr code
     *  get the data size, class (basic type), and naitve data type
     */
    if ((spaceid = H5Dget_space(ds_id->dat_id)) < 0) return 1;
    if ((ndim = H5Sget_simple_extent_dims(spaceid, dims, maxdims)) < 0)
        return 1;
    if (ndim > H5IO_MAXDIM) {
        printf("%s: # dimensions of dataset is > software maximum of %d\n",
                __FILE__, H5IO_MAXDIM);
        return 1;
    }
    if ((ds_datyp = H5Dget_type(ds_id->dat_id)) < 0) return 1;
    if ((d_class = H5Tget_class(ds_datyp)) == -1) return 1;
    if ((native_typ = H5Tget_native_type(ds_datyp, H5T_DIR_ASCEND)) < 0)
        return 1;

    nval = 1;
    idim = 0;
    if (ndim > 0) {
        while (idim < ndim)
            nval *= dims[idim++];
    }
    /*
     *  based on the class, read the data in
     */
    switch (d_class) {
    case H5T_STRING:
        if (nval == 0) {
            if (H5Dread(ds_id->dat_id, native_typ, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, data) < 0) return 1;
        } else {
            if (H5Dread(ds_id->dat_id, native_typ, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, data) < 0) return 1;
        }
        break;
    case H5T_INTEGER:
    case H5T_FLOAT:
        if (H5Dread(ds_id->dat_id, native_typ, H5S_ALL, H5S_ALL,
                H5P_DEFAULT, data) < 0) return 1;
        break;
    case H5T_TIME:
    case H5T_BITFIELD:
    case H5T_OPAQUE:
    case H5T_COMPOUND:
    case H5T_REFERENCE:
    case H5T_ENUM:
    case H5T_VLEN:
    case H5T_ARRAY:
    default:
        printf("Unable to handle the class for dataset\n");
        return 1;
        break;
    }
    /*
     *  close the allocated spaces before leaving
     */
    if (H5Sclose(spaceid) < 0) return 1;
    if (H5Tclose(ds_datyp) < 0) return 1;
    if (H5Tclose(native_typ) < 0) return 1;
    return 0;
}

int h5io_rd_ds_slice(h5io_str *ds_id, int *start, int *count, void *data)
/*******************************************************************

   h5io_rd_ds

   purpose: read a slice from a dataset into an array

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        ds_id            I      id for the dataset in the file
      int *             start            I      Start (0 origin) of slice to
                                                read, covers all array dims
      int *             count            I      how many items to read in
                                                each dimension
      void *            data             O      allocated space for the data 
                                                to be read into

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 25 Sep 2008     Original development

   Note that when called, this routine will set up the HDF 5 I/O needed
   and will save settings for subsequent calls.  If the slice size
   (count) stays the same, minimal changes are made to permit this.
   closing the dataset ID will unset these settings and free the space
   needed for the slice I/O (unlike h5io_rd_ds, which does a 1-time 
   read of the data needed)

 *******************************************************************/ {
    hsize_t dims[H5IO_MAXDIM], maxdims[H5IO_MAXDIM], dim_arr_mem[H5IO_MAXDIM];
    hsize_t lstart[H5IO_MAXDIM], lcount[H5IO_MAXDIM];
    hsize_t block[H5IO_MAXDIM], stride[H5IO_MAXDIM];
    H5T_class_t d_class;
    hid_t ds_datyp;
    int changed_count, indx;
    /*
     *  make sure the incoming id is OK - it must be a dataset
     */
    if (ds_id->type != H5IO_TYP_DAT_ID) {
        printf("%s: A dataset id must be passed in for this routine\n", __FILE__);
        return 1;
    }
    /*
     *  If this is the initial slice read, we need to set up for the read
     */
    if (ds_id->dat_rw_mode == H5IO_DAT_RWALL) {
        /*
         *  initial read of a slice do all the checks and set up as in h5io_rd_ds
         *  but additionally set up the range of the memory space (file space
         *  range is variable for each slice so set it further on
         */
        if ((ds_id->fil_space_id = H5Dget_space(ds_id->dat_id)) < 0) return 1;
        if ((ds_id->prev_ndim =
                H5Sget_simple_extent_dims(ds_id->fil_space_id, dims, maxdims)) < 0)
            return 1;
        if (ds_id->prev_ndim > H5IO_MAXDIM) {
            printf("%s: # dimensions of dataset is > software maximum of %d\n",
                    __FILE__, H5IO_MAXDIM);
            return 1;
        }
        if ((ds_datyp = H5Dget_type(ds_id->dat_id)) < 0) return 1;
        if ((d_class = H5Tget_class(ds_datyp)) == -1) return 1;
        if ((ds_id->native_typ =
                H5Tget_native_type(ds_datyp, H5T_DIR_ASCEND)) < 0)
            return 1;
        /*
         *  the data type is no longer of use
         */
        if (H5Tclose(ds_datyp) < 0) return 1;
        /*
         *  check that the class is a readable one
         */
        switch (d_class) {
        case H5T_INTEGER:
        case H5T_FLOAT:
            break;
        case H5T_STRING:
        case H5T_TIME:
        case H5T_BITFIELD:
        case H5T_OPAQUE:
        case H5T_COMPOUND:
        case H5T_REFERENCE:
        case H5T_ENUM:
        case H5T_VLEN:
        case H5T_ARRAY:
        default:
            printf("Unable to handle the class for dataset\n");
            return 1;
            break;
        }
        /*
         *  set space size for memory (will usually not change)
         */
        for (indx = 0; indx < ds_id->prev_ndim; indx++) {
            dim_arr_mem[indx] = *(count + indx);
            *(ds_id->prev_count + indx) = *(count + indx);
            lstart[indx] = 0;
            block[indx] = 1;
            stride[indx] = 1;
            lcount[indx] = *(count + indx);
        }
        if ((ds_id->mem_space_id =
                H5Screate_simple(ds_id->prev_ndim, dim_arr_mem, dim_arr_mem)) < 0)
            return 1;
        /*  */
        if (H5Sselect_hyperslab(ds_id->mem_space_id, H5S_SELECT_SET,
                lstart, stride, lcount, block) < 0) return 1;
        /*
         *  create space descriptor for file (set sub-area later)
         */
        if ((ds_id->fil_space_id = H5Dget_space(ds_id->dat_id)) < 0) return 1;
        /*
         *  initialization over, say so
         */
        ds_id->dat_rw_mode = H5IO_DAT_RWSLICE;
    }
    /*
     *
     *  If the count (slice size) was changed from the previous call, deal
     *  with that here
     */
    changed_count = 0;
    for (indx = 0; indx < ds_id->prev_ndim; indx++)
        if (*(count + indx) != *(ds_id->prev_count + indx))
            changed_count = 1;
    /*  */
    if (changed_count == 1) {
        /*
         *  since the slice size changed, re-set the space for memory
         */
        if (H5Sclose(ds_id->mem_space_id) < 0) return 1;
        /* */
        for (indx = 0; indx < ds_id->prev_ndim; indx++) {
            dim_arr_mem[indx] = *(count + indx);
            *(ds_id->prev_count + indx) = *(count + indx);
            lstart[indx] = 0;
            block[indx] = 1;
            stride[indx] = 1;
            lcount[indx] = *(count + indx);
        }
        if ((ds_id->mem_space_id =
                H5Screate_simple(ds_id->prev_ndim, dim_arr_mem, dim_arr_mem)) < 0)
            return 1;
        /*  */
        if (H5Sselect_hyperslab(ds_id->mem_space_id, H5S_SELECT_SET,
                lstart, stride, lcount, block) < 0) return 1;
    }
    /*
     *
     *  for each slice read, select the file slice to read
     */
    for (indx = 0; indx < ds_id->prev_ndim; indx++) {
        lstart[indx] = *(start + indx);
        ;
        block[indx] = 1;
        stride[indx] = 1;
        lcount[indx] = *(count + indx);
    }
    if (H5Sselect_hyperslab(ds_id->fil_space_id, H5S_SELECT_SET,
            lstart, stride, lcount, block) < 0) return 1;
    /*
     *  and read it
     */
    if (H5Dread(ds_id->dat_id, ds_id->native_typ, ds_id->mem_space_id,
            ds_id->fil_space_id, H5P_DEFAULT, (void *) data) < 0) return 1;
    /*
     *  and finish
     */
    return 0;
}

int h5io_mk_grp(h5io_str *id, char *grp_nam, h5io_str *grp_id)
/*******************************************************************

   h5io_mk_grp

   purpose: Create a group in an HDF 5 file

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        id               I      id for the opened file,
                                                or group
      char *            grp_nam          I      new group name to create
      h5io_str *        grp_id           O      id of the new group

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 29 Sep 2008     Original development

 *******************************************************************/ {
    /*
     *  make sure the incoming id is OK
     */
    if (id->type == H5IO_TYP_DAT_ID) {
        printf("%s: Cannot set group under a dataset\n", __FILE__);
        return 1;
    }
    /*
     *  just find the id of the group
     */
    if ((grp_id->grp_id = H5Gcreate1(id->grp_id, grp_nam, 0)) < 0)
        return 1;
    /*
     *  and set type for a group id
     */
    grp_id->type = H5IO_TYP_GRP_ID;
    return 0;
}

int h5io_mk_ds(h5io_str *id, char *ds_name, hid_t type, int ndim,
        int *dim_siz, h5io_str *ds_id)
/*******************************************************************

   h5io_mk_ds

   purpose: Create a dataset in a HDF file

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        id               I      id for the opened file,
                                                or group
      char *            ds_name          I      name of dataset to create
      hid_t             type             I      type of the data to be created
                                                (see below for more)
      int               ndim             I      number of dimensions of 
                                                data array to store
      int *             dim_siz          I      array of size ndim containing
                                                each dimension's length in
                                                elements
      h5io_str *        ds_id            O      id of the dataset created

   Notes
     More on type - I hate to bring this HDF 5 value up, but this is an 
     HDF 5 descriptor of one of the standard types used.  Ususlly, 
     pre-defined native computer types can be used.  other types are also 
     available.  All the types can be seen in the HDF5 Datatypes web page 
     (version 1.6)
     http://www.hdfgroup.org/HDF5/doc1.6/UG/UG_frame.html.

     The common native types are:
     Code                 Corresponding C Type
     -------------------  --------------------
     H5T_NATIVE_CHAR      char
     H5T_NATIVE_SCHAR     signed char
     H5T_NATIVE_UCHAR     unsigned char
     H5T_NATIVE_SHORT     short
     H5T_NATIVE_USHORT    unsigned short
     H5T_NATIVE_INT       int
     H5T_NATIVE_UINT      unsigned
     H5T_NATIVE_LONG      long
     H5T_NATIVE_ULONG     unsigned long
     H5T_NATIVE_LLONG     long long
     H5T_NATIVE_ULLONG    unsigned long long
     H5T_NATIVE_FLOAT     float
     H5T_NATIVE_DOUBLE    double
     H5T_NATIVE_LDOUBLE   long double

     Placing of strings is probably done in an entirely different way
     not handled here (this should be rare normally for datasets)

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 29 Sep 2008     Original development

 *******************************************************************/ {
    int idim;
    hsize_t ldim_siz[H5IO_MAXDIM];
    hid_t fil_space_id;
    /*
     *  make sure the incoming id is OK
     */
    if (id->type == H5IO_TYP_DAT_ID) {
        printf("%s: Cannot create a dataset under a dataset\n", __FILE__);
        return 1;
    }
    /*
     *  check inputs and transfer
     */
    if (ndim > H5IO_MAXDIM) {
        printf("%s: The created dataset cannot currently exceed %d dimensions\n",
                __FILE__, H5IO_MAXDIM);
        return 1;
    }
    for (idim = 0; idim < ndim; idim++)
        ldim_siz[idim] = dim_siz[idim];
    /*
     *  set up the file's data space
     */
    if ((fil_space_id = H5Screate_simple(ndim, ldim_siz, ldim_siz)) < 0)
        return 1;
    /*
     *  just find the id of the dataset
     */
    if ((ds_id->dat_id = H5Dcreate1(id->grp_id, ds_name, type,
            fil_space_id, H5P_DEFAULT)) < 0) return 1;
    /*
     *  and set type for a data id, default to entire data array read,
     *  and free the space id for now
     */
    ds_id->type = H5IO_TYP_DAT_ID;
    ds_id->dat_rw_mode = H5IO_DAT_RWALL;
    if (H5Sclose(fil_space_id) < 0) return 1;
    return 0;
}

int h5io_wr_attr(h5io_str *id, char *attr_name, hid_t out_type,
        int ndim, int *dim_siz, void *data)
/*******************************************************************

   h5io_wr_attr

   purpose: write data as an attribute to a hdf 5 file, group, or dataset
     Note that for strings, use h5io_wr_attr_str

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        id               I      id for the opened file,
                                                group, or dataset
      char *            attr_name        I      name of attribute to get
      hid_t             out_type         I      type of the data to be created
                                                (see h5io_mk_ds for more on
                                                type)
                                                The type being written is 
                                                assumed to be the same 
                                                general type (ie float, short)
      int               ndim             I      number of dimensions of
                                                data array to store if 0,
                                                store 1 value as a scalar
      int *             dim_siz          I      array of size ndim containing
                                                each dimension's length in
                                                elements (ignored if ndim = 0)
      void *            data             I      pointer to space made for 
                                                the data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 29 Sep 2008     Original development

 *******************************************************************/ {
    hid_t *base_grp, fil_spac_id, attr_id, native_type;
    hsize_t ldim_siz[H5IO_MAXDIM];
    int idim;
    /*
     *  make sure the # dims is OK
     */
    if (ndim > H5IO_MAXDIM) {
        printf("%s: The created dataset cannot currently exceed %d dimensions\n",
                __FILE__, H5IO_MAXDIM);
        return 1;
    }
    /*
     *  depending on array vs scalar, set the file space
     */
    if (ndim == 0) {
        if ((fil_spac_id = H5Screate(H5S_SCALAR)) < 0) return 1;
    } else {
        for (idim = 0; idim < ndim; idim++)
            ldim_siz[idim] = dim_siz[idim];
        if ((fil_spac_id = H5Screate_simple(ndim, ldim_siz, ldim_siz)) < 0)
            return 1;
    }
    /*
     *  create the attribute
     */
    base_grp = (id->type == H5IO_TYP_DAT_ID) ?
            &(id->dat_id) : &(id->grp_id);
    if ((attr_id = H5Acreate1(*base_grp, attr_name, out_type, fil_spac_id,
            H5P_DEFAULT)) < 0) return 1;

    if ((native_type = H5Tget_native_type(out_type, H5T_DIR_ASCEND)) < 0)
        return 1;
    /*
     *  write data out
     */
    if (H5Awrite(attr_id, native_type, data) < 0) return 1;
    /*
     *  close the allocated spaces before leaving
     */
    if (H5Aclose(attr_id) < 0) return 1;
    if (H5Sclose(fil_spac_id) < 0) return 1;
    if (H5Tclose(native_type) < 0) return 1;
    return 0;
}

int h5io_wr_attr_str(h5io_str *id, char *attr_name, int ndim, int *dim_siz, int str_len, char *data)
/*******************************************************************

   h5io_wr_attr_str

   purpose: write a string (or array of strings) to a hdf 5 attribute in
     a file, group, or dataset
     Note that for non-strings, use h5io_wr_attr

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        id               I      id for the opened file,
                                                group, or dataset
      char *            attr_name        I      name to give the attribute 
      int               ndim             I      number of dimensions of
                                                data array to store if 0,
                                                store 1 value as a scalar
      int *             dim_siz          I      array of size ndim containing
                                                each dimension's length in
                                                elements (ignored if ndim = 0)
      int               str_len          I      length of string(s) - this
                                                is required for an array, 
                                                but for ndim = 0 (scalar),
                                                if this is <= 0, the size
                                                will be computed with strlen()
      char *            data             O      pointer to space made for 
                                                the data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 29 Sep 2008     Original development

 *******************************************************************/ {
    hid_t *base_grp, fil_spac_id, attr_id, attr_typ;
    hsize_t ldim_siz[H5IO_MAXDIM];
    size_t lstr_len;
    int idim;
    /*
     *  make sure the # dims is OK
     */
    if (ndim > H5IO_MAXDIM) {
        printf("%s: The created dataset cannot currently exceed %d dimensions\n",
                __FILE__, H5IO_MAXDIM);
        return 1;
    }
    /*
     *  depending on array vs scalar, set the file space
     */
    if (ndim == 0) {
        if ((fil_spac_id = H5Screate(H5S_SCALAR)) < 0) return 1;
        lstr_len = (str_len <= 0) ? strlen(data) : str_len;
    } else {
        for (idim = 0; idim < ndim; idim++)
            ldim_siz[idim] = dim_siz[idim];
        lstr_len = str_len;
        if ((fil_spac_id = H5Screate_simple(ndim, ldim_siz, ldim_siz)) < 0)
            return 1;
    }
    /*
     *  create the attribute
     */
    base_grp = (id->type == H5IO_TYP_DAT_ID) ?
            &(id->dat_id) : &(id->grp_id);

    if ((attr_typ = H5Tcopy(H5T_C_S1)) < 0) return 1;
    if (H5Tset_size(attr_typ, lstr_len) < 0) return 1;
    if (H5Tset_strpad(attr_typ, H5T_STR_NULLTERM) < 0) return 1;

    if ((attr_id = H5Acreate1(*base_grp, attr_name, attr_typ, fil_spac_id,
            H5P_DEFAULT)) < 0) return 1;
    /*
     *  write data out
     */
    if (H5Awrite(attr_id, attr_typ, (void *) data) < 0) return 1;
    /*
     *  close the allocated spaces before leaving
     */
    if (H5Aclose(attr_id) < 0) return 1;
    if (H5Sclose(fil_spac_id) < 0) return 1;
    if (H5Tclose(attr_typ) < 0) return 1;
    return 0;
}

int h5io_wr_ds(h5io_str *ds_id, void *data)
/*******************************************************************

   h5io_wr_ds

   purpose: write an entire array into a HDF 5 dataset

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        ds_id            I      id for the dataset in the file
      void *            data             O      allocated space for the data 
                                                to be read into

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 30 Sep 2008     Original development

 *******************************************************************/ {
    hid_t spaceid, ds_datyp, native_typ;
    hsize_t dims[H5IO_MAXDIM], maxdims[H5IO_MAXDIM];
    H5T_class_t d_class;
    int ndim;
    /*
     *  make sure the incoming id is OK - it must be a dataset
     */
    if (ds_id->type != H5IO_TYP_DAT_ID) {
        printf("%s: A dataset id must be passed in for this routine\n", __FILE__);
        return 1;
    }
    /*
     *  just get the information on dimensions from the previously created space
     */
    if ((spaceid = H5Dget_space(ds_id->dat_id)) < 0) return 1;
    if ((ndim = H5Sget_simple_extent_dims(spaceid, dims, maxdims)) < 0)
        return 1;
    if (ndim > H5IO_MAXDIM) {
        printf("%s: # dimensions of dataset is > software maximum of %d\n",
                __FILE__, H5IO_MAXDIM);
        return 1;
    }
    if ((ds_datyp = H5Dget_type(ds_id->dat_id)) < 0) return 1;
    if ((d_class = H5Tget_class(ds_datyp)) == -1) return 1;
    if ((native_typ = H5Tget_native_type(ds_datyp, H5T_DIR_ASCEND)) < 0)
        return 1;
    /*
     *  check the class and write if int or float
     */
    switch (d_class) {
    case H5T_INTEGER:
    case H5T_FLOAT:
        if (H5Dwrite(ds_id->dat_id, native_typ, H5S_ALL, H5S_ALL,
                H5P_DEFAULT, data) < 0) return 1;
        break;
    case H5T_STRING:
    case H5T_TIME:
    case H5T_BITFIELD:
    case H5T_OPAQUE:
    case H5T_COMPOUND:
    case H5T_REFERENCE:
    case H5T_ENUM:
    case H5T_VLEN:
    case H5T_ARRAY:
    default:
        printf("Unable to handle the current class for dataset\n");
        return 1;
        break;
    }
    /*
     *  close the allocated spaces before leaving
     */
    if (H5Sclose(spaceid) < 0) return 1;
    if (H5Tclose(ds_datyp) < 0) return 1;
    if (H5Tclose(native_typ) < 0) return 1;
    return 0;
}

int h5io_wr_ds_slice(h5io_str *ds_id, int *start, int *count, void *data)
/*******************************************************************

   h5io_rd_ds

   purpose: write a slice to a dataset from an array
     similar to the read slice routine, h5io_rd_ds_slice

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        ds_id            I      id for the dataset in the file
      void *            data             O      data array to write

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 30 Sep 2008     Original development

   Note that when called, this routine will set up the HDF 5 I/O needed
   and will save settings for subsequent calls.  If the slice size
   (count) stays the same, minimal changes are made to permit this.
   closeing the dataset ID will unset these settings and free the space
   needed for the slice I/O (unlike h5io_rd_ds, which does a 1-time 
   read of the data needed)

 *******************************************************************/ {
    hsize_t dims[H5IO_MAXDIM], maxdims[H5IO_MAXDIM], dim_arr_mem[H5IO_MAXDIM];
    hsize_t lstart[H5IO_MAXDIM], lcount[H5IO_MAXDIM];
    hsize_t block[H5IO_MAXDIM], stride[H5IO_MAXDIM];
    H5T_class_t d_class;
    hid_t ds_datyp;
    int changed_count, indx;
    /*
     *  make sure the incoming id is OK - it must be a dataset
     */
    if (ds_id->type != H5IO_TYP_DAT_ID) {
        printf("%s: A dataset id must be passed in for this routine\n", __FILE__);
        return 1;
    }
    /*
     *
     *  If this is the initial slice write, we need to do the set up for it
     */
    if (ds_id->dat_rw_mode == H5IO_DAT_RWALL) {
        /*
         *  initial write of a slice - do all the checks and set up as in h5io_wr_ds
         *  but additionally set up the range of the memory space (file space
         *  range is variable for each slice so set it further on)
         */
        if ((ds_id->fil_space_id = H5Dget_space(ds_id->dat_id)) < 0) return 1;
        if ((ds_id->prev_ndim =
                H5Sget_simple_extent_dims(ds_id->fil_space_id, dims, maxdims)) < 0)
            return 1;
        if (ds_id->prev_ndim > H5IO_MAXDIM) {
            printf("%s: # dimensions of dataset is > software maximum of %d\n",
                    __FILE__, H5IO_MAXDIM);
            return 1;
        }
        if ((ds_datyp = H5Dget_type(ds_id->dat_id)) < 0) return 1;
        if ((d_class = H5Tget_class(ds_datyp)) == -1) return 1;
        if ((ds_id->native_typ =
                H5Tget_native_type(ds_datyp, H5T_DIR_ASCEND)) < 0)
            return 1;
        /*
         *  the data type is no longer of use
         */
        if (H5Tclose(ds_datyp) < 0) return 1;
        /*
         *  check that the class is a readable one
         */
        switch (d_class) {
        case H5T_INTEGER:
        case H5T_FLOAT:
            break;
        case H5T_STRING:
        case H5T_TIME:
        case H5T_BITFIELD:
        case H5T_OPAQUE:
        case H5T_COMPOUND:
        case H5T_REFERENCE:
        case H5T_ENUM:
        case H5T_VLEN:
        case H5T_ARRAY:
        default:
            printf("Unable to handle the class for dataset\n");
            return 1;
            break;
        }
        /*
         *  set space size for memory (will usually not change)
         */
        for (indx = 0; indx < ds_id->prev_ndim; indx++) {
            dim_arr_mem[indx] = *(count + indx);
            *(ds_id->prev_count + indx) = *(count + indx);
            lstart[indx] = 0;
            block[indx] = 1;
            stride[indx] = 1;
            lcount[indx] = *(count + indx);
        }
        if ((ds_id->mem_space_id =
                H5Screate_simple(ds_id->prev_ndim, dim_arr_mem, dim_arr_mem)) < 0)
            return 1;
        /*  */
        if (H5Sselect_hyperslab(ds_id->mem_space_id, H5S_SELECT_SET,
                lstart, stride, lcount, block) < 0) return 1;
        /*
         *  create space descriptor for file (set sub-area later)
         */
        if ((ds_id->fil_space_id = H5Dget_space(ds_id->dat_id)) < 0) return 1;
        /*
         *  initialization over, say so
         */
        ds_id->dat_rw_mode = H5IO_DAT_RWSLICE;
    }
    /*
     *
     *  If the count (slice size) was changed from the previous call, deal
     *  with that here
     */
    changed_count = 0;
    for (indx = 0; indx < ds_id->prev_ndim; indx++)
        if (*(count + indx) != *(ds_id->prev_count + indx))
            changed_count = 1;
    /*  */
    if (changed_count == 1) {
        /*
         *  since the slice size changed, re-set the space for memory
         */
        if (H5Sclose(ds_id->mem_space_id) < 0) return 1;
        /* */
        for (indx = 0; indx < ds_id->prev_ndim; indx++) {
            dim_arr_mem[indx] = *(count + indx);
            *(ds_id->prev_count + indx) = *(count + indx);
            lstart[indx] = 0;
            block[indx] = 1;
            stride[indx] = 1;
            lcount[indx] = *(count + indx);
        }
        if ((ds_id->mem_space_id =
                H5Screate_simple(ds_id->prev_ndim, dim_arr_mem, dim_arr_mem)) < 0)
            return 1;
        /*  */
        if (H5Sselect_hyperslab(ds_id->mem_space_id, H5S_SELECT_SET,
                lstart, stride, lcount, block) < 0) return 1;
    }
    /*
     *  select the file space for each slice written
     */
    for (indx = 0; indx < ds_id->prev_ndim; indx++) {
        lstart[indx] = *(start + indx);
        ;
        block[indx] = 1;
        stride[indx] = 1;
        lcount[indx] = *(count + indx);
    }
    if (H5Sselect_hyperslab(ds_id->fil_space_id, H5S_SELECT_SET,
            lstart, stride, lcount, block) < 0) return 1;
    /*
     *  and write it
     */
    if (H5Dwrite(ds_id->dat_id, ds_id->native_typ, ds_id->mem_space_id,
            ds_id->fil_space_id, H5P_DEFAULT, (void *) data) < 0) return 1;
    /*
     *  and finish
     */
    return 0;
}

int h5io_grab_ds(h5io_str *id, char *path_name, void *data)
/*******************************************************************

   h5io_grab_ds

   purpose: get the data from a dataset without keeping it open

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        id               I      id for the opened file,
                                                or group
      char *            path_name        I      path to set from current 
                                                location with the dataset name
                                                as the last part of the path
      void *            data             O      allocated space for the data
                                                to be read into

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 02 Oct 2009     Original development

 *******************************************************************/ {
    h5io_str ds_id;
    /*
     *  This is a convenience, built from h5io_set_ds, h5io_rd_ds, and h5io_close
     */
    if (h5io_set_ds(id, path_name, &ds_id) != 0) {
        printf("%s:%d - failed to open the %s dataset\n", __FILE__,
                __LINE__, path_name);
        return 1;
    }
    if (h5io_rd_ds(&ds_id, data) != 0) {
        printf("%s:%d - failed to read the %s dataset\n", __FILE__,
                __LINE__, path_name);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s:%d - failed to close the %s dataset\n", __FILE__,
                __LINE__, path_name);
        return 1;
    }
    return 0;
}

int h5io_grab_ds_slice(h5io_str *id, char *path_name, int *start,
        int *count, void *data)
/*******************************************************************

   h5io_grab_ds_slice

   purpose: get a slice of data from a dataset without keeping it open

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        id               I      id for the opened file,
                                                or group
      char *            path_name        I      path to set from current
                                                location with the dataset name
                                                as the last part of the path
      int *             start            I      Start (0 origin) of slice to
                                                read, covers all array dims
      int *             count            I      how many items to read in
                                                each dimension
      void *            data             O      allocated space for the data
                                                to be read into

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 09 Nov 2011     Original development

 *******************************************************************/ {
    h5io_str ds_id;
    /*
     *  This is a convenience, built from h5io_set_ds, h5io_rd_ds_slice, 
     *  and h5io_close
     */
    if (h5io_set_ds(id, path_name, &ds_id) != 0) {
        printf("%s:%d - failed to open the %s dataset\n", __FILE__,
                __LINE__, path_name);
        return 1;
    }
    if (h5io_rd_ds_slice(&ds_id, start, count, data) != 0) {
        printf("%s:%d - failed to read the %s dataset\n", __FILE__,
                __LINE__, path_name);
        return 1;
    }
    if (h5io_close(&ds_id) != 0) {
        printf("%s:%d - failed to close the %s dataset\n", __FILE__,
                __LINE__, path_name);
        return 1;
    }
    return 0;
}

int h5io_attr_exist(h5io_str *id, char *attr_name)
/*******************************************************************

   h5io_rd_attr

   purpose: read a hdf 5 attribute from the file, group, or dataset
     does an attribute exist?

   Returns type: int - return status: 0 is good
                        1 if no attribute by that name exists
                       -1 if an error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        id               I      id for the opened file,
                                                group, or dataset
      char *            attr_name        I      name of attribute to get

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 23 Sep 2008     Original development
      W. Robinson, SAIC 4 Oct 2016      add a H5Aclose to stem leaks

 *******************************************************************/ {
    int attr_exist;
    hid_t *base_grp, attr_id;
    H5E_auto_t old_func;
    void *old_client_data;
    /*
     *  set the error reporting off
     */
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    /*
     *  open the attribute
     */
    base_grp = (id->type == H5IO_TYP_DAT_ID) ?
            &(id->dat_id) : &(id->grp_id);
    attr_exist = 0;
    if ((attr_id = H5Aopen_name(*base_grp, attr_name)) < 0) {
        attr_exist = 1;
        H5Aclose(attr_id);
    }
    /*
     *  return to chatty error reporting
     */
    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);
    /*
     *  and return result
     */
    return attr_exist;
}

int h5io_grp_contents(h5io_str *id, int *n_obj, char ***o_names,
        int **types)
/*******************************************************************

   h5io_grp_contents

   purpose: get a list of object names and types stored in a group

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        id               I      id of group to use
      int *             n_obj            O      # of objects found
      char ***          o_names          O      pointer to string array of 
                                                names maintained in this 
                                                function
      int **            types            O      pointer to a list of object 
                                                types found
                                                values will be:
                                                H5G_GROUP - a group
                                                H5G_DATASET - a dataset

   NOTE - the o_names and types are local storage that will be re-used
    upon another call.  if the results are needed longer, they must 
    be transfered to other storage.

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 13 Oct 2010     Original development

 *******************************************************************/ {
    static char **lo_names;
    static int *ltypes, ln_obj = 0;
    H5G_info_t g_info;
    H5G_stat_t info;
    int i, len, len0;
    /*
     *  make sure the id is not for a dataset
     */
    if (id->type == H5IO_TYP_DAT_ID) {
        printf("%s, %d: Error: A dataset cannot be the input id\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  if previously used, free the storage (less memory leaks)
     */
    if (ln_obj != 0) {
        free(ltypes);
        for (i = 0; i < ln_obj; i++)
            free(lo_names[i]);
        free(lo_names);
        ln_obj = 0;
    }
    /*
     *  find the # of objects (attributes left out for some strange reason)
     */
    if (H5Gget_info(id->grp_id, &g_info) < 0) {
        return 1;
    }
    *n_obj = g_info.nlinks;
    ln_obj = *n_obj;
    lo_names = (char **) malloc(*n_obj * sizeof ( char *));
    ltypes = (int *) malloc(*n_obj * sizeof ( int));
    /*
     *  get the object names and types
     */
    for (i = 0; i < *n_obj; i++) {
        if ((len0 = H5Gget_objname_by_idx(id->grp_id, i, NULL, 500))
                < 0) {
            printf("%s, %d - H5Gget_objname_by_idx 1 failed\n",
                    __FILE__, __LINE__);
            return 1;
        }
        len0++;
        lo_names[i] = (char *) malloc(len0 * sizeof ( char));
        if ((len = H5Gget_objname_by_idx(id->grp_id, i, lo_names[i], len0))
                < 0) {
            printf("%s, %d - H5Gget_objname_by_idx failed\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if (H5Gget_objinfo(id->grp_id, lo_names[i], 0, &info) < 0) {
            printf("%s, %d - H5Gget_objinfo failed\n",
                    __FILE__, __LINE__);
            return 1;
        }
        *(ltypes + i) = info.type;
    }
    *o_names = lo_names;
    *types = ltypes;
    return 0;
}

int h5io_inq_path(h5io_str *id, char *path, int *found_path)
/*******************************************************************

   h5io_inq_path

   purpose: determine if a path defined by group names, exists 
            below the current location

   Returns type: int - return status: 0 is good
                        1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      h5io_str *        id               I      id of group to use
      char *            path             I      path of group names, like
                                                'Data/Products/May'
                                                Note that last part of path can
                                                also be a dataset name
      int *             found_path       O      indicates path exists if 1
                                                non-existant if 0

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 24 May 2013     Original development

 *******************************************************************/ {
    char *cptr, *tok, str_wrk[FILENAME_MAX], **grp_names;
    int32_t nelts, ielt, elt_num, ngid, nobj, iobj, igrp, found_elt;
    /* nelts # path elements, ngid # group IDs open */
    int *types;
    h5io_str *ar_gid, sgid;
    /*
     *  make sure the id is not for a dataset
     */
    if (id->type == H5IO_TYP_DAT_ID) {
        printf("%s, %d, E: h5io_inq_path, a dataset cannot be the input id\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  count the # elements in the path
     */
    nelts = 0;
    strcpy(str_wrk, path);
    cptr = str_wrk;
    while ((tok = strtok(cptr, "/")) != NULL) {
        nelts++;
        cptr = NULL;
    }

    ngid = nelts - 1;
    /*
     *  set up the list of group ids for the elements of the path
     *  and go through checking out the path sent in
     */
    ar_gid = (h5io_str *) malloc(ngid * sizeof ( h5io_str));
    strcpy(str_wrk, path);
    sgid = *id;
    cptr = str_wrk;
    for (ielt = 0; ielt < nelts; ielt++) {
        /* get the next path element = tok */
        if ((tok = strtok(cptr, "/")) == NULL) {
            printf("%s, %d, E: h5io_inq_path, tokens exhausted unexpectedly\n",
                    __FILE__, __LINE__);
            return 1;
        }
        /* get contents of current group */
        if (h5io_grp_contents(&sgid, &nobj, &grp_names, &types) != 0) {
            printf("%s, %d: h5io_inq_path problem, exiting\n", __FILE__, __LINE__);
            return 1;
        }
        /*  look through contents names and match up with path element */
        found_elt = 0;
        if (nobj > 0) {
            for (iobj = 0; iobj < nobj; iobj++) {
                if (strcmp(grp_names[iobj], tok) == 0) {
                    found_elt = 1;
                    elt_num = iobj;
                    break;
                }
            }
        }
        if (found_elt == 0) {
            *found_path = 0;
            if (ielt >= 1)
                for (igrp = 0; igrp < (ielt - 1); igrp++)
                    h5io_close(&(ar_gid[igrp]));
            free(ar_gid);
            return 0;
        }
        /*  make sure that all but the last path element is a group */
        if ((ielt < (nelts - 1)) && (types[elt_num] != H5G_GROUP)) {
            printf("%s, %d, E: h5io_inq_path\n", __FILE__, __LINE__);
            printf("    non-group found before the last element in the path\n");
            printf("    path: %s\n", path);
            printf("    non-group element: %s\n", grp_names[elt_num]);
            return 1;
        }
        /*
         *  set to the next element in the path and check
         */
        if (ielt < (nelts - 1)) {
            if (h5io_set_grp(&sgid, tok, &ar_gid[ielt]) != 0) {
                printf("%s, %d: h5io_inq_path problem, exiting\n",
                        __FILE__, __LINE__);
                return 1;
            }
        }
        sgid = ar_gid[ielt];
        cptr = NULL;
    }
    /*
     *  if we've arrived here, the path is there.  close all the ids
     *  and return success
     */
    for (igrp = 0; igrp < ngid; igrp++)
        h5io_close(&(ar_gid[igrp]));

    *found_path = 1;
    free(ar_gid);

    return 0;
}
