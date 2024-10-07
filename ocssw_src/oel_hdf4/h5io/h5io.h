/*******************************************************************

   h5io.h

   purpose: include file for the use of the hdf I/O routine aids

   Parameters: 
      Type              Name            I/O     Description
      ----              ----            ---     -----------

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 23-Sep-2008     Original development
      S. Bailey, Futuretech 01-Nov-2010	removed malloc.h header file
 *******************************************************************/

#ifndef __H5IO_H__
#define __H5IO_H__
/*
 *  Note that hdf5.h is needed for this
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>

#define H5IO_MAXDIM 6  /* the max # dimensions handled currently */

/* Identify the type of ID the h5io_str is */
#define H5IO_TYP_FILE_ID 0  /* this is a file ID */
#define H5IO_TYP_GRP_ID 1   /* this is a group ID */
#define H5IO_TYP_DAT_ID 2   /* this is a dataset ID */
#define H5IO_DAT_RWALL 0   /* When read/write dataset, do all */
#define H5IO_DAT_RWSLICE 1  /* When read/write dataset, do slices */

/*
 *  There are several HDF 5 IDs  handles needed for steps in the I/O and 
 *  different ones in this structure are set depending on the operation.
 *  The type will determine which is set and can be used to make sure the
 *  proper id is in use for the proper operation.
 */
struct h5io_str_d {
    int type;
    hid_t file_id; /* set with type H5IO_TYP_FILE_ID */
    hid_t grp_id; /* set with type H5IO_TYP_FILE_ID or H5IO_TYP_GRP_ID */
    hid_t dat_id; /* set with type H5IO_TYP_DAT_ID */
    int dat_rw_mode; /* data read / write mode: either all or in slices, */
    /* with several calls to do entire I/O  */
    /* defaults to all = H5IO_DAT_RWALL on open of dataset */
    hid_t fil_space_id; /* descriptors for slab I/O - these are only set*/
    hid_t mem_space_id; /* (and freed) if dat_rw_mode = H5IO_DAT_RWSLICE */
    hid_t native_typ; /* needed for repeated I/O */
    int prev_ndim; /* # dims of the dataset just needed for ref */
    hsize_t prev_count[ H5IO_MAXDIM ]; /* make sure slab read set right for 
                                          next read  */
};

typedef struct h5io_str_d h5io_str;

/*
 *  prototypes
 */
int h5io_openw(char *, int, h5io_str *);
int h5io_openr(char *, int, h5io_str *);
int h5io_rd_attr(h5io_str *, char *, void *);
int h5io_set_grp(h5io_str *, char *, h5io_str *);
int h5io_set_ds(h5io_str *, char *, h5io_str *);
int h5io_rd_ds(h5io_str *, void *);
int h5io_rd_ds_slice(h5io_str *, int *, int *, void *);
int h5io_close(h5io_str *);
int h5io_mk_grp(h5io_str *, char *, h5io_str *);
int h5io_mk_ds(h5io_str *, char *, hid_t, int, int *, h5io_str *);
int h5io_wr_attr(h5io_str *, char *, hid_t, int, int *, void *);
int h5io_wr_attr_str(h5io_str *, char *, int, int *, int, char *);
int h5io_wr_ds(h5io_str *, void *);
int h5io_wr_ds_slice(h5io_str *, int *, int *, void *);
int h5io_attr_exist(h5io_str *, char *);
/*  currently, h5io_info will do this job but keep as backup
int h5io_attr_info( h5io_str *, char *, H5T_class_t *, hid_t *, int*, 
  int *, int * );
 */
int h5io_info(h5io_str *, char *, H5T_class_t *, hid_t *, int*,
        int *, int *);
int h5io_grab_ds(h5io_str *, char *, void *);
int h5io_grab_ds_slice(h5io_str *, char *, int *, int *, void *);
int h5io_grp_contents(h5io_str *, int *, char ***, int **);
int h5io_inq_path(h5io_str *, char *, int *);

#endif // __H5IO_H__
