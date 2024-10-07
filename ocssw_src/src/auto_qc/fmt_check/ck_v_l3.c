#include "fmt_check.h"
#include <math.h>
#include <hdf4utils.h>

extern int fmt_status; /*  format checking status, see fmt_check */
extern int verbose; /* 0 - don't print info on each vgroup, 1 do  */
static l3_org_str l3_org; /* info on l3 organization */

void ck_v_l3(int32 file_id, fmt_str *fmt)
/*******************************************************************

   ck_v_l3

   purpose: Check the format of the Vdatas in Vgroup 'Level-3 Binned Data'
   for a level 3 dataset

   Returns type: none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             file_id          I      ID returned from Hopen
      fmt_str *         fmt              I      Format expected information

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       4-Apr-2005      Original development

 *******************************************************************/ {
    int32 nbins, tot_vgroup;
    int32 num_vgroups, dummy_i32, *ref_array, i, iret, ifound;

    /*
     * set up # bins locally
     */
    nbins = fmt->dim_id[ fmt->vg_nbin_indx ].dim_size;

    /*
     *  set up organization based on the resolve
     */
    l3_get_org(fmt->resolve, &l3_org);
    /* 
     *  open the thing for Vgroups 
     */
    printf("\n\nChecking the 'Level-3 Binned Data' Vgroup...\n\n");

    Vstart(file_id);

    /* get the ID for the top Vgroups
     *  ref man says to first use max_size = 0 to just get the #,
     *  so, we'll do this first too
     */
    if ((num_vgroups = Vlone(file_id, &dummy_i32, 0)) == -1) {
        printf("********************Program problem:\n");
        printf("Vgroup check: The Vlone #1 failed\n");
        fmt_status = fmt_status | 2;
    } else {
        if (verbose == 1)
            printf("the # of 'lone' Vgroups is %d\n", num_vgroups);
        /*  
         *  allocate space for the groups and get the IDs  
         */
        ref_array = (int32 *) malloc(num_vgroups * sizeof ( int32));
        if ((num_vgroups = Vlone(file_id, ref_array, num_vgroups)) == -1) {
            printf("********************Program problem:\n");
            printf("Vgroup check: The Vlone #2 failed\n");
            fmt_status = fmt_status | 2;
        } else {
            /*  
             *  next, loop thru the Vgroups and get class, name and # entries,
             *  = # of Vdatas in this case
             */
            ifound = 0;
            tot_vgroup = fmt->n_vgroup + 3;
            for (i = 0; i < num_vgroups; i++) {
                if ((iret = group_fnd(file_id, ref_array[i], nbins, tot_vgroup,
                        fmt->vg_info)) == 0) {
                    ifound = 1;
                }
            }
            if (ifound != 1) {
                printf("**** the Vgroup: 'Level-3 Binned Data' was not found\n");
                fmt_status = fmt_status | 2;
            }
        }
    }
    /* close the file  */
    Vend(file_id);
}

int32 group_fnd(int32 file_id, int32 vgroup_ref, int32 nbins, int32 n_vgroup,
        vg_info_str *vd_check)
/*******************************************************************

   group_fnd

   purpose: Check the format of the Vdatas in Vgroup 'Level-3 Binned Data'
   for a level 3 dataset, archive version 2.6

   Returns type: int32

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             file_id         I       file id from Hopen
      int32             vgroup_ref      I       current group ref #
      int32             nbins           I       # bins from attributes
                                                to check
      int32             n_vgroup        I       # of vgroups (3 + # prod)
      struct vd_check_str * vd_check    I       array of definitions of each 
                                                Vdata

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       16-Mar-1995     Original development
      W. Robinson, SAIC 6 Apr 2005      uogrades for table driven version

 *******************************************************************/ {
    int32 v_id, n_vgroup_entries;
    int32 vd_id, i, n_vd_fields, vdata_size, interlace;
    int32 n_rec, j, vd_typ, vd_ref;
    int32 *begin, *extent;
    char vgroup_class[VGNAMELENMAX], vgroup_name[VGNAMELENMAX];
    char vd_name[VSNAMELENMAX], vd_class[VSNAMELENMAX];
    char field_list[500];

    /*
     *  set up storage for start bin and # in each row
     */
    begin = (int32 *) malloc(l3_org.numrows * sizeof (int32));
    extent = (int32 *) malloc(l3_org.numrows * sizeof (int32));
    /* 
     *  attach to the Vgroup first
     */
    if ((v_id = Vattach(file_id, vgroup_ref, "r")) == -1) {
        printf("******************** Probable Program Error\n");
        printf("Problem on attaching to group with ref # %d\n", vgroup_ref);
        fmt_status = fmt_status | 2;
        return -1;
    }
    /*
     *  get the name and check with the desired name
     */
    if (Vinquire(v_id, &n_vgroup_entries, vgroup_name) == -1) {
        printf("******************** Probable Program Error\n");
        printf("Vinquire failed on group ref # %d\n",
                vgroup_ref);
        fmt_status = fmt_status | 2;
        return -1;
    }

    if (verbose == 1)
        printf("The group name being checked is:\n'%s'\n", vgroup_name);
    if (strcmp(vgroup_name, "Level-3 Binned Data") != 0) {
        return 1;
    }
    /*
     *  OK, if the name has passed, assume we have the right Vgroup
     *  and check the class and # entries
     */
    Vgetclass(v_id, vgroup_class);
    if (verbose == 1) {
        printf("For target Vgroup of '%s'\n# entries = %d and class = '%s'\n",
                vgroup_name, n_vgroup_entries, vgroup_class);
    }
    if (n_vgroup_entries != n_vgroup) {
        printf("**** Error, # entries (%d) is not = %d\n", n_vgroup_entries,
                n_vgroup);
        fmt_status = fmt_status | 2;
    }
    if (strcmp(vgroup_class, "PlanetaryGrid") != 0) {
        printf("**** Error, Vgroup class ('%s') is incorrect\n", vgroup_class);
        fmt_status = fmt_status | 2;
    }
    /*
     *  Go on to the Vdatas themselvs
     */
    printf("\n\nChecking Vdatas inside the Vgroup\n\n");

    for (i = 0; i < n_vgroup_entries; i++) {
        /*
         *  try to find each Vdata
         */
        if ((vd_ref = VSfind(file_id, vd_check[i].name)) < 0) {
            printf("**** Unable to locate Vdata with name '%s'\n",
                    vd_check[i].name);
            fmt_status = fmt_status | 2;
        } else {
            /*
             *  attach to this Vdata
             */
            if ((vd_id = VSattach(file_id, vd_ref, "r")) == -1) {
                printf(
                        "*************** Program Error: Unable to attach to Vdata name: '%s'\n",
                        vd_check[i].name);
                printf("          and vd_ref # %d\n", vd_ref);
                fmt_status = fmt_status | 2;
            } else {
                /*
                 *  get and check class
                 */
                VSgetclass(vd_id, vd_class);
                if (verbose == 1)
                    printf("Vdata: '%s' has class '%s'\n", vd_check[i].name, vd_class);
                if (strcmp(vd_class, vd_check[i].class) != 0) {
                    printf("****for Vdata: '%s', mismatch in the class\n",
                            vd_check[i].name);
                    printf("expected class: '%s'\n", vd_check[i].class);
                    printf("read class:     '%s'\n", vd_class);
                    fmt_status = fmt_status | 2;
                }
                /*
                 *  get and check the # records and field names
                 */
                if (VSinquire(vd_id, &n_rec, &interlace, field_list,
                        &vdata_size, vd_name) == -1) {
                    printf(
                            "******************* Program Error, VSinquire failed\n");
                    fmt_status = fmt_status | 2;
                } else {
                    if (verbose == 1)
                        printf("Vdata: '%s' has #records = %d, field list:\n'%s'\n",
                            vd_check[i].name, n_rec, field_list);
                    if (strcmp(field_list, vd_check[i].field_list) != 0) {
                        printf("****for Vdata: '%s', mismatch in the field list\n",
                                vd_check[i].name);
                        printf("expected list: '%s'\n", vd_check[i].field_list);
                        printf("read     list: '%s'\n", field_list);
                        fmt_status = fmt_status | 2;
                    }
                    if (vd_check[i].nrec == -1) {
                        /*
                         *  a note here, the # records in one of these Vdatas
                         *  can actually be >= the attribute 'Data Bins'.
                         *  Thus, the following check
                         */
                        if (n_rec < nbins) {
                            printf("****for Vdata: '%s', # records < # bins\n",
                                    vd_check[i].name);
                            printf("expected: %d, read: %d\n", nbins, n_rec);
                            fmt_status = fmt_status | 2;
                        }
                    } else {
                        if (n_rec != vd_check[i].nrec) {
                            printf("****for Vdata: '%s', mismatch in the # records\n",
                                    vd_check[i].name);
                            printf("expected: %d, read: %d\n",
                                    vd_check[i].nrec, n_rec);
                            fmt_status = fmt_status | 2;
                        }
                    }
                }
                /*
                 *  get and check # fields and their types
                 */
                if ((n_vd_fields = VSgetfields(vd_id, field_list)) == -1) {
                    printf("*************** Program error, VSgetfields failed.\n");
                    fmt_status = fmt_status | 2;
                } else {
                    if (verbose == 1)
                        printf("Vdata: '%s', has %d fields\n",
                            vd_check[i].name, n_vd_fields);
                    if (n_vd_fields != vd_check[i].n_fields) {
                        printf("****for Vdata: '%s', mismatch in the # fields\n",
                                vd_check[i].name);
                        printf("expected: %d, read: %d\n", vd_check[i].n_fields,
                                n_vd_fields);
                        fmt_status = fmt_status | 2;
                    } else {
                        for (j = 0; j < n_vd_fields; j++) {
                            vd_typ = VFfieldtype(vd_id, j);
                            if (verbose == 1)
                                printf("Vdata: '%s', field # %d has num type %d\n",
                                    vd_check[i].name, j, vd_typ);
                            /* WDR - don't check type of field 5 of bin list
                             if( vd_typ != vd_check[i].num_typs[j] )
                             */
                            if ((j != 5) && (i != 2) &&
                                    (vd_typ != vd_check[i].num_typs[j])) {
                                printf(
                                        "****for Vdata: '%s', field # %d has incorrect number type\n",
                                        vd_check[i].name, j);
                                printf("expected: %d, read: %d\n",
                                        vd_check[i].num_typs[j], vd_typ);
                                fmt_status = fmt_status | 2;
                            }
                        }
                    }
                }
                /*
                 *  Vdata has been attached to.  any extra checking can be added
                 *  here
                 */
                if (i == 0) /*  Check constant values in the SEAGrid */ {
                    chk_sea_grid(vd_id);
                }                    /* Check all the BinIndex field values except for "begin" and 
	  * "extent" fields.  Save the begin and extent buffers until
	  * Vdata BinList is read then cross verify the information
	  * written in these fields against "binno" of BinList vdata.
	  */
                else if (i == 1) {
                    chk_bin_index(vd_id, begin, extent);
                } else if (i == 2) /* Check binno of BinList    */ {
                    chk_bin_list(vd_id, nbins, begin, extent);
                } else /* look through the products*/ {
                    chk_l3_prod(vd_id, nbins, vd_check[i].name, field_list);
                }
                /*
                 *  remember to detach from the Vdata
                 */
                VSdetach(vd_id);
            }
        }
    }
    free(begin);
    free(extent);
    /*
     *  If we got here, all was well
     */
    return 0;
}

void chk_sea_grid(int32 vd_id)
/*******************************************************************

   chk_sea_grid

   purpose: Check the SEAGrid Vdata in Vgroup 'Level-3 Binned Data'
        to make sure its constant values are correct

   Returns type: none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             vd_id            I      ID of the Vdata

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       25-Apr-1995     Original development
      W. Robinson       1-Sep-2009      They finally fixed the value of 
                                        l3_org.bins_eq to be # bins around 
                                        equator - so remove fix in this check

 *******************************************************************/ {
    int32 n_vd_rec;

    /*
     *  This structure / union deserves some explaination.
     *  the Vdata record gets read in as a contiguous set of 3
     *  floats and 4 doubles.  the buffer for VSread should be a
     *  uint8, so it is 1 part of the union.  The sea_grid_str
     *  holds the actual values.  The array align makes sure that
     *  this whole thing starts at a double alignment.  The
     *  extra int32 dummy will assure that the float64 data
     *  starts at an aligned double boundary.  When reading in,
     *  the read starts at the location buf + 4 so that the data starts
     *  at the right place.
     */
    union sea_gr_u {

        struct sea_grid_str {
            int32 dummy;
            int32 registration;
            int32 straddle;
            int32 bins;
            float64 radius;
            float64 max_north;
            float64 max_south;
            float64 seam_lon;
        } str;
        uint8 buf[48];
        double align[6];
    } sea_grid;

    /*
     *  set to read the first record and all fields
     */
    if (VSsetfields(vd_id,
            "registration,straddle,bins,radius,max_north,max_south,seam_lon"
            ) == -1) {
        printf(
                "*************** Program Error: For Vdata: 'SEAGrid',\n"
                "Cannot do VSsetfields\n");
        fmt_status = fmt_status | 2;
        return;
    }

    if (VSseek(vd_id, 0) == -1) {
        printf(
                "*************** Program Error: For Vdata: 'SEAGrid',\n"
                "Cannot do VSseek\n");
        fmt_status = fmt_status | 2;
    }

    if ((n_vd_rec = VSread(vd_id, (sea_grid.buf + 4), 1, 0)) != 1) {
        printf(
                "**** For Vdata: 'SEAGrid', data read failed to read the 1 record \n"
                "      with return: %d\n", n_vd_rec);
        fmt_status = fmt_status | 2;
    } else {
        if (verbose == 1)
            printf("For Vdata 'SEAGrid', registration = %d\n",
                sea_grid.str.registration);
        if (sea_grid.str.registration != 5) {
            printf("**** For Vdata 'SEAGrid', mismatch in registration\n");
            printf("     expected: 5, read: %d\n", sea_grid.str.registration);
            fmt_status = fmt_status | 2;
        }
        if (verbose == 1)
            printf("For Vdata 'SEAGrid', straddle = %d\n",
                sea_grid.str.straddle);
        if (sea_grid.str.straddle != 0) {
            printf("**** For Vdata 'SEAGrid', mismatch in straddle\n");
            printf("     expected: 0, read: %d\n", sea_grid.str.straddle);
            fmt_status = fmt_status | 2;
        }
        if (verbose == 1)
            printf("For Vdata 'SEAGrid', bins = %d\n",
                sea_grid.str.bins);
        /*  The # bins around equator  was finally corrected in reprocessing 09
            if( sea_grid.str.bins != l3_org.numrows )
         */
        if (sea_grid.str.bins != l3_org.bins_eq) {
            printf("**** For Vdata 'SEAGrid', mismatch in bins\n");
            printf("     expected: %d, read: %d\n", l3_org.bins_eq,
                    sea_grid.str.bins);
            fmt_status = fmt_status | 2;
        }
        if (verbose == 1)
            printf("For Vdata 'SEAGrid', radius = %f\n",
                sea_grid.str.radius);
        if (sea_grid.str.radius != 6378.137) {
            printf("**** For Vdata 'SEAGrid', mismatch in radius\n");
            printf("     expected: 6378.137, read: %f\n", sea_grid.str.radius);
            fmt_status = fmt_status | 2;
        }
        if (verbose == 1)
            printf("For Vdata 'SEAGrid', max_north = %f\n",
                sea_grid.str.max_north);
        if (sea_grid.str.max_north != 90.0) {
            printf("**** For Vdata 'SEAGrid', mismatch in max_north/n");
            printf("     expected: 90.0, read: %f\n", sea_grid.str.max_north);
            fmt_status = fmt_status | 2;
        }
        if (verbose == 1)
            printf("For Vdata 'SEAGrid', max_south = %f\n",
                sea_grid.str.max_south);
        if (sea_grid.str.max_south != -90.0) {
            printf("**** For Vdata 'SEAGrid', mismatch in max_south/n");
            printf("     expected: -90.0, read: %f\n", sea_grid.str.max_south);
            fmt_status = fmt_status | 2;
        }
        if (verbose == 1)
            printf("For Vdata 'SEAGrid', seam_lon = %f\n",
                sea_grid.str.seam_lon);
        if (sea_grid.str.seam_lon != -180.0) {
            printf("**** For Vdata 'SEAGrid', mismatch in seam_lon/n");
            printf("     expected: -180.0, read: %f\n", sea_grid.str.seam_lon);
            fmt_status = fmt_status | 2;
        }
    }
}

void chk_bin_index(int32 vd_id, int32 *begin, int32 *extent)
/*******************************************************************

   chk_bin_index

   purpose: Check the BinIndex Vdata in Vgroup 'Level-3 Binned Data'
        to make sure its values are correct.  The values in fields 
        begin and extent will not be checked in this routine.  

   Returns type: none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             vd_id            I      ID of the Vdata
      int32*		begin		 O      # of 1st data containing bin
      int32*		extent		 O      # of bins stored in eh. row

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      L. Kumar   	23-Feb-1996     Original development
      W. Robinson       25 Jul 1996     correct hsize check so it fails on
                                        difference > 1/1000

 *******************************************************************/ {
    int32 i, status, row;
    int32 *int32buf;
    float64 *ft64buf, vsz;

    /*
     *  set up space needed
     */
    int32buf = (int32 *) malloc(l3_org.numrows * sizeof (int32));
    ft64buf = (float64 *) malloc(l3_org.numrows * sizeof (float64));

    vsz = l3_org.vsize; /* Note that all values pre-computed in l3_get_org */

    /*
     *  read and verify the row numbers from vdata BinIndex
     */
    if ((status = (rdvdata(vd_id, "row_num", 0, l3_org.numrows,
            (unsigned char *) &int32buf[0]))) >= 0) {
        if (verbose == 1)
            printf("For Vdata 'BinIndex', checking row_num \n");
        for (i = 0; i < l3_org.numrows; i++)
            if (int32buf[i] != i) {
                printf("***** For Vdata 'BinIndex', mismatch in row_num\n");
                printf("      expected: %d, read: %d\n", i, int32buf[i]);
                fmt_status = fmt_status | 2;
            }
    }

    /*
     *  read and verify the vsize from vdata BinIndex
     */
    if ((status = (rdvdata(vd_id, "vsize", 0, l3_org.numrows,
            (unsigned char *) ft64buf))) >= 0) {
        if (verbose == 1)
            printf("For Vdata 'BinIndex', checking vsize \n");
        for (i = 0; i < l3_org.numrows; i++)
            if (ft64buf[i] != vsz) {
                printf("***** For Vdata 'BinIndex', mismatch in vsize\n");
                printf("      expected: %15.10f, read: %15.10f\n", vsz, ft64buf[i]);
                fmt_status = fmt_status | 2;
            }
    }

    /*
     *  read and verify the hsize from vdata BinIndex
     */
    if ((status = (rdvdata(vd_id, "hsize", 0, l3_org.numrows,
            (unsigned char *) ft64buf))) >= 0) {
        if (verbose == 1)
            printf("For Vdata 'BinIndex', checking hsize \n");
        for (row = 0; row < l3_org.numrows; row++) {
            if ((fabs(ft64buf[row] - l3_org.hsize[row])) > 1.0 / 1000.0) {
                printf("***** For Vdata 'BinIndex', mismatch in hsize\n");
                printf("      expected: %10.8g, read: %10.8g\n", l3_org.hsize[row],
                        ft64buf[row]);
                printf("      for index %d\n", row);
                fmt_status = fmt_status | 2;
            }
        }
    }

    /*
     *  read and verify the start_num from vdata BinIndex
     */
    if ((status = (rdvdata(vd_id, "start_num", 0, l3_org.numrows,
            (unsigned char *) int32buf))) >= 0) {
        if (verbose == 1)
            printf("For Vdata 'BinIndex', checking start_num \n");
        for (i = 0; i < l3_org.numrows; i++) {
            if (int32buf[i] != l3_org.start_bin[i]) {
                printf("***** For Vdata 'BinIndex', mismatch in start_num\n");
                printf("      expected: %d, read: %d\n", l3_org.start_bin[i + 1],
                        int32buf[i]);
                fmt_status = fmt_status | 2;
            }
        }
    }

    /*
     *  read and verify the "max" (maximum no. of bins/row) from vdata BinIndex
     */
    if ((status = (rdvdata(vd_id, "max", 0, l3_org.numrows,
            (unsigned char *) int32buf))) >= 0) {
        if (verbose == 1)
            printf("For Vdata 'BinIndex', checking max \n");
        for (i = 0; i < l3_org.numrows; i++) {
            if (int32buf[i] != l3_org.max_bin[i]) {
                printf("***** For Vdata 'BinIndex', mismatch in max\n");
                printf("      expected: %d, read: %d\n", l3_org.max_bin[i],
                        int32buf[i]);
                fmt_status = fmt_status | 2;
            }
        }
    }

    /*
     *  read and verify the begin from vdata BinIndex
     */
    rdvdata(vd_id, "begin", 0, l3_org.numrows, (unsigned char *) begin);

    /*
     *  read and verify the extent from vdata BinIndex
     */
    rdvdata(vd_id, "extent", 0, l3_org.numrows, (unsigned char *) extent);

    free(int32buf);
    free(ft64buf);
}

void chk_bin_list(int32 vd_id, int32 nbins, int32 *begin, int32 *extent)
/*******************************************************************

   chk_bin_list

   purpose: Read the BinList Vdata in Vgroup 'Level-3 Binned Data' for
        field 'binno' to make sure begin and extent values stored in Vdata
        'BinIndex' are correct.  Also, verifies whether global attribute
        'Data Bins' is correct.

   Returns type: none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             vd_id            I      ID of the Vdata
      int32             nbins            I      # of bins
      int32*            begin            I      # of 1st data containing bin
      int32*            extent           I      # of bins stored in eh. row

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      L. Kumar          07-Mar-1996     Original development

 *******************************************************************/ {
    int32 i, binno[512], row = 0, rowcount = 0;
    int32 *loc_extent, *loc_begin, nrec;
    int32 start, nelts, rec, recsleft, lastbin;
    static int32 prev_bin = 0;

    /*
     *  set up storage needed
     */
    loc_extent = (int32 *) calloc(l3_org.numrows, sizeof ( int32));
    loc_begin = (int32 *) calloc(l3_org.numrows, sizeof ( int32));

    lastbin = l3_org.start_bin[0] + l3_org.max_bin[0] - 1;
    /*
     *  set # of records to read = 512
     */
    nelts = nrec = 512;
    /*
     *  Read binno field of vdata 'BinList', 512 bins at a time, 
     *  create begin and extent and compare the begin and extent with
     *  the input begin and extent
     */

    printf("\n nbins = %d", nbins);
    for (rec = 0; rec < nbins; rec += 512) {
        if ((recsleft = nbins - rec) < 512)
            nelts = nrec = recsleft;
        start = rec;
        if ((rdvdata(vd_id, "bin_num", start, nelts,
                (unsigned char *) binno)) >= 0) {
            /*
             * Check the bin #s to see they are in ascending order
             */
            if (prev_bin >= binno[0]) {
                printf("****For Vdata 'BinList', discrepancy in bin numbers\n");
                printf("      expected bin# greater than %d, but read: %d\n",
                        prev_bin, binno[0]);
                fmt_status = fmt_status | 2;
            }
            /*
             *  Generate begins and extents for the bin numbers read
             */
            for (i = 0; i < nrec; i++) {
                if (i < nrec - 1 && binno[i] > binno[i + 1]) {
                    printf("***For Vdata 'BinList', discrepancy in bin numbers");
                    printf("      expected bin# greater than %d, but read: %d\n",
                            prev_bin, binno[0]);
                    fmt_status = fmt_status | 2;
                }
                if (binno[i] > lastbin) { /* next row */
                    loc_extent[row] = rowcount;
                    while (binno[i] > lastbin) {
                        row++;
                        lastbin = l3_org.start_bin[row] + l3_org.max_bin[row] - 1;
                    }
                    rowcount = 1;
                    loc_begin[row] = binno[i];
                }
                else { /* add bin to rowcount */
                    rowcount++;
                }
            }
            prev_bin = binno[nrec - 1];
        }
    }

    /*  don't forget to count the last [partial] row of bins  */
    loc_extent[row] = rowcount;

    if (verbose == 1) {
        printf("\nFor Vdata 'BinIndex', checking begin ");
        printf("\nFor Vdata 'BinIndex', checking extent \n");
    }
    for (i = 0; i < l3_org.numrows; i++) {
        /*
                printf("\n[%d]\t%d\t%d", i, loc_begin[i+1], begin[i]);
                printf("\n[%d]\t%d\t%d", i, loc_extent[i+1], extent[i]);
         */
        if (loc_begin[i] != begin[i]) {
            printf("***** For Vdata 'BinIndex', mismatch in 'begin'\n");
            printf("      expected: %d for row %d, read: %d\n",
                    loc_begin[i], i, begin[i]);
            fmt_status = fmt_status | 2;
        }
        if (loc_extent[i] != extent[i]) {
            printf("***** For Vdata 'BinIndex', mismatch in 'extent'\n");
            printf("      expected: %d for row %d, read: %d\n",
                    loc_extent[i], i, extent[i]);
            fmt_status = fmt_status | 2;
        }
    }
    free(loc_extent);
    free(loc_begin);
}

void chk_l3_prod(int32 vd_id, int32 nbins, char *name, char *field_list)
/*******************************************************************

   ck_v_l3_v2p6

   purpose: Check the SEAGrid Vdata in Vgroup 'Level-3 Binned Data'
        to make sure its constant values are correct

   Returns type: none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             vd_id            I      ID of the Vdata
      int32             nbins            I      # of bins
      char *            name             I      name of Vdata
      char *            field_list       I      field list in Vdata

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       25-Apr-1995     Original development
      W. Robinson, SAIC 5 Apr 2005      adapt to check any product

 *******************************************************************/ {
    int32 n_vd_rec, i;
    float32 *p_data;
    float32 sum, sum_sq;

    /*
     *  make space for the sum and sum ^2 (such as they are)
     */
    p_data = (float32 *) malloc(2 * nbins * sizeof ( float32));

    /*
     *  set to read the sum and sum sq as INTERLACED
     */
    if (VSsetfields(vd_id, field_list) == -1) {
        printf(
                "*************** Program Error: For Vdata: '%s',\nCannot do VSsetfields\n",
                name);
        return;
    }

    if (VSseek(vd_id, 0) == -1) {
        printf(
                "*************** Program Error: For Vdata: '%s',\nCannot do VSseek\n",
                name);
    }
    if ((n_vd_rec = VSread(vd_id, (uint8 *) p_data, nbins, FULL_INTERLACE))
            != nbins) {
        printf(
                "**** For Vdata: '%s', data read failed to read the %d record \n"
                "      with return: %d\n", name, nbins, n_vd_rec);
    } else {
        /*
         *  get the mean of both ~sum and ~sum_sq
         */
        sum = 0;
        sum_sq = 0;
        for (i = 0; i < nbins; i++) {
            /*  this debug printout is voluminous
                     printf( "i: %d,  ~sum: %f,  ~sum_sq: %f\n", i,
             *( p_data + 2 * i ), *( p_data + 2 * i + 1 ) );
             */
            sum += *(p_data + 2 * i);
            sum_sq += *(p_data + 2 * i + 1);
        }
        sum = sum / nbins;
        sum_sq = sum_sq / nbins;

        if (verbose == 1)
            printf("For Vdata: '%s', mean of sum: %f, sum_sq: %f\n",
                name, sum, sum_sq);
    }

    free(p_data);
}


