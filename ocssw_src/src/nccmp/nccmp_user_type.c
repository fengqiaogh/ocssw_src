/*
 * nccmp_user_type.c
 *
 *  Created on: Sep 26, 2016
 *      Author: rhealy
 */
#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>
#include "nccmp_user_type.h"
//#define NCCMP_MAX_COMPOUND_FIELD_DIMS 4
//#define NCCMP_MAX_COMPOUND_FIELDS 20
// see https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-c/nc_005finq_005fcompound.html

nccmp_user_type_t* nccmp_load_group_usertype_array(int group_id, int *nuser_types) {
    int status;
    int i, typeid, *typeids = NULL, field_i;
    nccmp_user_type_t *comp_node;

    status = nc_inq_typeids(group_id, nuser_types, NULL);
    if (status != NC_NOERR) {
        fprintf(stderr, "Netcdf error status = %d\n", status);
        exit(status);
    }

    if (*nuser_types < 1) {
        *nuser_types = 0;
        return (NULL);
    }

    typeids = (int *) calloc(*nuser_types, sizeof (int));
    if (!typeids) {
        fprintf(stderr, "Error allocating typeids\n ");
        exit(1);
    }

    comp_node = (nccmp_user_type_t *) malloc((*nuser_types) * sizeof (nccmp_user_type_t));
    status = nc_inq_typeids(group_id, NULL, typeids);
    if (status != NC_NOERR) {
        fprintf(stderr, "Netcdf error status = %d\n", status);
        exit(status);
    }

    //results = nccmp_darray_create(*nuser_types);
    //using nc_inq_user_type for possible future implementation of other user types - rjh
    for (i = 0; i < *nuser_types; ++i) {
        typeid = typeids[i];
        comp_node[i].type_id = typeid;
        comp_node[i].group_id = group_id;
        status = nc_inq_user_type(group_id,
                typeid,
                comp_node[i].name,
                & comp_node[i].size,
                & comp_node[i].base_type,
                & comp_node[i].num_fields,
                & comp_node[i].user_class);
        if (status != NC_NOERR) {
            fprintf(stderr, "Netcdf error status = %d\n", status);
            exit(status);
        }
        if (comp_node[i].user_class != NC_COMPOUND) {
            fprintf(stderr, "Can only handle compound types at the moment\n ");
            exit(1);
        }

        if (comp_node[i].num_fields) {
            comp_node[i].root_size = comp_node[i].size; /* This node is a compound root. */
            comp_node[i].fields = (nccmp_user_type_t *) malloc(comp_node[i].num_fields * sizeof (nccmp_user_type_t));

            for (field_i = 0; field_i < comp_node[i].num_fields; ++field_i) {
                comp_node[i].fields[field_i].group_id = group_id;
                comp_node[i].fields[field_i].field_index = field_i;

                status = nc_inq_compound_field(group_id,
                        typeid,
                        field_i,
                        comp_node[i].fields[field_i].name,
                        & comp_node[i].fields[field_i].offset,
                        & comp_node[i].fields[field_i].type_id,
                        & comp_node[i].fields[field_i].num_dims,
                        comp_node[i].fields[field_i].dim_sizes);
                if (status != NC_NOERR) {
                    fprintf(stderr, "Netcdf error status = %d\n", status);
                    exit(status);
                }

            }
        }
    }

    free(typeids);

    return comp_node;
}
