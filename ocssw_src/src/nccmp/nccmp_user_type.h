/*
 * nccmp_user_type.h
 *
 *  Created on: Sep 27, 2016
 *      Author: rhealy
 */

#ifndef SRC_NCCMP_NCCMP_USER_TYPE_H_
#define SRC_NCCMP_NCCMP_USER_TYPE_H_
#include <stdint.h>
#define NCCMP_MAX_COMPOUND_FIELD_DIMS 4

typedef struct nccmp_user_type_t {
    nc_type base_type; /* Used by vlen and enum types for their primitive type. */
    int dim_sizes[NCCMP_MAX_COMPOUND_FIELD_DIMS];
    char field_index; /* Relative to parent compound. */
    int group_id; /* Group that defined the type. */
    int id; /* Unique instance id of a type in tree defined by nccmp. */
    char name[NC_MAX_NAME]; /* Base name like "myField". */
    int num_dims;
    size_t num_fields; /* Number of of fields in the compound type. */
    size_t offset; /* Byte offset within a compound. */
    struct nccmp_user_type_t *fields;
    size_t root_size; /* Total outermost root compound size. */
    size_t size; /* Number of bytes. */
    char *tree_name; /* Dotted name from root to child, like "myVar.myField1.myChild2". */
    nc_type type_id; /* Defined by file metadata. */
    nc_type user_class; /* Class like NC_BYTE or NC_COMPOUND. */
} nccmp_user_type_t;


#endif /* SRC_NCCMP_NCCMP_USER_TYPE_H_ */
