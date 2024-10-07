#ifndef FIND_VARIABLE_H
#define FIND_VARIABLE_H
#ifdef __cplusplus
extern "C" {
#endif
#include "stdlib.h"
#include "netcdf.h"

/**
 * @brief Recursively searches for a netcdf variable within an NC file or group. Search is recursive ( DFS tree traversal)
 * @param var_name variable name
 * @param ncid file or group id
 * @param netcdf_c_var_id output variable nc id
 * @param netcdf_c_grp_id output parent group nc id. Can be the file root
 * @return
 */
int32_t find_nc_variable_parent_grp_c_interface_id_ncid(const char* var_name, int32_t ncid,
                                                        int32_t* netcdf_c_var_id, int32_t* netcdf_c_grp_id);
int32_t find_nc_variable_possible_names(const char * possible_names[], int32_t number_of_names, int32_t ncid,
                                        int32_t* netcdf_c_var_id, int32_t* netcdf_c_grp_id, int32_t * file_index_name );
#ifdef __cplusplus
}
#endif

#endif
