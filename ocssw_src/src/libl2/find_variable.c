#include "find_variable.h"
#include "stdio.h"
#define MAXNUMNCGRP 100
int32_t find_nc_variable_parent_grp_c_interface_id_ncid(const char* var_name, int32_t ncid,
                                                        int32_t* netcdf_c_var_id, int32_t* netcdf_c_grp_id) {
    int32_t numgrps;
    int32_t ncids[MAXNUMNCGRP];
    int32_t status = nc_inq_varid(ncid, var_name, netcdf_c_var_id);
    if (status == NC_NOERR) {
        *netcdf_c_grp_id = ncid;
        return status;
    }

    status = nc_inq_grps(ncid, &numgrps, ncids);
    if (status != NC_NOERR)
        return status;
    else {
        for (int32_t i = 0; i < numgrps; i++) {
            status = find_nc_variable_parent_grp_c_interface_id_ncid(var_name, ncids[i], netcdf_c_var_id,
                                                                     netcdf_c_grp_id);
            if (status == NC_NOERR)
                return status;
        }
    }
    return NC_ENOTVAR;
};

int32_t find_nc_variable_possible_names(const char* possible_names[], int32_t number_of_names, int32_t ncid,
                                        int32_t* netcdf_c_var_id, int32_t* netcdf_c_grp_id,
                                        int32_t* file_index_name) {
    int32_t status = 0;
    if (number_of_names <= 0)
        return NC2_ERR;
    for (int32_t i = 0; i < number_of_names; i++) {
        status = find_nc_variable_parent_grp_c_interface_id_ncid(possible_names[i], ncid, netcdf_c_var_id,
                                                                 netcdf_c_grp_id);
        if (status == NC_NOERR) {
            *file_index_name = i;
            return status;
        }
    }
    return NC2_ERR;
}
