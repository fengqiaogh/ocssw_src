/*-----------------------------------------------------------------------------
    Function:  set_dim_names 

    Returns:   int32 (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        Set the dimension names for an SD

    Parameters: (in calling order)
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int32     sds_id        I   id of the SD to be set 
        char      *d0           I   name for the first dimension
        char      *d1           I   name for the second dimension
        char      *d2           I   name for the third dimension


    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Gary Fu        SAIC GSC       05/13/98  Original development 
----------------------------------------------------------------------------*/

#include "hdf.h"
#include "mfhdf.h"

int32 set_dim_names(int32 sds_id, char *d0, char *d1, char *d2) {
    int32 dim_id, iret;

    dim_id = SDgetdimid(sds_id, 0);
    if ((iret = SDsetdimname(dim_id, d0)) < 0)
        return FAIL;

    if (d1 != NULL && *d1 != 0) {
        dim_id = SDgetdimid(sds_id, 1);
        if ((iret = SDsetdimname(dim_id, d1)) < 0)
            return FAIL;
    }
    if (d2 != NULL && *d2 != 0) {
        dim_id = SDgetdimid(sds_id, 2);
        if ((iret = SDsetdimname(dim_id, d2)) < 0)
            return FAIL;
    }

    return SUCCEED;
}

