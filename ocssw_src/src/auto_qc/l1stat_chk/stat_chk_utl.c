#include "l1stat.h"
#include "l1stat_proto.h"

#define STMAX 255

extern char err_msg[1024];
extern int32 stat_status; /* status of statistical checking: 0 all good,
        		   	   1 program problem, 2 statistical problem, 
		 	   	   3 both problems  */

/*-----------------------------------------------------------------------------
    Function:  rdattr

    Returns:   int32 (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        The function rdattr reads the requested global attribute

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     sdfid        I    ID req to access HDF SDS interface
        char  *   attr_name    I    attribute name
        void  *   buf         I/O   pointer to data buffer

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     11/07/94  Original development

----------------------------------------------------------------------------*/
int32
#ifdef PROTOTYPE
rdattr(int32 sdfid, char *attr_name, void *buf)
#else
rdattr(sdfid, attr_name, buf)
int32 sdfid;
char *attr_name;
void *buf;
#endif
{
    int32 attrnum;

    if ((attrnum = SDfindattr(sdfid, attr_name)) < 0) {
        sprintf(err_msg, "rdattr: Failure in SDfindattr while trying \
			to read %s", attr_name);
        stat_status = stat_status | 1;
        return FAIL;
    }

    if ((SDreadattr(sdfid, attrnum, buf)) < 0) {
        sprintf(err_msg, "rdattr: Failure in SDreadattr while trying to \
		read %s", attr_name);
        stat_status = stat_status | 1;
        return FAIL;
    }

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  rdslice

    Returns:   int32 (sdsid)
        The return code is a negative value if any error occurs, otherwise,
        returns sdsid.

    Description:
        The function rdslice reads requested slice of data from the
        given named dataset

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     sdfid         I    ID req to access HDF SDS interface
        char      *name         I   SDS name
        int32     *start        I   start data dimension
        int32     *edge         I   no. of values to be read
        void      *buf          O   SDS data buffer

    NOTE:

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     11/02/94  Original development
----------------------------------------------------------------------------*/
int32
#ifdef PROTOTYPE
rdslice(int32 sdfid, char *name, int32 *start, int32 *edge, void *buf)
#else
rdslice(sdfid, name, start, edge, buf)
int32 sdfid, *start, *edge;
char *name;
void *buf;
#endif
{
    int32 index, sdsid, rank, num_type, nattrs;
    char sdsname[STMAX];
    /*
      clock_t  val, val2;
     */
    if ((index = SDnametoindex(sdfid, name)) < 0) {
        sprintf(err_msg, "rdslice: Failed to locate sds \"%s\" ", name);
        return FAIL;
    }
    if ((sdsid = SDselect(sdfid, index)) < 0) {
        sprintf(err_msg, "rdslice: SDselect failed for sds \"%s\" ", name);
        return FAIL;
    }

    if (edge[0] == 0 && edge[1] == 0 && edge[2] == 0)
        if ((SDgetinfo(sdsid, sdsname, &rank, edge, &num_type, &nattrs)) < 0)
            return FAIL;
    /*
      val = clock();
     */
    if ((SDreaddata(sdsid, start, NULL, edge, buf)) < 0) {
        sprintf(err_msg,
                "rdslice: SDreaddata error while reading \"%s\" ", name);
        return FAIL;
    }
    /*
      val2 = clock();
      printf("\ntime took for %s = %d\n", name, val2-val);
     */
    SDendaccess(sdsid);
    return SUCCEED;
}
