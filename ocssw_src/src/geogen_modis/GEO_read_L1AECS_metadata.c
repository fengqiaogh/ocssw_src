#include "PGS_MODIS_35251.h"
#include "GEO_input.h"
#include "GEO_product.h"
#include "L1a_data.h"
#include "smfio.h"

int  GEO_read_L1AECS_metadata(
               MODFILE * const l1a_file,
               ECS_metadata_struct * const ECS_metadata 
              )

/* 
!C**************************************************************************
!Description:   
        Routine in the input group of the Level-1A geolocation
        software to read ECS Inventory metadata from the 
        L1A product file.  

        Reads ECS Inventory Metadata:
                        RANGEBEGINNINGDATE
                        RANGEBEGINNINGTIME
                        RANGEENDINGDATE
                        RANGEENDINGTIME
                        DAYNIGHTFLAG
                        LOCALGRANULEID
                        PARAMETERVALUE.1 ("GRANULENUMBER")
                        ASSOCIATEDPLATFORMSHORTNAME.1
                        PRODUCTIONHISTORY

!Input Parameters:
        MODFILE *l1a_file - MAPI structure for Level 1A file

!Output Parameters:
                ECS_metadata - structure for ECS metadata inputs

Returns:
                SUCCEED         if all metadatea are retrieved
                FAIL            otherwise.

Global variables:
        None. 

Called by:
                GEO_prepare_l1a_data()

Routines Called:
                getMODISECSinfo - Reads an ECS metadata attribute from an
                        HDF global attribute.
                modsmf - writes error status messages to log

!Revision History:
 * $Log: GEO_read_L1AECS_metadata.c,v $
 * Revision 4.1  2003/02/21 22:49:42  kuyper
 * Corrected to use void* pointers with %p format code.
 *
 * Revision 3.1  2002/06/13 22:51:29  kuyper
 * Removed unnecessary NCSA acknowledgement.
 *
 * Revision 2.4  2001/01/17 16:06:37  vlin
 * PRODUCTIONHISTORY is added per CCR 507
 *
 * Revision 2.3  2000/08/12  23:02:24  fhliang
 * Added AssociatedPlatformShortName to the metadata list.
 * Changed return type to 'int'.
 * Cleaned up messages.
 * Changed 'rtval' to 'retval'.
 * Added static char filefunc and used it in calls to modsmf().
 * Removed an extra PDL paragraph.
 *
 * Revision 2.2  1998/03/04  00:48:51  jjb
 * Bug Fix MODxl00638: Initialized ECS_metadata ingest strings
 * 	to prevent writing uninitialized metadata strings
 * 	into the MOD03 product's metadata.
 *
 * Revision 2.1  1997/10/21  18:16:22  kuyper
 * Returned from ClearCase
 *
 * Revision /main/GEO_V2_DEV/3 1997/10/06 kuyper
 * Changed MECS_CORE to L1A_COREMETADATA.
 *
 * Revision /main/GEO_V2_DEV/2  1997/09/3  Ding
 * Modified the code to read V2 input data.
 * ding@ltpmail.gsfc.nasa.gov
 *
 * Revision 1.6.1.1  1997/07/21  22:20:17  kuyper
 * Merged in out-of-sequence changes.
 *
 * Revision 1.6  1997/07/21  16:24:34  kuyper
 * Baselined Version 1
 *
 *Parallel development:
 * Revision 1.5  1997/06/02  17:51:51  kuyper
 * Merged seed files.
 *
 * Revision 1.5  1997/03/26  18:12:56  fhliang
 * Initial revision of SDST delivery of GEO_read_L1AECS_metadata.c.
 *
	Revision 1.4  1997/03/11  22:05:34  fshaw
	*** empty log message ***
	
	Revision 1.3  1997/02/13  23:25:11  kuyper
	Fixed typo in error message.
	
        Revision 1.2  1997/01/13 19:41:15  kuyper
        Adjusted #include file list.
	James Kuyper Jr. <kuyper@ltpmail.gsfc.nasa.gov>

        Revision 1.1  1997/01/08 16:06:27  mikej
        Initial revision


Requirements:
                PR03-F-2.1-1
                PR03-F-2.1-2
                PR03-F-2.1-3


!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.
                
!END***************************************************************************
*/

{

#define N_ELEMENTS(foo) ((int) (sizeof(foo) / sizeof(foo[0])))
#define STRMAX 120
#define ARGMAX 80
    const char ECSINITSTRING[] = "NOT SET";  /* metadata output initialization*/

    /*
    Warning: The order and sequence of the initilizations can not
    be changed. We finish the initilization later and count on this order 
    */
    struct ECSmeta{
        char *attrNameStr; /* name of metadata */
        char *hdfAttrname; /* name of HDF global attribute */
        char datatype[DATATYPELENMAX]; /* metadata data type */
        long int n_elements; /* Number of metadata elements to retrieve */
        void * data; /* pointer to the data */
    } ECS_input_metadata[] = {
        {CORE_RANGE_BEG_DATE, L1A_COREMETADATA, TXT,
            (long int) sizeof(ECS_metadata->rangebeginningdate), NULL},
        {CORE_RANGE_BEG_TIME, L1A_COREMETADATA, TXT, 
            (long int) sizeof(ECS_metadata->rangebeginningtime), NULL},
        {CORE_RANGE_ENDING_DATE, L1A_COREMETADATA, TXT,
            (long int) sizeof(ECS_metadata->rangeendingdate), NULL},
        {CORE_RANGE_ENDING_TIME, L1A_COREMETADATA, TXT,
            (long int) sizeof(ECS_metadata->rangeendingtime), NULL},
        {CORE_DAYNIGHTFLAG, L1A_COREMETADATA, TXT,
            (long int) sizeof(ECS_metadata->operationmode), NULL},
        {CORE_LOCALGRANULEID, L1A_COREMETADATA, TXT,
            (long int) sizeof(ECS_metadata->localinputgranuleid), NULL},
        {CORE_PARAMETERVALUE ".1", L1A_COREMETADATA, TXT,
            (long int) sizeof(ECS_metadata->granulenumber), NULL},
        {CORE_ASSOCIATEDPLATFORMSHORTNAME ".1", L1A_COREMETADATA, TXT,
            (long int) sizeof(ECS_metadata->platformshortname), NULL},
        {MECS_PRODHISTORY, ARCHIVEMETADATA, TXT, (long int) sizeof( \
            (ECS_metadata->version_metadata).productionhistory), NULL}
    };

    int i, retval=SUCCEED;  /* loop index */
    int sameattr;
    char msgbuf[STRMAX];    /* scratch string buffer */

    static char filefunc[] = __FILE__", GEO_read_L1AECS_metadata";


  if((l1a_file == NULL) || (ECS_metadata == NULL) )
  {
      sprintf(msgbuf, " l1a_file = %p, ECS_metadata = %p",
	  (void*)l1a_file, (void*)ECS_metadata);
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

      return FAIL;
  } 
  
/*  Initialize ECS_metadata output */
    (void)memset(ECS_metadata,0,sizeof(ECS_metadata_struct));
    (void)strncpy(ECS_metadata->rangebeginningdate,ECSINITSTRING,
	   sizeof(ECS_metadata->rangebeginningdate)-1);
    (void)strncpy(ECS_metadata->rangebeginningtime,ECSINITSTRING,
	   sizeof(ECS_metadata->rangebeginningtime)-1);
    (void)strncpy(ECS_metadata->rangeendingdate,ECSINITSTRING,
	   sizeof(ECS_metadata->rangeendingdate)-1);
    (void)strncpy(ECS_metadata->rangeendingtime,ECSINITSTRING,
	   sizeof(ECS_metadata->rangeendingtime)-1);
    (void)strncpy(ECS_metadata->operationmode,ECSINITSTRING,
	   sizeof(ECS_metadata->operationmode)-1);
    (void)strncpy(ECS_metadata->localinputgranuleid,ECSINITSTRING,
	   sizeof(ECS_metadata->localinputgranuleid)-1);
    (void)strncpy(ECS_metadata->granulenumber,ECSINITSTRING,
	   sizeof(ECS_metadata->granulenumber)-1);
    (void)strncpy(ECS_metadata->platformshortname,ECSINITSTRING,
	   sizeof(ECS_metadata->platformshortname)-1);
    (void)strncpy((ECS_metadata->version_metadata).productionhistory,ECSINITSTRING,\
           sizeof((ECS_metadata->version_metadata).productionhistory)-1);

/*  Finish initializing the ECS_input_metadata array with the contents 
    defined above.
    
    See warning about sequence in declaration section 
*/
    ECS_input_metadata[0].data = &ECS_metadata->rangebeginningdate;
    ECS_input_metadata[1].data = &ECS_metadata->rangebeginningtime;
    ECS_input_metadata[2].data = &ECS_metadata->rangeendingdate;
    ECS_input_metadata[3].data = &ECS_metadata->rangeendingtime;
    ECS_input_metadata[4].data = &ECS_metadata->operationmode;
    ECS_input_metadata[5].data = &ECS_metadata->localinputgranuleid;
    ECS_input_metadata[6].data = &ECS_metadata->granulenumber;
    ECS_input_metadata[7].data = &ECS_metadata->platformshortname;
    ECS_input_metadata[8].data = &ECS_metadata->version_metadata.productionhistory;

    for (i = 0; i < N_ELEMENTS(ECS_input_metadata); i++)
    {
        sameattr = (i > 0 && !(strncmp(ECS_input_metadata[i].hdfAttrname,
                           ECS_input_metadata[i-1].hdfAttrname,
                          sizeof(ECS_input_metadata[i].hdfAttrname))));
        if (getMODISECSinfo(l1a_file, (sameattr ? NULL :
                   ECS_input_metadata[i].hdfAttrname), 
            ECS_input_metadata[i].attrNameStr, ECS_input_metadata[i].datatype,
            &ECS_input_metadata[i].n_elements, ECS_input_metadata[i].data) != MAPIOK && 
        strcmp(ECS_input_metadata[i].attrNameStr, MECS_PRODHISTORY) != 0) {
            sprintf(msgbuf, "getMODISECSinfo(%.*s)", ARGMAX,
	      ECS_input_metadata[i].attrNameStr);
            modsmf(MODIS_E_GEO, msgbuf, filefunc);

            retval = FAIL;
        }
    }
    
    return retval;
}
