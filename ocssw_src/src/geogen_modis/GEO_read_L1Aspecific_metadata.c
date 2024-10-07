/* File: GEO_read_L1Aspecific_metadata.c */

#include "PGS_MODIS_35251.h"
#include "GEO_input.h"
#include "L1a_data.h"
#include "smfio.h"
#include "granule_extract_metadata.h"

int GEO_read_L1Aspecific_metadata(
	MODFILE             * const l1a_file, 
        l1a_metadata_struct * const granule_metadata,
	int 		    * const number_of_scans
	)

/* 
!C**************************************************************************
 
!Description:   
		Routine in Input group of the Level-1A geolocation
                software to read L1A specific granule metadata from
		the input product.  These metadata are retrieved from
		individual global HDF attributes using MAPI calls.
		
!Input Parameters:
                MODFILE *l1a_file - the MAPI structure for the L1A product

!Output Parameters:
		granule_metadata - l1a specific granule metadata
		number_of_scans - the number of scans in the granule

Return parameter:
                SUCCESS if all the specified specific granule metadata are
			read from the l1a_file,
		FAIL otherwise.

Global variables: None

Call functions:
		int getMODISfileinfo - MAPI routine to read a global 
		  metadata attribute from an HDF file.
                modsmf(MODIS_X_MNEMONIC_STRING, "user message string", 
		  "function, GEO_write_parameters.c");


!Revision History:
 		$Log: GEO_read_L1Aspecific_metadata.c,v $
 		Revision 4.1  2003/02/21 22:44:48  kuyper
 		Corrected to use void* pointers with %p format code.

 		Revision 3.1  2002/06/13 22:54:16  kuyper
 		Removed unnecessary NCSA acknowledgement.

 		Revision 2.1  1997/10/21 18:16:22  kuyper
 		Returned from ClearCase

 * Revision /main/GEO_V2_DEV/1  1997/09/4  Ding
 * Revised to read Version 2 data.
 * ding@ltpmail.gsfc.nasa.gov
 *
 * Revision 1.7  1997/07/21  16:24:34  kuyper
 * Baselined Version 1
 *
 * Revision 1.7  1997/03/26  18:13:50  fhliang
 * Initial revision of SDST delivery of GEO_read_L1Aspecific_metadata.c.
 *
 		Revision 1.6  1997/02/13 19:53:38  kuyper
 		Merged seed files.

 		Revision 1.5  1997/01/30 20:20:50  kuyper
 		Backed out temporary changes made to match
 		the errors in the simulated data.

		Revision 1.4  1997/01/14  22:38:28  kuyper
		Changed to use M01 M-API macros indirectly.
		Temporary change to acommodate incorrectly simulated L1A files.
		
 		Revision 1.3  1997/01/14 19:30:25  kuyper
 		Adjusted header files.
		James Kuyper - kuyper@ltpmail.gsfc.nasa.gov

 		Revision 1.2  1997/01/14 18:44:30  mikej
 		Cosmetic typo fixed in comment.

		Revision 1.1  1997/01/14  18:20:58  mikej
		Initial revision
		

Requirements:
		PR03-F-2.1-1
		PR03-F-2.1-2

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

                
!END***************************************************************************
*/

{
/*   Initialize an array of structures, attribute_object's:
        attribute_object.name         text string of the attribute's name.
        attribute_object.data_type    DATATYPELENMAX long string of the 
                        	          attribute's M-API data type.
	attribute_object.n_elements   number of elements in the metadata
        attribute_object.data	      pointer to where to (eventually) store 
					the data from the attribute.

  Initialize int32 data type num_scans_in
    to retrieve the number of
    scans into the correct
    data type

    attribute_object initialization:

    .name        .data_type  	.n_elements	 .data
                                                          
================================================================================
               |  	    |		       | 
"Number of scans"      I32	1		 number_of_scans
"Max Earth Frames"     I32	1		 granule_metadata.
						   max_earth_frames
"Max SD Frames"	       I32	1		 granule_metadata.
						   max_sd_frames_in
"Max SV Frames"	       I32	1		 granule_metadata.max_sv_frames
"Incomplete Scans"     I32	1		 granule_metadata.
						   incomplete_scans
"Missing_packets"      I32	1		 granule_metadata.
						   missing_packets
"Packets with bad CRC" I32	1		 granule_metadata.
						   bad_CRC_packets
"Discarded Packets"    I32	1	         granule_metadata.
						   discarded_packets
	Note: data types of elements in the granule_metadata structure
		match the data types of the attributes in the L1A file.
*/

#define N_ELEMENTS(foo) ((int) (sizeof(foo) / sizeof(foo[0])))
#define STRMAX 120
#define NAMEMAX (STRMAX - 40)

    /*
    Warning: The order and sequence of the initilizations can not
    be changed. We finish the initilization later and count on this order 
    */
    struct L1Ameta{
        char *name; /* name of attribute */
        char data_type[DATATYPELENMAX]; /* M-API data type */
        long int n_elements; /* Number of metadata elements to retrieve */
        void * data; /* pointer to the data */
    } attribute_object[] = {
        {NUMBER_OF_SCANS, I32, 1L, NULL},
        {MAX_EARTH_FRAMES, I32, 1L, NULL},
        {MAX_SD_FRAMES, I32, 1L, NULL},
        {MAX_SV_FRAMES, I32, 1L, NULL},
        {INCOMPL_SCANS, I32, 1L, NULL},
        {MISSING_PACKETS, I32, 1L, NULL},
        {PACKTS_BAD_CRC,  I32, 1L, NULL},
        {DISCARD_PACKETS, I32, 1L, NULL}
    };
    int32 number_of_scans_in;
    int i, rtval=SUCCESS;  /* loop index */
    char msgbuf[STRMAX];    /* scratch string buffer */
 
/*   IF (l1a_file         = NULL)
  OR (granule_metadata = NULL)
  OR (number_of_scans  = NULL)
    BEGIN
      CALL modsmf function to report error in LogStatus:
       "GEO_read_L1Aspecific_metadata, GEO_read_L1Aspecific_metadata.c():
        MODIS_E_BAD_INPUT_ARG:XXXXXXXXX"
        timestamp "Invalid input argument l1a_file = " l1a_file
		  ", granule_metadata = " granule_metadata
		  ", number_of_scans = " number_of_scans

      RETURN FAIL
    END
*/
  if((l1a_file == NULL) || (granule_metadata == NULL) || 
    (number_of_scans == NULL )){
    sprintf(msgbuf, " l1a_file = %p, granule_metadata = %p, "
	"number_of_scans = %p",
	(void*)l1a_file, (void*)granule_metadata, (void*)number_of_scans);
    modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf,
      "GEO_read_L1Aspecific_metadata, GEO_read_L1Aspecific_metadata.c");
    return FAIL;
  } 
  
/*   Finish initializing the attribute_object array with the contents 
	defined above. 
    
    See warning about sequence in declaration section  
*/
    attribute_object[0].data = &number_of_scans_in;
    attribute_object[1].data = &granule_metadata->max_earth_frames;
    attribute_object[2].data = &granule_metadata->max_sd_frames;
    attribute_object[3].data = &granule_metadata->max_sv_frames;
    attribute_object[4].data = &granule_metadata->incomplete_scans;
    attribute_object[5].data = &granule_metadata->missing_packets;
    attribute_object[6].data = &granule_metadata->bad_CRC_packets;
    attribute_object[7].data = &granule_metadata->discarded_packets;
    
/*   FOR each metadata object, i, in the attribute_object array  
    BEGIN
      IF call to getMODISfileinfo(l1a_file,
                                  attribute_object[i].name,
                                  attribute_object[i].data_type,
                                  attribute_object[i].n_elements,
                                  attribute_object[i].data)
          to read the geolocation specific granule metadata 
	  from an HDF global attribute
        does NOT return MAPIOK
        BEGIN
          CALL modsmf function to report error in LogStatus:
           "putMODISfileinfo, GEO_read_L1Aspecific_metadata.c():
            MODIS_E_GEO:XXXXXXXXX"
            timestamp "Error returned by function reading" 
		attribute_object[i].name
	  RETURN FAIL
        END
    END FOR each attribute_object
*/
    for (i = 0; i < N_ELEMENTS(attribute_object); i++) {
        if (getMODISfileinfo(l1a_file, attribute_object[i].name, 
            attribute_object[i].data_type, &attribute_object[i].n_elements, 
            attribute_object[i].data) != MAPIOK) {
            sprintf(msgbuf, " reading %.*s", NAMEMAX,
	      attribute_object[i].name);
            modsmf(MODIS_E_GEO, msgbuf,
	      "getMODISfileinfo, GEO_read_L1Aspecific_metadata.c");
            rtval = FAIL;
        }
    }

    /* Read Extract Metadata */
    read_extract_metadata(l1a_file->sd_id,
                          &granule_metadata->Extract_Pixel_Offset,
                          &granule_metadata->Extract_Pixel_Count,
                          &granule_metadata->Extract_Line_Offset,
                          &granule_metadata->Extract_Line_Count);

/*   Perform explicit type conversion of number_of_scans_in 
	into number_of_scans.
*/
    *number_of_scans = (int) number_of_scans_in;
    
    return rtval;

/* END GEO_read_L1Aspecific_metadata
*/
}
