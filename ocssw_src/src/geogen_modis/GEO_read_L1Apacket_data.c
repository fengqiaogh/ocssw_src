#include "smfio.h"
#include "GEO_input.h"
#include "L1a_data.h"
#include "PGS_MODIS_35251.h"

PGSt_SMF_status GEO_read_L1Apacket_data(
	MODFILE		*l1a_file,
	int const	number_of_scans,
	uint16		earth_encoder_times[MAX_SCAN_NUMBER][ENCODER_LENGTH],
	sc_ancil_struct	sc_ancillary_data[2],
	int16		view_sector_start[MAX_SCAN_NUMBER][SECTOR_LENGTH],
	uint16		FRside[MAX_SCAN_NUMBER][ELEC_SIDES],
	uint16		SAside[MAX_SCAN_NUMBER][ELEC_SIDES],
	uint16		ss_cp_mode[]
)

/*
!C************************************************************************

!Description:   
               Subroutine in the input group of the Level-1A geolocation
               software to read mirror and spacecraft ancillary data
               from the L1A product file using the MAPI.
               PR03-CSC-07

!Input Parameters:
               number_of_scans - number of scans in the granule
               l1a_file - MAPI structure for L1A file

!Output Parameters:
       	earth_encoder_times -  array of Earth Encoder Time Data 
       	sc_ancillary_data - array of Spacecraft ancillary data.
       	view_sector_start - array of View Sector Start Encoder 
			Count and Vernier Count words
               FRside - Indicates which side of formatter electronics is
                       turned on
               SAside - Indicates which sides (sic) of the scan assembly
                       electronics are turned on.

Return parameter:
		MODIS_E_BAD_INPUT_ARG	If any argument has an invalid value.
		MODIS_E_GEO		If any subroutine failed.
		PGS_S_SUCCESS		Otherwise

Externally Defined:
               ATTITUDE_ANGLE_ROLL     "L1A_data.h"
               ATTITUDE_ANGLE_PITCH    "L1A_data.h"
               ATTITUDE_ANGLE_YAW      "L1A_data.h"
               ATTITUDE_RATE_ROLL      "L1A_data.c"
               ATTITUDE_RATE_PITCH     "L1A_data.h"
               ATTITUDE_RATE_YAW       "L1A_data.h"
               CR_FR_A_ON              "L1A_data.h"
               CR_FR_B_ON              "L1A_data.h"
               CR_SA_A_SCAN_ON         "L1A_data.h"
               CR_SA_B_SCAN_ON         "L1A_data.h"
               ELEC_SIDES              "GEO_geo.h"
               ENCODER_LENGTH          "GEO_geo.h" 
               EARTH_ENCODER_TIMES     "L1A_data.h"
               VIEW_SECTOR_START       "L1A_data.h"
               LAST_VALID_SCAN         "L1A_data.h"
               CURR_SC_ANCIL_DATA      "L1A_data.h"
               MAJCYC3COF7             "L1A_data.h"
               MAJCYC5BOF7             "L1A_data.h"
               PRIOR_SC_ANCIL_DATA     "L1A_data.h"
               MAPIOK                  "mapi.h"
               MAX_SCAN_NUMBER         "GEO_geo.h"
               MODIS_E_BAD_INPUT_ARG   "PGS_MODIS_35251.h"
               MODIS_E_GEO             "PGS_MODIS_35251.h"
               MODIS_W_GEO_MIRR_CHAN   "PGS_MODIS_35251.h"
               SC_PRIOR                "GEO_input.h"
               SC_CURRENT              "GEO_input.h"
               SC_POSITION_X           "L1A_data.h"
               SC_POSITION_Y           "L1A_data.h"
               SC_POSITION_Z           "L1A_data.h"
               SC_VELOCITY_X           "L1A_data.h"
               SC_VALOCITY_Y           "L1A_data.h"
               SC_VALOCITY_Z           "L1A_data.h"
               SECTOR_LENGTH           "GEO_input.h"

Called By:
               GEO_prepare_l1a_data

Call functions:
               modsmf - writes error status messages to log
               getMODISarray - MAPI routine to read from an array.
               getMODIStable - MAPI routine to read from a table.

Requirements: 
               PR03-F-2.1-1
               PR03-F-2.1-2
               PR03-I-3

!Revision History:
 * $Log: GEO_read_L1Apacket_data.c,v $
 * Revision 6.2  2010/06/01 20:33:30  kuyper
 * Changed to return status codes.
 * Made error messages more informative.
 *
 * Revision 6.1  2010/03/17 13:44:55  gbritz
 * Add ephemeris/attitude quality information to MOD03 file.
 *
 * Georgios Britzolakis georgios.britzolakis-1@nasa.gov
 *
 * Revision 5.1  2005/04/21 20:45:45  kuyper
 * Correct error message.
 *
 * Revision 4.1  2003/02/21 22:48:33  kuyper
 * Corrected to use void* pointers with %p format code.
 *
 * Revision 3.3  2002/06/13 22:52:21  kuyper
 * Removed unnecessary NCSA acknowledgement.
 *
 * Revision 3.2  2001/04/13 22:50:27  kuyper
 * Changed back to uint16 for earth_encoder_times. Re-written to portably
 *   decode the incorrectly typed data in the raw_mir_enc SDS.
 *
 * Revision 3.1  2001/04/02  20:34:21  kuyper
 * Corrected message mnemonic for getMODISarray message.
 *
 * Revision 2.8  2001/04/02  14:49:30  seaton
 * Changed
 *
 * Revision 2.7  2001/04/02  14:44:50  seaton
 * Entered walkthrough corrections.
 *
 * Revision 2.6  2001/03/14  17:59:00  seaton
 * Added electronics side information from the formatter to the set of data read in.
 *
 * Revision 2.5  1999/02/05  19:03:01  seaton
 * Made walkthrough modifications.
 *
 * Revision 2.4  1999/02/02  15:01:47  seaton
 * Changes to reflect S/C Ancillary Data being read from Vdatas instead of SDS from L1A product.
 *
 * Revision 2.3  1998/03/04  00:37:50  jjb
 * Changed earth_encoder_times and view_sector_start parameters
 * 	to uint16's.
 *
 * Revision 2.2  1997/11/06  21:39:34  kuyper
 * Make implicit conversion explicit.
 *
 * Revision 2.1  1997/10/21  18:16:22  kuyper
 * Returned from ClearCase
 *
 * Revision /main/GEO_V2_DEV/3 1997/09/26 kuyper
 * Fixed error messages.
 *
 * Revision /main/GEO_V2_DEV/2  1997/08/25  Ding
 * Added in mirror channel fields.
 * ding@ltpmail.gsfc.nasa.gov
 *
 * Revision 1.4.1.1  1997/07/21  22:20:17  kuyper
 * Merged in out-of-sequence changes.
 *
 * Revision 1.4  1997/07/21  16:24:34  kuyper
 * Baselined Version 1
 *
 * Revision 1.3  1997/03/26  18:13:15  fhliang
 * Initial revision SDST delivery of GEO_read_L1Apacket_data.c.
 *
 *Parallel development:
 * Revision 1.3  1997/06/02  17:55:35  kuyper
 * Merged seed files.
 *
		Revision 1.2  1997/01/14  23:35:44  kuyper
		Adjusted header file list.
		James Kuyper Jr. <kuyper@ltpmail.gsfc.nasa.gov>
		
		Revision 1.1  1996/12/20 16:15:39  mikej
		Initial revision

		
!Team-unique Header:
               This software is developed by the MODIS Science Data Support
               Team for the National Aeronautics and Space Administration,
               Goddard Space Flight Center, under contract NAS5-32373.
               
!END****************************************************************************
*/

{
#define N_ELEMENTS(foo) ((int) (sizeof(foo) / sizeof(foo[0])))
#define STRMAX 512
#define ARGMAX 80
#define MAX_SC_DATA_NAME 30
#define LAST_VALID_SCAN_INDEX 2

    /* Warning: The order and sequence of the packet_data initilizations can
    * not be changed. We finish the initilization later and count on this
    * order.
    */
    struct pack{
	char *name; /* name of SDS object in the L1A file */
	long int dims[2]; /* extent of locations in SDS */
	void * data; /* pointer to the data */
    }		packet_data[] = {
	{EARTH_ENCODER_TIMES, {MAX_SCAN_NUMBER, ENCODER_LENGTH}, NULL},
	{VIEW_SECTOR_START, {MAX_SCAN_NUMBER, SECTOR_LENGTH}, NULL}
    };
    long	buffsize;
    int		obj;
    char	msgbuf[STRMAX];
    static const char filefunc[] = __FILE__ ", GEO_read_L1Apacket_data";

    static char	SAfieldnm[] =
	CR_SA_A_SCAN_ON "," CR_SA_B_SCAN_ON "," LAST_VALID_SCAN;
    static char	FRfieldnm[] = CR_FR_A_ON "," CR_FR_B_ON "," LAST_VALID_SCAN;
    uint16	FRdata[MAX_SCAN_NUMBER][3];
    uint16	SAdata[MAX_SCAN_NUMBER][3];
    int16	raw_mir_enc[MAX_SCAN_NUMBER][ENCODER_LENGTH];
    /* SC Ancillary Data Field Structure */
    static struct{
       char *field_names;
       size_t offset;
       long int size;
    }		vdata_fields[] =
    {
       {TIME_STAMP,  offsetof(sc_ancil_struct, second_header),
	    (long)sizeof(sc_ancillary_data[0].second_header) },
       {SC_POSITION_X "," SC_POSITION_Y "," SC_POSITION_Z "," 
	    SC_VELOCITY_X "," SC_VELOCITY_Y "," SC_VELOCITY_Z,
	    offsetof(sc_ancil_struct, posvel),
	    (long)sizeof(sc_ancillary_data[0].posvel) },
       {ATTITUDE_ANGLE_ROLL "," ATTITUDE_ANGLE_PITCH "," ATTITUDE_ANGLE_YAW ","
	    ATTITUDE_RATE_ROLL "," ATTITUDE_RATE_PITCH "," ATTITUDE_RATE_YAW,
	    offsetof(sc_ancil_struct, attit_angvel),
	    (long)sizeof(sc_ancillary_data[0].attit_angvel) }
    };
    char	*ancil_data_name[2];
    int		phase;
    int		scan;

    if (number_of_scans == 0)
     return PGS_S_SUCCESS;

    /*   Checking input values */
    if(l1a_file == NULL || number_of_scans < 0 ||
	number_of_scans > MAX_SCAN_NUMBER || sc_ancillary_data == NULL ||
      SAside == NULL || FRside == NULL || ss_cp_mode == NULL)
    {
	sprintf(msgbuf, "l1a_file = %p, number_of_scans = %d, "
	    "sc_ancillary_data = %p, SAside = %p, FRside = %p, ss_cp_mode = %p",
	    (void*)l1a_file, number_of_scans, (void*)sc_ancillary_data,
	    (void*)SAside, (void*)FRside, (void*)ss_cp_mode);
       modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

       return MODIS_E_BAD_INPUT_ARG;
    }

    ancil_data_name[SC_PRIOR] = PRIOR_SC_ANCIL_DATA;
    ancil_data_name[SC_CURRENT] = CURR_SC_ANCIL_DATA;

/*   Finish initializing the packet_data array with the contents
	defined above. 
	
	See warning about sequence in declaration section 
*/
    packet_data[0].dims[0] = (long int) number_of_scans;
    packet_data[0].data = raw_mir_enc;
    packet_data[1].dims[0] = (long int) number_of_scans;
    packet_data[1].data = view_sector_start;


    for (obj = 0; obj < N_ELEMENTS(packet_data); obj++)
    {
	long start[2] = {0L, 0L};

	if (getMODISarray(l1a_file, packet_data[obj].name, L1A_ENGINEERING_GRP, 
	    start, packet_data[obj].dims, packet_data[obj].data) != MAPIOK)
	{
	    sprintf(msgbuf, "getMODISarray(\"%s\", \"%s\", \""
		L1A_ENGINEERING_GRP "\")", l1a_file->filename,
		packet_data[obj].name);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);

	    return MODIS_E_GEO;
	}
    }

    /* raw_mir_enc is incorrectly stored in the L1A file as an int16; it's
    * actually uint16. The decision has been made to not correct the MOD01
    * filespec. Therefore, the following loop emulates the 16-bit 2's
    * complement conversion that occured before the defect was detected, even
    * if on the current platform int16 is 1's complement, sign/magnitude, or
    * larger than 16 bits (which it will be, on platforms with no 16 bit
    * type).
    */
    for(scan = 0; scan < MAX_SCAN_NUMBER; scan++)
    {
	int enc;

	for(enc = 0; enc < ENCODER_LENGTH; enc++)
	{
	    int32 i32 = (int32)raw_mir_enc[scan][enc];
	    earth_encoder_times[scan][enc] =
		(uint16)(i32 <0 ? (int32)0x10000 + i32 : i32);
	}
    }



/*   FOR Each Prior and Current S/C Ancillary Data  */
    for (phase = SC_PRIOR; phase <= SC_CURRENT; phase++)
    {
	for (obj = 0; obj < N_ELEMENTS(vdata_fields); obj++)
	{
	    buffsize = vdata_fields[obj].size;
	    if (getMODIStable(l1a_file, ancil_data_name[phase],
		L1A_ENGINEERING_GRP, vdata_fields[obj].field_names, 0L,
		(long)number_of_scans, &buffsize,
		(unsigned char *)&sc_ancillary_data[phase] +
		vdata_fields[obj].offset) != MAPIOK)
	    {
		sprintf(msgbuf, "getMODIStable(\"%s\", \"%s\", \""
		    L1A_ENGINEERING_GRP "\", \"%s\", %ld)\n",
		    l1a_file->filename, ancil_data_name[phase],
		    vdata_fields[obj].field_names, (long)number_of_scans);
		modsmf(MODIS_E_GEO, msgbuf, filefunc);
		
		return MODIS_E_GEO;
	    } /* if getMODIStable */
	} /* end for each element */
    } /* end for phase */


/* retrieve mirr_chan data */
 
    buffsize=(long)sizeof(FRdata);
    if(getMODIStable(l1a_file, MAJCYC3COF7, L1A_ENGINEERING_GRP, FRfieldnm, 0L,
	(long)number_of_scans, &buffsize, (unsigned char *) FRdata) != MAPIOK)
    {
	sprintf(msgbuf, "getMODIStable(\"%s\", \"" MAJCYC3COF7 "\", \""
	    L1A_ENGINEERING_GRP "\", \"%s\", 0, %ld)\n",
	    l1a_file->filename, FRfieldnm, (long)number_of_scans);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	return MODIS_E_GEO;
    } 


    buffsize=(long)sizeof(SAdata);
    if(getMODIStable(l1a_file, MAJCYC5BOF7, L1A_ENGINEERING_GRP, SAfieldnm, 0L,
	(long)number_of_scans, &buffsize, (unsigned char *) SAdata) != MAPIOK)
    {
	sprintf(msgbuf, "getMODIStable(\"%s\", \"" MAJCYC5BOF7 "\", \""
	    L1A_ENGINEERING_GRP "\", \"%s\", 0, %ld)\n",
	    l1a_file->filename, SAfieldnm, (long)number_of_scans);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	return MODIS_E_GEO;
    }

    buffsize = (long)(number_of_scans*sizeof*ss_cp_mode);
    if(getMODIStable(l1a_file, M01MAJCYCALL1, L1A_ENGINEERING_GRP,
	M01SS_CP_MODE, 0L, (long)number_of_scans, &buffsize,
	(unsigned char *)ss_cp_mode) != MAPIOK)
    {
	sprintf(msgbuf, "getMODIStable(\"%s\", \"" M01MAJCYCALL1 "\", \""
	    L1A_ENGINEERING_GRP "\", \"" M01SS_CP_MODE "\", 0, %ld)\n",
	    l1a_file->filename, (long)number_of_scans);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
	
	return MODIS_E_GEO;
    }


    if((SAdata[number_of_scans-1][LAST_VALID_SCAN_INDEX])==65535)
    {
	modsmf(MODIS_W_GEO_MIRR_CHAN, "", filefunc);
	for(obj=number_of_scans-1; obj>=0; obj--)
	{
	    FRside[obj][CHAN_A]=0;
	    FRside[obj][CHAN_B]=0;
	    SAside[obj][CHAN_A]=0;
	    SAside[obj][CHAN_B]=0;
	}
    }
    else for(obj=number_of_scans-1; obj>=0; obj--)
    {
	if(FRdata[obj][LAST_VALID_SCAN_INDEX] != 65535)
	{
	    FRside[obj][CHAN_A] = FRdata[obj][CHAN_A];
	    FRside[obj][CHAN_B] = FRdata[obj][CHAN_B];
	}
	else
	{
	    if (obj == number_of_scans-1)
	    {
		FRside[obj][CHAN_A] = 0;
		FRside[obj][CHAN_B] = 0;
	    }
	    else
	    {
		FRside[obj][CHAN_A] = FRside[obj + 1][CHAN_A];
		FRside[obj][CHAN_B] = FRside[obj + 1][CHAN_B];
	    }
	}

	if(SAdata[obj][LAST_VALID_SCAN_INDEX] != 65535)
	{
	    SAside[obj][CHAN_A] = SAdata[obj][CHAN_A];
	    SAside[obj][CHAN_B] = SAdata[obj][CHAN_B];
	}
	else
	{
	    if(obj == number_of_scans - 1)
	    {
		SAside[obj][CHAN_A] = 0;
		SAside[obj][CHAN_B] = 0;
	    }
	    else
	    {
		SAside[obj][CHAN_A] = SAside[obj + 1][CHAN_A];
		SAside[obj][CHAN_B] = SAside[obj + 1][CHAN_B];
	    }
	} 
    }
    return PGS_S_SUCCESS;
}
