Algorithm - MODIS Level 1A (MOD_PR01)
Date      - 2015-08-21


1. README for MOD_PR01 v6.0.6 (Level 1A)


2. POINTS OF CONTACT:
  STM:

        James Kuyper
        NASA/GSFC MODIS Project - Science Data Support Team (SDST)
        7425 Forbes Blvd, Room 563
        Seabrook, MD 20706
        Telephone: (301) 867-2064
        E-mail: James.R.Kuyper@NASA.gov

  SSTG Contact:
	
	James Kuyper (see above)


3. BUILD
  To build MOD_PR01.exe:

    (1) Setup toolkit env vars (on moddev_c64)

          source ~/svn/COMMON/Scripts/Linux/64b/setup-5.2.17.csh.moddevC64 -f90
	
    (2) Use smfcompile to make the C error message include file:

	  $PGSBIN/smfcompile -f MODIS_35005.t -r

    (3) Build the executable by running the makefile:

	  make -f MOD_PR01.mk 

  To run MOD_PR01.exe:
    modify MOD_PR01.pcf to reflect the correct directories for input and 
    output data files and runtime parameters.  Type:

	MOD_PR01.exe


4.  EVIRONMENT VARIABLES used by the make file:
    See list in prologue of MOD_PR01.mk

5.  M-API version 6.0.1g was used in this software.


6.  HDF version 4.2r5 was used in this software.


7.  SDP TOOLKIT version 5.2.17 was used in this software.


8.  This code was developed on an Intel Xeon processor using Centos Linux.


9.  ANCILLARY DATA sets used in processing:
    None


10. MODIS PRODUCTS used as inputs:
    This software requires MOD000 or MOD00F (Level 0 partitioned data sets) as
    inputs.
    If a run requires only one L0 virtual data set, then that set is
    considered to be the "Current" L0 file and should be given the "Current"
    L0 file LUN (599002) in the PCF file.


11. OTHER 
    This software also requires an Engineering Data List Lookup
    table, ENG_DATA_LIST_TERRA or ENG_DATA_LIST_AQUA, as an input.

12. MODIS OUTPUT PRODUCTS
    MOD01 - MODIS/Terra Raw Radiances in Counts 5-Min L1A Swath
    MYD01 - MODIS/Aqua Raw Radiances in Counts 5-Min L1A Swath

13. PROBLEMS:
    None

14. DELIVERY SIZE:  1.25 MB


15. OUTPUT FILE SIZE:
    The approximate size of a Level 1A Day granule is 558.2 MB.
    The approximate size of a Level 1A Night granule is 186.5 MB.
    The approximate size of a Level 1A Empty granule is 0.5 MB.


16. TESTS PERFORMED:
    A 5 minute granule test was run on the L1A program for each of the
    following: day scan granule, night scan granule, and a mixed scan granule.


17. EXIT CODES:
    0	Indicates successful completion
    1	Indicates unsuccessful completion


18. ERROR LIST and operator actions:



MODIS_E_ARRAY_OUTPUT_ERR
        Error:          Unable to output data to SDS
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        handle_missing_scans.c write_pix_qual.c
                        write_scan.c write_scan_data.c write_scan_metadata.c

MODIS_E_ATTACHED_VDATAS
        Error:          Attempting to close the L1A output file with
                        Vdatas still attached
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        create_Vdata.c end_Vdata_access_to_file.c

MODIS_E_CAL_FR_CNT_EXC_LIM
        Error:          MODIS calibration frame count exceeds established limit
        Solution:       Possible corruption of L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        unpack_MODIS_header.c process_a_packet.c
                        unpack_packet_header.c

MODIS_E_CANT_PKT_PEEK
        Error:          Cannot peek at next packet in the L0 file
        Solution:       Internal error.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        process_a_packet.c 

MODIS_E_CHECKSUM_NOT_VALID
        Error:          Calculated pkt checksum did not match pkt checksum
        Solution:       Possible corruption of L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        check_checksum.c load_eng_data.c
                        process_a_packet.c process_a_scan.c
                        update_scan_metadata.c

MODIS_E_CLOSE_VDATA
        Error:          Unable to detach a Vdata from the L1A file
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.   
        Modules:        end_eng_data_access_to_file.c
        
MODIS_E_CORRUPT_PKT_HDR
        Error:          Could not successfully unpack packet header
        Solution:       Possible corruption of L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        load_eng_data.c
                        
MODIS_E_CREATE_L1A_GRANULE
        Error:          Error creating L1A granule
        Solution:       Notify SDST of fault. 
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        create_L1A_granule.c process_a_granule.c  
        
MODIS_E_CREATE_VDATA
        Error:          Unable to create a Vdata
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        create_Vdata.c init_L1A_HDF_vdatas.c
        
MODIS_E_CREATE_VDATA_FIELD
        Error:          Unable to create a Vdata field
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        create_Vdata.c
        
MODIS_E_DATATYPE_FAIL
        Error:          Datatype not compatible with Vdata name
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        create_Vdata_field.c
        
MODIS_E_EARTH_FR_CNT_EXC_LIM
        Error:          MODIS Earth frame count exceeds established limit
        Solution:       Possible corruption of L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        unpack_MODIS_header.c unpack_packet_header.c
        
MODIS_E_END_VDATA_FILE
        Error:          Could not cleanup Vdata
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        end_eng_data_access_to_file.c
        
MODIS_E_ENQUEUE
        Error:          Could not successfully put the failed packet on the
                        failed_pkt_queue
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.   
        Modules:        accumulate_failed_packets.c
                        
MODIS_E_FAILED_TIMECODE_CONV
        Error:          Failed to convert time code in MODIS packet   
        Solution:       Possible corruption of L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        unpack_packet_header.c unpack_secondary_header.c
                        
MODIS_E_FIELD_NAME_LIST_OVERRUN
        Error:          The Vdata field name list for VSsetfields exceeded the
                        size defined in VU_MAX_FIELD_NAME_LIST_SIZE
        Solution:       Increase the size of VU_MAX_FIELD_NAME_LIST_SIZE in the
                        header file VU_vdata_utilty.h (Code must be
                        re-built).
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.   
        Modules:        create_Vdata.c
                        
MODIS_E_FREE_QUEUE
        Error:          Attempting to free non-empty queue
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.   
        Modules:        free_queue.c
        
MODIS_E_GATTRIB_FAILED  
        Error:          Unable to write global attribute
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        write_global_+metadata.c
                        write_specific_granule_metadta.c
                        
MODIS_E_GET_VALID_L0_FILE
        Error:          Error in getting valid L0 file
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        get_valid_L0_file.c initialize_level1a.c
                        read_a_packet.c set_start_position.c
                        
MODIS_E_GETCONFIG_FAILED
        Error:          A non-success message was returned while trying to
                        retrieve data from the PCF file
        Solution:       Make sure all components in the .pcf file are
                        correct.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        get_pcf_config_data.c get_valid_L0_file.c
                        level1a.c parse_eng_data_list.c 
                        
MODIS_E_HANDLE_MISSING_SCANS
        Error:          Error handling missing scans  
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        handle_missing_scans.c process_a_granule.c
                        
MODIS_E_ID_TABLE_OVERFLOW
        Error:          Attempting to insert into a full Vdata_id_table
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        remember.c
                        
MODIS_E_INITIALIZATION_FAILED
        Error:          Program failed initialization
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        level1a.c
        
MODIS_E_INITIALIZE_GLOBAL_MD
        Error:          Not all global metadata values could be initialized
                        properly
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        process_a_granule.c
        
MODIS_E_INV_APID
        Error:          Packet APID does not fall within MODIS AM1 range
        Solution:       Possible corruption in input L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        process_a_packet.c unpack_packet_header.c
                        unpack_primary_header.c
                        
MODIS_E_INV_APID_TEST
        Error:          Packet APID indicates this is a Test packet
        Solution:       Possible curruption in input L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        process_a_packet.c unpack_packet_header.c
                        unpack_primary_header.c
        
MODIS_E_INV_PKT_SEQ_FLAG
        Error:          Group packet sequence flag is not valid
        Solution:       Possible curruption in input L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.   
        Modules:        process_a_packet.c unpack_packet_header.c
                        unpack_primary_header.c
                        
MODIS_E_INV_PKT_TIME            
       Error:           A packet has been detected with an invlaid
                        time code for the current scan
       Solutions:       Probable bit-flipped packet in L0 file.
                        No other actions required.
       Modules:         packet_of_scan.c, process_a_scan.c

MODIS_E_INV_PKT_TYPE    
        Error:          MODIS packet type is not normal
        Solution:       Possible curruption in input L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.   
        Modules:        process_a_packet.c unpack_packet_header.c
                        unpack_secondary_header.c
                        
MODIS_E_INV_QL_FLAG     
        Error:          Expedited data flag is not valid
        Solution:       Possible curruption in input L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.   
        Modules:        process_a_packet.c unpack_packet_header.c
                        unpack_secondary_header.c
                        
MODIS_E_INV_SEC_HDR_FLAG
        Error:          Packet secondary header flag is not valid
        Solution:       Possible curruption in input L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.   
        Modules:        process_a_packet.c unpack_packet_header.c
                        unpack_primary_header.c
                        
MODIS_E_INV_TYPE        
        Error:          Packet primary header type is not valid
        Solution:       Possible curruption in input L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.   
        Modules:        process_a_packet.c unpack_packet_header.c
                        unpack_primary_header.c
                        
MODIS_E_INV_VERSION     
        Error:          Packet primary herader version is not valid
        Solution:       Possible curruption in input L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.   
        Modules:        process_a_packet.c unpack_packet_header.c
                        unpack_primary_header.c
                        
MODIS_E_INVALID_VDATA_ORDER
        Error:          The order stored in the eng_data structure is invalid
        Solution:       Possible invalid value in the ENG_DATA_LIST.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.   
        Modules:        process_group2_packet1_vdata.c
        
MODIS_E_INVALID_VDATA_TYPE
        Error:          The datatype stored in the eng_data structure is
                        invalid
        Solution:       Possible invalid value in the ENG_DATA_LIST.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.   
        Modules:        init_L1A_HDF_vdatas.c
                        process_group2_packet1_vdata.c write_eng_data.c
                        
MODIS_E_LO_GETPACKET_FAILED
        Error:          PGS_IO_L0_GetPacket returned a non-success message
        Solution:       Possible curruption in input L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.   
        Modules:        NOT FOUND
        
MODIS_E_MALLOC_FAILED
        Error:          Unable to dynamically allocate memory
        Solution:       Possible out of memory on host machine.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        accumulate_failed_packets.cenqueue.c level1a.c   
                        make_queue.c
        
MODIS_E_MAX_VDATAS_EXCEEDED
        Error:          Attempting to attach more than VU_MAX_NUMBER_OF_VDATAS
                        Vdatas to the L1A file
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        remember.c
                        
MODIS_E_MET_SETATTR_FAILED
        Error:          Unable to set metadata value
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        write_ECS_metadata.c write_global_metadata.c
                        
MODIS_E_NO_SCANS_IN_PRODUCT
        Error:          Time of current scan is greater than expected end of
                        current L1A product
        Solution:       Possible curruption in input L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        level1a.c
                        
MODIS_E_NULL_POINTER
        Error;          A NULL Pointer was passed into procedure
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        create_Vdata.c create_eng_data_vdata_array.c
                        create_eng_data_vdata_array_field.c
                        init_L1A_HDF_vdatas.c parse_eng_data_list.c   
                        process_group2_packet1_vdata.c write_eng_data.c
                        
MODIS_E_PGS_IO_GEN_OPEN
        Error:          PGS_IO_Gen_Open failed to open the engineering data 
                        list file.
        Solution:       Check PCF file for LUN number of Eng data list
                        file aginst LUN number specified by the
                        PC_pcf_info.h file.
                        Check to make sure the ENG_DATA_LIST file exists
                        in the location specified in the .pcf file.
        Modules:        parse_eng_data_list.c

MODIS_E_PGS_IO_GEN_CLOSE
        Error:          PGS_IO_Gen_Close failed to close the engineering data
                        list file
        Solution:       Possible curruption in input L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        parse_eng_data_list.c
                        
MODIS_E_PGS_IO_LO_CLOSE
        Error:          PGS_IO_L0_Close failed to close the L0 file successfully
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        close_processing_run.c read_a_packet.c
                        set_start_position.c
                        
MODIS_E_RECALL_ID
        Error:          Error in recalling Vdata Id
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        close_Vdata.c write_Vdata.c
                        
MODIS_E_SDS_CREATE_FAILED
        Error:          Unable to initialize SDS
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        init_L1A_HDF_sdss.c init_L1A_pix_qual_HDF_sdss.c
                        init_L1A_scan_data_HDF_sdss.c
                        init_L1A_scan_meta_HDF_sdss.c
                        
MODIS_E_TAI_TO_UTC_FAILED
        Error:          PGE_TD_TAItoUTC returned a nun-success message Unable
                        to convert TAI time to UTC representation
        Solution:       Possible corruption in L0 file timecodes.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.   
        Modules:        initialize_global_metadata.c
        
MODIS_E_TOO_MANY_BITS   
        Error:          Function can extract up to 32 bits.
        Solution:       Extraction routine cannot extract over 32 bits.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        extr_bits.c
                        
MODIS_E_UTC_TO_TAI_FAILED
        Error:          PGS_TD_UTCtoTAI returned a non-success message Unable
                        to convert UTC time to TAI representation
        Solution:       Possible corruption in L0 file timecodes.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        get_pcf_config_data.c initialize_global_metadata.c
                        
MODIS_E_SCAN_PROCESS
        Error:          Unable to process scan successfully
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation. 
        Modules:        process_a_granule.c  
                        
MODIS_E_SCANCNT_NOT_VALID
        Error:          A Packet has been detected with an invalid
                        scan count field for the current scan
        Solutions:      Probable bit-flipped packet in L0 file.
                        No other actions required.
       Modules:         packet_of_scan.c, process_a_scan.c

MODIS_E_VDATA_BUFFER_OVERFLOW
        Error:          The S/C Ancillary Data has overflowed the data buffer
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        write_eng_data.c
        

MODIS_E_VSFDEFINE
        Error:          VSfdefine failed to define the Vdata field types 
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        create_Vdata_field.c
        
MODIS_E_VSSETFIELDS
        Error:          VSsetfields failed to set the Vdata field names  
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        create_Vdata.c
        
MODIS_E_WRITE_GLOBAL_METADATA
        Error:          Error writing global metadata
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        process_a_granule.c write_global_metadata.c
        
MODIS_E_WRITE_SCAN_FAIL
        Error:          Error writing a scan
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        process_a_granule.c write_scan.c
        
MODIS_E_WRITE_VDATA
        Error:          Failed to successfully write the Vdata to L1A file
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        write_Vdata.c write_eng_data.c
                        write_failed_packets.c
        
MODIS_F_GETREF_FAILED
        Error:          PGS_PC_GetReference returned a non-success message
        Solution:       Check PCF file for LUN number of variable in question
                        against LUN number specified by the PC_pcf_info.h file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        initialize_level1a.c
                        
MODIS_F_INVAL_L0_FILE_RET
        Error:          Cannot retrieve a valid L0 file
        Solution:       Check PCF file and/or LUN # specified in PC_pcf_data.h.
                        Possible file curruption.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        initialize_level1a.c
                        
MODIS_F_INVALID_ENG_DATA_LIST
        Error:          The engineering data list file is corrupted
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        initialize_level1a.c parse_eng_data_list.c
                        
MODIS_F_L0_HEADER_VAL_FAILED
        Error:          L0 file header failed validation
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        get_valid_L0_file.c set_start_position.c
                        validate_L0_header.c
                        
MODIS_F_L0_SETSTART_FAILED
        Error:          Unable to find an available packet with time after
                        start time of preloading engineering packets
        Solution:       Possible incorrect or invalid L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        initialize_level1a.c set_start_position.c
                        
MODIS_F_NO_PCF_CONFIG_DATA
        Error:          All configuration data could not be retrieved from the
                        PCF file
        Solution:       Check PCF file and/or LUN # specified in PC_pcf_data.h.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        initialize_level1a.c
                        
MODIS_F_PKT_READ_FAILED 
        Error:          Fatal error reading a packet
        Solution:       Possible curruption in input L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        initialize_level1a.c load_eng_data.c
                        process_a_granule.c process_a_scan.c read_a_packet.c
        
MODIS_F_PROCESSAGRANULE 
        Error:          Error processing a granule
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        level1a.c process_a_granule.c
                        
MODIS_F_UNABLE_TO_INIT_MCF
        Error:          Fatal error initializing a MCF file
        Solution:       Make sure the .mcf file is in the location
                        specified by the .pcf file.
                        If it is, Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        process_a_granule.c
                        
MODIS_F_WRITE_ENG_DATA_FAIL
        Error:          Fatal error writing engineering data
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        handle_missing_scans.c process_a_granule.c
                        write_eng_data.c write_scan.c

MODIS_M_PKT_NOT_IN_SCAN
        Error:          Packet not in current scan
        Solution:       Normal ending return code of the end of the current
                        scan. No action required.
        Modules:        packet_of_scan.c, process_a_scan.c
                        
MODIS_S_SUCCESS
        Error:          Routine ran successfully
        Solution:       No action required.
        Modules:        accumulate_failed_packets.c check_checksum.c
                        close_Vdata.c create_L1A_granule.c create_Vdata.c
                        create_Vdata_field.c end_vdata_access_to_file.c
                        end_eng_data_access_to_file.c enqueue.c extr_bits.c
                        failed_packets.c get_pcf_config.data.c
                        get_valid_L0_file.c handle_missing_scans.c
                        init_L1A_HDF_sdss.c init_L1A_HDF_vdatas.c
                        init_L1A_pix_qual_HDF_sdss.c
                        init_L1A_scan_qual_HDF_sdss.c
                        initialize_global_metadata.c initialize_level1a.c
                        load_eng_data.c parse_eng_data_list.c
                        process_a_granule.c process_a_packet.c
                        process_a_scan.c read_a_packet.c rmpath.c
                        set_start_position.c unpack_MODIS_header.c
                        unpack_packet_header.c unpack_primary_header.c
                        unpack_secondary_header.c update_scan_metadata.c
                        validate_L0_header.c write_ECS_metadata.c
                        write_Vdata.c write_eng_data.c
                        write_failed_packets.c write_global_metadata.c
                        write_pix_qual.c write_scan.c write_scan_data.c
                        write_scan_metadata.c
                        write_specific_granule_metadata.c
                        
MODIS_U_L1A_BEGIN
        Error:          The L1A process started
        Solution:       A LogStatus message stating the start time of L1A
                        procedure and the build date. Not an error,
                        for information only.
        Modules:        level1a.c
                        
MODIS_U_L1A_END
        Error:          The L1A process ended
        Solution:       A LogStatus message stating the stop time of the
                        L1A processing and the return code from the run.
        Modules:        level1a.c
                        
MODIS_W_EXCESS_DATA_LIST
        Error:          Warning: There are extra lines in the eng data list file
        Solution:       Extra lines were encountered at the end of the 
                        ENG_DATA_LIST. This is not an necessairly an error.
                        Check the ENG_DATA_LIST to make sure that the nature
                        of the lines at the end of the file are
                        understood.
        Modules:        parse_eng_data_list.c  
        
MODIS_W_ID_TABLE_NOT_INIT
        Error:          Warning: Vdata ID table has not been initialized
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        forget.c recall_id.c 
        
MODIS_W_INV_PKTLEN
        Error:          Packet length is not equal to a MODIS short or long
                        packet length
        Solution:       Possible curruption in input L0 file.
                        Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.   
        Modules:        process_a_packet.c unpack_packet_header.c
                        unpack_packet_header.c
                        
MODIS_W_NEGATIVE_VDATA_CTR
        Error:          Warning: Negative Vdata counter
        Solution:       Notify SDST of fault.
                        Retain LogStatus file and products generated when
                        error message occurred for SDST investigation.
        Modules:        forget.c
                        
MODIS_W_NO_MORE_PACKETS
        Error:          Warning: No more packets in Virtual file
        Solution:       This is stating that the L0 file ran out of
                        packets before the L1A process was complete. The L1A
                        process should try to open and execute the second L0
                        file (LUN 599002) if this happens during
                        processing on LUN 599001.
        Modules:        initialize_level1a.c load_eng_data.c
                        process_a_granule.c process_a_scan.c read_a_packet.c
                        
MODIS_W_PRELOAD_ENG_DATA
        Error:          Warning: There was no engineering data to preload for
                        the specified time given
        Solution:       This is just a warning that engineering data could
                        not be preloaded. This could be caused by not having a
                        70 second buffer of L0 data before the granule   
                        start time.
        Modules:        initialize_level1a.c load_eng_data.c
                        
MODIS_W_VDATA_NOT_IN_TABLE
        Error:          Warning: Vdata not found in table
        Solution:       The program is trying to access a Vdata name that
                        does not have a Vdata Id associated with it. Data for
                        this Vdata name will be lost.
        Modules:        forget.c
                        
        
18. PLATFORM NOTES:
a) The pcf file and the ENG_DATA_LIST file must correspond. 

Instrument	PCF file			ENG_DATA_LIST file
----------	--------			------------------
Terra		SatelliteInstrument = AM1M	ENG_DATA_LIST_TERRA
Aqua		SatelliteInstrument = PM1M	ENG_DATA_LIST_AQUA

If these values do not match, MOD_PR01 will exit 1 with a corresponding
LogStatus message.

NOTE: for MOD_PR01 there cannot be combined run. The above value must be 
either AM1M or PM1M. The AP1M combined value is not supported in MOD_PR01.

b) The output ESDT names will differ in the second character of the name as
documented in the MODIS PGE Modifications for Aqua paper. (MO for Terra
MY for Aqua).

c) The mcf files will be different for the two spacecraft. The mcf file for
Terra products is MOD01.mcf and the file for Aqua products is MYD01.mcf. The
mcf fields that are different are as follows:

				Terra		Aqua
                                -----		----
ShortName			MOD01		MYD01
AssociatedPlatformShortName	Terra		Aqua

d) The LocalGranuleID names will differ as the first 4 characters of the 
LocalGranuleID are the ESDT name. See #2 from above.

