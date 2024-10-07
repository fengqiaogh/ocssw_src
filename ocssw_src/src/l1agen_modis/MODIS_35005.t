%INSTR = MODIS
%LABEL = MODIS
%SEED = 35005

MODIS_E_ARRAY_OUTPUT_ERR        Unable to output data to SDS
MODIS_E_ATTACHED_VDATAS         Attempting to close the L1A output file with 
                                Vdatas still attached
MODIS_E_CAL_FR_CNT_INVALID      MODIS calibration frame count outside valid range
MODIS_E_CANT_PKT_PEEK           Cannot peek at next packet in the L0 file
MODIS_E_CHECKSUM_NOT_VALID      Calculated checksum
MODIS_E_CLOSE_VDATA             Unable to detach a Vdata from the L1A file
MODIS_E_CORRUPT_PKT_HDR         Could not successfully unpack packet
                                header
MODIS_E_CREATE_L1A_GRANULE      Error creating L1A granule
MODIS_E_CREATE_VDATA            Unable to create a Vdata
MODIS_E_CREATE_VDATA_FIELD      Unable to create a Vdata field
MODIS_E_DATATYPE_FAIL           Datatype not compatible with Vdata name
MODIS_E_EARTH_FR_CNT_INVALID    MODIS Earth frame count outside valid range.
MODIS_E_END_VDATA_FILE          Could not clean up Vdata
MODIS_E_ENQUEUE                 Could not successfully put the failed
                                packet on the failed_pkt_queue
MODIS_E_FAILED_TIMECODE_CONV    Failed to convert time code in MODIS
                                packet
MODIS_E_FIELD_NAME_LIST_OVERRUN The Vdata field name list for VSsetfields
                                exceeded the size defined in
                                VU_MAX_FIELD_NAME_LIST_SIZE 
MODIS_E_FREE_QUEUE              Attempting to free non-empty queue
MODIS_E_GATTRIB_FAILED          Unable to write global attribute
MODIS_E_GET_VALID_L0_FILE       Error in getting valid L0 file
MODIS_E_GETCONFIG_FAILED        A non-success message was returned while 
                                trying to retrieve data from the PCF file
MODIS_E_HANDLE_MISSING_SCANS    Error handling missing scans
MODIS_E_ID_TABLE_OVERFLOW       Attempting to insert into a full
                                Vdata_id_table
MODIS_E_INITIALIZATION_FAILED   Program failed initialization
MODIS_E_INITIALIZE_GLOBAL_MD    Not all global metadata values could be 
                                initialized properly
MODIS_E_INV_APID                Packet APID does not fall within MODIS AM1
                                range
MODIS_E_INV_APID_TEST           Packet APID indicates this is a Test
                                packet
MODIS_E_INV_PKT_SEQ_FLAG        Group packet sequence flag is not valid
MODIS_E_INV_PKT_TIME            A packet has been detected with an invlaid
                                time code for the current scan
MODIS_E_INV_PKT_TYPE            MODIS packet type is not normal
MODIS_E_INV_QL_FLAG             Expedited data flag is not valid
MODIS_E_INV_SEC_HDR_FLAG        Packet secondary header flag is not valid
MODIS_E_INV_SECTOR		Packet source id code inconsistent with packet type
MODIS_E_INV_TYPE                Packet primary header type is not valid
MODIS_E_INV_VERSION             Packet primary herader version is not
                                valid
MODIS_E_INVALID_VDATA_ORDER     The order stored in the eng_data structure
                                is invalid
MODIS_E_INVALID_VDATA_TYPE      The datatype stored in the eng_data
                                structure is invalid
MODIS_E_LO_GETPACKET_FAILED     PGS_IO_L0_GetPacket returned a non-success
                                message
MODIS_E_L1A			Error indication returned by function.
MODIS_E_MALLOC_FAILED           Unable to dynamically allocate memory
MODIS_E_MAX_VDATAS_EXCEEDED     Attempting to attach more than
                                VU_MAX_NUMBER_OF_VDATAS Vdatas to the L1A
                                file
MODIS_E_MET_SETATTR_FAILED      Unable to set metadata value
MODIS_E_NO_SCANS_IN_PRODUCT     Time of current scan is greater than
                                expected end of current L1A product 
MODIS_E_NULL_POINTER            A NULL Pointer was passed into procedure
MODIS_E_PGS_IO_GEN_OPEN         PGS_IO_Gen_Open failed to open the
                                engineering data list file.
MODIS_E_PGS_IO_GEN_CLOSE        PGS_IO_Gen_Close failed to close the
                                engineering data list file
MODIS_E_PGS_IO_LO_CLOSE         PGS_IO_L0_Close failed to close the L0 file
                                successfully
MODIS_E_RECALL_ID               Error in recalling Vdata Id
MODIS_E_SDS_CREATE_FAILED       Unable to initialize SDS
MODIS_E_TAI_TO_UTC_FAILED       PGE_TD_TAItoUTC returned a non-success
                                message Unable to convert TAI time to UTC
                                representation
MODIS_E_TOO_MANY_BITS           Function can extract up to 32 bits.
MODIS_E_UTC_TO_TAI_FAILED       PGS_TD_UTCtoTAI returned a non-success
                                message Unable to convert UTC time to TAI
                                representation
MODIS_E_SCAN_PROCESS            Unable to process scan successfully
MODIS_E_SCANCNT_NOT_VALID       Packet Number:
MODIS_E_VALIDATE_L0_HEADER      Corruption detected in the L0 File Header
MODIS_E_VDATA_BUFFER_OVERFLOW   The S/C Ancillary Data has overflowed the
                                data buffer
MODIS_E_VDATA_WRITE_FAILED      Could not write to Vdata
MODIS_E_VSFDEFINE               VSfdefine failed to define the Vdata field
                                types
MODIS_E_VSSETFIELDS             VSsetfields failed to set the Vdata field
                                names
MODIS_E_WRITE_GLOBAL_METADATA   Error writing global metadata
MODIS_E_WRITE_SCAN_FAIL         Error writing a scan
MODIS_E_WRITE_VDATA             Failed to successfully write the Vdata to
                                L1A file
MODIS_F_GETREF_FAILED           PGS_PC_GetReference returned a non-success
                                message
MODIS_F_INVAL_L0_FILE_RET       Cannot retrieve a valid L0 file
MODIS_F_INVALID_ENG_DATA_LIST   The engineering data list file is
                                corrupted
MODIS_F_L0_HEADER_VAL_FAILED    L0 file header failed validation
MODIS_F_L0_SETSTART_FAILED      Unable to find an available packet with
                                time after start time of preloading
                                engineering packets
MODIS_F_NO_PCF_CONFIG_DATA      All configuration data could not be
                                retrieved from the PCF file
MODIS_F_PKT_READ_FAILED         Fatal error reading a packet
MODIS_F_PROCESSAGRANULE         Error processing a granule
MODIS_F_UNABLE_TO_INIT_MCF      Fatal error initializing a MCF file
MODIS_F_WRITE_ENG_DATA_FAIL     Fatal error writing engineering data 
MODIS_M_PKT_NOT_IN_SCAN         Packet not in current scan
MODIS_S_SUCCESS                 Routine ran successfully
MODIS_U_L1A_BEGIN               The L1A process started   
MODIS_U_L1A_END                 The L1A process ended
MODIS_W_CREATE_L1A_GRANULE      Warning: No packets available to create L1A granule
MODIS_W_EXCESS_DATA_LIST        Warning: There are extra lines in the eng
                                data list file
MODIS_W_ID_TABLE_NOT_INIT       Warning: Vdata ID table has not been initialized
MODIS_W_NEGATIVE_VDATA_CTR      Warning: Negative Vdata counter  
MODIS_W_NO_MORE_PACKETS         Warning: No more packets in Virtual file
MODIS_W_INV_PKTLEN              Warning: Packet length is not equal to a MODIS
                                short or long packet length
MODIS_W_PRELOAD_ENG_DATA        Warning: There was no engineering data to
                                preload for the specified time given
MODIS_W_VDATA_NOT_IN_TABLE      Warning: Vdata not found in table
MODIS_W_TOO_MANY_SCANS          Warning: Too many scans of data in granule

