#ifndef L1A_PROTOTYPE_H
#define L1A_PROTOTYPE_H

/*
!C-INC************************************************************************

!Description:  This include file contains the prototypes for all L1A routines.

!Input Parameters: N/A

!Output Parameters: N/A

Externally Defined:  
               PGSt_SMF_status                 (PGS_SMF.h)
               PGSt_IO_L0_Packet               (PGS_IO.h)
               PGSt_IO_L0_VirtualDataSet       (PGS_IO.h)
               PH_PACKET_HEADER_t              (PH_pkt_hdr.h)
               FP_QUEUE_t                      (FP_failed_pkt_queue.h)
               PGSt_scTime                     (PGS_TD.h)
               PGSt_double                     (PGS_TD.h)

!Revision History:
 $Log: L1A_prototype.h,v $
 Revision 5.1  2007/01/24 22:36:31  kuyper
 Added gran_start_time parameter to handle_missing_scans().

 Revision 4.9  2003/11/12 21:04:52  kuyper
 Changed to define and use MAX_INPUTS macro.

 Revision 4.8  2003/09/30 18:49:29  kuyper
 Reinstated process_a_packet() prototype.

 Revision 4.7  2003/05/09 16:41:47  kuyper
 Changed process_next_packet to accept packet_header as a pointer to a pointer.

 Revision 4.6  2003/03/05 22:10:56  kuyper
 Added global_time_offset_array.

 Revision 4.5  2003/02/24 21:43:15  kuyper
 Changed pointer parameters of log_fmt_msg() to be pointers to
   const.

 Revision 4.4  2002/10/21 19:10:29  vlin
  global variables "global_first_gran_start_time" & "global_last_gran_stop_time" added

 Revision 4.2  2002/09/19 14:31:58  vlin
 prototype for function "initialize_level1a" updated.

 Revision 4.1  2002/08/28 18:44:43  vlin
 Routines initialize_global_metadata, get_pcf_config_data, process_a_granule,
 initialize_scan_data, process_next_packet, and process_a_packet are updated.
 global_input_pointer added.
 vlin@saicmodis.com

               Revision 3.0 1998/10/16   09:40:00 EDT
               John Seaton/SAIC/GSC (seaton@ltpmail.gsfc.nasa.gov)
               Added new procedure process_group2_packet1_vdata.

               Revision 2.0  1997/07/10  11:00
               Tom Johnson/SAIC/GSC (johnson@ltpmail.gsfc.nasa.gov)
               Added include files for packet processing

               Revision 1.0  1997/09/22  10:00
               Qi Huang (qhuang@ltpmail.gsfc.nasa.gov)
               Original include file

!Team-unique Header:

     This software is developed by the MODIS Science Data Support Team 
     for the National Aeronautics and Space Administration, 
     Goddard Space Flight Center, under contract NAS5-32373.

 Design Notes: 
     The ".h" below was specifically written for development in C.
     Any other language choice may require reworking of the ".h"
     before coding can begin.

!END**********************************************************************
*/

#include "hdf.h"
#include "mapi.h" 
#include "PGS_IO.h"
#include "PGS_SMF.h"
#include "EN_eng_data.h"
#include "FP_failed_pkt_queue.h"
#include "MD_metadata.h"
#include "PD_pkt_data.h"
#include "PH_pkt_hdr.h"
#include "SC_scan.h"
#include "VU_vdata_utility.h"
#include "PC_pcf_info.h"


/************************* L1A Routine Prototype *****************************/

PGSt_SMF_status   accumulate_failed_packets (PGSt_IO_L0_Packet  *pkt, 
                                             FP_QUEUE_t          failed_pkts);

int16  attached_Vdata_counter (int16 action);

PGSt_SMF_status   check_checksum (PH_PACKET_HEADER_t  pkt_header,
                                  uint16              *pkt_contents);

void   close_processing_run (PGSt_IO_L0_VirtualDataSet L0_file);

PGSt_SMF_status   close_Vdata (char *Vdata_name);

void   compute_SD_start_time (PH_PACKET_HEADER_t    *pkt_header,
                             PGSt_double           *SD_start_time);

void   compute_global_time_offsets  (PGSt_double  scan_rate);

PGSt_SMF_status   create_L1A_granule (EN_VDATA_TYPE_t  *eng_data, 
                                      int16            nscans,
                                      MODFILE          **L1A_file_ptr);

PGSt_SMF_status   create_Vdata(char *Vdata_name, 
                               char field_names[][VU_MAX_NAME_LENGTH], 
                               char data_types[][VU_MAX_DATA_TYPE_STRING_LENGTH],
                               int16 num_fields,
                               uint16 order[]);

PGSt_SMF_status   create_Vdata_field (char    *Vdata_name,
                                      int32    Vdata_id,
                                      char    *field_name,
                                      char    *data_type,
                                      int32    order);

void   create_eng_data_vdata_array (char               *eng_data_name,
                                    EN_VDATA_TYPE_t    *eng_data,
                                    uint16             *curr_eng_data_index,
                                    uint16             *curr_field_index);

void   create_eng_data_vdata_array_field (char              *field_name,
                                          uint16            num_bits,
                                          uint16            start_bit_pos,
                                          uint16            order,
                                          uint16            type,
                                          EN_VDATA_TYPE_t   *eng_data,
                                          uint16            curr_eng_data_index,
                                          uint16            *curr_field_index);

void   create_missing_scans (int16            prev_scan_num,
                             PGSt_double      scan_rate,
                             PGSt_double      *SD_start_time,
                             MD_SCAN_MET_t    *scan_metadata,
                             SC_PIXEL_QUALITY_DATA_t   *pixel_qual_data);

char   *dequeue (FP_QUEUE_t  Q);

PGSt_SMF_status   end_Vdata_access_to_file (MODFILE  *L1A_file);

PGSt_SMF_status   end_eng_data_access_to_file (MODFILE          *L1A_file,
                                               EN_VDATA_TYPE_t  *eng_data);

PGSt_SMF_status   enqueue (char      *item,  
                           FP_QUEUE_t Q);

int16   equal_strings (char  *a,
                       char  *b);

uint32   extr_bits (uint8  *a, 
                    int    start_byte,
                    int    start_bit,
                    int    num_bits);

void   finalize_pixel_qual_data (SC_PIXEL_QUALITY_DATA_t   *scan_pixel,
                                 MD_SCAN_MET_t             *scan_meta);

void   finalize_scan_metadata (MD_SCAN_MET_t   *scan_meta,
                               int16           num_packets);

void   forget (char  *Vdata_name);

void   free_queue (FP_QUEUE_t  Q);

int16  get_empty_slot (void);

int16  get_index (char  *Vdata_name);

#define MAX_INPUTS 3
extern char   global_input_pointer[MAX_INPUTS][PGSd_PC_VALUE_LENGTH_MAX];
                         /* Names of the input files: List file name, 
                            prior L0 file name, current L0 file name  */  

#define SECTOR_TIME_OFFSETS 5
extern	PGSt_double  global_time_offset_array[SECTOR_TIME_OFFSETS];
			/* Aarray containing the time offsets for each sector */
			/* from the beginning of the scan (i.e. SD sector) */

extern  PGSt_double  global_first_gran_start_time;
                         /* Start time for the first L1A output granule */

extern  PGSt_double  global_last_gran_stop_time;
                         /* Stop time for the last L1A output granule */

int16  get_number_of_attached_Vdatas(void);

PGSt_SMF_status   get_pcf_config_data(PCF_CONFIG_t  *pcf_config);

PGSt_SMF_status   get_valid_L0_file (PGSt_tag                    spacecraft_tag,
                                     PGSt_IO_L0_VirtualDataSet  *L0_file,
                                     PGSt_double                *start_time, 
                                     PGSt_double                *stop_time);

PGSt_SMF_status   handle_missing_scans (MODFILE                   *L1A_file_ptr,
                                        PGSt_double                SD_start_time, 
                                        PGSt_double               *pre_SD_time,
                                        PGSt_double                scan_rate,
                                        int                       *scan_number,
                                        MD_ECS_GRA_INV_MET_t      *ecs_gran_meta,
                                        MD_L1A_SPECIFIC_MET_t     *L1A_specific_meta,
                                        PGSt_SMF_boolean          *gran_start_time_used,
					PGSt_double		   gran_start_time,
                                        EN_VDATA_TYPE_t           *eng_data);

PGSt_SMF_status   init_L1A_HDF_sdss (MODFILE  *L1A_file_ptr,
                                     int16    nscans);

PGSt_SMF_status   init_L1A_HDF_vdatas (EN_VDATA_TYPE_t  *eng_data,
                                       MODFILE          *L1A_file);

PGSt_SMF_status   init_L1A_pix_qual_HDF_sdss  (MODFILE  *mfile, 
                                               int      nscans);

PGSt_SMF_status   init_L1A_scan_data_HDF_sdss (MODFILE  *mfile, 
                                               int      nscans);

PGSt_SMF_status   init_L1A_scan_meta_HDF_sdss (MODFILE  *mfile, 
                                               int      nscans);

PGSt_SMF_status   initialize_global_metadata (PGSt_double             gran_start_time,
                                              PGSt_double             gran_end_time, 
                                              int                     nscans,
                                              PCF_CONFIG_t            *pcf_config,
                                              MD_ECS_GRA_INV_MET_t    *ecs_gra_inv_met,
                                              MD_L1A_SPECIFIC_MET_t   *l1a_specific_met);

void   initialize_id_table  (void);

PGSt_SMF_status   initialize_level1a (PCF_CONFIG_t               *pcf_config,
                                      EN_VDATA_TYPE_t            *eng_data,
                                      PGSt_IO_L0_Packet          *pkt,
                                      PH_PACKET_HEADER_t         *pkt_header,
                                      PGSt_IO_L0_VirtualDataSet  *L0_file);

void   initialize_pixel_qual_data (SC_PIXEL_QUALITY_DATA_t    *scan_pixel);

void   initialize_scan (SC_SCAN_DATA_t           *L1A_scan,
                        SC_PIXEL_QUALITY_DATA_t  *scan_pixel,
                        MD_SCAN_MET_t            *scan_meta);

void   initialize_scan_metadata (MD_SCAN_MET_t    *scan_meta);

void   initialize_scan_data (SC_SCAN_DATA_t    *L1A_scan);

int32   L1A_datatype_to_DFNT(char *datatype);

PGSt_SMF_status   load_eng_data (PGSt_double                scan_rate,
                                 EN_VDATA_TYPE_t           *eng_data,
                                 PGSt_IO_L0_Packet         *pkt,
                                 PH_PACKET_HEADER_t        *pkt_header,
                                 PGSt_IO_L0_VirtualDataSet  L0_file);

void   log_fmt_msg (PGSt_SMF_status  code, 
                    const char *routine, 
                    const char *msg_fmt, 
                    ...);

FP_QUEUE_t   make_queue (void);

void   output_daymode_data_to_scan (PH_PACKET_HEADER_t  *pkt_header,
                                    uint16              *pkt_contents,
                                    SC_SCAN_DATA_t      *L1A_scan);

void   output_eng1_pkt1_to_scan (PGSt_IO_L0_Packet  *pkt,
                                 SC_SCAN_DATA_t     *L1A_scan);

void   output_eng1_pkt2_to_scan (PGSt_IO_L0_Packet  *pkt,
                                 SC_SCAN_DATA_t     *L1A_scan);

void   output_eng2_pkt1_to_scan (PGSt_IO_L0_Packet  *pkt,
                                 SC_SCAN_DATA_t     *L1A_scan);

void   output_eng2_pkt2_to_scan (PGSt_IO_L0_Packet  *pkt,
                                 SC_SCAN_DATA_t     *L1A_scan);

void   output_eng_data_to_scan (PH_PACKET_HEADER_t  *pkt_header,
                                PGSt_IO_L0_Packet   *pkt,
                                SC_SCAN_DATA_t      *L1A_scan);

void   output_nightmode_data_to_scan (PH_PACKET_HEADER_t  *pkt_header,
                                      uint16              *pkt_contents,
                                      SC_EV_1KM_NIGHT      EV_1km_night);

PGSt_SMF_status packet_of_scan ( PH_PACKET_HEADER_t *pkt_header,
                                 PGSt_double next_scan_start_time,
                                 int8           *previous_scan_count,
                                 SC_SCAN_PROC_STATE_t scan_proc_state[5],
                                 PGSt_IO_L0_VirtualDataSet  *L0_file );

PGSt_SMF_status   parse_eng_data_list (EN_VDATA_TYPE_t  *eng_data);

PGSt_SMF_status   process_a_granule (PGSt_IO_L0_VirtualDataSet   L0_file,
                                     PGSt_double                 gran_start_time,
                                     PGSt_double                 gran_end_time,
                                     PCF_CONFIG_t               *pcf_config,
                                     EN_VDATA_TYPE_t            *eng_data, 
                                     PH_PACKET_HEADER_t         *pkt_header, 
                                     PGSt_IO_L0_Packet          *pkt,
                                     FP_QUEUE_t                 *failed_pkts);

PGSt_SMF_status   process_a_packet (PGSt_IO_L0_VirtualDataSet  *L0_file,
				    PGSt_IO_L0_Packet          *pkt,
				    PH_PACKET_HEADER_t         *packet_header,
				    uint16                     *packet_cont);

PGSt_SMF_status   process_a_scan (int                        *scan_number,
                                  PGSt_IO_L0_Packet          *pkt,
                                  PGSt_double                *scan_rate, 
                                  PGSt_double                *scan_time, 
                                  SC_SCAN_DATA_t             *L1A_scan,
                                  MD_SCAN_MET_t              *scan_meta,
                                  EN_VDATA_TYPE_t            *eng_data,
                                  FP_QUEUE_t                 *failed_pkts,
                                  PH_PACKET_HEADER_t         *pkt_header,
                                  SC_PIXEL_QUALITY_DATA_t    *scan_pixel,
                                  PGSt_IO_L0_VirtualDataSet  *L0_file);

void   process_cp_hk_tlmy (EN_VDATA_TYPE_t    *eng_data,
                           PGSt_IO_L0_Packet  *eng_pkt_2_1,
                           uint16             scan_number);

void   process_eng_packet (EN_VDATA_TYPE_t       *eng_data, 
                           int                    scan_number, 
                           PH_PACKET_HEADER_t    *pkt_header,
                           PGSt_IO_L0_Packet     *pkt);

void   process_group2_packet1_vdata (PGSt_IO_L0_Packet *pkt, EN_VDATA_TYPE_t *eng_data); 

PGSt_SMF_status  process_next_packet(
	PGSt_IO_L0_VirtualDataSet	*L0_file,
	PH_PACKET_HEADER_t		**packet_header);

void   process_sci_eng_data (EN_VDATA_TYPE_t    *eng_data,
                             PGSt_IO_L0_Packet  *eng_pkt_1_2,
                             uint16             scan_number);

void   put_cal_data_in_scan (PH_PACKET_HEADER_t   *pkt_header,
                             uint16               *pkt_cont,
                             SC_CAL_250M          SC_250m,
                             SC_CAL_500M          SC_500m,
                             SC_CAL_1KM_DAY       SC_1km_day, 
                             SC_CAL_1KM_NIGHT     SC_1km_night);

void   put_earth_data_in_scan (PH_PACKET_HEADER_t   *pkt_header,
                               uint16               *pkt_cont,
                               SC_EV_250M           SC_250m,
                               SC_EV_500M           SC_500m,
                               SC_EV_1KM_DAY        SC_1km_day, 
                               SC_EV_1KM_NIGHT      SC_1km_night);

void   put_pkt_cont_in_scan (PH_PACKET_HEADER_t   pkt_header,
                             PGSt_IO_L0_Packet    *pkt,
                             uint16               *pkt_contents,
                             SC_SCAN_DATA_t       *L1A_scan);

PGSt_SMF_status   read_a_packet (PGSt_IO_L0_VirtualDataSet  *L0_file,
                                 PGSt_IO_L0_Packet          *pkt);

int32  recall_id (char  *Vdata_name);

void   remember (char    *Vdata_name,
                 int32    Vdata_id);

void   reset_last_valid_scan (EN_VDATA_TYPE_t  *eng_data);

PGSt_SMF_status   set_start_position (PGSt_IO_L0_VirtualDataSet  L0_file,
                                      PGSt_double                *preload_start_time,
                                      PGSt_double                start_time,
                                      PGSt_double                stop_time);

PGSt_SMF_status   unpack_MODIS_header (PGSt_IO_L0_Packet   *pkt,
                                       PH_PACKET_HEADER_t  *packet_header);

void   unpack_packet_contents (PGSt_IO_L0_Packet   *pkt,
                               PH_PACKET_HEADER_t  *pkt_header,
                               uint16              *pkt_contents);

PGSt_SMF_status   unpack_packet_header (PGSt_IO_L0_Packet   *pkt,
                                        PH_PACKET_HEADER_t  *packet_header);

PGSt_SMF_status   unpack_primary_header (PGSt_IO_L0_Packet   *pkt,
                                         PH_PACKET_HEADER_t  *packet_header);

PGSt_SMF_status   unpack_secondary_header (PGSt_IO_L0_Packet   *pkt,
                                           PH_PACKET_HEADER_t  *packet_header);

void   update_eng_data (uint16              index,
                        PGSt_IO_L0_Packet   *eng_packet,
                        uint16              scan_number,
                        EN_VDATA_TYPE_t     *eng_data, 
                        int                 use_cp_prior_offset);

void   update_eng_data_for_maj_cycle_n (uint16             major_cycle,
                                        PGSt_IO_L0_Packet  *eng_pkt_2_1,
                                        uint16             scan_number,
                                        EN_VDATA_TYPE_t    *eng_data,
                                        int                is_cp_hk_prior_section);

void   update_global_metadata (MD_SCAN_MET_t            *scan_meta,
                               MD_ECS_GRA_INV_MET_t     *ecs_gra_inv_met,
                               MD_L1A_SPECIFIC_MET_t    *l1a_specific_met);

void   update_pixel_qual_data (PH_PACKET_HEADER_t       pkt_header,
                               int16                    qual_value,           
                               SC_PIXEL_QUALITY_DATA_t  *scan_pix);

void   update_scan_metadata (PH_PACKET_HEADER_t  packet_header,
                             PGSt_SMF_status     packetStatus,
                             MD_SCAN_MET_t       *scan_meta,
                             int16               *qual_value);

int   validate_L0_header (PGSt_IO_L0_VirtualDataSet L0_file);

PGSt_SMF_status   write_ECS_metadata (PGSt_MET_all_handles  md_handles,
                                      MD_ECS_GRA_INV_MET_t  *ecs_gra_inv_met);

PGSt_SMF_status   write_Vdata (char            *Vdata_name,
                               unsigned char   *data,
                               int32           num_records);

PGSt_SMF_status   write_eng_data (EN_VDATA_TYPE_t  *eng_data);

PGSt_SMF_status   write_failed_packets (FP_QUEUE_t failed_pkts);

PGSt_SMF_status   write_global_metadata (MODFILE                *mfile,
                                         PGSt_MET_all_handles   md_handles,
                                         MD_ECS_GRA_INV_MET_t   *ecs_gra_inv_met,
                                         MD_L1A_SPECIFIC_MET_t  *l1a_specific_met);

PGSt_SMF_status   write_pix_qual (MODFILE                   *L1A_file_ptr, 
                                  SC_PIXEL_QUALITY_DATA_t   *pix_qual,
                                  int16                     scan_num);

PGSt_SMF_status   write_scan (MODFILE                  *L1A_file_ptr, 
                             SC_SCAN_DATA_t           *L1A_scan, 
                             MD_SCAN_MET_t            *scan_meta, 
                             SC_PIXEL_QUALITY_DATA_t  *pix_qual, 
                             FP_QUEUE_t                failed_pkts, 
                             EN_VDATA_TYPE_t          *eng_data);

PGSt_SMF_status   write_scan_data (MODFILE         *L1A_file_ptr, 
                                   SC_SCAN_DATA_t  *L1A_scan,
                                   int16           scan_num,
                                   char            *scan_type);

PGSt_SMF_status   write_scan_metadata (MODFILE        *L1A_file_ptr, 
                                       MD_SCAN_MET_t  *scan_meta);

PGSt_SMF_status   write_specific_granule_metadata (MODFILE                *mfile,
                                                   MD_L1A_SPECIFIC_MET_t  *l1a_specific_met);

#endif

