
# l1agen_modis Changelog

## <VERSION> - 2015-03-12

 - isolated MODIS extract metadata code
 - removed extra libs from modis CMakeLists.txt
 - Collection "6.1" of MODIS L1 processing code: l1agen_modis v6.0.6, geogen_modis v6.1.0, l1bgen_modisa v6.2.1, l1bgen_modist v6.2.2.
 - Cleaned up our changes to MODIS L1 code
 - Updated MODIS L1 code to be compatible with opt (aka lib3) r14227, and simplified changes wrt SDST delivery so as to make future integration easier.
 - moving MODIS L1 files out of subdirectory
  
### Source Changed

  * accumulate_failed_packets.c
  * attached_Vdata_counter.c
  * Changelog.md
  * check_checksum.c
  * close_processing_run.c
  * close_Vdata.c
  * CMakeLists.txt
  * compute_global_time_offsets.c
  * compute_SD_start_time.c
  * create_eng_data_vdata_array.c
  * create_eng_data_vdata_array_field.c
  * create_L1A_granule.c
  * create_missing_scans.c
  * create_Vdata.c
  * create_Vdata_field.c
  * dequeue.c
  * end_eng_data_access_to_file.c
  * end_Vdata_access_to_file.c
  * EN_eng_data.h
  * enqueue.c
  * equal_strings.c
  * extr_bits.c
  * finalize_pixel_qual_data.c
  * finalize_scan_metadata.c
  * forget.c
  * FP_failed_pkt_queue.h
  * free_queue.c
  * from_SDST/HISTORY.txt
  * from_SDST/MOD01.mcf
  * from_SDST/MOD_PR01.mk
  * from_SDST/MOD_PR01.pcf
  * from_SDST/MOD_PR01_pr.txt
  * from_SDST/MYD01.mcf
  * from_SDST/MYD_PR01.pcf
  * from_SDST/PACKINGLIST.txt
  * from_SDST/README.txt
  * get_empty_slot.c
  * get_index.c
  * get_number_of_attached_Vdatas.c
  * get_pcf_config_data.c
  * get_valid_L0_file.c
  * handle_missing_scans.c
  * HISTORY.txt
  * initialize_global_metadata.c
  * initialize_id_table.c
  * initialize_level1a.c
  * initialize_pixel_qual_data.c
  * initialize_scan.c
  * initialize_scan_data.c
  * initialize_scan_metadata.c
  * init_L1A_HDF_sdss.c
  * init_L1A_HDF_vdatas.c
  * init_L1A_pix_qual_HDF_sdss.c
  * init_L1A_scan_data_HDF_sdss.c
  * init_L1A_scan_meta_HDF_sdss.c
  * L1A_datatype_to_DFNT.c
  * L1A_prototype.h
  * level1a.c
  * load_eng_data.c
  * log_fmt_msg.c
  * make_queue.c
  * MD_metadata.h
  * MOD01.mcf
  * MODIS_35005.t
  * MOD_PR01.mk
  * MOD_PR01.pcf
  * MOD_PR01_pr.txt
  * MS_misc.h
  * MYD01.mcf
  * MYD_PR01.pcf
  * output_daymode_data_to_scan.c
  * output_eng1_pkt1_to_scan.c
  * output_eng1_pkt2_to_scan.c
  * output_eng2_pkt1_to_scan.c
  * output_eng2_pkt2_to_scan.c
  * output_eng_data_to_scan.c
  * output_nightmode_data_to_scan.c
  * packet_of_scan.c
  * packet_stats.h
  * PACKINGLIST.txt
  * parse_eng_data_list.c
  * PC_pcf_info.h
  * PD_pkt_data.h
  * PGS_MODIS_35005.h
  * PH_pkt_hdr.h
  * print_stats.c
  * process_a_granule.c
  * process_a_packet.c
  * process_a_scan.c
  * process_cp_hk_tlmy.c
  * process_eng_packet.c
  * process_group2_packet1_vdata.c
  * process_sci_eng_data.c
  * put_cal_data_in_scan.c
  * put_earth_data_in_scan.c
  * put_pkt_cont_in_scan.c
  * read_a_packet.c
  * README.txt
  * recall_id.c
  * remember.c
  * reset_last_valid_scan.c
  * SC_scan.h
  * set_start_position.c
  * unpack_MODIS_header.c
  * unpack_packet_contents.c
  * unpack_packet_header.c
  * unpack_primary_header.c
  * unpack_secondary_header.c
  * update_eng_data.c
  * update_eng_data_for_maj_cycle_n.c
  * update_global_metadata.c
  * update_pixel_qual_data.c
  * update_scan_metadata.c
  * validate_L0_header.c
  * version.h
  * VU_vdata_utility.h
  * write_ECS_metadata.c
  * write_eng_data.c
  * write_failed_packets.c
  * write_global_metadata.c
  * write_pix_qual.c
  * write_scan.c
  * write_scan_data.c
  * write_scan_metadata.c
  * write_specific_granule_metadata.c
  * write_Vdata.c
## <VERSION> - 2015-03-12
### Added
  * Changelog.md
