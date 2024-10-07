cmake_minimum_required(VERSION 3.12)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
set(LIBS ${HDF4_LIBRARIES})

find_package(Mapi REQUIRED)
include_directories(${Mapi_INCLUDE_DIRS})
list(APPEND LIBS ${Mapi_LIBRARIES})

find_package(SDST REQUIRED)
include_directories(${SDST_INCLUDE_DIRS})
list(APPEND LIBS ${SDST_LIBRARIES})

add_executable(l1agen_modis
  accumulate_failed_packets.c
  attached_Vdata_counter.c
  check_checksum.c
  close_processing_run.c
  close_Vdata.c
  compute_SD_start_time.c
  compute_global_time_offsets.c
  create_L1A_granule.c
  create_Vdata.c
  create_Vdata_field.c
  create_eng_data_vdata_array.c
  create_eng_data_vdata_array_field.c
  create_missing_scans.c
  dequeue.c
  end_Vdata_access_to_file.c
  end_eng_data_access_to_file.c
  enqueue.c
  equal_strings.c
  extr_bits.c
  finalize_pixel_qual_data.c
  finalize_scan_metadata.c
  forget.c
  free_queue.c
  get_empty_slot.c
  get_index.c
  get_number_of_attached_Vdatas.c
  get_pcf_config_data.c
  get_valid_L0_file.c
  handle_missing_scans.c
  init_L1A_HDF_sdss.c
  init_L1A_HDF_vdatas.c
  init_L1A_pix_qual_HDF_sdss.c
  init_L1A_scan_data_HDF_sdss.c
  init_L1A_scan_meta_HDF_sdss.c
  initialize_global_metadata.c
  initialize_id_table.c
  initialize_level1a.c
  initialize_pixel_qual_data.c
  initialize_scan.c
  initialize_scan_data.c
  initialize_scan_metadata.c
  L1A_datatype_to_DFNT.c
  level1a.c
  load_eng_data.c
  log_fmt_msg.c
  make_queue.c
  output_daymode_data_to_scan.c
  output_eng1_pkt1_to_scan.c
  output_eng1_pkt2_to_scan.c
  output_eng2_pkt1_to_scan.c
  output_eng2_pkt2_to_scan.c
  output_eng_data_to_scan.c
  output_nightmode_data_to_scan.c
  packet_of_scan.c
  parse_eng_data_list.c
  print_stats.c
  process_a_granule.c
  process_a_packet.c
  process_a_scan.c
  process_cp_hk_tlmy.c
  process_eng_packet.c
  process_group2_packet1_vdata.c
  process_sci_eng_data.c
  put_cal_data_in_scan.c
  put_earth_data_in_scan.c
  put_pkt_cont_in_scan.c
  read_a_packet.c
  recall_id.c
  remember.c
  reset_last_valid_scan.c
  set_start_position.c
  unpack_MODIS_header.c
  unpack_packet_contents.c
  unpack_packet_header.c
  unpack_primary_header.c
  unpack_secondary_header.c
  update_eng_data.c
  update_eng_data_for_maj_cycle_n.c
  update_global_metadata.c
  update_pixel_qual_data.c
  update_scan_metadata.c
  validate_L0_header.c
  write_ECS_metadata.c
  write_Vdata.c
  write_eng_data.c
  write_failed_packets.c
  write_global_metadata.c
  write_pix_qual.c
  write_scan.c
  write_scan_data.c
  write_scan_metadata.c
  write_specific_granule_metadata.c
)

target_link_libraries(l1agen_modis
  ${LIBS}
)

install (TARGETS l1agen_modis DESTINATION bin)
