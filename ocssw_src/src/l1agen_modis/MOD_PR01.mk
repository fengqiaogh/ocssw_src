.POSIX:
	
###########################################################################
# !make
#
# !Makefile Name: MOD_PR01.mk
#
# !Description:  Makefile for MOD_PR01 Version 2 Level 1A code.
#                This makefile is used to generate MOD_PR01.exe
#
# !Variables:
#
#      VARIABLES     DESCRIPTION
#      ~~~~~~~~~     ~~~~~~~~~~~
#      TARGET        Program executable name
#      ADD_CFLAGS    Additional compiler options
#      LIB           Libraries
#      MATH          Math Libraries
#      SRCS          Source files
#      OBJ           Object files
#      INC           Include files
#
# !Env Variables
#
#      ENV VARIABLES DESCRIPTION
#      ~~~~~~~~~~~~~ ~~~~~~~~~~~
#      ~~~~~~~~~~~~~ ~~~~~~~~~~~
#      CFLAGS        Compiler flags which are set by the script
#                    (e.g. v2_n32_f77)
#      CC            The compiler set by the script
#      API_INC       Include directory of MODIS Appl Prog Interface set
#                    by the script
#      SDST_INC      Include directory of MODIS SDST Toolkit set by the 
#                    script
#      PGSINC        Include directory of PGS toolkit set by the script
#      HDFINC        Include directory of HDF 4 set by the script
#      HDF5INC       Include directory of HDF 5 set by the script (optional)
#      HDFEOS_INC    Include directory of HDFEOS 4 set by the script
#      API_LIB       Library of MODIS API set by the script
#      SDST_LIB      Library of MODIS SDST Toolkit set by the script
#      PGSLIB        Library of PGS toolkit set by the script
#      HDFLIB        Library of HDF 4 set by the script
#      HDF5LIB	     Library of HDF 5 set by the script (optional)
#      HDFEOS_LIB    Library of HDFEOS 4 set by the script
#
# !Team-Unique Header:
#      This software was developed by MODIS Science Data Support Team
#      for the National Aeronautics and Space Administration,
#      Goddard Space Flight Center, under contract NAS5-32373.
#
# !Revision History:
# $Log: MOD_PR01.mk,v $
# Revision 6.3  2015/04/14 20:45:20  jkuyper
# Dropped LIB_HDF5OPT, which was only needed when we were linking the shared-object version.
#
# Revision 6.2  2011/11/02 21:47:46  kuyper
# Changes to enforce POSIX conformance.
#
# Revision 6.1  2010/08/25 16:27:32  kuyper
# Updated for compatibilty with SDP Toolkit 5.2.16, while retaining
#   compatibility with older versions.
# Changed to link to HDF5 library statically, rather than dynamically.
# Changed to turn on POSIX features, thereby indirectly turning off XOPEN
#   extensions that are inappropriately turned on by the MODIS project's use
#   of -xansi option.
# Added print_stats.
#
# Revision 5.1  2004/10/04 19:30:15  seaton
# Revised for Collection 5.
# seaton@saicmodis.com
#
# Revision 4.4  2003/05/12 22:30:30  kuyper
# Corrected to match the macro name restrictions in the Software Delivery
#   Guide.
#
# Revision 4.2  2003/03/11 19:31:09  kuyper
# Changed to optionally make use of HDF 5.
#
# James Kuyper Jr. (kuyper@saicmodis.com)
#
# revision 1.0 1997/12/09   9:45:00
# Tom Johnson/GSC  (johnson@ltpmail.gsfc.nasa.gov)
# Original makefile
#
# !Note:
#       1) must follow MAPI, SDP TK, HDF, and other header/libs order
#       2) HDFLIB must be explicit and in order of -lmfhdf, -ldf,-ljpeg, -lz.
#       3) system environments, such as APIINC, PGSINC, HDFINC, PGSLIB,
#          API_LIB etc. should be defined before executing this makefile.
#
# !END
##############################################################################

# Define executable name

TARGET = MOD_PR01.exe

# Combine the includes using pre-defined includes and your includes

INC_HDF5OPT = -I$(HDF5INC)
INC =      -I$(API_INC) \
           -I$(SDST_INC) \
	   $(INC_HDF5OPT:-I=) \
           -I$(PGSINC) \
           -I$(HDFINC)


# Combine the libraries using pre-defined libraries
LIB_DL = $(BRAND:linux=-ldl)$(BRAND:sgi64=-lc)
LIB =      -L$(API_LIB) -lmapi \
           -L$(SDST_LIB) -lsdst \
           -L$(PGSLIB) -lPGSTK -lhdf5 -L$(HDFLIB) -lmfhdf -ldf	\
	   $${SZIPHOME:+-lsz} $${SZIPHOME:+$(LIB_DL:$(BRAND)=)} -ljpeg -lz -lm


# Define additional C flags

         # Flags for optimization run
         ADD_CFLAGS = -D_POSIX_C_SOURCE=1 $(BRAND:linux=-ffloat-store)

         # Flags for debugging
         # Need to also remove the O3 from the .c.o section below. 
         # ADD_CFLAGS = -woff 1209  -O0 -fullwarn -g

# Define object files

OBJ = \
      accumulate_failed_packets.o           \
      attached_Vdata_counter.o              \
      check_checksum.o                      \
      close_processing_run.o                \
      close_Vdata.o                         \
      compute_SD_start_time.o               \
      compute_global_time_offsets.o         \
      create_L1A_granule.o                  \
      create_Vdata.o                        \
      create_Vdata_field.o                  \
      create_eng_data_vdata_array.o         \
      create_eng_data_vdata_array_field.o   \
      create_missing_scans.o                \
      dequeue.o                             \
      end_Vdata_access_to_file.o            \
      end_eng_data_access_to_file.o         \
      enqueue.o                             \
      equal_strings.o                       \
      extr_bits.o                           \
      finalize_pixel_qual_data.o            \
      finalize_scan_metadata.o              \
      forget.o                              \
      free_queue.o                          \
      get_empty_slot.o                      \
      get_index.o                           \
      get_number_of_attached_Vdatas.o       \
      get_pcf_config_data.o                 \
      get_valid_L0_file.o                   \
      handle_missing_scans.o                \
      init_L1A_HDF_sdss.o                   \
      init_L1A_HDF_vdatas.o                 \
      init_L1A_pix_qual_HDF_sdss.o          \
      init_L1A_scan_data_HDF_sdss.o         \
      init_L1A_scan_meta_HDF_sdss.o         \
      initialize_global_metadata.o          \
      initialize_id_table.o                 \
      initialize_level1a.o                  \
      initialize_pixel_qual_data.o          \
      initialize_scan.o                     \
      initialize_scan_data.o                \
      initialize_scan_metadata.o            \
      L1A_datatype_to_DFNT.o                \
      level1a.o                             \
      load_eng_data.o                       \
      log_fmt_msg.o                         \
      make_queue.o                          \
      output_daymode_data_to_scan.o         \
      output_eng1_pkt1_to_scan.o            \
      output_eng1_pkt2_to_scan.o            \
      output_eng2_pkt1_to_scan.o            \
      output_eng2_pkt2_to_scan.o            \
      output_eng_data_to_scan.o             \
      output_nightmode_data_to_scan.o       \
      packet_of_scan.o                      \
      parse_eng_data_list.o                 \
      print_stats.o			    \
      process_a_granule.o                   \
      process_a_packet.o                    \
      process_a_scan.o                      \
      process_cp_hk_tlmy.o                  \
      process_eng_packet.o                  \
      process_group2_packet1_vdata.o        \
      process_sci_eng_data.o                \
      put_cal_data_in_scan.o                \
      put_earth_data_in_scan.o              \
      put_pkt_cont_in_scan.o                \
      read_a_packet.o                       \
      recall_id.o                           \
      remember.o                            \
      reset_last_valid_scan.o               \
      set_start_position.o                  \
      unpack_MODIS_header.o                 \
      unpack_packet_contents.o              \
      unpack_packet_header.o                \
      unpack_primary_header.o               \
      unpack_secondary_header.o             \
      update_eng_data.o                     \
      update_eng_data_for_maj_cycle_n.o     \
      update_global_metadata.o              \
      update_pixel_qual_data.o              \
      update_scan_metadata.o                \
      validate_L0_header.o                  \
      write_ECS_metadata.o                  \
      write_Vdata.o                         \
      write_eng_data.o                      \
      write_failed_packets.o                \
      write_global_metadata.o               \
      write_pix_qual.o                      \
      write_scan.o                          \
      write_scan_data.o                     \
      write_scan_metadata.o                 \
      write_specific_granule_metadata.o


# Make the process

$(TARGET)      :       $(OBJ)
	$(CC)   $(CFLAGS) $(OBJ) $(LIB) -o $(TARGET)

.c.o   :
	$(CC) -O3  $(CFLAGS) $(ADD_CFLAGS:$(BRAND)=) $(INC) -c  $<



#******************** End of make file ********************
