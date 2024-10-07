cmake_minimum_required(VERSION 3.12)

set(LIBS timeutils)

add_executable(l0qc_viirs
    l0qc_viirs.c
    init_packet_read.c
    read_packet.c
    get_viirs_packet_time.c
    read_viirs_scan_packets.c
)

add_executable(scqc_viirs
    scqc_viirs.c
)

add_executable(snpp_diary_to_llh
	snpp_diary_to_llh.c
	convert_diary.c
	orb2lla.c
)

target_link_libraries(l0qc_viirs ${LIBS})

target_link_libraries(scqc_viirs ${LIBS})

target_link_libraries(snpp_diary_to_llh ${LIBS})


install (TARGETS l0qc_viirs scqc_viirs snpp_diary_to_llh DESTINATION bin)

