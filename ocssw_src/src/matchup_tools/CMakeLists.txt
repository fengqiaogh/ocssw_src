cmake_minimum_required(VERSION 3.12)

install(PROGRAMS fd_matchup.py DESTINATION bin RENAME find_matchup)
install(PROGRAMS mk_matchup.py DESTINATION bin RENAME make_matchup)
install(PROGRAMS mk_matchup_sst.py DESTINATION bin RENAME make_sst_matchup)
install(PROGRAMS seabass2L1B.py DESTINATION bin)
install(DIRECTORY seabass DESTINATION bin)
