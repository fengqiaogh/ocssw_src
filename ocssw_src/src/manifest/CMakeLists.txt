cmake_minimum_required(VERSION 3.12)

install(PROGRAMS manifest.py update_manifest_dirs build_new_tag fix_mac_rpath.py DESTINATION bin)
install(PROGRAMS manifest_ocssw.py  RENAME manifest_ocssw DESTINATION bin)
install(PROGRAMS install_ocssw.py  RENAME install_ocssw DESTINATION bin)
install(PROGRAMS build_manifest_ocssw.py RENAME build_manifest_ocssw DESTINATION bin)
