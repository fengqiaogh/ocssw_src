#
# GetGitRevisionDescription
#
# Derive version information using the SHA1 from the latest git commit
#
include(GetGitRevisionDescription)
git_describe(VERSION --tags)

#parse the git version to get just the latest SHA
string (REGEX REPLACE  "^(.*-g)(.*$)" "\\2" GITSHA "${VERSION}")
set(VERSION_SHORT "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

configure_file(${CMAKE_CURRENT_LIST_DIR}/version.h.in
                ${CMAKE_CURRENT_SOURCE_DIR}/version.h)
set(version_file "${CMAKE_CURRENT_SOURCE_DIR}/version.h")

