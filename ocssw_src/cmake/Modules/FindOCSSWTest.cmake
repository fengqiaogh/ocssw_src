# - Try to find the OCSSW test directory
# Once done this will define
#  OCSSWTestDir_FOUND - has the test tree

find_path(OCSSWTestDir_FOUND README.OCSSW_TESTS
  HINTS $ENV{OCSSWROOT}/test
)
