
macro(SET_OCSSW_POLICY)
  # set the <package>_ROOT policy
  if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
  endif()
endmacro(SET_OCSSW_POLICY)

