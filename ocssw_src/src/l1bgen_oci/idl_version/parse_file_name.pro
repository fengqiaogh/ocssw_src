  function parse_file_name,file_path_name
; IDL function to parse a file name from a path/name string

  fpstr = strsplit(file_path_name, '/', /extract)
  npf = n_elements(fpstr)
  return, fpstr[npf-1]
  end
