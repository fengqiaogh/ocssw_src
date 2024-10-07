int read_extract_metadata(int32 sd_id,
                          int32 *extract_pixel_offset,
                          int32 *extract_pixel_count,
                          int32 *extract_line_offset,
                          int32 *extract_line_count);
int write_extract_metadata(int32 sd_id,
                           int32 extract_pixel_offset,
                           int32 extract_pixel_count,
                           int32 extract_line_offset,
                           int32 extract_line_count);
