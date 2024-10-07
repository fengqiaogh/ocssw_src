#include "mfhdf.h"

int read_extract_metadata(int32 sd_id,
                          int32 *extract_pixel_offset,
                          int32 *extract_pixel_count,
                          int32 *extract_line_offset,
                          int32 *extract_line_count) {

    int32 attr_index;
    *extract_pixel_offset = -1;
    *extract_pixel_count = -1;
    *extract_line_offset = -1;
    *extract_line_count = -1;

    attr_index = SDfindattr(sd_id, "Extract Pixel Offset");
    if (attr_index != FAIL)
        SDreadattr(sd_id, attr_index, (VOIDP) extract_pixel_offset);

    attr_index = SDfindattr(sd_id, "Extract Pixel Count");
    if (attr_index != FAIL)
        SDreadattr(sd_id, attr_index, (VOIDP) extract_pixel_count);

    attr_index = SDfindattr(sd_id, "Extract Line Offset");
    if (attr_index != FAIL)
        SDreadattr(sd_id, attr_index, (VOIDP) extract_line_offset);

    attr_index = SDfindattr(sd_id, "Extract Line Count");
    if (attr_index != FAIL)
        SDreadattr(sd_id, attr_index, (VOIDP) extract_line_count);

    return SUCCEED;
}

int write_extract_metadata(int32 sd_id,
                           int32 extract_pixel_offset,
                           int32 extract_pixel_count,
                           int32 extract_line_offset,
                           int32 extract_line_count) {

    if (extract_pixel_offset != -1)
        SDsetattr(sd_id, "Extract Pixel Offset",
                  DFNT_INT32, 1,
                  (VOIDP) &extract_pixel_offset);

    if (extract_pixel_count != -1)
        SDsetattr(sd_id, "Extract Pixel Count",
                  DFNT_INT32, 1,
                  (VOIDP) &extract_pixel_count);

    if (extract_line_offset != -1)
        SDsetattr(sd_id, "Extract Line Offset",
                  DFNT_INT32, 1,
                  (VOIDP) &extract_line_offset);

    if (extract_line_count != -1)
        SDsetattr(sd_id, "Extract Line Count",
                  DFNT_INT32, 1,
                  (VOIDP) &extract_line_count);

    return SUCCEED;
}
