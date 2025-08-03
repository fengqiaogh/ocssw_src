#ifndef GEO_REGION_H
#define GEO_REGION_H

#ifdef __cplusplus
extern "C" {
#endif

void set_georegion_filename(const char* filename);
int get_georegion(float lat, float lon);
void close_georegion_file();

#ifdef __cplusplus
}

// C++ section of header
#include <string>

void set_georegion_filename(const std::string& filename);

#endif

#endif
