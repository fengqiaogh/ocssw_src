#ifndef OCSSW_CACHESIZE_H
#define OCSSW_CACHESIZE_H
#include "OutFile.h"
#include "L3File.h"
#include "Wave3DParsing.hpp"

size_t get_chunk_along_lines();

bool set_cache_size_chunk_size_read(const std::string &fileName, bool verbose = false);
bool set_cache_size_chunk_size_write(l3::L3File *l3file, OutFile *file, OutFile *file2,  bool verbose = false);

#endif  // OCSSW_CACHESIZE_H
