
#ifndef L1B_OCI_OPTIONS_H
#define L1B_OCI_OPTIONS_H

#include "device.hpp"
#include <clo.h>

namespace oel {
class L1bOptions {
   public:
    L1bOptions(int argc, char *argv[]);
    ~L1bOptions();

    std::string l1aFilename;  // Input filename
    std::string l1bFilename;  // Output filename
    std::string calibrationLutFilename;
    std::string geolocationLutFilename;
    std::string demFile;  // Digital Elevation Model
    std::string digitalObjectId;
    std::string processingVersion;
    std::string ephFile;           // Definititive ephemeris file, used for geolocation
    std::string xtalkLutFilename;  // Cross-talk LUT filename

    bool radianceGenerationEnabled;
    bool disableGeolocation;
    bool enableCrosstalk;  // Used in tandem with xtalkLutFilename. Latter has default value if not set
                           // explicitly, so never empty.

    size_t startingLine;  // start line (first line = 0)
    int32_t endingLine;   // ending line inclusive (first line = 0)

    int deflate; // Compression factor

    clo_optionList_t *optionList;
};
}  // namespace oel

#endif
