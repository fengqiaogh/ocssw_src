
#ifndef __OCI_DEVICE_H__
#define __OCI_DEVICE_H__

#include <string>

enum Device {
    BLUE,
    RED,
    SWIR,
    CCD,              // Either red or blue
    ENUM_DEVICE_SIZE  // Sentinel value just in case the size of this enum is needed
};

/**
 * @brief Indicates which celestial body is being observed
 *
 * @note Defined in device.hpp because it's needed by multiple files; doesn't pertain to a device specifically
 */
enum LocatingContext { GEO, HELIO, SELENO, NONE };

std::string determineColor(Device device);

#endif