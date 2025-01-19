
#ifndef __OCI_DEVICE_H__
#define __OCI_DEVICE_H__

#include <string>

enum Device {
    BLUE,
    RED,
    SWIR,
    ENUM_DEVICE_SIZE  // Sentinel value just in case the size of this enum is needed
};

std::string determineColor(Device device);

#endif