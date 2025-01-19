#include "device.hpp"
#include <stdexcept>

std::string determineColor(Device device) {
    switch (device) {
        case BLUE:
            return "blue";
        case RED:
            return "red";
        case SWIR:
            return "SWIR";
        default:
            throw std::invalid_argument("Device doesn't exist");
    }
}