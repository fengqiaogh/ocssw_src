#ifndef _INSTLM_H
#define _INSTLM_H

#include "swl0_types.h"
#include "genutils.h"

static int instlm_list[3][5] = {

    {1, 0, 1, 1, 0}, /* mnf #1 */
    {0, 1, 1, 0, 1}, /* mnf #2 */
    {1, 1, 0, 1, 1} /* mnf #3 */
};


static float instlm_limits[32][4] = {

    { -1.334, 66.667, 5.0, 35.0}, /* Band 1/2 FPA Temperature */
    { -1.334, 66.667, 5.0, 35.0}, /* Band 3/4 FPA Temperature */
    { -1.334, 66.667, 5.0, 35.0}, /* Band 5/6 FPA Temperature */
    { -1.334, 66.667, 5.0, 35.0}, /* Band 7/8 FPA Temperature */
    { -1.334, 66.667, 5.0, 35.0}, /* Telescope Motor Temperature */
    { -1.334, 66.667, 4.0, 35.0}, /* Tilt Base Temperature */
    { -1.334, 66.667, 5.0, 35.0}, /* Tilt Platform temperature */
    { -1.334, 66.667, 5.0, 35.0}, /* Half-Angle Motor Temperature */
    { 0.26, 5.36, 1.0, 3.2}, /* Power Supply-A Input Current */
    { 0.26, 5.36, 1.0, 3.2}, /* Power Supply-B Input Current */
    { 0.0, 19.125, 14.0, 16.0}, /* +15 V Analog Power Voltage */
    {-19.125, 0.0, -16.0, -14.0}, /* -15 V Analog Power Voltage */
    { 0.0, 6.375, 4.9, 5.7}, /* +5 V Logic Power Voltage */
    { -1.334, 66.667, 5.0, 50.0}, /* Power Supply Temperature */
    { -1.334, 66.667, 5.0, 35.0}, /* B1/B2 post-Amplifier Temperature */
    { -1.334, 66.667, 5.0, 50.0}, /* Servo Driver Temperature */
    { 0.0, 38.25, 28.5, 31.0}, /* +30 V Servo Power Voltage */
    { 0.0, 26.622, 19.5, 22.0}, /* +21 V Servo Power Voltage */
    {-26.622, 0.0, -22.0, -19.5}, /* -21 V Servo Power Voltage */
    { 0.0, 6.375, 4.9, 5.4}, /* +5 V Servo Power Voltage */
    { -377.0, 1795.6, 1165.0, 1325.0}, /* Angular Momentum Speed */
    { 0.0, 367.2, 40.0, 240.0}, /* Tilt Platform Position */
    { 0.0, 367.2, 40.0, 240.0}, /* Tilt Base Position */
    { 0.0, 35.7, 24.0, 32.0}, /* +28 V Heater Power */
    { 0.0, 0.612, 0.0, 0.4}, /* Telescope-A Motor Current */
    { 0.0, 0.612, 0.0, 0.4}, /* Telescope-B Motor Current */
    { 0.0, 0.612, 0.0, 0.3}, /* Half-Angle-A Motor Current */
    { 0.0, 0.612, 0.0, 0.3}, /* Half-Angle-B  Motor Current */
    { -1.25, 1.25, -1.3, 1.3}, /* Servo-A Phase Error */
    { -1.25, 1.25, -1.3, 1.3}, /* Servo-B  Phase Error */
    { 0.0, 4.08, 0.0, 0.6}, /* Ang-Mom-Comp-A Motor Current */
    { 0.0, 4.08, 0.0, 0.6} /* Ang-Mom-Comp-B Motor Current */
};


#endif

