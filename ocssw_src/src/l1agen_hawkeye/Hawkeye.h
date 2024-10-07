/*

        Hawkeye.h

        Contains the Constants, Enumerations and Structure definitions
        for communicating with the Hawkeye Camera.
        
*/

#ifndef _HAWKEYE_
#define _HAWKEYE_

#include <stdint.h>

#define HAWKEYE_CANDC_PORT                              5625                    // TCP/IP port for Command and Control
#define HAWKEYE_IMAGE_PORT                              5626                    // TCP/IP port for Image Data

#define FINDERSCOPE_VBINNING                    1                               // set to 1 (unbinned) or 2 for 2:1 in-chip vertical binning
#define FINDERSCOPE_FPGA_BINNING                0                               // when the above is 2, set this to 1 to bin in the FPGA (sum 2 rows),
                                                                                                                //  or set it to 0 to bin in the Micron Sensor (weird, part average, part min value)

#define FINDERSCOPE_LIGHT_WIDTH                 752                             // Width of the Finderscope Sensor
#define FINDERSCOPE_LIGHT_HEIGHT                (480/FINDERSCOPE_VBINNING)                      // Height of the Finderscope Sensor with 2:1 vertical binning
#define FINDERSCOPE_CAPTURE_DARK_WIDTH  36                              // Width of Captured Dark at start of each Micron row
#define FINDERSCOPE_DOWNLOAD_DARK_WIDTH 4                               // Width of Downloaded Dark at start of each Micron row
#define FINDERSCOPE_CAPTURE_WIDTH               (FINDERSCOPE_CAPTURE_DARK_WIDTH + FINDERSCOPE_LIGHT_WIDTH)
#define FINDERSCOPE_DOWNLOAD_WIDTH              (FINDERSCOPE_DOWNLOAD_DARK_WIDTH + FINDERSCOPE_LIGHT_WIDTH)
#define FINDERSCOPE_MAX_IMAGES                  25                            // Maximum number of Finderscope Images per Exposure
#if FINDERSCOPE_VBINNING == 1
  #define FINDERSCOPE_MAX_IMAGES_UGA    250                             //  "        " for UGA
  #define FINDERSCOPE_MIN_PERIOD_MS             18                              // Minimum finderscope period in ms
#else
  #define FINDERSCOPE_MAX_IMAGES_UGA    500                             //  "        " for UGA
  #if FINDERSCOPE_FPGA_BINNING
    #define FINDERSCOPE_MIN_PERIOD_MS   18                              // Minimum finderscope period in ms
  #else
    #define FINDERSCOPE_MIN_PERIOD_MS   15                              // Minimum finderscope period in ms
  #endif
#endif
#define FINDERSCOPE_MAX_PERIOD_MS               80                              // Maximum finderscope period in ms

#define BAND_IMAGE_HEIGHT                               6000                    // Maximum Height of Spectral  Band Images
#define BAND_IMAGE_WIDTH                                1800                    // Width of Spectral Band Images in Light Pixels
#define BAND_DARK_WIDTH                                 16                              // No Dark Pixels at start of row
#define BAND_TOTAL_WIDTH                                (BAND_DARK_WIDTH + BAND_IMAGE_WIDTH)

#define BAND_DARK_TOP                                   0                                       // enumeration constant for Darks at top of Spectral Images
#define BAND_DARK_BOTTOM                                1                                       // enumeration constant for Darks at bottom of Spectral Images
#define BAND_DARK_POSITION                              BAND_DARK_BOTTOM        // Compile time option for position of Spectral Darks
#define BAND_DARK_BOTTOM_TRANSITION             20                                      // When darks at the bottom this indicates how many rows are in transition
                                                                                                                        // from light to dark ||| Ask Alan (9+9+2 safety)

#define BAND_US_PER_EXPOSURE_COUNT      10                                      // scale factor between the asked for interval value and time in microseconds
                                                                                                                //  1 : API Interval of 100 => 100us interval in FPGA
                                                                                                                // 10 : API Interval of 100 => 10*100 = 1000us interval in FPGA
                                                                                                                // Also scales exposures (log signals)
#define BAND_DARK_MAX_RECORDS                   12                              // Maximum number of Average Spectral Dark Pixels rows to encode


#define TELEMETRY_MAX_RECORDS                   12                              // Maximum number of Telemetry Records per Exposure

#define MISSION_LOG_MAX_RECORDS                 3                               // Maximum number of Mission Log Records per Exposure
#define MISSION_LOG_LENGTH                              10000                   // Maximum Length of the Mission Log File

#define MIN_SECTOR_SIZE                                 30                              // Sectors must be at least this long (bytes)
#define MAX_SECTOR_SIZE                                 1024                    // And no longer than this

#define ACK 0x06                // Response to command - acknowledged
#define NAK 0x15                // Response to command - bad checksum
#define DC1 0x11                // Response to command - resource busy
#define DC2 0x12                // Response to command - resource not available
#define CAN 0x18                // Response to command - malformed data payload

/*
 *              Type definitions for the commands and various
 *              enumerated fields.
 */

typedef enum { /*   0 -   5 */ HC_NULL, HC_PING, HC_SET_POWER_STATE, HC_GET_POWER_STATE, HC_GET_TELEMETRY, HC_SET_EXPOSURE_PARAMETERS,
/*   6 -  10 */ HC_GET_EXPOSURE_PARAMETERS, HC_START_EXPOSURE, HC_END_EXPOSURE, HC_GET_EXPOSURE_STATE, HC_GENERATE_TEST_EXPOSURE,
/*  11 -  15 */ HC_SET_IMAGE_SECTOR_SIZE, HC_SET_SPECTRAL_COMPRESSION_PARAMETERS, HC_GET_SPECTRAL_COMPRESSION_PARAMETERS, HC_POST_PROCESS_IMAGE, HC_GET_MISSION_TIME,
/*  16 -  20 */ HC_SET_MISSION_TIME, HC_GET_IMAGE_SECTORS, HC_FLASH_OPERATIONS, HC_LAST_COMMAND,
/* 100 - 105 */ HC_UGA_NULL = 100, HC_UGA_SET_EXPOSURE_PARAMETERS, HC_UGA_GET_EXPOSURE_PARAMETERS, HC_UGA_POST_PROCESS_IMAGE, HC_UGA_GET_EXPOSURE_STATE, HC_UGA_GENERATE_TEST_EXPOSURE,
                                HC_UGA_LAST_COMMAND,
/* 200 - 205 */ HC_DB_NULL = 200, HC_DB_DOWNLOAD, HC_DB_MICRON_GRAB, HC_DB_KLI_GRAB, HC_DB_KLI_STREAM, HC_DB_ECHO,
                                HC_DB_CAMERA_INFO, HC_DB_LAST_COMMAND,
} HAWKEYE_COMMAND;

typedef enum { PT_SHORT, PT_MEDIUM, PT_LONG } PING_TYPE;
typedef enum { PS_POWERED_OFF, PS_LOW_POWER, PS_FULL_POWER, PS_START_SHUTDOWN, PS_SHUTDOWN_COMPLETE } POWER_STATE;
typedef enum { TC_SOFTWARE_VERSION, TC_FPGA_VERSION, TC_INDEX, TC_CCD1_TEMP, TC_CCD2_TEMP,
                           TC_CCD3_TEMP, TC_CCD4_TEMP, TC_FPGA_TEMP, TC_FPGA_VAUX, TC_FPGA_VINT, TC_FPGA_VNVP, TC_CCD_VDD_OC,
                           TC_AD7490_CH01, TC_AD7490_CH02, TC_AD7490_CH03, TC_AD7490_CH04, TC_AD7490_CH05, TC_AD7490_CH06, TC_AD7490_CH07, TC_AD7490_CH08,
                           TC_AD7490_CH09, TC_AD7490_CH10, TC_AD7490_CH11, TC_AD7490_CH12, TC_AD7490_CH13, TC_AD7490_CH14, TC_AD7490_CH15, TC_AD7490_CH16,
                           TC_NO_CHANNELS } TELEMETRY_CHANNELS;
typedef enum { TI_NOT_INTERPRETED, TI_NOMINAL, TI_LOW, TI_HIGH } TELEMETRY_INTERPRETATION;
typedef enum { ES_IDLE, ES_ACTIVE, ES_POST_PROCESSING_IMAGE, ES_POST_PROCESSING_COMPLETE } EXPOSURE_STATE;
typedef enum { TI_GRADIENT, TI_SQUARE_BULLS_EYE } TEST_IMAGE;
typedef enum { DC_UNCOMPRESSED, DC_PACKED, DC_DELTA } DATA_COMPRESSION;
typedef enum { FO_ERASE, FO_SAVE, FO_RESTORE } FLASH_OPERATIONS_COMMAND;
typedef enum { SS_NONE, SS_SOLENOID1, SS_SOLENOID2 } SHUTTER_SOLENOID;          // used in decoding Type 3 Block of Image info
typedef enum { RO_GREEN_FIRST, RO_BLUE_FIRST } READOUT_ORDER;                           //  "        "

typedef enum { UGF_DARK_PIXELS = 1, UGF_USE_SIMULATOR = 2 } DB_MICRON_GRAB_FLAGS;
typedef enum { DLB_MICRON, DLB_KLI4104, DLB_STREAM} DB_DOWNLOAD_BUFFER;
typedef enum { RR_VALID, RR_ACK, RR_NAK, RR_CAN, RR_DC1, RR_DC2, RR_BAD_CHECKSUM, RR_BAD_START, RR_BAD_COMMAND, RR_BAD_LENGTH, RR_UNKNOWN } VALIDATE_RESULT;
typedef enum { DBE_CANDC, DBE_QSPI } DB_ECHO_CHANNEL;
/*
 *              Command and Response Struct Definitions
 *
 *              Each Hawkeye Command has a struct defined for passing
 *              parameters to the command and for returning the response.
 *
 *              These struct are byte packed so there are no extra bytes
 *              between the various variables in the structs.
 *
 */
#pragma pack(push)
#pragma pack(1)

/*
 *      Ping Command
 */
typedef struct {
        uint16_t pingType;
} PingParams;

typedef struct {
        uint8_t ack;
} PingShortResponse;

typedef struct {
        uint8_t ping[10];
} PingMediumResponse;

typedef struct {
        uint8_t ping[256];
} PingLongResponse;

/*
 *      Set / Get Power State
 */
typedef struct {
        uint16_t powerState;
} SetPowerStateParams;

typedef struct {
        uint16_t powerStateIn;
        uint16_t powerStateOut;
} GetPowerStateResponse;

/*
 *      Get Telemetry
 */
typedef struct {
        uint16_t channelValue;
        uint16_t channelInterp;
} TelemetryPair;

typedef struct {
        TelemetryPair channel[TC_NO_CHANNELS];
} GetTelemetryResponse;

/*
 *      Set / Get Mission Time
 */

typedef struct {
        uint8_t epoch_time[5];
} SetMissionTimeParams;

typedef SetMissionTimeParams GetMissionTimeResponse;

/*
 *      Set / Get Image Parameters
 */

typedef struct {
        uint16_t channelBitfield;
        uint16_t ccd1Exposure;
        uint16_t ccd2Exposure;
        uint16_t ccd3Exposure;
        uint16_t ccd4Exposure;
        uint16_t height;
        uint16_t interval;
        uint16_t oversampling;
        uint16_t darkHeight;
        uint16_t finderscopeExposure;
        uint16_t noFinderscopeImages;
} SetExposureParametersParams;

typedef SetExposureParametersParams GetExposureParametersResponse;

/*
 *      Start Exposure
 */
typedef struct {
        uint32_t exposureID;
} StartExposureParams;

/*
 *      End Exposure
 */
typedef struct {
        uint8_t epoch_time[5];
} EndExposureParams;

/*
 *      Get Exposure State
 */
typedef struct {
        uint32_t exposureID;
        uint16_t exposureState;
        uint16_t rowsAcquired;
        uint16_t errorCode;
} GetExposureStateResponse;

/*
 *      Generate Test Image
 */
typedef struct {
        uint16_t testImage;
        uint16_t height;
} GenerateTestExposureParams;

/*
 *      Post Process Image
 */
typedef struct {
        uint16_t cropLeft;
        uint16_t cropTop;
        uint16_t cropWidth;
        uint16_t cropHeight;
        uint16_t darkSubtraction;
        uint16_t imageBinning;
        uint16_t finderscopeBinning;
        uint16_t compression;
        uint32_t imageID;
} PostProcessImageParams;

/*
 *      Set Image Sector Size
 */
typedef struct {
        uint16_t sectorSize;
} SetImageSectorSizeParams;

typedef struct {
        uint32_t imageSize;
        uint32_t numberSectors;
} SetImageSectorSizeResponse;

/*
 *      Set Spectral Compression Parameters
 */
typedef struct {
        uint16_t slope1;
        uint16_t slope2;
        uint16_t knee;
} BandCompression;

typedef struct {
        uint16_t ccd1Gain;
        uint16_t ccd2Gain;
        uint16_t ccd3Gain;
        uint16_t ccd4Gain;
        BandCompression bandCompression[8];
} SetSpectralCompressionParametersParams;

typedef SetSpectralCompressionParametersParams GetSpectralCompressionParametersResponse;

/*
 *      Get Image Sectors
 */
typedef struct {
        uint32_t startingSector;
        uint32_t noSectors;
} GetImageSectorsParams;

typedef struct {
        uint32_t sectorNumber;
        //
        // The following field is repeated N times
        // where N = (Payload Length - 4) / 2
        uint16_t pixel;
} GetImageSectorsResponse;

/*
 *      Flash Operations
 */
typedef struct {
        uint16_t command;
} FlashOperationsParams;

/*
 *              University of Georgia additions to the
 *              Hawkeye Commands
 */

/*
 *      UGA Set / Get Exposure Parameters
 */

typedef struct {
        uint16_t finderscopeExposure;
        uint16_t noFinderscopeImages;
        uint16_t finderscopeVerticalBinning;
        uint16_t framePeriodMS;
} UGASetExposureParametersParams;

typedef UGASetExposureParametersParams UGAGetExposureParametersResponse;

/*
 *      UGA Post Process Image
 */
typedef struct {
        uint16_t finderscopeBinning;
        uint16_t compression;
        uint32_t imageID;
} UGAPostProcessImageParams;

/*
 *      UGA Get Exposure State
 */
typedef struct {
        uint32_t exposureID;
        uint16_t exposureState;
        uint16_t imagesAcquired;
        uint16_t errorCode;
} UGAGetExposureStateResponse;

/*
 *      UGA Generate Test Image
 */
typedef struct {
        uint16_t testImage;
        uint16_t noImages;
} UGAGenerateTestExposureParams;

/*
 *              Debug Commands
 *
 *              The following structs ar for the Debug Commands
 *              which are used for testing the system during
 *              development.
 *
 *              All Debug structs and enums start with DB
 */

typedef struct {
        uint16_t imageBuffer;
        uint32_t startPixel;
        uint32_t pixelLength;
} DBDownloadParams;

typedef struct {
        uint32_t exposure;
        uint16_t flags;
        uint16_t height;
        uint16_t width;
} DBMicronGrabParams;

typedef struct {
        uint16_t error;
        uint16_t height;
        uint16_t width;
        uint16_t min;
        uint16_t max;
        uint32_t ave100X;
} DBMicronGrabResponse;

typedef struct {
        uint32_t exposure;
        uint32_t interval;
        uint16_t oversampling;
        uint16_t channel;
        uint16_t height;
        uint16_t width;
        uint16_t left;
} DBKLIGrabParams;

typedef DBMicronGrabResponse DBKLIGrabResponse;

typedef struct {
        uint16_t channel;
} DBKLIStreamParams;

typedef DBMicronGrabResponse DBKLIStreamResponse;

typedef struct {
        uint16_t channel;                 // enum DB_ECHO_CHANNEL
        uint8_t data[200];
} DBEchoParams;

typedef DBEchoParams DBEchoResponse;

typedef struct {
        uint16_t firmwareVersion;         // BCD xx.xx
        uint16_t fpgaVersion;                     //  "    "
        uint16_t finderscopeHeight;       // vertical pixels in micron sensor
        uint16_t finderscopeWidth;        // horizontal pixels in micron sensor
        uint16_t finderscopeVBinning;     // vertical binning in micron sensor
} DBCameraInfoResponse;

#pragma pack(pop)

#ifdef __cplusplus
        extern "C" uint16_t Checksum(uint8_t *src, int len);
        extern "C" uint16_t Swap2(uint16_t us);
        extern "C" uint32_t Swap4(uint32_t ul);
        extern "C" void Swap2Copy(uint16_t *dest, uint16_t *scr, int len);
        extern "C" int HawkeyeBuildCommand(uint8_t *dest, int command, void *pParams, uint16_t paramsLen);
        extern "C" VALIDATE_RESULT HawkeyeValidateResponse(uint8_t *src, int command, void* pResponse, uint16_t responseLen);
#else
        extern uint16_t Checksum(uint8_t *src, int len);
        extern uint16_t Swap2(uint16_t us);
        extern uint32_t Swap4(uint32_t ul);
        extern void Swap2Copy(uint16_t *dest, uint16_t *scr, int len);
        extern int HawkeyeBuildCommand(uint8_t *dest, int command, void *pParams, uint16_t paramsLen);
        extern VALIDATE_RESULT HawkeyeValidateResponse(uint8_t *src, int command, void* pResponse, uint16_t responseLen);
#endif

#endif // _HAWKEYE_

