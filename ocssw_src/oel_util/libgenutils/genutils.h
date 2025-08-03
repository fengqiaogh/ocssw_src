#ifndef _GENUTILS_H
#define _GENUTILS_H

#include <stdint.h>
#include <stdio.h>

#include "passthebuck.h"

// circumference of the earth at the equator in meters
#define EARTH_CIRCUMFERENCE 40075016.6856
#define OEL_PI 3.14159265358979323846264338327950288 // Copy of the current C++ definition of Pi
#define OEL_RADEG (180.0 / OEL_PI) // Radians to degrees using C++ definition of Pi
#define OEL_DEGRAD (OEL_PI / 180.0) // Degrees to radians using C++ definition of Pi

#ifdef __cplusplus
extern "C" {
#endif

/**\file 
   This is a bunch of useful general utilities.
 */

/** Define a few useful BAD values
 */
#define BAD_FLT  -32767.0
#define BAD_INT    -32767
#define BAD_UINT    65535
#define BAD_BYTE     -128
#define BAD_UBYTE     255

/** Do we want to print out extra info to stdout */
extern int want_verbose;

/** Is this machine little endian
    \return 0 for BIG_ENDIAN machines, 1 for LITTLE_ENDIAN mahcines 
 */
int endianess(void);

void parse_file_name(const char *inpath, char *outpath);

/** Swap bytes in place 
    \param in pointer to memory to byte swap
    \param nbyte size of object to reverse
    \param ntime number of objects to swap
    \return 0 always
 */
int swapc_bytes(char *in, int nbyte, int ntime);

/** Swap bytes from in to out.  in and out should not be overlapping memory
    \param in pointer to source memory to byte swap
    \param out pointer to destination memory
    \param nbyte size of object to reverse
    \param ntime number of objects to swap
    \return 0 always
 */
int swapc_bytes2(const char *in, char *out, int nbyte, int ntime);

/** read a binary file and swap bytes if necessary
    \param little_endian set to 0 for big endian files, 1 for little endian files
    \param ptr memory to read into
    \param size size of object to read in bytes
    \param nmemb number of objects to read
    \param stream file to read from
    \return number of objects read
 */
size_t fread_swap(int little_endian, void *ptr, size_t size, size_t nmemb, FILE *stream);

/** write a binary file and swap bytes if necessary
    \param little_endian set to 0 for big endian files, 1 for little endian files
    \param ptr memory to write
    \param size size of object to write in bytes
    \param nmemb number of objects to write
    \param stream file to write to
    \return number of objects written
 */
size_t fwrite_swap(int little_endian, const void *ptr, size_t size, size_t nmemb, FILE *stream);

void spline(float [], float [], int, float, float, float []);
void splint(float [], float [], float [], int, float, float *);
void lspline(float xin [], float yin [], int32_t nin,
        float xout[], float yout[], int32_t nout);
float linterp(float xin [], float yin [], int32_t nin, float xout);
float bioBandShift(float xin [], float yin [], float chl, int32_t nin, float xout);

int32_t filesize_(const char *filename);
int32_t filesize(const char *filename);

int getlut_file(char *lut_file, short *rlut, short *glut, short * blut);

char *lowcase(char *instr);
char *upcase(char *instr);

int isValidInt(const char* str);

void* allocateMemory(size_t numBytes, const char* name);

void trimBlanks(char* str);
char* trimBlanksDup(const char* str);

int getFileFormatIndex(const char* str);
const char* getFileFormatName(const char* str);
const char* getFileFormatExtension(const char* str);

char* replace_ocroots(const char* inStr);

/**
 * @brief Pixel resolution string to meters
 * 
 * @param resolutionStr resolution string
 * @return double in meters
 */
double str2resolution(char const * resolutionStr);

/**
 * @brief Pixel resolution meters to string
 * 
 * @param resolution in meters
 * @return const char* resolution as a string with units
 */
const char* resolution2str(double resolution);

/**
 * @brief Pixel resolution meters to degrees
 * 
 * @param resolution in meters
 * @return double in degrees
 */
double resolution2degrees(double resolution);

/**
 * @brief Pixel resolution degrees to meters
 * 
 * @param degrees resolution in degrees
 * @return double in meters
 */
double degrees2resolution(double degrees);

/**
 * @brief convert l2bin style resolve string to number line in L3bin file
 * 
 * @param resolve l2bin style string
 * @return int32_t number of lines in L3bin file
 */
int32_t resolve2binRows(const char *resolve);

/**
 * @brief convert l2bin style resolve string to pixel resolution in meters
 * 
 * @param resolve l2bin style string
 * @return double resolution in meters
 */
double resolve2resolution(const char *resolve);

/**
 * @brief chek if the file is url
 * @param file input file
 * @return 1, if it's an url. 0 if it's not
 */
int check_url(const char *file);

#ifdef __cplusplus
}

// C++ functions
#include <string>
#include <vector>
void replaceOCroots(std::string& str);

/**
 * @brief Pixel resolution string to meters
 * 
 * @param resolutionStr pixel resolution
 * @return double resolution in meters
 */
double string2resolution(std::string resolutionStr);

/**
 * @brief Pixel resolution from meters to resolution string
 * 
 * @param resolution in meters
 * @return std::string resolution string
 */
std::string resolution2string(double resolution);

/**
 * @brief Text file reader 
 * 
 * @param fileName 
 * @return vector<string>
 */
std::vector<std::string> readFileList(std::string fileName);

/**
 * @brief convert an l2bin resolve string in L3bin file num rows and resolution
 * 
 * @param resolve l2bin stype resolve string
 * @param nrows number of rows for L3bin file
 * @param resolution resolution in meters
 */
void resolve2binRows(std::string resolve, int32_t &nrows, double &resolution);


#endif // C++

#endif
