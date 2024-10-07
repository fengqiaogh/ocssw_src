/**
 *
 * @authors Jakob C. Lindo (SSAI)
 * @date Apr 2024
 * @brief A class to simulate a single stream object over multiple OCI L0 files
 *
 */

#ifndef _L0STREAM_OCI_
#define _L0STREAM_OCI_

#include <genutils.h>
#include <iostream>
#include <fstream>

class L0Stream {
   public:
    typedef struct L0File {
        std::fstream *fileStream;  // Open L0 stream
        bool hasHeader;            // `true` if `fileStream` has an OCI header, `false` otherwise
        size_t fileSize;           // Number of bytes in `fileStream`
        size_t currPos;            // Current position of read head in `fileStream`
        size_t bytesLeft;          // Amount of bytes still in this stream
    } L0File;

    L0Stream(std::vector<std::string> &l0FileNames);
    ~L0Stream();

    /**
     * @brief Read from the current stream, switching streams as necessary.
     * @param buffer The data structure into which bytes will be read.
     * @param numBytes The size of the block to be read.
     */
    void read(char *buffer, size_t num);

    /**
     * @brief Move the read/write head a number of bytes in a particular direction, switching streams as
     * necessary.
     * @param num The number of bytes to seek. May be negative, indicating a backwards seek.
     * @param dir The direction to seek.
     */
    void seekg(std::streamoff num, std::ios_base::seekdir dir);

    /**
     * @brief Indicate if the fstream::badbit or fstream::failbit is set
     * @return `true` if fstream::badbit or fstream::failbit is set, `false` otherwise
     */
    bool fail();

    /**
     * @brief Indicate if fstream::eof is true
     * @return fstream::eof()
     */
    bool eof();

    /**
     * @brief Obtain the position of the read/write head from the beginning of the current filestream.
     * @return -1 if `fail()` or `eof()` is true, otherwise the position of the read/write head
     */
    int tellg();

    static void strip_oci_header(std::fstream &fileStream);

    /**
     * @brief Strips the PACE OCI header off of `fileStream` if it exists, and returns whether it does so or
     * not
     * @param fileStream An open fstream
     * @return `true` if a header exists in the file, `false` if not.
     */
    bool hasHeader(std::fstream &fileStream);

    /**
     * @brief Obtain the size in bytes of `fileStream`. This method does **not** guarantee that the read/write
     * head returns to its position before being called.
     * @param fileStream An open fstream.
     * @return The number of bytes in `fileStream`
     */
    size_t getFileSize(std::fstream &fileStream);

    /**
     * @brief Obtain the number of bytes between the read/write head and the end of the current stream.
     * @param fileStream The stream object
     * @return The number of bytes between the read/write head and the end of the current stream.
     */
    size_t getBytesLeft();

    /**
     * @brief For use when the program has reached a bad state
     * @param message A helpful message to print to cerr indicating what went wrong
     */
    void exit_unhappy(std::string message);

   private:
    std::vector<L0File *> l0Files;
    size_t idxCurrStream;
};

#endif  // _L0STREAM_OCI_