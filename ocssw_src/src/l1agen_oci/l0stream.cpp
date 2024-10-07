/**
 *
 * @authors Jakob C. Lindo (SSAI)
 * @date Apr 2024
 * @brief The implementation of a class to simulate a single stream object over multiple OCI L0 files
 *
 */

#include "l0stream.hpp"
#include <cstring>

#define PACE_OCI_HDR_LEN 64
#define EXP_CONTENT_TYPE "cFE1"
#define EXP_CONTENT_SUBTYPE " SCI"
#define EXP_SPACECRAFT_ID "PACE"
#define EXP_PROCESSOR_ID " OCI"
#define CURR_L0_FILE this->l0Files.at(idxCurrStream)

using namespace std;

L0Stream::L0Stream(vector<string> &l0FileNames) {
    for (auto filename : l0FileNames) {
        idxCurrStream = 0;

        L0File *l0File = new L0File;

        l0File->fileStream = new fstream(filename, fstream::in | fstream::binary);
        l0File->fileSize = getFileSize(*l0File->fileStream);
        l0File->hasHeader = hasHeader(*l0File->fileStream);
        l0File->bytesLeft = l0File->fileSize;

        l0Files.push_back(l0File);
    }
}

L0Stream::~L0Stream() {
    for (auto l0File : l0Files) {
        l0File->fileStream->close();
        delete(l0File->fileStream);
        delete(l0File);
    }
}

void L0Stream::exit_unhappy(string message) {
    cerr << "-E- " << message << endl;
    exit(EXIT_FAILURE);
}

void L0Stream::strip_oci_header(fstream &fileStream) {
    /**
     * Typical OCI header
     * 00000000  63 46 45 31 20 53 43 49  00 00 00 40 50 41 43 45  |cFE1 SCI...@PACE|
     * 00000010  00 00 00 03 20 4f 43 49  7c 8c a0 d3 52 07 30 00  |.... OCI|...R.0.|
     * 00000020  4f 43 49 30 30 30 30 32  31 32 34 39 2e 6f 63 69  |OCI000021249.oci|
     * 00000030  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  |................|
     */

    char hdr[PACE_OCI_HDR_LEN];  // Full .oci file header
    fileStream.read(hdr, PACE_OCI_HDR_LEN);

    if (string(hdr + 0x00, 4) != EXP_CONTENT_TYPE || string(hdr + 0x04, 4) != EXP_CONTENT_SUBTYPE ||
        string(hdr + 0x0c, 4) != EXP_SPACECRAFT_ID || string(hdr + 0x14, 4) != EXP_PROCESSOR_ID)
        fileStream.seekg(0, fileStream.beg);
}

bool L0Stream::hasHeader(fstream &fileStream) {
    L0Stream::strip_oci_header(fileStream);
    return fileStream.tellg() == 0 ? false : true;
}

size_t L0Stream::getFileSize(std::fstream &fileStream) {
    fileStream.seekg(0, fileStream.end);
    size_t fileSize = fileStream.tellg();

    fileStream.seekg(0, fileStream.beg);
    return fileSize;
}

size_t L0Stream::getBytesLeft() {
    return CURR_L0_FILE->fileSize - CURR_L0_FILE->currPos;
}

void L0Stream::read(char *buffer, size_t numBytes) {
    fstream *currStream = CURR_L0_FILE->fileStream;
    size_t bytesLeft = CURR_L0_FILE->bytesLeft;

    if (bytesLeft <= numBytes) {
        int numBytesAfterCarryover = numBytes - bytesLeft;

        char *carryover = (char *)calloc(bytesLeft, sizeof(char));
        currStream->read(carryover, bytesLeft);

        memcpy(buffer, carryover, bytesLeft);
        free(carryover);

        idxCurrStream++;  // Switch streams

        if (idxCurrStream >= l0Files.size())  // No more files to read
            return;

        currStream->read(buffer + bytesLeft, numBytesAfterCarryover);  // Read remaining bytes
        CURR_L0_FILE->currPos = CURR_L0_FILE->fileStream->tellg();
        CURR_L0_FILE->bytesLeft = getBytesLeft();
        return;
    }

    CURR_L0_FILE->fileStream->read(buffer, numBytes);
    CURR_L0_FILE->currPos = CURR_L0_FILE->fileStream->tellg();
    CURR_L0_FILE->bytesLeft = getBytesLeft();

    if (CURR_L0_FILE->currPos == CURR_L0_FILE->fileSize)
        idxCurrStream++;
}

void L0Stream::seekg(std::streamoff num, std::ios_base::seekdir dir) {
    int currPos = CURR_L0_FILE->currPos;
    int abs_num = abs(num);  // Absolute value of num

    if (num < 0 && currPos < abs_num) {  // Seeking toward beginning, into the previous file

        if (idxCurrStream <= 0)
            exit_unhappy("Tried to seek into a file that doesn't exist");

        CURR_L0_FILE->fileStream->seekg(0);  // Reset pos of this stream to beginning
        idxCurrStream--;
        CURR_L0_FILE->fileStream->seekg(abs_num - currPos, CURR_L0_FILE->fileStream->end);
    } else {
        CURR_L0_FILE->fileStream->seekg(num, dir);
    }
}

bool L0Stream::fail() {
    if (idxCurrStream < l0Files.size())
        return CURR_L0_FILE->fileStream->fail();
    else
        return true;
}

bool L0Stream::eof() {
    if (idxCurrStream < l0Files.size())
        return CURR_L0_FILE->fileStream->eof();
    else
        return true;
}

int L0Stream::tellg() {
    if (idxCurrStream < l0Files.size())
        return CURR_L0_FILE->fileStream->tellg();
    else
        return -1;
}
