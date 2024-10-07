/* loghdr.h - contains nebula file header structure definition */

#include <stdint.h>

struct LOGHDR {
    unsigned short int headerVersion;
    unsigned short int timeMode;
    uint32_t blockSize;
    uint32_t dataSize;
    unsigned short int dataOffset;
    unsigned short int headerMode;
    uint32_t dataId;
    uint32_t timeStamp;
    unsigned short int timeXStamp;
    unsigned short int statusWord;
    unsigned short int dataMode;
    unsigned short int groupNumber;
    uint32_t localTimeStamp;
    unsigned short int localTimeXStamp;
    unsigned short int reserved;
};


typedef unsigned char BYTE;
typedef unsigned short int WORD;
typedef uint32_t DWORD;
typedef float IEEE_FLOAT;
typedef double IEEE_DOUBLE;

#define LOBYTE(w)        ((BYTE) (w & 0xFF))
#define HIBYTE(w)        ((BYTE) (((WORD)(w) >> 8) & 0xFF))
#define LOWORD(l)        ((WORD)(DWORD) (l & 0xFFFF))
#define HIWORD(l)        ((WORD) (((DWORD)(l) >> 16) & 0xFFFF))

#define SWAP_WORD(w)     ((WORD) ((WORD)HIBYTE(w)|(((WORD)LOBYTE(w))<<8)))
#define SWAP_LONG(l)     ((DWORD) ((DWORD)SWAP_WORD(HIWORD(l))) \
                                        | ((DWORD)(SWAP_WORD(LOWORD(l)))<<16))

#define HDR_WORD(w)      ((swap)?SWAP_WORD(w):w)
#define HDR_LONG(l)      ((swap)?SWAP_LONG(l):l)


