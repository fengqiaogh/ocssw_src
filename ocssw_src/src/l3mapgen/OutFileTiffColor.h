#ifndef OUTFILETIFFCOLOR_H
#define OUTFILETIFFCOLOR_H

#include "OutFileTiff.h"
#include <memory>

class OutFileTiffColor : public OutFileTiff {
    uint8_t* fileData;

   public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~OutFileTiffColor();

    /**
     * @brief Set the size of the TIFF Color output file.
     * @param width  Width of the output file.
     * @param height Height of the output file.
     */
    virtual void setSize(int32_t width, int32_t height);

    /**
     * @brief Write a line of data to the TIFF Color output file.
     */
    virtual void writeLine();

    /**
     * @brief Set the color configuration for the TIFF output file.
     */
    virtual void setTiffColor();
};

#endif // OUTFILETIFFCOLOR_H
