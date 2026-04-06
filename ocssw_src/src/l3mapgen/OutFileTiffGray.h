#ifndef OUTFILETIFFGRAY_H
#define OUTFILETIFFGRAY_H

#include "OutFileTiff.h"

class OutFileTiffGray : public OutFileTiff {
    float* fileData = nullptr;

   public:
    /**
     * @brief Destructor for OutFileTiffGray class.
     */
    virtual ~OutFileTiffGray();

    /**
     * @brief Set the size of the TIFF Gray output file.
     * @param width  Width of the output file.
     * @param height Height of the output file.
     */
    virtual void setSize(int32_t width, int32_t height);

    /**
     * @brief Write a line of data to the TIFF Gray output file.
     */
    virtual void writeLine();

    /**
     * @brief Set color configuration for the TIFF Gray output file.
     */
    virtual void setTiffColor();
};

#endif // OUTFILETIFFGRAY_H