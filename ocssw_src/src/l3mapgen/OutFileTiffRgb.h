#ifndef OUTFILETIFFRGB_H
#define OUTFILETIFFRGB_H

#include "OutFileTiff.h"
#include "OutFile.h"

class OutFileTiffRgb : public OutFileTiff {
    uint8_t* fileData = nullptr;

   public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~OutFileTiffRgb();

    /**
     * @brief Set the size of the TIFF RGB output file.
     * @param width  Width of the output file.
     * @param height Height of the output file.
     */
    virtual void setSize(int32_t width, int32_t height);

    /**
     * @brief Write a line of data to the TIFF RGB output file.
     */
    virtual void writeLine();

    /**
     * @brief Set the color configuration for the TIFF RGB output file.
     */
    virtual void setTiffColor();

    /**
     * @brief Set a pixel value in the TIFF RGB output file.
     * @param x     Pixel x-coordinate.
     * @param val   Pixel value.
     * @param prod  Product index.
     */
    virtual void setPixel(int32_t x, double val, int32_t prod = 0);

    /**
     * @brief Set an RGB pixel value in the TIFF RGB output file.
     * @param x     Pixel x-coordinate.
     * @param red   Red value.
     * @param green Green value.
     * @param blue  Blue value.
     */
    virtual void setPixelRGB(int32_t x, float red, float green, float blue);
    virtual void landPixel(int32_t x);
    virtual void fillPixel(int32_t x);
    virtual void missingPixel(int32_t x);
};

#endif // OUTFILETIFFRGB_H