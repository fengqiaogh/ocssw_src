
/**
 * @file OutFilePngRgb.h
 * @brief Desgines the OutFilePngRgb class for handling RGB PNG file output.
 */
#ifndef OUTFILEPNGRGB_H
#define OUTFILEPNGRGB_H

#include "OutFile.h"

/**
 * @class OutFilePngRgb
 * @brief A class for creating/managing RGB PNG file output.
 * 
 * OutFilePngRgb inherits from OutFile and provides functionality
 * specific to RGB PNG file format.
 */
class OutFilePngRgb : public OutFile {
    FILE* outfp; ///< File pointer for output
    uint8_t* fileData; ///< Buffer for file data
    png_structp pngPtr; ///< PNG structure pointer
    png_infop infoPtr; ///< PNG info structure pointer
    uint numText;

   public:
    /**
     * @brief Default constructor.
     */
    OutFilePngRgb();

    /**
     * @brief Virtual destructor.
     */
    virtual ~OutFilePngRgb();


    /**
     * @brief Set the size of the PNG RGB output file.
     * @param width  Width of the output file.
     * @param height Height of the output file.
     */
    virtual void setSize(int32_t width, int32_t height);

    /**
     * @brief Open the PNG RGB output file and write the header.
     * @return True if the file is successfully opened, else False.
     */
    virtual bool open();

    /**
     * @brief Close the PNG RGB output file.
     * @return True if the file is successfully closed, else false.
     */
    virtual bool close();

    /**
     * @brief Set a pixel value in the PNG RGB output file.
     * @param x     Pixel x-coordinate.
     * @param val   Pixel value to set.
     * @param prod  Product index (default is 0).
     */
    virtual void setPixel(int32_t x, double val, int32_t prod = 0);

    /**
     * @brief Set the RGB values of a pixel in the PNG RGB output file.
     * @param x      Pixel x-coordinate.
     * @param red    Red component value of the pixel.
     * @param green  Green component value.
     * @param blue   Blue component value.
     */
    virtual void setPixelRGB(int32_t x, float red, float green, float blue);

    /**
     * @brief Set a pixel as a land pixel.
     * @param x Pixel x-coordinate.
     */
    virtual void landPixel(int32_t x);

    /**
     * @brief Fill a pixel in the PNG RGB output file.
     * @param x Pixel x-coordinate.
     */
    virtual void fillPixel(int32_t x);

    /**
     * @brief Mark a pixel as a missing pixel in the PNG RGB output file.
     * @param x Pixel x-coordinate.
     */
    virtual void missingPixel(int32_t x);

    /**
     * @brief Write a line of data to the PNG RGB output file.
     */
    virtual void writeLine();
};

#endif // OUTFILEPNGRGB