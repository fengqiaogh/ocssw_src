/**
 * @file OutFilePpmRgb.h
 * @brief Defines the OutFilePpmRgb class for handling RGB PPM file output
 */

#ifndef OUTFILEPPMRGB_H
#define OUTFILEPPMRGB_H

#include "OutFile.h"

/**
 * @class OutFilePpmRgb
 * @brief A class for creating a managing RGB PPM file output
 * 
 * OutFilePpmRgb inherits from OutFile and provides functionality
 * specific to RGB PPM (portable pixmap) file format
 */
class OutFilePpmRgb : public OutFile {
    FILE* outfp; ///< File pointer for output operations
    uint8_t* fileData; ///< Buffer for file data

   public:
    /**
     * @brief Default constructor.
     */
    OutFilePpmRgb();

    /**
     * @brief Virtual destructor.
     */
    virtual ~OutFilePpmRgb();

    /**
     * @brief Set the size of the PPM RGB output file.
     * @param width  Width of the output file in pixels.
     * @param height Height of the output file in pixels.
     */
    virtual void setSize(int32_t width, int32_t height) override;

    /**
     * @brief Open the PPM RGB output file and write the header.
     * @return true if the file is successfully opened, else false.
     */
    virtual bool open() override;

    /**
     * @brief Close the PPM RGB output file.
     * @return true if the file is successfully closed, false otherwise.
     */
    virtual bool close() override;

    /**
     * @brief Set a pixel value in the PPM RGB output file.
     * @param x     Pixel x-coordinate.
     * @param val   Pixel value.
     * @param prod  Product index.
     */
    virtual void setPixel(int32_t x, double val, int32_t prod = 0) override;

    /**
     * @brief Set the RGB values of a pixel in the PPM RGB output file.
     * @param x      Pixel x-coordinate.
     * @param red    Red component value.
     * @param green  Green component value.
     * @param blue   Blue component value.
     */
    virtual void setPixelRGB(int32_t x, float red, float green, float blue) override;

    /**
     * @brief Set a pixel to land value in the PPM RGB output file.
     * @param x Pixel x-coordinate.
     */
    virtual void landPixel(int32_t x) override;

    /**
     * @brief Set pixel as fill pixel in the PPM RGB output file.
     * @param x Pixel x-coordinate.
     */
    virtual void fillPixel(int32_t x) override;

    /**
     * @brief Mark a pixel as missing in the PPM RGB output file (if within range).
     * @param x Pixel x-coordinate.
     */
    virtual void missingPixel(int32_t x) override;

    /**
     * @brief Write a line of data to the PPM RGB output file.
     */
    virtual void writeLine() override;
};


#endif // OUTFILEPPMRGB_H