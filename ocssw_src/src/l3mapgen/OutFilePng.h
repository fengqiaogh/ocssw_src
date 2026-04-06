/**
 * @file OutFilePng.h
 * @brief Defines the OutFilePng class for handling PNG file output
 */

#ifndef OUTFILEPNG_H
#define OUTFILEPNG_H

#include "OutFile.h"

/**
 * @class OutFilePng
 * @brief A class for creating/managing PNG file output.
 * 
 * OutFilePng inherits from OutFile and provides functionality specific 
 * to PNG (portable network graphics) file format
 */
class OutFilePng : public OutFile {
    FILE* outfp; ///< File pointer for output operations
    uint8_t* fileData; ///< Buffer for file data
    bool isColor; ///< Flag indicating if the output is color
    png_structp pngPtr; ///< PNG structure pointer
    png_infop infoPtr; ///< PNG ingo structure pointer
    uint numText; ///< Number of text chunks

   public:
    /**
     * @brief Constructor for OutFilePng class.
     * @param color Whether the output PNG file is in color.
     */
    OutFilePng(bool color);

    /**
     * @brief Virtual destructor.
     */
    virtual ~OutFilePng();

    /**
     * @brief Set the size of the PNG output file.
     * @param width  Width of image in pixels.
     * @param height Height of the image in pixels.
     */
    virtual void setSize(int32_t width, int32_t height);

    /**
     * @brief Open the PNG output file and write the header.
     * @return True if the file is successfully opened, else false.
     */
    virtual bool open();

    /**
     * @brief Close the PNG output file.
     * @return True if the file is successfully closed, else false.
     */
    virtual bool close();

    /**
     * @brief Write a line of data to the PNG output file.
     */
    virtual void writeLine();
};

#endif // OUTFILEPNG_H