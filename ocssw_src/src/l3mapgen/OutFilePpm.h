/**
 * @file OutFilePpm.h
 * @brief Defines the OutFilePpm class for handling PPM file output
 */

#ifndef OUTFILEPPM_H
#define OUTFILEPPM_H

#include "OutFilePgm.h"

/**
 * @class OutFilePpm
 * @brief A class for creating/managing PPM file output.
 * 
 * OutFilePpm inherits from OutFilePgm and provides functionality specific to PPM (portable pixmap)
 * file format.
 */
class OutFilePpm : public OutFilePgm {
   public:
    /**
     * @brief Set the size of the PPM output file.
     * @param width  Width of the output file.
     * @param height Height of the output file.
     */
    virtual void setSize(int32_t width, int32_t height);

    /**
     * @brief Open the PPM output file and write the header.
     * @return True if the file is successfully opened, else false.
     */
    virtual bool open();

    /**
     * @brief Write a line of data to the PPM output file.
     */
    virtual void writeLine();
};

#endif // OUTFILEPPM_H