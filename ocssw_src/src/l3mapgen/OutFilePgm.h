/**
 * @file OutFilePgm.h
 * @brief Defines the OutFilePgm class for handling PGM file output
 */

#ifndef OUTFILEPGM_H
#define OUTFILEPGM_H


#include "OutFile.h"

/**
 * @class OutFilePgm
 * @brief A class for creating/managing PGM file output.
 * 
 * OutFilePgm inherits from OutFile and provides functionality specific to PGM (portable graymap)
 * file format.
 */
class OutFilePgm : public OutFile {
   protected:
    FILE* outfp;
    uint8_t* fileData;

   public:
    /**
     * @brief Default constructor.
     */
    OutFilePgm();

    /**
     * @brief Virtual destructor.
     */
    virtual ~OutFilePgm();

    /**
     * @brief Set the size of the PGM output file.
     * @param width  Width of the output file in pixels.
     * @param height Height of the output file in pixels.
     */
    virtual void setSize(int32_t width, int32_t height) override;

    /**
     * @brief Open the PGM output file and write the header.
     * @return true if the file is successfully opened, else false.
     */
    virtual bool open() override;

    /**
     * @brief Close the PGM output file.
     * @return True if the file is successfully closed.
     */
    virtual bool close() override;

    /**
     * @brief Write a line of data to the PGM output file.
     */
    virtual void writeLine() override;
};

#endif // OUTFILEPGM_H