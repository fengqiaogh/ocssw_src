/**
 * @file OutFileTiff.h
 * @brief Defines OutFileTiff class for handling TIFF file output
 */

#ifndef OUTFILETIFF_H
#define OUTFILETIFF_H

#include "OutFile.h"
#include <xtiffio.h>
#include <geotiffio.h>
#include <geotiff.h>
#include <geo_tiffp.h>

class OutFileTiff : public OutFile {
   protected:
    TIFF* tiff = nullptr;
    GTIF* gtif = nullptr;

   public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~OutFileTiff();

    /**
     * @brief Open the TIFF output file and write the header.
     * @return True if the file is successfully opened, else false.
     */
    virtual bool open();

    /**
     * @brief Close the TIFF output file.
     * @return true if the file is successfully closed, else false
     */
    virtual bool close();


    virtual void setTiffColor() = 0;
};

#endif // OUTFILETIFF_H