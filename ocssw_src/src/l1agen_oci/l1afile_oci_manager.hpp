#ifndef _L1AFILE_OCI_MANAGER_
#define _L1AFILE_OCI_MANAGER_

#include <map>
#include "l1afile_oci.hpp"
#include "dimension_shape.hpp"

// Manages writing to multiple L1A files for OCI. When there's a datatype change in-between OCI science data packets, 
// create a new l1a file and write to it. When it ends, resume writing into the science l1a file.
//
// It wraps around L1aFile class and manages multiple versions of it if there's a datatype change.

class L1aFileManager {
    public:
        L1aFileManager();
        ~L1aFileManager();
        
        // reference multiple L1a files based on data type id
        std::map<short, L1aFile> fileCabinet; 

        // tracks the UNLIMITED dimensions for each l1a file so when file changes, you know where
        // to properly append the data
        std::map<short, DimensionShape> fileDimShapes;


        /**
         * @brief create a l1a file for dataType
         * @param dataType 
         * @return 
         */
        L1aFile* createL1aOutputFile(short dataType);
        
        /**
         * @brief checks if an L1A file has been made for the datatype
         * @param dataType 
         * @return 
         */
        bool contains(short dataType);

        L1aFile* getL1aFile(short dataType);

        DimensionShape* getCurrL1aFileDimShape(short dataType);

        bool incrementNumScansShape(short dataType, size_t incrementBy);

        bool incrementNumMceScanShape(short dataType, size_t incrementBy);

        bool incrementNumScaScanShape(short dataType, size_t incrementBy);

        bool incrementAttRecordsShape(short dataType, size_t incrementBy);

        bool incrementOrbRecordsShape(short dataType, size_t incrementBy);

        bool incrementTlmPacketsShape(short dataType, size_t incrementBy);

        bool incrementTiltSampleShape(short dataType, size_t incrementBy);

        void closeAllL1aFiles();

        void closeAndRemoveFile(short datatype);

};





#endif