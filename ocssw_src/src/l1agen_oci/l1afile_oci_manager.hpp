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


        // outlistData =  {
        //      1 : { 
        //          "data_type": "SOL",
        //          "file_name": "PACE_OCI_SOL.nc",
        //          "start": "2025-08-07T02:38:24.639",
        //          "end": "2025-08-07T03:49:34.639",
        //          "flag": 0,
        //          
        //      }
        //  }
        // key == DATA_TYPE id 
        std::map<short, std::map<std::string, std::string>> outlistData;

        // files that do not need to be merged. For now: Lunar
        std::vector<std::string> noMergeOutlistStrings;

        // track the last data type added so if files cannot be merged due to a large gap
        // we know what to write out right away
        short lastDataTypeSeen = -1;
    

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

        /**
         * @brief notes down all files generated and track their time, complete flag and data type
         * @param filename 
         * @param startTime 
         * @param endTime 
         * @param completeFlag 
         * @param dataType 
         */
        void addFileToOutlistBuffer(std::string filename, std::string startTime,
                                std::string endTime, int completeFlag, 
                                std::string dataType, short dataTypeId
        );

        /**
         * @brief write the buffer information into the provided file name
         * @param outlist output file
         */
        void dumpOutlistBuffer(std::string outlist);
        
        /**
         * @brief Note the last data type that was added to the buffer
         * @param dataType 
         */
        void updateLastDataTypeSeen(short dataType);

        /**
         * @brief if the next file's time gap is too large, move the last
         *         file out of the buffer and into the noMergeOutlist
         */
        void processPrevFile();

        bool filesWereGenerated();

};





#endif