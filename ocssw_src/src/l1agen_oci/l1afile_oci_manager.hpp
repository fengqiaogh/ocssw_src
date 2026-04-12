#ifndef _L1AFILE_OCI_MANAGER_
#define _L1AFILE_OCI_MANAGER_

#include <map>
#include "l1afile_oci.hpp"
#include "dimension_shape.hpp"
#include "navigation_timeframe.hpp"

// Manages writing to multiple L1A files for OCI. When there's a datatype change in-between OCI science data packets, 
// create a new l1a file and write to it. When it ends, resume writing into the science l1a file.
//
// It wraps around L1aFile class and manages multiple versions of it if there's a datatype change.

class L1aFileManager {
    public:
        L1aFileManager();
        ~L1aFileManager();

        // most recent ancillarty packet time (in unix time) for each data type
        // each index == time for the data type.
        // ie: science data index == 1, so lastSeenPacketTime[1] is the last time seen for
        // science data
        std::vector<double> prevAncillaryPktTimes;
        


        //********************************************************* */
        //*************** File Tracking Variables ***************** */
        //********************************************************* */

        
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

    

        /**
         * @brief update the previous time tracking for dataType with new time
         * @param dataType 
         * @param time 
         */
        void updatePrevAncillaryPktTimeFor(short dataType, double time);


        /**
         * @brief grab the last ancillary packet time for dataType file
         * @param dataType 
         * @return 
         */
        double getPrevAncillaryPktTimeFor(short dataType);
    


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

        /**
         * @brief grabs and returns a reference to the L1A file with dataType
         * @param dataType 
         * @return 
         */
        L1aFile* getL1aFile(short dataType);

        /**
         * @brief grabs and returns a reference to the dimension shape for 
         *          the L1A file with dataType 
         * @param dataType 
         * @return 
         */
        DimensionShape* getCurrL1aFileDimShape(short dataType);

        bool incrementNumScansShape(short dataType, size_t incrementBy);

        bool incrementNumMceScanShape(short dataType, size_t incrementBy);

        bool incrementNumScaScanShape(short dataType, size_t incrementBy);

        bool incrementAttRecordsShape(short dataType, size_t incrementBy);

        bool incrementOrbRecordsShape(short dataType, size_t incrementBy);

        bool incrementTlmPacketsShape(short dataType, size_t incrementBy);

        bool incrementTiltSampleShape(short dataType, size_t incrementBy);

        /**
         * @brief goes through all currently opened L1A files and closes them
         */
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
         * @param returnStatus update return status if any errors happen during dumping
         */
        void dumpOutlistBuffer(std::string outlist);
        
        /**
         * @brief Note the last data type that was added to the buffer
         * @param dataType 
         */
        void updateLastDataTypeSeen(short dataType);

        /**
         * @brief Generate outlist string for the data type in the format:
         *          "fileName startTime endTime completeFlag dataType"
         * @param dataType 
         * @return std::string 
         */
        std::string processOutlistStringFor(short dataType);

        /**
         * @brief if the next file's time gap is too large, move the last
         *         file out of the buffer and into the noMergeOutlist
         * @param dataType file type to process
         */
        void processOutlistFor(short dataType);

        bool filesWereGenerated();


        //********************************************************* */
        //****** navigation coverage tracking Variables *********** */
        //********************************************************* */
        // Only Scinece files will be using this dataType == 1
        // if multiple due to file time gaps, then the previous one should be closed
        // and reset this navTimeFrame

        NavigationTimeFrame navTimeFrame;

        /**
         * @brief pointer to the NavigationTimeFrame object so it can be passed into a function
         *          and get updated when processing scan times
         * @return 
         */
        NavigationTimeFrame* getNavTimeFrameRef();

};





#endif