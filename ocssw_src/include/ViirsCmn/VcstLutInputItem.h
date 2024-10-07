/*******************************************************************************
 *
 * NAME: VcstLutInputItem.h
 *
 * DESCRIPTION: Base class for all Geolocation and Calibration LUT Input Item
 * object classes.  Provides basic functionality to read binary files and manage
 * data associated with characteristics of the files.
 *
 * Based on various IDPS object classes published by Raytheon Company that support
 * similar functionality.
 *
 *******************************************************************************/

#ifndef VcstLutInputItem_h
#define VcstLutInputItem_h

#include <string>

class VcstLutInputItem {
public:

    /**
     * Constructor
     *
     * @param lutName The LUT name that identifies the data buffer.
     */

    VcstLutInputItem(const std::string& groupName, size_t size);

    /**
     * Destructor
     */

    virtual ~VcstLutInputItem();

    /**
     * Retrieve instance
     */

    static VcstLutInputItem* getInstance() {
        return new VcstLutInputItem;
    }

    /**
     * Handles the details of retrieving the data item from disk.
     *
     * @return PRO_SUCCESS or PRO_FAIL
     */

    virtual int getData();

    inline void* getDataPtr() const {
        return dataPtr_;
    }

    inline bool dataPtrUnknown() const {
        return !dataPtrEstablished();
    }

    const std::string getShortName() const;

    void setShortName(const std::string& shortName);

    const std::string getFilePath() const;

    void setFilePath(const std::string& filePath);

    const std::string getFileName() const;

    inline void setDataPtr(void* ptr) {
        dataPtr_ = ptr;
    }

protected:

    static const std::size_t MAX_IO_SIZE = 1073741824;

    /**
     * Checks attributes of this data item to determine if it is able to be
     * read from disk
     */

    virtual bool isItemOkForIO();

    /**
     * Inline gets
     */

    inline const std::string getGroupName() const {
        return groupName_;
    }

    inline const size_t getDataSize() const {
        return dataSize_;
    }

    inline const bool dataPtrEstablished() const {
        return (dataPtr_ != 0);
    }

    /**
     * Inline sets
     */

    inline void setGroupName(std::string strName) {
        groupName_ = strName;
    }

    inline void setDataSize(size_t size) {
        dataSize_ = size;
    }

    inline void setIsAllocated(bool bVal) {
        isAllocated_ = bVal;
    }

    /**
     * Determine file size on disk
     */

    size_t determineFileSize(const std::string& fileName);

    /**
     * Allocate memory
     */

    int allocateMemory(size_t size);

    /**
     * Release allocated memory
     */

    void releaseMemory();

    /**
     * Convert endianness of LUT file.
     */

    virtual int convertEndianness();

    /**
     * Converts the endianness of the parameter by performing the appropriate
     * byte swapping.
     *
     */

    template<typename T>
    static void byteSwap(T& aValue);

    /**
     * This data item's group name
     */

    std::string groupName_;

private:

    /**
     * Default constructor.
     */

    VcstLutInputItem();

    /**
     * Copy constructor.
     */

    VcstLutInputItem(const VcstLutInputItem& right);

    /**
     * Assignment operator.
     */

    VcstLutInputItem& operator=(const VcstLutInputItem& right);

    /**
     * Inserts file metadata
     */

    int readMetadata(std::string filepath);

    /**
     * Pointer to the data buffer
     */

    void* dataPtr_;

    /**
     * Size of the data buffer (bytes) pointed to by the data pointer.
     */

    size_t dataSize_;

    /**
     * Set to indicate memory has been allocated
     */

    bool isAllocated_;

    /**
     * Set to indicate big endian file requiring conversion
     */

    bool bLittleEndian_;

    /**
     * The short name.
     */

    std::string shortName_;

    /**
     * The file path.
     */

    std::string filePath_;

};

//---------------------------------------------------------------------------
// Converts the endianness of the parameter by performing the appropriate
// byte swapping.
//---------------------------------------------------------------------------

template<typename T>
void VcstLutInputItem::byteSwap(T& aValue) {
    T tempValue = aValue; // a temporary copy of the value

    // Pointers to the first byte of both variables
    unsigned char* aValuePtr = reinterpret_cast<unsigned char*> (&aValue);
    unsigned char* tempValuePtr = reinterpret_cast<unsigned char*> (&tempValue);

    // Swap the byte order
    for (unsigned int byte = 0; byte < sizeof (aValue); ++byte) {
        aValuePtr[byte] = tempValuePtr[(sizeof (aValue) - 1) - byte];
    }
}

#endif
