/*******************************************************************************
 *
 * NAME: VcstPolarWanderItem.h
 *
 * DESCRIPTION: This class is derived from the AdlCmnInputItem class. It represents
 * a HDF5 Polar Wander data item capable of converting an IDPS Polar Wander AUX
 * file into a format that is compatible with the interface requirements of the
 * common geolocation algorithm.
 *
 * It is adapted from a utility application developed by
 *
 * Author:      Richard P. Cember
 *              Computational Physics, Inc. and Integrated Program Office
 *
 * Modified by: Paul E. Meade
 *
 *  Description: A source file to convert the USNO polar motion "finals" file
 *               into the corresponding binary file that is read by the SDR processes,
 *               with shortName USNO-POLARWANDER-UT1 and fitting into struct type
 *               CMNGEO_POLAR_UT1_TYPE.  This file is a highly modified version of
 *               Build 1.5.0.18's ING/MSD/EarthOrientation/src/IngMsdEarth-
 *               Orientation_Converter.cpp and its corresponding header (.h) file.
 *
 *******************************************************************************/

#ifndef VcstPolarWanderItem_h
#define VcstPolarWanderItem_h

#include <string>
#include <VcstCmnGeoStructs.h>
#include <VcstLutInputItem.h>

/**
 * This class is derived from the AdlCmnInputItem class. It represents a
 * HDF5 Polar Wander data item.
 */

class VcstPolarWanderItem : public VcstLutInputItem {
public:

    /**
     * Constructor
     *
     * @param groupName The group name that identifies the data buffer.
     * @param client The DMS client pointer.  Null for ADL, but present
     *        to retain the proper data item interface.
     * @param caller who instantiated this object
     * @param collectionShortName The CSN for this input item
     */
    VcstPolarWanderItem();

    /**
     * Destructor
     */
    virtual ~VcstPolarWanderItem();

    /**
     * Retrieve instance
     */
    static VcstPolarWanderItem* getInstance() {
        return new VcstPolarWanderItem;
    }

    /**
     * Handles the details of retrieving the data item from disk.
     *
     * @return PRO_SUCCESS or PRO_FAIL
     */
    virtual int getData();

protected:

    /**
     * Copy constructor.
     */
    VcstPolarWanderItem(const VcstPolarWanderItem& right);

    /**
     * Assignment operator.
     */
    VcstPolarWanderItem& operator=(const VcstPolarWanderItem& right);

    /**
     * Inserts file metadata
     */
    int readMetadata(std::string filepath);

    /**
     * Method to read an ASCII polar wander file into memory
     */
    int ReadAsciiPolarWanderFile(std::string filepath, char* bufPtr);

    /**
     * Method to read an HDF5 polar wander file into memory
     */
    int ReadH5PolarWanderFile(std::string filepath, char* bufPtr);

    /**
     * Method to copy the static anc out of its native format into an
     * internal format
     */
    int fillInternalBuffer(char* bufPtr);

    /**
     * Extracts a field from the source row and places it into a string stream
     * as a null terminated character array.
     *
     * @param converter Stream used to hold extracted value
     *
     * @param source Source buffer for the data being extracted
     *
     * @param offset Offset into the source buffer for the data being
     *               extracted
     *
     * @param length Length of the data being extracted
     */

    void extractField(std::stringstream &converter, unsigned char* source,
            const unsigned int offset, const unsigned int length);

    /**
     * Interim structure used to hold an entry from a finals file before
     * converting it into the internal format.
     */
    struct finalsEntry {
        /**
         * Year of the entry
         */
        int year;

        /**
         * Month of the entry
         */
        int month;

        /**
         * Day of the entry
         */
        int day;

        /**
         * Modified julian date of the entry
         */
        double modifiedJulianDate;

        /**
         * Polar flag of the entry
         */
        char polarFlag;

        /**
         * X wander of the entry
         */
        double x_wan;

        /**
         * Error in the above for the entry
         */
        double x_err;

        /**
         * Y wander of the entry
         */
        double y_wan;

        /**
         * Error in the above for the entry
         */
        double y_err;

        /**
         * Difference flag of the entry
         */
        char diffFlag;

        /**
         * UT1mUTC value of the entry.
         */
        double UT1mUTC;
    };

    unsigned char* nativeBuffer_;
    int nativeBufferSize_;

};

#endif
