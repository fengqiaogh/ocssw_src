/*******************************************************************************
 *
 * NAME: VcstGeoLutInputItem
 *
 *
 *******************************************************************************/

#ifndef VcstGeoLutInputItem_h
#define VcstGeoLutInputItem_h

#include <VcstLutInputItem.h>

class VcstGeoLutInputItem : public VcstLutInputItem {
public:

    /**
     * Constructor
     *
     * @param lutName The LUT name that identifies the data buffer.
     */
    VcstGeoLutInputItem(const std::string& groupName, size_t size);

    /**
     * Destructor
     */
    virtual ~VcstGeoLutInputItem();

protected:

    /**
     * Convert endianness of LUT file.
     */
    int convertEndianness();

    int convert_GEO_PARAM();

    int convert_CMNGEO_SAA_COEFF();

};



#endif
