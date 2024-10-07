/**************************************************************************
 *
 * NAME: VcstGeoAncInputItem
 *
 * DESCRIPTION: See description preceding class declaration below
 *
 * REFERENCES:
 * none
 *
 * LIMITATIONS:
 * none
 *
 * NOTES (MISCELLANEOUS) SECTION:
 * none
 *
 * --------- -------- -----------------------------------
 * <header review information>
 *
 **************************************************************************/

#ifndef _VcstGeoAncInputItem_h_
#define _VcstGeoAncInputItem_h_

#include <string>
#include <map>


#include <VcstLutInputItem.h>

#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;

/**
 * Derived class of AdlCmnInputItem.  This class represents a single
 * SRTM30Plus tile.  The difference between this class and an AdlCmnInput item
 * is that this class must use a tile ID to retrieve the appropriate file from
 * disk.
 */
class VcstGeoAncInputItem : public VcstLutInputItem {
public:

    /**
     * Constructor
     *
     * @param groupName The group name of this data item
     * @param size The size of this data item's buffer
     * @param tileId The tile ID of this item
     */
    VcstGeoAncInputItem(const std::string& groupName,
            size_t size,
            const std::string& tileId);

    /**
     * Destructor
     */
    virtual ~VcstGeoAncInputItem();


    /**
     * Determines the CSN the config guide.
     *
     * @param env Pointer to the environment that this data item is part of.
     *
     * @return PRO_SUCCESS or PRO_FAIL
     */
    virtual int setCSN();

    /**
     * Accessors
     */
    const std::string& getTileId() const;

    /**
     * Mutators
     */
    void setTileId(const std::string& tileId);

    /**
     * Handles the details of retrieving the data item from disk.
     *
     * @return PRO_SUCCESS or PRO_FAIL
     */
    virtual int getData();

protected:

    /**
     * Checks attributes of this data item to determine if it is able to be
     * read from disk
     */
    virtual bool isItemOkForIO();


    /**
     * Convert endianness of LUT file.
     */
    virtual int convertEndianness();

    /**
     * Read netCDF DEM LUT file.
     * @return PRO_SUCCESS or PRO_FAIL
     */
    virtual int readNetCDF();

private:

    /**
     * Default constructor.
     */
    VcstGeoAncInputItem();

    /**
     * Copy constructor.
     */
    VcstGeoAncInputItem(const VcstGeoAncInputItem& right);

    /**
     * Assignment operator.
     */
    VcstGeoAncInputItem& operator=(const VcstGeoAncInputItem& right);

    /**
     * Inserts file metadata
     */
    int readMetadata(std::string filepath);

    //Map to store terrain DEM tiles.
    std::map<std::string, std::string> ecoTileMap;

    int setTerrainTileMap();
    bool demBigEndian;
    bool demNetcdf;

    /**
     * The tile ID of the tile that this class represents.
     */
    std::string tileId_;

};


#endif
