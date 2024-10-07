/**************************************************************************
 * 
 * NAME: VcstPolarTileDatabase
 *
 **************************************************************************/

#ifndef _VcstPolarTileDatabase_H_
#define _VcstPolarTileDatabase_H_

#include <VcstTileDatabase.h>   // base class
#include <VcstMapDataSet.h>     // Hemisphere
#include <string>

// forward declaration to be SSPM compliant
// #include is in .cpp
class VcstPolarStereographicDataSet;

/**
 * Derived from ProCmnTileDatabase.  Polar Stereographic projection
 * information, including NORTHERN and SOUTHERN hemisphere map data sets.
 */
class VcstPolarTileDatabase : public VcstTileDatabase {
public:

    /**
     * Constructor
     */
    VcstPolarTileDatabase();

    /*
     * Destructor
     */
    /* virtual */
    ~VcstPolarTileDatabase();

    /**
     * Access the requested Hemisphere Mds
     *
     * @param hem NORTHERN or SOUTHERN
     * @retval requested mds, or 0 if invalid (EQUATOR)
     */
    VcstPolarStereographicDataSet* getHemisphereMds(Hemisphere hem);

    /**
     * Return the string representation of the tile specified by the
     * row and column in the overall database (grid).
     *
     * @param hem NORTHERN or SOUTHERN
     * @param aDatabaseRow Row in the overall database
     * @param aDatabaseCol Column in the overall database
     * @retval string representation of the tile
     */
    std::string calculateSrc(Hemisphere hem,
            double aDatabaseRow,
            double aDatabaseCol) const;

    /**
     * Return the string representation of the tile specified by the row
     * and column in the overall database (grid) in the stringstream parameter.
     *
     * @param hem NORTHERN or SOUTHERN
     * @param aDatabaseRow Row in the overall database
     * @param aDatabaseCol Column in the overall database
     * @param aStr stringstream representation of the tile
     */
    void calculateSrc(Hemisphere hem,
            double aDatabaseRow,
            double aDatabaseCol,
            std::ostringstream& aStr) const;

    /**
     * Return the string representation of the tile specified by the
     * row and column in the overall database (grid).
     *
     * @param hem NORTHERN or SOUTHERN
     * @param aTileRow Row for the given tile
     * @param aTileCol Column for the given tile
     * @retval string representation of the tile
     */
    std::string calculateSrc(Hemisphere hem,
            int aTileRow,
            int aTileCol) const;

    /**
     * Return the string representation of the tile specified by the
     * row and column in the overall database (grid).
     *
     * @param hem NORTHERN or SOUTHERN
     * @param aTileRow Row for the given tile
     * @param aTileCol Column for the given tile
     * @param aStr stringstream representation of the tile
     */
    void calculateSrc(Hemisphere hem,
            int aTileRow,
            int aTileCol,
            std::ostringstream& aStr) const;

    /**
     * Return the string representation of the tile specified by the
     * row and column in the overall database (grid).
     *
     * @param hem NORTHERN or SOUTHERN
     * @param aTileNum Overall tile number
     * @retval string representation of the tile id
     */
    std::string calculateSrc(Hemisphere hem,
            int aTileNum) const;

    /**
     * Return the string representation of the tile id specified by the
     * row and column in the overall database (grid).
     *
     * @param hem NORTHERN or SOUTHERN
     * @param aTileNum Overall tile number
     * @param aStr stringstream representation of the tile id
     */
    void calculateSrc(Hemisphere hem,
            int aTileNum,
            std::ostringstream& aStr) const;

protected:

    /**
     * Northern Mds class container
     */
    VcstPolarStereographicDataSet* theNorthernMds_;

    /**
     * Southern Mds class container
     */
    VcstPolarStereographicDataSet* theSouthernMds_;

private:

    /**
     * Privitized copy constructor
     * 
     * @param rhs this class to copy
     */
    VcstPolarTileDatabase(const VcstPolarTileDatabase& rhs);

    /**
     * Privitized assignment operator
     * 
     * @param rhs this class to assign
     */
    VcstPolarTileDatabase& operator=(const VcstPolarTileDatabase& rhs);

};

#endif  // _VcstPolarTileDatabase_H_
