/**************************************************************************
 * 
 * NAME: VcstTileDatabase
 *
 **************************************************************************/

#ifndef VcstTileDatabase_h
#define VcstTileDatabase_h

/**
 * Abstract Base Class.  Derived classes exist for Polar Stereographic and
 * Sinusoidal projections.  Use these tile database classes to obtain
 * dimensions including row and column.  Also obtain tile number.
 */
class VcstTileDatabase {
public:

    /*
     * Destructor
     *
     * Enables safe deleting of derived class object stored in base 
     * class pointer -- otherwise memory leak.
     */
    virtual ~VcstTileDatabase();

    /*
     * Returns the number of tiles in a row
     */
    int getNumOfTileRows() const;

    /*
     * Returns the number of tiles in a column
     */
    int getNumOfTileCols() const;

    /*
     * Returns the number of rows in a tile
     */
    int getNumOfRowsPerTile() const;

    /*
     * Returns the number of columns in a tile
     */
    int getNumOfColsPerTile() const;

    /*
     * Returns the total number of rows in the database (grid)
     */
    int getNumOfDbRows() const;

    /*
     * Returns the total number of columns in the database (grid)
     */
    int getNumOfDbCols() const;

    /*
     * Maps from database (grid) row to row offset inside tile
     * 
     * @param aDatabaseRow Overall row in database (grid)
     * @retval row offset inside tile from database row
     */
    virtual int calculateRowInTile(double aDatabaseRow) const;

    /*
     * Maps from database (grid) column to column offset inside tile
     * 
     * @param aDatabaseCol Overall column in database (grid)
     * @retval column offset inside tile from database column
     */
    virtual int calculateColInTile(double aDatabaseCol) const;

    /*
     * Maps from database (grid) row to tile row
     * 
     * @param aDatabaseRow Overall row in database (grid)
     * @retval tile row from database grid
     */
    virtual int calculateTileRow(double aDatabaseRow) const;

    /*
     * Maps from database (grid) column to tile column
     * 
     * @param aDatabaseCol Overall column in database (grid)
     * @retval tile column from database grid
     */
    virtual int calculateTileCol(double aDatabaseCol) const;

    /*
     * Maps from database (grid) row and column to tile number
     *
     * @param aDatabaseRow Overall row in database (grid)
     * @param aDatabaseCol Overall column in database (grid)
     * @retval tile number from database row and column
     */
    int calculateTileNum(double aDatabaseRow,
            double aDatabaseCol) const;

    /*
     * Maps from tile row and column to tile number
     *
     * @param aTileRow Tile row
     * @param aTileCol Tile column
     * @retval tile number from mapped row and column
     */
    int calculateTileNum(int aTileRow,
            int aTileCol) const;

    /*
     * Decompose the tile number into tile row/column.
     * @param aTileNum (in) a tile number to decompose.
     * @param aTileRow (out) place to put the tile row.
     * @param aTileCol (out) place to put the tile column.
     */
    virtual void decomposeTileNum(int aTileNum, int& aTileRow, int& aTileCol);

protected:

    /*
     * Constructor is protected forcing into abstract class
     *
     * @param aNumOfTileRows Total number of tile rows
     * @param aNumOfTileCols Total number of tile columns
     * @param aNumOfRowsPerTile Total number of rows in a tile
     * @param aNumOfColsPerTile Total number of columns in a tile
     */
    VcstTileDatabase(int aNumOfTileRows,
            int aNumOfTileCols,
            int aNumOfRowsPerTile,
            int aNumOfColsPerTile);

    /**
     * The number of tile rows
     */
    int theNumOfTileRows_;
    /**
     * The number of tile columns
     */
    int theNumOfTileCols_;
    /**
     * The number of rows per tile
     */
    int theNumOfRowsPerTile_;
    /**
     * The number of columns per tile
     */
    int theNumOfColsPerTile_;
    /**
     * The number of database rows
     */
    int theNumOfDbRows_;
    /**
     * The number of database columns
     */
    int theNumOfDbCols_;

private:

    /**
     * Privitized copy constructor
     *
     * @param rhs this class to copy
     */
    VcstTileDatabase(const VcstTileDatabase& rhs);
    /**
     * Privitized assignment operator
     *
     * @param rhs this class to assign
     */
    VcstTileDatabase& operator=(const VcstTileDatabase& rhs);

};


//------------------------------------------------------------------------------

inline int
VcstTileDatabase::getNumOfTileRows() const {
    return theNumOfTileRows_;
}

//------------------------------------------------------------------------------

inline int
VcstTileDatabase::getNumOfTileCols() const {
    return theNumOfTileCols_;
}

//------------------------------------------------------------------------------

inline int
VcstTileDatabase::getNumOfRowsPerTile() const {
    return theNumOfRowsPerTile_;
}

//------------------------------------------------------------------------------

inline int
VcstTileDatabase::getNumOfColsPerTile() const {
    return theNumOfColsPerTile_;
}

//------------------------------------------------------------------------------

inline int
VcstTileDatabase::getNumOfDbRows() const {
    return theNumOfDbRows_;
}

//------------------------------------------------------------------------------

inline int
VcstTileDatabase::getNumOfDbCols() const {
    return theNumOfDbCols_;
}

//------------------------------------------------------------------------------

inline int
VcstTileDatabase::calculateTileNum(double aDatabaseRow,
        double aDatabaseCol) const {
    return calculateTileNum(calculateTileRow(aDatabaseRow),
            calculateTileCol(aDatabaseCol));
}

//------------------------------------------------------------------------------

inline int
VcstTileDatabase::calculateTileNum(int aTileRow,
        int aTileCol) const {
    return (aTileRow * getNumOfTileCols() + aTileCol);
}

#endif  // VcstTileDatabase_h
