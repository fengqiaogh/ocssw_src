/**************************************************************************
 *
 * NAME: VcstCalculateGridCenter.h
 *
 * DESCRIPTION: contains function declarations for functions used in the
 *  ProSdrViirs geolocation routines.
 *
 **************************************************************************/

#include <VcstGeoRctnglStruct.h>
#include <VcstPolarStereographicDataSet.h>

template <typename T>

void calculateGridCenter
(
        int scanIdx,
        // const int inActScans,
        const int inDet,
        const int inMidRow,
        const ViirsGeoRctnglType& inRect,
        VcstPolarStereographicDataSet* inPsds,
        T* gridPtr
        ) {
    // Fill the map data set in the grid
    inPsds->fillMds(&gridPtr->gmds);

    int row;
    int col;
    int width;

    //  for(int scanIdx=0; scanIdx < inActScans; scanIdx++)
    // {
    // Compute the row number of the middle row of a scan
    //  row = (scanIdx * inDet) + inMidRow;
    row = inMidRow;

    // Calculate the cntr grid row/col for each rectangle in the scan
    for (int rectIdx = 0; rectIdx < inRect.numRctngl; rectIdx++) {
        // Determine the width of the rectangle
        width = inRect.rcol[rectIdx] - inRect.lcol[rectIdx] + 1;

        // Capture the middle column of this rectangle
        col = inRect.mcol[rectIdx];

        // No matter what the width of the rectangle is, computing
        // the cntr grid row/col will require at least 2 pixels in
        // the middle column of the rectangle.
        // Ensure that these two pixels are not fill.
        if ((gridPtr->grow[row][col] < FLOAT32_FILL_TEST) ||
                (gridPtr->grow[row + 1][col] < FLOAT32_FILL_TEST)) {
            gridPtr->ctr_grow[scanIdx][rectIdx] = ERR_FLOAT64_FILL;
            gridPtr->ctr_gcol[scanIdx][rectIdx] = ERR_FLOAT64_FILL;
        } else {
            // If the width of the rectangle is an odd number, then
            // there is an exact middle column.  Use the 2 grid
            // row/col values to calculate the cntr grid row/col
            if ((width % 2) != 0) {
                gridPtr->ctr_grow[scanIdx][rectIdx] =
                        (gridPtr->grow[row][col] +
                        gridPtr->grow[row + 1][col]) / 2.e0;

                gridPtr->ctr_gcol[scanIdx][rectIdx] =
                        (gridPtr->gcol[row][col] +
                        gridPtr->gcol[row + 1][col]) / 2.e0;
            } else {
                // The rectangle width is an even number, therefore there
                // is no exact middle column.  Use the 4 surrounding grid
                // row/col values to calculate the cntr grid row/col
                if ((gridPtr->grow[row][col + 1] < FLOAT32_FILL_TEST) ||
                        (gridPtr->grow[row + 1][col + 1] < FLOAT32_FILL_TEST)) {
                    gridPtr->ctr_grow[scanIdx][rectIdx] = ERR_FLOAT64_FILL;
                    gridPtr->ctr_gcol[scanIdx][rectIdx] = ERR_FLOAT64_FILL;
                } else {
                    gridPtr->ctr_grow[scanIdx][rectIdx] =
                            (gridPtr->grow[row][col] +
                            gridPtr->grow[row][col + 1] +
                            gridPtr->grow[row + 1][col] +
                            gridPtr->grow[row + 1][col + 1]) / 4.e0;

                    gridPtr->ctr_gcol[scanIdx][rectIdx] =
                            (gridPtr->gcol[row][col] +
                            gridPtr->gcol[row][col + 1] +
                            gridPtr->gcol[row + 1][col] +
                            gridPtr->gcol[row + 1][col + 1]) / 4.e0;
                }
            }
        }
    }
    // }
}


