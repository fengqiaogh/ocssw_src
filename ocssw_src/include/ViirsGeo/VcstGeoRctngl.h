/******************************************************************************
 * 
 * NAME: VcstGeoRctngl
 *
 * DESCRIPTION: Defines functions to interpolate geolocation data
 *
 *
 ******************************************************************************/

#ifndef VCST_GEO_RCTNGL_H
#define VCST_GEO_RCTNGL_H

#include <VcstViirsGeo.h>
#include <VcstGeoDataStructs.h>
#include <VcstPolarStereographicDataSet.h>
#include <VcstGeoRctnglStruct.h>
#include <VcstGeoParameters.h>

class VcstCmnGeo;

// Geolocates the decimated rectangle pointes
void geolocateDecim(ViirsGeoRctnglType* inRec,
        viirsSdrGeoPtrs* inPtrs,
        VcstCmnGeo *geoPtr,
        ViirsGeoDecimType* outDecim);

// Geolocates all points in the granule by interpolating to the full
// strucure using the decimated geolocation data
int geolocateFullFromDecim(int iscan,
        ViirsGeoRctnglType* inRec,
        viirsSdrGeoPtrs* inPtrs,
        ViirsGeoDecimType* inDecim,
        VcstCmnGeo* geoPtr,
        ViirsGeoFullType* outFull,
        int extractPixelLimits[2]);

// Interpolates grow/gcol and lat/lon values for full locations using 
// the decimated geolocation data
int interpLocations(int iscan,
        ViirsGeoRctnglType* inRec,
        viirsSdrGeoPtrs* inPtrs,
        ViirsGeoDecimType* inDecim,
        VcstPolarStereographicDataSet* psds,
        ViirsGeoFullType* outFull,
        int extractPixelLimits[2]);

// Interpolates satellite, sun, and moon angles for full locations using
// the decimated geolocation data
void interpAngles(int iscan,
        ViirsGeoRctnglType* inRec,
        viirsSdrGeoPtrs* inPtrs,
        ViirsGeoDecimType* inDecim,
        ViirsGeoFullType* outFull,
        int extractPixelLimits[2]);

// Performs quadratic interpolation of the decimated data to the full
// data structure
void quadInterp(const double* Yval_in,
        double Xmid,
        double Xend,
        int idxlim,
        double* Yval_out);

// Improves the satellite angles near nadir
int fixNadirSatAngles(int iscan,
        viirsSdrGeoPtrs* inPtrs,
        VcstCmnGeo *geoPtr,
        ViirsGeoFullType* outFull,
        int extractPixelLimits[2]);

// Improves the satellite angles near the pole
int fixPoleSatAngles(int iscan,
        viirsSdrGeoPtrs* inPtrs,
        ViirsGeoRctnglType* inRec,
        VcstCmnGeo *geoPtr,
        ViirsGeoFullType* outFull,
        int extractPixelLimits[2]);

#endif

