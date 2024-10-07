#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A script to convert PRISM data from the CORAL project 
to netCDF proxy data in the format of PACE-OCI
written by J.Scott on 2017/10/03 (joel.scott@nasa.gov)
"""

def main():

    import argparse
    import os
    from PRISM_module import readL1B, readL2, readL3
    from netCDF4 import Dataset
    from datetime import datetime, timedelta
    import numpy as np

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description='''\
      This script converts PRISM data from the CORAL project to netCDF
      Outputs:
        A single netCDF4 OB.DAAC L2 file of L1, L2, and L3 PRISM data products from the CORAL project
      Inputs:
        The argument-list is a set of --keyword value pairs.
      Example usage:
        prismcoral2nc.py --l1dir=$HOME/prm20160908t234616_rdn_v1p1 --l2dir=$HOME/prm20160908t234616_rb_v1p1 --l3dir=$HOME/prm20160908t234616_bc_v1p1
      Caveats:
        Compatibility: This script was developed for Python 3.6
      License:
        /*=====================================================================*/
                         NASA Goddard Space Flight Center (GSFC) 
                 Software distribution policy for Public Domain Software

         The fd_matchup.py code is in the public domain, available without fee for 
         educational, research, non-commercial and commercial purposes. Users may 
         distribute this code to third parties provided that this statement appears
         on all copies and that no charge is made for such copies.

         NASA GSFC MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THE SOFTWARE
         FOR ANY PURPOSE. IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED
         WARRANTY. NEITHER NASA GSFC NOR THE U.S. GOVERNMENT SHALL BE LIABLE FOR
         ANY DAMAGE SUFFERED BY THE USER OF THIS SOFTWARE.
        /*=====================================================================*/
      ''',add_help=True)

    parser.add_argument('--l1dir', nargs=1, required=True, type=str, help='''\
      String specifier for PRISM L1(B) data dir 
      (i.e. - folder where L1 tar.gz bundle was unzipped to)
      ''')

    parser.add_argument('--l2dir', nargs=1, required=True, type=str, help='''\
      String specifier for PRISM L2 data dir 
      (i.e. - folder where L2 tar.gz bundle was unzipped to)
      ''')

    parser.add_argument('--l3dir', nargs=1, required=True, type=str, help='''\
      String specifier for PRISM L3 data dir 
      (i.e. - folder where L3 tar.gz bundle was unzipped to)
      ''')

    parser.add_argument('--odir', nargs=1, type=str, default=(['./']), help=('''\
      OPTIONAL: output directory
      If none provided, netCDF data will be saved to the current working directory
      '''))

    args=parser.parse_args()

    if not args.l1dir or not args.l2dir or not args.l3dir:
        parser.error("you must specify l1dir, l2dir, and l3dir")
    else:
        dict_args=vars(args)

    #remove trailing slashes
    if '/' in dict_args['l1dir'][0][-1]:
        dict_args['l1dir'][0] = dict_args['l1dir'][0][:-1]

    if '/' in dict_args['l2dir'][0][-1]:
        dict_args['l2dir'][0] = dict_args['l2dir'][0][:-1]

    if '/' in dict_args['l3dir'][0][-1]:
        dict_args['l3dir'][0] = dict_args['l3dir'][0][:-1]

    if '/' in dict_args['odir'][0][-1]:
        dict_args['odir'][0] = dict_args['odir'][0][:-1]

    #load L1 data
    print('Reading PRISM L1 files from:',dict_args['l1dir'][0])
    dsL1 = readL1B(dict_args['l1dir'][0])
    #load L2 data
    print('Reading PRISM L2 files from:',dict_args['l2dir'][0])
    dsL2 =  readL2(dict_args['l2dir'][0])
    #load L3 data
    print('Reading PRISM L3 files from:',dict_args['l3dir'][0])
    dsL3 =  readL3(dict_args['l3dir'][0])

    #mk netcdf filename
    dt_start = dsL1.hdr['rdn_img']['start datetime']
    nc_fname = dt_start.strftime('PRM%Y%j%H%M%S.L2.nc')

    #mk netCDF4 file
    nc_fullpath = dict_args['odir'][0] + '/' + nc_fname
    print('Creating netCDF4 file',nc_fullpath)
    nc_fid = Dataset(nc_fullpath, 'w', clobber=True, format='NETCDF4') #clobber=True will over-write file, if it exists

    #mk global attributes
    print('Creating global attributes')
    nc_fid.title = "PRISM Level-2 Data"
    nc_fid.sensor = "PRISM (Portable Remote Imaging SpectroMeter)"
    nc_fid.instrument = "PRISM"
    nc_fid.platform = "Tempus Applied Solutions Gulfstream-IV"
    nc_fid.product_name = nc_fname
    nc_fid.processing_version = "V1.0"
    nc_fid.prism_version = dsL1.hdr['rdn_img']['prism_version']
    nc_fid.Conventions = "CF-1.6"
    nc_fid.institution = "NASA Jet Propulsion Laboratory/California Institute of Technology"
    nc_fid.license = "https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/"
    nc_fid.naming_authority = "gov.nasa.gsfc.sci.oceandata"
    dt_now = datetime.now()
    nc_fid.date_created = dt_now.strftime("%Y-%m-%dT%H:%M:%SZ")
    nc_fid.keywords_vocabulary = "NASA Global Change Master Directory (GCMD) Science Keywords"
    nc_fid.stdname_vocabulary = "NetCDF Climate and Forecast (CF) Metadata Convention"
    nc_fid.creator_name = "Jet Propulsion Laboratory/California Institute of Technology"
    nc_fid.creator_email = "sarah.r.lundeen@jpl.nasa.gov"
    nc_fid.creator_url = "https://prism.jpl.nasa.gov/index.html"
    nc_fid.project = "CORAL EVS-2"
    nc_fid.project_url = "https://coral.jpl.nasa.gov/"
    nc_fid.publisher_name = "NASA Goddard Space Flight Center (GSFC), OB.DAAC"
    nc_fid.publisher_url = "https://oceancolor.gsfc.nasa.gov"
    nc_fid.publisher_email = "data@oceancolor.gsfc.nasa.gov"
    nc_fid.identifier_product_doi_authority = "http://dx.doi.org"
    nc_fid.identifier_product_doi = "10.5067/PRISM/CORAL/L2/RBEN/V1"
    nc_fid.processing_level = "L2"
    nc_fid.spatialResolution = dsL1.hdr['rdn_img']['map info'][5] + " m"
    nc_fid.flight_line = dt_start.strftime("prm%Y%j%H%M%S")
    nc_fid.history = "netCDF4 file created by OCSSW's prismcoral2nc.py using Level-1, Level-2, and Level-3 flight line data from the Portable Remote Imaging SpectroMeter (PRISM) instrument mounted on Tempus Applied Solutions Gulfstream-IV (G-IV) aircraft flying at a nominal operating altitude of 8.5 km from the COral Reef Airborne Laboratory (CORAL) Earth Venture Suborbital-2 (EVS-2) project. Input files were " + dict_args['l1dir'][0].split('/')[-1] + ', ' + dict_args['l2dir'][0].split('/')[-1] + ', and ' + dict_args['l3dir'][0].split('/')[-1];
    nc_fid.time_coverage_start = dt_start.strftime("%Y-%m-%dT%H:%M:%SZ")
    dt_end = dt_start + timedelta(microseconds=np.nanmax(dsL1.data['utc time']))
    nc_fid.time_coverage_end = dt_end.strftime("%Y-%m-%dT%H:%M:%SZ")
    nc_fid.northernmost_latitude = float(np.nanmax(dsL1.data['latitude (wgs-84)']))
    nc_fid.southernmost_latitude = float(np.nanmin(dsL1.data['latitude (wgs-84)']))
    nc_fid.easternmost_longitude = float(np.nanmax(dsL1.data['longitude (wgs-84)']))
    nc_fid.westernmost_longitude = float(np.nanmin(dsL1.data['longitude (wgs-84)']))
    nc_fid.geospatial_lat_units = "degrees_north" ;
    nc_fid.geospatial_lon_units = "degrees_east" ;
    nc_fid.geospatial_lat_max = float(np.nanmax(dsL1.data['latitude (wgs-84)']))
    nc_fid.geospatial_lat_min = float(np.nanmin(dsL1.data['latitude (wgs-84)']))
    nc_fid.geospatial_lon_max = float(np.nanmax(dsL1.data['longitude (wgs-84)']))
    nc_fid.geospatial_lon_min = float(np.nanmin(dsL1.data['longitude (wgs-84)']))
    nc_fid.day_night_flag = "Day"
    nc_fid.earth_sun_distance_correction = float(np.nanmean(dsL1.data['earth-sun distance (au)']))

    #mk dimensions
    print('Creating dimensions')
    nc_fid.createDimension('number_of_scans', dsL1.data['utc time'].shape[0])
    nc_fid.createDimension('pixels_per_scan', dsL1.data['utc time'].shape[1])
    nc_fid.createDimension('number_of_bands', len(dsL1.hdr['rdn_img']['wavelength']))
    nc_fid.createDimension('number_of_benthic_bands', len(dsL2.hdr['rb_img']['wavelength']))

    #mk groups
    print('Creating groups')
    nc_fid.createGroup('/sensor_band_parameters')
    nc_fid.createGroup('/scan_line_attributes')
    nc_fid.createGroup('/navigation_data')
    nc_fid.createGroup('/observation_data')
    nc_fid.createGroup('/derived_data')

    #mk sbp-group vars
    print('Creating and writing sensor_band_parameters')
    nc_wavelength = nc_fid.createVariable('/sensor_band_parameters/wavelength', 
                                          "f4", 
                                          dimensions=('number_of_bands',), 
                                          zlib=True, 
                                          complevel=5, 
                                          chunksizes=(nc_fid.dimensions['number_of_bands'].size,), 
                                          fill_value=-999.)
    nc_wavelength.long_name = 'Band center wavelengths'
    nc_wavelength.units = 'nm'
    nc_wavelength.valid_min = 350.0
    nc_wavelength.valid_max = 1050.0
    nc_wavelength[:] = [float(d) for d in dsL1.hdr['rdn_img']['wavelength']]

    nc_fwhm = nc_fid.createVariable('/sensor_band_parameters/fwhm', 
                                    "f4", 
                                    dimensions=('number_of_bands',), 
                                    zlib=True, 
                                    complevel=5, 
                                    chunksizes=(nc_fid.dimensions['number_of_bands'].size,), 
                                    fill_value=-999.)
    nc_fwhm.long_name = 'Band full-width half-maximums'
    nc_fwhm.units = 'nm'
    nc_fwhm.valid_min = 0.0
    nc_fwhm.valid_max = 100.0
    nc_fwhm[:] = [float(d) for d in dsL1.hdr['rdn_img']['fwhm']][:nc_fid.dimensions['number_of_bands'].size]

    nc_rbwavelength = nc_fid.createVariable('/sensor_band_parameters/benthic_wavelength', 
                                          "f4", 
                                          dimensions=('number_of_benthic_bands',), 
                                          zlib=True, 
                                          complevel=5, 
                                          chunksizes=(nc_fid.dimensions['number_of_benthic_bands'].size,), 
                                          fill_value=-999.)
    nc_rbwavelength.long_name = 'Band center benthic reflectance wavelengths'
    nc_rbwavelength.units = 'nm'
    nc_rbwavelength.valid_min = 350.0
    nc_rbwavelength.valid_max = 1035.0
    nc_rbwavelength[:] = [float(d) for d in dsL2.hdr['rb_img']['wavelength']]

    #mk sla-group vars
    print('Creating and writing scan_line_attributes')
    nc_scan_start_time = nc_fid.createVariable('/scan_line_attributes/scan_start_time', 
                                    "f8", 
                                    dimensions=('number_of_scans',), 
                                    zlib=True, 
                                    complevel=5, 
                                    chunksizes=(nc_fid.dimensions['number_of_scans'].size,), 
                                    fill_value=-999.)
    nc_scan_start_time.long_name = 'Scan start time (UTC)'
    nc_scan_start_time.units = 'seconds'
    nc_scan_start_time.valid_min = 0.
    nc_scan_start_time.valid_max = 2000000000.
    dsL1.data['scan_start_time'][np.isnan(dsL1.data['scan_start_time'])] = nc_scan_start_time._FillValue
    nc_scan_start_time[:] = dsL1.data['scan_start_time']

    nc_scan_end_time = nc_fid.createVariable('/scan_line_attributes/scan_end_time', 
                                    "f8", 
                                    dimensions=('number_of_scans',), 
                                    zlib=True, 
                                    complevel=5, 
                                    chunksizes=(nc_fid.dimensions['number_of_scans'].size,), 
                                    fill_value=-999.)
    nc_scan_end_time.long_name = 'Scan end time (UTC)'
    nc_scan_end_time.units = 'seconds'
    nc_scan_end_time.valid_min = 0.
    nc_scan_end_time.valid_max = 2000000000.
    dsL1.data['scan_end_time'][np.isnan(dsL1.data['scan_end_time'])] = nc_scan_end_time._FillValue
    nc_scan_end_time[:] = dsL1.data['scan_end_time']

    #mk nd-group vars
    print('Creating and writing navigation_data')
    nc_lon = nc_fid.createVariable('/navigation_data/longitude', 
                                   "f4", 
                                   dimensions=('number_of_scans','pixels_per_scan'), 
                                   zlib=True, 
                                   complevel=5, 
                                   chunksizes=(1,nc_fid.dimensions['pixels_per_scan'].size), 
                                   fill_value=-999.)
    nc_lon.long_name = 'Longitudes of pixel locations'
    nc_lon.units = 'degrees_east'
    nc_lon.valid_min = -180.
    nc_lon.valid_max = 180.
    dsL1.data['longitude (wgs-84)'][np.isnan(dsL1.data['longitude (wgs-84)'])] = nc_lon._FillValue
    nc_lon[:,:] = dsL1.data['longitude (wgs-84)']

    nc_lat = nc_fid.createVariable('/navigation_data/latitude', 
                                   "f4", 
                                   dimensions=('number_of_scans','pixels_per_scan'), 
                                   zlib=True, 
                                   complevel=5, 
                                   chunksizes=(1,nc_fid.dimensions['pixels_per_scan'].size), 
                                   fill_value=-999.)
    nc_lat.long_name = 'Latitudes of pixel locations'
    nc_lat.units = 'degrees_north'
    nc_lat.valid_min = -90.
    nc_lat.valid_max = 90.
    dsL1.data['latitude (wgs-84)'][np.isnan(dsL1.data['latitude (wgs-84)'])] = nc_lat._FillValue
    nc_lat[:,:] = dsL1.data['latitude (wgs-84)']

    nc_altitude = nc_fid.createVariable('/navigation_data/height', 
                                           "f4", 
                                           dimensions=('number_of_scans','pixels_per_scan'), 
                                           zlib=True, 
                                           complevel=5, 
                                           chunksizes=(1,nc_fid.dimensions['pixels_per_scan'].size), 
                                           fill_value=-999.)
    nc_altitude.long_name = 'Terrain height at pixel locations'
    nc_altitude.units = 'meters'
    nc_altitude.valid_min = -1000.
    nc_altitude.valid_max = 10000.
    dsL1.data['elevation (m)'][np.isnan(dsL1.data['elevation (m)'])] = nc_altitude._FillValue
    nc_altitude[:,:] = dsL1.data['elevation (m)']

    nc_range = nc_fid.createVariable('/navigation_data/range', 
                                     "i2", 
                                     dimensions=('number_of_scans','pixels_per_scan'), 
                                     zlib=True, 
                                     complevel=5, 
                                     chunksizes=(1,nc_fid.dimensions['pixels_per_scan'].size), 
                                     fill_value=-999.)
    nc_range.long_name = 'Aircraft-to-pixel range'
    nc_range.units = 'meters'
    nc_range.valid_min = 0.
    nc_range.valid_max = 25000.
    dsL1.data['path length (m)'][np.isnan(dsL1.data['path length (m)'])] = nc_range._FillValue
    nc_range[:,:] = dsL1.data['path length (m)']

    nc_senz = nc_fid.createVariable('/navigation_data/sensor_zenith', 
                                    "f4", 
                                    dimensions=('number_of_scans','pixels_per_scan'), 
                                    zlib=True, 
                                    complevel=5, 
                                    chunksizes=(1,nc_fid.dimensions['pixels_per_scan'].size), 
                                    fill_value=-999.)
    nc_senz.long_name = 'Sensor zenith angle at pixel locations'
    nc_senz.units = 'degrees'
    nc_senz.valid_min = 0.
    nc_senz.valid_max = 90.
    dsL1.data['to-sensor zenith (0 to 90 degrees from zenith)'][np.isnan(dsL1.data['to-sensor zenith (0 to 90 degrees from zenith)'])] = nc_senz._FillValue
    nc_senz[:,:] = dsL1.data['to-sensor zenith (0 to 90 degrees from zenith)']

    nc_sena = nc_fid.createVariable('/navigation_data/sensor_azimuth', 
                                    "f4", 
                                    dimensions=('number_of_scans','pixels_per_scan'), 
                                    zlib=True, 
                                    complevel=5, 
                                    chunksizes=(1,nc_fid.dimensions['pixels_per_scan'].size), 
                                    fill_value=-999.)
    nc_sena.long_name = 'Sensor azimuth angle at pixel locations'
    nc_sena.units = 'degrees'
    nc_sena.valid_min = 0.
    nc_sena.valid_max = 360.
    dsL1.data['to-sensor azimuth (0 to 360 degrees cw from n)'][np.isnan(dsL1.data['to-sensor azimuth (0 to 360 degrees cw from n)'])] = nc_sena._FillValue
    nc_sena[:,:] = dsL1.data['to-sensor azimuth (0 to 360 degrees cw from n)']

    nc_solz = nc_fid.createVariable('/navigation_data/solar_zenith', 
                                    "f4", 
                                    dimensions=('number_of_scans','pixels_per_scan'), 
                                    zlib=True, 
                                    complevel=5, 
                                    chunksizes=(1,nc_fid.dimensions['pixels_per_scan'].size), 
                                    fill_value=-999.)
    nc_solz.long_name = 'Solar zenith angle at pixel locations'
    nc_solz.units = 'degrees'
    nc_solz.valid_min = 0.
    nc_solz.valid_max = 90.
    dsL1.data['to-sun zenith (0 to 90 degrees from zenith)'][np.isnan(dsL1.data['to-sun zenith (0 to 90 degrees from zenith)'])] = nc_solz._FillValue
    nc_solz[:,:] = dsL1.data['to-sun zenith (0 to 90 degrees from zenith)']

    nc_sola = nc_fid.createVariable('/navigation_data/solar_azimuth', 
                                    "f4", 
                                    dimensions=('number_of_scans','pixels_per_scan'), 
                                    zlib=True, 
                                    complevel=5, 
                                    chunksizes=(1,nc_fid.dimensions['pixels_per_scan'].size), 
                                    fill_value=-999.)
    nc_sola.long_name = 'Solar azimuth angle at pixel locations'
    nc_sola.units = 'degrees'
    nc_sola.valid_min = 0.
    nc_sola.valid_max = 360.
    dsL1.data['to-sun azimuth (0 to 360 degrees cw from n)'][np.isnan(dsL1.data['to-sun azimuth (0 to 360 degrees cw from n)'])] = nc_sola._FillValue
    nc_sola[:,:] = dsL1.data['to-sun azimuth (0 to 360 degrees cw from n)']

    #mk od-group vars
    print('Creating and writing observation_data')
    nc_lt = nc_fid.createVariable('/observation_data/Lt', 
                                  "f4", 
                                  dimensions=('number_of_bands','number_of_scans','pixels_per_scan'), 
                                  zlib=True, 
                                  complevel=5, 
                                  chunksizes=(nc_fid.dimensions['number_of_bands'].size,1,nc_fid.dimensions['pixels_per_scan'].size), 
                                  fill_value=-999.)
    nc_lt.long_name = 'Top of atmosphere radiance'
    nc_lt.units = 'uW cm-2 nm-1 sr-1'
    nc_lt.valid_min = 0.
    nc_lt.valid_max = 800.
    dsL1.data['lt'][np.isnan(dsL1.data['lt'])] = nc_lt._FillValue
    nc_lt[:,:,:] = dsL1.data['lt']

    #mk dd-group vars
    print('Creating and writing derived_data')
    nc_rrs = nc_fid.createVariable('/derived_data/Rrs', 
                                   "f4", 
                                   dimensions=('number_of_bands','number_of_scans','pixels_per_scan'), 
                                   zlib=True, 
                                   complevel=5, 
                                   chunksizes=(nc_fid.dimensions['number_of_bands'].size,1,nc_fid.dimensions['pixels_per_scan'].size), 
                                   fill_value=-999.)
    nc_rrs.long_name = 'Remote sensing reflectance'
    nc_rrs.units = 'sr-1'
    nc_rrs.valid_min = 0.
    nc_rrs.valid_max = 1000.
    nc_rrs.comment = 'Scaled estimates of the water-leaving reflectance directly above the water surface for each measured wavelength, scaled by Pi with the assumption of a Lambertian surface. This is comparable to the Hemispherical Directed Reflectance Function (Shaepman-Strub et al., 2006).  It is calculated using a variant of the ATREM atmospheric correction method (Gao et al., 1993; Gao et al., 2009).  The modifications are described in texts by Thompson et al. (2015a and 2015b).'
    nc_rrs.reference = 'Gao, B. C., Heidebrecht, K. B., & Goetz, A. F. (1993). Derivation of scaled surface reflectances from AVIRIS data. Remote sensing of Environment, 44(2-3), 165-178.; Gao, B. C., Montes, M. J., Davis, C. O., & Goetz, A. F. (2009). Atmospheric correction algorithms for hyperspectral remote sensing data of land and ocean. Remote Sensing of Environment, 113, S17-S24.; Schaepman-Strub, G., Schaepman, M. E., Painter, T. H., Dangel, S., & Martonchik, J. V. (2006). Reflectance quantities in optical remote sensing — Definitions and case studies. Remote sensing of environment, 103(1), 27-42.; Thompson, D. R., Gao, B. C., Green, R. O., Roberts, D. A., Dennison, P. E., & Lundeen, S. R. (2015a). Atmospheric correction for global mapping spectroscopy: ATREM advances for the HyspIRI preparatory campaign. Remote Sensing of Environment, 167, 64-77.; Thompson, D. R., Seidel, F. C., Gao, B. C., Gierach, M. M., Green, R. O., Kudela, R. M., & Mouroulis, P. (2015b). Optimizing irradiance estimates for coastal and inland water imaging spectroscopy. Geophysical Research Letters, 42(10), 4116-4123.'
    dsL2.data['rrs'][np.isnan(dsL2.data['rrs'])] = nc_rrs._FillValue
    nc_rrs[:,:,:] = dsL2.data['rrs']

    nc_rb = nc_fid.createVariable('/derived_data/Rb', 
                                   "f4", 
                                   dimensions=('number_of_benthic_bands','number_of_scans','pixels_per_scan'), 
                                   zlib=True, 
                                   complevel=5, 
                                   chunksizes=(nc_fid.dimensions['number_of_benthic_bands'].size,1,nc_fid.dimensions['pixels_per_scan'].size), 
                                   fill_value=-999.)
    nc_rb.long_name = 'Benthic reflectance'
    nc_rb.units = 'sr-1'
    nc_rb.valid_min = 0.
    nc_rb.valid_max = 1000.
    nc_rb.comment = 'Estimates of the diffuse bottom reflectance (i.e. the ratio of downwelling to upwelling flux) for each wavelength, formed by estimating the apparent optical properties AOPs) of the water column.  The calculation uses the relationship for shallow-water reflectance described by Maritorena et al. (1994). The bottom reflectance is modeled using a linear nonnegative combination of a set of one or more basis endmembers from a suitable library of bottom reflectances.  This process provides smooth-looking spectra across all wavelengths, but we note that values for longer wavelengths becomes less certain as the depth of the water increases.'
    nc_rb.reference = 'Maritorena, S., Morel, A., & Gentili, B. (1994). Diffuse reflectance of oceanic shallow waters: Influence of water depth and bottom albedo. Limnology and oceanography, 39(7), 1689-1703.'
    dsL2.data['rb'][np.isnan(dsL2.data['rb'])] = nc_rb._FillValue
    nc_rb[:,:,:] = dsL2.data['rb']

    nc_depth = nc_fid.createVariable('/derived_data/depth', 
                                   "f4", 
                                   dimensions=('number_of_scans','pixels_per_scan'), 
                                   zlib=True, 
                                   complevel=5, 
                                   chunksizes=(1,nc_fid.dimensions['pixels_per_scan'].size), 
                                   fill_value=-999.)
    nc_depth.long_name = 'water column depth'
    nc_depth.units = 'meters'
    nc_depth.valid_min = 0.
    nc_depth.valid_max = 10000.
    nc_depth.comment = 'Estimated depth of the water column used to derive benthic reflectance.'
    dsL2.data['depth (m)'][np.isnan(dsL2.data['depth (m)'])] = nc_depth._FillValue
    nc_depth[:,:] = dsL2.data['depth (m)']

    nc_coral = nc_fid.createVariable('/derived_data/coral', 
                                   "f4", 
                                   dimensions=('number_of_scans','pixels_per_scan'), 
                                   zlib=True, 
                                   complevel=5, 
                                   chunksizes=(1,nc_fid.dimensions['pixels_per_scan'].size), 
                                   fill_value=-999.)
    nc_coral.long_name = 'Probability of coral benthic cover type'
    nc_coral.units = 'unitless'
    nc_coral.valid_min = 0.
    nc_coral.valid_max = 1.
    nc_coral.comment = 'Derived probability associated with coral as the benthic cover classification for each seafloor pixel.'
    dsL3.data['coral'][np.isnan(dsL3.data['coral'])] = nc_coral._FillValue
    nc_coral[:,:] = dsL3.data['coral']

    nc_sand = nc_fid.createVariable('/derived_data/sand', 
                                   "f4", 
                                   dimensions=('number_of_scans','pixels_per_scan'), 
                                   zlib=True, 
                                   complevel=5, 
                                   chunksizes=(1,nc_fid.dimensions['pixels_per_scan'].size), 
                                   fill_value=-999.)
    nc_sand.long_name = 'Probability of sand benthic cover type'
    nc_sand.units = 'unitless'
    nc_sand.valid_min = 0.
    nc_sand.valid_max = 1.
    nc_sand.comment = 'Derived probability associated with sand as the benthic cover classification for each seafloor pixel.'
    dsL3.data['sand'][np.isnan(dsL3.data['sand'])] = nc_sand._FillValue
    nc_sand[:,:] = dsL3.data['sand']

    nc_algae = nc_fid.createVariable('/derived_data/algae', 
                                   "f4", 
                                   dimensions=('number_of_scans','pixels_per_scan'), 
                                   zlib=True, 
                                   complevel=5, 
                                   chunksizes=(1,nc_fid.dimensions['pixels_per_scan'].size), 
                                   fill_value=-999.)
    nc_algae.long_name = 'Probability of algae benthic cover type'
    nc_algae.units = 'unitless'
    nc_algae.valid_min = 0.
    nc_algae.valid_max = 1.
    nc_algae.comment = 'Derived probability associated with algae as the benthic cover classification for each seafloor pixel.'
    dsL3.data['algae'][np.isnan(dsL3.data['algae'])] = nc_algae._FillValue
    nc_algae[:,:] = dsL3.data['algae']

    #close nc file
    nc_fid.close()
    print('Successfully created', nc_fullpath)

    return


if __name__ == "__main__": main()
