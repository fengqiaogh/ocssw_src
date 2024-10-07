#!/usr/bin/env python3
# coding: utf-8

"""
A script to perform searches of the EarthData Common Metadata Repository (CMR)
for satellite granule names and download links.
written by J.Scott on 2016/12/12 (joel.scott@nasa.gov)

Modified by Inia Soto on 2020/08/17 (inia.m.sotoramos@nasa.gov)
Replaced wget function with python requests
"""

def main():

    import argparse
    import os
    from datetime import timedelta
    from math import isnan
    from collections import OrderedDict

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description='''\
      This program perform searches of the EarthData Search (https://search.earthdata.nasa.gov/search) Common Metadata
      Repository (CMR) for satellite granule names given an OB.DAAC satellite/instrument and lat/lon/time point or range.
      
      Outputs:
         1) a list of OB.DAAC L2 satellite file granule names that contain the input criteria, per the CMR's records.
         2) a list of public download links to fetch the matching satellite file granules, per the CMR's records.
      
      Inputs:
        The argument-list is a set of --keyword value pairs.

      Example usage calls:
         fd_matchup.py --sat=modist --slat=23.0 --slon=170.0 --stime=2015-11-16T09:00:00Z --max_time_diff=8
         fd_matchup.py --sat=modist --stime=2015-11-15T09:00:00Z --etime=2015-11-17T09:00:00Z --slat=23.0 --elat=25.0 --slon=170.0 --elon=175.0
         fd_matchup.py --sat=modist --max_time_diff=4 --seabass_file=[your SB file name].sb
         fd_matchup.py --sat=modist --slat=23.0 --slon=170.0 --stime=2015-11-16T09:00:00Z --max_time_diff=8 --get_data=[Your path]

      Caveats:
        * This script is designed to work with files that have been properly
          formatted according to SeaBASS guidelines (i.e. Files that passed FCHECK).
          Some error checking is performed, but improperly formatted input files
          could cause this script to error or behave unexpectedly. Files
          downloaded from the SeaBASS database should already be properly formatted, 
          however, please email seabass@seabass.gsfc.nasa.gov and/or the contact listed
          in the metadata header if you identify problems with specific files.

        * It is always HIGHLY recommended that you check for and read any metadata
          header comments and/or documentation accompanying data files. Information 
          from those sources could impact your analysis.

        * Compatibility: This script was developed for Python 3.5.

        * Requires a valid .netrc file in the user home ($HOME), e.g.:
          machine urs.earthdata.nasa.gov login USERNAME password PASSWD

      ''',add_help=True)

    parser.add_argument('--sat', nargs=1, required=True, type=str, choices=['modisa','modist','viirsn','viirsj1','goci','meris','czcs','octs','seawifs'], help='''\
      String specifier for satellite platform/instrument
      
      Valid options are:
      -----------------
      modisa  = MODIS on AQUA
      modist  = MODIS on TERRA
      viirsn  = VIIRS on Suomi-NPP
      viirsj1 = VIIRS on JPSS-1/NOAA-20
      meris   = MERIS on ENVISAT
      goci    = GOCI on COMS
      czcs    = CZCS on Nimbus-7
      seawifs = SeaWiFS on OrbView-2
      octs    = OCTS on ADEOS-I
      ''')

    parser.add_argument('--data_type', nargs=1, type=str, default=(['*']), choices=['oc','iop','sst'], help='''\
      OPTIONAL: String specifier for satellite data type
      Default behavior returns all product suites
      
      Valid options are:
      -----------------
      oc   = Returns OC (ocean color) product suite
      iop  = Returns IOP (inherent optical properties) product suite
      sst  = Returns SST product suite (including SST4 where applicable)
      ''')

    parser.add_argument('--slat', nargs=1, type=float, help=('''\
      Starting latitude, south-most boundary
      If used with --seabass_file, will override lats in the file
      Valid values: (-90,90N)
      '''))

    parser.add_argument('--elat', nargs=1, type=float, help=('''\
      Ending latitude, north-most boundary
      If used with --seabass_file and --slat, will override lats in the file
      Valid values: (-90,90N)
      '''))

    parser.add_argument('--slon', nargs=1, type=float, help=('''\
      Starting longitude, west-most boundary
      If used with --seabass_file, will override lons in the file
      Valid values: (-180,180E)
      '''))

    parser.add_argument('--elon', nargs=1, type=float, help=('''\
      Ending longitude, east-most boundary
      If used with --seabass_file and --slon, will override lons in the file
      Valid values: (-180,180E)
      '''))

    parser.add_argument('--stime', nargs=1, type=str, help='''\
      Time (point) of interest in UTC
      Default behavior: returns matches within +/- MAX_TIME_DIFF (default +/-3 hours) about this given time
      If used with ETIME, this creates a search time window, between STIME and ETIME.
      Valid format: string of the form: yyyy-mm-ddThh:mm:ssZ
      OPTIONALLY: Use with --max_time_diff or --etime
      ''')

    parser.add_argument('--max_time_diff', nargs=1, type=float, default=([3]), help=('''\
      Maximum time difference between satellite and in situ point
      OPTIONAL: default value +/-3 hours
      Valid values: decimal number of hours (0-36)
      Use with --seabass_file OR --stime
      '''))

    parser.add_argument('--etime', nargs=1, type=str, help='''\
      Maximum time (range) of interest in UTC
      Valid format: string of the form: yyyy-mm-ddThh:mm:ssZ
      Use with --stime
     ''')

    parser.add_argument('--seabass_file', nargs='+', type=argparse.FileType('r'), help='''\
      Valid SeaBASS file name or list of file names
      File must contain latitude, longitude, and date-time information as fields.
      ''')

    parser.add_argument('--get_data', nargs=1, type=str, help='''\
      Flag to download all identified satellite granules.
      Requires the use of an HTTP request.
      Set to the desired output directory.
      ''')

    parser.add_argument('--verbose', default=False, action='store_true', help=('''\
      OPTIONAL: Displays HTTP requests for each Earthdata CMR query.
      '''))

    args=parser.parse_args()

    if not args.sat:
        parser.error("you must specify an satellite string to conduct a search")
    else:
        dict_args=vars(args)
        sat = dict_args['sat'][0]

    #dictionary of lists of CMR platform, instrument, collection names
    dict_plat = {}
    dict_plat['modisa']  = ['MODIS','AQUA','MODISA_L2_']
    dict_plat['modist']  = ['MODIS','TERRA','MODIST_L2_']
    dict_plat['viirsn']  = ['VIIRS','Suomi-NPP','VIIRSN_L2_']
    dict_plat['viirsj1'] = ['VIIRS','NOAA-20','VIIRSJ1_L2_']
    dict_plat['meris']   = ['MERIS','ENVISAT','MERIS_L2_']
    dict_plat['goci']    = ['GOCI','COMS','GOCI_L2_']
    dict_plat['czcs']    = ['CZCS','Nimbus-7','CZCS_L2_']
    dict_plat['seawifs'] = ['SeaWiFS','OrbView-2','SeaWiFS_L2_']
    dict_plat['octs']    = ['OCTS','ADEOS-I','OCTS_L2_']

    if sat not in dict_plat:
        parser.error('you provided an invalid satellite string specifier. Use -h flag to see a list of valid options for --sat')

    if args.get_data:
        if not dict_args['get_data'][0] or not os.path.exists(dict_args['get_data'][0]):
            parser.error('invalid --get_data target download directory provided.')

    if dict_args['max_time_diff'][0] < 0 or dict_args['max_time_diff'][0] > 36:
        parser.error('invalid --max_time_diff value provided. Please specify an value between 0 and 36 hours. Received --max_time_diff = ' + str(dict_args['max_time_diff'][0]))
    else:
        twin_Hmin = -1 * int(dict_args['max_time_diff'][0])
        twin_Mmin = -60 * (dict_args['max_time_diff'][0] - int(dict_args['max_time_diff'][0]))
        twin_Hmax = 1 * int(dict_args['max_time_diff'][0])
        twin_Mmax = 60 * (dict_args['max_time_diff'][0] - int(dict_args['max_time_diff'][0]));

    #beginning of file/loop if-condition
    if args.seabass_file:
    
        for filein_sb in dict_args['seabass_file']:

            ds = check_SBfile(parser, filein_sb.name)
            
            granlinks = OrderedDict()
            hits = 0

            if args.slat and args.slon and args.elat and args.elon:
                check_lon(parser, dict_args['slon'][0])
                check_lon(parser, dict_args['elon'][0])

                check_lat(parser, dict_args['slat'][0])
                check_lat(parser, dict_args['elat'][0])

                check_lat_relative(parser, dict_args['slat'][0], dict_args['elat'][0])
                check_lon_relative(parser, dict_args['slon'][0], dict_args['elon'][0])

                #loop through times in file
                for dt in ds.datetime:
                    tim_min = dt + timedelta(hours=twin_Hmin,minutes=twin_Mmin) #use as: tim_min.strftime('%Y-%m-%dT%H:%M:%SZ')
                    tim_max = dt + timedelta(hours=twin_Hmax,minutes=twin_Mmax) #use as: tim_max.strftime('%Y-%m-%dT%H:%M:%SZ')

                    url = 'https://cmr.earthdata.nasa.gov/search/granules.json?page_size=2000' + \
                                '&provider=OB_DAAC' + \
                                '&bounding_box=' + str(dict_args['slon'][0]) + ',' + str(dict_args['slat'][0]) + ',' + \
                                                   str(dict_args['elon'][0]) + ',' + str(dict_args['elat'][0]) + \
                                '&instrument=' + dict_plat[sat][0] + \
                                '&platform=' + dict_plat[sat][1] + \
                                '&short_name=' + dict_plat[sat][2] + dict_args['data_type'][0] + \
                                '&options[short_name][pattern]=true' + \
                                '&temporal=' + tim_min.strftime('%Y-%m-%dT%H:%M:%SZ') + ',' + tim_max.strftime('%Y-%m-%dT%H:%M:%SZ') + \
                                '&sort_key=short_name'

                    if dict_args['verbose']:
                        print(url)

                    content = send_CMRreq(url)

                    [hits, granlinks] = process_CMRreq(content, hits, granlinks)

            elif args.slat and args.slon and not args.elat and not args.elon:
                check_lon(parser, dict_args['slon'][0])
                check_lat(parser, dict_args['slat'][0])

                #loop through times in file
                for dt in ds.datetime:
                    tim_min = dt + timedelta(hours=twin_Hmin,minutes=twin_Mmin) #use as: tim_min.strftime('%Y-%m-%dT%H:%M:%SZ')
                    tim_max = dt + timedelta(hours=twin_Hmax,minutes=twin_Mmax) #use as: tim_max.strftime('%Y-%m-%dT%H:%M:%SZ')

                    url = 'https://cmr.earthdata.nasa.gov/search/granules.json?page_size=2000' + \
                                '&provider=OB_DAAC' + \
                                '&point=' + str(dict_args['slon'][0]) + ',' + str(dict_args['slat'][0]) + \
                                '&instrument=' + dict_plat[sat][0] + \
                                '&platform=' + dict_plat[sat][1] + \
                                '&short_name=' + dict_plat[sat][2] + dict_args['data_type'][0] + \
                                '&options[short_name][pattern]=true' + \
                                '&temporal=' + tim_min.strftime('%Y-%m-%dT%H:%M:%SZ') + ',' + tim_max.strftime('%Y-%m-%dT%H:%M:%SZ') + \
                                '&sort_key=short_name'

                    if dict_args['verbose']:
                        print(url)

                    content = send_CMRreq(url)

                    [hits, granlinks] = process_CMRreq(content, hits, granlinks)

            else:
                ds = check_SBfile_latlon(parser, ds)

                for lat,lon,dt in zip(ds.lat,ds.lon,ds.datetime):
                    if isnan(lat) or isnan(lon):
                        continue
                    check_lon(parser, lon)
                    check_lat(parser, lat)

                    tim_min = dt + timedelta(hours=twin_Hmin,minutes=twin_Mmin) #use as: tim_min.strftime('%Y-%m-%dT%H:%M:%SZ')
                    tim_max = dt + timedelta(hours=twin_Hmax,minutes=twin_Mmax) #use as: tim_max.strftime('%Y-%m-%dT%H:%M:%SZ')

                    url = 'https://cmr.earthdata.nasa.gov/search/granules.json?page_size=2000' + \
                                '&provider=OB_DAAC' + \
                                '&point=' + str(lon) + ',' + str(lat) + \
                                '&instrument=' + dict_plat[sat][0] + \
                                '&platform=' + dict_plat[sat][1] + \
                                '&short_name=' + dict_plat[sat][2] + dict_args['data_type'][0] + \
                                '&options[short_name][pattern]=true' + \
                                '&temporal=' + tim_min.strftime('%Y-%m-%dT%H:%M:%SZ') + ',' + tim_max.strftime('%Y-%m-%dT%H:%M:%SZ') + \
                                '&sort_key=short_name'

                    if dict_args['verbose']:
                        print(url)

                    content = send_CMRreq(url)

                    [hits, granlinks] = process_CMRreq(content, hits, granlinks)

            print_CMRreq(hits, granlinks, dict_plat[sat], args, dict_args, filein_sb.name)

    #end of file/loop if-condition
    #beginning of lat/lon/time if-condition
    else:
    
        granlinks = OrderedDict()
        
        #Define time vars from input
        if args.stime and not args.etime:
            dt = check_time(parser, dict_args['stime'][0])

            tim_min = dt + timedelta(hours=twin_Hmin,minutes=twin_Mmin) #use as: tim_min.strftime('%Y-%m-%dT%H:%M:%SZ')
            tim_max = dt + timedelta(hours=twin_Hmax,minutes=twin_Mmax) #use as: tim_max.strftime('%Y-%m-%dT%H:%M:%SZ')

        elif args.stime and args.etime:
            tim_min = check_time(parser, dict_args['stime'][0])
            tim_max = check_time(parser, dict_args['etime'][0])
            
            check_time_relative(parser, tim_min, tim_max)

        else:
            parser.error('invalid time: All time inputs MUST be in UTC. Must receive --stime=YYYY-MM-DDTHH:MM:SSZ')

        #Define lat vars from input and call search query
        if args.slat and args.slon and not args.elat and not args.elon:
            check_lon(parser, dict_args['slon'][0])
            check_lat(parser, dict_args['slat'][0])

            url = 'https://cmr.earthdata.nasa.gov/search/granules.json?page_size=2000' + \
                            '&provider=OB_DAAC' + \
                            '&point=' + str(dict_args['slon'][0]) + ',' + str(dict_args['slat'][0]) + \
                            '&instrument=' + dict_plat[sat][0] + \
                            '&platform=' + dict_plat[sat][1] + \
                            '&short_name=' + dict_plat[sat][2] + dict_args['data_type'][0] + \
                            '&options[short_name][pattern]=true' + \
                            '&temporal=' + tim_min.strftime('%Y-%m-%dT%H:%M:%SZ') + ',' + tim_max.strftime('%Y-%m-%dT%H:%M:%SZ') + \
                            '&sort_key=short_name'

            if dict_args['verbose']:
                print(url)

            content = send_CMRreq(url)

        elif args.slat and args.elat and args.slon and args.elon:
            check_lon(parser, dict_args['slon'][0])
            check_lon(parser, dict_args['elon'][0])

            check_lat(parser, dict_args['slat'][0])
            check_lat(parser, dict_args['elat'][0])

            check_lat_relative(parser, dict_args['slat'][0], dict_args['elat'][0])
            check_lon_relative(parser, dict_args['slon'][0], dict_args['elon'][0])

            url = 'https://cmr.earthdata.nasa.gov/search/granules.json?page_size=2000' + \
                            '&provider=OB_DAAC' + \
                            '&bounding_box=' + str(dict_args['slon'][0]) + ',' + str(dict_args['slat'][0]) + ',' + \
                                                   str(dict_args['elon'][0]) + ',' + str(dict_args['elat'][0]) + \
                            '&instrument=' + dict_plat[sat][0] + \
                            '&platform=' + dict_plat[sat][1] + \
                            '&short_name=' + dict_plat[sat][2] + dict_args['data_type'][0] + \
                            '&options[short_name][pattern]=true' + \
                            '&temporal=' + tim_min.strftime('%Y-%m-%dT%H:%M:%SZ') + ',' + tim_max.strftime('%Y-%m-%dT%H:%M:%SZ') + \
                            '&sort_key=short_name'

            if dict_args['verbose']:
                print(url)

            content = send_CMRreq(url)

        else:
            parser.error('invalid combination of --slat and --slon OR --slat, --elat, --slon, and --elon arguments provided. All latitude inputs MUST be between -90/90N deg. All longitude inputs MUST be between -180/180E deg.')

        #Parse json return for the lat/lon/time if-condition
        processANDprint_CMRreq(content, granlinks, dict_plat[sat], args, dict_args, tim_min, tim_max)

    return


def check_SBfile(parser, file_sb):
    """ function to verify SB file exists, is valid, and has correct fields; returns data structure """
    import os
    from seabass.SB_support import readSB

    if os.path.isfile(file_sb):
        ds = readSB(filename=file_sb, 
                    mask_missing=True, 
                    mask_above_detection_limit=True, 
                    mask_below_detection_limit=True, 
                    no_warn=True)
    else:
        parser.error('ERROR: invalid --seabass_file specified. Does: ' + file_sb + ' exist?')

    ds.datetime = ds.fd_datetime()
    if not ds.datetime:
        parser.error('missing fields in SeaBASS file. File must contain date/time, date/hour/minute/second, year/month/day/time, OR year/month/day/hour/minute/second')
    
    return ds


def check_SBfile_latlon(parser, ds):
    """ function to verify lat/lon exist in SB file's data structure """

    import re
    from numpy import mean

    ds.lon = []
    ds.lat = []

    if 'lat' in ds.data and 'lon' in ds.data:

        for lat,lon in zip(ds.data['lat'],ds.data['lon']):
            check_lat(parser, float(lat))
            check_lon(parser, float(lon))
            ds.lat.append(float(lat))
            ds.lon.append(float(lon))

    elif 'north_latitude' in ds.headers and \
        'south_latitude' in ds.headers and \
        'west_longitude' in ds.headers and \
        'east_longitude' in ds.headers:

        lat_n = re.search("([+|-]?\d*\.?\d*)\[(deg|DEG)\]", ds.headers['north_latitude'])
        lat_s = re.search("([+|-]?\d*\.?\d*)\[(deg|DEG)\]", ds.headers['south_latitude'])
        lon_w = re.search("([+|-]?\d*\.?\d*)\[(deg|DEG)\]", ds.headers['west_longitude'])
        lon_e = re.search("([+|-]?\d*\.?\d*)\[(deg|DEG)\]", ds.headers['east_longitude'])

        try:

            #for i in range(ds.length):

            ds.lat.append(mean([float(lat_n.group(1)), float(lat_s.group(1))]))
            ds.lon.append(mean([float(lon_w.group(1)), float(lon_e.group(1))])) #won't handle a dateline crossing well

        except:

            parser.error('/north_latitude, /south_latitude, /west_longitude, or /east_longitude headers not formatted correctly; unable to parse in file: {:}'.format(ds.filename))

    else:

        parser.error('missing headers/fields in SeaBASS file. File must contain lat,lon information')

    return ds


def check_lat(parser, lat):
    """ function to verify lat range """
    if abs(lat) > 90.0:
        parser.error('invalid latitude: all LAT values MUST be between -90/90N deg. Received: ' + str(lat))
    return


def check_lon(parser, lon):
    """ function to verify lon range """
    if abs(lon) > 180.0:
        parser.error('invalid longitude: all LON values MUST be between -180/180E deg. Received: ' + str(lon))
    return

def check_lat_relative(parser, slat, elat):
    """ function to verify two lats relative to each other """
    if slat > elat:
        parser.error('invalid latitude: --slat MUST be less than --elat. Received --slat = ' + str(slat) + ' and --elat = ' + str(elat))
    return


def check_lon_relative(parser, slon, elon):
    """ function to verify two lons relative to each other """
    if slon > elon:
        parser.error('invalid longitude: --slon MUST be less than --elon. Received --slon = ' + str(slon) + ' and --elon = ' + str(elon))
    return


def check_time(parser, tim):
    """ function to verify time """
    import re
    from datetime import datetime

    try:
        tims = re.search("(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})Z", tim);
        dt = datetime(year=int(tims.group(1)), \
                      month=int(tims.group(2)), \
                      day=int(tims.group(3)), \
                      hour=int(tims.group(4)), \
                      minute=int(tims.group(5)), \
                      second=int(tims.group(6)))
    except:
        parser.error('invalid time: All time inputs MUST be in UTC in the form: YYYY-MM-DDTHH:MM:SSZ Received: ' + tim)
    return dt


def check_time_relative(parser, tim_min, tim_max):
    """ function to verify two times relative to each other """
    if tim_min > tim_max:
        parser.error('invalid time: --stime MUST be less than --etime. Received --stime = ' + \
                     tim_min.strftime('%Y-%m-%dT%H:%M:%SZ') + ' and --etime = ' + \
                     tim_max.strftime('%Y-%m-%dT%H:%M:%SZ'))
    return


def send_CMRreq(url):
    """ function to submit a given URL request to the CMR; return JSON output """
    import requests

    req = requests.get(url)
    content = req.json()

    return content


#def send_CMRreq(url):
#    """ function to submit a given URL request to the CMR; return JSON output """
#    from urllib import request
#    import json
#
#    req = request.Request(url)
#    req.add_header('Accept', 'application/json')
#    content = json.loads(request.urlopen(req).read().decode('utf-8'))
#
#    return content


def process_CMRreq(content, hits, granlinks):
    """ function to process the return from a single CMR JSON return """

    try:
        hits = hits + len(content['feed']['entry'])
        for entry in content['feed']['entry']:
            granid = entry['producer_granule_id']
            granlinks[granid] = entry['links'][0]['href']
    except:
        print('WARNING: No matching granules found for a row. Continuing to search for granules from the rest of the input file...')

    return hits, granlinks


def download_file(url, out_dir):
    '''
    download_file downloads a file
    given URL and out_dir strings
    syntax fname_local = download_file(url, out_dir)
    '''
    import os,subprocess
    import requests
    from pathlib import Path

    out_path = Path(out_dir).expanduser() #even if already a Path object, this convert should work
    path_local = out_path / url.split('/')[-1]

    print('Downloading',url.split('/')[-1],'to',out_path)

    #Inia M. Soto Ramos @ 2020.08.13
    #Remove wget and use request to download the file and prevent wget issues
    #cmd1 = 'wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies ' + \
    #       '--auth-no-challenge=on --keep-session-cookies --content-disposition ' + \
    #       f'{url} -O {path_local.as_posix()}'

    if not out_path.exists():
        os.makedirs(str(out_dir))

    if path_local.exists():

        print(f'{path_local.name} found - deleting and attempting to re-download...')
        path_local.unlink()

    else:

        try:
            print(f'downloading {path_local.name}')
            #s = subprocess.run(cmd1, encoding='ascii', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #Added by Inia Soto 2020.08.
            req2 = requests.get(url)
            with open(path_local.as_posix(), 'wb') as f:
                f.write(req2.content)

        except Exception as e:
            if 'HTTPError' in e:
                print('Error in download_file:',e)
                return -1

            else:
                print('Unknown error:', e)

        if not path_local.exists():
            print('Error in download_file: local copy of file not found; download unsuccessful')
            return -2

        print('Successfully downloaded',url.split('/')[-1],'to',out_path)

    return str(path_local)


def print_CMRreq(hits, granlinks, plat_ls, args, dict_args, sbfile):
    """" function to print the CMR results from a SB file """

    if hits > 0:
        unique_hits = 0
        for granid in granlinks:
            unique_hits = unique_hits + 1

            print('Matching ' + plat_ls[1] + '/' + plat_ls[0] + ' granule (' + granid + ') found for: ' + sbfile)

            if args.get_data and dict_args['get_data'][0]:
                fname_loc = download_file(granlinks[granid], dict_args['get_data'][0])
            else:
                print('Download link: ' + granlinks[granid])
            
            print(' ')
        
        print('Number of granules found: ' + str(unique_hits))
        print(' ')
        
    else:
        print('WARNING: No granules found for ' + plat_ls[1] + '/' + plat_ls[0] + ' and any lat/lon/time inputs.')

    return


def processANDprint_CMRreq(content, granlinks, plat_ls, args, dict_args, tim_min, tim_max):
    """ function to process AND print the return from a single CMR JSON return """

    try:
        hits = len(content['feed']['entry'])
        for entry in content['feed']['entry']:
            granid = entry['producer_granule_id']
            granlinks[granid] = entry['links'][0]['href']

            print('Matching ' + plat_ls[1] + '/' + plat_ls[0] + ' granule (' + granid + ') found.')

            if args.get_data and dict_args['get_data'][0]:
                fname_loc = download_file(granlinks[granid], dict_args['get_data'][0])
            else:
                print('Download link: ' + granlinks[granid])

            print(' ')
        
        print('Number of granules found: ' + str(hits))
        print(' ')

    except:
        print('WARNING: No matching granules found for ' + plat_ls[1] + '/' + plat_ls[0] + \
              ' containing the requested lat/lon area during the ' + \
              str(dict_args['max_time_diff'][0]) + '-hr window of ' + \
              tim_min.strftime('%Y-%m-%dT%H:%M:%SZ') + ' to ' + tim_max.strftime('%Y-%m-%dT%H:%M:%SZ'))

    return


if __name__ == "__main__": main()
