#!/usr/bin/env python3
import argparse
from pathlib import Path
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter, LatitudeLocator

def get_palette(palfile):

    palpath = None
    if os.getenv("OCDATAROOT"):
        dataroot = Path(os.getenv("OCDATAROOT"))

        palpath = dataroot / 'common' / 'palette' / palfile
        print(palpath)
        if not Path.exists(palpath):
            err_msg = "Cannot find palette %s" % palfile
            sys.exit("Cannot find palette {palfile}".format(palfile=str(palfile)))
        palpath = str(palpath)
    else:
        sys.exit("$OCDATAROOT needs to be defined")
    colors = []

    fp = open(palpath, "r")
    for line in fp:
        r,g,b = line.split()
        colors.append((int(r)/255.,int(g)/255.,int(b)/255.))


    cmap = mpl.colors.LinearSegmentedColormap.from_list("palfile",colors, gamma=1.0)
    return cmap

def get_crs_from_proj(projStr):
    crs = None
    #proj=+proj=aea +ellps=WGS84 +datum=WGS84 +lon_0=-80.000000 +lat_0=27.500000 +lat_1=15.000000 +lat_2=40.000000


    if 'EPSG' in projStr:
        print("Since I have to trust querying https://epsg.io/ (and need an internet connection to do so)...things might go wrong...but we'll give it a try :)")
        key,epsgcode = projStr.strip().split(':')
        crs = ccrs.epsg(epsgcode)
    else:
        proj={}
        proj['utmSouth'] = False
        for element in projStr.split():
            if '+south' in element:
                proj['utmSouth'] = True
            if '=' in element:
                key, value = element.split('=')
                proj[key[1:]] = value 
        if 'proj' in proj:
            if proj['proj'] == 'eqc':
                crs = ccrs.PlateCarree()
            if proj['proj'] == 'aea':
                crs = ccrs.AlbersEqualArea(central_longitude=float(proj['lon_0']), central_latitude=float(proj['lat_0']), standard_parallels=(float(proj['lat_1']), float(proj['lat_2'])))
            if proj['proj'] == 'moll':
                crs = ccrs.Mollweide(central_longitude=float(proj['lon_0']))
            if proj['proj'] == 'lcc':
                crs = ccrs.LambertConformal(central_longitude=float(proj['lon_0']), central_latitude=float(proj['lat_0']), standard_parallels=(float(proj['lat_1']), float(proj['lat_2'])))
            if proj['proj'] == 'merc':
                crs = ccrs.Mercator(central_longitude=float(proj['lon_0']))
            if proj['proj'] == 'tmerc':
                crs = ccrs.TransverseMercator(central_longitude=float(proj['lon_0']),central_latitude=float(proj['lat_0']), approx=False)
            if proj['proj'] == 'utm':
                crs = ccrs.UTM(proj['zone'], southern_hemisphere=proj['utmSouth'])
            if proj['proj'] == 'stere':
                crs = ccrs.Stereographic(central_longitude=float(proj['lon_0']), central_latitude=float(proj['lat_0']), true_scale_latitude=float(proj['lat_ts']))
            if proj['proj'] == 'ortho':
                crs = ccrs.Orthographic(central_longitude=float(proj['lon_0']), central_latitude=float(proj['lat_0']))

    return crs

def normalize_lon(lons):
    idx = np.where(lons - 180 < -360)
    if idx:
        lons[idx] += 360
    idx = np.where(lons + 180 > 360)
    if idx:
        lons[idx] -= 360
    return lons

def calculate_gridline_increments(minval,maxval):
    if maxval < minval:
        maxval += 360

    if ((maxval-minval) < 15): # If you're small, just fix it at an increment of 10
        grid = np.round(np.linspace(np.floor(minval-0.5),np.ceil(maxval+0.5),10),decimals=1)
    else:
        increment = int(np.floor((maxval-minval)/4)+1)
        grid = np.round(np.linspace(np.floor(minval-increment/2),np.ceil(maxval+increment/2),increment),decimals=0)

    return(grid)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="add overlays to mapped PNG images output from mapgen (or l3mapgen)")
    parser.add_argument('ifile', type=str, help='input mapped image file')
    parser.add_argument('--projinfo', '-p', type=str, default=None, help='input projection info file; default: <ifile>.projtxt')
    parser.add_argument('--ofile', '-o', type=str, default=None, help='output filename; default: <ifile basename>.overlay.png')
    parser.add_argument('--coastline', '-c', action='store_true', help='add coastline overlay')
    parser.add_argument('--gridline', '-g', action='store_true', help='add gridline overlay')
    parser.add_argument('--label', type=str, default=None, help='add label string (e.g. source file for the image)')
    parser.add_argument('--nogridlabels', action='store_true', help='do not lable gridlines')
    parser.add_argument('--colorbar', action='store_true',help='add a colorbar; requires datamin, datamax; optionally: scaletype, palfile, cbartitle')
    parser.add_argument('--datamin', type=float, default=None, help='minimum value for colorbar scale')
    parser.add_argument('--datamax', type=float, default=None, help='maximum value for colorbar scale')
    parser.add_argument('--scale_type', choices=['linear','log'], default=None, help='colorbar scaling method; default: linear')
    parser.add_argument('--palfile', type=str, default='default-black0.pal', help='colorbar palette name; see $OCDATAROOT/common/palette/')
    parser.add_argument('--cbartitle', type=str,default=None, help='colorbar title string')
    parser.add_argument('--minX', type=float, help='minimum longitude extent in meters')
    parser.add_argument('--maxX', type=float, help='maximum longitude extent in meters')
    parser.add_argument('--minY', type=float, help='minimum latitude extent in meters')
    parser.add_argument('--maxY', type=float, help='maximum latitude extent in meters')
    parser.add_argument('--north', '-n', type=float, help='maximum latitude extent in degrees')
    parser.add_argument('--south', '-s', type=float, help='minimum latitude extent in degrees')
    parser.add_argument('--west', '-w', type=float, help='minimum longitude extent in degrees')
    parser.add_argument('--east', '-e', type=float, help='maximum longitude extent in degrees')
    parser.add_argument('--proj', type=str, help='PRØJ-formatted projection string')
    parser.add_argument('--use_transparency',action='store_true', help='output transparent image')

    args = parser.parse_args()

    ofile = args.ofile
    if not ofile:
        ofilePath = Path(args.ifile)
        ofile = "{dir}/{stem}.{ext}".format(dir=ofilePath.parent,stem=ofilePath.stem, ext='overlay.png')

    pinfofile = args.projinfo
    if not pinfofile:
        pinfofile = "{ifile}.projtxt".format(ifile=args.ifile)

    pinfo = {}
    asString = ['proj','scale_type']
    with open(pinfofile,'r') as pj:
        for line in pj:
            if line.startswith('#'):
                continue
            key, value = line.split('=',1)
            if key in asString:
                pinfo[key] = value.rstrip()
            else:
                pinfo[key] = float(value.rstrip())

    pinfokeys = ['minX','maxX','minY','maxY','north','south','east','west','datamin','datamax','scale_type','proj']
    for key,value in (vars(args)).items():
        if key in pinfokeys and value:
            pinfo[key] = value

    img_extent_meters = (pinfo['minX'],pinfo['maxX'],pinfo['minY'],pinfo['maxY'])
    img_extent = (pinfo['west'],pinfo['east'],pinfo['south'],pinfo['north'])

    img = plt.imread(args.ifile)

    # set the figure with the size from the input image dimensions
    dpi = 72

    # get image size in inches and create figure with buffer around the image
    imgwidth = img.shape[1]/float(dpi)
    imgheight = img.shape[0]/float(dpi)
    # add padding
    padheight = np.ceil(imgwidth*0.3)
    padwidth = np.ceil(imgheight*0.3)

    figwidth = np.round(imgwidth + padwidth, decimals=2)
    figheight = np.round(imgheight + padheight, decimals=2)
    figrows = int(np.ceil(4*figheight))

    fig = plt.figure(figsize=(figwidth, figheight), dpi=dpi)

    # set the plot projection
    crs = get_crs_from_proj(pinfo['proj'])
    if not crs:
        sys.exit("Projection not supported: {proj}".format(proj=pinfo['proj']))
    bottom = int(np.floor(0.025*figrows))
    ax = plt.subplot2grid((figrows,7),(0,0),rowspan=figrows-bottom, colspan=7, fig=fig,projection=crs)

    # Image Extent
    ax.set_extent(img_extent_meters,crs=crs)

    # Create the image
    ax.imshow(img, extent=img_extent_meters, origin='upper', transform=crs)

    # Add image name
    if args.label:
        ax.set_title(args.label, fontsize=18)

    # Add coastlines
    if args.coastline:
        ax.coastlines(resolution='50m', color='white', linewidth=1)

    # Add gridlines
    labelsize = 14 * (figheight / 14.0)
    if figwidth < 8:
        labelsize = 10

    labelgrid = True
    if args.nogridlabels:
        labelgrid = False

    gl = ax.gridlines(draw_labels=labelgrid, dms=False,
                    x_inline=False, y_inline=False,
                    linewidth=2, color='gray', alpha=0.5)
    gl.bottom_labels = False
    if not args.gridline:
        gl.xlines = False
        gl.ylines = False
    gl.xlabel_style = {'size': labelsize}
    gl.ylabel_style = {'size': labelsize}

    lons = calculate_gridline_increments(pinfo['west'],pinfo['east'])
    lats = calculate_gridline_increments(pinfo['south'],pinfo['north'])

    gl.xlocator = mticker.FixedLocator(normalize_lon(lons))
    gl.ylocator = mticker.FixedLocator(lats)

    # add colorbar
    if args.colorbar:
        cbax = plt.subplot2grid((figrows,7),(figrows-bottom,2),fig=fig,rowspan=bottom,colspan=3)

        cmap = get_palette(args.palfile)

        norm = mpl.colors.Normalize(vmin=pinfo['datamin'], vmax=pinfo['datamax'])
        formatter = None

        if (pinfo['scale_type']).lower() == 'log':
            if pinfo['datamin'] <= 0:
                sys.exit("Log scales don't like zero (or less than zero)")
            formatter = mpl.ticker.LogFormatter(10, labelOnlyBase=False) 
            norm = mpl.colors.LogNorm(vmin=pinfo['datamin'], vmax=pinfo['datamax'])

        cb = mpl.colorbar.ColorbarBase(cbax, cmap=cmap,
                                norm=norm,
                                orientation='horizontal',
                                format=formatter)

        if args.cbartitle:
            cb.set_label(args.cbartitle,size=labelsize)
        cbax.tick_params(labelsize=labelsize)
        cbax.minorticks_on()

    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)
    # Save the image
    plt.savefig(ofile, dpi=dpi, bbox_inches='tight', pad_inches=0.2,transparent=args.use_transparency)