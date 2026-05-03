#! /usr/bin/env python3
import argparse
from pathlib import Path
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LogNorm, ListedColormap, Normalize
from matplotlib.cm import ScalarMappable
from PIL import Image
import os

def _crs_from_proj4(proj_string):
    # Parse a Proj.4 string and return an appropriate cartopy CRS subclass
    params = {}
    for token in proj_string.split():
        if '=' in token:
            key, val = token.lstrip('+').split('=', 1)
            params[key] = val
        else:
            params[token.lstrip('+')] = True

    proj   = params.get('proj', 'merc')
    lon_0  = float(params.get('lon_0', 0))
    lat_0  = float(params.get('lat_0', 0))

    ellps  = params.get('ellps', params.get('datum', None))
    globe  = ccrs.Globe(ellipse=ellps) if ellps else None
    gkw    = {'globe': globe} if globe else {}

    if proj == 'merc':
        return ccrs.Mercator(central_longitude=lon_0, **gkw)
    elif proj in ('eqc', 'longlat', 'latlong'):
        return ccrs.PlateCarree(central_longitude=lon_0)
    elif proj == 'lcc':
        lat_1 = float(params.get('lat_1', lat_0))
        lat_2 = float(params.get('lat_2', lat_0))
        return ccrs.LambertConformal(central_longitude=lon_0, central_latitude=lat_0,
                                     standard_parallels=(lat_1, lat_2), **gkw)
    elif proj == 'stere':
        return ccrs.Stereographic(central_longitude=lon_0, central_latitude=lat_0, **gkw)
    elif proj == 'ortho':
        return ccrs.Orthographic(central_longitude=lon_0, central_latitude=lat_0)
    elif proj == 'moll':
        return ccrs.Mollweide(central_longitude=lon_0)
    elif proj == 'robin':
        return ccrs.Robinson(central_longitude=lon_0)
    elif proj == 'sinu':
        return ccrs.Sinusoidal(central_longitude=lon_0)
    else:
        print(f"Warning: unrecognized projection '{proj}', falling back to Mercator")
        return ccrs.Mercator(central_longitude=lon_0, **gkw)

def annotate_l3mapgen(input_file, output_file=None, label_text=None):
    # Set path for Cartopy to find its downloaded shapefiles
    # Make sure the path is expanded *before* Cartopy tries to load anything
    cartopy_data_root = os.environ['OCDATAROOT'] + '/common/cartopy'
    cartopy.config['data_dir'] = os.path.expanduser(cartopy_data_root)
    cartopy.config['pre_existing_data_dir'] = os.path.expanduser(cartopy_data_root)

    input_path = Path(input_file)
    is_png = input_path.suffix.lower() == '.png'

    if is_png:
        img = Image.open(input_file)
        meta = img.info

        proj_string = meta.get('projString')
        if proj_string is None:
            raise ValueError("No 'projString' metadata found in PNG file")
        
        projection = _crs_from_proj4(proj_string)
        proj_name = proj_string.split('+proj=')[1].split()[0].title() if '+proj=' in proj_string else 'Unknown'
        print(f"Using projection from PNG metadata {proj_string}")

        resolution = float(meta['resolution'])
        min_x = float(meta['minX'])
        max_y = float(meta['maxY'])
        max_x = min_x + int(meta['width']) * resolution
        min_y = max_y - int(meta['height']) * resolution

        is_true_color = img.mode in ('RGB', 'RGBA')

        if img.mode == 'P':
            raw_indices = np.array(img)
            fill_idx = int(meta.get('fill_value_index', 0))

        img_data = np.array(img.convert('RGBA'))

        if img.mode == 'P':
            img_data[raw_indices == fill_idx] = [255, 255, 255, 0]

        data_var = input_path.stem

    else:        
            
        ds = nc.Dataset(input_file)

        # Determine projection from file; default to Mercator
        proj_attr = getattr(ds, "map_projection", None)
        if proj_attr:
            proj_attr = proj_attr.lower()
            if proj_attr in ["equidistant cylindrical", "platecarree", "plate carrée", "eqc"]:
                projection = ccrs.PlateCarree()
                proj_name = "Equidistant Cylindrical"
            else:
                projection = ccrs.Mercator()
                proj_name = "Mercator"
            print(f"Using projection from file: {proj_attr}")
        else:
            projection = ccrs.Mercator()
            proj_name = "Mercator"
            print("No projection found in file. Defaulting to Mercator")

        # Find primary data variable
        data_var = next((v for v in ds.variables if v.lower() not in ['latitude', 'longitude', 'lat', 'lon']), None)
        if data_var is None:
            raise ValueError("No suitable data variable found in NetCDF file.")

        var = ds.variables[data_var]
        lats = ds.variables['lat'][:]
        lons = ds.variables['lon'][:]
        data = var[:]

        # Determine whether to use log scale
        use_log = np.all(data > 0)
        units = getattr(var, "units", "")

    fig = plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=projection)

    if is_png:
        ax.imshow(
            img_data,
            origin='upper',
            extent=[min_x, max_x, min_y, max_y],
            transform=projection
        )

        if not is_true_color:
            raw_palette = img.getpalette()
            n = len(raw_palette) // 3
            colors = np.array(raw_palette).reshape(n, 3) / 255.0
            cmap = ListedColormap(colors)

            vmin = float(meta.get('data_minimum', 0))
            vmax = float(meta.get('data_maximum', n - 1))
            sm = ScalarMappable(cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax))
            sm.set_array([])

            cbar = fig.colorbar(sm, ax=ax, orientation='horizontal',
                                pad=0.05, shrink=0.5, aspect=40)
            cbar.set_label(data_var, fontsize=10)
            print(f"PNG mode '{img.mode}': adding palette colorbar ({n} colors, "
                  f"range {vmin} - {vmax})")
        else:
            print(f"PNG mode '{img.mode}': true color, no colorbar added")

    else:
        # Plot data with automatic log scaling if appropriate
        if use_log:
            mesh = ax.pcolormesh(
                lons, lats, data,
                transform=ccrs.PlateCarree(),
                cmap='viridis',
                norm=LogNorm(vmin=np.nanmin(data[data>0]), vmax=np.nanmax(data))
            )
        else:
            mesh = ax.pcolormesh(
                lons, lats, data,
                transform=ccrs.PlateCarree(),
                cmap='viridis'
            )

        # Units and colorbar label
        cbar_label = f"{data_var} ({units})" if units else data_var
        cbar = fig.colorbar(mesh, ax=ax, orientation='horizontal', pad=0.05, aspect=40)
        cbar.set_label(cbar_label, fontsize=10)

    # Add map features
    ax.coastlines(resolution='50m', linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='white', linewidth=0.5)
    if not is_png:
        ax.add_feature(cfeature.LAND, facecolor='lightgray')
        ax.add_feature(cfeature.OCEAN, facecolor='lightblue')

    # Gridlines with labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='white', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10}
    gl.ylabel_style = {'size': 10}

    # Axis labels
    ax.set_xlabel("Longitude", fontsize=12)
    ax.set_ylabel("Latitude", fontsize=12)

    # Descriptive title
    title = f"{data_var} Map ({proj_name})"
    if not is_png and units:
        title += f" — {units}"
    if label_text:
        title += f" — {label_text}"
    plt.title(title, fontsize=14, fontweight='bold')

    # Save output
    if output_file is None:
        output_file = input_path.with_name(input_path.stem + f'_{proj_name}').with_suffix('png')
    else:
        output_file = Path(output_file)

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved annotated map to: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate L3 NetCDF or PNG map output with graticules, coastlines, and optional labels")
    parser.add_argument("-i", "--ifile", required=True, help="Input L3 NetCDF or PNG file")
    parser.add_argument("-o", "--ofile", help="Output image file")
    parser.add_argument("--label", help="Optional label to annotate map title")

    args = parser.parse_args()
    label_text = args.label if args.label else Path(args.ifile).stem
    annotate_l3mapgen(args.ifile, args.ofile, label_text)