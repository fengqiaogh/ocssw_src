#! /usr/bin/env python3
import argparse
from pathlib import Path
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LogNorm
import os

def annotate_l3mapgen(input_file, output_file=None, label_text=None):
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

    # Set path for Cartopy to find its downloaded shapefiles
    # Make sure the path is expanded *before* Cartopy tries to load anything
    cartopy_data_root = os.environ['OCDATAROOT'] + '/common/cartopy'
    cartopy.config['data_dir'] = os.path.expanduser(cartopy_data_root)
    cartopy.config['pre_existing_data_dir'] = os.path.expanduser(cartopy_data_root)

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

    fig = plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=projection)

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

    # Add map features
    ax.coastlines(resolution='50m', linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')

    # Gridlines with labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10}
    gl.ylabel_style = {'size': 10}

    # Axis labels
    ax.set_xlabel("Longitude", fontsize=12)
    ax.set_ylabel("Latitude", fontsize=12)

    # Units and colorbar label
    units = getattr(var, "units", "")
    cbar_label = f"{data_var} ({units})" if units else data_var

    # Colorbar
    cbar = fig.colorbar(mesh, ax=ax, orientation='horizontal', pad=0.05, aspect=40)
    cbar.set_label(cbar_label, fontsize=10)

    # Descriptive title
    title = f"{data_var} Map ({proj_name})"
    if units:
        title += f" — {units}"
    if label_text:
        title += f" — {label_text}"
    plt.title(title, fontsize=14, fontweight='bold')

    # Save output
    if output_file is None:
        output_file = Path(input_file + f'_{proj_name}').with_suffix('.png')
    else:
        output_file = Path(output_file)

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved annotated map to: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate L3 NetCDF output with graticules, coastlines, and optional labels")
    parser.add_argument("-i", "--ifile", required=True, help="Input L3 NetCDF file")
    parser.add_argument("-o", "--ofile", help="Output image file")
    parser.add_argument("--label", help="Optional label to annotate map title")

    args = parser.parse_args()
    label_text = args.label if args.label else Path(args.ifile).stem
    annotate_l3mapgen(args.ifile, args.ofile, label_text)