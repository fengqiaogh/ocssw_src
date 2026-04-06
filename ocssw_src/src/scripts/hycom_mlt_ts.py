#!/usr/bin/env python3
# coding: utf-8

# ## HYCOM Mixed Layer Thickness Calculator 
# 
# Point of Contact: Ivona Cetinic & Emerson Sirk (ivona.cetinic@nasa.gov, emerson.a.sirk@nasa.gov)
# Institution: NASA GSFC Ocean Ecology Lab
# 
# This program:
# 
# * Calculates Mixed Layer Thickness (MLT) using GSW library
# * MLT is calculated using the 0.03 kg/m³ criterium, surface density is calculated for top 10 m 
# * MLT  defaults to being capped at 1000 m max, to mimic the initial approach that OSU had
# * Outputs streamlined MLT results without intermediate variables
# * Missing values are -32767.0
# 
# Main function is asking for a salinity file, temperature file, output file, and latitude and longitude range.  
# If the whole lat/lon resolution of the file needs to be run, then last two parameters should stay empty.  
# For example, for north pacific one would call the function as:  
#     *success = mlt_gen_flexible(temp_ifile, sal_ifile, mlt_ofile, lat_range = (0, 65), lon_range = (120, 260))*  
# but to run the whole file, no spatial subseting is:  
#     *success = mlt_gen_flexible(temp_ifile, sal_ifile, mlt_ofile)*
# 
# ### Input files
# #### Temperature should be .nc hycom file, as given example (test_t.nc)
# * Dimensions: ['time', 'depth', 'lat', 'lon']
# * Variables: ['water_temp', 'time', 'depth', 'lat', 'lon']
# 
# #### Salinity should be .nc hycom file, as given example (test_s.nc) 
# * Dimensions: ['time', 'depth', 'lat', 'lon']
# * Variables: ['salinity', 'time', 'depth', 'lat', 'lon']
# 
# 
# ### Output file: ofile, mimics original mixed layer thickness .nc file 
# Output file structure:
# * Dimensions: ['time', 'lat', 'lon']
# * Variables: ['time', 'lat', 'lon', 'mixed_layer_thickness']

import os
import sys
import warnings

def check_required_modules():
    """Check if all required modules are available before proceeding"""
    required_modules = {
        'numpy': 'numpy',
        'netCDF4': 'netCDF4', 
        'gsw': 'gsw',
        'argparse': 'argparse',
        'xarray': 'xarray'
    }
    
    missing_modules = []
    
    for module_name, import_name in required_modules.items():
        try:
            __import__(import_name)
        except ImportError:
            missing_modules.append(module_name)
    
    if missing_modules:
        print("\n" + "="*60)
        print("  HYCOM Mixed Layer Thickness Calculator v1.0")
        print("  NASA GSFC Ocean Ecology Lab")
        print("="*60)
        print("\nERROR: Missing required Python modules:")
        for module in missing_modules:
            print(f"  - {module}")

        print("\nTo install missing modules, run:")
        if missing_modules:
            modules_str = ' '.join(missing_modules)
            print(f"  pip install {modules_str}")

        print("\nNote: GSW (Gibbs SeaWater) is essential for oceanographic calculations.")
        print()
        sys.exit(1)

    return True

# Check modules first, before other imports
check_required_modules()

# Now import the modules
import argparse
from datetime import date, datetime, timedelta
import numpy as np
import xarray as xr
from netCDF4 import Dataset
import gsw
warnings.filterwarnings("ignore")

def mlt_gen_flexible(temp_ifile, sal_ifile, ofile, lat_range=None, lon_range=None, debug=False, max_depth=1000):
    try:
        temp_data = Dataset(temp_ifile, "r")
        sal_data = Dataset(sal_ifile, "r")

        print(f" Processing Temperature: {os.path.basename(temp_ifile)}")
        print(f" Processing Salinity: {os.path.basename(sal_ifile)}")
        if debug:
            # Additional debugging information  # <- This line needs proper indentation
            print(f" Full Temperature file path: {temp_ifile}")
            print(f" Full Salinity file path: {sal_ifile}")
            print(f" Temperature file size: {os.path.getsize(temp_ifile) / (1024*1024):.2f} MB")
            print(f" Salinity file size: {os.path.getsize(sal_ifile) / (1024*1024):.2f} MB")
            print(f" Maximum depth limit: {max_depth}m")
            
        if lat_range:
            print(f" Latitude subset: {lat_range[0]} to {lat_range[1]}")
        if lon_range:
            print(f" Longitude subset: {lon_range[0]} to {lon_range[1]}")

        # Step 1: Find coordinate variables (use temperature file as reference)
        coord_mapping = {}
        input_data = temp_data  # Use temp file for coordinates

        # Find latitude
        lat_candidates = [v for v in input_data.variables if "lat" in v.lower()]
        if lat_candidates:
            coord_mapping["lat"] = lat_candidates[0]
            lat_data = input_data.variables[coord_mapping["lat"]]
            print(f" Latitude: {coord_mapping['lat']} {lat_data.shape}")
        else:
            raise ValueError("No latitude variable found")

        # Find longitude
        lon_candidates = [v for v in input_data.variables if "lon" in v.lower()]
        if lon_candidates:
            coord_mapping["lon"] = lon_candidates[0]
            lon_data = input_data.variables[coord_mapping["lon"]]
            print(f" Longitude: {coord_mapping['lon']} {lon_data.shape}")
        else:
            raise ValueError("No longitude variable found")

        # Find depth
        depth_candidates = [
            v
            for v in input_data.variables
            if any(d in v.lower() for d in ["depth", "lev", "z"])
        ]
        if depth_candidates:
            coord_mapping["depth"] = depth_candidates[0]
            depth_data = input_data.variables[coord_mapping["depth"]]
            print(f" Depth: {coord_mapping['depth']} {depth_data.shape}")
        else:
            raise ValueError("No depth variable found")

        # Find time (optional) - PRESERVE TIME INFO
        time_value = None
        time_units = "days since 1900-01-01 00:00:00"
        time_candidates = [v for v in input_data.variables if "time" in v.lower()]
        if time_candidates:
            coord_mapping["time"] = time_candidates[0]
            time_data = input_data.variables[coord_mapping["time"]]
            print(f" Time: {coord_mapping['time']} {time_data.shape}")
            # Extract time value and units for output
            if hasattr(time_data, 'units'):
                time_units = time_data.units
            time_value = time_data[0] if len(time_data.shape) > 0 else time_data[:]

        # Extract coordinate arrays
        lat_array = lat_data[:]
        lon_array = lon_data[:]
        depth_array = depth_data[:]

        # Ensure 1D coordinate arrays
        if lat_array.ndim > 1:
            lat_array = lat_array.flatten()
        if lon_array.ndim > 1:
            lon_array = lon_array.flatten()
        if depth_array.ndim > 1:
            depth_array = depth_array.flatten()

        # Apply lat/lon subsetting if specified
        lat_indices = None
        lon_indices = None
        
        if lat_range:
            lat_min, lat_max = lat_range
            lat_indices = np.where((lat_array >= lat_min) & (lat_array <= lat_max))[0]
            if len(lat_indices) == 0:
                raise ValueError(f"No latitude points found in range {lat_range}")
            lat_array = lat_array[lat_indices]
            print(f" Latitude subset: {len(lat_indices)} points from {lat_min} to {lat_max}")
            
        if lon_range:
            lon_min, lon_max = lon_range
            lon_indices = np.where((lon_array >= lon_min) & (lon_array <= lon_max))[0]
            if len(lon_indices) == 0:
                raise ValueError(f"No longitude points found in range {lon_range}")
            lon_array = lon_array[lon_indices]
            print(f" Longitude subset: {len(lon_indices)} points from {lon_min} to {lon_max}")

        print(f" Coordinate arrays after subsetting:")
        print(
            f"   Lat: {len(lat_array)} points ({lat_array.min():.2f} to {lat_array.max():.2f})"
        )
        print(
            f"   Lon: {len(lon_array)} points ({lon_array.min():.2f} to {lon_array.max():.2f})"
        )
        print(
            f"   Depth: {len(depth_array)} levels ({depth_array.min():.1f} to {depth_array.max():.1f}m)"
        )

        if debug and lat_range:
            print(f" Latitude range covers {lat_range[1] - lat_range[0]:.1f} degrees")
        if debug and lon_range:
            print(f" Longitude range covers {lon_range[1] - lon_range[0]:.1f} degrees")
    
        temp_var = None
        for var_name, var in temp_data.variables.items():
            if len(var.shape) >= 3:  # Must be at least 3D
                var_lower = var_name.lower()
                long_name = getattr(var, "long_name", "").lower()
                standard_name = getattr(var, "standard_name", "").lower()

                # Check for temperature
                if (
                    any(
                        temp_word in var_lower + long_name + standard_name
                        for temp_word in ["temp", "temperature"]
                    )
                    and not temp_var
                ):
                    temp_var = var_name
                    break

        # Step 3: Find salinity variable in salinity file
        sal_var = None
        for var_name, var in sal_data.variables.items():
            if len(var.shape) >= 3:  # Must be at least 3D
                var_lower = var_name.lower()
                long_name = getattr(var, "long_name", "").lower()
                standard_name = getattr(var, "standard_name", "").lower()

                # Check for salinity
                if (
                    any(
                        sal_word in var_lower + long_name + standard_name
                        for sal_word in ["sal", "salinity"]
                    )
                    and not sal_var
                ):
                    sal_var = var_name
                    break

        if not temp_var or not sal_var:
            print(
                f" Missing variables - Temperature: {temp_var}, Salinity: {sal_var}"
            )
            return False

        temp_var_data = temp_data.variables[temp_var]
        sal_var_data = sal_data.variables[sal_var]

        print(f" Temperature: {temp_var} {temp_var_data.shape} {temp_var_data.dimensions}")
        print(f" Salinity: {sal_var} {sal_var_data.shape} {sal_var_data.dimensions}")

        # Function to extract and subset 3D array from either file
        def extract_and_subset_3d_array(data_var, var_name, file_type):
            """Extract 3D array from data variable and apply subsetting"""
            shape = data_var.shape
            
            if len(shape) == 4:  # Likely (time, depth, lat, lon) or similar
                dims = data_var.dimensions
                print(f"   {file_type} dimensions: {dims}")

                # Find axis positions
                time_axis = None
                depth_axis = None
                lat_axis = None
                lon_axis = None

                for i, dim in enumerate(dims):
                    if "time" in dim.lower():
                        time_axis = i
                    elif any(d in dim.lower() for d in ["depth", "lev", "z"]):
                        depth_axis = i
                    elif "lat" in dim.lower():
                        lat_axis = i
                    elif "lon" in dim.lower():
                        lon_axis = i

                print(f"   {file_type} axes - Time: {time_axis}, Depth: {depth_axis}, Lat: {lat_axis}, Lon: {lon_axis}")

                # Create slicing object based on dimensions
                slice_obj = [slice(None)] * len(shape)
                
                if time_axis is not None:
                    slice_obj[time_axis] = 0  # First time step
                
                if lat_axis is not None and lat_indices is not None:
                    slice_obj[lat_axis] = lat_indices
                    
                if lon_axis is not None and lon_indices is not None:
                    slice_obj[lon_axis] = lon_indices
                
                # Extract and subset in one operation
                array_3d = data_var[tuple(slice_obj)]

            elif len(shape) == 3:  # Already 3D, need to determine dimension order
                dims = data_var.dimensions
                print(f"   {file_type} dimensions (3D): {dims}")
                
                # Find axis positions
                depth_axis = None
                lat_axis = None
                lon_axis = None

                for i, dim in enumerate(dims):
                    if any(d in dim.lower() for d in ["depth", "lev", "z"]):
                        depth_axis = i
                    elif "lat" in dim.lower():
                        lat_axis = i
                    elif "lon" in dim.lower():
                        lon_axis = i
                
                # Create slicing object based on dimensions
                slice_obj = [slice(None)] * 3
                
                if lat_axis is not None and lat_indices is not None:
                    slice_obj[lat_axis] = lat_indices
                    
                if lon_axis is not None and lon_indices is not None:
                    slice_obj[lon_axis] = lon_indices
                
                # Extract and subset in one operation
                array_3d = data_var[tuple(slice_obj)]
            else:
                raise ValueError(f"Unexpected {file_type} data shape: {shape}")

            return array_3d

        # Extract 3D arrays from both files with subsetting
        temp_array = extract_and_subset_3d_array(temp_var_data, temp_var, "Temperature")
        sal_array = extract_and_subset_3d_array(sal_var_data, sal_var, "Salinity")

        print(f"   Extracted arrays: temp {temp_array.shape}, sal {sal_array.shape}")
    
        if debug:
            print(f"   Memory usage for arrays: ~{temp_array.nbytes/(1024*1024):.1f} MB each")

        # Ensure correct order (depth, lat, lon) for both arrays
        expected_shape = (len(depth_array), len(lat_array), len(lon_array))
        print(f"   Expected shape: {expected_shape}")

        def reorder_array(array, array_name, expected_shape):
            """Reorder array dimensions to match expected shape"""
            if array.shape != expected_shape:
                print(f" {array_name} shape mismatch, attempting to transpose...")
                
                if array.shape == (len(lat_array), len(lon_array), len(depth_array)):
                    array = array.transpose(2, 0, 1)  # (lat,lon,depth) -> (depth,lat,lon)
                    print(f"   {array_name} transposed to: {array.shape}")
                elif array.shape == (len(depth_array), len(lon_array), len(lat_array)):
                    array = array.transpose(0, 2, 1)  # (depth,lon,lat) -> (depth,lat,lon)
                    print(f"   {array_name} transposed to: {array.shape}")
                elif array.shape[0] == len(depth_array) and array.shape[1]*array.shape[2] == len(lat_array)*len(lon_array):
                    # Special case for when lat/lon indices have been applied but dimensions need reordering
                    print(f" {array_name} dimensions need special handling after subsetting")
                    # This requires dataset-specific handling
            
            return array

        temp_array = reorder_array(temp_array, "Temperature", expected_shape)
        sal_array = reorder_array(sal_array, "Salinity", expected_shape)

        # Verify both arrays have the same final shape
        if temp_array.shape != sal_array.shape:
            print(f" Warning: Temperature shape {temp_array.shape} and salinity shape {sal_array.shape} differ")
            # Try to reconcile shapes if they're close but not identical
            # This might be needed if subsetting worked slightly differently for each file

        # Limit to mmax_depth
        depth_mask = depth_array <= max_depth
        if np.any(~depth_mask):
            n_depths_original = len(depth_array)
            depth_array = depth_array[depth_mask]
            temp_array = temp_array[depth_mask, :, :]
            sal_array = sal_array[depth_mask, :, :]
            print(
                f"\n Limited to top {max_depth}m: {len(depth_array)} levels (was {n_depths_original})"
            )

        # Step 8: Calculate MLT using the corrected arrays
        print("\n Calculating pressure and density...")

        def apply_enhanced_qc(temp_array, sal_array, depth_array):
            """Apply enhanced oceanographic quality control checks"""
            original_temp_valid = np.sum(~np.isnan(temp_array))
            original_sal_valid = np.sum(~np.isnan(sal_array))

            final_temp_valid = np.sum(~np.isnan(temp_array))
            final_sal_valid = np.sum(~np.isnan(sal_array))

            print(
                f"   Temperature: {original_temp_valid}  {final_temp_valid} valid points"
            )
            print(f"   Salinity: {original_sal_valid} {final_sal_valid} valid points")

            return temp_array, sal_array

        # Apply enhanced quality control
        temp_array, sal_array = apply_enhanced_qc(temp_array, sal_array, depth_array)

        # Calculate pressure
        pressure_array = np.zeros_like(temp_array)

        for i, depth_val in enumerate(depth_array):
            if i % max(1, len(depth_array) // 10) == 0:
                progress = (i + 1) / len(depth_array) * 100
                print(f"   Progress: {progress:.0f}% (depth {depth_val:.1f}m)")

            for j, lat_val in enumerate(lat_array):
                pressure_array[i, j, :] = gsw.p_from_z(-depth_val, lat_val)

        # Continue with GSW calculations...
        print("   Converting to GSW variables...")

        lon_grid, lat_grid = np.meshgrid(lon_array, lat_array)

        # Process in chunks to save memory
        chunk_size = min(10, len(depth_array))
        density_array = np.zeros_like(temp_array)

        for start_idx in range(0, len(depth_array), chunk_size):
            end_idx = min(start_idx + chunk_size, len(depth_array))

            temp_chunk = temp_array[start_idx:end_idx, :, :]
            sal_chunk = sal_array[start_idx:end_idx, :, :]
            pres_chunk = pressure_array[start_idx:end_idx, :, :]

            abs_sal_chunk = np.zeros_like(sal_chunk)
            cons_temp_chunk = np.zeros_like(temp_chunk)

            for k in range(end_idx - start_idx):
                abs_sal_chunk[k, :, :] = gsw.SA_from_SP(
                    sal_chunk[k, :, :], pres_chunk[k, :, :], lon_grid, lat_grid
                )
                cons_temp_chunk[k, :, :] = gsw.CT_from_t(
                    abs_sal_chunk[k, :, :], temp_chunk[k, :, :], pres_chunk[k, :, :]
                )

            density_array[start_idx:end_idx, :, :] = gsw.rho(
                abs_sal_chunk, cons_temp_chunk, 0
            )

        # ENHANCED MLT CALCULATION WITH INTERPOLATION PHILOSOPHY
        print("\n Calculating Mixed Layer Thickness with interpolation...")

        def calculate_robust_surface_reference(
            density_array, depth_array, max_surface_depth=10
        ):
            """Calculate surface reference density with better handling"""
            surface_mask = depth_array <= max_surface_depth
            surface_indices = np.where(surface_mask)[0]

            if len(surface_indices) == 0:
                # If no data within 10m, use shallowest available
                surface_indices = [0]
                print(f"\n No data <= {max_surface_depth}m, using shallowest level")

            # Use maximum instead of mean/median for surface density reference
            surface_density = np.nanmean(density_array[surface_indices, :, :], axis=0)

            # Additional QC on surface reference
            surface_density[(surface_density < 990) | (surface_density > 1060)] = (
                np.nan
            )

            return surface_density

        def calculate_mlt_with_interpolation(
            density_array, depth_array, surface_density, threshold=0.03
        ):
            """Calculate MLT using interpolation philosophy similar to reference code"""
            nlat, nlon = surface_density.shape
            mlt_result = np.full((nlat, nlon), -32767.0)

            # Create threshold array
            threshold_density = surface_density + threshold

            print(f"   Processing {nlat}x{nlon} grid points...")

            # Process each water column
            for j in range(nlat):
                if j % max(1, nlat // 10) == 0:
                    progress = j / nlat * 100
                    print(f"   Grid progress: {progress:.0f}% (row {j}/{nlat})")

                for i in range(nlon):
                    if not np.isnan(threshold_density[j, i]):
                        # Extract this profile
                        profile = density_array[:, j, i]

                        # Remove NaNs and corresponding depths
                        valid_mask = ~np.isnan(profile)
                        if np.sum(valid_mask) < 3:  # Need at least 3 points
                            continue

                        clean_depth = depth_array[valid_mask]
                        clean_density = profile[valid_mask]

                        # Find depths > 10m for MLT search 
                        search_mask = clean_depth > 10.0
                        if not np.any(search_mask):
                            continue

                        search_depths = clean_depth[search_mask]
                        search_density = clean_density[search_mask]

                        # DENSE INTERPOLATION 
                        # Create fine depth grid with 2m resolution
                        depth_fine = np.arange(
                            search_depths.min(), min(search_depths.max(), max_depth) + 2, 2.0
                        )

                        if len(depth_fine) < 2:
                            continue

                        # Interpolate density to fine grid
                        density_fine = np.interp(
                            depth_fine, search_depths, search_density
                        )

                        # Find threshold crossing 
                        exceeds_idx = np.where(density_fine > threshold_density[j, i])[
                            0
                        ]

                        if len(exceeds_idx) > 0:
                            # Found the MLT - use interpolated value
                            mlt_result[j, i] = depth_fine[exceeds_idx[0]]
                        else:
                            # Well-mixed to search limit (capped at max_depth)
                            mlt_result[j, i] = min(search_depths.max(), max_depth)

            return mlt_result

        # Calculate robust surface reference
        surface_density = calculate_robust_surface_reference(density_array, depth_array)

        # Calculate MLT using interpolation philosophy
        mlt_result = calculate_mlt_with_interpolation(
            density_array, depth_array, surface_density
        )

        # Apply final cap of Max_depth
        mlt_result[mlt_result > max_depth] = max_depth

        # Save output with REQUIRED STRUCTURE
        output = Dataset(ofile, mode="w")

        # Add more detailed metadata
        output.title = "Mixed Layer Thickness from separate HYCOM Temperature and Salinity Data (Enhanced with Interpolation)"
        output.source_files = f"Temperature: {os.path.basename(temp_ifile)}, Salinity: {os.path.basename(sal_ifile)}"
        output.date_created = date.today().strftime("%Y-%m-%d")
        output.method = "Density threshold (0.03 kg/m³) with 2m interpolation"
        output.max_depth_limit = f"{max_depth}m" 
        if lat_range:
            output.latitude_range = f"{lat_range[0]} to {lat_range[1]}"
        if lon_range:
            output.longitude_range = f"{lon_range[0]} to {lon_range[1]}"

        # CREATE REQUIRED DIMENSIONS: time, lat, lon
        output.createDimension("time", 1)  # Single time step
        output.createDimension("lat", len(lat_array))
        output.createDimension("lon", len(lon_array))

        # CREATE COORDINATE VARIABLES
        time_out = output.createVariable("time", "f8", ("time",))
        lat_out = output.createVariable("lat", "f8", ("lat",))
        lon_out = output.createVariable("lon", "f8", ("lon",))
        
        # Set coordinate data
        if time_value is not None:
            time_out[:] = [time_value]
            time_out.units = time_units
        else:
            # Create a default time value if none found
            time_out[:] = [0]
            time_out.units = "days since 1900-01-01 00:00:00"
        
        time_out.long_name = "time"
        lat_out[:] = lat_array
        lat_out.units = "degrees_north"
        lat_out.long_name = "latitude"
        lon_out[:] = lon_array  
        lon_out.units = "degrees_east"
        lon_out.long_name = "longitude"

        # CREATE REQUIRED VARIABLES WITH TIME DIMENSION
        mixed_layer_var = output.createVariable(
            "mixed_layer_thickness", "f4", ("time", "lat", "lon"), fill_value=-32767.0
        )
        
        # Set the MLT data (add time dimension by expanding)
        mlt_data_3d = np.expand_dims(mlt_result, axis=0)  # Add time dimension
        mixed_layer_var[:] = mlt_data_3d 
        
        # Set attributes
        mixed_layer_var.units = "meters"
        mixed_layer_var.long_name = f"Mixed Layer Thickness (interpolated, capped at {max_depth}m)"
        mixed_layer_var.method = "Density threshold with 2m interpolation"

        # Close all files
        temp_data.close()
        sal_data.close()
        output.close()

        valid_mlt = mlt_result[mlt_result != -32767.0]
        print(f"\n Success! Valid MLT points: {len(valid_mlt):,}")
        if len(valid_mlt) > 0:
            print(f"   Range: {np.min(valid_mlt):.1f} - {np.max(valid_mlt):.1f} m")
            print(f"   Mean: {np.mean(valid_mlt):.1f} m")
            print(f"   Median: {np.median(valid_mlt):.1f} m")
            capped_count = np.sum(valid_mlt == max_depth)
            if capped_count > 0:
                print(
                    f"   Capped at {max_depth}m: {capped_count:,} points ({capped_count/len(valid_mlt)*100:.1f}%)"
                )
        print(f" Output: {ofile}")
        print(" Structure: Dimensions=['time', 'lat', 'lon'], Variables=['mixed_layer_thickness', 'time', 'lat', 'lon']")

        return True

    except Exception as e:
        print(f" Error: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
# Print program name and version
    print(f"{os.path.basename(sys.argv[0])} Ver 1.0\n")

    parser = argparse.ArgumentParser(
        description='Calculate Mixed Layer Thickness from HYCOM Temperature and Salinity data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Examples:
          # Process entire domain:
            python mlt_calculator.py --temp-ifile temp.nc --sal-ifile sal.nc --mlt-ofile output_mlt.nc
  
          # Process North Pacific subset:
            python mlt_calculator.py --temp-ifile temp.nc --sal-ifile sal.nc --mlt-ofile output_mlt.nc --lat-range 0 65 --lon-range 120 260
  
          # Process Atlantic subset:
            python mlt_calculator.py --temp-ifile temp.nc --sal-ifile sal.nc --mlt-ofile atlantic_mlt.nc --lat-range 0 70 --lon-range -80 20
        """
    )
    
    # Required arguments
    parser.add_argument('--temp-ifile', required=True, help='Input temperature NetCDF file path')
    parser.add_argument('--sal-ifile', required=True, help='Input salinity NetCDF file path')
    parser.add_argument('--mlt-ofile', required=True, help='Output MLT NetCDF file path')
    
    # Optional spatial subsetting
    parser.add_argument('--lat-range', nargs=2, type=float, metavar=('MIN', 'MAX'),
                       help='Latitude range for spatial subsetting (degrees North)')
    parser.add_argument('--lon-range', nargs=2, type=float, metavar=('MIN', 'MAX'),
                       help='Longitude range for spatial subsetting (degrees East)')
    # Optional depth limit
    parser.add_argument('--max-depth', type=float, default=1000, metavar='DEPTH',
                       help='Maximum depth limit in meters (default: 1000)')
    # Optional flags
    parser.add_argument('--debug', '-d', action='store_true', 
                   help='Enable debugging output')
    parser.add_argument('--version','-v', action='version', version='HYCOM MLT Calculator v1.0 NASA GSFC Ocean Ecology Lab')
    parser.add_argument('--logfile','-l', help='Write log output to file', metavar='LOGFILE')

    args = parser.parse_args()
    # check if we are writing an optional log
    if args.logfile:
        import logging
        logging.basicConfig(
            filename=args.logfile, 
            level=logging.INFO,
            format='%(asctime)s - %(message)s',
            filemode='w'
        ) 

    # print Python module search paths if debug mode
    if args.debug:
        print("\n Python module search paths (sys.path):")
        for i, path in enumerate(sys.path):
            print(f" {i+1}. {path}")
        print()
    # Validate input files exist
    for file_path in [args.temp_ifile, args.sal_ifile]:
        if not os.path.exists(file_path):
            print(f" Error: Input file does not exist: {file_path}")
            sys.exit(1)
    
    # Check output directory exists
    output_dir = os.path.dirname(args.mlt_ofile)
    if output_dir and not os.path.exists(output_dir):
        print(f" Error: Output directory does not exist: {output_dir}")
        sys.exit(1)
    
    # Validate latitude range
    if args.lat_range:
        if not (-90 <= args.lat_range[0] <= 90 and -90 <= args.lat_range[1] <= 90):
            print(" Error: Latitude values must be between -90 and 90")
            sys.exit(1)
        if args.lat_range[0] >= args.lat_range[1]:
            print(" Error: Latitude minimum must be less than maximum")
            sys.exit(1)
    
    # Validate longitude range  
    if args.lon_range:
        if args.lon_range[0] >= args.lon_range[1]:
            print(" Error: Longitude minimum must be less than maximum")
            sys.exit(1)
    
    # Convert to tuples if provided
    lat_range = tuple(args.lat_range) if args.lat_range else None
    lon_range = tuple(args.lon_range) if args.lon_range else None
    
    # Validate max depth
    if args.max_depth <= 0:
        print(" Error: Maximum depth must be positive")
        sys.exit(1)
    if args.max_depth > 2000:  # Reasonable upper limit
        print(" Warning: Maximum depth is very large (>2km)")

    # Print processing info
    print(" HYCOM Mixed Layer Thickness Calculator")
    print(f" Temperature file: {args.temp_ifile}")
    print(f" Salinity file: {args.sal_ifile}")
    print(f" Output file: {args.mlt_ofile}")
    print(f" Maximum depth limit: {args.max_depth}m")
    if lat_range:
        print(f" Latitude range: {lat_range[0]}° to {lat_range[1]}°N")
    if lon_range:
        print(f" Longitude range: {lon_range[0]}° to {lon_range[1]}°E")
    print()
    
# Call mlt_gen_flexible
    success = mlt_gen_flexible(args.temp_ifile, args.sal_ifile, args.mlt_ofile,
                          lat_range=lat_range, lon_range=lon_range, 
                          debug=args.debug, max_depth=args.max_depth)
    
    if success:
        print("\n MLT calculation completed successfully!")
        sys.exit(0)
    else:
        print("\n MLT calculation failed!")
        sys.exit(1)
if __name__ == "__main__":
    main()

