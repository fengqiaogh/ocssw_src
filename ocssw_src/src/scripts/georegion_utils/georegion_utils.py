import sys
from shapely import wkt
from shapely.geometry import Point
import pyproj
from shapely.ops import transform
from shapely.geometry import LineString
import numpy as np
import rasterio
from rasterio.transform import xy
import geopandas
import numpy as np
import magic
def check_file_type(filepath: str):
    mime = magic.Magic(mime=True)
    file_type = mime.from_file(filepath)

    if file_type == "text/plain" or file_type == "text/csv":
        return "ASCII text"
    elif file_type == "application/x-qgis" or "application/octet-stream" in file_type:
        # Shapefiles often show as octet-stream, need additional check
        if filepath.endswith(".shp"):
            return "Shapefile"
    elif file_type == "application/x-hdf5" or file_type == "application/x-netcdf":
        return "NetCDF"
    return file_type


def determine_projection(wkt_polygon):
    polygon = wkt.loads(wkt_polygon)
    proj_type = "EPSG:3031"
    project_transform = pyproj.Transformer.from_crs(
        "EPSG:4326", proj_type, always_xy=True
    ).transform
    polygon_projected = transform(project_transform, polygon)
    south_pole = Point(0, -90)
    south_pole_proj = transform(project_transform, south_pole)
    contains_south_pole = polygon_projected.contains(south_pole_proj)
    if contains_south_pole:
        area_south = polygon_projected.area
        proj_type = "EPSG:3995"
        project_transform = pyproj.Transformer.from_crs(
            "EPSG:4326", proj_type, always_xy=True
        ).transform
        polygon_projected = transform(project_transform, polygon)
        area_north = polygon_projected.area
        proj_type = "EPSG:3395"
        project_transform = pyproj.Transformer.from_crs(
            "EPSG:4326", proj_type, always_xy=True
        ).transform
        polygon_projected = transform(project_transform, polygon)
        equator = LineString([(-180, 0), (180, 0)])
        equator_proj = transform(project_transform, equator)
        crosses_equator = polygon_projected.intersects(equator_proj)
        if crosses_equator:
            print(
                "-E-: shape {} that cross both equator and poles are not supported".format(
                    wkt_polygon
                )
            )
            sys.exit(1)
        if area_south < area_north:
            return "EPSG:3031"
        else:
            return "EPSG:3995"
    return "EPSG:3395"


def get_bounds(wkt_polygon, proj_type, resolve=10000):
    polygon = wkt.loads(wkt_polygon)
    project_transform = pyproj.Transformer.from_crs(
        "EPSG:4326", proj_type, always_xy=True
    ).transform
    project_invert_transform = pyproj.Transformer.from_crs(
        proj_type, "EPSG:4326", always_xy=True
    ).transform
    polygon_projected = transform(project_transform, polygon)
    max_y = polygon_projected.bounds[3]
    min_y = polygon_projected.bounds[1]
    max_x = polygon_projected.bounds[2]
    min_x = polygon_projected.bounds[0]
    scale = resolve  # meter
    # check if scale at least five smaller than bounding box
    if (max_x - min_x) / scale < 5:
        scale = (max_x - min_x) / 5
        print(f"  -W-: adjusted scale to {scale} for x dimension")
    if (max_y - min_y) / scale < 5:
        scale = (max_y - min_y) / 5
        print(f"  -W-: adjusted scale to {scale} for y dimension")
    xsize = int((max_x - min_x) / scale) + 1
    ysize = int((max_y - min_y) / scale) + 1
    trs = rasterio.transform.from_origin(min_x, max_y, scale, scale)
    img_shape = (ysize, xsize)
    rows, cols = np.indices(img_shape)

    # Convert pixel indices to polar projection coordinates
    x_coords, y_coords = xy(trs, rows, cols)

    lon, lat = project_invert_transform(x_coords.flatten(), y_coords.flatten())
    lon = np.reshape(lon, img_shape)
    lat = np.reshape(lat, img_shape)
    south_pole = Point(0, -90)
    north_pole = Point(0, 90)
    south_pole_proj = transform(project_transform, south_pole)
    north_pole_proj = transform(project_transform, north_pole)
    R = 6378137.0  # Earth's radius in meters

    delta_lon = np.min(scale / 2 / (R * np.cos(np.pi * lat / 180))) * 180 / np.pi
    delta_lat = scale / 2 / R * 180 / np.pi
    lat_max = float(lat.max()) + float(
        delta_lat
    )  # add small buffer to avoid edge issues
    lat_min = float(lat.min()) - float(delta_lat)
    lon_max = float(lon.max()) + float(delta_lon)
    lon_min = float(lon.min()) - float(delta_lon)
    contains_south_pole = polygon_projected.contains(south_pole_proj)
    contains_north_pole = polygon_projected.contains(north_pole_proj)
    if proj_type == "EPSG:3031" and contains_south_pole:
        lat_min = -90
        lon_max = 180
        lon_min = -180
    elif proj_type == "EPSG:3995" and contains_north_pole:
        lat_max = 90
        lon_max = 180
        lon_min = -180
    lon_max = min(180, lon_max)
    lon_min = max(-180, lon_min)
    lat_max = min(90, lat_max)
    lat_min = max(-90, lat_min)
    return lon_max, lon_min, lat_max, lat_min


def read_wkt_file(path_wkt: str):
    with open(path_wkt) as f:
        wkts = f.read().splitlines()
    return geopandas.GeoSeries.from_wkt(wkts)