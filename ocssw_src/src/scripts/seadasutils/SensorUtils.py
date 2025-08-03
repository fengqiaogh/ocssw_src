import os
import json
from typing import Dict, Optional

"""
Utility functions for determining directories for each sensor.
"""

sensor_data = {}


def load_sensorlist(filename: Optional[str] = None) -> Dict[str, Dict[str, str]]:
    """
    Read list of sensor definitions from JSON file
    """
    global sensor_data
    sensor_data = {}
    try:
        if not filename:
            filename = os.path.join(
                os.getenv("OCDATAROOT"), "common", "SensorInfo.json"
            )
        with open(filename, "r") as infile:
            sensor_data = json.load(infile)

        # Ensure all required keys are present in each sensor definition
        required_keys = ["sensor", "instrument", "platform", "dir"]
        for _, sensor_info in sensor_data.items():
            for key in required_keys:
                if key not in sensor_info:
                    sensor_info[key] = ""  # Add a key if it's missing
    except Exception as e:
        print(f"Error loading sensor list: {e}")


def by_sensor(name: str) -> Optional[Dict[str, Dict[str, str]]]:
    """
    Get sensor defs from unique sensor name

    This function performs a case-insensitive search over the sensor information file for a matching sensor name.

    Args:
        name: The name of the sensor to look for

    Returns:
        A dictionary containing sensor information

    Example:
        >>> oci_info = by_sensor("oci")
        >>> oci_info
        ... {'sensor': 'OCI', 'instrument': 'OCI', 'platform': 'PACE', 'dir': 'oci'}
    """
    global sensor_data
    if len(sensor_data) == 0:
        load_sensorlist()
    try:
        return next(
            sensor_data[data]
            for data in sensor_data.keys()
            if sensor_data[data]["sensor"].lower() == name.lower()
        )
    except StopIteration:
        return None


def by_instplat(inst: str, plat: str) -> Optional[Dict[str, Dict[str, str]]]:
    """
    Retrieve sensor definitions based on the given instrument and platform.

    This function searches the global sensor_data dictionary for a matching
    instrument and platform combination. The search is case-insensitive.

    Args:
        inst: The instrument name to search for.
        plat: The platform name to search for.

    Returns:
        A dictionary containing sensor
        definitions if a match is found, None otherwise.

    Example:
        >>> inst_plat = by_instplat("MODIS", "Aqua")
        >>> inst_plat
        ... {'sensor': 'MODISA', 'instrument': 'MODIS', 'platform': 'Aqua', 'dir': 'modis/aqua'}

    Note:
        This function loads the sensor data if it hasn't been loaded already.
    """
    global sensor_data
    if len(sensor_data) == 0:
        load_sensorlist()
    try:
        return next(
            sensor_data[data]
            for data in sensor_data
            if (sensor_data[data]["instrument"].lower() == inst.lower())
            and (sensor_data[data]["platform"].lower() == plat.lower())
        )
    except StopIteration:
        return None


def by_desc(name: str) -> Optional[Dict[str, Dict[str, str]]]:
    """
    Get sensor definitions based on a descriptive string.

    Searches for sensor information using a case-insensitive string that matches
    entries in the sensor info file.

    Args:
        name: A case-insensitive string to search for sensor information.

    Returns:
        A nested dictionary where the outer keys are sensor identifiers and
        the inner dictionaries contain sensor-specific information fields.
        Returns None if no match is found.

    Example:
        >>> result = by_desc("MODIS")
        >>> if result:
        ...     print(result.keys())
        ['MODIS_AQUA', 'MODIS_TERRA']
    """
    global sensor_data
    if len(sensor_data) == 0:
        load_sensorlist()
    try:
        return next(
            sensor_data[data]
            for data in sensor_data
            if name.lower()
            in (
                sensor_data[data]["sensor"].lower(),
                sensor_data[data]["platform"].lower(),
                sensor_data[data].get("subdir"),  # Might not exist, so we use get
                sensor_data[data]["dir"].lower(),
                sensor_data[data]["instrument"].lower(),
            )
        )
    except StopIteration:
        return None


# end of class SensorUtils

# test routines below
if __name__ == "__main__":
    if len(sensor_data) == 0:

        print("\nby sensor:")
        namelist = ["modisa", "seawifs", "bogus"]
        for sensor in namelist:
            print(sensor, ":\t", by_sensor(sensor))

        print("\nby inst/plat:")
        instlist = ["modis", "VIIRS", "modis", "bogus", "Aquarius"]
        platlist = ["Terra", "JPSS-1", "bogus", "Aqua", "SAC-D"]
        for inst, plat in zip(instlist, platlist):
            print(inst, plat, ":\t", by_instplat(inst, plat))

        print("\nby any name:")
        namelist = [
            "seawifs",
            "aquarius",
            "modisa",
            "modist",
            "viirsn",
            "viirsj1",
            "aqua",
            "terra",
            "npp",
            "j1",
            "bogus",
        ]
        for name in namelist:
            print(name, ":\t", by_desc(name))
