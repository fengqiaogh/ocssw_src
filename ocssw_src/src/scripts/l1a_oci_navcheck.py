#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset as NetCDF

__version__ = "1.0.0 2026-03-16"

# This script checks if an L1A OCI file has good navigation coverage based on the number of attitude records compared to
# the expected size of the granule in seconds.
#
# Run: l1a_oci_navcheck ifile
# Return Code:
#   0:  Good navigation coverage (enough attitude records)
#   1:  Bad file because it is missing what L1A files should have (e.g., missing att records, missing scans, etc.)
#  110: Regenerate L1A file using latest l1agen_oci

def main():
    print(f"l1a_oci_navcheck version {__version__}")

    parser = argparse.ArgumentParser(description="Check if L1A OCI file has good navigation coverage.")
    parser.add_argument("ifile", type=str, help="Path to the input file")

    args = parser.parse_args()
    fileName = args.ifile

    # open the file and grab att_records and number_of_scans and compute 
    try:     
        ncFile = NetCDF(fileName, "r", "NETCDF4")
        attRecords = len(ncFile.dimensions["att_records"])
        numScans = len(ncFile.dimensions["number_of_scans"])
        granuleSizeInSeconds = numScans / 5.7  # OCI's scan rate is 5.7 scans per second

        status = 0 if attRecords >= (granuleSizeInSeconds) else 110
        if status == 0:
            print(f"Good navigation coverage. Attitude records: {attRecords}, Granule Time: {granuleSizeInSeconds}")
        else:
            print(f"Need Regeneration for {fileName}. Attitude records: {attRecords}, Granule Time: {granuleSizeInSeconds}")
        return status

    except Exception as e:
        print(f"Error reading file {fileName}: {e}")
        return 1

if __name__ == "__main__":
    exit(main())