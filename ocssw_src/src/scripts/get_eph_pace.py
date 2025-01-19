#!/usr/bin/env python3

import argparse
import os
from datetime import datetime, timezone, time
import fnmatch
from pathlib import Path
import shutil

BASE_PATH = "/glusteruser/sdsmo/PACE/data/EPH"
VERBOSE = False


def get_day_of_year(date=None, verbose=False):
    if date is None:
        date = datetime.now(timezone.utc)
    first_day = datetime(date.year, 1, 1, tzinfo=date.tzinfo)

    day_of_year = (date - first_day).days  # Yo dawg ...
    # If the time is past 1600, we have to look at the next EPH file, not the one from above
    if date.time() >= time(16, 0, 0):
        if verbose:
            print("Time past 1600, adding one to julian day")
        day_of_year += 1
    return day_of_year


def parse_datetime(datetime_string):
    try:
        return datetime.strptime(datetime_string, "%Y%m%dT%H%M%S").replace(
            tzinfo=timezone.utc
        )
    except ValueError:
        raise argparse.ArgumentTypeError(
            "Invalid datetime format. Please use YYYYMMDDTHHMMSS format."
        )


def find_files(julian_day, year, verbose=False):
    dir_path = f"{BASE_PATH}/{year}"
    pattern = f"PACE_EPH_DEF_{year}{julian_day:03d}_0*.oem"

    if verbose:
        print(f'Searching {dir_path} using pattern "{pattern}"')

    matching_files = []

    try:
        for filename in os.listdir(dir_path):
            if fnmatch.fnmatch(filename, pattern):
                full_path = os.path.join(dir_path, filename)
                matching_files.append(full_path)
                if verbose:
                    print(f"    Found matching file: {full_path}")

    except Exception as e:
        print(f"Error while searching for files: {e}")

    if not matching_files and verbose:
        print("No matching files found.")

    return matching_files


def parse_file_datetime(datetime_string):
    date_time, fractional = datetime_string.split(".")
    dt = datetime.strptime(date_time, "%Y-%m-%dT%H:%M:%S")
    microseconds = int(fractional[:6].ljust(6, "0"))
    dt = dt.replace(microsecond=microseconds, tzinfo=timezone.utc)
    return dt


def verify_file(filename, target_datetime, verbose=False):
    with open(filename, "r") as file:
        start_time = None
        stop_time = None
        for line in file:
            if "START_TIME" in line:
                start_time = line.split(" = ")[1].strip()
            elif "STOP_TIME" in line:
                stop_time = line.split(" = ")[1].strip()
            if start_time and stop_time:
                break

        if start_time and stop_time:
            if verbose:
                print(
                    f"Verifying {target_datetime} is between {start_time} and {stop_time}"
                )
            start = parse_file_datetime(start_time)
            stop = parse_file_datetime(stop_time)
            return start <= target_datetime <= stop
    return False


def copy_file(source_path, destination_path, copy_metadata, verbose):
    if verbose:
        print(f"Attempting to copy {source_path} to {destination_path}")

    try:
        if copy_metadata:
            shutil.copy2(source_path, destination_path)
        else:
            shutil.copy(source_path, destination_path)
        if verbose:
            print("Copy successful")
    except IOError as ioe:
        print("Unable to copy file; ", ioe)
    except:
        print("Unexpected error: ")


def is_leap_year(year):
    return (year % 4 == 0) and (year % 100 != 0 or year % 400 == 0)


def add_args(parser):
    parser.add_argument(
        "-d",
        "--datetime",
        type=parse_datetime,
        help="Specify a date-time in YYYYMMDDTHHMMSS format",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Display more information"
    )
    parser.add_argument(
        "-o",
        "--odir",
        help="Full path to desired output directory. \nDefaults to current working directory (%s)"
        % Path.cwd(),
        default=Path.cwd(),
    )
    parser.add_argument(
        "-c",
        "--copy_metadata",
        help="Copy the metadata of the ephemeris file",
        action="store_true",
    )


def main():
    parser = argparse.ArgumentParser(
        description="Find and download definitive ephemeris files based on Julian day"
    )
    add_args(parser)
    args = parser.parse_args()

    if args.datetime:
        date = args.datetime
        if args.verbose:
            print("Datetime specified, using ", date)
    else:
        date = datetime.now(timezone.utc)
        if args.verbose:
            print("Datetime not specified, using ", date)

    julian_day = get_day_of_year(date)
    year = date.year

    if args.verbose:
        print(f"Date-time: {date.strftime('%Y-%m-%d %H:%M:%S')} UTC")
        print(f"Julian day (0-indexed): {julian_day}")

    files = find_files(julian_day, year, args.verbose)

    if not files:
        if args.verbose:
            print(
                f"find_files() was unable to find files with day {julian_day} in year {year}"
            )
        return 1

    verified_files = [f for f in files if verify_file(f, date, args.verbose)]

    if verified_files:
        if args.verbose:
            print("File verification successful; ")
            for file in verified_files:
                print("   ", file)
        copy_file(verified_files[0], args.odir, args.copy_metadata, args.verbose)
        return 0
    else:
        """
        The definitive ephemeris files cover some 32 hours, and there are 3 overlapping files per day, offset by 8 hours.
        If the first try doesn't find at least one appropriate file, the next day might have a correct time span. This happens commonly
        when the date lines up with the final of the aforementioned 3 files, but the time is off the end of that third file.

        Note that Julian days are denoted by 12:00:00 GMT, not 00:00:00
        """
        if args.verbose:
            print("No files match the specified datetime, trying next day")

        if (
            is_leap_year(year)
            and julian_day == 366
            or not is_leap_year(year)
            and julian_day == 365
        ):
            if args.verbose:
                print("Rolling over the end of the year")
            next_day = 0
            year += 1

        next_day = julian_day + 1
        files = find_files(next_day, year, args.verbose)

        if not files:
            if args.verbose:
                print(
                    f"find_files() was unable to find files with day {next_day} in year {year}"
                )
            return 1

        verified_files = [f for f in files if verify_file(f, date, args.verbose)]

        if verified_files:
            if args.verbose:
                print("File verification successful; ")
                for file in verified_files:
                    print("   ", file)
            copy_file(verified_files[0], args.odir, args.copy_metadata, args.verbose)
            return 0
        else:
            if args.verbose:
                print("No files match the next day either")
            return 1


if __name__ == "__main__":
    if main() != 0:
        print("Unable to find definitive ephemeris files for given datetime")
    else:
        print("done")
