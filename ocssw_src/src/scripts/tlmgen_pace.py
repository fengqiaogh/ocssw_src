#! /usr/bin/env python3

import argparse
import glob
import logging
import os.path
import sys
from io import BytesIO

import pandas as pd
from telemetry import ccsdspy
from telemetry.PacketUtils import *

__version__ = "1.1.0 (2024-05-20)"


def main():
    print("tlmgen_pace", __version__)

    # Read command line options
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Convert S-band housekeeping telemetry to CSV for specified fields",
        epilog="""
EXIT Status:
0   : Success
1   : Fatal error
101 : Non-fatal file open error
102 : Invalid file or instrument from CFE header
103 : Invalid packet header
104 : Invalid packet [header/datatype/length]
110 : No valid packets found
120 : Multiple warnings; see log
""",
    )
    parser.add_argument(
        "ifile",
        type=str,
        help="path to S-band telemetry (HSK) file OR list of input files, one per line",
    )
    parser.add_argument(
        "-o", "--ofile", type=str, help="output CSV file; defaults to ifile.csv"
    )
    parser.add_argument(
        "--packetdir",
        type=str,
        help="path to directory containing packet structures for desired mnemonics.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="print status messages",
    )
    args = parser.parse_args()
    status = 0

    loglevel = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    # locate CSV directory
    if args.packetdir is None:
        pktDir = os.path.join(
            os.getenv("OCDATAROOT"), "telemetry", "pace", "monitoring"
        )
    else:
        pktDir = args.packetdir
    if not os.path.exists(pktDir):
        logging.error(f"ERROR: The directory {pktDir} does not exist.")
        return 1

    # read conversions
    csvfile = os.path.join(pktDir, "LinearConverters.csv")
    if os.path.exists(csvfile):
        conversions = pd.read_csv(csvfile)
    for column in ["slope", "intercept"]:
        if isinstance(conversions[column][0], str):
            conversions[column] = [eval(v) for v in conversions[column]]

    # read packet definitions
    packetDef = {}
    for csvfile in glob.glob(os.path.join(pktDir, "APID[0-9]*.*.csv")):
        apid = int(
            os.path.basename(csvfile).split(".", maxsplit=1)[0].replace("APID", "")
        )
        packetDef[apid] = ccsdspy.FixedLength.from_file(csvfile)

        # attach conversion, if present
        for name in [f._name for f in packetDef[apid]._fields]:
            row = conversions.loc[conversions["mnemonic"] == name]
            if len(row) > 0:
                packetDef[apid].add_converted_field(
                    name,
                    name,
                    ccsdspy.converters.LinearConverter(
                        slope=row.slope.values[0], intercept=row.intercept.values[0]
                    ),
                )

    ofile = f"{args.ifile}.csv" if args.ofile is None else args.ofile

    # Is input file tlm or list of files?
    filelist = []
    infile = os.path.expandvars(args.ifile)

    try:
        with open(infile, mode="rt") as flist:  # try to read as list of files
            try:
                input_list = True
                for ifile in flist:
                    filelist.append(os.path.expandvars(ifile.rstrip()))
            except UnicodeDecodeError:
                input_list = False  # contains invalid byte - infile is binary

    except IOError as e:
        logging.error(f"{e}; exiting")
        return 1

    if not len(filelist):  # if that didn't work, it's a single HSK file
        filelist.append(infile)

    writeheader = True

    # Step through all input files
    for filename in filelist:
        logging.info(f"Reading {filename}")
        fname = os.path.basename(filename)
        dictList = []

        try:
            ifile = open(filename, mode="rb")
        except IOError as e:
            status = 120 if status > 0 else 101
            logging.warning(f"{e}; continuing")
            continue  # input_list errors already caught above

        # Read any file header(s)
        filehdr = readFileHeader(ifile)
        if filehdr:
            logging.info(filehdr)
            logging.info("")

            # Is it the right kind of file?
            desired = (
                filehdr["subtype"] == 101
                and filehdr["length"] == 64
                and filehdr["SCID"] == b"PACE"
                and filehdr["processorid"] in (1, 2, 30)
            )

        if not filehdr or not desired:
            status = 120 if status > 0 else 102
            if input_list:
                logging.warning(f"File {filename} has invalid header; continuing")
                continue  # go to next file
            else:
                logging.error(f"File {filename} has invalid header; returning")
                return status

        # read CCSDS packets
        for packet in ccsdspy.utils.iter_packet_bytes(
            ifile, include_primary_header=True
        ):
            data = packet[6:]
            header = ccsdspy.utils.read_primary_headers(BytesIO(packet))
            for k, v in header.items():
                header[k] = v[0]  # remove outer array
            logging.info(header)

            # check for invalid header
            if header["CCSDS_PACKET_LENGTH"] > 16378:
                status = 120 if status > 0 else 103
                logging.warning(
                    f"File {filename} contains invalid CCSDS packet header: {header}"
                )
                continue  # go to next packet

            # check for truncated data
            if len(data) < header["CCSDS_PACKET_LENGTH"] + 1:
                status = 120 if status > 0 else 104
                logging.warning(
                    f"File {filename} has unexpected EOF: expected"
                    f" {header['CCSDS_PACKET_LENGTH']+1} more bytes, got {len(data)}"
                )
                break  # done with this file

            # check for invalid timestamp
            if (
                header["CCSDS_SECONDARY_FLAG"] == 1
                and data[0] > 112  # after  2017-07-18T05:49:15Z
                and data[0] < 192  # before 2060-01-28T16:50:35Z
            ):
                timestamp = tai58_as_datetime(readTimestamp(data))
            else:
                continue  # done with this packet

            # parse defined mnemonics
            apid = header["CCSDS_APID"]
            if apid in packetDef.keys():
                myDict = (packetDef[apid]).load(BytesIO(packet))
                for key, val in myDict.items():
                    if not key.endswith("time"):
                        outdict = {}
                        outdict["filename"] = fname
                        outdict["time_val"] = timestamp.strftime(
                            "%Y-%m-%d %H:%M:%S.%f"
                        )[:-3]
                        outdict["var"] = key
                        outdict["value"] = val[0]
                        outdict["alert_type"] = ""  # populated later
                        outdict["recorded"] = ""  # database ingest datetime
                        dictList.append(outdict)

        # end (for packet in packets)

        # close input file
        ifile.close()

        # write new dataframe
        if len(dictList) > 0:
            logging.info(f"Writing {len(dictList)} records from {fname} to {ofile}")
            df = pd.DataFrame(dictList)
            df.to_csv(ofile, sep=",", index=False, header=writeheader, mode="a")
            writeheader = False

    # end (for filename in filelist)

    if writeheader:  # true if no records written
        logging.warning(f"No requested packets found.")
        status = 120 if status > 0 else 110

    if status:
        logging.warning(f"Exiting with status code {status}")
    return status


if __name__ == "__main__":
    sys.exit(main())
