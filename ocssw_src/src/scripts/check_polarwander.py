#! /usr/bin/env python3
"""
Verify validity of USNO polar wander file
Format description: http://maia.usno.navy.mil/ser7/readme.finals2000A
"""

import re
import sys


def isValidLine(line):
    valid = True
    ipflag = re.compile("^[IP]$")  # IERS (I) or Prediction (P) flag

    try:

        int(line[0:2])     # year
        int(line[2:4])     # month number
        int(line[4:6])     # day of month
        float(line[6:16])  # fractional Modified Julian Date (MJD UTC)

        valid &= ipflag.match(line[16]) is not None
        # IP flag for Bull. A polar motion values
        float(line[17:27])  # Bull. A PM-x (sec. of arc)
        float(line[27:36])  # error in PM-x (sec. of arc)
        float(line[36:46])  # Bull. A PM-y (sec. of arc)
        float(line[46:57])  # error in PM-y (sec. of arc)

        valid &= ipflag.match(line[57]) is not None
        # IP flag for Bull. A UT1-UTC values
        float(line[58:68])  # Bull. A UT1-UTC (sec. of time)
        float(line[68:78])  # error in UT1-UTC (sec. of time)

        if len(line[78:95].split()) > 0:  # optional section
            float(line[78:86])   # Bull. A LOD (msec. of time)
            float(line[86:95])   # error in LOD (msec. of time)

        if len(line[95:134].split()) > 0:  # optional section
            valid &= ipflag.match(line[95]) is not None
            # IP flag for Bull. A nutation values
            float(line[96:106])   # Bull. A dX wrt IAU2000A Nutation (msec. of arc)
            float(line[106:115])  # error in dX (msec. of arc)
            float(line[115:125])  # Bull. A dY wrt IAU2000A Nutation (msec. of arc)
            float(line[125:134])  # error in dY (msec. of arc)

        if len(line[134:].split()) > 0:  # optional section
            float(line[134:144])  # Bull. B PM-x (sec. of arc)
            float(line[144:154])  # Bull. B PM-y (sec. of arc)
            float(line[154:165])  # Bull. B UT1-UTC (sec. of time)
            float(line[165:175])  # Bull. B dX wrt IAU2000A Nutation (msec. of arc)
            float(line[175:185])  # Bull. B dY wrt IAU2000A Nutation (msec. of arc)

    except (ValueError, IndexError):
        return False
    except Exception as e:
        print('Unexpected Exception: {:}'.format(e))
        return False

    return valid


if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('Checks a USNO Polar Wander file for validity.')
        print('Exits with status = 1 upon first invalid line.')
        callseq = sys.argv[0] + ' filename'
        print('Usage:', callseq)
        sys.exit(1)

    infile = open(sys.argv[1], 'r')
    for line in infile:
        if not isValidLine(line):
            print('Format error in line:')
            print(line)
            infile.close()
            sys.exit(1)

    infile.close()
    sys.exit(0)
