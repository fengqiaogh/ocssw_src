# Changelog for l1agen_oci

Follow the format below when adding modifications to the code. Put most recent change on top. <br>
You can find the older modification history in a code block below. This is before this markdown file was made and all modification history were comments in the code,.
### [version] - YYYY-MM-DD
- Things you changed, or comments.
- ...
- ...

-----------------------

### [01.19.00] - 2024-07-12
- Code clean up version.
- Removed `l0info_oci` because there's a python version.
- Deleted comment in CMakeLists that builds the code.
- Formatted all code to follow C/C++ style




<br><br><br><br><br>

## Code log up until July 2024 for `l1agen_oci.cpp`
```
//  Modification history:
//  Programmer     Organization   Date     Ver     Description of change
//  ----------     ------------   ----     ---     ---------------------
//  Joel Gales     FutureTech     09/20/18 0.10    Original development
//                                                 based on IDL routines
//                                                 developed by F. Patt
//  Joel Gales     SAIC           10/26/18 0.20    Complete initial alpha version
//  Joel Gales     SAIC           12/20/18 0.30    Complete initial version
//                                                 of SWIR bands
//  Joel Gales     SAIC           04/26/19 0.40    Implement code changes from
//                                                 F. Patt
//  Joel Gales     SAIC           05/20/19 0.50    Implement code changes from
//                                                 F. Patt
//  Joel Gales     SAIC           06/04/19 0.60    Implement code changes from
//                                                 F. Patt
//  Joel Gales     SAIC           07/01/19 0.70    Add support for outlist
//  Joel Gales     SAIC           07/23/19 0.75    Implement code changes from
//                                                 F. Patt
//  Joel Gales     SAIC           08/02/19 0.80    Flag and ignore packets
//                                                 greater than 1200 bytes.
//                                                 Remove common routines and
//                                                 Place in common.cpp
//  Joel Gales     SAIC           08/09/19 0.81    Fix memory overwrite bug in
//                                                 unpack_ccd_packet() in the
//                                                 ossdata array.  Initialize
//                                                 "lines" and "bands" arrays.
//                                                 Add support for granules with
//                                                 missing blue/red bands.
//  Joel Gales     SAIC           10/25/19 0.82    Exit if EOF before finding
//                                                 good packet
//  Joel Gales     SAIC           11/20/19 0.85    Implement code changes from
//                                                 F. Patt
//  Joel Gales     SAIC           11/25/19 0.86    Change 60 to maxsc fpr dspn
//                                                 comparison
//  Joel Gales     SAIC           11/27/19 0.87    Change L1A name to standard
//                                                 Trap granules with no ancillary
//                                                 granules
//  Joel Gales     SAIC           12/05/19 0.88    Add nametage command line
//                                                 option
//  Joel Gales     SAIC           01/03/20 0.90    Update and add additional
//                                                 telemetry fields
//  Joel Gales     SAIC           01/22/20 0.91    Add scomp/maxgap
//  Joel Gales     SAIC           01/23/20 0.92    Check for EOF when checking
//                                                 for zero science pixels
//                                                 Zero out sci/dark fields
//  Joel Gales     SAIC           01/28/20 0.93    Bug fixes and updates
//  Joel Gales     SAIC           02/07/20 0.94    Implement changes from F.Patt
//                                                 (020720)
//  Joel Gales     SAIC           02/20/20 0.95    Fixed bugs with bagerr, ragerr,
//                                                 digerr and dautemp
//  Joel Gales     SAIC           02/25/20 0.96    Fixed bug in ccsds..
//                                                 Writing 6 bytes into uint32_t
//  Joel Gales     SAIC           03/04/20 0.97    Implement changes from F.Patt
//                                                 (022820)
//  Joel Gales     SAIC           03/24/20 0.98    Fix HAM_side issue
//  Joel Gales     SAIC           04/02/20 0.99    Implement SWIR mode changes
//  Liang Hong     SAIC           04/28/20 0.9901  Fixed last scan_line_attributes
//                                                 records in last output granule
//  Liang Hong     SAIC           05/20/20 0.9902  Handle dark zone science data
//  Liang Hong     SAIC           08/11/20 0.9903  Handle 0 pix 0 band sci/cal data in l1a
//                                                 "no data" and order variation in input
//  Liang Hong     SAIC           08/24/20 0.9904  fixed a bug to close file with 0-scan
//                                                 ; demand non-negative indices
//  Liang Hong     SAIC           09/02/20 0.9905  HKT and science data overlap check; order
//                                                 SWIR bands in ascending wavelength
//  Liang Hong     SAIC           10/28/20 0.9906  handle fill values in SWIR
//  Liang Hong     SAIC           10/29/20 0.9907  APID for the MCE HK packet changed to 713
//  Liang Hong     SAIC           11/23/20 0.9908  fixed rare start/end time error; SWIR band
//                                                 -specific pixel shifts as input option;
//                                                 run with optional granule start time;
//                                                 fixed science packet sequence error flag
//                                                 fixed HAM_side value after index nmce
//  Liang Hong     SAIC           12/01/20 0.99.00 fixed number_of_filled_scans in metadata
//  Liang Hong     SAIC           04/22/21 0.99.10 fixed no ancillary data exit conditions
//  Liang Hong     SAIC           06/17/21 0.99.20 generate 1 telemetry entry when no data
//                                                 bug fix in duplicated reading of last tlm
//  Liang Hong     SAIC           07/27/21 0.99.21 return 110 for No ancillary packets
//  Liang Hong     SAIC           01/07/22 1.00.00 temperature fields update in HKT packets
//                                                 SWAP_4 in common.h updated
//  Liang Hong     SAIC           01/11/22 1.00.01 OCI SWIR Raw mode reading correction
//  Liang Hong     SAIC           01/25/22 1.00.11 blue and red spectral mode order correction
//  Liang Hong     SAIC           03/11/22 1.01.00 added telemetry for the solar calibrator;
//                                                 references in metadata; write ancillary_tlm
//  Liang Hong     SAIC           03/30/22 1.02.00 SPCA and lunar stare data types added
//  Liang Hong     SAIC           04/14/22 1.03.00 update error checking and handling
//  Liang Hong     SAIC           04/21/22 1.03.01 updated packet # threshold; mode table check
//  Liang Hong     SAIC           04/24/22 1.03.02 exit read packets when endfile
//  Liang Hong     SAIC           05/11/22 1.04.00 update of APID list and check all for spin #
//                                                 added navigation data from hkt input
//  Liang Hong     SAIC           05/27/22 1.04.01 limit packet reading to maxpkts
//  Liang Hong     SAIC           06/24/22 1.05.00 added CCD masks; fixed bugs in nav data
//  Liang Hong     SAIC           11/18/22 1.06.00 noSPW option set for OCI test data only
//  Liang Hong     SAIC           11/22/22 1.07.00 granule starts with specified start time
//  Liang Hong     SAIC           11/28/22 1.08.01 updated start time, outfile and --noSPW options
//  Liang Hong     SAIC           12/05/22 1.08.03 clear outlist if no L1A generated, return 111
//  Liang Hong     SAIC           01/10/23 1.09.04 added data type to last column of outlist
//  Gwyn Fireman   SAIC           02/02/23 1.09.05 CF-compliant output; new CDL file; history attribute
//  Liang Hong     SAIC           02/28/23 1.10.05 update OCI solar cal datatypes
//  Liang Hong 	   SAIC           05/01/23 1.11.00 navigation data read and fill update
//  Liang Hong     SAIC           05/15/23 1.11.01 bug fix in PACE navigation data read
//  Liang Hong     SAIC           05/30/23 1.12.01 added usage of tilt flag
//  Liang Hong     SAIC           06/23/23 1.13.00 Fill value update; fixed cross date data search issue
//  Liang Hong     SAIC           06/26/23 1.14.00 Added maxgap as an option allowed missing scans/file
//  Gwyn Fireman   SAIC           07/10/23 1.14.01  Read global metadata from json file
//  Liang Hong     SAIC           09/12/23 1.15.00 Added SCA diffuser; 1 granule given start, granule_len
//  Wei Jiang      SAIC           03/15/24 1.16.00 maxtime (mtime) does not update if there's a data type change
//  Wei Jiang      SAIC           03/18/24 1.16.00 handles single .oci file that have 64 byte header intact
```
<br><br>
## Code log up until July 2024 for `common.cpp`
```
//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     FutureTech     08/02/19 0.10  Original development
//  Liang Hong     SAIC           04/14/22 0.20  l1agen_oci v1.03.00
//  Liang Hong     SAIC           04/24/22 0.21  l1agen_oci v1.03.02
//  Liang Hong	   SAIC           05/05/22 0.30  l1agen_oci v1.04.00
//  Liang Hong     SAIC           11/18/22 0.40  l1agen_oci v1.06.00
//  Jakob Lindo    SSAI           03/16/24 0.50  Added aggregation indication
```