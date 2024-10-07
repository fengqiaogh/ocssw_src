1. README for MOD_PR02AQUA (AQUA) 
   Version: 6.2.1
   Deliver to SDST: August 15, 2017 


2. POINT OF CONTACT:  

   Developers:
   ----------
     Xu Geng
     SSAI
     MODIS Characterization Support Team
     EMAIL: xu.geng@ssaihq.com
     TEL: 301-867-2057
     10210 Greenbelt Road, Suite 600
     Lanham, Maryland 20706

   SDST Contacts:
   -------------
     James Kuyper
     Global Science and Technology
     MODIS Science Data Support Team
     EMAIL: James.R.Kuyper@nasa.gov
     TEL: 301-867-2064
     10210 Greenbelt Road, Suite 600
     Lanham, Maryland 20706

3.BUILD:
 a.) Create the file PGS_MODIS_36110.h in the $PGSINC directories, and the
     file PGS_36000 in the $PGSMSG directories (every toolkit directory tree
     you want to compile with). This requires setting the PGS environment
     variables first as noted below. (Note: if the files already exist, this
     step should be skipped unless the files require updating, ie the
     MODIS_36100.t file has changed. To accomplish this command, write
     permission is required for $PGSINC and $PGSMSG directories.)

     on moddev-c7, the following commands would be executed:
            source /data/common/util/bin/setup-5.2.19.csh
            $PGSBIN/smfcompile -f MODIS_36100.t -r -i

 b.) Edit the Process Control Files to set the file paths to the current
     working directory. Check the L1B output file names and path, make sure
     you have write privileges. The L1A input and L1B output file names
     should be changed whenever new input/output is desired.
     NOTE: the value of PGEVersion is now set in the PCF file for PGE02.

     See other comments about the PCF file in the "Satellite Notes" and
     "Reprocessing Notes" Sections later in this README file.

 c.) Enviroment variables for PGS, HDF, HDFEOS, PGSMSG, and SGI_MODE must be 
     set before the makefile is used to compile the code.  First set up 
     command aliases:

     on moddev-c7, the following commands would be executed:
            source /data/common/util/bin/setup-5.2.19.csh

 d.) PCF file notes:
     Two run-time parameters in the PCF are used to set the current code/LUT 
     version and to turn production of 250m and 500m resolution data off when 
     all scans of a granule are in night mode:

     PGE02 Version Number
     ===================
     In the PCF file, LUN 800500 identifies the PGE02 Version number, which is
     checked against code macro "PGE02_VERSION" defined in MOD_PR02 code. The
     value defined after the last "|" defines the PGE02 version.

     Example:
     800500|PGE02 Version|6.2.1

     MCST Version Number
     ===================
     In the PCF file, LUN 800610 identifies the MCST Version number, which is      
     checked against "MCST_Version" as identified within the LUT hdf files 
     defined elsewhere in the PCF. The value defined after the last "|" defines 
     the MCST version.

     Example:
     800610|MCSTLUTVersion|6.2.1.0_Aqua

     Create High Resolution Night Mode Data Output Off/On
     ====================================================
     In the PCF file, LUN 800615 controls whether the 250m and 500m data sets
     will be written out when all scans in a granule are NIGHT mode.  To enable 
     writing out night mode High Resolution data, set the flag to 1; to disable, 
     set to 0.

     For disabling creating 250m and 500m output data sets when in night mode:
     800615|Write_Night_Mode_HiRes_Data|0

     For enabling creating 250m and 500m output data sets when in night mode:
     800615|Write_Night_Mode_HiRes_Data|1

     Satellite Notes
     ===============
     To install the PGE for Aqua, use the following files and set their names 
     in the PCF file appropriately:

     MYD021KM.mcf
     MYD02HKM.mcf
     MYD02QKM.mcf
     MYD02OBC.mcf
     MYD02_Emissive_LUTs.hdf
     MYD02_QA_LUTs.hdf
     MYD02_Reflective_LUTs.hdf

     In the PCF file, LUN 800510 identifies the platform (Terra or Aqua). The
     value defined after the last "|" defines the satellite.

     For Aqua:
     800510|Satellite; AM1M=Terra, PM1M=Aqua|PM1M

     The value defined for this LUN should be consistent with the MCF and LUT 
     files installed with the PGE and with their names as defined elsewhere in 
     the PCF.

     Reprocessing Notes
     ==================
     Two run-time parameters in the PCF are used to set the values of the ECS 
     core metadata fields "ReprocessingPlanned" and "ReprocessingActual":

     800600|ReprocessingPlanned|further update is anticipated
     800605|ReprocessingActual|processed once

     For forward processing, set the value of "ReprocessingActual" to 
     "processed once".  Set the value of "ReprocessingPlanned" to "further 
     update is anticipated".

     Processing Environment
     ======================
     LUN 800550 has been identified by SDST as the LUN to use if the metadata
     value "ProcessingEnvironment" is to be read from the PCF file.  To avoid
     many "error" messages when this value is left blank (even though that is
     permissible, the SDP Toolkit issues error messages when it cannot find
     the variable), MCST has chosen to have L1B determine the
     "ProcessingEnvironment" variable by using the POSIX function "getenv" to
     mimic the action of the command-line query "uname -a" and fill the
     variable in at all times.  Hence the LUN is commented out in the PCF file.

     Processing Center
     ======================
     LUN 800620 has been identified by SDST as the LUN to use for the metadata
     value "ProcessingCenter", which is to be read from the PCF file.  Set 
     this value to the location the data is processed, e.g.

     800620|ProcessingCenter|MODAPS

 e.) Use the makefile to create the executable "MOD_PR02AQUA.exe".
     Optionally, remove old object files and executables:
          make clean -f MOD_PR02AQUA.mk

     Compile and link:
          make -f MOD_PR02AQUA.mk

 f.) Running the code:
     The name of the PCF file must be defined by the environment variable:
          PGS_PC_INFO_FILE  

     Enter the command:
          MOD_PR02AQUA.exe


4. ENVIRONMENT VARIABLE:
   PGSINC       defines directory for SDP toolkit headers
   HDFINC       defines directory for HDF headers
   HDFEOS_INC   defines directory for HDF-EOS headers
   CFLAGS       defines compiler flags
   CC           defines the C compiler
   PGSLIB       defines directory for SDP toolkit library
   HDFLIB       defines directory for HDF library
   HDFEOS_LIB   defines directory for HDFEOS library
   HDF5INC      defines directory for HDF5 headers
   HDF5LIB      defines directory for HDF5 library
   SZIPHOME     defines directory for SZIP compression library
   BRAND        defines the operating system


5. MAPI: 
   No M-API used.


6. HDF: 
   HDF version 4.2r10 and HDF-EOS version 2.19


7. SDP TOOLKIT: 
   SDP-Toolkit versions 5.2.19


8. This code was developed on the MCST computer server "luna",
   OS Linux version 3.10.0-514.6.1.el7.x86_64,
   16 Intel(R) Xeon(R) CPU E7330  @ 2.40GHz
   gcc version 4.8.5 20150623 (Red Hat 4.8.5-11)


9. ANCILLARY DATA:
   This code uses no ancillary data.


10.MODIS INPUT PRODUCTS:
   Three L1A (MOD_PR01) files, one Geolocation (MOD_PR03) file. See 
   "MOD_PR02AQUA_pr.txt" for more details.


11.OTHER INPUTS:
   Three lookup table files (MYD02_Emissive_LUTs.hdf, MYD02_Reflective_
   LUTs.hdf, and MYD02_QA_LUTs.hdf). See "MOD_PR02AQUA_pr.txt" for more 
   details.


12.MODIS OUTPUT PRODUCTS:
   Four L1B (MOD_PR02AQUA) files (MYD021KM, MYD02HKM, MYD02QKM, and MYD02OBC).


13.PROBLEMS:
   None.


14.DELIVERY SIZE:
   Total size of untarred delivery is approximately 3526 MB.
   See PACKINGLIST for the information of files delivered.


15.OUTPUT FILE SIZE:
   Total size of expected output files (203 scans) is approximately 1838 MB.
   The stack size required is approximately 68 MB.


16.TESTS PERFORMED:

   The code has been tested with MODIS data collected on the Aqua satellite 
   which was produced by PGE01 at SDST.  It runs to completion, producing 
   the 4 L1B output files and the outputs look reasonable.


17.EXIT CODES:

   0:	normal end; 
   1: 	any fatal errors; and 
   233: number of scans less or equal to 0


18.ERROR LIST:
   See the file "error_messages.txt" for a description of the MOD_PR02AQUA
   error messages.

