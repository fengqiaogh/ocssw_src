#****************************************************************************
#MODIS_36100.t
#
#Description: Level 1B MODIS SMF SEED file
#
#Revision History:
#$Log: MODIS_36100.t,v $
#Revision 1.2  2004-12-23 12:08:22-05  ltan
#Terra MOD_PR02 code and LUT update to V5.0.0 (Collection 5)
#
#Revision 1.1  2002/02/06 18:53:20  alice
#Initial revision
#
# Revision 1.1  2002/02/06  18:43:38  alice
# Initial revision
#
# Revision 1.1  2000/08/18  17:37:23  rogers
# Initial revision
#
#
#  Revision 01.01  1997/08/27  
#  Revision to replace LABEL field MCST with MODIS.
#  Zhidong Hao (hao@barebackride.gsfc.nasa.gov)
#
#  Revision 01.00  1996/02/05  09:13:00 
#  Revision to replace LABEL field MODIS with MCST.
#  Joan Baden(baden@highwire.gsfc.nasa.gov)
#
#Team-unique Header:
#  This software is developed by the MODIS Characterization Support
#  Team for the National Aeronautics and Space Administration,
#  Goddard Space Flight Center, under contract NAS5-32373.
#
#References and Credits
#  Written by Joan R. Baden 
#	      
#             Research and Data Systems Corporation (RDC)
#	      SAIC/GSC MODIS Support Office
#             7501 Forbes Blvd, Suite 103
#             Seabrook MD 20706  
#
#	      301-352-2128
#
#	      baden@highwire.gsfc.nasa.gov
#	   
#
#Design Notes
#
#  Externals: none
#
#END
#*****************************************************************************
%INSTR = MODIS
%LABEL = MODIS
%SEED  = 36100

#Dynamic timestamp message
MODIS_F_OUT_OF_MEMORY	    Out of memory running...%sTIME:%s
MODIS_F_MEM_FREE_FAIL        Memory Free failed running...%sTIME:%s
MODIS_F_FILE_NOT_FOUND	    File not found running...%sTIME:%s  
MODIS_F_READ_ERROR	    Read Error running...%sTIME:%s
MODIS_F_WRITE_ERROR          Write Error running...%sTIME:%s
MODIS_F_OUT_OF_RANGE	    Out of range running...%sTIME:%s
MODIS_W_OUT_OF_RANGE	    Warning: Out of range running...%sTIME:%s
MODIS_W_TIME_INCORRECT       Warning:Time incorrect running...%sTIME:%s
MODIS_F_NOK		    Fatal error(NOK) running...%sTIME:%s
MODIS_F_HDF_ERROR            Fatal error (HDF) running...%sTIME:%s
MODIS_F_NO_MORE              Fatal Error: No more data available (EOF) running...%sTIME:%s
MODIS_S_NO_MORE              No more data available (EOF) running...%sTIME:%s
MODIS_F_FILE_NOT_OPENED	    File not opened running...%sTIME:%s
MODIS_F_FILE_NOT_CREATED     File not created running...%sTIME:%s
MODIS_F_INVALID_ARGUMENT     Invalid argument running...%sTIME:%s
MODIS_E_TESTING		    Testing error running...%sTIME:%s
MODIS_S_OK      	            Completed Processing Successfully...%sTIME:%s

#End of MODIS_36100.t
