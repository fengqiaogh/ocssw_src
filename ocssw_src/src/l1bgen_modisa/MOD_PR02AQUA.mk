###########################################################################
# !make
#
# !Makefile Name: MOD_PR02AQUA.mk
#
# !Description:
#	This makefile used to generate MOD_PR02AQUA.exe
#
# !Variables:
#
# 	VARIABLES     DESCRIPTION
# 	~~~~~~~~~     ~~~~~~~~~~~
# 	TARGET        Program executable name
# 	ADD_CFLAGS    Additional compiler options
# 	LIB           Libraries
# 	OBJ           Object files
# 	INC           Include files
# 	INC_FILES     Additional include files
# 	INC_HDF5OPT   Additional HDF5 include files
# 	LIB_HDF5OPT   Additional HDF5 libraries
# 	LIB_DL        Libraries supporting dynamic loading
#
# !Env Variables
#
# 	ENV VARIABLES DESCRIPTION
# 	~~~~~~~~~~~~~ ~~~~~~~~~~~
# 	CFLAGS        Compiler flags which are set by the script 
#		      (e.g. v2_n32_f77)
# 	CC            The compiler set by the script
# 	PGSINC        Include directory of PGS toolkit set by the script
# 	HDFINC        Include directory of HDF set by the script
# 	HDFEOS_INC    Include directory of HDFEOS set by the script
# 	PGSLIB        Library of PGS toolkit set by the script
# 	HDFLIB        Library of HDF set by the script
# 	HDFEOS_LIB    Library of HDFEOS set by the script
# 	BRAND         Operating system set by the script
# 	HDF5INC       Include directory of HDF5 set by the script 
# 	HDF5LIB       Library of HDF5 set by the script
# 	SZIPHOME      Library of SZIP compression set by the script 
#
#
# !Team-Unique Header:
#    	This software was developed by MODIS Science Data Support Team
#    	for the National Aeronautics and Space Administration,
#    	Goddard Space Flight Center, under contract NAS5-32373.
#
# !Revision History:
#
#  DATE                AUTHOR                  DESCRIPTION
#  ----------          -----------             -----------
#  01/21/2011          Liqin Tan               Updated for compatibility with SDP Toolkit 5.2.16.
#  11/02/2004          Liqin Tan               Added optional support for HDF5.
#  10/26/2004          Liqin Tan               Changed to use -ffloat-store
#                                              option when $(BRAND) is linux.
#  2/7/2002            Alice Isaacman          Changed executable name to
#                                              MOD_PR02AQUA.
#  3/23/2000           Jim Rogers              Added -DALLOW_MISSING_GRANULES
#                                              to implement production rules
#                                              change (CCR 519).
#  5/26/98             Gang Ye                 Updated for DAAC's requirement
#  4/23/98             David Catozzi           Updated for Version 2.1
#  4/15/97             Zhidong Hao             Updated for version 2 delivery 
#  Nov. 1996           Neal Devine             Initial revision
#
# !END
##############################################################################
# Define executable name
TARGET = MOD_PR02AQUA.exe

# Combine the includes using pre-defined includes and your includes
INC_HDF5OPT = -I$(HDF5INC)
INC =	-I$(PGSINC) \
	-I$(HDFINC) \
	-I$(HDFEOS_INC) \
	$(INC_HDF5OPT:-I=)

# Combine the libraries using pre-defined libraries and your library
LIB_HDF5OPT = /$(HDF5LIB)/libhdf5.a
LIB_DL = $(BRAND:linux=-ldl)$(BRAND:sgi64=-lc)
LIB =	-L$(PGSLIB) -lPGSTK \
	-L$(HDFLIB) -lmfhdf -ldf -ljpeg -lz  \
	-L$(HDFEOS_LIB) -lhdfeos -lGctp -lm \
	$(LIB_HDF5OPT://libhdf5.a=) \
	$${SZIPHOME:+-lsz} $${SZIPHOME:+$(LIB_DL:$(BRAND)=)} -ljpeg -lz -lm

# Define additional C flags
ADD_CFLAGS = -O3 -DALLOW_MISSING_GRANULES -D_POSIX_C_SOURCE=1 $(BRAND:linux=-ffloat-store)

# Additional include files
INC_FILES = FNames.h Granule.h GranuleP.h HDF_Lib.h L1B_Setup.h\
            L1B_SetupP.h L1B_Tables.h Metadata.h MetadataP.h\
            PGS_Error_Codes.h Preprocess.h PreprocessP.h\
            Reflective_Cal.h Reflective_CalP.h PGS_MODIS_36110.h\
            Emissive_Cal.h

# Define object files
OBJ =	L1B.o Preprocess.o L1B_Tables.o Granule.o HDF_Lib.o Metadata.o\
	Reflective_Cal.o Emissive_Cal.o L1B_Setup.o 

# Make the process
$(TARGET) :	$(OBJ)
	$(CC)	$(CFLAGS) $(ADD_CFLAGS:sgi64=) $(OBJ) $(LIB) -o $(TARGET)

.c.o :	$(INC_FILES)
	$(CC)	$(CFLAGS) $(ADD_CFLAGS:sgi64=) $(INC) -c  $< -o $@

# Delete object files
clean:
		/bin/rm -f $(OBJ)

# Delete the executable file $TARGET
clean_exe:
		/bin/rm -f $(TARGET)

#******************** End of make file ********************

