# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

MAKEFILE = $(EXE).mk

# Progam to make
EXE	= ice2hdf

# Object modules for EXE
OBJ	= ice2hdf.o


# Include file locations
INCLUDE	= -I$(LIB3_INC) -I$(INCDIR)/swfinc -I$(INCDIR)/utils

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -lhdfutils -lmfhdf -ldf -ljpeg -lz 


include $(MAKEFILE_APP_TEMPLATE)
