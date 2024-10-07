# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

MAKEFILE = Makefile

# Progam to make
EXE	= solar_xray

# Object modules for EXE
OBJ	= \
solar_xray.o

# Include file locations
INCLUDE	= -I$(INCDIR)/swfinc -I$(LIB3_INC)

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -lgenutils -lhdf5 -lz

include $(MAKEFILE_APP_TEMPLATE)
