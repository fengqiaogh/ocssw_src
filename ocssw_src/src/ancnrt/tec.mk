# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

MAKEFILE = Makefile

# Progam to make
EXE	= tec

# Object modules for EXE
OBJ	= \
tec.o

# Include file locations
INCLUDE	= -I$(INCDIR)/swfinc -I$(LIB3_INC)

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -lhdf5 -lz

include $(MAKEFILE_APP_TEMPLATE)
