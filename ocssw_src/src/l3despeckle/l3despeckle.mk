# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

MAKEFILE = $(EXE).mk

# Progam to make
EXE	= l3despeckle

# Object modules for EXE
OBJ	= l3despeckle.o \
          l3despeckle_input.o

# Include file locations
INCLUDE = -I$(INCDIR)/swfinc -I$(LIB3_INC)

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -lbin++ -lgenutils -lbin -lhdfutils -lmfhdf \
	  -ldf -lhdf5 -ljpeg -lz

LD = $(CXX)

include $(MAKEFILE_APP_TEMPLATE)






