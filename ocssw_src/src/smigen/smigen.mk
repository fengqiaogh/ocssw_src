# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

MAKEFILE = Makefile

# Progam to make
EXE	= smigen

# Object modules for EXE
OBJ	= smigen.o \
          smigen_input.o \
	  put_smi.o

# Include file locations
INCLUDE	= -I$(INCDIR)/utils -I$(INCDIR)/swfinc -I$(LIB3_INC)

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

LLIBS 	= -lbin++ -lbin -lhdfutils -lgenutils -lmfhdf -ldf -lhdf5 -ljpeg -lz


LOCAL_CXXFLAGS = -DHDF5

LD = $(CXX)

include $(MAKEFILE_APP_TEMPLATE)

ifdef ppcDarwinArchitecture
  LOCAL_CFLAGS = -O0 -DHDF5
endif
