# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

MAKEFILE = Makefile

# Progam to make
EXE	= ancnrt

# Object modules for EXE
OBJ	= ancnrt.o \
	  readgrib.o \
	  fillenv.o \
	  ANCroutines.o \
	  rdattr.o \
	  countann.o \
	  rdfiles.o \
	  rdsdsid.o \
	  pexit.o \
	  julian.o \
	  lnstrg.o \
	  leap.o \
	  resize_oz.o

# Include file locations
INCLUDE	= -I$(INCDIR)/swfinc -I$(LIB3_INC)

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -lmfhdf -ldf -ljpeg -lz

include $(MAKEFILE_APP_TEMPLATE)
