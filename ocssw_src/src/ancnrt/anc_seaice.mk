# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

MAKEFILE = Makefile

# Progam to make
EXE	= anc_seaice

# Object modules for EXE
OBJ	= anc_seaice.o \
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
	  resize_oz.o \
	  mk_smooth_ice_map.o

# Include file locations
INCLUDE	= -I$(INCDIR)/swfinc -I$(LIB3_INC)

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -laquarius -lmfhdf -ldf -ljpeg -lz

include $(MAKEFILE_APP_TEMPLATE)
