# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

MAKEFILE = $(EXE).mk

# Progam to make
EXE	= l1agen_seawifs

# Object modules for EXE
OBJ	= main_l1agen.o \
	  swl0_utils.o \
	  swl1_hdf.o \
	  getl0indx.o \
	  printindx.o \
          temporal_anomaly.o \
	  getl0scene.o \
	  getl0scene_nav.o \
	  printscene.o \
          getorbdata.o \
          getnavdata.o \
          printnav.o \
          getl1rec.o \
          getorbnum.o \
          instlm.o \
          mkmeta.o \
          getmetrics.o

# Include file locations
INCLUDE	= -I$(INCDIR)/utils -I$(INCDIR)/swfinc -I$(LIB3_INC)

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -lseawifs -lswfnav -lnav \
          -lgenutils -lmfhdf -ldf -ljpeg -lz

include $(MAKEFILE_APP_TEMPLATE)
