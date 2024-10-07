# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

MAKEFILE = Makefile

# Progam to make
EXE	= o3nrt_omi

# Object modules for EXE
OBJ	= \
o3nrt_omi.o \
ANCroutines.o \
rdgrid.o \
countann.o \
fillenv.o \
pexit.o \
gregor.o \
leap.o \
julian.o \
getfn.o \
lnstrg.o \
epochbk.o \
julday.o \
rdattr.o \
rdfiles.o \
rdsdsid.o \
world_avg.o

# Include file locations
INCLUDE	= -I$(INCDIR)/swfinc -I$(LIB3_INC)

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -lgenutils -lmfhdf -ldf -ljpeg -lz

include $(MAKEFILE_APP_TEMPLATE)
