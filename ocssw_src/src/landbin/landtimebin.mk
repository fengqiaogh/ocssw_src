# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

# Makefile name
MAKEFILE= landtimebin.mk

# Progam to make
EXE	= landtimebin

# Object modules for EXE
OBJ	= landtimebin.o 

# Include file locations
INCLUDE	= -I$(LIB3_INC)

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -lmfhdf -ldf -ljpeg -lz

# Program specific additions to compiler flags
LOCAL_CFLAGS = -DHDF41r3

include $(MAKEFILE_APP_TEMPLATE)
