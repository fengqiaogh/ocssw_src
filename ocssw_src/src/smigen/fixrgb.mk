# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

MAKEFILE = Makefile

# Progam to make
EXE	= fixrgb

# Object modules for EXE
OBJ	=  fixrgb.o

# Include file locations
INCLUDE	= 

# Library locations
LIBS    = -L$(LIB3_LIB)

include $(MAKEFILE_APP_TEMPLATE)

ifdef ppcDarwinArchitecture
  LOCAL_CFLAGS = -O0
endif
