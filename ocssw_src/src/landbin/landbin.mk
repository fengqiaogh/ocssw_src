# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

# Makefile name
MAKEFILE= landbin.mk

# Progam to make
EXE	= landbin

# Object modules for EXE
OBJ	= make_L3_v1.1.o \
          calib_calibrate_l1a.o \
	  calib_get_cal.o \
	  calib_get_cal_misc.o \
	  chand.o \
	  compute_l1b.o \
	  csalbr.o \
	  get_attributes.o \
	  get_calib_sds.o \
	  get_navig_sds.o \
	  make_psurf.o \
	  proj_cproj.o \
	  proj_hamfor.o \
	  proj_haminv.o \
	  proj_molwfor.o \
	  proj_molwinv.o \
	  proj_report.o \
	  proj_robfor.o \
	  proj_robinv.o \
	  read_write.o

# Include file locations
INCLUDE	= -I$(LIB3_INC)

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -lmfhdf -ldf -ljpeg -lz

# Program specific additions to compiler flags
LOCAL_CFLAGS = -DHDF41r3

include $(MAKEFILE_APP_TEMPLATE)

