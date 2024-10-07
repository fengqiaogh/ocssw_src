# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

MAKEFILE = $(EXE).mk

# Progam to make
EXE	= l1agen_czcs

# Object modules for EXE
OBJ	=  \
l1czcs.o \
read_crtt.o \
czcs_l1_write.o \
hdr_2_gattr.o \
czcs_ctl_pt.o \
lladjust.o \
cz_ll_upd.o \
wrt_czcs_sla.o \
time_str.o \
wrt_czcs_qual.o \
cz_clean.o \
cz_dat_alloc.o \
cz_sd_set.o \
set_czcs_ctl_data.o \
create_sds.o \
set_dim_names.o \
fill_orb_dat.o \
time_utl.o \
day_to_ofile.o \
rd_smmr_orb.o \
lonlat.o \
hdfio.o \
asap_int2.o

# Include file locations
INCLUDE	= -I$(INCDIR)/swfinc -I$(INCDIR)/utils -I$(LIB3_INC)

LOCAL_CFLAGS = -D__DEC_MACHINE

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -lczcs -lseawifs -lnav -lgenutils \
          -lmfhdf -ldf -ljpeg -lz 


include $(MAKEFILE_APP_TEMPLATE)
#include MakeApp.linux.tpl
