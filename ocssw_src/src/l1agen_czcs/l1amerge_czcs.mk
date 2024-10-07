# ===========================================================================
#	Begin Program Specific Block
# ===========================================================================

MAKEFILE = $(EXE).mk

# Progam to make
EXE	= l1amerge_czcs

# Object modules for EXE
OBJ	=  \
czl1merge.o \
cz_dat_alloc.o \
cz_l1_read.o \
cz_mov_scn.o \
cztimqual.o \
fill_mstr.o \
hdfio.o \
olap_resolve.o \
czcs_l1_write.o \
wrt_czcs_qual.o \
read_file_list.o \
set_czcs_ctl_data.o \
wrt_czcs_sla.o \
create_sds.o \
set_dim_names.o \
cz_meta_adj.o \
cz_ll_upd.o

# Include file locations
INCLUDE	= -I$(INCDIR)/swfinc -I$(INCDIR)/utils -I$(LIB3_INC)

LOCAL_CFLAGS = -D__DEC_MACHINE

# Library locations
LIBS 	= -L$(LIBDIR) -L$(LIB3_LIB)

# Libraries required to link
LLIBS 	= -lczcs -lseawifs -lnav -lgenutils \
          -lmfhdf -ldf -ljpeg -lz 


include $(MAKEFILE_APP_TEMPLATE)
