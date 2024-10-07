.POSIX:
	
###########################################################################
# !make
#
# makefile name :     MOD_PR03.mk
#
# !Description: This makefile is used to generate MOD_PR03.exe
#
# !Variables:
#   INC			Include directories
#   LIB			Library directories
#   OBJ			Object files
#   TARGET		Target executable
#
# !Env Variables
#   AC_INCLUDE		Include directory for IMSL functions, set by the
#			iptsetup script (optional)
#   API_LIB		Library of MODIS Application Progra set by the script
#   API_INC		Include directory of MODIS Application Program
#			Interface set by the script
#   CC			Compiler definition for C code set by the script
#   CFLAGS		compiler flags which is set by the script
#   HDFINC		Include directory of HDF 4 set by the script
#   HDF5INC		Include directory of HDF 5 set by the script (optional)
#   HDFLIB		Library of HDF 4 set by the script
#   HDF5LIB		Library of HDF 5 set by the script (optional)
#   HDFEOS_INC		Include directory of HDF-EOS set by the script
#   HDFEOS_LIB		Library of HDF-EOS set by the script
#   LINK_CNL_STATIC	Libraries for IMSL C, set by iptsetup script (optional)
#   PGSINC		Include directory of PGS toolkit set by the script
#   PGSLIB		Library of PGS toolkit set by the script
#   SDST_INC		Include directory of SDST toolkit set by the script
#   SDST_LIB		Library of SDST toolkit set by the script
#
# !Team-unique Header:
#    This software was developed by MODIS Science Data Support Team for
#    the National Aeronautics and Space Administration, Goddard Space
#    Flight Center, under contract NAS5-32373.
#
#
# !Revision History:
#MOD_PR03.mk,v
#Revision 6.7  2015/05/27 21:25:08  jkuyper
#Corrected hdf5 library link.
#
#Revision 6.6  2015/04/14 20:46:43  jkuyper
#Dropped LIB_HDF5OPT, which was only needed when we were linking the shared-object version.
#
#Revision 6.5  2014/08/01 21:45:40  jkuyper
#Corrected POSIX macro name.
#
#Revision 6.4  2011/11/02 21:48:51  kuyper
#Changes to enforce POSIX conformance.
#
#Revision 6.3  2010/04/23 14:46:17  kuyper
#Removed obsolete modules.
#Updated for compatibility with SDP Toolkit 5.2.16.
#
#Revision 6.2  2010/03/24 20:18:14  kuyper
#Resolved Bug 2827 by #defining POSIX_C_SOURCE=1, thereby turning off some of
#  the SGI extensions to C that are turned on by the fact that CFLAGS
#  contains '-xansi' when code is built on our IRIX machines.
#Added GEO_get_ephatt_inputs.o to resolve Bug 2471.
#Changed to link in HDF5 library (which is not actually used) statically,
#  rather than dynamically.
#
# James Kuyper	james.kuyper@sigmaspace.com
#Revision 6.1  2009/05/29 19:11:50  kuyper
#Added new modules for Collection 6.
#
#Revision 5.3  2004/10/14 20:42:49  kuyper
#Explained about Revision keyword, corrected use of it.
#
#Revision 5.2  2004/08/27 21:08:37  kuyper
#Changed spelling of maneuver code file.
#
#Revision 5.1  2004/07/16 23:22:44  kuyper
#Changed to use -ffloat-store option when $(BRAND) is linux.
#
#Revision 4.3  2003/08/27 17:38:25  kuyper
#Removed header.depend stuff
#
#Revision 4.2  2003/08/20 14:52:01  kuyper
#Removed duplicate object file.
#Changed to make inclusion in unit test make files work better.
#
#Revision 4.1  2003/07/27 17:56:44  kuyper
#Added optional support for HDF5.
#Corrected handling of optional IMSL use.
#Added new modules for 4.1.0.
#
#Revision 3.2  2002/09/24 22:55:48  kuyper
#Made IMSL optional.
#Removed chebyshev modules.
#Added GEO_check_ea_headers.o.
#Moved other rules so a simple 'make' command will make $(TARGET).
#
#Revision 3.1  2002/06/12 14:23:51  kuyper
#Dropped GEO_INC.
#Simplified header.depend building.
#Changed to put message header files in current directory.
#
#Revision 2.9  2001/03/19 16:23:11  pliu
#Added rules for building GEO_get_utcpole_metadata and GEO_get_version_metadata.
#
# Revision 2.8  2000/09/05  18:09:42  kuyper
# With version 5.2.6 of toolkit, we need -ljpeg to link to HDF1.1r3.
#
# Revision 2.7  1999/03/12  18:25:35  kuyper
# Alphabetized environment variable list
# Split up INC and LIB lines.
# Increased optimization level to -O3
#
# Revision 2.6  1999/02/02  18:13:47  kuyper
# Replaced smf.c with SDST library.
# Merged GEO_prepare_ancil_data.c and GEO_interp_ephemeris_attitude.c into
#   GEO_ephem_attit.c.
# Added GEO_get_geoid.c
#
# Revision 2.5  1998/03/27  21:02:25  kuyper
# Modified header.* and auto.depend dependencies.
#
# Revision 2.5  1998/03/27  20:24:25  kuyper
# Removed header.sed dependencies; they don't work.
#
# Revision 2.4  1997/11/17  19:18:59  kuyper
# Added LINK_CNL_STATIC to prolog, dropped -lgen
#
# Revision 2.3  1997/10/24  13:36:01  kuyper
# Reinstated automatic dependency checking. Made more robust.
#
# Revision 2.2  1997/10/22  16:02:15  raj
# Removed GEO_get_orbit_node.o, added GEO_create_swath.o, GEO_landsea_mask.o,
#   GKgetpoint.o modules.
#
# Revision 2.1  1997/10/21  18:16:22  kuyper
# Returned from ClearCase
#
# Note:
#            1) must follow MAPI, SDP TK, HDF, and other header/libs order
#            2) HDFLIB must be explicit and in order of -lmfhdf, -ldf,
#               -ljpeg, and -lz.
#            3) system environments, such as API_INC, PGSINC, HDFINC,
#               API_INC, APILIB etc. should be defined before executing
#               this makefile (see Item 3 in README.txt).
# !END
###############################################################################

# Set up include file directories

#Optional include directories:
INC_IMSL = -I$(AC_INCLUDE)
INC_HDF5OPT = -I$(HDF5INC)

# The Revision keyword gets substituted by RCS each time this file is checked
# out. Because the DAAC checks in our code using RCS, or some similar tool, but
# without respecting our RCS revision numbers, this file has been set to be
# automatically checked out with the "-kv" option, removing the Revision
# keyword. Therefore, to check it out with a lock you must use the "-kkv"
# option.
# The double quotes are interpreted by C, making the expansion of the
# MAKEFILE_REVISION macro a string literal. The single quotes prevent the
# double quotes from being removed by the shell. Whew!
INC =	-I$(HDFEOS_INC)\
	-I$(SDST_INC)\
	-I$(HDFINC)\
	-I$(PGSINC)\
	-I$(API_INC)\
	$(INC_IMSL:-I=-DNOIMSL) \
	$(INC_HDF5OPT:-I=) \
	'-DMAKEFILE_REVISION="6.7"'
ADD_CFLAGS = -D_POSIX_C_SOURCE=1 $(BRAND:linux=-ffloat-store)

TARGET = MOD_PR03.exe

# location of the mapi directory, HDF directory, and other library link directories
LIB_DL = $(BRAND:linux=-ldl)$(BRAND:sgi64=-lc)
LIB =	-L$(API_LIB) -lmapi\
	-L$(SDST_LIB) -lsdst\
	-L$(PGSLIB) -lPGSTK \
	-L$(HDFEOS_LIB) -lhdfeos -lGctp -lhdf5 -L$(HDFLIB) -lmfhdf -ldf \
	$${SZIPHOME:+-lsz} $${SZIPHOME:+$(LIB_DL:$(BRAND)=)} -ljpeg -lz -lm

# Instrumentation
OBJI =	GEO_get_inst_mirr_normal.o \
	GEO_get_sample_time.o \
	GEO_get_view_vec.o \
	GEO_interp_mirr_ang.o \
	GEO_interp_mirr_enc.o

# Geolocation
OBJL =	GEO_check_ea_headers.o \
	GEO_earth_location.o \
	GEO_ellip_position.o \
	GEO_ephem_attit.o \
	GEO_get_ephatt_inputs.o \
	GEO_get_T_inst2ecr.o \
	GEO_aggregate.o \
	GEO_hires.o \
	GEO_solar_and_lunar_vectors.o \
	GEO_interp_ECR.o

# Terrain information
OBJT =	GEO_DEM.o \
	GEO_get_geoid.o \
	GEO_landsea_mask.o \
	GEO_terrain_correct.o

# Preparation
OBJN =	GEO_maneuver.o \
	GEO_prepare_l1a_data.o \
       	GEO_prepare_mirr_data.o \
	GEO_read_param_file.o

# Metadata preparation
OBJO =	GEO_initialize_product.o \
	GEO_derived_products.o \
	GEO_update_L1A_metadata.o \
	GEO_get_GRing_points.o \
	GEO_get_bounding_coords.o \
	GEO_create_swath.o \
	GEO_get_utcpole_metadata.o \
	GEO_get_version_metadata.o

# Main process
OBJM = GEO_location_main.o GEO_locate_one_granule.o  GEO_locate_one_scan.o

# Metadata reading
OBJR =	GEO_read_L1Aspecific_metadata.o \
	GEO_read_L1Ascan_metadata.o \
	GEO_read_L1Apacket_data.o \
	GEO_read_L1AECS_metadata.o \
	GEO_read_L1Atemp_data.o

# Utility modules
OBJU =  GEO_mat_vec_mul3.o \
	GEO_poly_fit.o \
	GEO_vec_unit3.o \
	GEO_vec_length3.o \
	GEO_vec_mul3.o \
        GEO_poly_coef1.o \
	GEO_vec_prod3.o

# Data validation
OBJV =	GEO_validate_derived_products.o \
	GEO_validate_earth_location.o \
	GEO_abs_limit_check.o \
	GEO_del_limit_check.o \
	GEO_find_next_flag.o

# Metadata writing
OBJW =	GEO_write_ECS_metadata.o \
	GEO_write_geospecific_metadata.o \
	GEO_write_one_scan.o \
	GEO_write_scan_data.o \
	GEO_write_granule_metadata.o \
	GEO_write_parameters.o \
	GEO_write_scan_metadata.o \
	GEO_write_input_metadata.o


OBJIMSL = \
	imsl_d_lin_sol_gen.o \
	imsl_d_spline_interp.o \
	imsl_d_spline_value.o \
	imsl_error.o

# If LINK_CNL_STATIC is defined normally, the following two lines result in
# removing the '-lm' at the end of it, which is unacceptable in the dependency
# list for $(TARGET). However, if it is not defined, then these two lines result
# in the addition of 'pseudo_imsl.a' to the dependency list.
CNL_OBJ = $(LINK_CNL_STATIC)pseudo_imsl.a
OBJ = $(OBJM) $(OBJI) $(OBJU) $(OBJT) $(OBJL) $(OBJN) $(OBJO) $(OBJR) $(OBJV) $(OBJW) $(CNL_OBJ:-lmpseudo_imsl.a=)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LIB) -o $@

pseudo_imsl.a: $(OBJIMSL)
	$(AR) $(ARFLAGS) $@ $?

clean:
	-rm -f *.o *.exe

PGS_MODIS_35251.h: 35251.t
	smfcompile -f 35251.t -r

# Explicit -o allows for making the target file in the specified directory,
# rather than the current one.
.c.o:
	@echo "\nCompiling" $*
	$(CC) -O3 $(CFLAGS) $(ADD_CFLAGS:$(BRAND)=) $(INC) -c $< -o $@

