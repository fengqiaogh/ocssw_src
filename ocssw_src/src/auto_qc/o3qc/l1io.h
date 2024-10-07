/*******************************************************************

   l1io.h

   purpose: include file for the use of the level 1 direct I/O routines

   Parameters: 
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      unsigned char *   l1a_path        I       path for level 1A file

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       15-Feb-1993     Original development

 *******************************************************************/

/*
 *  Note that hdf.h is needed for this and navigation defs
 */
#include <stdint.h>
#include "nav_l1io.h"
#include "ancil.h"

/*
 *  the l1info_struct structure has file ids...
 */
typedef struct l1info_struct_d {
    int32_t fid;
    int32_t sdfid;
} l1info_struct;

/*
 *  prototypes
 */
void pexit(char *string);
int pwarning(char *string);

int32_t open_hdf(char *fname, l1info_struct *l1info);
int32_t read_g_attr(l1info_struct l1info, char *name, int32_t *n_type,
        int32_t *count, void *data);
int day2mday(int year, int day_of_year, int *month, int *day_of_month);


/* prototypes from ANCroutines.c */
int startHDF(char *outfile, int32_t *sdfid, int32_t *fid, int32_t mode);
int32_t setupGrid(int32_t fid, char *grpname);
int32_t gridToGrid(int32_t outergridid, int32_t innergridid);
int32_t writeGeom(int32_t fid, int32_t gridid, char *geomname, int32_t bin_meth,
        int32_t registration, float vsize, float hsize,
        float max_north, float max_south, float max_west,
        float max_east);
int32_t findGeomId(int32_t fid, char *geomname);
int32_t linkGeom(int32_t gridid, int32_t geomid);
int32_t detachGeom(int32_t geomid);
int addAttr(int32_t sdsid, char *dataattr, int32_t datatype, char *dataunit);
int setSDSref(int32_t sdsid, int32_t gridid);
int deattachHDFgrid(int32_t gridid);
int closeHDFstructs(int32_t sdfid, int32_t fid);
int32_t wrtsds(int32_t sdfid, int rank, int32_t *shape, int32_t datatype,
        char *datalabel, void *data);
int32_t rewrtsds(int32_t sdsid, int32_t *shape, void *data);
int rdsds(char *filename, char *vgname, char *sdsname, int32_t *dimsizes,
        void *inData);
int wrtattr(int32_t dfile, struct annotation *annot, int numannarr);


int32_t l1io_open(char *, l1info_struct*, int32_t *, int32_t *);
int32_t l1io_read(l1info_struct, int, int16 *,
        navblockType *);
void l1io_close(l1info_struct);
