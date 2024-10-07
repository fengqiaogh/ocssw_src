/* 
 * File:   ancnrt_proto.h
 * Author: dshea
 *
 * Created on February 16, 2011, 11:49 AM
 */

#ifndef ANCNRT_PROTO_H
#define ANCNRT_PROTO_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* fortran  prototypes */
void rdgrid_();
void gregor_(int *julday, int *year, int *month, int *day);

int readgrib2(char *file, int npix, int nlin, int rgmode,
        float *data);
int grib2_t(char *grib2_t_str, int *year, int *doy, int *hour, int *npix,
        int *nlin, int *h_fcst);
int readgrib2_3d(char *file, char *grib2_t_str, int npix, int nlin,
        float *data, int nprfl, int *year, int *month, int *day, int *hour);
int count_annot(char *filename);
void shift_180(float *in, int npix, int nlin, float fact, float *out);


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
int32_t SDSinFile(char *sdsname, char *longname, char *units, char *datafmt,
        int32_t datatype, int32_t sdfid, int32_t rank, int32_t *shape, void *data, int32_t gridid);


void pexit(char *string);
int pwarning(char *string);


#ifdef __cplusplus
}
#endif

#endif /* ANCNRT_PROTO_H */

