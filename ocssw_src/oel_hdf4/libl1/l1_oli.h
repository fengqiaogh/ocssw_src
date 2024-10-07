#ifndef L1_OLI_H
#define L1_OLI_H

int openl1_oli(filehandle *l1file);
int readl1_oli(filehandle *l1file, int recnum, l1str *l1rec, int lonlat);
int closel1_oli(filehandle *l1file);
int get_oli_nom_angles(char *meta_filename, int32_t npix, int32_t nscan, int32_t iscan,
        float *solz, float *sola, float *senz, float *sena, short *sca_num, short *det_num,
        int32_t isGeomPerBand);
int get_oli_angles( char *emeta_filename, int32_t npix, int32_t nscan, int32_t iscan,
        float *solz, float *sola, float *senz, float *sena,
        int32_t isGeomPerBand);
#endif
