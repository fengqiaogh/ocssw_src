#include "filehandle.h"

#include <string.h>
#include <stdlib.h>

float* calloc_nbandsf(int32_t nbands, float *nbarray, float init_val) {
    int i;

    if ((nbarray = (float *) calloc(nbands, sizeof (float))) == NULL) {
        printf("-E- : Error allocating float memory in alloc_nbandsf\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < nbands; i++)
        nbarray[i] = init_val;

    return nbarray;
}

int32_t* calloc_nbandsi32t(int32_t nbands, int32_t *nbarray, int32_t init_val) {
    int i;
    if ((nbarray = (int32_t *) calloc(nbands, sizeof (int32_t))) == NULL) {
        printf("-E- : Error allocating float memory in alloc_nbandsi\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < nbands; i++)
        nbarray[i] = init_val;

    return nbarray;
}

int* calloc_nbandsi(int32_t nbands, int *nbarray, int init_val) {
    int i;
    if ((nbarray = (int *) calloc(nbands, sizeof (int))) == NULL) {
        printf("-E- : Error allocating float memory in alloc_nbandsi\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < nbands; i++)
        nbarray[i] = init_val;

    return nbarray;
}

void filehandle_init(filehandle *file) {
    int32_t i;

    strcpy(file->name, "");
    file->format = -1;
    file->sensorID = -1;
    file->subsensorID = -1;
    strcpy(file->spatialResolution, "");

    file->length = 0;
    file->spix = 0;
    file->epix = -1;
    file->npix = 0;
    file->nscan = 0;
    file->nbands = 0;
    file->nbandsir = 0;
    file->nlvl = 42;  /* fixed 42 GMAO FP-IT levels now */
    file->n_refl_loc = 10;  /* reflectance location 3rd dim count */
    file->n_cloud_phase = 2;  /* # cloud phases we retrieve */
    file->bindx = NULL;
    file->ndets = 1;
    file->mode = READ;
    strcpy(file->l2prod, "");
    strcpy(file->def_l2prod, "");
    file->sd_id = 0;
    file->tot_prod = 0;
    for (i = 0; i < L1_MAXPROD; i++)
        strcpy(file->l2_prod_names[i], "");
    file->prodptr = NULL;
    file->productInfos = NULL;

    file->geofile = NULL;
    file->orbit_node_lon = -999.0;
    file->orbit_number = 0;
    file->node_crossing_time = 0;
    memset(file->flag_cnt, 0, L1_NFLAGS * sizeof (int32_t));
    file->terrain_corrected = 0;
    file->sv_with_moon = 0;
    for(i=0; i<8; i++)
        file->grp_id[i] = -1;
    
    file->iwave = NULL;
    file->fwave = NULL;
    file->fwhm = NULL;
    file->Fobar = NULL;
    file->Fonom = NULL;
    file->Tau_r = NULL;
    file->k_oz = NULL;
    file->k_no2 = NULL;
    file->aw = NULL;
    file->bbw = NULL;

    file->private_data = NULL;
    
    return;
}
