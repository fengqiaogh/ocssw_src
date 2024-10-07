#include <stdlib.h>
#include <assert.h>
#include "l12_proto.h"

static int lutLoaded = 0;

// scale factors from lut file
static float scale_551;
static float offset_551;
static float scale_488;
static float offset_488;
static float scale_sst;
static float offset_sst;

// actual lut file
static int32_t lut_num_rows = 0;
static float *lut_nlw_551;
static float *lut_nlw_488;
static float *lut_sst;
static int16_t *lut_class_ward;
static int16_t *lut_class_k;
static int16_t *lut_class_34k_w;

// products for current scan line
static int32_t current_iscan = -1;
static float *class_ward;
static float *class_k;
static float *class_34k_w;

static int ib551 = -1;
static int ib488 = -1;

void init_owmc(l2str *l2rec) {
    char titleStr[256] = "";

    int32 sd_id;
    int32 status;
    int32 dims[H4_MAX_VAR_DIMS];

    int32_t ip;
    int32_t ipb;
    int32_t row;

    float delta;
    float best_delta;
    int32_t best_lut_row;

    float raw_551;
    float raw_488;
    float raw_sst;

    float cooked_551;
    float cooked_488;
    float cooked_sst;

    l1str *l1rec = l2rec->l1rec;

    if (!lutLoaded) {

        printf("Loading OWMC lut file %s.\n", input->owmcfile);

        /* Open the file and initiate the SD interface */
        sd_id = SDstart(input->owmcfile, DFACC_RDONLY);
        if (sd_id == -1) {
            fprintf(stderr, "-E- %s line %d: Could not open OWMC lut file %s.\n",
                    __FILE__, __LINE__, input->owmcfile);
            exit(1);
        }

        // Read the attribute data.
        READ_GLBL_ATTR_E("Title", titleStr);

        // read title OK
        if (strcmp(titleStr, "OWMC LUT") != 0) {
            fprintf(stderr, "-E- %s line %d: Global attribute \"Title\" from %s should be \"OWMC LUT\".\n",
                    __FILE__, __LINE__, input->owmcfile);
            exit(1);
        }

        // Read the attribute data.
        READ_GLBL_ATTR_E("scale_551", &scale_551);
        READ_GLBL_ATTR_E("offset_551", &offset_551);
        READ_GLBL_ATTR_E("scale_488", &scale_488);
        READ_GLBL_ATTR_E("offset_488", &offset_488);
        READ_GLBL_ATTR_E("scale_sst", &scale_sst);
        READ_GLBL_ATTR_E("offset_sst", &offset_sst);

        status = getDims(sd_id, "nlw_551", dims);
        if (status) {
            fprintf(stderr, "-E- %s line %d: Could get dimensions for \"nlw_551\" from, %s.\n",
                    __FILE__, __LINE__, input->owmcfile);
            exit(1);
        }
        lut_num_rows = dims[0];

        // allocate the memory for the lut
        lut_nlw_551 = (float*) malloc(lut_num_rows * sizeof (float));
        assert(lut_nlw_551);
        lut_nlw_488 = (float*) malloc(lut_num_rows * sizeof (float));
        assert(lut_nlw_488);
        lut_sst = (float*) malloc(lut_num_rows * sizeof (float));
        assert(lut_sst);

        lut_class_ward = (int16_t *) malloc(lut_num_rows * sizeof (int16_t));
        assert(lut_class_ward);
        lut_class_k = (int16_t *) malloc(lut_num_rows * sizeof (int16_t));
        assert(lut_class_k);
        lut_class_34k_w = (int16_t *) malloc(lut_num_rows * sizeof (int16_t));
        assert(lut_class_34k_w);

        // allocate memory for the current line
        class_ward = (float *) malloc(l1rec->npix * sizeof (float));
        assert(class_ward);
        class_k = (float *) malloc(l1rec->npix * sizeof (float));
        assert(class_k);
        class_34k_w = (float *) malloc(l1rec->npix * sizeof (float));
        assert(class_34k_w);

        // read in the lut data
        READ_SDS_E("nlw_551", lut_nlw_551, 0, 0, 0, lut_num_rows, 1, 1);
        READ_SDS_E("nlw_488", lut_nlw_488, 0, 0, 0, lut_num_rows, 1, 1);
        READ_SDS_E("sst", lut_sst, 0, 0, 0, lut_num_rows, 1, 1);
        READ_SDS_E("class_ward", lut_class_ward, 0, 0, 0, lut_num_rows, 1, 1);
        READ_SDS_E("class_k", lut_class_k, 0, 0, 0, lut_num_rows, 1, 1);
        READ_SDS_E("class_34k_w", lut_class_34k_w, 0, 0, 0, lut_num_rows, 1, 1);

        /* Close the file */
        status = SDend(sd_id);
        if (status == FAIL) {
            fprintf(stderr, "-E- %s line %d: Could not close OWMC lut file %s.\n",
                    __FILE__, __LINE__, input->owmcfile);
            exit(1);
        }

        // setup the band indexes for the sensor
        if ((ib551 = bindex_get(545)) < 0)
            if ((ib551 = bindex_get(550)) < 0)
                if ((ib551 = bindex_get(555)) < 0)
                    if ((ib551 = bindex_get(560)) < 0) {
                        printf("-E- %s line %d: can't find 551 band\n", __FILE__, __LINE__);
                        exit(1);
                    }

        if ((ib488 = bindex_get(488)) < 0) {
            printf("-E- %s line %d: can't find 488 band\n", __FILE__, __LINE__);
            exit(1);
        }

        lutLoaded = 1;
    }

    // if the current scan line changed
    // load up the products for the current line
    if (l1rec->iscan != current_iscan) {
        current_iscan = l1rec->iscan;

        // loop pixels
        for (ip = 0; ip < l1rec->npix; ip++) {
            ipb = ip * l1rec->l1file->nbands;

            raw_551 = l2rec->nLw[ipb + ib551];
            raw_488 = l2rec->nLw[ipb + ib488];

            // check for missing values
            if ((raw_551 != BAD_FLT) && (raw_551 != BAD_FLT)) {
                if (l2rec->sst && (l2rec->sst[ip] > -2.0))
                    raw_sst = l2rec->sst[ip];
                else
                    raw_sst = l1rec->sstref[ip];

                cooked_551 = (raw_551 - offset_551) / scale_551;
                cooked_488 = (raw_488 - offset_488) / scale_488;
                cooked_sst = (raw_sst - offset_sst) / scale_sst;

                best_delta = fabsf(cooked_551 - lut_nlw_551[0]) +
                        fabsf(cooked_488 - lut_nlw_488[0]) +
                        fabsf(cooked_sst - lut_sst[0]);
                best_lut_row = 0;

                // loop through lut
                for (row = 1; row < lut_num_rows; row++) {
                    delta = fabsf(cooked_551 - lut_nlw_551[row]) +
                            fabsf(cooked_488 - lut_nlw_488[row]) +
                            fabsf(cooked_sst - lut_sst[row]);
                    if (delta < best_delta) {
                        best_delta = delta;
                        best_lut_row = row;
                    }
                } // for lut row

                class_ward[ip] = lut_class_ward[best_lut_row];
                class_k[ip] = lut_class_k[best_lut_row];
                class_34k_w[ip] = lut_class_34k_w[best_lut_row];
            } else {
                class_ward[ip] = BAD_FLT;
                class_k[ip] = BAD_FLT;
                class_34k_w[ip] = BAD_FLT;
            }

        } // for pixel
    } // if scan line changed
}

float *get_class_ward_owmc(l2str *l2rec) {
    init_owmc(l2rec);
    return class_ward;
}

float *get_class_k_owmc(l2str *l2rec) {
    init_owmc(l2rec);
    return class_k;
}

float *get_class_34k_w_owmc(l2str *l2rec) {
    init_owmc(l2rec);
    return class_34k_w;
}
