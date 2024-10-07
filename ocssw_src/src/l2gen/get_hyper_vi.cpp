/**
 * @name get_hyper_vi.c - Vegetation indices using exclusively hyperspectral instruments
 * @brief A collection of scientific algorithms to produce vegetation index products.
 * @authors Jakob Lindo (SSAI)
 * @date Feb 2024,
 */

#include "l12_proto.h"
#include "vegetation_indices.h"

/**
 * @name get_pri
 * @brief Calculate the Photochemical Reflectance Index
 * @cite Grace, J., Nichol, C., Disney, M., Lewis, P., Quaife, T., and Bowyer, P.,  Can we measure terrestrial
 * photosynthesis from space directly, using spectral reflectance and fluorescence? Global Change Biology,
 * 13(7), 1484-1497, https://doi.org/10.1111/j.1365-2486.2007.01352.x, (2007)
 * @param l1rec A level one file
 * @param pri An array into which the result of calculation will be stored
 */
void get_pri(l1str *l1rec, float pri[]) {
    int32_t num_pixels = l1rec->npix;
    int32_t pixel_idx = -0;
    int32_t band_idx = -0;
    double rhos_530 = -0.0;
    double rhos_570 = -0.0;

    static uint idx_530 = bindex_get(530);
    static uint idx_570 = bindex_get(570);
    std::vector<double> rhos_values;

    for (pixel_idx = 0; pixel_idx < num_pixels; pixel_idx++) {
        pri[pixel_idx] = UNDEFINED;
        band_idx = l1rec->l1file->nbands * pixel_idx;

        rhos_530 = l1rec->rhos[band_idx + idx_530];
        rhos_570 = l1rec->rhos[band_idx + idx_570];

        // Pixel attributes and metadata
        if (!rhos_values.empty())
            rhos_values.clear();
        rhos_values.push_back(rhos_530);
        rhos_values.push_back(rhos_570);
        double pixel_elevation = l1rec->dem[pixel_idx];
        double pixel_mask = l1rec->flags[pixel_idx] & LAND_MASK;

        if (invalid_pixel(pixel_elevation, pixel_mask, rhos_values.data(), rhos_values.size())) {
            l1rec->flags[pixel_idx] |= PRODFAIL;
            continue;
        }

        float pri_val = (rhos_530 - rhos_570) / (rhos_530 + rhos_570);
        pri[pixel_idx] = clamp(pri_val, VI_MINVAL, VI_MAXVAL);
    }
}

/**
 * @name get_cire
 * @brief Calculate Chlorophyll Index Red Edge
 * @cite Gitelson, A.A., Viña, A., Rundquist, D.C., Ciganda, V., Arkebauer, T.J., Remote Estimation of Canopy
 * Chlorophyll Content in Crops, Geophysical Research Letters, 32(8), L08403, doi:10.1029/2005GL022688 (2005)
 * @param l1rec A level one file
 * @param cire An array into which the result of calculation will be stored
 */
void get_cire(l1str *l1rec, float cire[]) {
    int32_t num_pixels = l1rec->npix;
    int32_t pixel_idx = -0;
    int32_t band_idx = -0;
    double rhos_800 = -0.0;
    double rhos_705 = -0.0;

    static uint idx_800 = bindex_get(800);
    static uint idx_705 = bindex_get(705);
    std::vector<double> rhos_values;

    for (pixel_idx = 0; pixel_idx < num_pixels; pixel_idx++) {
        cire[pixel_idx] = UNDEFINED;
        band_idx = l1rec->l1file->nbands * pixel_idx;

        rhos_800 = l1rec->rhos[band_idx + idx_800];
        rhos_705 = l1rec->rhos[band_idx + idx_705];

        // Pixel attributes and metadata
        if (!rhos_values.empty())
            rhos_values.clear();

        rhos_values.push_back(rhos_800);
        rhos_values.push_back(rhos_705);
        double pixel_elevation = l1rec->dem[pixel_idx];
        double pixel_mask = l1rec->flags[pixel_idx] & LAND_MASK;

        if (invalid_pixel(pixel_elevation, pixel_mask, rhos_values.data(), rhos_values.size())) {
            l1rec->flags[pixel_idx] |= PRODFAIL;
            continue;
        }

        float cire_val = (rhos_800 / rhos_705) - 1;
        cire[pixel_idx] = clamp(cire_val, VI_MINVAL, VI_MAXVAL);
    }
}

/**
 * @name get_car
 * @brief Calculate a chlorophyll-carotenoid index at 495 705 and 800 nm.
 * @cite Gitelson, A.A., Zur, Y., Chivkunova, O.B., and Merzlyak, M.N., Assessing Carotenoid Content in Plant
 * Leaves with Reflectance Spectroscopy, Photochemistry and Photobiology, 75(3), 272-281,
 * https://doi.org/10.1562/0031-8655(2002)0750272ACCIPL2.0.CO2 (2002)
 * @param l1rec A level one file
 * @param car An array into which the result of calculation will be stored
 */
void get_car(l1str *l1rec, float car[]) {
    int32_t num_pixels = l1rec->npix;
    int32_t pixel_idx = -0;
    int32_t band_idx = -0;
    double rhos_495 = -0.0;
    double rhos_705 = -0.0;
    double rhos_800 = -0.0;

    static uint idx_495 = bindex_get(495);
    static uint idx_705 = bindex_get(705);
    static uint idx_800 = bindex_get(800);
    std::vector<double> rhos_values;

    for (pixel_idx = 0; pixel_idx < num_pixels; pixel_idx++) {
        car[pixel_idx] = UNDEFINED;
        band_idx = l1rec->l1file->nbands * pixel_idx;

        rhos_495 = l1rec->rhos[band_idx + idx_495];
        rhos_705 = l1rec->rhos[band_idx + idx_705];
        rhos_800 = l1rec->rhos[band_idx + idx_800];

        // Pixel attributes and metadata
        if (!rhos_values.empty())
            rhos_values.clear();

        rhos_values.push_back(rhos_495);
        rhos_values.push_back(rhos_705);
        rhos_values.push_back(rhos_800);
        double pixel_elevation = l1rec->dem[pixel_idx];
        double pixel_mask = l1rec->flags[pixel_idx] & LAND_MASK;

        if (invalid_pixel(pixel_elevation, pixel_mask, rhos_values.data(), rhos_values.size())) {
            l1rec->flags[pixel_idx] |= PRODFAIL;
            continue;
        }

        float car_val = ((1 / rhos_495) - (1 / rhos_705)) * rhos_800;
        car[pixel_idx] = clamp(car_val, VI_MINVAL, VI_MAXVAL);
    }
}

/**
 * @name get_ari
 * @brief Calculate the Anthocyanin Reflectance Index
 * @cite Gitelson, A.A., Merzlyak, M.N., and Chivkunova, O.B., Optical properties and non-destructive
 * estimation of anthocyanin content in plant leaves, Photochemistry and Photobiology, 74(1), 38-45,
 * https://doi.org/10.1562/0031-8655(2001)0740038OPANEO2.0.CO2 (2001)
 * @param l1rec A level one file
 * @param ari An array into which the result of calculation will be stored
 */
void get_ari(l1str *l1rec, float ari[]) {
    int32_t num_pixels = l1rec->npix;
    int32_t pixel_idx = -0;
    int32_t band_idx = -0;
    double rhos_550 = -0.0;
    double rhos_705 = -0.0;
    double rhos_800 = -0.0;

    static uint idx_550 = bindex_get(550);
    static uint idx_705 = bindex_get(705);
    static uint idx_800 = bindex_get(800);

    std::vector<double> rhos_values;

    for (pixel_idx = 0; pixel_idx < num_pixels; pixel_idx++) {
        ari[pixel_idx] = UNDEFINED;
        band_idx = l1rec->l1file->nbands * pixel_idx;

        rhos_550 = l1rec->rhos[band_idx + idx_550];
        rhos_705 = l1rec->rhos[band_idx + idx_705];
        rhos_800 = l1rec->rhos[band_idx + idx_800];

        // Pixel attributes and metadata
        if (!rhos_values.empty())
            rhos_values.clear();

        rhos_values.push_back(rhos_550);
        rhos_values.push_back(rhos_705);
        rhos_values.push_back(rhos_800);
        double pixel_elevation = l1rec->dem[pixel_idx];
        double pixel_mask = l1rec->flags[pixel_idx] & LAND_MASK;

        if (invalid_pixel(pixel_elevation, pixel_mask, rhos_values.data(), rhos_values.size())) {
            l1rec->flags[pixel_idx] |= PRODFAIL;
            continue;
        }

        float ari_val = ((1 / rhos_550) - (1 / rhos_705)) * rhos_800;
        ari[pixel_idx] = clamp(ari_val, VI_MINVAL, VI_MAXVAL);
    }
}

/**
 * @name get_hyper_vi(3) - Main entry point for getting hyperspectral vegetation indices
 * @brief Takes in a level 1 file, product number/catalogue index, and an array in which the result will be
 * stored.
 * @param l1rec A level one file.
 * @param product_number A catalogue index indicating which hyperspectral VI the caller wants.
 * @param product The array, passed by reference, in which the resultant product will be stored.
 * @return 0 if successful, non-zero if unsuccessful.
 */
int get_hyper_vi(l1str *l1rec, int product_number, float product[]) {
    switch (product_number) {
        case CAT_pri:
            get_pri(l1rec, product);
            break;
        case CAT_cire:
            get_cire(l1rec, product);
            break;
        case CAT_car:
            get_car(l1rec, product);
            break;
        case CAT_ari:
            get_ari(l1rec, product);
            break;
        default:
            printf("Error in %s; Unknown product specifier: %d\n", __FILE__, product_number);
            exit(FATAL_ERROR);
    }
    return (EXIT_SUCCESS);
}