/******************************************************************************
 *
 *			Photosynthetic Available Radiation Utilities
 *				Sensor-Independent Shared Subroutine
 *					by Robert J. Lossing
 *						July 30, 2013
 *				Ocean Color Processing Group
 *
 *******************************************************************************/

/******************************************************************************
 *
 *      Retrieve Phase function and Omega, given input scattering angle
 *        and Angstrom coefficient.
 *
 *      Authors: Robert Frouin <RFrouin@ucsd.edu>, scientific algorithms
 *              John McPherson <JMcPherson@ucsd.edu>, program structure
 *
 *      Note: Angstrom exponents in tables aren't in increasing order,
 *        need to determine proper order to use:
 *   			Model  omega(6)       alpha(6,8)     true_order
 *      			1   1.0000000E+00 -8.9877717E-02  2nd
 *     			2   9.9999797E-01 -9.1842167E-02  1st *min
 *     			3   9.8501903E-01  5.1216149E-01  8
 *      			4   9.8840600E-01  4.0654629E-01  7
 *     			5   9.9607497E-01  1.9749035E-01  4
 *      			6   9.9880499E-01  8.1851527E-02  3rd
 *      			7   9.7776997E-01  7.7916741E-01  9
 *     			8   9.9351001E-01  3.9742991E-01  6
 *    			9   9.9790698E-01  2.1781628E-01  5
 *     			10  9.5734900E-01  1.6471386E+00  12th *max
 *     			11  9.8160601E-01  1.4880587E+00  11th
 *     			12  9.9180502E-01  1.2930322E+00  10th
 *
 *******************************************************************************/

#include "par_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NMODEL 12
#define NS 15
#define NT 9

static float *tablewavelengths;
static float *tablephaseangles;
static float *tablealphas;  // angstrom exponent
static float *tableomegas;  // single scattering albedo
static float *tableaerphasefunc;

#define FILE_ERROR(str, name)  {        if (str.fid == -1) { \
printf("-Error--: could not open \"%s\" for reading.\n, See line %d in %s. Exiting...", name, __LINE__, __FILE__); \
exit(EXIT_FAILURE);\
}}

inline float kasten_equation(float solz) {
    float CosSZ = cos(solz * OEL_PI / 180.0);
    return (1.003198 * CosSZ + 0.101632) /
           (CosSZ * CosSZ + 0.090560 * CosSZ + 0.003198);
}

size_t search(const float *arr, size_t s, size_t e, float val, size_t *i_s,
              size_t *i_e) {
    const bool acsending = arr[s] < arr[e];
    if (acsending) {
        if (val >= arr[e]) {
            *i_s = e;
            *i_e = e;
            return e;
        }
        if (val <= arr[s]) {
            *i_s = s;
            *i_e = s;
            return s;
        }
    } else {
        if (val <= arr[e]) {
            *i_s = e;
            *i_e = e;
            return e;
        }
        if (val >= arr[s]) {
            *i_s = s;
            *i_e = s;
            return s;
        }
    }
    while (e - s > 1) {
        size_t m = (s + e) / 2;  // compute median
        if (acsending) {
            if (arr[m] <= val) {
                s = m;
            } else {
                e = m;
            }
        } else {
            if (arr[m] <= val) {
                e = m;
            } else {
                s = m;
            }
        }
    }
    {
        *i_s = s;
        *i_e = s + 1;
    }
    return s;
}

size_t get_reduced_product(size_t n_dims, const size_t *dims) {
    size_t result = 1;
    for (size_t i = 0; i < n_dims; i++) {
        result *= dims[i];
    }
    return result;
}

size_t get_hyper_cube_index(size_t n_dims, const size_t *dims,
                            const size_t *indexes) {
    if (n_dims > 1)
        return indexes[0] * get_reduced_product(n_dims - 1, dims + 1) +
               get_hyper_cube_index(n_dims - 1, dims + 1, indexes + 1);
    else
        return indexes[0];
}

// 2^n_dims, float width, float diff
float hyper_cube_linear_interp(size_t n_dims, float *hyper_cube, float *width) {
    if (n_dims == 1)  // hyper_cube[2]
    {
        float ans = (hyper_cube[1] - hyper_cube[0]) * width[0] + hyper_cube[0];
        return ans;
    } else {
        size_t total_index = pow(2, n_dims - 1);
        float *small_hypercube = (float *) calloc(total_index, sizeof(float));
        for (size_t ind = 0; ind < total_index; ind++) {
            small_hypercube[ind] =
                    (hyper_cube[ind + total_index] - hyper_cube[ind]) * width[0] +
                    hyper_cube[ind];
            // {
            //     std::cout << "small hypercube " << n_dims << " " << ind << "
            //     "
            //               << small_hypercube[ind] << " "
            //               << hyper_cube[ind + total_index] << " "
            //               << hyper_cube[ind] << " " << width[0] << std::endl;
            // }
        }
        float ans =
                hyper_cube_linear_interp(n_dims - 1, small_hypercube, width + 1);
        free(small_hypercube);
        return ans;
    }
}

float interpnd(size_t n_dims, const size_t *dims, const float *point,
               float **grid, const float *lut) {
    size_t *st_pts = (size_t *) calloc(n_dims, sizeof(size_t));
    size_t *end_pts = (size_t *) calloc(n_dims, sizeof(size_t));;
    for (size_t i_dim = 0; i_dim < n_dims; i_dim++) {
        size_t dim_size = dims[i_dim];
        size_t st = 0;
        size_t end = dim_size - 1;
        search(grid[i_dim], st, end, point[i_dim], &st, &end);
        st_pts[i_dim] = st;
        end_pts[i_dim] = end;
    }
    size_t total_index = pow(2, n_dims);
    float *hypercube = (float *) calloc(total_index, sizeof(float));
    float *width = (float *) calloc(n_dims, sizeof(float));
    size_t *indexes = (size_t *) calloc(n_dims, sizeof(size_t));
    for (size_t i_dim = 0; i_dim < n_dims; i_dim++) {
        size_t st = st_pts[i_dim];
        size_t end = end_pts[i_dim];
        if (st != end)
            width[i_dim] = (point[i_dim] - grid[i_dim][st]) /
                           (grid[i_dim][end] - grid[i_dim][st]);
        else
            width[i_dim] = 0.0;
        // { std::cout << "width[i_dim]" << width[i_dim] << std::endl; }
    }

    for (size_t i = 0; i < total_index; i++) {
        size_t hyper_cube_index = 0;
        for (size_t i_dim = 0; i_dim < n_dims; i_dim++) {
            size_t bit = (i >> i_dim) & 1;
            if (bit) {
                indexes[i_dim] = end_pts[i_dim];
                hyper_cube_index += pow(2, n_dims - i_dim - 1);
            } else {
                indexes[i_dim] = st_pts[i_dim];
            }
            // {
            //     std::cout << "index " << i_dim << " i_dim " << indexes[i_dim]
            //               << std::endl;
            // }
        }
        size_t index_lut = get_hyper_cube_index(n_dims, dims, indexes);
        hypercube[hyper_cube_index] = lut[index_lut];
        // {
        //     std::cout << i << " " << hyper_cube_index << " " << index_lut <<
        //     " "
        //               << hypercube[hyper_cube_index] << std::endl;
        // }
    }
    float ans = hyper_cube_linear_interp(n_dims, hypercube, width);
    free(st_pts);
    free(end_pts);
    free(hypercube);
    free(width);
    free(indexes);
    return ans;
}

void get_nc_dim(const idDS *struct_nc, const char *name, size_t *len) {
    int dim_id;
    int status = nc_inq_dimid((*struct_nc).fid, name, &dim_id);
    if (status == 0) status = nc_inq_dimlen((*struct_nc).fid, dim_id, len);
    if (status != 0) {
        printf("--Error--:  Dimension %s could not be read. See line %d in %s. Exiting...", name, __LINE__, __FILE__);
        exit(EXIT_FAILURE);
    }
}

void get_nc_var(const idDS *struct_nc, const char *name, float *var,
                const size_t *start, const size_t *count) {
    int var_id;
    int status = nc_inq_varid((*struct_nc).fid, name, &var_id);
    if (status == 0)
        status = nc_get_vara_float((*struct_nc).fid, var_id, start, count, var);
    if (status != 0) {
        printf("--Error--: Variable %s could not be read. See line %d in %s. Exiting...", name, __LINE__, __FILE__);
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief Get the luts data object
 * For now it includes the MERRA data. Should be moved after the merra data is
 * read.
 */
void get_luts_data(l2str *l2rec, luts_par *luts_data) {
    const char *dataroot = getenv("OCDATAROOT");
    char lut_rho_path[FILENAME_MAX] = "";
    char lut_tg_path[FILENAME_MAX] = "";
    char lut_td_path[FILENAME_MAX] = "";
    char lut_dobson_path[FILENAME_MAX] = "";
    char lut_watvap_path[FILENAME_MAX] = "";
    char lut_surfoceanalbedo[FILENAME_MAX] = "";
    char luts_scalar_par_path[FILENAME_MAX] = "";
    char luts_scalar_inst_par_path[FILENAME_MAX] = "";
    sprintf(lut_rho_path, "%s/common/LUT_rho.nc", dataroot);
    sprintf(lut_td_path, "%s/common/LUT_td.nc", dataroot);
    sprintf(lut_tg_path, "%s/common/LUT_tg.nc", dataroot);
    sprintf(lut_dobson_path, "%s/common/LUT_dobson.nc", dataroot);
    sprintf(lut_watvap_path, "%s/common/LUT_watvap.nc", dataroot);
    sprintf(lut_surfoceanalbedo, "%s/common/LUT_ocean_surface_albedo.nc",
            dataroot);
    sprintf(luts_scalar_par_path, "%s/common/LUT_scalar_PAR.nc", dataroot);
    sprintf(luts_scalar_inst_par_path, "%s/common/LUT_scalar_inst_PAR.nc", dataroot);
    int16_t year, month, mday, doy;
    double sec;
    unix2yds(l2rec->l1rec->scantime, &year, &doy, &sec);
    unix2ymds(l2rec->l1rec->scantime, &year, &month, &mday, &sec);
    // Reading rho lut
    {
        idDS luts_rho_struct = openDS(lut_rho_path);
        FILE_ERROR(luts_rho_struct, lut_rho_path)
        size_t solar_zenith_angle, view_zenith_angle, relative_azimuth_angle,
                angstrom_coefficients, optical_depth_at_550_nm, wavelength;
        get_nc_dim(&luts_rho_struct, "solar_zenith_angle", &solar_zenith_angle);
        get_nc_dim(&luts_rho_struct, "view_zenith_angle", &view_zenith_angle);
        get_nc_dim(&luts_rho_struct, "relative_azimuth_angle",
                   &relative_azimuth_angle);
        get_nc_dim(&luts_rho_struct, "angstrom_coefficients",
                   &angstrom_coefficients);
        get_nc_dim(&luts_rho_struct, "optical_depth_at_550_nm",
                   &optical_depth_at_550_nm);
        get_nc_dim(&luts_rho_struct, "wavelength", &wavelength);
        size_t start[] = {0, 0, 0, 0, 0, 0};
        size_t count[] = {solar_zenith_angle, view_zenith_angle,
                          relative_azimuth_angle, angstrom_coefficients,
                          optical_depth_at_550_nm, wavelength};
        const size_t size_all = solar_zenith_angle * view_zenith_angle *
                                relative_azimuth_angle * angstrom_coefficients *
                                optical_depth_at_550_nm * wavelength;
        luts_data->lut_rho = calloc(size_all, sizeof(float));

        luts_data->rhodims.dimsolar_zenith_angle = solar_zenith_angle;
        luts_data->rhodims.dimview_zenith_angle = view_zenith_angle;
        luts_data->rhodims.dimrelative_azimuth_angle = relative_azimuth_angle;
        luts_data->rhodims.dimangstrom_coefficients = angstrom_coefficients;
        luts_data->rhodims.dimoptical_depth_at_550_nm = optical_depth_at_550_nm;
        luts_data->rhodims.dimwavelength = wavelength;

        luts_data->rhodims.solar_zenith_angle =
                calloc(solar_zenith_angle, sizeof(float));
        luts_data->rhodims.view_zenith_angle =
                calloc(view_zenith_angle, sizeof(float));
        luts_data->rhodims.relative_azimuth_angle =
                calloc(relative_azimuth_angle, sizeof(float));
        luts_data->rhodims.angstrom_coefficients =
                calloc(angstrom_coefficients, sizeof(float));
        luts_data->rhodims.optical_depth_at_550_nm =
                calloc(optical_depth_at_550_nm, sizeof(float));
        luts_data->rhodims.wavelength = calloc(wavelength, sizeof(float));
        get_nc_var(&luts_rho_struct, "lut_rho", luts_data->lut_rho, start,
                   count);
        size_t start_dim[] = {0};
        size_t count_dim[] = {0};
        count_dim[0] = solar_zenith_angle;
        get_nc_var(&luts_rho_struct, "solar_zenith_angle",
                   luts_data->rhodims.solar_zenith_angle, start_dim, count_dim);

        count_dim[0] = view_zenith_angle;
        get_nc_var(&luts_rho_struct, "view_zenith_angle",
                   luts_data->rhodims.view_zenith_angle, start_dim, count_dim);

        count_dim[0] = relative_azimuth_angle;
        get_nc_var(&luts_rho_struct, "relative_azimuth_angle",
                   luts_data->rhodims.relative_azimuth_angle, start_dim,
                   count_dim);

        count_dim[0] = angstrom_coefficients;
        get_nc_var(&luts_rho_struct, "angstrom_coefficients",
                   luts_data->rhodims.angstrom_coefficients, start_dim,
                   count_dim);

        count_dim[0] = optical_depth_at_550_nm;
        get_nc_var(&luts_rho_struct, "optical_depth_at_550_nm",
                   luts_data->rhodims.optical_depth_at_550_nm, start_dim,
                   count_dim);

        count_dim[0] = wavelength;
        get_nc_var(&luts_rho_struct, "wavelength",
                   luts_data->rhodims.wavelength, start_dim, count_dim);

        endDS(luts_rho_struct);
    }
    // Reading tg lut
    {
        idDS luts_tg_struct = openDS(lut_tg_path);
        FILE_ERROR(luts_tg_struct, lut_tg_path)
        size_t wavelength, air_mass, ozone_concentration, water_vapor_pressure;
        get_nc_dim(&luts_tg_struct, "wavelength", &wavelength);
        get_nc_dim(&luts_tg_struct, "air_mass", &air_mass);
        get_nc_dim(&luts_tg_struct, "ozone_concentration",
                   &ozone_concentration);
        get_nc_dim(&luts_tg_struct, "water_vapor_pressure",
                   &water_vapor_pressure);
        size_t start[] = {0, 0, 0, 0};
        size_t count[] = {wavelength, air_mass, water_vapor_pressure,
                          ozone_concentration};
        const size_t size_all =
                wavelength * air_mass * ozone_concentration * water_vapor_pressure;
        luts_data->tgdims.dimair_mass = air_mass;
        luts_data->tgdims.dimozone_concentration = ozone_concentration;
        luts_data->tgdims.dimwavelength = wavelength;
        luts_data->tgdims.dimwater_vapor_pressure = water_vapor_pressure;
        luts_data->lut_tg = calloc(size_all, sizeof(float));
        luts_data->tgdims.air_mass = calloc(air_mass, sizeof(float));
        luts_data->tgdims.wavelength = calloc(wavelength, sizeof(float));
        luts_data->tgdims.ozone_concentration =
                calloc(ozone_concentration, sizeof(float));
        luts_data->tgdims.water_vapor_pressure =
                calloc(water_vapor_pressure, sizeof(float));

        get_nc_var(&luts_tg_struct, "lut_tg", luts_data->lut_tg, start, count);

        size_t start_dim[] = {0};
        size_t count_dim[] = {0};

        count_dim[0] = wavelength;
        get_nc_var(&luts_tg_struct, "wavelength", luts_data->tgdims.wavelength,
                   start_dim, count_dim);

        count_dim[0] = water_vapor_pressure;
        get_nc_var(&luts_tg_struct, "water_vapor_pressure",
                   luts_data->tgdims.water_vapor_pressure, start_dim,
                   count_dim);

        count_dim[0] = ozone_concentration;
        get_nc_var(&luts_tg_struct, "ozone_concentration",
                   luts_data->tgdims.ozone_concentration, start_dim, count_dim);

        count_dim[0] = air_mass;
        get_nc_var(&luts_tg_struct, "air_mass", luts_data->tgdims.air_mass,
                   start_dim, count_dim);
        endDS(luts_tg_struct);
    }
    // Reading td lut
    {
        idDS luts_td_struct = openDS(lut_td_path);
        FILE_ERROR(luts_td_struct, lut_td_path)
        size_t wavelength, air_mass, angstrom_coefficients,
                optical_depth_at_550_nm;
        get_nc_dim(&luts_td_struct, "wavelength", &wavelength);
        get_nc_dim(&luts_td_struct, "air_mass", &air_mass);
        get_nc_dim(&luts_td_struct, "angstrom_coefficients",
                   &angstrom_coefficients);
        get_nc_dim(&luts_td_struct, "optical_depth_at_550_nm",
                   &optical_depth_at_550_nm);

        size_t start[] = {0, 0, 0, 0};
        size_t count[] = {wavelength, air_mass, angstrom_coefficients,
                          optical_depth_at_550_nm};
        const size_t size_all = wavelength * air_mass * angstrom_coefficients *
                                optical_depth_at_550_nm;

        luts_data->tddims.dimair_mass = air_mass;
        luts_data->tddims.dimoptical_depth_at_550_nm = optical_depth_at_550_nm;
        luts_data->tddims.dimwavelength = wavelength;
        luts_data->tddims.dimangstrom_coefficients = angstrom_coefficients;
        luts_data->lut_td = calloc(size_all, sizeof(float));

        luts_data->tddims.air_mass = calloc(air_mass, sizeof(float));
        luts_data->tddims.wavelength = calloc(wavelength, sizeof(float));
        luts_data->tddims.angstrom_coefficients =
                calloc(angstrom_coefficients, sizeof(float));
        luts_data->tddims.optical_depth_at_550_nm =
                calloc(optical_depth_at_550_nm, sizeof(float));
        get_nc_var(&luts_td_struct, "lut_td", luts_data->lut_td, start, count);

        size_t start_dim[] = {0};
        size_t count_dim[] = {0};

        count_dim[0] = wavelength;
        get_nc_var(&luts_td_struct, "wavelength", luts_data->tddims.wavelength,
                   start_dim, count_dim);

        count_dim[0] = angstrom_coefficients;
        get_nc_var(&luts_td_struct, "angstrom_coefficients",
                   luts_data->tddims.angstrom_coefficients, start_dim,
                   count_dim);

        count_dim[0] = optical_depth_at_550_nm;
        get_nc_var(&luts_td_struct, "optical_depth_at_550_nm",
                   luts_data->tddims.optical_depth_at_550_nm, start_dim,
                   count_dim);

        count_dim[0] = air_mass;
        get_nc_var(&luts_td_struct, "air_mass", luts_data->tddims.air_mass,
                   start_dim, count_dim);
        endDS(luts_td_struct);
    }

    // Reading Ozone Conc. Lut
    {
        idDS luts_dobson_struct = openDS(lut_dobson_path);
        FILE_ERROR(luts_dobson_struct, lut_dobson_path)
        size_t days, latitude;
        get_nc_dim(&luts_dobson_struct, "days", &days);
        get_nc_dim(&luts_dobson_struct, "latitude", &latitude);
        luts_data->ozonedims.dimlatitude = latitude;
        luts_data->ozonedims.dimdays = days;
        size_t start[] = {
                0,
                0,
        };
        size_t count[] = {days, latitude};
        const size_t size_all = days * latitude;
        luts_data->lut_dobson = calloc(size_all, sizeof(float));
        luts_data->ozonedims.latitude = calloc(latitude, sizeof(float));
        luts_data->ozonedims.days = calloc(days, sizeof(float));
        get_nc_var(&luts_dobson_struct, "Dobson", luts_data->lut_dobson, start,
                   count);
        size_t start_dim[] = {0};
        size_t count_dim[] = {0};
        count_dim[0] = latitude;
        get_nc_var(&luts_dobson_struct, "latitude",
                   luts_data->ozonedims.latitude, start_dim, count_dim);
        count_dim[0] = days;
        if (year % 4 == 0)
            get_nc_var(&luts_dobson_struct, "days_leap",
                       luts_data->ozonedims.days, start_dim, count_dim);
        else
            get_nc_var(&luts_dobson_struct, "days_no_leap",
                       luts_data->ozonedims.days, start_dim, count_dim);
        endDS(luts_dobson_struct);
    }
    // Reading Wat. Vapor. Lut.
    {
        idDS luts_watvap_struct = openDS(lut_watvap_path);
        FILE_ERROR(luts_watvap_struct, lut_watvap_path)
        size_t days, latitude;
        get_nc_dim(&luts_watvap_struct, "days", &days);
        get_nc_dim(&luts_watvap_struct, "latitude", &latitude);
        luts_data->watvapdims.dimlatitude = latitude;
        luts_data->watvapdims.dimdays = days;
        size_t start[] = {
                0,
                0,
        };
        size_t count[] = {days, latitude};
        const size_t size_all = days * latitude;
        luts_data->lut_watvap = calloc(size_all, sizeof(float));
        luts_data->watvapdims.latitude = calloc(latitude, sizeof(float));
        luts_data->watvapdims.days = calloc(days, sizeof(float));
        get_nc_var(&luts_watvap_struct, "WatVap", luts_data->lut_watvap, start,
                   count);
        size_t start_dim[] = {0};
        size_t count_dim[] = {0};
        count_dim[0] = latitude;
        get_nc_var(&luts_watvap_struct, "latitude",
                   luts_data->watvapdims.latitude, start_dim, count_dim);
        count_dim[0] = days;
        if (year % 4 == 0)
            get_nc_var(&luts_watvap_struct, "days_leap",
                       luts_data->watvapdims.days, start_dim, count_dim);
        else
            get_nc_var(&luts_watvap_struct, "days_no_leap",
                       luts_data->watvapdims.days, start_dim, count_dim);
        endDS(luts_watvap_struct);
    }
    // Reading Surface Ocean Albedo Lut
    {
        idDS luts_soa_struct = openDS(lut_surfoceanalbedo);
        FILE_ERROR(luts_soa_struct, lut_surfoceanalbedo)
        size_t dim_jwvl, dim_wv1, dim_wv3, dim_wlt_reft;
        get_nc_dim(&luts_soa_struct, "dim_jwvl", &dim_jwvl);
        get_nc_dim(&luts_soa_struct, "dim_wv1", &dim_wv1);
        get_nc_dim(&luts_soa_struct, "dim_wv3", &dim_wv3);
        get_nc_dim(&luts_soa_struct, "dim_wlt_reft", &dim_wlt_reft);

        luts_data->soa_lut.dim_jwvl = dim_jwvl;
        luts_data->soa_lut.dim_wv1 = dim_wv1;
        luts_data->soa_lut.dim_wv3 = dim_wv3;
        luts_data->soa_lut.dim_wlt_reft = dim_wlt_reft;

        luts_data->soa_lut.jwl = (float *) calloc(dim_jwvl, sizeof(float));
        luts_data->soa_lut.achl = (float *) calloc(dim_jwvl, sizeof(float));

        luts_data->soa_lut.wv1 = (float *) calloc(dim_wv1, sizeof(float));
        luts_data->soa_lut.bw1 = (float *) calloc(dim_wv1, sizeof(float));

        luts_data->soa_lut.wv3 = (float *) calloc(dim_wv3, sizeof(float));
        luts_data->soa_lut.aw3 = (float *) calloc(dim_wv3, sizeof(float));

        luts_data->soa_lut.wlt_reft =
                (float *) calloc(dim_wlt_reft, sizeof(float));
        luts_data->soa_lut.reft = (float *) calloc(dim_wlt_reft, sizeof(float));

        size_t start_dim[] = {0};
        size_t count_dim[] = {0};

        count_dim[0] = dim_jwvl;
        get_nc_var(&luts_soa_struct, "jwvl", luts_data->soa_lut.jwl, start_dim,
                   count_dim);
        get_nc_var(&luts_soa_struct, "achl", luts_data->soa_lut.achl, start_dim,
                   count_dim);

        count_dim[0] = dim_wv1;
        get_nc_var(&luts_soa_struct, "wv1", luts_data->soa_lut.wv1, start_dim,
                   count_dim);
        get_nc_var(&luts_soa_struct, "bw1", luts_data->soa_lut.bw1, start_dim,
                   count_dim);

        count_dim[0] = dim_wv3;
        get_nc_var(&luts_soa_struct, "wv3", luts_data->soa_lut.wv3, start_dim,
                   count_dim);
        get_nc_var(&luts_soa_struct, "aw3", luts_data->soa_lut.aw3, start_dim,
                   count_dim);

        count_dim[0] = dim_wlt_reft;
        get_nc_var(&luts_soa_struct, "wlt_reft", luts_data->soa_lut.wlt_reft,
                   start_dim, count_dim);
        get_nc_var(&luts_soa_struct, "reft", luts_data->soa_lut.reft, start_dim,
                   count_dim);
        endDS(luts_soa_struct);
    }
    // Reading Scalar Par and Mean Cosine Luts
    {
        idDS luts_scalar_par_struct = openDS(luts_scalar_par_path);
        FILE_ERROR(luts_scalar_par_struct, luts_scalar_par_path)
        size_t dim_wind_speed, dim_doy, dim_latitude;
        get_nc_dim(&luts_scalar_par_struct, "wind_speed", &dim_wind_speed);
        get_nc_dim(&luts_scalar_par_struct, "doy", &dim_doy);
        get_nc_dim(&luts_scalar_par_struct, "latitude", &dim_latitude);
        luts_data->scalar_luts.dim_latitude = dim_latitude;
        luts_data->scalar_luts.dim_wind_speed = dim_wind_speed;
        luts_data->scalar_luts.dim_doy = dim_doy;
        luts_data->scalar_luts.latitude =
                (float *) calloc(dim_latitude, sizeof(float));
        luts_data->scalar_luts.doy = (float *) calloc(dim_doy, sizeof(float));
        luts_data->scalar_luts.wind_speed =
                (float *) calloc(dim_wind_speed, sizeof(float));
        size_t start_dim[] = {0};
        size_t count_dim[] = {0};
        count_dim[0] = dim_latitude;
        get_nc_var(&luts_scalar_par_struct, "latitude",
                   luts_data->scalar_luts.latitude, start_dim, count_dim);
        count_dim[0] = dim_wind_speed;
        get_nc_var(&luts_scalar_par_struct, "wind_speed",
                   luts_data->scalar_luts.wind_speed, start_dim, count_dim);
        count_dim[0] = dim_doy;
        get_nc_var(&luts_scalar_par_struct, "doy", luts_data->scalar_luts.doy,
                   start_dim, count_dim);

        size_t start[] = {0, 0, 0};
        size_t count[] = {dim_wind_speed, dim_doy, dim_latitude};
        const size_t size_all = dim_wind_speed * dim_doy * dim_latitude;
        luts_data->scalar_luts.CF = (float *) calloc(size_all, sizeof(float));
        luts_data->scalar_luts.mu_cosine =
                (float *) calloc(size_all, sizeof(float));
        luts_data->scalar_luts.mu_cosine_c =
                (float *) calloc(size_all, sizeof(float));
        luts_data->scalar_luts.PARo = (float *) calloc(size_all, sizeof(float));
        luts_data->scalar_luts.PARo_c =
                (float *) calloc(size_all, sizeof(float));
        luts_data->scalar_luts.PARc_above =
                (float *) calloc(size_all, sizeof(float));
        get_nc_var(&luts_scalar_par_struct, "CF", luts_data->scalar_luts.CF,
                   start, count);
        get_nc_var(&luts_scalar_par_struct, "mu_cosine",
                   luts_data->scalar_luts.mu_cosine, start, count);
        get_nc_var(&luts_scalar_par_struct, "mu_cosine_c",
                   luts_data->scalar_luts.mu_cosine_c, start, count);
        get_nc_var(&luts_scalar_par_struct, "PARo", luts_data->scalar_luts.PARo,
                   start, count);
        get_nc_var(&luts_scalar_par_struct, "PARo_c",
                   luts_data->scalar_luts.PARo_c, start, count);
        get_nc_var(&luts_scalar_par_struct, "PARc_above",
                   luts_data->scalar_luts.PARc_above, start, count);
                   
    }
    // Reading Scalar Inst Par LUT
    {
        idDS luts_scalar_par_inst_struct = openDS(luts_scalar_inst_par_path);
        FILE_ERROR(luts_scalar_par_inst_struct, luts_scalar_inst_par_path)
        size_t dim_wind_speed, dim_solz, dim_cot;
        get_nc_dim(&luts_scalar_par_inst_struct, "wind_speed", &dim_wind_speed);
        get_nc_dim(&luts_scalar_par_inst_struct, "solz", &dim_solz);
        get_nc_dim(&luts_scalar_par_inst_struct, "cot", &dim_cot); // cot is for cf grid
        luts_data->scalar_inst_luts.dim_wind_speed = dim_wind_speed;
        luts_data->scalar_inst_luts.dim_solz = dim_solz;
        luts_data->scalar_inst_luts.dim_cot = dim_cot;
        luts_data->scalar_inst_luts.wind_speed = (float *) calloc(dim_wind_speed, sizeof(float));
        luts_data->scalar_inst_luts.cos_solz = (float *) calloc(dim_solz, sizeof(float));
        luts_data->scalar_inst_luts.cot = (float *) calloc(dim_cot, sizeof(float));
        luts_data->scalar_inst_luts.pard_m_cs = (float *) calloc(dim_wind_speed * dim_solz, sizeof(float));
        luts_data->scalar_inst_luts.pard_p_cs = (float *) calloc(dim_wind_speed * dim_solz, sizeof(float));
        luts_data->scalar_inst_luts.pard_m_oc = (float *) calloc(dim_wind_speed * dim_solz, sizeof(float));
        luts_data->scalar_inst_luts.cf_pard_m = (float *) calloc(dim_wind_speed * dim_solz * dim_cot, sizeof(float));
        luts_data->scalar_inst_luts.cf_pard_p = (float *) calloc(dim_wind_speed * dim_solz * dim_cot, sizeof(float));

        // reading 1d vars
        {
            size_t start[] = {0};
            size_t count[] = {dim_wind_speed};
            get_nc_var(&luts_scalar_par_inst_struct, "wind_speed", luts_data->scalar_inst_luts.wind_speed,
                       start, count);
            count[0] = dim_solz;
            get_nc_var(&luts_scalar_par_inst_struct, "cos_solz", luts_data->scalar_inst_luts.cos_solz,
                       start, count);
            count[0] = dim_cot;
            get_nc_var(&luts_scalar_par_inst_struct, "cot", luts_data->scalar_inst_luts.cot,
                       start, count);
        }
        // reading 2d and 3d vars
        {
            size_t start[] = {0, 0, 0};
            size_t count[] = {dim_wind_speed, dim_solz, dim_cot};
            get_nc_var(&luts_scalar_par_inst_struct, "PARo_below_clear_sky", luts_data->scalar_inst_luts.pard_m_cs,
                       start, count);
            get_nc_var(&luts_scalar_par_inst_struct, "PARd_above_clear_sky", luts_data->scalar_inst_luts.pard_p_cs,
                       start, count);
            get_nc_var(&luts_scalar_par_inst_struct, "PARo_below_overcast", luts_data->scalar_inst_luts.pard_m_oc,
                       start, count);

            get_nc_var(&luts_scalar_par_inst_struct, "CF_PARo_below", luts_data->scalar_inst_luts.cf_pard_m,
                       start, count);
            get_nc_var(&luts_scalar_par_inst_struct, "CF_PARd_above", luts_data->scalar_inst_luts.cf_pard_p,
                       start, count);
        }
    }
}

void GetAerPhase(l2str *l2rec, int ip, int32_t nbands, float angstrom,
                 float *aerphasefunc, float *omega, float *modelAngstrom) {
    float phaseangle = l2rec->l1rec->scattang[ip];
    float distx, dista, p1, p2;
    int i, j, ix1 = 0, ix2, ia1 = 0, ia2, widx;

    static int firstCall = TRUE;

    // The data table is not in monotonically-increasing order, so here are the
    // indices in proper order:

    int ii[12] = {1, 0, 5, 4, 8, 7, 3, 2, 6, 11, 10, 9};

    // Read the aerosol data file:

    if (firstCall) {
        firstCall = FALSE;
        tablewavelengths =
                (float *) allocateMemory(nbands * sizeof(float), "tablewavelengths");
        tablephaseangles =
                (float *) allocateMemory(NPHASE * sizeof(float), "tablephaseangles");
        tablealphas = (float *) allocateMemory(nbands * NMODEL * sizeof(float),
                                               "tablealphas");
        tableomegas = (float *) allocateMemory(nbands * NMODEL * sizeof(float),
                                               "tableomegas");
        tableaerphasefunc = (float *) allocateMemory(
                NPHASE * nbands * NMODEL * sizeof(float), "tableaerphasefunc");
        read_aerosol_par(l2rec, nbands, tablewavelengths, tablephaseangles,
                         tablealphas, tableomegas, tableaerphasefunc);
    }

    // Find out which tablephaseangles-values in the table will be used:

    if (phaseangle < tablephaseangles[0]) {
        //        printf("Input phase angle too small, setting to 0. Was %f",
        //        phaseangle);
        phaseangle = tablephaseangles[0];
        ix1 = 0;
    } else if (phaseangle >= tablephaseangles[NPHASE - 1]) {
        //        printf("Input phase angle too large, setting to 180. Was %f",
        //                phaseangle);
        phaseangle = tablephaseangles[NPHASE - 1];
        ix1 = NPHASE - 2;
    } else {
        for (i = NPHASE - 2; i > 0; i--) {
            if (phaseangle >= tablephaseangles[i]) {
                ix1 = i;
                break;
            }
        }
    }
    ix2 = ix1 + 1;

    distx = (tablephaseangles[ix2] - phaseangle) /
            (tablephaseangles[ix2] - tablephaseangles[ix1]);

    // Find out which models in the table to use:

    widx = windex(670.0, tablewavelengths, nbands);

    if (angstrom < tablealphas[widx * NMODEL + ii[0]]) {
        angstrom = tablealphas[widx * NMODEL + ii[0]];
        ia1 = 0;
    } else if (angstrom >= tablealphas[widx * NMODEL + ii[NMODEL - 1]]) {
        angstrom = tablealphas[widx * NMODEL + ii[NMODEL - 1]];
        ia1 = NMODEL - 2;
    } else {
        for (j = NMODEL - 2; j > 0; j--) {
            if (angstrom >= tablealphas[widx * NMODEL + ii[j]]) {
                ia1 = j;
                break;
            }
        }
    }
    ia2 = ia1 + 1;
    dista = (tablealphas[widx * NMODEL + ii[ia2]] - angstrom) /
            (tablealphas[widx * NMODEL + ii[ia2]] -
             tablealphas[widx * NMODEL + ii[ia1]]);

    // Loop through the wavelengths and perform interpolations:

    for (i = 0; i < nbands; i++) {
        omega[i] = (tableomegas[i * NMODEL + ii[ia1]] * dista) +
                   (tableomegas[i * NMODEL + ii[ia2]] * (1.0 - dista));

        modelAngstrom[i] = (tablealphas[i * NMODEL + ii[ia1]] * dista) +
                           (tablealphas[i * NMODEL + ii[ia2]] * (1.0 - dista));

        p1 = (tableaerphasefunc[i * NMODEL * NPHASE + ii[ia1] * NPHASE + ix1] *
              distx) +
             (tableaerphasefunc[i * NMODEL * NPHASE + ii[ia1] * NPHASE + ix2] *
              (1.0 - distx));

        p2 = (tableaerphasefunc[i * NMODEL * NPHASE + ii[ia2] * NPHASE + ix1] *
              distx) +
             (tableaerphasefunc[i * NMODEL * NPHASE + ii[ia2] * NPHASE + ix2] *
              (1.0 - distx));

        aerphasefunc[i] = (p1 * dista) + (p2 * (1.0 - dista));
    }

    return;
}

/********************************************************************************
 *
 *	Description:
 *       This is to read the SeaWiFS aerosol data file for the PAR.
 *                               Menghua Wang 6/22/99
 *	Parameters:
 *       angl, R, ---- scattering angles.
 *       omega0, R, --- single scattering albedo.
 *       alpha, R, --- Angstrom coefficient alpha(lambda,865).
 *       s11, R, --- aerosol phase function.
 *
 *	SeaWiFS aerosol models:
 *       1-2: Oceanic RH=90 and 99%
 *       3-6: Maritime RH=50%, 70, 90, and 99%
 *       7-9: Coastal RH=50, 90, and 99%
 *       10-12: Tropospheric RH=50, 90, and 99%
 *
 *******************************************************************************/
void read_aerosol_par(l2str *l2rec, int32_t nbands, float *tablewavelengths,
                      float *tablephaseangles, float *tablealphas,
                      float *tableomegas, float *tableaerphasefunc) {
    int imodel[NMODEL];

    char aerosolfilename[FILENAME_MAX] = "";
    char *dataroot;
    FILE *aerosoldata;
    // Get OCDATAROOT path from user:
    if ((dataroot = getenv("OCDATAROOT")) == NULL) {
        printf("OCDATAROOT environment variable is not defined.\n");
        return;
    }
    sprintf(aerosolfilename, "%s/%s/%s_aerosol_par.dat", dataroot,
            sensorId2SensorDir(l2rec->l1rec->l1file->sensorID),
            sensorId2SensorDir(l2rec->l1rec->l1file->sensorID));

    // Opens the aerosol data file for reading.
    if ((aerosoldata = fopen(aerosolfilename, "r")) == NULL) {
        printf("Error reading aerosol table for PAR from %s.\n",
               aerosolfilename);
        exit(EXIT_FAILURE);
    }
    printf("Loading aerosol properties for PAR from %s.\n", aerosolfilename);
    int j, jj, iph;
    int num_model;
    for (jj = 0; jj < NPHASE; jj++) {
        if ((fscanf(aerosoldata, "%f", &tablephaseangles[jj])) == 0) {
            printf("Problem reading phase angles from %s.\n", aerosolfilename);
            exit(EXIT_FAILURE);
        }
    }
    for (j = 0; j < nbands; j++) {
        for (iph = 0; iph < NMODEL; iph++) {
            if ((fscanf(aerosoldata, "%d %d %f", &imodel[iph], &num_model,
                        &tablewavelengths[j])) == 0) {
                printf("Problem reading model number and wavelength from %s.\n",
                       aerosolfilename);
                exit(EXIT_FAILURE);
            }
            if ((fscanf(aerosoldata, "%f %f", &tableomegas[j * NMODEL + iph],
                        &tablealphas[j * NMODEL + iph])) == 0) {
                printf("Problem reading SSA and Angstrom from %s.\n",
                       aerosolfilename);
                exit(EXIT_FAILURE);
            }
            for (jj = 0; jj < NPHASE; jj++) {
                if ((fscanf(aerosoldata, "%f",
                            &tableaerphasefunc[j * NMODEL * NPHASE +
                                               iph * NPHASE + jj])) == 0) {
                    printf("Problem reading phase function values from %s.\n",
                           aerosolfilename);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    fclose(aerosoldata);
}

/********************************************************************************
 *
 *	  Estimate Dobson units from climatology, given the month and
 *      latitude.  Table has 12 columns, one per month, and 35 rows,
 *      from 85 N to -85, in steps of 5 deg lat.
 *      Perform 2d interpolation (latitude,month_day) to get estimate.
 *
 *     Author: John McPherson <JMcPherson@ucsd.edu>
 *
 *      Dobson look-up table from Keating et al., Ozone reference models
 *      for the middle atmosphere, Handbook for MAP, Vol. 31, 1-36, 1989
 *
 */

float EstimateDobson(int32_t year, int32_t month, int32_t day, float lat) {
    int i1, i2, m1, m2;
    float daymid, d1, d2;
    float t, t1, t2, fac, fac2, difflat;
    float dobson;

    int tabdobson[12][35] = {
            {395, 395, 395, 395, 395, 392, 390, 387, 376, 354, 322, 292,
                    269, 254, 248, 246, 247, 251, 255, 260, 266, 271, 277, 286,
                    295, 306, 319, 334, 344, 344, 338, 331, 324, 320, 316},

            {433, 433, 433, 436, 432, 428, 426, 418, 402, 374, 338, 303,
                    278, 261, 251, 246, 248, 250, 254, 258, 262, 265, 270, 278,
                    286, 294, 303, 313, 322, 325, 324, 317, 306, 299, 294},

            {467, 470, 460, 459, 451, 441, 433, 420, 401, 377, 347, 316,
                    291, 271, 260, 254, 254, 255, 257, 259, 261, 264, 269, 277,
                    284, 289, 296, 305, 312, 315, 317, 312, 305, 299, 295},

            {467, 465, 462, 455, 444, 431, 421, 410, 395, 373, 348, 325,
                    304, 287, 275, 267, 261, 259, 258, 259, 260, 263, 271, 278,
                    284, 289, 297, 306, 314, 318, 319, 313, 302, 302, 302},

            {411, 414, 416, 415, 410, 406, 402, 394, 382, 363, 342, 324,
                    307, 291, 279, 271, 264, 260, 258, 257, 258, 264, 271, 281,
                    291, 303, 312, 318, 322, 323, 322, 322, 322, 322, 322},

            {371, 371, 370, 368, 367, 372, 375, 372, 360, 341, 323, 311,
                    301, 290, 282, 275, 268, 263, 259, 256, 258, 264, 273, 289,
                    306, 319, 327, 328, 328, 337, 337, 337, 337, 337, 337},

            {333, 332, 332, 334, 338, 346, 350, 346, 335, 321, 310, 302,
                    296, 289, 284, 280, 274, 268, 262, 259, 261, 268, 279, 295,
                    315, 331, 340, 342, 338, 344, 340, 340, 340, 340, 340},

            {311, 308, 308, 313, 320, 327, 330, 326, 319, 310, 303, 298,
                    291, 286, 283, 281, 277, 273, 268, 264, 266, 274, 288, 306,
                    327, 343, 353, 355, 351, 339, 325, 307, 294, 294, 294},

            {283, 291, 302, 308, 312, 317, 318, 313, 307, 300, 295, 290,
                    284, 279, 279, 279, 278, 276, 272, 270, 273, 282, 295, 313,
                    333, 348, 360, 367, 368, 353, 324, 291, 267, 253, 230},

            {299, 299, 299, 309, 315, 317, 317, 312, 302, 291, 283, 280,
                    275, 270, 268, 267, 263, 263, 265, 269, 277, 287, 301, 317,
                    336, 354, 371, 387, 402, 402, 374, 333, 294, 274, 259},

            {314, 314, 314, 314, 332, 332, 327, 322, 311, 297, 284, 276,
                    270, 263, 261, 260, 258, 259, 264, 270, 278, 286, 298, 311,
                    323, 335, 350, 366, 381, 390, 388, 376, 357, 346, 341},

            {358, 358, 358, 358, 358, 358, 353, 349, 338, 320, 299, 281,
                    267, 256, 252, 251, 251, 253, 257, 264, 272, 279, 287, 297,
                    307, 318, 332, 347, 358, 365, 366, 364, 358, 356, 353}};

    int days[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    //   Set up for the time interpolation:

    if (year % 4 == 0) days[1] = days[1] + 1;

    daymid = days[month - 1] / 2.0;

    if (day >= daymid) {
        m1 = month - 1;
        m2 = m1 + 1;
        if (m2 > 12) m2 = 1;
        t = day;
    } else {
        m2 = month - 1;
        m1 = m2 - 1;
        if (m1 < 1) m1 = 12;
        t = days[m1] + day;
    }

    t1 = days[m1] / 2.0;
    t2 = days[m1] + (days[m2] / 2.0);

    // Perform the lat. interpolation, and further set up for time interp.:

    if (lat >= 85.0) {
        d1 = tabdobson[m1][1];
        d2 = tabdobson[m2][1];
    } else if (lat <= -85.0) {
        d1 = tabdobson[m1][34];
        d2 = tabdobson[m2][34];
    } else {
        i1 = 17 - (int) (lat / 5.0);
        i2 = i1 + 1;
        fac = (tabdobson[month][i2] - tabdobson[month][i1]) / (-5.0);
        difflat = lat - (90.0 - (i1 * 5.0));
        d1 = tabdobson[m1][i1] + (fac * difflat);
        d2 = tabdobson[m2][i1] + (fac * difflat);
    }

    // Complete the time interpolation to get final estimate:

    fac2 = (d2 - d1) / (t2 - t1);

    dobson = d1 + (t - t1) * fac2;
    return dobson;
}

/* ****************************************************************************
 *
 *       Interpolate water vapor table, as for dobson (ozone).
 *       Table has 12 columns, one per month, and 35 rows,
 *       from 85 N to 85 S, in steps of 5 deg lat.
 *       winter peak (Jan 15), summer peak (July 15), north hemis
 *       Perform 2d interpolation (latitude,month_day) to get estimate.
 *
 *      Tropical: if ( ABS( Lat ) .LE. 20.0 ) WatVap = 4.12
 *
 *      Mid-latitude case: if ( ABS( Lat ) .LE. 60.0 ) Then
 *        if ( (North.And.Summer) .Or. (South.And.Winter) ) WatVap = 2.93
 *        else WatVap = 0.85
 *
 *      Sub-arctic case:
 *         if ( (North.And.Summer) .Or. (South.And.Winter) ) WatVap = 2.10
 *         else WatVap = 0.42
 *
 *      after McClatchey, R.A., et al, _Optical Properties of the Atmosphere_,
 *      AFCRL 71-0279 354, 91 pp., 1971
 *
 */

float EstimateWatVap(int32_t year, int32_t month, int32_t day, float lat) {
    int i1, i2, m1, m2;
    month -= 1;
    float daymid, d1, d2;
    float t, t1, t2, fac, fac2, difflat;
    float watvap;

    float tabwv[12][35] = {
            {0.42, 0.47, 0.52, 0.56, 0.61, 0.66, 0.71, 0.75, 0.80, 0.85, 1.26, 1.67,
                    2.08, 2.48, 2.89, 3.30, 3.71, 4.12, 3.97, 3.82, 3.67, 3.53, 3.38, 3.23,
                    3.08, 2.93, 2.84, 2.75, 2.65, 2.56, 2.47, 2.38, 2.28, 2.19, 2.10},

            {0.70, 0.76, 0.81, 0.87, 0.92, 0.98, 1.03, 1.09, 1.14, 1.20, 1.56, 1.93,
                    2.29, 2.66, 3.02, 3.39, 3.75, 4.12, 3.93, 3.74, 3.54, 3.35, 3.16, 2.97,
                    2.78, 2.58, 2.50, 2.41, 2.33, 2.24, 2.16, 2.07, 1.99, 1.90, 1.82},

            {0.98, 1.04, 1.11, 1.17, 1.23, 1.29, 1.36, 1.42, 1.48, 1.54, 1.87, 2.19,
                    2.51, 2.83, 3.15, 3.48, 3.80, 4.12, 3.88, 3.65, 3.41, 3.18, 2.94, 2.71,
                    2.47, 2.24, 2.16, 2.08, 2.00, 1.93, 1.85, 1.77, 1.69, 1.62, 1.54},

            {1.26, 1.33, 1.40, 1.47, 1.54, 1.61, 1.68, 1.75, 1.82, 1.89, 2.17, 2.45,
                    2.73, 3.01, 3.28, 3.56, 3.84, 4.12, 3.84, 3.56, 3.28, 3.01, 2.73, 2.45,
                    2.17, 1.89, 1.82, 1.75, 1.68, 1.61, 1.54, 1.47, 1.40, 1.33, 1.26},

            {1.54, 1.62, 1.69, 1.77, 1.85, 1.93, 2.00, 2.08, 2.16, 2.24, 2.47, 2.71,
                    2.94, 3.18, 3.41, 3.65, 3.88, 4.12, 3.80, 3.48, 3.15, 2.83, 2.51, 2.19,
                    1.87, 1.54, 1.48, 1.42, 1.36, 1.29, 1.23, 1.17, 1.11, 1.04, 0.98},

            {1.82, 1.90, 1.99, 2.07, 2.16, 2.24, 2.33, 2.41, 2.50, 2.58, 2.78, 2.97,
                    3.16, 3.35, 3.54, 3.74, 3.93, 4.12, 3.75, 3.39, 3.02, 2.66, 2.29, 1.93,
                    1.56, 1.20, 1.14, 1.09, 1.03, 0.98, 0.92, 0.87, 0.81, 0.76, 0.70},

            {2.10, 2.19, 2.28, 2.38, 2.47, 2.56, 2.65, 2.75, 2.84, 2.93, 3.08, 3.23,
                    3.38, 3.53, 3.67, 3.82, 3.97, 4.12, 3.71, 3.30, 2.89, 2.49, 2.08, 1.67,
                    1.26, 0.85, 0.80, 0.75, 0.71, 0.66, 0.61, 0.56, 0.52, 0.47, 0.42},

            {1.82, 1.90, 1.99, 2.07, 2.16, 2.24, 2.33, 2.41, 2.50, 2.58, 2.78, 2.97,
                    3.16, 3.35, 3.54, 3.74, 3.93, 4.12, 3.75, 3.39, 3.02, 2.66, 2.29, 1.93,
                    1.56, 1.20, 1.14, 1.09, 1.03, 0.98, 0.92, 0.87, 0.81, 0.76, 0.70},

            {1.54, 1.62, 1.69, 1.77, 1.85, 1.93, 2.00, 2.08, 2.16, 2.24, 2.47, 2.71,
                    2.94, 3.18, 3.41, 3.65, 3.88, 4.12, 3.80, 3.48, 3.15, 2.83, 2.51, 2.19,
                    1.87, 1.54, 1.48, 1.42, 1.36, 1.29, 1.23, 1.17, 1.11, 1.04, 0.98},

            {1.26, 1.33, 1.40, 1.47, 1.54, 1.61, 1.68, 1.75, 1.82, 1.89, 2.17, 2.45,
                    2.73, 3.01, 3.28, 3.56, 3.84, 4.12, 3.84, 3.56, 3.28, 3.01, 2.73, 2.45,
                    2.17, 1.89, 1.82, 1.75, 1.68, 1.61, 1.54, 1.47, 1.40, 1.33, 1.26},

            {0.98, 1.04, 1.11, 1.17, 1.23, 1.29, 1.36, 1.42, 1.48, 1.54, 1.87, 2.19,
                    2.51, 2.83, 3.15, 3.48, 3.80, 4.12, 3.88, 3.65, 3.41, 3.18, 2.94, 2.71,
                    2.47, 2.24, 2.16, 2.08, 2.00, 1.93, 1.85, 1.77, 1.69, 1.62, 1.54},

            {0.70, 0.76, 0.81, 0.87, 0.92, 0.98, 1.03, 1.09, 1.14, 1.20, 1.56, 1.93,
                    2.29, 2.66, 3.02, 3.39, 3.75, 4.12, 3.93, 3.74, 3.54, 3.35, 3.16, 2.97,
                    2.78, 2.58, 2.50, 2.41, 2.33, 2.24, 2.16, 2.07, 1.99, 1.90, 1.82}};

    int days[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    // Set up for the time interpolation:

    if (year % 4 == 0) days[1] = days[1] + 1;

    daymid = days[month] / 2.0;

    if (day >= daymid) {
        m1 = month;
        m2 = m1 + 1;
        if (m2 > 11) m2 = 0;
        t = day;
    } else {
        m2 = month;
        m1 = m2 - 1;
        if (m1 < 0) m1 = 11;
        t = days[m1] + day;
    }

    t1 = days[m1] / 2.0;
    t2 = days[m1] + (days[m2] / 2.0);

    // Perform the lat. interpolation, and further set up for time interp.:

    if (lat >= 85.0) {
        d1 = tabwv[m1][0];
        d2 = tabwv[m2][0];
    } else if (lat <= (-85.0)) {
        d1 = tabwv[m1][34];
        d2 = tabwv[m2][34];
    } else {
        i1 = 16 - (int) (lat / 5.0);
        i2 = i1 + 1;
        fac = (tabwv[month][i2] - tabwv[month][i1]) / (-5.0);
        difflat = lat - (90.0 - (i1 * 5.0));
        d1 = tabwv[m1][i1] + (fac * difflat);
        d2 = tabwv[m2][i1] + (fac * difflat);
    }

    // Complete the time interpolation to get final estimate:

    fac2 = (d2 - d1) / (t2 - t1);

    watvap = d1 + (t - t1) * fac2;

    return watvap;
}

/********************************************************************************
 *       		Get Earth-Sun distance correction	(VarSol)
 *
 *            Calculation of the variability of the solar constant during the
 *year. jday is the number of the day in the month dsol is a multiplicative
 *factor to apply to the mean value of solar constant.
 *
 ********************************************************************************/

float varsol(int32_t jday) {
    float om;
    float dsol;

    om = (0.9856f * (float) (jday - 4)) * OEL_PI / 180.f;
    dsol = 1.f / powf((1.f - .01673f * cosf(om)), 2.f);

    return dsol;
}

/*************************************************************************************
 *        Compute the rise and set time for the given Julian day of year
 **
 *************************************************************************************/

void triseset(int32_t jday, float xlon, float xlat, float *trise, float *tset) {
    float fac, xlo, xla, xj, a1, a2, et, a3, delta, cosah, ah1, ah2, tsv1, tsv2,
            tsm1, tsm2;

    fac = (float) OEL_PI / 180.f;

    xlo = xlon * fac;
    xla = xlat * fac;
    xj = (float) (jday - 1);

    a1 = (1.00554f * xj - 6.28306f) * fac;
    a2 = (1.93946f * xj + 23.35089f) * fac;
    et = -7.67825f * sinf(a1) - 10.09176f * sinf(a2);
    a3 = (0.9683f * xj - 78.00878f) * fac;

    delta = 23.4856f * sinf(a3) * fac;
    cosah = (-(sinf(xla) * sinf(delta))) / (cosf(xla) * cosf(delta));

    if ((cosah < -1.f) || (cosah > 1.f)) {
        *trise = 0.0f;
        *tset = 24.0f;
        return;
    }

    ah1 = acosf(cosah);
    ah2 = -ah1;

    tsv1 = ah1 / (15.f * fac);
    tsv2 = ah2 / (15.f * fac);
    tsm1 = tsv1 + 12.f - et / 60.f;
    tsm2 = tsv2 + 12.f - et / 60.f;

    *trise = tsm2 - (xlo / (15.f * fac));
    *tset = tsm1 - (xlo / (15.f * fac));
    return;
}

/********************************************************************************
 *     Get Sun position
 ********************************************************************************/

float get_solz(int32_t jday, float time, float lon, float lat) {
    float asol;
    double tsm, xla, xj, tet, a1, a2, a3, a4, a5, et, tsv, ah, b1, b2, b3, b4,
            b5, b6, b7, delta, amuzero, elev, az, caz, azim, pi2;

    //     solar position (zenithal angle asol,azimuthal angle phi0 in degrees)
    //    j is the day number in the year
    //
    //   mean solar time (heure decimale) (aka decimal time)
    double fac = OEL_PI / 180.;
    tsm = time + lon / 15.;
    xla = lat * fac;
    xj = (float) (jday - 1);
    tet = 2. * OEL_PI * xj / 365.;

    //    time equation (in mn.dec)
    a1 = 0.000075;
    a2 = 0.001868;
    a3 = 0.032077;
    a4 = 0.014615;
    a5 = 0.040849;
    et = a1 + a2 * cos(tet) - a3 * sin(tet) - a4 * cos(2. * tet) -
         a5 * sin(2. * tet);
    et = et * 12. * 60. / OEL_PI;

    //     true solar time

    tsv = tsm + et / 60.;
    tsv -= 12.;

    //     hour angle

    ah = tsv * 15. * fac;

    //     solar declination (in radian)

    b1 = 0.006918;
    b2 = 0.399912;
    b3 = 0.070257;
    b4 = 0.006758;
    b5 = 0.000907;
    b6 = 0.002697;
    b7 = 0 - .001480;
    delta = b1 - b2 * cos(tet) + b3 * sin(tet) - b4 * cos(2. * tet) +
            b5 * sin(2. * tet) - b6 * cos(3. * tet) + b7 * sin(3. * tet);

    //     elevation,azimuth

    amuzero = sin(xla) * sin(delta) + cos(xla) * cos(delta) * cos(ah);
    if ((fabs(amuzero) - 1.000) > 0.00000) amuzero = copysign(amuzero, 1.);
    elev = asin(amuzero);
    az = cos(delta) * sin(ah) / cos(elev);
    if ((fabs(az) - 1.000) > 0.00000) az = copysign(az, 1.);

    caz =
            (-cos(xla) * sin(delta) + sin(xla) * cos(delta) * cos(ah)) / cos(elev);
    azim = asin(az);

    if (caz <= 0.) azim = OEL_PI - azim;

    if (caz > 0. && az <= 0.) azim = 2. * OEL_PI + azim;

    azim = azim + OEL_PI;
    pi2 = 2.f * OEL_PI;

    if (azim > pi2) azim -= pi2;

    elev /= fac;

    //    conversion in degrees

    asol = (float) (90. - elev);
    //	phi0 = azim / fac;

    return asol;
}

/******************************************************************************
 *      Given inputs cos(SZ) and tauA500, interpolate the appropriate * surface
 *albedo from the table. This is for the case where tau < 1.0	  * Inputs are
 *assumed to be in the valid ranges (the associated solar 	  * zenith angle
 *should be in range 0.0 - 87.134 degrees, and tau		  * should be in
 *the range 0.0 - 0.99)
 **
 ******************************************************************************/

float interp_as_taulow(float csz, float tau) {
    int s1, s2, t1, t2, i;

    float stot, sdist, ttot, tdist, slope, alb1[2], as;

    float cszt[NS] = {0.05, 0.09, 0.15, 0.21, 0.27, 0.33, 0.39, 0.45,
                      0.52, 0.60, 0.68, 0.76, 0.84, 0.92, 1.00};

    float taut[NT] = {0.00, 0.05, 0.10, 0.16, 0.24, 0.35, 0.50, 0.70, 0.99};

    //       Store the table here rather than read a new file. Fortran stores
    //      2-d arrays internally by columns (col1, col2, col3 ...)

    float ast[NS][NT] = {{0.1592253, 0.1258933, 0.1088253, 0.0968594, 0.0880966,
                                 0.0815492, 0.0766266, 0.0721790, 0.0679746},
                         {0.1944972, 0.1637040, 0.1432702, 0.1240363, 0.1066057,
                                 0.0922638, 0.0826043, 0.0757601, 0.0700167},
                         {0.1978650, 0.1785707, 0.1643212, 0.1482773, 0.1304572,
                                 0.1119496, 0.0963537, 0.0839010, 0.0740948},
                         {0.1785221, 0.1673987, 0.1589092, 0.1486409, 0.1359193,
                                 0.1209508, 0.1061981, 0.0920763, 0.0793196},
                         {0.1531512, 0.1473233, 0.1426633, 0.1366985, 0.1289257,
                                 0.1188203, 0.1078114, 0.0957847, 0.0831087},
                         {0.1281988, 0.1255739, 0.1234628, 0.1204020, 0.1161201,
                                 0.1101539, 0.1030820, 0.0943791, 0.0841410},
                         {0.1061354, 0.1052812, 0.1048179, 0.1036617, 0.1016918,
                                 0.0986726, 0.0948040, 0.0894514, 0.0822926},
                         {0.0877530, 0.0881279, 0.0883850, 0.0884469, 0.0880587,
                                 0.0870923, 0.0854952, 0.0826607, 0.0782380},
                         {0.0708031, 0.0720173, 0.0727886, 0.0736250, 0.0741367,
                                 0.0746769, 0.0747173, 0.0741294, 0.0722523},
                         {0.0567974, 0.0582123, 0.0593260, 0.0604172, 0.0615151,
                                 0.0627740, 0.0639047, 0.0647193, 0.0650727},
                         {0.0472026, 0.0485713, 0.0495830, 0.0507434, 0.0519943,
                                 0.0536504, 0.0551176, 0.0566950, 0.0581553},
                         {0.0406631, 0.0419177, 0.0429259, 0.0440614, 0.0451823,
                                 0.0467889, 0.0483670, 0.0501810, 0.0522433},
                         {0.0366926, 0.0377500, 0.0384530, 0.0393431, 0.0404503,
                                 0.0419569, 0.0434430, 0.0451987, 0.0473454},
                         {0.0343793, 0.0352937, 0.0358116, 0.0364891, 0.0374246,
                                 0.0385732, 0.0398924, 0.0414707, 0.0435983},
                         {0.0331561, 0.0337733, 0.0341567, 0.0346916, 0.0354239,
                                 0.0364011, 0.0374280, 0.0387133, 0.0405543}};

    //     Set up interpolation:

    //      Locate subcells in array to use:

    s1 = NS - 2;
    s2 = NS - 1;
    for (i = 0; i < NS - 1; i++) {
        if (csz < cszt[i + 1]) {
            s1 = i;
            s2 = i + 1;
            break;
        }
    }
    t1 = NT - 2;
    t2 = NT - 1;
    for (i = 0; i < NT - 1; i++) {
        if (tau < taut[i + 1]) {
            t1 = i;
            t2 = i + 1;
            break;
        }
    }

    stot = cszt[s2] - cszt[s1];
    sdist = csz - cszt[s1];
    ttot = taut[t2] - taut[t1];
    tdist = tau - taut[t1];

    slope = (ast[s2][t1] - ast[s1][t1]) / stot;
    alb1[0] = ast[s1][t1] + (slope * sdist);
    slope = (ast[s2][t2] - ast[s1][t2]) / stot;
    alb1[1] = ast[s1][t2] + (slope * sdist);

    slope = (alb1[1] - alb1[0]) / ttot;
    as = alb1[0] + (slope * tdist);

    return as;
}

/********************************************************************************
 *      Given input cos(SZ), interpolate the appropriate surface albedo
 ** from the table. This is for the case where tau > 1.0
 ** Input is assumed to be in the valid range (the associated solar
 ** zenith angle should be in range 0.0 - 87.134 degrees)
 **
 ********************************************************************************/

float interp_as_tauhigh(float csz) {
    int s1, s2, i;

    float stot, sdist, slope, as;

    float cszt[NS] = {0.05, 0.09, 0.15, 0.21, 0.27, 0.33, 0.39, 0.45,
                      0.52, 0.60, 0.68, 0.76, 0.84, 0.92, 1.00};

    float ast[NS] = {0.0623166, 0.0630070, 0.0640739, 0.0650397, 0.0657760,
                     0.0661690, 0.0659415, 0.0651271, 0.0635101, 0.0610652,
                     0.0583470, 0.0556146, 0.0530646, 0.0509498, 0.0490149};

    //     Set up interpolation:

    //      Locate subcells in array to use:

    s1 = NS - 2;
    s2 = NS - 1;
    for (i = 0; i < NS; i++) {
        if (csz < cszt[i + 1]) {
            s1 = i;
            s2 = i + 1;
            break;
        }
    }

    stot = cszt[s2] - cszt[s1];
    sdist = csz - cszt[s1];

    slope = (ast[s2] - ast[s1]) / stot;
    as = ast[s1] + (slope * sdist);

    return as;
}

void getcldalbe(float *TauCld, float *CF, float cosSZ, float t_obs,
                float *t_range, float *albe_obs, float *TauCld_obs,
                float *CF_obs, size_t t_step, float *wl, size_t bands_num) {
    float wave[] = {370., 440., 560., 650., 720.};
    float a1[] = {0.78, 0.65, 0.53, 0.49, 0.51};
    float b1[] = {0.31, 0.59, 0.79, 0.88, 0.82};
    float b2[] = {-0.0782, -0.1212, -0.1054, -0.1276, -0.1054};
    float b3[] = {-0.9218, -0.8788, -0.8946, -0.8724, -0.8946};
    float k1[] = {1.3, 1.3, 1.3, 1.3, 1.3};
    float k2[] = {0.2000, 0.0229, 0.0409, 0.0436, 0.0575};
    float k3[] = {0.8000, 1.1754, 1.2436, 1.1836, 1.3637};
    float cc[] = {0.0835, 0.1008, 0.1153, 0.1215, 0.1307};
    float dd[] = {0.0759, 0.0954, 0.1126, 0.1192, 0.1218};

    size_t s1, s2;

    s1 = t_step - 2;
    s2 = t_step - 1;

    float Stot, Sdist, slope;

    search(t_range, 0, t_step - 1, t_obs, &s1, &s2);

    Stot = t_range[s2] - t_range[s1];
    Sdist = t_obs - t_range[s1];

    if (t_obs <= t_range[0]) Sdist = 0;
    if (t_obs >= t_range[t_step - 1]) Sdist = Stot;
    if (s1 != s2)
        slope = (TauCld[s2] - TauCld[s1]) / Stot;
    else
        slope = 0.0;
    *TauCld_obs = TauCld[s1] + (slope * Sdist);
    if (s1 != s2)
        slope = (CF[s2] - CF[s1]) / Stot;
    else
        slope = 0.0;
    *CF_obs = CF[s1] + (slope * Sdist);
    // if(isnan(*CF_obs))
    // {
    //     printf("CF_OBS IS NAN,%f, %f, %f ",CF[s1],slope,Sdist);
    //     exit(EXIT_FAILURE);
    // }
    float aa, bb;

    float trc[5];

    float alpha = 0;  //  # alpha - surface albedo

    for (size_t i = 0; i < 5; i++) {
        aa = a1[i] + (1 - a1[i]) * exp(-*TauCld_obs * k1[i]);
        bb = b1[i] * (1 + b2[i] * exp(-*TauCld_obs * k2[i]) +
                      b3[i] * exp(-*TauCld_obs * k3[i]));

        trc[i] =
                (aa + cosSZ * bb) / (1 + *TauCld_obs * (cc[i] - dd[i] * alpha));
    }
    for (size_t band = 0; band < bands_num; band++) {
        search(wave, 0, 4, wl[band], &s1, &s2);
        Stot = wave[s2] - wave[s1];
        float Sdist = wl[band] - wave[s1];
        float slope = (trc[s2] - trc[s1]) / Stot;
        if (s1 == s2) slope = 0.0;
        float trc_obs = trc[s1] + (slope * Sdist);
        albe_obs[band] = (1 - trc_obs) * (*CF_obs);
    }
}

// ##############################################################################
// # new parameterization of As, from Jin et al. (2011)
// # Jin, Z., Y. Qiao, Y. Wang, Y. Fang, and W. Yi, 2011: A new parameterization
// of # spectral and broadband ocean surface albedo. Opt. Express, 19, no. 27,
// 26429-26443, # doi:10.1364/OE.19.026429.
// ##############################################################################

// # Regression function (Equation 4)

float ff(float u, float v) {
    float p[] = {0.0152, -1.7873, 6.8972, -8.5778, 4.071, -7.6446,
                 0.1643, -7.8409, -3.5639, -2.3588, 10.0538};
    const float zmod =
            (p[0] + p[1] * u + p[2] * pow(u, 2) + p[3] * pow(u, 3) + p[4] * v +
             p[5] * u * v) *
            exp(p[6] + p[7] * u + p[8] * pow(u, 2) + p[9] * v + p[10] * u * v);
    return zmod;
}

/// @brief
/// @param wl
/// @param luts_data
/// @return
float getrefm(float wl, const luts_par *luts_data) {
    size_t dims[] = {luts_data->soa_lut.dim_wlt_reft};
    float *grid[1];
    grid[0] = luts_data->soa_lut.wlt_reft;
    float point[] = {wl};
    float refm = interp1d(dims, point, grid, luts_data->soa_lut.reft);
    if (refm > 1.34) refm = 1.34;
    return refm;
}

float getr0(float wl, float chl, float u0, const luts_par *luts_data) {
    float *grid[1];
    float point[] = {wl};
    float chlabs, bw, aw;
    // get chlor
    {
        size_t dims[] = {luts_data->soa_lut.dim_jwvl};
        grid[0] = luts_data->soa_lut.jwl;
        chlabs = interp1d(dims, point, grid, luts_data->soa_lut.achl);
    }

    // get bw
    {
        size_t dims[] = {luts_data->soa_lut.dim_wv1};
        grid[0] = luts_data->soa_lut.wv1;
        bw = interp1d(dims, point, grid, luts_data->soa_lut.bw1);
    }

    // get aw
    {
        size_t dims[] = {luts_data->soa_lut.dim_wv3};
        grid[0] = luts_data->soa_lut.wv3;
        aw = interp1d(dims, point, grid, luts_data->soa_lut.aw3);
    }
    const float ylmd = exp(0.014 * (440.0 - wl));
    const float aw440 = 0.00635;
    const float ap = 0.06 * chlabs * pow(chl, 0.65) +
                     0.2 * (aw440 + 0.06 * pow(chl, 0.65)) * ylmd;

    const float bp550 = 0.416 * pow(chl, 0.766);

    float v, bbp;

    if (chl >= 2) {
        v = 0.;
        bbp = (0.002 + 0.01 * (0.5 - 0.25 * log10(chl))) * bp550;
    } else {
        if (chl > 0.02) {
            v = 0.5 * (log10(chl) - 0.3);
            bbp =
                    (0.002 + 0.01 * (0.5 - 0.25 * log10(chl)) * pow(wl / 550., v)) *
                    bp550;
        } else {
            bbp = 0.019 * (550. / wl) * bp550;
        }
    }

    // # hb=h/(h+2*bbpf*(1.-h))	;Morel-Gentili(1991), Eq (12)
    float hb = 0.5 * bw / (0.5 * bw + bbp);

    // #;		use Morel 91 to get f

    float f =
            0.6279 - 0.2227 * hb - 0.0513 * hb * hb + (-0.3119 + 0.2465 * hb) * u0;

    float r0 = f * (0.5 * bw + bbp) / (aw + ap);
    return r0;
}

float getosa(float wl, float sza, float wind, float chl, float fr,
             const luts_par *luts_data) {  // # wl	wavelength (nm), set wl=0 for broadband
    // # chl	chlorophyll (g/m3)
    // # sza	solar zenith angle (degree)
    // # wind	wind speed (m/s)
    // # fr	fraction of diffuse incidence
    const float sig = sqrt(0.003 + 0.00512 * wind);  //       #convert to Sigma

    const float refm = getrefm(wl, luts_data);

    const float csz = cos(sza * OEL_PI / 180.);

    const float xx2 = sqrt(1.0 - (1.0 - csz * csz) / refm / refm);

    float rr0 = (0.5 * (pow(((xx2 - csz * refm) / (xx2 + csz * refm)), 2) +
                        pow(((csz - refm * xx2) / (csz + refm * xx2)), 2)));

    float rrr = 0.5 * (pow(((xx2 - 1.34 * csz) / (xx2 + 1.34 * csz)), 2) +
                       pow(((csz - 1.34 * xx2) / (csz + 1.34 * xx2)), 2));

    const float r11 = rr0 - ff(csz, sig) * rr0 / rrr;  // # direct surface
    // albedo

    // # Water volume scattering
    const float r22 = 0.48168549 - 0.014894708 * sig - 0.20703885 * sig * sig;

    float r00 = getr0(wl, chl, csz, luts_data);

    const float rw = r00 * (1. - r22) * (1. - r11) /
                     (1. - r00 * r22);  // #direct water-leaving albedo

    const float ue =
            0.676;  //               # the equivalent u_unif for diffuse incidence
    const float ue2 = sqrt(1.0 - (1.0 - ue * ue) / refm / refm);

    rr0 = 0.5 * (pow(((ue2 - refm * ue) / (ue2 + refm * ue)), 2) +
                 pow(((ue - refm * ue2) / (ue + refm * ue2)), 2));
    rrr = 0.5 * (pow(((ue2 - 1.34 * ue) / (ue2 + 1.34 * ue)), 2) +
                 pow(((ue - 1.34 * ue2) / (ue + 1.34 * ue2)), 2));

    float r11df = rr0 - ff(ue, sig) * rr0 / rrr;

    r00 = getr0(wl, chl, ue, luts_data);

    float rwdf = r00 * (1. - r22) * (1. - r11df) /
                 (1. - r00 * r22);  // # diffuse water-leaving albedo

    float pc[] = {-0.1482, -0.012, 0.1609, -0.0244};

    float rdf;
    if (fr > 0.99)
        rdf = -0.1479 + 0.1502 * refm - 0.0176 * sig * refm;
    else
        rdf = pc[0] + pc[1] * sig + pc[2] * refm +
              pc[3] * sig * refm;  // # surface diffuse (Eq 5a-5b)

    // float rwt =
    //     (1. - fr) * rw +
    //     fr * rwdf;  //          #total parameterized water-leaving albedo

    float albt = (1. - fr) * (r11 + rw) +
                 fr * (rdf + rwdf);  //          # total parameterized albedo

    // # Foam correction (Eq 16-17)
    float fwc = (2.95e-6) * pow(wind, 3.52);

    float osa = (1. - fwc) * albt + fwc * 0.55;
    return osa;
}

float SunGlint(float sz, float vz, float ra, float ws) {
    double radeg = OEL_RADEG;
    const float ssz = sin(sz * radeg);
    const float svz = sin(vz * radeg);
    const float csz = cos(sz * radeg);
    const float cvz = cos(vz * radeg);
    const float cra = cos(ra * radeg);

    // # Determine omega
    const float cos_2omega = csz * cvz + ssz * svz * cra;
    const float omega = acos(cos_2omega) / 2.0;  // # in radians
    const float sOmega = sin(omega);
    const float cOmega = cos(omega);
    const float omegaP = asin(sOmega / 1.34);  //  # in radians, Snell's law
    const float omegaDp = omega + omegaP;
    const float omegaDm = omega - omegaP;

    const float term1 = sin(omegaDm) / sin(omegaDp);
    const float term2 = tan(omegaDm) / tan(omegaDp);
    float r_omega = 0.5 * (term1 * term1 + term2 * term2);
    if (omega < OEL_DEGRAD) r_omega = 0.0211119;

    const float beta = acos((csz + cvz) / (2.0 * cOmega));  //   #  in radians
    const float cBeta = cos(beta);
    const float tBeta = tan(beta);
    const float sig2 = 0.003 + 0.00512 * ws;
    const float P_V_beta =
            (1.0 / (2.0 * sig2 * OEL_PI)) * exp(-tBeta * tBeta / (2.0 * sig2));
    const float rho_g =
            OEL_PI * P_V_beta * r_omega / (4.0 * csz * cvz * pow((cBeta), 4));
    return rho_g;
}

void set_searh_grid(size_t pts[][2], float width[], size_t ndims,
                    const size_t *dims, const float *point, float **grid) {
    for (size_t i_dim = 0; i_dim < ndims; i_dim++) {
        size_t dim_size = dims[i_dim];
        size_t st = 0;
        size_t end = dim_size - 1;
        search(grid[i_dim], st, end, point[i_dim], &st, &end);
        pts[i_dim][0] = st;
        pts[i_dim][1] = end;
    }
    for (size_t i_dim = 0; i_dim < ndims; i_dim++) {
        size_t st = pts[i_dim][0];
        size_t end = pts[i_dim][1];
        if (st != end)
            width[i_dim] = (point[i_dim] - grid[i_dim][st]) /
                           (grid[i_dim][end] - grid[i_dim][st]);
        else
            width[i_dim] = 0.0;
        // { std::cout << "width[i_dim]" << width[i_dim] << std::endl; }
    }
}

size_t get6dindex(const size_t *point, const size_t *dims) {
    return point[0] * dims[1] * dims[2] * dims[3] * dims[4] * dims[5] +
           point[1] * dims[2] * dims[3] * dims[4] * dims[5] +
           point[2] * dims[3] * dims[4] * dims[5] +
           point[3] * dims[4] * dims[5] + point[4] * dims[5] + point[5];
}

float interp6d(const size_t *dims, const float *point, float **grid,
               const float *lut) {
    const size_t ndims = 6;
    size_t pts[ndims][2];
    float width[ndims];
    set_searh_grid(pts, width, ndims, dims, point, grid);
    float temp0[2][2][2][2][2][2];
    for (size_t i0 = 0; i0 < 2; i0++)
        for (size_t i1 = 0; i1 < 2; i1++)
            for (size_t i2 = 0; i2 < 2; i2++)
                for (size_t i3 = 0; i3 < 2; i3++)
                    for (size_t i4 = 0; i4 < 2; i4++)
                        for (size_t i5 = 0; i5 < 2; i5++) {
                            size_t point_pos[] = {pts[0][i0], pts[1][i1],
                                                  pts[2][i2], pts[3][i3],
                                                  pts[4][i4], pts[5][i5]};
                            temp0[i0][i1][i2][i3][i4][i5] =
                                    lut[get6dindex(point_pos, dims)];
                        }
    float temp1[2][2][2][2][2];
    for (size_t i1 = 0; i1 < 2; i1++)
        for (size_t i2 = 0; i2 < 2; i2++)
            for (size_t i3 = 0; i3 < 2; i3++)
                for (size_t i4 = 0; i4 < 2; i4++)
                    for (size_t i5 = 0; i5 < 2; i5++) {
                        temp1[i1][i2][i3][i4][i5] =
                                (temp0[1][i1][i2][i3][i4][i5] -
                                 temp0[0][i1][i2][i3][i4][i5]) *
                                width[0] +
                                temp0[0][i1][i2][i3][i4][i5];
                    }

    float temp2[2][2][2][2];
    for (size_t i2 = 0; i2 < 2; i2++)
        for (size_t i3 = 0; i3 < 2; i3++)
            for (size_t i4 = 0; i4 < 2; i4++)
                for (size_t i5 = 0; i5 < 2; i5++) {
                    temp2[i2][i3][i4][i5] =
                            (temp1[1][i2][i3][i4][i5] - temp1[0][i2][i3][i4][i5]) *
                            width[1] +
                            temp1[0][i2][i3][i4][i5];
                }

    float temp3[2][2][2];
    for (size_t i3 = 0; i3 < 2; i3++)
        for (size_t i4 = 0; i4 < 2; i4++)
            for (size_t i5 = 0; i5 < 2; i5++) {
                temp3[i3][i4][i5] =
                        (temp2[1][i3][i4][i5] - temp2[0][i3][i4][i5]) * width[2] +
                        temp2[0][i3][i4][i5];
            }
    float temp4[2][2];
    for (size_t i4 = 0; i4 < 2; i4++)
        for (size_t i5 = 0; i5 < 2; i5++) {
            temp4[i4][i5] = (temp3[1][i4][i5] - temp3[0][i4][i5]) * width[3] +
                            temp3[0][i4][i5];
        }

    float temp5[2];
    for (size_t i5 = 0; i5 < 2; i5++) {
        temp5[i5] = (temp4[1][i5] - temp4[0][i5]) * width[4] + temp4[0][i5];
    }

    float temp6;
    temp6 = (temp5[1] - temp5[0]) * width[5] + temp5[0];

    return temp6;
}

size_t get4dindex(const size_t *point, const size_t *dims) {
    return point[0] * dims[1] * dims[2] * dims[3] +
           point[1] * dims[2] * dims[3] + point[2] * dims[3] + point[3];
}

float interp4d(const size_t *dims, const float *point, float **grid,
               const float *lut) {
    const size_t ndims = 4;
    size_t pts[ndims][2];
    float width[ndims];
    set_searh_grid(pts, width, ndims, dims, point, grid);
    float temp0[2][2][2][2];
    for (size_t i0 = 0; i0 < 2; i0++)
        for (size_t i1 = 0; i1 < 2; i1++)
            for (size_t i2 = 0; i2 < 2; i2++)
                for (size_t i3 = 0; i3 < 2; i3++) {
                    size_t point_pos[] = {pts[0][i0], pts[1][i1], pts[2][i2],
                                          pts[3][i3]};
                    temp0[i0][i1][i2][i3] = lut[get4dindex(point_pos, dims)];
                }
    float temp1[2][2][2];
    for (size_t i1 = 0; i1 < 2; i1++)
        for (size_t i2 = 0; i2 < 2; i2++)
            for (size_t i3 = 0; i3 < 2; i3++) {
                temp1[i1][i2][i3] =
                        (temp0[1][i1][i2][i3] - temp0[0][i1][i2][i3]) * width[0] +
                        temp0[0][i1][i2][i3];
            }

    float temp2[2][2];
    for (size_t i2 = 0; i2 < 2; i2++)
        for (size_t i3 = 0; i3 < 2; i3++) {
            temp2[i2][i3] = (temp1[1][i2][i3] - temp1[0][i2][i3]) * width[1] +
                            temp1[0][i2][i3];
        }

    float temp3[2];
    for (size_t i3 = 0; i3 < 2; i3++) {
        temp3[i3] = (temp2[1][i3] - temp2[0][i3]) * width[2] + temp2[0][i3];
    }

    float temp4;
    temp4 = (temp3[1] - temp3[0]) * width[3] + temp3[0];

    return temp4;
}

size_t get3dindex(const size_t *point, const size_t *dims) {
    return point[0] * dims[1] * dims[2] + point[1] * dims[2] + point[2];
}

float interp3d(const size_t *dims, const float *point, float **grid,
               const float *lut) {
    const size_t ndims = 3;
    size_t pts[ndims][2];
    float width[ndims];
    set_searh_grid(pts, width, ndims, dims, point, grid);
    float temp0[2][2][2];
    for (size_t i0 = 0; i0 < 2; i0++)
        for (size_t i1 = 0; i1 < 2; i1++)
            for (size_t i2 = 0; i2 < 2; i2++) {
                size_t point_pos[] = {pts[0][i0], pts[1][i1], pts[2][i2]};
                temp0[i0][i1][i2] = lut[get3dindex(point_pos, dims)];
            }
    float temp1[2][2];
    for (size_t i1 = 0; i1 < 2; i1++)
        for (size_t i2 = 0; i2 < 2; i2++) {
            temp1[i1][i2] = (temp0[1][i1][i2] - temp0[0][i1][i2]) * width[0] +
                            temp0[0][i1][i2];
        }

    float temp2[2];
    for (size_t i2 = 0; i2 < 2; i2++) {
        temp2[i2] = (temp1[1][i2] - temp1[0][i2]) * width[1] + temp1[0][i2];
    }

    float temp3 = (temp2[1] - temp2[0]) * width[2] + temp2[0];

    return temp3;
}

size_t get1dindex(const size_t *point, const size_t *dims) { return point[0]; }

float interp1d(const size_t *dims, const float *point, float **grid,
               const float *lut) {
    const size_t ndims = 1;
    size_t pts[ndims][2];
    float width[ndims];
    set_searh_grid(pts, width, ndims, dims, point, grid);
    float temp0[2];
    for (size_t i0 = 0; i0 < 2; i0++) {
        size_t point_pos[] = {pts[0][i0]};
        temp0[i0] = lut[get1dindex(point_pos, dims)];
    }

    float temp1;
    temp1 = (temp0[1] - temp0[0]) * width[0] + temp0[0];

    return temp1;
}
