#include "argpar.h"
#include "olog.h"
#include "olog/loader.h"
#include "shash.h"
#include "val_extract.h"
#include "vincenty.h"

#include <errno.h>
#include <float.h>
#include <string.h>

#define OFILE_DEFAULT_EXT ".qc"

#ifdef DEBUG
#define dprintf(...) do { printf(__VA_ARGS__); } while(0)
#else
#define dprintf(...) do {} while (0)
#endif

#define eprintf(...) do { fprintf(stderr, __VA_ARGS__); } while(0)

#define str(s) #s
#define expanded_str(s) str(s)

typedef struct sstval_extract_main_input {
    char *ofile;
    unsigned product_count;
    shash *products_printed;
    FILE *base_ofile_h;
    char *qual_check;
    double qual_check_distance;
    nc_var **qual_check_var_copies;
    double **qual_check_data_copies;
    int box_size;
    val_extract_arguments *val_extract_arguments;
    double slat_orig, slon_orig;
} sstval_extract_main_input;

static argpar_option options[] = {
        { "ofile", 'o', "FILE", 0, "base output file path" },
        { "qual_check", 'q', "VAR", 0, "quality variable (lower = better)" },
        { "qual_check_distance", 'd', "Km", OPTION_DBL, "max distance to find better quality" },
        { 0, 0, 0, 0, "Output files:", -2 },
        { 0, 0, 0, OPTION_DOC,
                "Any number of products may be specified.  If none are given, every product that appears to be "
                        "geospatial, having scan line and pixel as the only dimensions, are processed, outputting statistics "
                        "such as min, max, and mean." },
        { 0, 0, 0, OPTION_DOC, "Two sets of output files are created, one base file describing the extract as a whole, "
                "located at <ofile> or, if not set, <ifile>" OFILE_DEFAULT_EXT ", and one file for each product requested, "
        "located at <base>.<product>.  If the file is ignored by the L2QC, the product files are omitted." },
        { 0, 0, 0, 0, "Return values:", -1 },
        { str(VALEXTRACT_ERR_NONE) "=" expanded_str(VALEXTRACT_ERR_NONE), 0, 0, OPTION_DOC,
                "Extract successfully processed" },
        { str(VALEXTRACT_ERR_NCFILE_ERR) "=" expanded_str(VALEXTRACT_ERR_NCFILE_ERR), 0, 0, OPTION_DOC,
                "NetCDF file is malformed" },
        { str(VALEXTRACT_ERR_NCFILE_INVALID) "=" expanded_str(VALEXTRACT_ERR_NCFILE_INVALID), 0, 0, OPTION_DOC,
                "NetCDF file is fine but isn't in the format expected (doesn't have geospatial dimensions, etc)" },
        { str(VALEXTRACT_ERR_FLAG) "=" expanded_str(VALEXTRACT_ERR_FLAG), 0, 0, OPTION_DOC,
                "Error processing/finding flags" },
        { str(VALEXTRACT_ERR_VARIABLE) "=" expanded_str(VALEXTRACT_ERR_VARIABLE), 0, 0, OPTION_DOC,
                "Error processing/finding a product" },
        { str(VALEXTRACT_ERR_INPUT) "=" expanded_str(VALEXTRACT_ERR_INPUT), 0, 0, OPTION_DOC,
                "Bad or no input given" },
        { str(VALEXTRACT_ERR_L2QC) "=" expanded_str(VALEXTRACT_ERR_L2QC), 0, 0, OPTION_DOC,
                "L2QC flag over threshold" },
        { str(VALEXTRACT_ERR_UNKNOWN) "=" expanded_str(VALEXTRACT_ERR_UNKNOWN), 0, 0, OPTION_DOC,
                "Unknown error, such as a malloc failure or permissions problem." },
        { 0 } // tell argpar to stop checking options
};

static const char doc[] =
        "val_extract is a utility designed to process a small section of a Level-2 "
                "NetCDF file given either a bounding box or a center point and a box size.  "
                "A product list can be given, as well, to limit processing targets.  sstval_extract "
                "is a proxy for this library that can also find higher quality pixels in order to "
                "re-center the extract.";

static const char args_doc[] = "ifile=<file> <box definition> [products...]";

const char *argpar_program_name = "sstval_extract";
const char *argpar_program_version = "1.3.11";

static int parse_options(int key, char *argv, struct argpar_state *state) {
    int ret = 0;

    sstval_extract_main_input *arguments = state->input;
    switch (key) {
    case 'o':
        arguments->ofile = (void*) argv;
        break;
    case 'q':
        arguments->qual_check = (void*) argv;
        break;
    case 'd':
        if (state->argv_as_dbl_err) {
            eprintf("qual_check_distance must be valid number\n");
            return ARGPAR_ERR_ABORT;
        } else {
            arguments->qual_check_distance = state->argv_as_dbl;
        }
        break;
    case ARGPAR_KEY_ARG:
        arguments->val_extract_arguments->products[arguments->val_extract_arguments->product_count++] = argv;
        break;
    case ARGPAR_KEY_INIT:
        state->child_inputs[0] = arguments->val_extract_arguments;
        arguments->val_extract_arguments->product_count = 0;
        arguments->qual_check = NULL;
        arguments->qual_check_distance = 0;

        arguments->val_extract_arguments->products = malloc(state->argc * sizeof(char*));
        if (arguments->val_extract_arguments->products == NULL) {
            return ARGPAR_ERR_UNKNOWN;
        }
        break;
    }
    return ret;
}

static const argpar_child argpar_children[] = { { &val_extract_argpar }, { 0 } };
static argpar sstval_extract_main_argpar = { options, parse_options, args_doc, doc, argpar_children };

static const char *TOTAL_PIXEL_STR = "pixel_count";
static const char *TOTAL_PIXEL_NOT_FLAGGED_STR = "unflagged_pixel_count";
static const char *TOTAL_PIXEL_FLAGGED_STR = "flagged_pixel_count";
static const char *ASCII_TIME_STR = "time";
static const char *PROD_NAME_STR = "name";
static const char *MAX_PIXEL_STR = "max";
static const char *MIN_PIXEL_STR = "min";
static const char *MEAN_PIXEL_STR = "mean";
static const char *STDDEV_STR = "stddev";
static const char *RMS_STR = "rms";
static const char *VALID_PIXEL_STR = "valid_pixel_count";
static const char *MEDIAN_PIXEL_STR = "median";
static const char *CENTER_PIXEL_STR = "center_value";

static void print_stats_to_file(FILE *h, const char *prefix, nc_var_stats *stats) {
    fprintf(h, "%s%s=%d\n", prefix, VALID_PIXEL_STR, stats->count);
    fprintf(h, "%s%s=%f\n", prefix, MAX_PIXEL_STR, stats->max);
    fprintf(h, "%s%s=%f\n", prefix, MIN_PIXEL_STR, stats->min);
    fprintf(h, "%s%s=%f\n", prefix, MEAN_PIXEL_STR, stats->mean);
    fprintf(h, "%s%s=%f\n", prefix, MEDIAN_PIXEL_STR, stats->median);
    fprintf(h, "%s%s=%f\n", prefix, STDDEV_STR, stats->stddev);
    fprintf(h, "%s%s=%f\n", prefix, RMS_STR, stats->rms);
}

static unsigned count_char(const char *haystack, char needle) {
    unsigned i, count;
    for (i = 0, count = 0; haystack[i]; i++) {
        count += (haystack[i] == needle);
    }
    return count;
}
static int strcmp_p(const void *p1, const void *p2) {
    unsigned p1_slashes = count_char(*(char * const *) p1, '/');
    unsigned p2_slashes = count_char(*(char * const *) p2, '/');
    if (p1_slashes < p2_slashes) {
        return -1;
    } else if (p1_slashes > p2_slashes) {
        return 1;
    }
    return strcmp(*(char * const *) p1, *(char * const *) p2);
}

static int save_extract(int key, void *nc_input, void *user_input) {
    sstval_extract_main_input *input = user_input;
    switch (key) {
    case VALEXTRACT_KEY_INIT:
        break;
    case VALEXTRACT_KEY_SUCCESS:
        break;
    case VALEXTRACT_KEY_ERROR:
        break;
    case VALEXTRACT_KEY_FINI:
        if (input->base_ofile_h) {
            fprintf(input->base_ofile_h, "\n");
            fprintf(input->base_ofile_h, "versions=val_extract=%s sstval_extract-main=%s\n", val_extract_version(), argpar_program_version);
            fclose(input->base_ofile_h);
        }
        if (input->products_printed) {
            shash_destroy(input->products_printed);
        }
        break;
    case VALEXTRACT_KEY_FILE: {
        nc_region *region = nc_input;
        nc_file *file = region->file;
        FILE *ofile_h = NULL;
        if ((ofile_h = fopen(input->ofile, "w")) == NULL) {
            eprintf("Error opening output file %s, %s\n", input->ofile, strerror(errno));
            return -1;
        }

        fprintf(ofile_h, "%s=%s\n", ASCII_TIME_STR, region->ascii_time);
        fprintf(ofile_h, "utime=%zd\n", region->utime);
        fprintf(ofile_h, "dim_length_lines=%d\n", file->dim_lengths[file->line_dimid]);
        fprintf(ofile_h, "dim_length_pixels=%d\n", file->dim_lengths[file->pixel_dimid]);
        if (input->val_extract_arguments->start_lat != VALEXTRACT_UNSET) {
            if (input->val_extract_arguments->box_size == VALEXTRACT_UNSET) {
                fprintf(ofile_h, "slat=%f\n", input->val_extract_arguments->start_lat);
                fprintf(ofile_h, "slon=%f\n", input->val_extract_arguments->start_lon);
                fprintf(ofile_h, "elat=%f\n", input->val_extract_arguments->end_lat);
                fprintf(ofile_h, "elon=%f\n", input->val_extract_arguments->end_lon);
            } else {
                fprintf(ofile_h, "slat=%f\n", input->slat_orig);
                fprintf(ofile_h, "slon=%f\n", input->slon_orig);
                fprintf(ofile_h, "qual_check=%s\n", input->qual_check);
                fprintf(ofile_h, "qual_check_distance=%f\n", input->qual_check_distance);
                fprintf(ofile_h, "qual_check_distance_actual=%f\n", vincenty_distance(input->slat_orig, input->slon_orig, input->val_extract_arguments->start_lat, input->val_extract_arguments->start_lon) / 1000.);
            }
        } else {
            fprintf(ofile_h, "sline=%d\n", input->val_extract_arguments->start_line);
            fprintf(ofile_h, "spixl=%d\n", input->val_extract_arguments->start_pixel);
            if (input->val_extract_arguments->box_size == VALEXTRACT_UNSET) {
                fprintf(ofile_h, "eline=%d\n", input->val_extract_arguments->end_line);
                fprintf(ofile_h, "epixl=%d\n", input->val_extract_arguments->end_pixel);
            }
        }
        if (input->val_extract_arguments->box_size != VALEXTRACT_UNSET) {
            fprintf(ofile_h, "boxsize=%d\n", input->val_extract_arguments->box_size);
        }
        fprintf(ofile_h, "center_lat=%f\n", region->lat);
        fprintf(ofile_h, "center_lon=%f\n", region->lon);
        fprintf(ofile_h, "center_line=%zd\n", region->center.line);
        fprintf(ofile_h, "center_pixel=%zd\n", region->center.pixel);
        fprintf(ofile_h, "%s=%d\n", TOTAL_PIXEL_STR, region->pixel_count);
        fprintf(ofile_h, "%s=%d\n", TOTAL_PIXEL_NOT_FLAGGED_STR, region->unflagged_pixel_count);
        fprintf(ofile_h, "%s=%d\n", TOTAL_PIXEL_FLAGGED_STR, region->flagged_pixel_count);
        fprintf(ofile_h, "flag_counts=");
        int flag_i;
        for (flag_i = 0; flag_i < file->flag_count; flag_i++) {
            if (strcmp("SPARE", file->flag_meanings[flag_i])) {
                fprintf(ofile_h, "%*s%s=%d", flag_i ? 1 : 0, "", file->flag_meanings[flag_i], region->flag_counts[flag_i]);
            }
        }
        fprintf(ofile_h, "\nflag_percents=");
        for (flag_i = 0; flag_i < file->flag_count; flag_i++) {
            if (strcmp("SPARE", file->flag_meanings[flag_i])) {
                double flag_percent = (double) region->flag_counts[flag_i] / region->pixel_count;
                fprintf(ofile_h, "%*s%s=%f", flag_i ? 1 : 0, "", file->flag_meanings[flag_i], flag_percent);
            }
        }
        fprintf(ofile_h, "\nvariables=");
        input->base_ofile_h = ofile_h;
        input->product_count = 0;
        if ((input->products_printed = shash_create(0)) == NULL) {
            eprintf("Error creating product hash\n");
            return -1;
        }

        if (file->attributes) {
            char attr_filename[255];
            strcpy(attr_filename, input->ofile);
            strcat(attr_filename, ".global_attrs");
            FILE *attr_file_h = NULL;
            if ((attr_file_h = fopen(attr_filename, "w")) == NULL) {
                eprintf("Error opening output file %s, %s\n", attr_filename, strerror(errno));
                return -1;
            }
            int attr_count = shash_rewind(file->attributes), att_i = 0;
            const char *attribute_names[attr_count], *v;
            while (!shash_next(file->attributes, &attribute_names[att_i++], &v))
                ;
            qsort(attribute_names, attr_count, sizeof(char*), strcmp_p);
            for (att_i = 0; att_i < attr_count; att_i++) {
                fprintf(attr_file_h, "%s=%s\n", attribute_names[att_i], shash_get(file->attributes, attribute_names[att_i]));
            }
            fclose(attr_file_h);
        }
    }
        break;
    case VALEXTRACT_KEY_VAR: {
        nc_var *var = nc_input;
        nc_region *region = var->region;
        nc_file *file = var->file;

        size_t center = SIZE_MAX;

        if (!var->is_geospatial) {
            if (var->dim_ids[0] == file->line_dimid) {
                center = region->center.line;
            } else if (var->dim_ids[0] == file->pixel_dimid) {
                center = region->center.pixel;
            }
        }

        if (var->is_geospatial || center != SIZE_MAX) {
            if (input->product_count++) {
                fprintf(input->base_ofile_h, ",");
            }
            char output_filename[255];
            strcpy(output_filename, input->ofile);
            strcat(output_filename, ".");
            if (shash_get(input->products_printed, var->name) == NULL) {
                shash_set(input->products_printed, var->name, "");
                strcat(output_filename, var->name);
                fprintf(input->base_ofile_h, "%s", var->name);
            } else {
                char var_filename[strlen(var->name) + 4];
                sprintf(var_filename, "%s_", var->name);
                char *where_to_write = var_filename + strlen(var->name) + 1;
                unsigned i = 1;
                do {
                    sprintf(where_to_write, "%u", ++i);
                } while (shash_get(input->products_printed, var_filename) != NULL);
                shash_set(input->products_printed, var_filename, "");
                strcat(output_filename, var_filename);
                fprintf(input->base_ofile_h, "%s", var_filename);
            }

            FILE *output_h = NULL;
            if ((output_h = fopen(output_filename, "w")) == NULL) {
                eprintf("Error opening output file %s, %s\n", output_filename, strerror(errno));
                return -1;
            }

            fprintf(output_h, "%s=%s\n", PROD_NAME_STR, var->name);
            fprintf(output_h, "group_name=%s\n", var->group_name);

            if (var->attributes) {
                int attr_count = shash_rewind(var->attributes), att_i = 0;
                const char *attribute_names[attr_count], *v;
                while (!shash_next(var->attributes, &attribute_names[att_i++], &v))
                    ;
                qsort(attribute_names, attr_count, sizeof(char*), strcmp_p);
                for (att_i = 0; att_i < attr_count; att_i++) {
                    fprintf(output_h, "%s=%s\n", attribute_names[att_i], shash_get(var->attributes, attribute_names[att_i]));
                }
            }

            if (var->is_geospatial) {
                fprintf(output_h, "%s=%f\n", CENTER_PIXEL_STR, var->stats.center_value);
                if (var->valid_data){
                    print_stats_to_file(output_h, "", &var->stats);
                    print_stats_to_file(output_h, "filtered_", &var->filtered_stats);
                    print_stats_to_file(output_h, "iqr_", &var->iqr_stats);
                }
            } else {
                if (center != SIZE_MAX) {
                    double v = 0;
                    int nc_ret;
                    if ((nc_ret = nc_get_var1_double(var->gid, var->varid, &center, &v)) != NC_NOERR) {
                        printf("Error getting center for %s, error %d\n", var->name, nc_ret);
                        return -1;
                    }
                    fprintf(output_h, "%s=%f\n", CENTER_PIXEL_STR, v);
                }
            }
            fclose(output_h);
        }
    }
        break;
    }
    return 0;
}

static int find_best_quality(int key, void *nc_input, void *user_input) {
    sstval_extract_main_input *input = user_input;
    switch (key) {
    case VALEXTRACT_KEY_INIT:
        break;
    case VALEXTRACT_KEY_SUCCESS: {
        nc_var *var = input->qual_check_var_copies[2];
        nc_file *file = var->file;

        const int line_length = file->dim_lengths[file->line_dimid];
        const int pixl_length = file->dim_lengths[file->pixel_dimid];

        const int border_bottom = (input->box_size >> 1);
        const int line_border_top = line_length - border_bottom;
        const int pixl_border_top = pixl_length - border_bottom;

        const double max_distance = input->qual_check_distance * 1000;

        double *lat = input->qual_check_data_copies[0];
        double *lon = input->qual_check_data_copies[1];
        double *qual = input->qual_check_data_copies[2];

        int best_i = -1;
        double best_quality = DBL_MAX, closest_best_quality = DBL_MAX;

        int line = 0, pixel = 0;

        for (int i = 0; i < var->region->pixel_count; i++) {
            if (lat[i] != -999 && qual[i] >= 0 && (var->region->pixel_flags[i] & input->val_extract_arguments->ignore_flags_mask) == 0) {
                if (!(line < border_bottom || pixel < border_bottom || line > line_border_top || pixel > pixl_border_top)) {
                    const double distance = vincenty_distance(input->slat_orig, input->slon_orig, lat[i], lon[i]);
                    if (distance <= max_distance && (qual[i] < best_quality || (qual[i] == best_quality && distance < closest_best_quality))) {
                        best_quality = qual[i];
                        closest_best_quality = distance;
                        best_i = i;
                    }
                }
            }
            ++pixel;
            if (pixel >= pixl_length) {
                pixel = 0;
                ++line;
            }
        }

        if (best_i >= 0) {
            input->val_extract_arguments->start_lat = lat[best_i];
            input->val_extract_arguments->start_lon = lon[best_i];
        } else {
            input->val_extract_arguments->start_lat = DBL_MAX;
            input->val_extract_arguments->start_lon = DBL_MAX;
        }
//			printf("find_best: %f, %f (%f)\n", input->val_extract_arguments->start_lat, input->val_extract_arguments->start_lon, best_quality);
//			printf("find_best distance: %f\n", vincenty_distance(input->val_extract_arguments->start_lat, input->val_extract_arguments->start_lon, var->region->lat, var->region->lon));
    }
        break;
    case VALEXTRACT_KEY_ERROR:
        break;
    case VALEXTRACT_KEY_FINI:
        break;
    case VALEXTRACT_KEY_FILE: {
//		nc_region *region = nc_input;
    }
        break;
    case VALEXTRACT_KEY_VAR: {
        nc_var *var = nc_input;
        int copy_location = 2;
        if (!strcmp(var->name, "latitude")) {
            copy_location = 0;
        } else if (!strcmp(var->name, "longitude")) {
            copy_location = 1;
        }
        input->qual_check_var_copies[copy_location] = malloc(sizeof(nc_var));
        memcpy(input->qual_check_var_copies[copy_location], var, sizeof(nc_var));
        input->qual_check_data_copies[copy_location] = malloc(sizeof(double) * var->region->pixel_count);
        memcpy(input->qual_check_data_copies[copy_location], var->data, sizeof(double) * var->region->pixel_count);
    }
        break;
    }
    return 0;
}

int main(int argc, char *argv[]) {

    if (argc <= 1) {
        argpar_usage_default(&sstval_extract_main_argpar);
        return VALEXTRACT_ERR_INPUT;
    }
    sstval_extract_main_input input = { 0 };
    val_extract_arguments arguments = {
            .geospatial_only = 1, .geospatial_to_double = 1,
            .val_extract_parser = (val_extract_parser) save_extract,
            .user_input = (void*) &input
    };
    input.val_extract_arguments = &arguments;

    olog *olog_default = olog_load_default();

    int ret = VALEXTRACT_ERR_NONE;
    if (argpar_parse_args(&sstval_extract_main_argpar, argc, argv, 0, NULL, &input)) {
        ret = VALEXTRACT_ERR_INPUT;
    } else {
        input.slat_orig = arguments.start_lat;
        input.slon_orig = arguments.start_lon;

        if (input.qual_check) {
            double box_size_km_save = arguments.box_size_km;
            int box_size_save = arguments.box_size;
            int product_count_save = arguments.product_count;
            int ignore_flag_count_save = arguments.ignore_flag_count;
            int count_flag_count_save = arguments.count_flag_count;
            int l2qc_threshold_count_save = arguments.l2qc_threshold_count;
            int valid_range_count_save = arguments.valid_range_count;
            char **products_save = arguments.products;
            char **ignore_flags_save = arguments.ignore_flags;
            char **count_flags_save = arguments.count_flags;
            char **l2qc_flags_save = arguments.l2qc_flags;
            double *l2qc_thresholds_save = arguments.l2qc_thresholds;
            const char *ifile_save = arguments.ifile;
            bool global_atts_save = arguments.global_atts;
            bool variable_atts_save = arguments.variable_atts;
            val_extract_valid_range *valid_ranges_save = arguments.valid_ranges;

            char *qual_products[] = { "latitude", "longitude", input.qual_check };
            arguments.product_count = 3;
            arguments.products = qual_products;
//			arguments.box_size_km = input.qual_check_distance;
//			arguments.box_size = -1;
            arguments.start_lat = -90;
            arguments.end_lat = 90;
            arguments.start_lon = -180;
            arguments.end_lon = 180;
//			arguments.ignore_flags = NULL;
            arguments.count_flags = NULL;
            arguments.l2qc_flags = NULL;
            arguments.l2qc_thresholds = NULL;
            arguments.valid_ranges = NULL;
            arguments.global_atts = false;
            arguments.variable_atts = false;

//			arguments.ignore_flag_count = 0;
            arguments.count_flag_count = 0;
            arguments.l2qc_threshold_count = 0;
            arguments.valid_range_count = 0;

            arguments.val_extract_parser = (val_extract_parser) find_best_quality;

            input.box_size = box_size_save;
            input.qual_check_var_copies = malloc(sizeof(nc_var*) * 3);
            input.qual_check_data_copies = calloc(sizeof(double*), 3);

            ret = val_extract(&arguments);
//			arguments.products = NULL;
//			val_extract_clean(&arguments);

            for (int i = 0; i < 3; i++) {
                if (input.qual_check_data_copies[i]) {
                    free(input.qual_check_data_copies[i]);
                    free(input.qual_check_var_copies[i]);
                }
            }
            free(input.qual_check_data_copies);
            free(input.qual_check_var_copies);

            double start_lat_save = arguments.start_lat;
            double start_lon_save = arguments.start_lon;

            val_extract_clear_args(&arguments);

            arguments.products = products_save;
            arguments.ignore_flags = ignore_flags_save;
            arguments.count_flags = count_flags_save;
            arguments.l2qc_flags = l2qc_flags_save;
            arguments.valid_ranges = valid_ranges_save;
            arguments.l2qc_thresholds = l2qc_thresholds_save;
            arguments.ifile = ifile_save;
            arguments.lat_and_lon = true;
            arguments.product_count = product_count_save;
            arguments.box_size = box_size_save;
            arguments.box_size_km = box_size_km_save;
            arguments.ignore_flag_count = ignore_flag_count_save;
            arguments.count_flag_count = count_flag_count_save;
            arguments.l2qc_threshold_count = l2qc_threshold_count_save;
            arguments.valid_range_count = valid_range_count_save;
            arguments.start_lat = start_lat_save;
            arguments.start_lon = start_lon_save;
            arguments.global_atts = global_atts_save;
            arguments.variable_atts = variable_atts_save;
            arguments.val_extract_parser = (val_extract_parser) save_extract;

            if (ret != VALEXTRACT_ERR_NONE || start_lat_save == DBL_MAX) {
                olog_crit(olog_default, "Error finding closest best-quality pixel\n");
                if (ret == VALEXTRACT_ERR_NONE) {
                    olog_crit(olog_default, "val-extract succeeded, a close-enough point was likely not found in file\n");
                    ret = VALEXTRACT_ERR_POINT_NOT_FOUND;
                }
            }
        }
        if (ret == VALEXTRACT_ERR_NONE) {
            char output_filename[255];
            if (!input.ofile) {
                strcpy(output_filename, arguments.ifile);
                strcat(output_filename, OFILE_DEFAULT_EXT);
                input.ofile = output_filename;
            }
            ret = val_extract(&arguments);
        }
    }
    val_extract_clean(&arguments);
    argpar_clean(&sstval_extract_main_argpar);

    olog_destroy(olog_default);

    return ret;
}
