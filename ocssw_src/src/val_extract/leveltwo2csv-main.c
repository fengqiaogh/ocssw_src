#include "argpar.h"
#include "val_extract.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define OFILE_DEFAULT_EXT ".csv"

#ifdef DEBUG
#define dprintf(...) do { printf(__VA_ARGS__); } while(0)
#else
#define dprintf(...) do {} while (0)
#endif

#define eprintf(...) do { fprintf(stderr, __VA_ARGS__); } while(0)

#define str(s) #s
#define expanded_str(s) str(s)

typedef struct leveltwo2csv_main_input {
    char *ofile;
    char *product;
    val_extract_arguments *val_extract_arguments;
} leveltwo2csv_main_input;

static argpar_option options[] = {
        { "ofile", 'o', "FILE", 0, "output file path" },
        { "product", 'p', "FILE", 0, "product to dump" },
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
        "A wrapper around val-extract to dump a product of an L2 file to a CSV";

static const char args_doc[] = "ifile=<file> [<box definition>] [products...]";

static int parse_options(int key, char *argv, struct argpar_state *state) {
    int ret = 0;

    leveltwo2csv_main_input *arguments = state->input;
    switch (key) {
    case 'o':
        arguments->ofile = (void*) argv;
        break;
    case 'p':
        arguments->product = (void*) argv;
        break;
    case ARGPAR_KEY_SUCCESS:
        arguments->val_extract_arguments->products = malloc((arguments->product ? 3 : 2) * sizeof(char*));
        if (arguments->val_extract_arguments->products == NULL) {
            return ARGPAR_ERR_UNKNOWN;
        }
        arguments->val_extract_arguments->products[0] = "latitude";
        arguments->val_extract_arguments->products[1] = "longitude";
        if (arguments->product) {
            arguments->val_extract_arguments->products[2] = arguments->product;
            arguments->val_extract_arguments->product_count = 3;
        } else {
            arguments->val_extract_arguments->product_count = 2;
        }
        break;
    case ARGPAR_KEY_INIT:
        state->child_inputs[0] = arguments->val_extract_arguments;
        arguments->val_extract_arguments->product_count = 0;
        arguments->val_extract_arguments->box_size = 0;
        arguments->val_extract_arguments->start_lat = -90;
        arguments->val_extract_arguments->end_lat = 90;
        arguments->val_extract_arguments->start_lon = -180;
        arguments->val_extract_arguments->end_lon = 180;
        break;
    }
    return ret;
}

static const argpar_child argpar_children[] = { { &val_extract_argpar }, { 0 } };
static argpar leveltwo2csv_main_argpar = { options, parse_options, args_doc, doc, argpar_children };

const char *argpar_program_name = "leveltwo2csv";
const char *argpar_program_version = "0.0.1";

static int save_extract(int key, void *nc_input, void *user_input) {
    leveltwo2csv_main_input *input = user_input;
    switch (key) {
    case VALEXTRACT_KEY_INIT:
        break;
    case VALEXTRACT_KEY_SUCCESS: {
        FILE *ofile_h = NULL;
        if ((ofile_h = fopen(input->ofile, "w")) == NULL) {
            eprintf("Error opening output file %s, %s\n", input->ofile, strerror(errno));
            return -1;
        }
        FILE *ivar_h[input->val_extract_arguments->product_count];
        int i;
        char input_filename[3][255];
        for (i = 0; i < input->val_extract_arguments->product_count; i++) {
            strcpy(input_filename[i], input->ofile);
            strcat(input_filename[i], ".");
            strcat(input_filename[i], input->val_extract_arguments->products[i]);
            if ((ivar_h[i] = fopen(input_filename[i], "r")) == NULL) {
                eprintf("Error opening output file %s, %s\n", input_filename[i], strerror(errno));
                return -1;
            }
            if (i) {
                fprintf(ofile_h, ",");
            }
            fprintf(ofile_h, "%s", input->val_extract_arguments->products[i]);
        }

        int file_ret = 0;
        while (file_ret != EOF) {
            fprintf(ofile_h, "\n");
            for (i = 0; i < input->val_extract_arguments->product_count && file_ret != EOF; i++) {
                double v = NAN;
                file_ret = fscanf(ivar_h[i], "%lf\n", &v);
                if (file_ret != EOF) {
                    if (i) {
                        fprintf(ofile_h, ",");
                    }
                    fprintf(ofile_h, "%f", v);
                }
            }
        }
        fprintf(ofile_h, "\n");

        for (i = 0; i < input->val_extract_arguments->product_count; i++) {
            fclose(ivar_h[i]);
            remove(input_filename[i]);
        }
        fclose(ofile_h);
    }
        break;
    case VALEXTRACT_KEY_ERROR:
        break;
    case VALEXTRACT_KEY_FINI:
        break;
    case VALEXTRACT_KEY_VAR: {
        nc_var *var = nc_input;
        if (!var->is_geospatial) {
            return 0;
        }
        char output_filename[255];
        strcpy(output_filename, input->ofile);
        strcat(output_filename, ".");
        strcat(output_filename, var->name);

        FILE *output_h = NULL;
        if ((output_h = fopen(output_filename, "w")) == NULL) {
            eprintf("Error opening output file %s, %s\n", output_filename, strerror(errno));
            return -1;
        }

        int i;
        double *values = (double*) var->data;
        for (i = 0; i < var->region->pixel_count; i++) {
            if (var->has_fill && values[i] == var->fill) {
                fprintf(output_h, "%s\n", "NULL");
            } else {
                double v = values[i];
                if (var->has_scale) {
                    v *= var->scale;
                }
                if (var->has_offset) {
                    v += var->offset;
                }
                fprintf(output_h, "%f\n", v);
            }
        }
        fclose(output_h);
    }
        break;
    }
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc <= 1) {
        argpar_usage_default(&leveltwo2csv_main_argpar);
        return VALEXTRACT_ERR_INPUT;
    }
    leveltwo2csv_main_input input = { 0 };
    val_extract_arguments arguments = {
            .geospatial_only = 1, .geospatial_to_double = 1,
            .val_extract_parser = (val_extract_parser) save_extract,
            .user_input = (void*) &input
    };
    input.val_extract_arguments = &arguments;

    int ret = EXIT_SUCCESS;
    if (argpar_parse_args(&leveltwo2csv_main_argpar, argc, argv, 0, NULL, &input)) {
        ret = VALEXTRACT_ERR_INPUT;
    } else {
        char output_filename[255];
        printf("%s\n", arguments.ifile);
        if (!input.ofile) {
            strcpy(output_filename, arguments.ifile);
            strcat(output_filename, OFILE_DEFAULT_EXT);
            input.ofile = output_filename;
        }
        ret = val_extract(&arguments);
    }
    val_extract_clean(&arguments);
    argpar_clean(&leveltwo2csv_main_argpar);
    return ret;
}

