#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "smigen_input.h"
#include "smiinc.h"
#include "palette.h"
#include "genutils.h"

int get_item(char *arg, instr *input_str);
void set_param_string(instr *input_str);
int par_file(char *pfile, instr *input_str);
int get_def_l2prod(char *def_file, instr *input_str);

int input_init(instr *input_str) {
    input_str->ifile[0] = '\0';
    input_str->ofile[0] = '\0';
    input_str->pfile[0] = '\0';
    input_str->prod [0] = '\0';
    input_str->parms[0] = '\0';

    input_str->stype = 0;
    input_str->meas = 1;
    input_str->datamin = 0.0;
    input_str->datamax = 0.0;

    input_str->lonwest = -180.0;
    input_str->loneast = +180.0;
    input_str->latsouth = -90.0;
    input_str->latnorth = +90.0;

    strcpy(input_str->pversion, "Unspecified");
    strcpy(input_str->projection, "RECT");
    strcpy(input_str->resolution, "9km");
    strcpy(input_str->palfile, "DEFAULT");

    input_str->gap_fill = 0;
    input_str->seam_lon = -180;

    input_str->proddesc[0] = '\0';
    input_str->units[0] = '\0';
    input_str->precision[0] = '\0';

    input_str->minobs = 0;
    input_str->deflate = 0;
    input_str->oformat[0] = '\0';
    return 0;
}

/*-----------------------------------------------------------------------------
    Function:  smigen_input

    Returns:   int (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        Convert the arguments from the command line into a structure input
        variable.

    Parameters: (in calling order)
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       argc          I   number of arguments
        char      **argv        I   list of arguments
        instr     input         O   structure variable for inputs

----------------------------------------------------------------------------*/

int smigen_input(int argc, char **argv, instr *input) {
    int i;

    if (argc == 1) return (-1);

    /*                                                                  */
    /* Set input values to defaults                                     */
    /*                                                                  */
    if (input_init(input) != 0) {
        printf("-E- %s: Error initializing input structure.\n", __FILE__);
        return (-1);
    }

    /*                                                                  */
    /* Loop through command arguments and update input struct           */
    /*                                                                  */
    for (i = 1; i < argc; i++)
        if (get_item(argv[i], input) != 0)
            return -1;


#if 0
    /* Define default palette */
    /* ---------------------- */
    if (strcmp(input->prod, "NDVI") == 0 ||
            strcmp(input->prod, "EVI") == 0)
        memcpy(input->palette, ndvi_palette, 768);
    else if (strcmp(input->prod, "sst") == 0 ||
            strcmp(input->prod, "sst4") == 0)
        memcpy(input->palette, sst_palette, 768);
    else
        memcpy(input->palette, default_palette, 768);


    /*                                                                  */
    /* Add any input-dependent defaults                                 */
    /*                                                                  */
    if (strcmp(input->palfile, "DEFAULT") != 0) {
        short *r, *g, *b;
        if (!(r = (short *) calloc(256, sizeof (short)))) {
            fprintf(stderr, "smigen: Error allocating space for red.\n");
            return -1;
        };
        if (!(g = (short *) calloc(256, sizeof (short)))) {
            fprintf(stderr, "smigen: Error allocating space for green.\n");
            return -1;
        };
        if (!(b = (short *) calloc(256, sizeof (short)))) {
            fprintf(stderr, "smigen: Error allocating space for blue.\n");
            return -1;
        };
        if (getlut_file(input->palfile, r, g, b)) {
            fprintf(stderr, "Error reading palette file %s\n", input->pfile);
            return -1;
        }
        for (i = 0; i < 256; i++) {
            input->palette[i * 3] = r[i];
            input->palette[i * 3 + 1] = g[i];
            input->palette[i * 3 + 2] = b[i];
        }
        free(r);
        free(g);
        free(b);
    }
#endif
    return 0;
}

void set_param_string(instr *input) {

    char str_buf[4096];
    /*                                                                  */
    /* Build string of parameters for metadata                          */
    /*                                                                  */
    sprintf(str_buf, "ifile = %s|", input->ifile);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "ofile = %s|", input->ofile);
    strcat(input->parms, str_buf);
    if (input->pfile[0] != '\0') {
        sprintf(str_buf, "pfile = %s|", input->pfile);
        strcat(input->parms, str_buf);
    }
    sprintf(str_buf, "prod = %s|", input->prod);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "palfile = %s|", input->palfile);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "processing version = %s|", input->pversion);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "meas = %d|", input->meas);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "stype = %d|", input->stype);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "datamin = %f|", input->datamin);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "datamax = %f|", input->datamax);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "lonwest = %f|", input->lonwest);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "loneast = %f|", input->loneast);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "latsouth = %f|", input->latsouth);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "latnorth = %f|", input->latnorth);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "resolution = %s|", input->resolution);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "projection = %s|", input->projection);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "gap_fill = %d|", input->gap_fill);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "seam_lon = %f|", input->seam_lon);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "minobs = %d|", input->minobs);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "deflate = %d|", input->deflate);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "oformat = %s|", input->oformat);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "precision = %s|", input->precision);
    strcat(input->parms, str_buf);
}

int get_item(char *arg, instr *input) {
    char *tmp_str;
    char keyword [20];
    char parm_str[4096];
    char tmp_file[FILENAME_MAX];
    int ilen1, iret;
    static int lon_limit_used = 0;
    static int seam_used = 0;

    if ((tmp_str = strchr(arg, '=')) == NULL) {
        printf("Invalid argument \"%s\"\n", arg);
        return -1;
    }

    ilen1 = tmp_str - arg;
    strncpy(keyword, arg, ilen1);
    keyword[ilen1] = '\0';
    strcpy(parm_str, tmp_str + 1);

    /* change keyword to lower case */
    tmp_str = keyword;
    while (*tmp_str != '\0') {
        if (isupper(*tmp_str)) *tmp_str = tolower(*tmp_str);
        tmp_str++;
    }

    if (strcmp(keyword, "par") == 0) {
        iret = par_file(parm_str, input);
        return iret;

    } else if (strcmp(keyword, "ifile") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->ifile, tmp_file);

    } else if (strcmp(keyword, "ofile") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->ofile, tmp_file);

    } else if (strcmp(keyword, "prod") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->prod, tmp_file);

    } else if (strcmp(keyword, "palfile") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->palfile, tmp_file);

    } else if (strcmp(keyword, "pversion") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->pversion, tmp_file);

    } else if (strcmp(keyword, "stype") == 0) {
        input->stype = atoi(parm_str);

    } else if (strcmp(keyword, "meas") == 0) {
        input->meas = atoi(parm_str);
        if (input->meas == 4) strcpy(input->prod, "pixels");
        if (input->meas == 5) strcpy(input->prod, "scenes");

    } else if (strcmp(keyword, "datamin") == 0) {
        input->datamin = atof(parm_str);

    } else if (strcmp(keyword, "datamax") == 0) {
        input->datamax = atof(parm_str);

    } else if (strcmp(keyword, "lonwest") == 0) {
        input->lonwest = atof(parm_str);
        lon_limit_used = 1;

    } else if (strcmp(keyword, "loneast") == 0) {
        input->loneast = atof(parm_str);
        lon_limit_used = 1;

    } else if (strcmp(keyword, "latsouth") == 0) {
        input->latsouth = atof(parm_str);

    } else if (strcmp(keyword, "latnorth") == 0) {
        input->latnorth = atof(parm_str);

    } else if (strcmp(keyword, "resolution") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->resolution, tmp_file);

    } else if (strcmp(keyword, "projection") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->projection, tmp_file);

    } else if (strcmp(keyword, "gap_fill") == 0) {
        input->gap_fill = atoi(parm_str);

    } else if (strcmp(keyword, "seam_lon") == 0) {
        input->seam_lon = atof(parm_str);
        seam_used = 1;

    } else if (strcmp(keyword, "proddesc") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->proddesc, tmp_file);

    } else if (strcmp(keyword, "units") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->units, tmp_file);

    } else if (strcmp(keyword, "precision") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->precision, tmp_file);

    } else if (strcmp(keyword, "minobs") == 0) {
        input->minobs = atol(parm_str);

    } else if (strcmp(keyword, "deflate") == 0) {
        input->deflate = atoi(parm_str);

    } else if (strcmp(keyword, "oformat") == 0) {
        const char* tmpStr = getFileFormatName(parm_str);
        strcpy(input->oformat, tmpStr);

    } else {
        goto Invalid_return;
    }

    if (lon_limit_used == 1 && seam_used == 1) {
        printf("LONEAST/LONWEST cannot be used with SEAM_LON.\n");
        exit(1);
    }

    return 0;


Invalid_return:
    printf("Invalid argument \"%s\"\n", arg);
    return -1;
}

int par_file(char *pfile, instr *input) {
    FILE *fp;
    char arg[2048];
    int iret;

    if ((fp = fopen(pfile, "r")) == NULL) {
        printf("Error on opening the parameter file - %s\n", pfile);
        return -1;
    }

    while ((fgets(arg, 2047, fp)) != NULL) {
        /* skip the comment or blank line */
        if (arg[0] == '#' || arg[0] == ';' || arg[0] == ' ' || arg[0] == '\0' ||
                arg[0] == '\n')
            continue;

        arg[strlen(arg) - 1] = '\0'; /* replace the last char new_line to NULL */
        iret = get_item(arg, input);
        if (iret != 0) {
            fclose(fp);
            return -1;
        }
    }

    fclose(fp);
    return 0;
}
