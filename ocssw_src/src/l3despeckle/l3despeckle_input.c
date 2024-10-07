#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "l3despeckle_input.h"

int input_init(instr *input_str) {
    input_str->infile[0] = '\0';
    input_str->ofile[0] = '\0';
    input_str->pfile[0] = '\0';

    strcpy(input_str->out_parm, "DEFAULT");
    strcpy(input_str->repimg, "ORIGINAL");

    input_str->tflag = 'O';

    input_str->syear = 9999;
    input_str->sday = 999;

    input_str->eyear = 9999;
    input_str->eday = 999;

    input_str->sorbit = -1;
    input_str->eorbit = -1;

    input_str->reduce_fac = 1;

    input_str->noext = 0;

    input_str->merged[0] = '\0';

    input_str->loneast = +180;
    input_str->lonwest = -180;
    input_str->latnorth = +90;
    input_str->latsouth = -90;

    input_str->verbose = 0;

    return 0;
}

/*-----------------------------------------------------------------------------
    Function:  l3despeckle_input

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

int l3despeckle_input(int argc, char **argv, instr *input) {
    int i;
    char str_buf[4096];

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


    /*                                                                  */
    /* Build string of parameters for metadata                          */
    /*                                                                  */
    sprintf(str_buf, "INFILE = %s|", input->infile);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "OFILE = %s|", input->ofile);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "PFILE = %s|", input->ofile);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "SYEAR = %ld|", input->syear);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "EYEAR = %ld|", input->eyear);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "SDAY = %ld|", input->sday);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "EDAY = %ld|", input->eday);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "SORBIT = %ld|", input->sorbit);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "EORBIT = %ld|", input->eorbit);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "OUT_PARM = %s|", input->out_parm);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "TFLAG = %c|", input->tflag);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "REPIMG = %s|", input->repimg);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "REDUCE_FAC = %ld|", input->reduce_fac);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "MERGED = %s|", input->merged);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "LONEAST = %f|", input->loneast);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "LONWEST = %f|", input->lonwest);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "LATNORTH = %f|", input->latnorth);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "LATSOUTH = %f|", input->latsouth);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "VERBOSE = %ld|", input->verbose);
    strcat(input->parms, str_buf);

    return 0;
}

int get_item(char *arg, instr *input) {
    char *tmp_str;
    char keyword [20];
    char parm_str[4096];
    char tmp_file[FILENAME_MAX];
    int ilen1, iret;

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

    if (strncmp(keyword, "parfile", 3) == 0) {
        iret = par_file(parm_str, input);
        return iret;

    } else if (strcmp(keyword, "in") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->infile, tmp_file);

    } else if (strcmp(keyword, "out") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->ofile, tmp_file);

    } else if (strcmp(keyword, "infile") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->infile, tmp_file);

    } else if (strcmp(keyword, "ofile") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->ofile, tmp_file);

    } else if (strcmp(keyword, "pfile") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->pfile, tmp_file);

    } else if (strcmp(keyword, "rep") == 0) {
        strcpy(input->repimg, parm_str);

    } else if (strcmp(keyword, "tflag") == 0) {
        memcpy(&input->tflag, parm_str, 1);

    } else if (strcmp(keyword, "syear") == 0) {
        input->syear = atoi(parm_str);

    } else if (strcmp(keyword, "eyear") == 0) {
        input->eyear = atoi(parm_str);

    } else if (strcmp(keyword, "sday") == 0) {
        input->sday = atoi(parm_str);

    } else if (strcmp(keyword, "eday") == 0) {
        input->eday = atoi(parm_str);

    } else if (strcmp(keyword, "orbit1") == 0) {
        input->sorbit = atoi(parm_str);

    } else if (strcmp(keyword, "orbit2") == 0) {
        input->eorbit = atoi(parm_str);

    } else if (strcmp(keyword, "gsfcqual") == 0) {

    } else if (strcmp(keyword, "out_parm") == 0) {
        strcpy(input->out_parm, ":");
        strcat(input->out_parm, parm_str);
        strcat(input->out_parm, ":");

    } else if (strcmp(keyword, "reduce_fac") == 0) {
        input->reduce_fac = atoi(parm_str);

    } else if (strcmp(keyword, "noext") == 0) {
        input->noext = atoi(parm_str);

    } else if (strcmp(keyword, "merged") == 0) {
        parse_file_name(parm_str, tmp_file);
        strcpy(input->merged, tmp_file);

    } else if (strcmp(keyword, "loneast") == 0) {
        input->loneast = atof(parm_str);

    } else if (strcmp(keyword, "lonwest") == 0) {
        input->lonwest = atof(parm_str);

    } else if (strcmp(keyword, "latnorth") == 0) {
        input->latnorth = atof(parm_str);

    } else if (strcmp(keyword, "latsouth") == 0) {
        input->latsouth = atof(parm_str);

    } else if (strcmp(keyword, "verbose") == 0) {
        input->verbose = atoi(parm_str);

    } else {
        goto Invalid_return;

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
