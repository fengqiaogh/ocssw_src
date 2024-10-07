#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hdf.h"
#include "mfhdf.h"

#define MAXNVDATA 128

int32 wr_vdata(char *outname, int32 fileid_w, int32 vgid, char *name,
        char *class1, int32 n_flds, int32 n_recs_to_write,
        char *fldname[], int32 type[], int32 noext, uint8 *data,
        int32 verbose) {

    int32 i, j;
    int32 status;
    int32 n_write;

    static int32 prod_count = 0;

    static int32 n_vdata_open = 0;
    static int32 vdata_id_w[MAXNVDATA];
    static int32 aid[MAXNVDATA];
    static char *vdata_id_name[MAXNVDATA];
    static char empty[8] = {"_EMPTY_"};

    char buffer[1000];


    /* Closing Vdata */
    /* ------------- */
    /*  if (n_recs_to_write == 0) { */
    if (data == NULL) {
        for (j = 0; j < n_vdata_open; j++) {
            if (strcmp(vdata_id_name[j], name) == 0) {
                printf("Detaching: %s\n", vdata_id_name[j]);
                free(vdata_id_name[j]);
                vdata_id_name[j] = &empty[0];
                VSdetach(vdata_id_w[j]);
                if (aid[j] != -1) Hendaccess(aid[j]);
                return 0;
            }
        }
    }


    /* If existing vdata goto write */
    /* ---------------------------- */
    for (j = 0; j < n_vdata_open; j++) {
        if (strcmp(vdata_id_name[j], name) == 0) {
            goto WRITE_ONLY;
        }
    }


    /* Setup new vdata */
    /* --------------- */
    j = n_vdata_open;
    vdata_id_w[n_vdata_open] = VSattach(fileid_w, -1, "w");
    aid[n_vdata_open] = -1;
    vdata_id_name[n_vdata_open] = (char *) malloc(strlen(name) + 1);
    strcpy(vdata_id_name[n_vdata_open], name);
    printf("Setting Up: %s %d\n", name, vdata_id_w[n_vdata_open]);

    for (i = 0; i < n_flds; i++) {
        status = VSfdefine(vdata_id_w[n_vdata_open], fldname[i], type[i], 1);
        if (status != 0) {
            printf("Error defining \"%s\".\n", fldname[i]);
        }
    }

    /* Set fieldnames */
    /* -------------- */
    strcpy(buffer, fldname[0]);
    for (i = 1; i < n_flds; i++) {
        strcat(buffer, ",");
        strcat(buffer, fldname[i]);
    }
    VSsetfields(vdata_id_w[j], buffer);

    Vinsert(vgid, vdata_id_w[n_vdata_open]);
    VSsetname(vdata_id_w[n_vdata_open], name);
    VSsetclass(vdata_id_w[n_vdata_open], class1);

    VSsetblocksize(vdata_id_w[n_vdata_open], 4096 * 6);


    /* Setup External Data Vdatas */
    /* -------------------------- */
    if (strcmp(class1, "DataSubordinate") == 0 && noext == 0) {

        static FILE *sfile;

        memset(buffer, 0, sizeof (buffer));
        status = VSwrite(vdata_id_w[j], (uint8 *) buffer, 1, FULL_INTERLACE);
        VSdetach(vdata_id_w[j]);
        vdata_id_w[j] = VSattach(fileid_w, VSfind(fileid_w, name), "w");

        /*  synthesize file name  */
        strcpy(buffer, outname);
        sprintf(&buffer[strlen(buffer)], ".x%02d", prod_count);

        /*  open file and write name of main file  **
         **  in first 512 bytes                     */
        sfile = fopen(buffer, "w");

        if (strrchr(outname, '/') == NULL)
            strcpy(buffer, outname);
        else
            strcpy(buffer, strrchr(outname, '/') + 1);
        status = fwrite(buffer, 512, 1, sfile);
        fclose(sfile);

        /*  convert to external element  */
        strcpy(buffer, outname);
        sprintf(&buffer[strlen(buffer)], ".x%02d", prod_count);
        aid[j] = HXcreate(fileid_w, DFTAG_VS, (uint16) VSfind(fileid_w, name),
                buffer, 512, 0);

        VSdetach(vdata_id_w[j]);
        vdata_id_w[j] = VSattach(fileid_w, VSfind(fileid_w, name), "w");
        VSseek(vdata_id_w[j], 0);

        prod_count++;
    }

    n_vdata_open++;


WRITE_ONLY:

    if (n_recs_to_write > 0) {
        n_write = VSwrite(vdata_id_w[j], data, n_recs_to_write, FULL_INTERLACE);

        if (n_write != n_recs_to_write) {
            printf("n_write: %d n_recs_to_write: %d \n", n_write, n_recs_to_write);
            printf("Error writing to vdata: \"%s\" (%d)\n", name, vdata_id_w[j]);
            exit(1);
        } else {
            if (verbose == 1)
                printf("Writing %d records to %s (%d)\n",
                    n_write, name, VSelts(vdata_id_w[j]));
        }
    }

    return 0;
}

