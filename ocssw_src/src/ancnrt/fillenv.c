/********************************************************************* 
 * fillenv_annot
 *
 * fillenv_annot - fills array with DAAC-type annotations read from file
 * 
 * Contains sscanf format for reading two double-quoted (") strings on
 * a line.
 * Author: Brian D. Schieber GSC/SAIC 5/93
 *********************************************************************/

#include <mfhdf.h>
#include "ancil.h"
#include "ancnrt_proto.h"

struct annotation *fillenv_annot(char *filename) {

    int cnt;
    FILE *fp;
    char s[MAXDESCLEN];
    char Tlabel[MAXLABLEN];
    char Ttype[MAXLABLEN];
    char Tdescr[MAXDESCLEN];
    struct annotation *tannot;

    if ((tannot = (struct annotation *)
            malloc(sizeof (struct annotation) * 100))
            == NULL) pexit("malloc Annotation");

    if ((fp = fopen(filename, "r")) == NULL) {
        printf("Error opening %s\n", filename);
        exit(1);
    }

    cnt = 0;
    while (fgets(s, 200, fp) != NULL) {
        if (!strncmp(&s[0], "#", 1)) continue;

        /* sscanf(&s[1], "%[^'\"] \" \"  %[^'\"]", Tlabel, Tdescr);  */

        sscanf(&s[1], "%[^'\"] \" \"  %[^'\"] \" \"  %[^'\"]",
                Tlabel, Ttype, Tdescr);

        /* printf("Count [%d]: [%s] [%s] [%s]\n", cnt, Tlabel, Ttype, Tdescr); */

        strcpy(tannot[cnt].label, Tlabel);
        if (!strcmp(Ttype, "DFNT_CHAR")) tannot[cnt].type = DFNT_CHAR;
        else if (!strcmp(Ttype, "DFNT_INT16")) tannot[cnt].type = DFNT_INT16;
        else if (!strcmp(Ttype, "DFNT_INT32")) tannot[cnt].type = DFNT_INT32;
        else if (!strcmp(Ttype, "DFNT_FLOAT32")) tannot[cnt].type = DFNT_FLOAT32;
        else {
            printf("Error finding datatype in fillenv.c\n");
            break;
        }
        strcpy(tannot[cnt].descr, Tdescr);

        strcpy(Tlabel, "");
        strcpy(Ttype, "");
        strcpy(Tdescr, "");

        cnt++;
    }

    if (fclose(fp) != 0) {
        printf("Error closing %s\n", filename);
        exit(1);
    }
    return (tannot);

} /* fillenv_annot() */
