/********************************************************************* 
 * countann
 *
 * countann - count number of annotations in an HDF annotation file
 * Brian Schieber SAIC/GSC 5/93
 *********************************************************************/

#include "ancil.h"
#include "ancnrt_proto.h"

int count_annot(char *filename) {
    int cnt;
    FILE *fp;
    char s[MAXDESCLEN];

    if ((fp = fopen(filename, "r")) == NULL) {
        printf("Error opening %s\n", filename);
        return (-1);
    }

    cnt = 0;
    while (fgets(s, 200, fp) != NULL) {
        if (!strncmp(&s[0], "#", 1)) continue;
        cnt++;
    }

    fclose(fp);
    return (cnt);

} /* count_annot() */
