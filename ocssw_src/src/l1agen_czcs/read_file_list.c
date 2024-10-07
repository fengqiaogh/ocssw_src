#include <stdio.h>
#include <string.h>

int read_file_list(char *filename, char files[][FILENAME_MAX], int maxfiles)
/*******************************************************************

   read_file_list

   purpose: read a list of file names (or any strings) from a file

   Returns type: int - count of # files (lines) read

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            filename         I      name of file with list
      char[][]          files            O      returned file names
      int               maxfiles         I      max # files to read

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      B. Franz, SAIC    13 Sep 2002      Original development
      W. Robinson, SAIC  9 Aug 2004      make a stand-alone routine

 *******************************************************************/ {
    FILE *fp = NULL;
    int nfiles = 0;

    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stderr,
                "-W- %s line %d: unable to open %s for reading\n",
                __FILE__, __LINE__, filename);
        return (-1);
    }

    while (nfiles < maxfiles && (fscanf(fp, "%s\n", files[nfiles]) != EOF))
        if (strlen(files[nfiles]) > 1)
            nfiles++;

    fclose(fp);

    return (nfiles);
}
