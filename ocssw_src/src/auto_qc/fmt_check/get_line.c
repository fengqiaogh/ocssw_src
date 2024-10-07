#include <string.h>
#include <stdio.h>

char *get_line(char* str, int nchar, FILE *stream, int enc)
/*******************************************************************

   s_parse

   purpose: similar to fgets, read a line in which does not begin 
            with a special comment character, remove any newline and 
            return the line

   Returns type: char * the pointer to the line storage or NULL
            if an error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char*             str              O      array to fill with line
      int               nchar            I      limit on # char to read
      FILE *            stream           I      input file id
      int               enc              I      comment character such as
                                                # ...

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       14-Mar-1995     Original development

 *******************************************************************/ {
    char *ptr;

    /*
     *  read in lines until the start character is not the comment char
     */
    do {
        if (fgets(str, nchar, stream) == NULL)
            return NULL;
    } while (*str == enc);

    /*
     *  remove any ending newline that can be left by the fgets routine
     */
    if ((ptr = strchr(str, '\n')) != NULL)
        *ptr = 0;

    /*
     *  and exit
     */
    return str;
}
