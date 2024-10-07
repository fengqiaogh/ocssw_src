#include <string.h>
#include "fmt_check.h"

char *s_parse(char* str, int enc)
/*******************************************************************

   s_parse

   purpose: along the lines of strtok, s_parse parses a string
            using the space as a delimiter as well as an enclosing
            character like " or ' to return substrings

   Returns type: char * the pointer to the next substring or NULL
            if string is exhausted

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char*             str              I      input string on first call
                                                NULL for subsequent calls
      int               enc              I      enclosing character such as
                                                ' or "

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       10-Mar-1995     Original development

 *******************************************************************/ {
    static char *ptr;
    char *s_ptr;
    int start, ienc;

    /* get the location of this string locally  */
    if (str != NULL)
        ptr = str;
    /*  1 - get to a character which is non blank or enclosing */
    start = 0;
    ienc = 0;
    do {
        if (*ptr == 0)
            return NULL;
        else if (*ptr == enc) {
            ienc = 1;
            s_ptr = ptr + 1;
            start = (*(ptr + 1) == 0) ? 0 : 1;
        } else if (*ptr != ' ') {
            start = 1;
            s_ptr = ptr;
        }
        ptr++;
    } while (start == 0);

    /*  find the next blank or, if a quote the ending quote */

    start = 0;
    do {
        if (*ptr == 0)
            return s_ptr;
        else if (ienc == 1 && *ptr == enc) {
            *ptr++ = 0;
            return s_ptr;
        } else if (ienc == 0 && *ptr == ' ') {
            *ptr++ = 0;
            return s_ptr;
        } else
            ptr++;
    } while (start == 0);

    return NULL;

}








