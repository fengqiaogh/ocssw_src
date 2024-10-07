/***********************************************************************
 *
 *     This subroutine converts uppercase letters in INSTR to lowercase.
 *
 *     B.D.Schieber, GSC, 3/94
 *
 ************************************************************************/
#include <ctype.h>

char *lowcase(char *instr) {
    char *strptr = instr;

    while (*instr != '\0') {
        if (isupper(*instr)) *instr = tolower(*instr);
        instr++;
    }

    return (strptr);
}

