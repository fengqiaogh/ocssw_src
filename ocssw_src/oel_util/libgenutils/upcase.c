/***********************************************************************
 *
 *     This subroutine converts lowercase letters in INSTR to uppercase.
 *
 *     B.D.Schieber, GSC, 2/93
 *
 ************************************************************************/
#include <ctype.h>

char *upcase(char *instr) {
    char *strptr = instr;

    while (*instr != '\0') {
        if (islower(*instr)) *instr = toupper(*instr);
        instr++;
    }

    return (strptr);
}

