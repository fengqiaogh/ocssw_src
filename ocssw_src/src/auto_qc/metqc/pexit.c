/***************************************************************** 
 *
 * print error message and exit
 *
 * Brian Schieber SAIC/GSC 5/93
 *****************************************************************/
#include <stdio.h>
#include <stdlib.h>

void pexit(char *string) {
    printf("FATAL Error in routine... %s\n", string);
    exit(-1);
}

/***************************************************************** 
 *
 * print warning message and return error code 
 *
 *****************************************************************/

void pwarning(char *string) {
    printf("WARNING: Error in routine... %s\n", string);
}


