#include <genutils.h>

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

/** trim white space off of beginning and end of string.
    \param str string to modify
 */
void trimBlanks(char* str) {

    // bail if the string is null
    if (str == NULL)
        return;

    // bail if empty string
    int i = strlen(str) - 1;
    if (i < 0)
        return;

    // find last non-blank char
    while (i >= 0) {
        if (isspace(str[i]))
            i--;
        else
            break;
    }
    int length = i + 1;

    // get rid of beginning spaces
    i = 0;
    while (i < length) {
        if (isspace(str[i]))
            i++;
        else
            break;
    }

    if (i > 0) {
        int j = 0;
        while (i < length) {
            str[j++] = str[i++];
        }
        // NULL terminate the string
        str[j] = 0;

    } else {
        // NULL terminate the string
        str[length] = 0;
    }

}

/** Allocate a string without beginning and ending white space.
    \param str string to use as source (left unmodified)
    \return pointer to newly allocated string without whitespace on ends
 */
char* trimBlanksDup(const char* str) {
    // bail if the string is null
    if (str == NULL)
        return strdup("");

    // bail if empty string
    int i = strlen(str) - 1;
    if (i < 0)
        return strdup("");

    // find last non-blank char
    while (i >= 0) {
        if (isspace(str[i]))
            i--;
        else
            break;
    }
    int length = i + 1;

    // get rid of beginning spaces
    i = 0;
    while (i < length) {
        if (isspace(str[i]))
            i++;
        else
            break;
    }

    char* str2 = (char*) malloc(length - i + 1);

    int j = 0;
    while (i < length) {
        str2[j++] = str[i++];
    }
    // NULL terminate the string
    str2[j] = 0;

    return str2;
}

