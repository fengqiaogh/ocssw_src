#include <genutils.h>

#include <ctype.h>

/** check if a string is a valid integer
    @param str string to check
    @return 1 if valid int, 0 if not
 */
int isValidInt(const char* str) {

    // Handle leading sign character.
    //
    if (*str == '-')
        str++;
    else if (*str == '+')
        str++;

    // Handle empty string or just "-" or "+"
    //
    if (!*str)
        return 0;

    // Check for non-digit chars in the rest of the stirng.
    //
    while (*str) {
        if (isdigit(*str))
            str++;
        else
            return 0;
    }

    return 1;
}

