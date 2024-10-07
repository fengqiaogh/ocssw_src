
/* determine if year is a leap year or not */
int isleap(int year) {
    if (((year % 400) == 0) || (((year % 4) == 0) && ((year % 100) != 0)))
        return 1;
    else
        return 0;
}
