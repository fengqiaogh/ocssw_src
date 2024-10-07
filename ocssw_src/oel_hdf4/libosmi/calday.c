static int days[] = {-1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

int daysInYear(int year) {
    if ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0)
        return 366;
    else return 365;
}

int daysInMonth(int year, int month) {
    if ((month == 2) && (daysInYear(year) == 366)) return 29;
    else return days[month];
}

int getCalDay(int year, int jday, int *month, int *day) {
    int i, sub_days, mdays;
    sub_days = jday;
    for (i = 1; i < 13; i++) {
        mdays = daysInMonth(year, i);
        sub_days -= mdays;
        if (sub_days <= 0)
            break;
    }
    if (sub_days == 0) {
        *month = i;
        *day = daysInMonth(year, i);
    } else {
        *month = i;
        *day = daysInMonth(year, i) + sub_days;
    }
    return 0;
}

int getJulianDay(int year, int month, int day) {
    int i;
    int julday = 0;

    if (month < 1 || month > 12)
        return -1;
    if (day < 1 || day > daysInMonth(year, month))
        return -1;
    for (i = 1; i < month; i++)
        julday += daysInMonth(year, i);
    julday += day;
    return julday;
}
