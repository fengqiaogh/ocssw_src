#include <timeutils.h>

double zulu2unix(char *zulu) {
    int yy = 0, dd = 0;
    int hh = 0, mn = 0, sc = 0, cc = 0, mc = 0;
    int32_t year, day;
    double sec;
    int status = 0;
    int32_t len = 0;
    len = strlen(zulu);

    switch (len) {
    case 7:
        status = sscanf(zulu, "%4d%3d", &yy, &dd);
        break;
    case 9:
        status = sscanf(zulu, "%4d%3d%2d", &yy, &dd, &hh);
        break;
    case 11:
        status = sscanf(zulu, "%4d%3d%2d%2d", &yy, &dd, &hh, &mn);
        break;
    case 13:
        status = sscanf(zulu, "%4d%3d%2d%2d%2d", &yy, &dd, &hh, &mn, &sc);
        break;
    case 15:
        status = sscanf(zulu, "%4d%3d%2d%2d%2d%2d", &yy, &dd, &hh, &mn, &sc, &cc);
        break;
    case 16:
    case 19:
        status = sscanf(zulu, "%4d%3d%2d%2d%2d%3d", &yy, &dd, &hh, &mn, &sc, &mc);
        break;
    default:
        printf("%s : unrecognized zulu time format\n", __FILE__);
        exit(1);
    }
    if (status < 2) {
        printf("%s : problem interpreting zulu time format\n", __FILE__);
        exit(1);
    }
    year = yy;
    day = dd;
    sec = (double) (hh * 3600 + mn * 60 + sc + cc / 100. + mc / 1000.);

    return (yds2unix(year, day, sec));
}

