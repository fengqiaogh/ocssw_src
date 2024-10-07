#include <timeutils.h>

#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[]) {
    // test this leap second record
    // 2015 JUL  1 =JD 2457204.5  TAI-UTC=  36.0       S + (MJD - 41317.) X 0.0      S
   
    double unixTime; // unix time = UTC epoch 1970
    double tai58Time;  // TAI58 time
    double tai93Time;  // TAI93 time
    double leapSeconds;
    
    double leapSecondsAt93 = 27;
    double tai58at93 = 1104537600.0 + leapSecondsAt93;
    double unixAt93 = ymds2unix(1993, 1, 1, 0.0);
    
    unixTime = ymds2unix(1958, 1, 1, 0.0);
    tai93Time = unix_to_tai93(unixTime);
    printf("%s unix=%.1f, TAI93=%.1f\n", unix2isodate(unixTime, 'G'), unixTime, tai93Time);

    unixTime = ymds2unix(1993, 1, 1, 0.0);
    tai93Time = unix_to_tai93(unixTime);
    printf("%s unix=%.1f, TAI93=%.1f\n", unix2isodate(unixTime, 'G'), unixTime, tai93Time);

    unixTime = ymds2unix(2015, 7, 1, 0.0);
    tai93Time = unix_to_tai93(unixTime);
    printf("%s unix=%.1f, TAI93=%.1f\n", unix2isodate(unixTime, 'G'), unixTime, tai93Time);

    unixTime = tai93_to_unix(tai93Time);
    printf("%s unix=%.1f, TAI93=%.1f\n", unix2isodate(unixTime, 'G'), unixTime, tai93Time);
    
    double u93 = ymds2unix(1993, 1, 1, 0.0);
    double u58 = ymds2unix(1958, 1, 1, 0.0);
    printf("unix93 - unix58 = %.1f\n", u93-u58);
   
    
    printf("\n\n");
    unixTime = ymds2unix(2015, 7, 1, 0.0);
    unixTime -= 11;
    
    // loop over the leap second
    for(int i=0; i<22; i++) {
        tai58Time = unix_to_tai58(unixTime);
        tai93Time = unix_to_tai93(unixTime);
        leapSeconds = tai93Time - (unixTime - unixAt93) + leapSecondsAt93;
        printf("%s unix=%.1f, TAI58=%.1f, TAI93=%.1f, leap=%.1f\n", unix2isodate(unixTime, 'G'), 
                unixTime, tai58Time, tai93Time, leapSeconds);
        unixTime += 1;
    }
    

    printf("\n\n");
    unixTime = ymds2unix(2015, 7, 1, 0.0);
    tai93Time = unix_to_tai93(unixTime);
    tai93Time -= 11;

    for(int i=0; i<23; i++) {
        unixTime = tai93_to_unix(tai93Time);
        tai58Time = tai93Time + tai58at93;
        leapSeconds = tai93Time - (unixTime - unixAt93) + leapSecondsAt93;
        printf("%s unix=%.1f, TAI58=%.1f, TAI93=%.1f, leap=%.1f\n", unix2isodate(unixTime, 'G'), 
                unixTime, tai58Time, tai93Time, leapSeconds);
        tai93Time += 1;
    }
    

    // test the values around the leap second 2015-7-1
    int returnVal = EXIT_SUCCESS;
    int numPoints = 4;
    double unix_times[]  = {1435708799, 1435708799, 1435708800, 1435708801};
    double tai58[] = {1814400034, 1814400035, 1814400036, 1814400037};
    double tai93[] = { 709862407,  709862408,  709862409,  709862410};

    // test unix -> tai58
    // and  unix -> tai93
    for(int i=0; i<numPoints; i++) {
        if(i==1) // we need to skip the second point since UTC can not represent it
            i++;
        tai58Time = unix_to_tai58(unix_times[i]);
        if(tai58Time != tai58[i]) {
            printf("unix_to_tai58(%.1f) Failed, %.1f != %.1f\n", unix_times[i], tai58Time, tai58[i]);
            returnVal = EXIT_FAILURE;
        }
        tai93Time = unix_to_tai93(unix_times[i]);
        if(tai93Time != tai93[i]) {
            printf("unix_to_tai93(%.1f) Failed, %.1f != %.1f\n", unix_times[i], tai58Time, tai58[i]);
            returnVal = EXIT_FAILURE;
        }
    }

    // test tai58 -> unix
    // and  tai93 -> unix
    for(int i=0; i<numPoints; i++) {
        unixTime = tai58_to_unix(tai58[i]);
        if(unixTime != unix_times[i]) {
            printf("tai58_to_unix(%.1f) Failed, %.1f != %.1f\n", tai58[i], unixTime, unix_times[i]);
            returnVal = EXIT_FAILURE;
        }
        unixTime = tai93_to_unix(tai93[i]);
        if(unixTime != unix_times[i]) {
            printf("tai93_to_unix(%.1f) Failed, %.1f != %.1f\n", tai93[i], unixTime, unix_times[i]);
            returnVal = EXIT_FAILURE;
        }
    }

    return returnVal;
}

