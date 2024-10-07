/* ---------------------------------------------------------------------- */
/*                 runl1a version of getstationinfo()                     */
/* ---------------------------------------------------------------------- */

/*
This function fills a structure with information about the local
HRPT station.  The information is read from a disk file on the
local system.

The disk file should be a text file containing lines of the form:

NAME: VALUE

where NAME is one of the strings listed below in the #define statements,
and VALUE is the appropriate string for the local station.  The NAME
field is not case sensitive, but the VALUE field is.  Long VALUE fields
may be continued on subsequent lines by beginning each continuation
line with whitespace.  Only the first colon is significant for separating
the NAME and VALUE fields.  All whitespace is trimmed from the beginning
and end of the NAME and VALUE fields, and adjacent whitespace characters
are collapsed into a single space.

Norman Kuring		28-Oct-1997

Integrated with SEADAS version, 20-May-2003, BAF.
 */

#include <GetStationInfo.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <passthebuck.h>

#define CALLOC(ptr,typ,num) {                                           \
  (ptr) = (typ *)calloc((num) , sizeof(typ));                           \
  if((ptr) == NULL){                                                    \
    fprintf(stderr,"-E- %s line %d: Memory allocation failure.\n",      \
    __FILE__,__LINE__);                                                 \
    return(MEMORY_ALLOCATION_ERROR);                                   \
  }                                                                     \
}

#define ENVIRONMENT_VAR  "HRPT_STATION_IDENTIFICATION_FILE"
#define CODE_STR  "code"
#define DATA_CENTER_STR  "data center"
#define STATION_NAME_STR "station name"
#define STATION_LAT_STR  "station latitude"
#define STATION_LON_STR  "station longitude"

int GetStationInfo(char *stationInfoFile, StationInfo *stationInfo) {

    char *bytes, *string, *bp, *next, **sd, *sp;
    int gotcolon, word;
    struct stat stat_struct;
    FILE *fh;
    char *code, *data_center, *station_name;
    char *station_latitude, *station_longitude;

    code = data_center = station_name = station_latitude = station_longitude = NULL;

    /* Get the station ID filename */
    if (stationInfoFile == NULL) {
        stationInfoFile = getenv(ENVIRONMENT_VAR);
        if (stationInfoFile == NULL) {
            fprintf(stderr, "-W- %s line %d: ", __FILE__, __LINE__);
            fprintf(stderr, "Environment variable, \"%s\", not set.\n", ENVIRONMENT_VAR);
            return (MISSING_ENVVAR);
        }
    }
    printf("Station information derived from %s\n", stationInfoFile);

    /* Get the file size and make sure it's greater than zero. */
    if (stat(stationInfoFile, &stat_struct) != 0) {
        fprintf(stderr, "-W- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "stat() failed for file, %s.  ", stationInfoFile);
        perror("Reason");
        return (STAT_FAILURE);
    }
    if (stat_struct.st_size < 1) {
        fprintf(stderr, "-W- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "File, %s, is empty.\n", stationInfoFile);
        return (FILE_IS_EMPTY);
    }

    /* Allocate memory to hold the file contents. */
    CALLOC(bytes, char, stat_struct.st_size + 1);
    CALLOC(string, char, stat_struct.st_size + 1);

    /* Open the file and read its contents. */
    if ((fh = fopen(stationInfoFile, "r")) == NULL) {
        fprintf(stderr, "-W- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "Could not open file, %s .  ", stationInfoFile);
        perror("Reason");
        return (FOPEN_FAILURE);
    }
    if (fread(bytes, 1, stat_struct.st_size, fh) != stat_struct.st_size) {
        fprintf(stderr, "-W- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "Error reading file, %s .  ", stationInfoFile);
        perror("Reason");
        return (FREAD_FAILURE);
    }
    fclose(fh);

    /* Merge continuation lines (i.e. lines that begin with whitespace). */
    bp = bytes;
    next = bp + 1;
    while (*next) {
        if (*bp == '\n' && isspace(*next)) *bp = ' ';
        bp++;
        next++;
    }

    /* Parse the file contents for name=value pairs and fill a structure
     *  with the desirable ones (see the #define lines above). */
    bp = bytes;
    sp = string;
    gotcolon = 0;
    word = 0;
    while (*bp) {
        while (isspace(*bp) && *bp != '\n') bp++;
        if (*bp == '\n') {
            gotcolon = 0;
            *sp = 0;
            if (sd != NULL) {
                if (*string)
                    sp--;
                *sp = 0;
                *sd = strdup(string);
                sd = NULL;
            }
            sp = string;
            bp++;
        }
        while (*bp && !isspace(*bp) && (gotcolon || *bp != ':')) {
            word = 1;
            if (gotcolon)
                *sp++ = *bp;
            else
                *sp++ = tolower(*bp);
            bp++;
        }
        if (word) {
            *sp = ' ';
            sp++;
            word = 0;
        }
        if (*bp == ':' && !gotcolon) {
            gotcolon = 1;
            bp++;
            sp--;
            *sp = 0;
            if (strcmp(string, CODE_STR) == 0) {
                sd = &code; /* sd = string destination */
            } else if (strcmp(string, DATA_CENTER_STR) == 0) {
                sd = &data_center;
            } else if (strcmp(string, STATION_NAME_STR) == 0) {
                sd = &station_name;
            } else if (strcmp(string, STATION_LAT_STR) == 0) {
                sd = &station_latitude;
            } else if (strcmp(string, STATION_LON_STR) == 0) {
                sd = &station_longitude;
            } else {
                sd = NULL;
            }
            sp = string;
        }
    }
    /*
    Account for the case of a line that is terminated by the
    end of the file rather than by a newline (0x0a) character.
     */
    if (sd != NULL && *string) {
        *sd = strdup(string);
    }

    if (code != NULL) {
        /* Copy at most 4 chars into stationInfo->code then terminate it */
        strncpy(stationInfo->code, code, 4);
        stationInfo->code[4] = (char) 0;

        /***
         **** This does not work with the 4-char code logic in swl1_hdf.c:DataTypeString():1346.
         **** JGW 2/24/2006
         ****
        if(code != NULL && strlen(code) >= 3){
          if(isalpha(code[0]) && isalpha(code[1]) && isalpha(code[2])){
            code[3] = 0;
            strcpy(stationInfo->code, code);
          }
         ****
         ****/
        free(code);
        code = NULL;
    }
    if (data_center != NULL) {
        stationInfo->data_center = data_center;
    }
    if (station_name != NULL) {
        stationInfo->station_name = station_name;
    }
    if (station_latitude != NULL) {
        stationInfo->station_latitude = (float) atof(station_latitude);
        free(station_latitude);
        station_latitude = NULL;
    }
    if (station_longitude != NULL) {
        stationInfo->station_longitude = (float) atof(station_longitude);
        free(station_longitude);
        station_longitude = NULL;
    }

    return (LIFE_IS_GOOD);
}



