#include <stdio.h>
#include <genutils.h>
#include <string.h>
#include <strings.h>  // strcasecmp
#include <stdlib.h>

#define THE_STRING_MAX 128
#define FILE_FORMAT_FILE "$OCDATAROOT/common/file_formats.txt"

static int isFormatsLoaded = 0;
static int numFormats = 0;
static char** formatIndex = NULL;
static char** formatName = NULL;
static char** formatExtension = NULL;

/*
This is a sample format file.  Actually the current one.


# The colon-separated columns are:
# 1) index number
# 2) format_name
# 3) default extension
#
1:HDF4:
2:netCDF4:nc
3:HDF5:h5
4:BIN:bin
5:PNG:png
6:TIFF:tiff
7:PPM:ppm
8:JPG:jpg

 */

/**
 * read the format file into the static format arrays
 */
void readFileFormats() {
    int index;
    char line[THE_STRING_MAX];
    char str1[THE_STRING_MAX];
    char fileName[FILENAME_MAX];
    int count;
    FILE* fp;
    int fileLineNum;

    isFormatsLoaded = 1;

    // open the file
    parse_file_name(FILE_FORMAT_FILE, fileName);
    fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("Error - Could not open file format definition file %s\n", fileName);
        exit(1);
    }

    // count lines in the file
    fileLineNum = 0;
    while (fgets(line, THE_STRING_MAX, fp)) {
        fileLineNum++;
        trimBlanks(line);
        if (line[0] && line[0] != '#') {
            count = sscanf(line, "%d:%s", &index, str1);
            if (count != 2) {
                printf("-E- %s Line %d: error reading file type record on line %d of file %s.\n",
                        __FILE__, __LINE__, fileLineNum, fileName);
                exit(1);
            }
            numFormats++;
        }
    }

    // rewind file
    fseek(fp, 0, SEEK_SET);

    // allocate storage
    formatIndex = (char**) calloc(numFormats, sizeof (char*));
    formatName = (char**) calloc(numFormats, sizeof (char*));
    formatExtension = (char**) calloc(numFormats, sizeof (char*));
    if (!formatIndex || !formatName || !formatExtension) {
        printf("Error - Could not allocate space for the file format definition file %s\n", fileName);
        exit(1);
    }

    fileLineNum = 0;
    index = 0;
    while (fgets(line, THE_STRING_MAX, fp)) {
        fileLineNum++;
        trimBlanks(line);
        if (line[0] && line[0] != '#') {
            char *s1, *s2;
            s1 = line;

            // get file type index
            s2 = strchr(s1, ':');
            if (s2 == NULL) {
                printf("-E- %s Line %d: did not fine first ':' on line %d of file %s.\n",
                        __FILE__, __LINE__, fileLineNum, fileName);
                exit(1);
            }
            *s2 = '\0';
            formatIndex[index] = trimBlanksDup(s1);
            s1 = s2 + 1;

            // get file type name
            s2 = strchr(s1, ':');
            if (s2 == NULL) {
                printf("-E- %s Line %d: did not find second ':' on line %d of file %s.\n",
                        __FILE__, __LINE__, fileLineNum, fileName);
                exit(1);
            }
            *s2 = '\0';
            formatName[index] = trimBlanksDup(s1);
            s1 = s2 + 1;

            // get file type extension
            formatExtension[index] = trimBlanksDup(s1);

            index++;
            if (index > numFormats) {
                printf("-E- %s Line %d: Read more lines the second time through on line %d of file %s.\n",
                        __FILE__, __LINE__, fileLineNum, fileName);
                exit(1);
            }
        }
    }
    fclose(fp);
}

/**
 * get the internal index into the format arrays
 * @param str format string to look up
 * @return index into arrays if found, else -1
 */
int getInternalIndex(const char* str) {
    int i;
    int result = -1;

    // bail if NULL is passed in
    if (str == NULL)
        return (result);

    char* inStr = trimBlanksDup(str);

    // bail if empty string
    if (inStr[0] == 0)
        return (result);

    if (!isFormatsLoaded)
        readFileFormats();

    // loop through the arrays
    for (i = 0; i < numFormats; i++) {

        // look in the index array
        if (strcmp(inStr, formatIndex[i]) == 0) {
            result = i;
            break;
        }

        // look in the format array
        if (strcasecmp(inStr, formatName[i]) == 0) {
            result = i;
            break;
        }

        // look in the extension array
        if (strcasecmp(inStr, formatExtension[i]) == 0) {
            result = i;
            break;
        }

    }

    free(inStr);
    return result;
}

/**
 * convert the input string to the file format index.
 * Input string can be index, name or extension
 * Formats are defined in $OCDATAROOT/common/file_formats.txt
 *
 * @param str input string
 * @return index associated with the file format
 */
int getFileFormatIndex(const char* str) {
    int i = getInternalIndex(str);
    if (i == -1) {
        return -1;
    }
    return atoi(formatIndex[i]);
}

/**
 * convert the input string to the normalized file format name.
 * Input string can be index, name or extension
 * Formats are defined in $OCDATAROOT/common/file_formats.txt
 *
 * @param str input string
 * @return pointer to internal string holding the normalized format
 *         string, or NULL
 */
const char* getFileFormatName(const char* str) {
    int i = getInternalIndex(str);
    if (i == -1) {
        return NULL;
    }
    return formatName[i];
}

/**
 * convert the input string to the file format extension.
 * Input string can be index, name or extension
 * Formats are defined in $OCDATAROOT/common/file_formats.txt
 *
 * @param str input string
 * @return pointer to internal string holding the file format extension
 */
const char* getFileFormatExtension(const char* str) {
    int i = getInternalIndex(str);
    if (i == -1) {
        return NULL;
    }
    return formatExtension[i];
}

