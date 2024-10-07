#include <clo.h>

#include <string.h>             /* strdup,strcmp */
#include <strings.h>            /* strcasecmp */
#include <stdlib.h>             /* malloc,free,getenv */
#include <assert.h>             /* assert */
#include <ctype.h>              /* tolower */
#include <stdio.h>              /* printf */
#include <errno.h>              /* errno */

#include <genutils.h>

// how many spaces needed to indent each XML nesting level
#define XML_INDENT_STRING "    "

static char *clo_versionString = NULL;
static char *clo_helpString = NULL;
static int clo_ignoreKeyCase = 1;
static int clo_enableDumpOptions = 1;
static clo_programMetadataList_t clo_programMetadataList;
static char **selectOptionKeys = NULL;

// should reading the command line or par file add options to the list
static int clo_enableExtraOptions = 0;

// are command line position options allowed
static int clo_enablePositionOptions = 0;

// is par file processing enabled
static int clo_enableParOption = 1;

/** CLO callback function for the "-help" option. */
void clo_helpOptionCb(struct clo_option_t *option) {
    // check for version
    if (clo_getOptionBool(option)) {
        clo_printUsage((clo_optionList_t*) option->cb_data);
        exit(0);
    }
}

/** CLO callback function for the "-version" option. */
void clo_versionOptionCb(struct clo_option_t *option) {
    // check for help
    if (clo_getOptionBool(option)) {
        clo_printVersion();
        exit(0);
    }
}

/**
 * Callback function for "par" option
 *
 */
void clo_parOptionCb(clo_option_t *option) {
    if (option->cb_flag)
        clo_readFile((clo_optionList_t*) option->cb_data, option->valStr);
}

/**
 * function to manage an array that grows as you add elements
 *
 * @param array pointer to the array
 * @param storageSize current internal storage size of the array
 * @param count current size of the array
 * @param ptr pointer to add to the array
 */
void clo_addToArray(void ***array, int *storageSize, int *count, void *ptr) {
    assert(array);
    assert(storageSize);
    assert(count);
    assert(ptr);

    if (*count >= *storageSize) {
        int i = 0;
        void **oldArray = *array;
        *storageSize += CLO_CHUNK_SIZE;
        *array = (void**) malloc(sizeof (void*) * (*storageSize));
        assert(*array);
        if (oldArray) {
            for (; i < *count; i++)
                (*array)[i] = oldArray[i];
            free(oldArray);
        }
        for (; i < *storageSize; i++)
            (*array)[i] = NULL;
    }
    (*array)[*count] = ptr;
    (*count)++;
}

/**
 * Expand any environment variables referenced in the string
 * Here are some examples:
 *
 * $HOME/otherDir
 * text$(HOME)moreText
 * text${HOME}moreText
 * $HOME/more$(HOME)thanOne
 *
 * @param str string to substitute environment variables into
 * @return pointer to a newly allocated string that can to be freed with free(3)
 */
char* clo_envExpandString(const char* str) {
    char outStr[SILLYSTRING];
    char envStr[512];
    const char* inPtr;
    char* outPtr;
    char* envPtr;
    const char* envStart;
    char* envVal;
    int i;

    inPtr = str;
    outPtr = outStr;
    while (*inPtr) {

        if (*inPtr == '$') { /* found an environment variable */
            envStart = inPtr;
            inPtr++;
            if (*inPtr == '(' || *inPtr == '{')
                inPtr++;
            envPtr = envStr;
            while (isalpha(*inPtr) || isdigit(*inPtr) || *inPtr == '_') {
                *envPtr = *inPtr;
                inPtr++;
                envPtr++;
            }
            if (*inPtr == ')' || *inPtr == '}')
                inPtr++;
            *envPtr = '\0';
            envPtr++;

            envVal = getenv(envStr);
            if (envVal) {
                i = strlen(envVal);
                memcpy(outPtr, envVal, i);
                outPtr += i;
            } else {
                i = inPtr - envStart;
                memcpy(outPtr, envStart, i);
                outPtr += i;
            }

        } else { /* copy char into outStr */
            *outPtr = *inPtr;
            inPtr++;
            outPtr++;
        }
    } // while ptr1

    *outPtr = '\0';
    return strdup(outStr);
}

/**
 * Return the string representation of the passed data type
 *
 * @param dataType type to decode
 * @return pointer a static string representation of the type
 */
char* clo_dataTypeToString(enum clo_dataType_t dataType) {
    char* dataTypeStr;
    switch (dataType) {
    case CLO_TYPE_BOOL:
        dataTypeStr = "boolean";
        break;
    case CLO_TYPE_INT:
        dataTypeStr = "int";
        break;
    case CLO_TYPE_INT64:
        dataTypeStr = "int64";
        break;
    case CLO_TYPE_FLOAT:
        dataTypeStr = "float";
        break;
    case CLO_TYPE_DOUBLE:
        dataTypeStr = "double";
        break;
    case CLO_TYPE_STRING:
        dataTypeStr = "string";
        break;
    case CLO_TYPE_IFILE:
        dataTypeStr = "ifile";
        break;
    case CLO_TYPE_OFILE:
        dataTypeStr = "ofile";
        break;
    case CLO_TYPE_HELP:
        dataTypeStr = "helpString";
        break;
    case CLO_TYPE_POSITION:
        dataTypeStr = "position";
        break;
    default:
        dataTypeStr = "unknownType";
        break;
    }
    return dataTypeStr;
}

/**
 * divide strVal (or defaultVal) into an array of strings and set count
 *
 * @param option structure to operate on
 */
void clo_parseOptionString(clo_option_t *option) {
    assert(option);
    int i;
    char *str;
    char *ptr;

    // delete old string array
    if (option->strArray) {
        for (i = 0; i < option->count; i++) {
            if (option->strArray[i])
                free(option->strArray[i]);
        }
        free(option->strArray);
        option->strArray = NULL;
    }
    option->count = 0;
    int storage = 0;

    // parse the new one, using the valStr or the defaultStr
    if (option->valStr)
        str = strdup(option->valStr);
    else if (option->defaultVal)
        str = strdup(option->defaultVal);
    else {
        return;
    }

    ptr = strtok(str, CLO_ARRAY_DELIMITER);
    while (ptr) {
        clo_addToArray((void***) &(option->strArray), &storage,
                &(option->count), (void*) strdup(ptr));

        ptr = strtok(NULL, CLO_ARRAY_DELIMITER);
    }

    free(str);
}

/**
 * Convert strVal into an array of strings into strArray, then fill up the
 * valArray based on data type.
 *
 * @param option structure to work on
 */
void clo_parseOption(clo_option_t *option) {
    assert(option);

    if (option->valArray) {
        free(option->valArray);
        option->valArray = NULL;
    }
    clo_parseOptionString(option);

    if (option->count == 0) {
        return;
    }

    switch (option->dataType) {
    case CLO_TYPE_BOOL:
    {
        int i;
        int *data;
        char *str;
        char *ptr;
        data = (int*) malloc(sizeof (int) * option->count);
        option->valArray = (void*) data;
        for (i = 0; i < option->count; i++) {
            ptr = str = strdup(option->strArray[i]);

            // make the string all lower case
            while (*ptr) {
                *ptr = tolower(*ptr);
                ptr++;
            }
            if (strcmp(str, "0") == 0 || strcmp(str, "f") == 0
                    || strcmp(str, "false") == 0 || strcmp(str, "n") == 0
                    || strcmp(str, "no") == 0 || strcmp(str, "off") == 0) {
                data[i] = 0;
            } else if (strcmp(str, "1") == 0 || strcmp(str, "t") == 0
                    || strcmp(str, "true") == 0 || strcmp(str, "y") == 0
                    || strcmp(str, "yes") == 0 || strcmp(str, "on") == 0) {
                data[i] = 1;
            } else {
                fprintf(stderr,
                        "-E- clo_parseOption: Invalid boolean value=%s for option key=%s\n",
                        option->strArray[i], option->key);
                exit(1);
            }
            free(str);
        }
    }
        break;
    case CLO_TYPE_INT:
    {
        int i;
        char *ptr;
        int *data = (int*) malloc(sizeof (int) * option->count);
        option->valArray = (void*) data;
        for (i = 0; i < option->count; i++) {
            errno = 0;
            data[i] = strtol(option->strArray[i], &ptr, 0);
            if ((errno != 0) || (ptr == option->strArray[i])
                    || (*ptr != '\0')) {
                fprintf(stderr,
                        "-E- clo_parseOption: Invalid integer value=%s for option key=%s\n",
                        option->strArray[i], option->key);
                exit(1);
            }
        }
    }
        break;
    case CLO_TYPE_INT64:
    {
        int i;
        char *ptr;
        int64_t *data = (int64_t*) malloc(sizeof (int64_t) * option->count);
        option->valArray = (void*) data;
        for (i = 0; i < option->count; i++) {
            errno = 0;
            data[i] = strtol(option->strArray[i], &ptr, 0);
            if ((errno != 0) || (ptr == option->strArray[i])
                    || (*ptr != '\0')) {
                fprintf(stderr,
                        "-E- clo_parseOption: Invalid integer value=%s for option key=%s\n",
                        option->strArray[i], option->key);
                exit(1);
            }
        }
    }
        break;
    case CLO_TYPE_FLOAT:
    {
        int i;
        char *ptr;
        float *data = (float*) malloc(sizeof (float) * option->count);
        option->valArray = (void*) data;
        for (i = 0; i < option->count; i++) {
            errno = 0;
            data[i] = strtod(option->strArray[i], &ptr);
            if ((errno != 0) || (ptr == option->strArray[i])
                    || (*ptr != '\0')) {
                fprintf(stderr,
                        "-E- clo_parseOption: Invalid float value=%s for option key=%s\n",
                        option->strArray[i], option->key);
                exit(1);
            }
        }
    }
        break;
    case CLO_TYPE_DOUBLE:
    {
        int i;
        char *ptr;
        double *data = (double*) malloc(sizeof (double) * option->count);
        option->valArray = (void*) data;
        for (i = 0; i < option->count; i++) {
            errno = 0;
            data[i] = strtod(option->strArray[i], &ptr);
            if ((errno != 0) || (ptr == option->strArray[i])
                    || (*ptr != '\0')) {
                fprintf(stderr,
                        "-E- clo_parseOption: Invalid double value=%s for option key=%s\n",
                        option->strArray[i], option->key);
                exit(1);
            }
        }
    }
        break;
    case CLO_TYPE_STRING: /* do nothing */
    case CLO_TYPE_IFILE: /* do nothing */
    case CLO_TYPE_OFILE: /* do nothing */
    case CLO_TYPE_HELP: /* do nothing */
    case CLO_TYPE_POSITION: /* do nothing */
    default:
        break;
    } // switch dataType
}

/**
 * Tell the library to ignore the case of the option key or not
 *
 * @param val 1 = case insensitive, 0 case sensitive
 */
void clo_setIgnoreKeyCase(int val) {
    if (val)
        clo_ignoreKeyCase = 1;
    else
        clo_ignoreKeyCase = 0;
}

/**
 * get the state of option key case sensitivity
 *
 * @return 1 = case insensitive, 0 case sensitive
 */
int clo_getIgnoreKeyCase() {
    return clo_ignoreKeyCase;
}

/**
 * set whether or not the library looks for the "-dump_options" key
 *
 * @param val 1 = enable dump, 0 = disable
 */
void clo_setEnableDumpOptions(int val) {
    clo_enableDumpOptions = val;
}

/**
 * is the library currently looking for the "-dump_options" key
 *
 * @return 1 = enable dump, 0 = disable
 */
int clo_getEnableDumpOptions() {
    return clo_enableDumpOptions;
}

/**
 * set whether or not the library will automatically add options to
 * the list when reading the command line or par file.
 *
 * @param val 1 = enable adding, 0 = disable
 */
void clo_setEnableExtraOptions(int val) {
    clo_enableExtraOptions = val;
}

/**
 * is the library set to automatically add options to
 * the list when reading the command line or par file.
 *
 * @return 1 = enable adding, 0 = disable
 */
int clo_getEnableExtraOptions() {
    return clo_enableExtraOptions;
}

/**
 * set the version string for the program
 *
 * @param str string to use for the version string
 */
void clo_setVersion(const char* str) {
    if (clo_versionString)
        free(clo_versionString);
    if (str)
        clo_versionString = strdup(str);
    else
        clo_versionString = NULL;
}

/**
 * Set the version string using the program name and version string.
 * The string is set to: "progName versionStr (__DATE__ __TIME__)"
 *
 * @param programName name of the program
 * @param versionStr string to use for the version number
 */
void clo_setVersion2(const char* programName, const char* versionStr) {
    int size = 30 + strlen(programName) + strlen(versionStr);
    char* str = (char*) malloc(size);
    sprintf(str, "%s %s (%s %s)", programName, versionStr, __DATE__, __TIME__);
    if (clo_versionString)
        free(clo_versionString);
    clo_versionString = str;
}

/**
 * get the version string for the program
 *
 * @return the version string
 */
char* clo_getVersion() {
    return clo_versionString;
}

/**
 * set the help string for the program
 *
 * @param str string to use for help
 */
void clo_setHelpStr(const char *str) {
    if (clo_helpString)
        free(clo_helpString);
    if (str)
        clo_helpString = strdup(str);
    else
        clo_helpString = NULL;
}

/**
 * get the help string for the program
 *
 * @return the help string
 */
char* clo_getHelpStr() {
    return clo_helpString;
}

/**
 * set the list of option keys that will be printed by help and XML
 * output.  The last pointer in the array must be NULL.  If keys is
 * NULL then all options will be used.  The strings are not copied,
 * so don't delete the array or the strings.
 *
 * @param keys NULL terminated array of strings containing the
 *             key names
 */
void clo_setSelectOptionKeys(char **keys) {
    selectOptionKeys = keys;
}

/**
 * get the array of select option keys
 *
 * @return NULL terminated array of key names
 */
char **clo_getSelectOptionKeys() {
    return selectOptionKeys;
}

/**
 * create an empty option list container
 *
 * @return pointer to newly allocated and initialized structure
 */
clo_optionList_t* clo_createList() {
    clo_optionList_t* list;
    clo_option_t* option;

    list = (clo_optionList_t*) malloc(sizeof (clo_optionList_t));
    list->storageSize = 0;
    list->numOptions = 0;
    list->options = NULL;

    list->positionStorageSize = 0;
    list->positionNumOptions = 0;
    list->positionOptions = NULL;

    option = clo_addOption(list, "help", CLO_TYPE_BOOL, NULL,
            "print usage information");
    clo_addOptionAlias(option, "h");
    option->cb_data = list;
    option->cb = clo_helpOptionCb;

    option = clo_addOption(list, "version", CLO_TYPE_BOOL, NULL,
            "print the version\n        information");
    option->cb_data = list;
    option->cb = clo_versionOptionCb;

    option = clo_addOption(list, "dump_options", CLO_TYPE_BOOL, NULL,
            "print\n        information about each option");

    option = clo_addOption(list, "dump_options_paramfile", CLO_TYPE_OFILE,
            NULL, "print\n        information about each option to paramfile");

    option = clo_addOption(list, "dump_options_xmlfile", CLO_TYPE_OFILE, NULL,
            "print\n        information about each option to XML file");

    if(clo_enableParOption) {
        option = clo_addOption(list, "par", CLO_TYPE_IFILE, NULL,
                "input parameter file");
        option->cb_data = list;
        option->cb = clo_parOptionCb;
        option->cb_flag = 1;
    }
        
    return list;
}

/**
 * allocate and init an option structure
 *
 * @param key name of the option
 * @param dataType data storage type for the option
 * @param defaultVal initial value
 * @param desc long description of this option
 * @return newly allocated and initialized option structure
 */
clo_option_t* clo_createOption(const char *key, enum clo_dataType_t dataType,
        const char *defaultVal, const char *desc) {
    clo_option_t* opt;

    assert(key);

    opt = (clo_option_t*) malloc(sizeof (clo_option_t));
    opt->key = clo_trimDashesDup(key);
    opt->dataType = dataType;
    if (defaultVal)
        opt->defaultVal = strdup(defaultVal);
    else
        opt->defaultVal = NULL;
    if (desc)
        opt->desc = strdup(desc);
    else
        opt->desc = NULL;
    opt->valStr = NULL;
    opt->source = strdup("default");
    opt->cb = NULL;
    opt->cb_data = NULL;
    opt->strArray = NULL;
    opt->valArray = NULL;
    opt->count = 0;

    if (dataType == CLO_TYPE_BOOL && defaultVal == NULL) {
        opt->defaultVal = strdup("false");
    }

    opt->aliases = NULL;
    opt->numAliases = 0;
    opt->aliasStorageSize = 0;

    opt->position = -1;

    if (defaultVal)
        clo_parseOption(opt);

    return opt;
}

/**
 * add and alias for key to the option
 *
 * @param option option to add the alias to
 * @param alias new key alias for the option
 */
void clo_addOptionAlias(clo_option_t *option, const char* alias) {
    assert(option);
    assert(alias);
    clo_addToArray((void***) &(option->aliases), &(option->aliasStorageSize),
            &(option->numAliases), (void*) clo_trimDashesDup(alias));
}

/**
 * add and alias for key to the option
 *
 * @param list to add the alias to
 * @param key in the list for the getting the option
 * @param alias new key alias for the option
 */
void clo_addAlias(clo_optionList_t *list, const char *key, const char* alias) {
    assert(list);
    assert(key);
    assert(alias);
    clo_option_t *option = clo_findOption(list, key);
    if (!option) {
        fprintf(stderr,
                "-E- clo_addAlias: key=%s does not exist.\n", key);
        exit(1);
    }
    clo_addToArray((void***) &(option->aliases), &(option->aliasStorageSize),
            &(option->numAliases), (void*) clo_trimDashesDup(alias));
}

/**
 * add the given option to the option list
 *
 * @param list option list to add the option into
 * @param option pointer to option structure to add
 */
void clo_insertOption(clo_optionList_t *list, clo_option_t *option) {
    assert(list);
    assert(option);

    clo_addToArray((void***) &(list->options), &(list->storageSize),
            &(list->numOptions), (void*) option);
}

/**
 * create a new option and add it to the list
 *
 * @param list list to add the option to
 * @param key name of the new option
 * @param dataType data type for this option
 * @param defaultVal default value
 * @param desc long description of the option
 * @return the newly created option
 */
clo_option_t* clo_addOption(clo_optionList_t *list, const char *key,
        enum clo_dataType_t dataType, const char *defaultVal,
        const char *desc) {
    clo_option_t* option;
    assert(list);
    assert(key);

    option = clo_findOption(list, key);
    if (option) {
        fprintf(stderr,
                "-E- clo_addOption: option with key=%s already exists.\n", key);
        exit(1);
    }
    option = clo_createOption(key, dataType, defaultVal, desc);
    clo_insertOption(list, option);
    return option;
}

/**
 * create a new option and add it to the position list
 *
 * @param list list to add the option to
 * @return the newly created option
 */
clo_option_t* clo_addPositionOption(clo_optionList_t *list) {
    clo_option_t* option;
    assert(list);

    option = clo_createOption("pos", CLO_TYPE_POSITION, NULL, "position parameter");
    option->position = list->positionNumOptions;
    clo_addToArray((void***) &(list->positionOptions), &(list->positionStorageSize),
            &(list->positionNumOptions), (void*) option);
    return option;
}

/**
 * delete all options in the position list
 *
 * @param list list to add the option to
 */
void clo_clearPositionOptions(clo_optionList_t *list) {
    assert(list);

    int i;
    for (i = 0; i < list->positionNumOptions; i++) {
        if (list->positionOptions[i])
            clo_deleteOption(list->positionOptions[i]);
    }
    list->positionNumOptions = 0;
}

/**
 * allocate a new option and make a deep copy of the given option
 *
 * @param option option to copy
 * @return newly allocated and initialized option
 */
clo_option_t* clo_copyOption(clo_option_t *option) {
    int i;
    clo_option_t* opt;

    assert(option);
    if (!option->key) {
        fprintf(stderr,
                "-E- clo_copyOption: can not copy an option with a NULL key\n");
        exit(1);
    }

    opt = (clo_option_t*) malloc(sizeof (clo_option_t));
    opt->key = strdup(option->key);
    opt->dataType = option->dataType;
    if (option->defaultVal)
        opt->defaultVal = strdup(option->defaultVal);
    else
        opt->defaultVal = NULL;
    if (option->desc)
        opt->desc = strdup(option->desc);
    else
        opt->desc = NULL;
    if (option->valStr)
        opt->valStr = strdup(option->valStr);
    else
        opt->valStr = NULL;
    if (option->source)
        opt->source = strdup(option->source);
    else
        opt->source = NULL;
    opt->cb = option->cb;
    opt->cb_data = option->cb_data;

    opt->aliases = NULL;
    opt->numAliases = 0;
    opt->aliasStorageSize = 0;

    for (i = 0; i < option->numAliases; i++)
        clo_addToArray((void***) &(opt->aliases), &(opt->aliasStorageSize),
            &(opt->numAliases), (void*) strdup(option->aliases[i]));

    // go ahead and make the code re-parse the option for now, so
    // I don't have to figure out how to copy the rest of the data
    opt->strArray = NULL;
    opt->valArray = NULL;
    opt->count = 0;
    opt->position = option->position;
    clo_parseOption(opt);

    return opt;
}

/**
 * allocate a new list and make a deep copy of all of the options
 *
 * @param list list to copy
 * @return newly allocated and copied list
 */
clo_optionList_t* clo_copyList(clo_optionList_t *list) {
    int i;
    clo_option_t *option;

    assert(list);
    clo_optionList_t* newList = clo_createList();
    for (i = 0; i < list->numOptions; i++) {
        option = clo_copyOption(list->options[i]);
        clo_insertOption(newList, option);
    }
    for (i = 0; i < list->positionNumOptions; i++) {
        option = clo_copyOption(list->positionOptions[i]);
        clo_addToArray((void***) &(list->positionOptions), &(list->positionStorageSize),
                &(list->positionNumOptions), (void*) option);
    }

    return newList;
}

/**
 * delete the option and all of its internal storage
 *
 * @param option option to delete
 */
void clo_deleteOption(clo_option_t *option) {
    int i;
    assert(option);
    if (option->key) {
        free(option->key);
        /* set to NULL in case someone tries to access a deleted option */
        option->key = NULL;
    }
    if (option->defaultVal) {
        free(option->defaultVal);
        option->defaultVal = NULL;
    }
    if (option->desc) {
        free(option->desc);
        option->desc = NULL;
    }
    if (option->valStr) {
        free(option->valStr);
        option->valStr = NULL;
    }
    if (option->source) {
        free(option->source);
        option->source = NULL;
    }
    if (option->strArray) {
        for (i = 0; i < option->count; i++) {
            if (option->strArray[i])
                free(option->strArray[i]);
        }
        free(option->strArray);
        option->strArray = NULL;
    } // if strArray
    option->count = 0;
    if (option->valArray) {
        free(option->valArray);
        option->valArray = NULL;
    }
    for (i = 0; i < option->numAliases; i++)
        free(option->aliases[i]);
    free(option->aliases);
    option->aliases = NULL;
    option->numAliases = 0;
    option->aliasStorageSize = 0;

    free(option);
}

/**
 * delete the list and all of its internal options
 *
 * @param list list to delete
 */
void clo_deleteList(clo_optionList_t *list) {
    int i;
    assert(list);
    if (list->options) {
        for (i = 0; i < list->numOptions; i++) {
            if (list->options[i])
                clo_deleteOption(list->options[i]);
        }
        free(list->options);
        list->options = NULL;
        list->numOptions = 0;
    }

    if (list->positionOptions) {
        for (i = 0; i < list->positionNumOptions; i++) {
            if (list->positionOptions[i])
                clo_deleteOption(list->positionOptions[i]);
        }
        free(list->positionOptions);
        list->positionOptions = NULL;
        list->positionNumOptions = 0;
    }

    free(list);
}

/**
 * get the option at index i
 *
 * @param list list to search through
 * @param i index into the list
 * @return option pointer at i, or NULL if not found
 */
clo_option_t* clo_getOption(clo_optionList_t *list, int i) {
    assert(list);
    if (i < 0 || i >= list->numOptions)
        return NULL;
    return list->options[i];
}

/** Remove the dashes from beginning of the string.
    \param str string to remove the dashes from
 */
void clo_trimDashes(char* str) {

    // bail if the string is null
    if (str == NULL)
        return;

    // first trim the spaces
    trimBlanks(str);

    // bail if empty string
    int length = strlen(str);
    if (length <= 0)
        return;

    // get rid of leading dashes
    int i = 0;
    while (i < length && (str[i] == '-' || isspace(str[i])))
        i++;

    // bail if there were no dashes
    if (i == 0)
        return;

    int j = 0;
    while (i < length) {
        str[j++] = str[i++];
    }

    // NULL terminate the string
    str[j] = 0;
}

/** Allocate a string without beginning dashes.
    \param str string to use as source (left unmodified)
    \return pointer to newly allocated string without whitespace on ends or dashes
 */
char* clo_trimDashesDup(const char* str) {
    char* str2 = trimBlanksDup(str);
    clo_trimDashes(str2);
    return str2;
}

/**
 * search list for the option with the given key
 *
 * @param list list to search
 * @param key key to look for
 * @return option pointer or NULL if not found
 */
clo_option_t* clo_findOption(clo_optionList_t *list, const char *key) {
    assert(list);
    assert(key);

    int j;
    int i;
    clo_option_t *option;
    char *keyStr;
    clo_option_t* result = NULL;

    keyStr = clo_trimDashesDup(key);
    // search the defined options
    i = 0;
    while (i < list->numOptions) {
        option = list->options[i];
        if (clo_ignoreKeyCase) {
            if (strcasecmp(option->key, keyStr) == 0) {
                result = option;
                break;
            }
            for (j = 0; j < option->numAliases; j++)
                if (strcasecmp(option->aliases[j], keyStr) == 0) {
                    result = option;
                    break;
                }
        } else {
            if (strcmp(option->key, keyStr) == 0) {
                result = option;
                break;
            }
            for (j = 0; j < option->numAliases; j++)
                if (strcmp(option->aliases[j], keyStr) == 0) {
                    result = option;
                    break;
                }
        }

        i++;
    }

    free(keyStr);

    return result;
}

/**
 * get the number of options in the list
 * @param list list to get the count for
 * @return number of options
 */
int clo_getNumOptions(clo_optionList_t *list) {
    assert(list);
    return list->numOptions;
}

/**
 * return the option's raw string. If the option has not been set the
 * default value is returned.  Exit is called if the option has not
 * been set and the default value is NULL.
 *
 * @param option option to use
 * @return raw string value for the option
 */
char* clo_getOptionRawString(clo_option_t *option) {

    if (option->valStr)
        return option->valStr;

    if (option->defaultVal)
        return option->defaultVal;

    fprintf(stderr, "-E- clo_getOptionRawString: option=%s needs to be set\n",
            option->key);
    exit(1);
}

/**
 * return the first string value from option. If the option has not been
 * set the default value is returned.  Option must have only 1 value.
 *
 * @param option option to use
 * @return first value of the option
 */
char* clo_getOptionString(clo_option_t *option) {
    if (option->dataType != CLO_TYPE_STRING
            && option->dataType != CLO_TYPE_IFILE
            && option->dataType != CLO_TYPE_OFILE
            && option->dataType != CLO_TYPE_HELP) {
        fprintf(stderr,
                "-E- clo_getOptionString: option=%s is not a string option\n",
                option->key);
        exit(1);
    }
    if (option->count == 0) {
        if (option->valStr) {
            return option->valStr;
        } else {
            fprintf(stderr,
                    "-E- clo_getOptionString: option=%s needs to be set\n",
                    option->key);
            exit(1);
        }
    }

    if (option->count != 1) {
        fprintf(stderr,
                "-E- clo_getOptionString: option=%s needs to be set with only one value\n",
                option->key);
        exit(1);
    }
    return option->strArray[0];
}

/**
 * return the first boolean value from option. If the option has not been
 * set the default value is returned.  Option must have only 1 value.
 *
 * @param option option to use
 * @return first value of the option
 */
int clo_getOptionBool(clo_option_t *option) {
    if (option->dataType != CLO_TYPE_BOOL) {
        fprintf(stderr,
                "-E- clo_getOptionBool: option=%s is not a boolean option\n",
                option->key);
        exit(1);
    }
    if (option->count <= 0) {
        return 0;
    }
    if (option->count != 1) {
        fprintf(stderr,
                "-E- clo_getOptionBool: option=%s count needs to be 1\n",
                option->key);
        exit(1);
    }
    return ((int*) (option->valArray))[0];
}

/**
 * return the first integer value from option. If the option has not been
 * set the default value is returned.  Option must have only 1 value.
 *
 * @param option option to use
 * @return first value of the option
 */
int clo_getOptionInt(clo_option_t *option) {
    if (option->dataType != CLO_TYPE_INT) {
        fprintf(stderr,
                "-E- clo_getOptionInt: option=%s is not an integer option\n",
                option->key);
        exit(1);
    }
    if (option->count <= 0) {
        fprintf(stderr, "-E- clo_getOptionInt: option=%s needs to be set\n",
                option->key);
        exit(1);
    }
    if (option->count != 1) {
        fprintf(stderr, "-E- clo_getOptionInt: option=%s count needs to be 1\n",
                option->key);
        exit(1);
    }
    return ((int*) (option->valArray))[0];
}

/**
 * return the first integer value from option. If the option has not been
 * set the default value is returned.  Option must have only 1 value.
 *
 * @param option option to use
 * @return first value of the option
 */
int64_t clo_getOptionInt64(clo_option_t *option) {
    if (option->dataType != CLO_TYPE_INT64) {
        fprintf(stderr,
                "-E- clo_getOptionInt64: option=%s is not a 64 bit integer option\n",
                option->key);
        exit(1);
    }
    if (option->count <= 0) {
        fprintf(stderr, "-E- clo_getOptionInt64: option=%s needs to be set\n",
                option->key);
        exit(1);
    }
    if (option->count != 1) {
        fprintf(stderr, "-E- clo_getOptionInt: option=%s count needs to be 1\n",
                option->key);
        exit(1);
    }
    return ((int64_t*) (option->valArray))[0];
}

/**
 * return the first float value from option. If the option has not been
 * set the default value is returned.  Option must have only 1 value.
 *
 * @param option option to use
 * @return first value of the option
 */
float clo_getOptionFloat(clo_option_t *option) {
    if (option->dataType != CLO_TYPE_FLOAT) {
        fprintf(stderr,
                "-E- clo_getOptionFloat: option=%s is not a float option\n",
                option->key);
        exit(1);
    }
    if (option->count <= 0) {
        fprintf(stderr, "-E- clo_getOptionFloat: option=%s needs to be set\n",
                option->key);
        exit(1);
    }
    if (option->count != 1) {
        fprintf(stderr,
                "-E- clo_getOptionFloat: option=%s count needs to be 1\n",
                option->key);
        exit(1);
    }
    return ((float*) (option->valArray))[0];
}

/**
 * return the first double value from option. If the option has not been
 * set the default value is returned.  Option must have only 1 value.
 *
 * @param option option to use
 * @return first value of the option
 */
double clo_getOptionDouble(clo_option_t *option) {
    if (option->dataType != CLO_TYPE_DOUBLE) {
        fprintf(stderr,
                "-E- clo_getOptionDouble: option=%s is not a double option\n",
                option->key);
        exit(1);
    }
    if (option->count <= 0) {
        fprintf(stderr, "-E- clo_getOptionDouble: option=%s needs to be set\n",
                option->key);
        exit(1);
    }
    if (option->count != 1) {
        fprintf(stderr,
                "-E- clo_getOptionDouble: option=%s count needs to be 1\n",
                option->key);
        exit(1);
    }
    return ((double*) (option->valArray))[0];
}

/**
 * Return the count and string array from the given option. If the option
 * has not been set the default value is returned.  The pointer returned
 * references internally allocated memory stored in the option structure,
 * so it should be valid until the value of option is changed.

 * @param option option to get the information from
 * @param count [out] pointer where number of strings in array is written
 * @return array of strings
 */
char** clo_getOptionStrings(clo_option_t *option, int *count) {
    if (option->dataType != CLO_TYPE_STRING
            && option->dataType != CLO_TYPE_IFILE
            && option->dataType != CLO_TYPE_OFILE
            && option->dataType != CLO_TYPE_HELP) {
        fprintf(stderr,
                "-E- clo_getOptionStrings: option=%s is not a string option\n",
                option->key);
        exit(1);
    }
    assert(count);
    *count = option->count;
    return option->strArray;
}

/**
 * Return the count and boolean array from the given option. If the option
 * has not been set the default value is returned.  The pointer returned
 * references internally allocated memory stored in the option structure,
 * so it should be valid until the value of option is changed.

 * @param option option to get the information from
 * @param count [out] pointer where number of booleans in array is written
 * @return array of booleans
 */
int* clo_getOptionBools(clo_option_t *option, int *count) {
    if (option->dataType != CLO_TYPE_BOOL) {
        fprintf(stderr,
                "-E- clo_getOptionBools: option=%s is not a boolean option\n",
                option->key);
        exit(1);
    }
    assert(count);
    *count = option->count;
    return (int*) option->valArray;
}

/**
 * Return the count and int array from the given option. If the option
 * has not been set the default value is returned.  The pointer returned
 * references internally allocated memory stored in the option structure,
 * so it should be valid until the value of option is changed.

 * @param option option to get the information from
 * @param count [out] pointer where number of ints in array is written
 * @return array of ints
 */
int* clo_getOptionInts(clo_option_t *option, int *count) {
    if (option->dataType != CLO_TYPE_INT) {
        fprintf(stderr,
                "-E- clo_getOptionInts: option=%s is not an integer option\n",
                option->key);
        exit(1);
    }
    assert(count);
    *count = option->count;
    return (int*) option->valArray;
}

/**
 * Return the count and float array from the given option. If the option
 * has not been set the default value is returned.  The pointer returned
 * references internally allocated memory stored in the option structure,
 * so it should be valid until the value of option is changed.

 * @param option option to get the information from
 * @param count pointer where number of floats in array is written
 * @return array of floats
 */
float* clo_getOptionFloats(clo_option_t *option, int *count) {
    if (option->dataType != CLO_TYPE_FLOAT) {
        fprintf(stderr,
                "-E- clo_getOptionFloats: option=%s is not a float option\n",
                option->key);
        exit(1);
    }
    assert(count);
    *count = option->count;
    return (float*) option->valArray;
}

/**
 * Return the count and double array from the given option. If the option
 * has not been set the default value is returned.  The pointer returned
 * references internally allocated memory stored in the option structure,
 * so it should be valid until the value of option is changed.

 * @param option option to get the information from
 * @param count [out] pointer where number of doubles in array is written
 * @return array of doubles
 */
double* clo_getOptionDoubles(clo_option_t *option, int *count) {
    if (option->dataType != CLO_TYPE_DOUBLE) {
        fprintf(stderr,
                "-E- clo_getOptionDoubles: option=%s is not a double option\n",
                option->key);
        exit(1);
    }
    assert(count);
    *count = option->count;
    return (double*) option->valArray;
}

/**
 * Find the option "key" and return the raw string value. Exit is called
 * if "key" is not found.  If the option has not been set the default
 * value is returned.  Exit is called if the option has not been set
 * and the default value is NULL.
 *
 * @param list list to search through
 * @param key option key to find
 * @return internal pointer to the raw string value
 */
char* clo_getRawString(clo_optionList_t *list, const char *key) {
    clo_option_t* option = clo_findOption(list, key);
    if (option == NULL) {
        fprintf(stderr, "-E- clo_getRawString: option=%s not found\n", key);
        exit(1);
    }
    return clo_getOptionRawString(option);
}

/**
 * Find the option "key" and return the first value. exit is called if
 * "key" is not found.  If the option has not been set the default
 * value is returned.
 *
 * @param list list of options to search
 * @param key option name to search for
 * @return internal pointer to the first string value
 */
char* clo_getString(clo_optionList_t *list, const char *key) {
    clo_option_t* option = clo_findOption(list, key);
    if (option == NULL) {
        fprintf(stderr, "-E- clo_getString: option=%s not found\n", key);
        exit(1);
    }
    return clo_getOptionString(option);
}

/**
 * Find the option "key" and return the first value. exit is called if
 * "key" is not found.  If the option has not been set the default
 * value is returned.
 *
 * @param list list of options to search
 * @param key option name to search for
 * @return internal pointer to the first boolean value
 */
int clo_getBool(clo_optionList_t *list, const char *key) {
    clo_option_t* option = clo_findOption(list, key);
    if (option == NULL) {
        fprintf(stderr, "-E- clo_getBool: option=%s not found\n", key);
        exit(1);
    }
    return clo_getOptionBool(option);
}

/**
 * Find the option "key" and return the first value. exit is called if
 * "key" is not found.  If the option has not been set the default
 * value is returned.
 *
 * @param list list of options to search
 * @param key option name to search for
 * @return internal pointer to the first int value
 */
int clo_getInt(clo_optionList_t *list, const char *key) {
    clo_option_t* option = clo_findOption(list, key);
    if (option == NULL) {
        fprintf(stderr, "-E- clo_getInt: option=%s not found\n", key);
        exit(1);
    }
    return clo_getOptionInt(option);
}

/**
 * Find the option "key" and return the first value. exit is called if
 * "key" is not found.  If the option has not been set the default
 * value is returned.
 *
 * @param list list of options to search
 * @param key option name to search for
 * @return internal pointer to the first int value
 */
int64_t clo_getInt64(clo_optionList_t *list, const char *key) {
    clo_option_t* option = clo_findOption(list, key);
    if (option == NULL) {
        fprintf(stderr, "-E- clo_getInt: option=%s not found\n", key);
        exit(1);
    }
    return clo_getOptionInt64(option);
}

/**
 * Find the option "key" and return the first value. exit is called if
 * "key" is not found.  If the option has not been set the default
 * value is returned.
 *
 * @param list list of options to search
 * @param key option name to search for
 * @return internal pointer to the first float value
 */
float clo_getFloat(clo_optionList_t *list, const char *key) {
    clo_option_t* option = clo_findOption(list, key);
    if (option == NULL) {
        fprintf(stderr, "-E- clo_getFloat: option=%s not found\n", key);
        exit(1);
    }
    return clo_getOptionFloat(option);
}

/**
 * Find the option "key" and return the first value. exit is called if
 * "key" is not found.  If the option has not been set the default
 * value is returned.
 *
 * @param list list of options to search
 * @param key option name to search for
 * @return internal pointer to the first double value
 */
double clo_getDouble(clo_optionList_t *list, const char *key) {
    clo_option_t* option = clo_findOption(list, key);
    if (option == NULL) {
        fprintf(stderr, "-E- clo_getDouble: option=%s not found\n", key);
        exit(1);
    }
    return clo_getOptionDouble(option);
}

/**
 * Find the option "key" and return the count and string array. exit is
 * called if "key" is not found.  If the option has not been set the
 * default value is returned.  The pointer returned references
 * internally allocated memory stored in the option structure, so it
 * should be valid until the value of option is changed.
 *
 * @param list list of options to search through
 * @param key option name to search for
 * @param count [out] number of items in the array
 * @return an array of strings.  This pointer is pointing to internal
 *         memory stored in the option structure.
 */
char** clo_getStrings(clo_optionList_t *list, const char *key, int *count) {
    assert(count);
    clo_option_t* option = clo_findOption(list, key);
    if (option == NULL) {
        fprintf(stderr, "-E- clo_getStrings: option=%s not found\n", key);
        exit(1);
    }
    return clo_getOptionStrings(option, count);
}

/**
 * Find the option "key" and return the count and boolean array. exit is
 * called if "key" is not found.  If the option has not been set the
 * default value is returned.  The pointer returned references
 * internally allocated memory stored in the option structure, so it
 * should be valid until the value of option is changed.
 *
 * @param list list of options to search through
 * @param key option name to search for
 * @param count [out] number of items in the array
 * @return an array of booleans.  This pointer is pointing to internal
 *         memory stored in the option structure.
 */
int* clo_getBools(clo_optionList_t *list, const char *key, int *count) {
    assert(count);
    clo_option_t* option = clo_findOption(list, key);
    if (option == NULL) {
        fprintf(stderr, "-E- clo_getBools: option=%s not found\n", key);
        exit(1);
    }
    return clo_getOptionBools(option, count);
}

/**
 * Find the option "key" and return the count and int array. exit is
 * called if "key" is not found.  If the option has not been set the
 * default value is returned.  The pointer returned references
 * internally allocated memory stored in the option structure, so it
 * should be valid until the value of option is changed.
 *
 * @param list list of options to search through
 * @param key option name to search for
 * @param count [out] number of items in the array
 * @return an array of ints.  This pointer is pointing to internal
 *         memory stored in the option structure.
 */
int* clo_getInts(clo_optionList_t *list, const char *key, int *count) {
    assert(count);
    clo_option_t* option = clo_findOption(list, key);
    if (option == NULL) {
        fprintf(stderr, "-E- clo_getInts: option=%s not found\n", key);
        exit(1);
    }
    return clo_getOptionInts(option, count);
}

/**
 * Find the option "key" and return the count and float array. exit is
 * called if "key" is not found.  If the option has not been set the
 * default value is returned.  The pointer returned references
 * internally allocated memory stored in the option structure, so it
 * should be valid until the value of option is changed.
 *
 * @param list list of options to search through
 * @param key option name to search for
 * @param count [out] number of items in the array
 * @return an array of floats.  This pointer is pointing to internal
 *         memory stored in the option structure.
 */
float* clo_getFloats(clo_optionList_t *list, const char *key, int *count) {
    assert(count);
    clo_option_t* option = clo_findOption(list, key);
    if (option == NULL) {
        fprintf(stderr, "-E- clo_getFloats: option=%s not found\n", key);
        exit(1);
    }
    return clo_getOptionFloats(option, count);
}

/**
 * Find the option "key" and return the count and float array. exit is
 * called if "key" is not found.  If the option has not been set the
 * default value is returned.  The pointer returned references
 * internally allocated memory stored in the option structure, so it
 * should be valid until the value of option is changed.
 *
 * @param list list of options to search through
 * @param key option name to search for
 * @param count [out] number of items in the array
 * @return an array of doubles.  This pointer is pointing to internal
 *         memory stored in the option structure.
 */
double* clo_getDoubles(clo_optionList_t *list, const char *key, int *count) {
    assert(count);
    clo_option_t* option = clo_findOption(list, key);
    if (option == NULL) {
        fprintf(stderr, "-E- clo_getDoubles: option=%s not found\n", key);
        exit(1);
    }
    return clo_getOptionDoubles(option, count);
}

/**
 * Set the string for this option and run any callback
 *
 * @param option structure to modify
 * @param val string to enter into the structure to have parsed
 * @param source where did this value come from
 * @return 0 if OK, -1 if error
 */
int clo_setOptionString(clo_option_t *option, const char *val,
        const char *source) {
    assert(option);

    char *ptr;
    char *deleteAfter;
    int i;

    // delete the passed in string after copying it
    deleteAfter = option->valStr;
    option->valStr = NULL;

    // delete the old parsed values
    if (option->strArray) {
        for (i = 0; i < option->count; i++) {
            if (option->strArray[i])
                free(option->strArray[i]);
        }
        free(option->strArray);
        option->strArray = NULL;
    } // if strArray
    option->count = 0;
    if (option->valArray) {
        free(option->valArray);
        option->valArray = NULL;
    }

    // assign the new value string
    if (val == NULL) {
        if (option->dataType == CLO_TYPE_BOOL) {
            option->valStr = strdup("true");
        } else {
            fprintf(stderr,
                    "-E- clo_setOptionString: No value given for option key=%s\n",
                    option->key);
            exit(1);
        }
    } else {
        // trim spaces and quotes off the front of value
        while (isspace(*val) || (*val == '"'))
            val++;
        ptr = strdup(val);
        i = strlen(ptr) - 1;
        while (i >= 0) {
            if (isspace(ptr[i]) || (ptr[i] == '"'))
                ptr[i] = '\0';
            else
                break;
            i--;
        }
        option->valStr = ptr;
    }

    // delete old value string
    if (deleteAfter)
        free(deleteAfter);

    // assign the source so we can track the options source
    deleteAfter = option->source;
    if (source)
        option->source = strdup(source);
    else
        option->source = strdup("undefined");

    // delete old source string
    if (deleteAfter)
        free(deleteAfter);

    // fill up the strArray, valArray and count
    clo_parseOption(option);

    // call the callback routine if it has one
    if (option->cb) {
        (option->cb)(option);
    }

    return 0;
}

/**
 * Set the string for this option and run any callback
 *
 * @param list option list to search
 * @param key option key to search for
 * @param val string to enter into the structure to have parsed
 * @param source where did this value come from
 * @return 0 if OK, -1 if error
 */
int clo_setString(clo_optionList_t *list, const char* key,
        const char *val, const char *source) {
    assert(list);
    assert(key);
    assert(val);
    clo_option_t* option = clo_findOption(list, key);
    if (option == NULL) {
        fprintf(stderr, "-E- clo_setString: option=%s not found\n", key);
        return -1;
    }
    return clo_setOptionString(option, val, source);
}

/**
 * set whether or not the library will accept positional options
 *
 * @param val 1 = enable processing, 0 = disable
 */
void clo_setEnablePositionOptions(int val) {
    clo_enablePositionOptions = val;
}

/**
 * is the library set to automatically add options to
 * the list when reading the command line or par file.
 *
 * @return 1 = enable adding, 0 = disable
 */
int clo_getEnablePositionOptions() {
    return clo_enablePositionOptions;
}

/**
 * get the number of position options in the list
 * @param list list to get the count for
 * @return number of options
 */
int clo_getPositionNumOptions(clo_optionList_t *list) {
    assert(list);
    return list->positionNumOptions;
}

/**
 * get the option at index i
 *
 * @param list list to search through
 * @param i index into the list
 * @return option pointer at i, or NULL if not found
 */
clo_option_t* clo_getPositionOption(clo_optionList_t *list, int i) {
    assert(list);
    if (i < 0 || i >= list->positionNumOptions)
        return NULL;
    return list->positionOptions[i];
}

char* clo_getPositionString(clo_optionList_t *list, int pos) {
    clo_option_t* option = clo_getPositionOption(list, pos);
    if (option == NULL) {
        fprintf(stderr, "-E- clo_getPositionString: option %d not found on command line\n", pos);
        exit(1);
    }
    return clo_getOptionRawString(option);
}

/**
 * set whether or not the library will add the par file option
 *     this must be set before calling clo_createList()
 *
 * @param val 1 = enable processing, 0 = disable
 */
void clo_setEnableParOption(int val) {
    clo_enableParOption = val;
}

/**
 * is the library set to automatically add the par option
 *
 * @return 1 = enable, 0 = disable
 */
int clo_getEnableParOption() {
    return clo_enableParOption;
}

/**
 * print the current value of the option
 *
 * @param option option pointer to print
 */
void clo_printOptionVal(clo_option_t *option) {
    int i;

    assert(option);
    if (option->count > 1)
        printf("[");

    switch (option->dataType) {
    case CLO_TYPE_BOOL:
    {
        int* vals = (int*) option->valArray;
        if (option->count <= 0)
            printf("<noValue>");
        else {
            if (vals[0])
                printf("true");
            else
                printf("false");
            for (i = 1; i < option->count; i++) {
                if (vals[i])
                    printf(",true");
                else
                    printf(",false");
            } // for
        }
    }
        break;
    case CLO_TYPE_INT:
    {
        int* vals = (int*) option->valArray;
        if (option->count <= 0)
            printf("<noValue>");
        else {
            printf("%d", vals[0]);
            for (i = 1; i < option->count; i++)
                printf(",%d", vals[i]);
        }
    }
        break;
    case CLO_TYPE_FLOAT:
    {
        float* vals = (float*) option->valArray;
        if (option->count <= 0)
            printf("<noValue>");
        else {
            printf("%f", vals[0]);
            for (i = 1; i < option->count; i++)
                printf(",%f", vals[i]);
        }
    }
        break;
    case CLO_TYPE_DOUBLE:
    {
        double* vals = (double*) option->valArray;
        if (option->count <= 0)
            printf("<noValue>");
        else {
            printf("%f", vals[0]);
            for (i = 1; i < option->count; i++)
                printf(",%f", vals[i]);
        }
    }
        break;
    case CLO_TYPE_STRING:
    case CLO_TYPE_IFILE:
    case CLO_TYPE_OFILE:
        if (option->count <= 0)
            printf("<noValue>");
        else {
            printf("%s", option->strArray[0]);
            for (i = 1; i < option->count; i++)
                printf(",%s", option->strArray[i]);
        }
        break;
    case CLO_TYPE_HELP:
        break;
    default:
        printf("unknown type");
        break;
    } // switch

    if (option->count > 1)
        printf("]");
}

/**
 * print the options that have been set in a terse manner
 * @param list option list to print
 */
void clo_printVals(clo_optionList_t *list) {
    int i;

    for (i = 0; i < list->numOptions; i++) {
        if (list->options[i]->dataType == CLO_TYPE_HELP)
            continue;
        if (list->options[i]->valStr) {
            printf("%s=", list->options[i]->key);
            clo_printOptionVal(list->options[i]);
            printf("\n");
        }
    } // for options

}

/**
 * print the full description of the option
 *
 * @param option option to print
 */
void clo_printOption(clo_option_t *option) {
    int i;
    assert(option);
    char* dataTypeStr;

    if (option->dataType == CLO_TYPE_HELP) {
        if (option->desc)
            printf("   %s", option->desc);
        return;
    }

    dataTypeStr = clo_dataTypeToString(option->dataType);

    printf("   %s (%s)", option->key, dataTypeStr);
    if (option->numAliases > 0) {
        printf(" (alias=%s", option->aliases[0]);
        for (i = 1; i < option->numAliases; i++)
            printf(",%s", option->aliases[i]);
        printf(")");
    }
    if (option->defaultVal)
        printf(" (default=%s)", option->defaultVal);
    if (option->desc)
        printf(" = %s", option->desc);
}

/**
 * print the full description and value of an option
 *
 * @param option option to dump
 */
void clo_dumpOption(clo_option_t *option) {
    assert(option);

    int i;
    char* dataTypeStr;

    if (option->dataType == CLO_TYPE_HELP)
        return;

    dataTypeStr = clo_dataTypeToString(option->dataType);

    printf("   %s (%s)", option->key, dataTypeStr);
    if (option->numAliases > 0) {
        printf(" (alias=%s", option->aliases[0]);
        for (i = 1; i < option->numAliases; i++)
            printf(",%s", option->aliases[i]);
        printf(")");
    }
    if (option->defaultVal)
        printf(" (default=%s)", option->defaultVal);

    if (option->valStr) {
        printf(" (current=");
        clo_printOptionVal(option);
        printf(")");
    }
    if (option->source)
        printf(" (source=%s)", option->source);
}

/**
 * print the full description and value of all of the options
 *
 * @param list option list to print
 */
void clo_printOptions(clo_optionList_t *list) {
    int i;
    char *key;
    clo_option_t *option;

    assert(list);
    if (selectOptionKeys) {
        i = 0;
        key = selectOptionKeys[i];
        while (key) {
            option = clo_findOption(list, key);
            if (option) {
                if (option->desc) { /* only print if description */
                    clo_printOption(option);
                    printf("\n");
                } // description exists
            } else {
                fprintf(stderr,
                        "clo_printOptions - Could not find option \"%s\" in option list.\n",
                        key);
            }
            i++;
            key = selectOptionKeys[i];
        } // while key
    } else {
        for (i = 0; i < list->numOptions; i++) {
            if (list->options[i]->desc) { /* only print if description */
                clo_printOption(list->options[i]);
                printf("\n");
            }
        }
    } // for options
}

/**
 * print the version information
 */
void clo_printVersion() {
    if (clo_versionString)
        printf("%s\n", clo_versionString);
    else
        printf("Version String not set.\n");
}

/**
 * print the help string
 */
void clo_printHelpString() {
    if (clo_helpString)
        printf("%s\n", clo_helpString);
}

/**
 * print the program usage information
 *
 * @param list option list to print out
 */
void clo_printUsage(clo_optionList_t *list) {
    clo_printVersion();
    clo_printHelpString();
    clo_printOptions(list);
}

/**
 * print version string and dump all options in the list
 *
 * @param list option list to print out
 */
void clo_dumpOptions(clo_optionList_t *list) {
    int i;

    assert(list);
    clo_printVersion();
    for (i = 0; i < list->numOptions; i++) {
        clo_dumpOption(list->options[i]);
        printf("\n");
    } // for options
}

/**
 * Read a string assuming it is one command line option.  Lookup the key
 * in list and set the value and the source.
 *
 * @param list option list to search
 * @param str string containing "key=value"
 * @param source source of the string.  Either "command line" or
 *               the file name
 */
void clo_readString(clo_optionList_t *list, const char *str, const char *source) {
    assert(list);
    assert(str);
    assert(source);
    char* newStr = strdup(str);
    char* keyStr;
    char* valStr;

    // see if there is an '=' in the str
    char* equalPtr = strchr(newStr, '=');
    if (equalPtr) { /* equal found */
        *equalPtr = '\0';
        keyStr = newStr;
        valStr = equalPtr + 1;
    } else { /* equal not found */
        trimBlanks(newStr);
        if (newStr[0] == '-')
            valStr = "true";
        else
            valStr = NULL;
        keyStr = newStr;
    }

    keyStr = clo_envExpandString(keyStr);
    if (valStr) {
        valStr = clo_envExpandString(valStr);
    }

    // find the option 
    clo_option_t* option = clo_findOption(list, keyStr);

    // if option not found and enableExtra Options is set,
    // add it and set the description field to "UNDEFINED_OPTION"
    if (option == NULL) {
        if (valStr) {
            if (clo_enableExtraOptions) {
                option = clo_addOption(list, keyStr, CLO_TYPE_STRING,
                        NULL, "UNDEFINED_OPTION");
            } else {
                fprintf(stderr, "-E- clo_readString: unknown option \"%s\" from %s\n",
                        str, source);
                exit(1);
            }
        } else {
            if (strcmp(source, "command line") == 0) {
                if (clo_enablePositionOptions) {
                    option = clo_addPositionOption(list);
                    valStr = keyStr;
                    keyStr = NULL;
                } else {
                    fprintf(stderr, "-E- clo_readString: unknown option \"%s\" from %s\n",
                            str, source);
                    exit(1);
                }
            } else {
                if (clo_enableExtraOptions) {
                    option = clo_addOption(list, keyStr, CLO_TYPE_BOOL,
                            NULL, "UNDEFINED_OPTION");
                } else {
                    fprintf(stderr, "-E- clo_readString: unknown option \"%s\" from %s\n",
                            str, source);
                    exit(1);
                }
            }
        }
    }

    clo_setOptionString(option, valStr, source);

    // clean up the memory
    free(newStr);
    free(keyStr);
    free(valStr);
}

/**
 * Read command line arguments handling function.
 * First call to readArgsPar with file descending enabled, then disabled to allow Command Line options override
 *
 * @param list option list to search
 * @param argc number of elements in argv
 * @param argv array of command line argument strings
 */

void clo_readArgs(clo_optionList_t *list, int argc, char *argv[]) {

    clo_readArgsPar(list, argc, argv, 1);
    clo_readArgsPar(list, argc, argv, 0);

}

/**
 * Read command line arguments.  Lookup the key in list and set the
 * value.
 *
 * @param list option list to search
 * @param argc number of elements in argv
 * @param argv array of command line argument strings
 * @param enableFileDecending enable file descending
 */
void clo_readArgsPar(clo_optionList_t *list, int argc, char *argv[], int enableFileDescending) {
    int i;

    // first clear out the position params if there are any
    clo_clearPositionOptions(list);

    clo_option_t* option;
    //Find the par option and set enableFileDecending
    option = clo_findOption(list, "par");
    if(option)
        option->cb_flag = enableFileDescending;
    for (i = 1; i < argc; i++) {
        clo_readString(list, argv[i], "command line");
    }

    // check for dump
    if (clo_enableDumpOptions) {
        int dumpExit = 0;

        if (clo_getBool(list, "dump_options")) {
            clo_dumpOptions(list);
            dumpExit = 1;
        }

        option = clo_findOption(list, "dump_options_paramfile");
        if (clo_isOptionSet(option)) {
            printf("writing options param file to %s\n",
                    clo_getOptionString(option));
            clo_writeParameterFile(list, clo_getOptionString(option));
            dumpExit = 1;
        }

        option = clo_findOption(list, "dump_options_xmlfile");
        if (clo_isOptionSet(option)) {
            printf("writing options XML file to %s\n",
                    clo_getOptionString(option));
            clo_writeXmlFile(list, clo_getOptionString(option));
            dumpExit = 1;
        }

        if (dumpExit) {
            exit(0);
        }

    }

}

/**
 * copy (overwrite) options from readList that have been set
 *
 * @param list option list to copy to
 * @param readList option list to copy from
 */
void clo_readOptions(clo_optionList_t *list, clo_optionList_t *readList) {
    int i;
    clo_option_t *option;
    char *str;

    assert(list);
    assert(readList);

    for (i = 0; i < readList->numOptions; i++) {
        option = readList->options[i];
        if (option->valStr) {
            str = (char*) malloc(
                    strlen(option->key) + strlen(option->valStr) + 5);
            sprintf(str, "%s=%s", option->key, option->valStr);
            clo_readString(list, str, option->source);
            free(str);
        }
    }

    // delete position options first
    for (i = 0; i < list->positionNumOptions; i++) {
        clo_deleteOption(list->positionOptions[i]);
    }
    list->positionNumOptions = 0;

    for (i = 0; i < readList->positionNumOptions; i++) {
        option = clo_copyOption(readList->positionOptions[i]);
        clo_addToArray((void***) &(list->positionOptions), &(list->positionStorageSize),
                &(list->positionNumOptions), option);
    }
}

/**
 * read a par file
 * @param list option list to write to
 * @param fileName parameter file to read
 */
void clo_readFile(clo_optionList_t *list, const char *fileName) {
    assert(list);
    assert(fileName);

    //    if(clo_checkFileRecursion(fileName))
    //        return;

    FILE *fp;
    char line[SILLYSTRING];
    char *ptr;
    int lineNumber = 0;
    char *sourceStr;

    // copy the string since the pointer passed in might get freed
    // while we are processing options from the file.
    sourceStr = strdup(fileName);

    if ((fp = fopen(fileName, "r")) == NULL) {
        fprintf(stderr, "-E- clo_readFile: Can't open parameter file - %s\n",
                fileName);
        exit(1);
    }

    while ((fgets(line, SILLYSTRING, fp)) != NULL) {
        lineNumber++;

        /* skip the comment or blank line */
        ptr = line;
        while (isspace(*ptr))
            ptr++;
        if (*ptr == '#' || *ptr == ';' || *ptr == '\0')
            continue;

        clo_readString(list, ptr, sourceStr);

    } // while line

    free(sourceStr);
    fclose(fp);
}

/**
 * Was the option explicitly set.
 *
 * @param option option to check
 * @return 1 if option was set else 0.  0 is returned if option not found.
 */
int clo_isOptionSet(clo_option_t *option) {
    if (option->valStr)
        return 1;
    return 0;
}

/**
 * Was the option explicitly set.
 *
 * @param list list of options to search through
 * @param key name of option to search for
 * @return 1 if option was set else 0.  0 is returned if option not found.
 */
int clo_isSet(clo_optionList_t *list, const char *key) {
    clo_option_t* option;

    option = clo_findOption(list, key);
    if (option)
        return clo_isOptionSet(option);

    return 0;
}

/**
 * Write the list to a parameter file
 *
 * @param list option list to write
 * @param filename file name to write to
 * @return 0 if OK else -1
 */
int clo_writeParameterFile(clo_optionList_t *list, const char *filename) {
    clo_option_t* option;
    int i;
    FILE *fout;
    char *key;

    assert(list);
    assert(filename);

    // open file
    fout = fopen(filename, "w");
    if (fout == NULL) {
        fprintf(stderr,
                "clo_writeParameterFile - Could not open \"%s\" for writing.\n",
                filename);
        return -1;
    }

    if (selectOptionKeys) {
        i = 0;
        key = selectOptionKeys[i];
        while (key) {
            option = clo_findOption(list, key);
            if (option) {
                if (clo_isOptionSet(option)
                        && strcmp(option->key, "dump_options_paramfile")
                        && strcmp(option->key, "par")) {
                    fprintf(fout, "%s = %s\n", option->key, option->valStr);
                } // description exists
            } else {
                fprintf(stderr,
                        "clo_writeParameterFile - Could not find option \"%s\" in option list.\n",
                        key);
            }
            i++;
            key = selectOptionKeys[i];
        } // while key
    } else {
        for (i = 0; i < list->numOptions; i++) {
            option = list->options[i];
            if (clo_isOptionSet(option)
                    && strcmp(option->key, "dump_options_paramfile")
                    && strcmp(option->key, "par")) {
                fprintf(fout, "%s = %s\n", option->key, option->valStr);
            }
        }
    } // for options

    fclose(fout);

    return 0;
}

/**
 * add program metadata XML entries
 * @param tag XML tag to add
 * @param value XML value to add
 */
void clo_addXmlProgramMetadata(const char* tag, const char* value) {
    clo_programMetadata_t* metadata = (clo_programMetadata_t*) malloc(
            sizeof (clo_programMetadata_t));
    assert(metadata);
    metadata->tag = strdup(tag);
    assert(metadata->tag);
    metadata->value = strdup(value);
    assert(metadata->value);
    clo_addToArray((void***) &(clo_programMetadataList.entries),
            &(clo_programMetadataList.storageSize),
            &(clo_programMetadataList.numEntries), (void*) metadata);
}

/**
 * clear all of the program metadata XML entries
 */
void clo_clearXmlProgramMetadata() {
    int i;

    for (i = 0; i < clo_programMetadataList.numEntries; i++) {
        free(clo_programMetadataList.entries[i]->tag);
        free(clo_programMetadataList.entries[i]->value);
    }
    free(clo_programMetadataList.entries);
    clo_programMetadataList.entries = NULL;
    clo_programMetadataList.numEntries = 0;
    clo_programMetadataList.storageSize = 0;
}

/**
 * write an XML start tag
 *
 * @param fout file handle to write to
 * @param level indent level
 * @param tag XML tag string
 */
void clo_writeXmlStartTag(FILE *fout, int level, const char *tag) {
    int i;
    for (i = 0; i < level; i++)
        fprintf(fout, XML_INDENT_STRING);
    fprintf(fout, "<%s>\n", tag);
}

/**
 * write an XML end tag
 *
 * @param fout file handle to write to
 * @param level indent level
 * @param tag XML tag string
 */
void clo_writeXmlEndTag(FILE *fout, int level, const char *tag) {
    int i;
    for (i = 0; i < level; i++)
        fprintf(fout, XML_INDENT_STRING);
    fprintf(fout, "</%s>\n", tag);
}

/**
 * write an XML tag and escape the value with CDATA is necessary
 *
 * @param fout file handle to write to
 * @param level indent level
 * @param tag XML tag string
 * @param value XML value string
 */
void clo_writeXmlTag(FILE *fout, int level, const char *tag, const char *value) {
    int i;
    for (i = 0; i < level; i++)
        fprintf(fout, XML_INDENT_STRING);
    if (value == NULL || strlen(value) == 0) {
        fprintf(fout, "<%s/>\n", tag);
    } else if (strpbrk(value, "<>&"))
        fprintf(fout, "<%s><![CDATA[%s]]></%s>\n", tag, value, tag);
    else
        fprintf(fout, "<%s>%s</%s>\n", tag, value, tag);
}

/**
 * write valid value XML tag
 *
 * @param fout file handle to write to
 * @param level indent level
 * @param value value to write in the tag
 * @param description description of this valid value
 */
static void writeXmlValidValue(FILE *fout, int level, const char *value,
        const char *description) {
    clo_writeXmlStartTag(fout, level, "validValue");
    clo_writeXmlTag(fout, level + 1, "value", value);
    clo_writeXmlTag(fout, level + 1, "description", description);
    clo_writeXmlEndTag(fout, level, "validValue");
}

/**
 * write the description tag
 *
 * @param fout file handle to write to
 * @param level indent level
 * @param str long string to write
 */
static void writeXmlDescription(FILE *fout, int level, const char *str) {
    char *buff, *buff2;
    char *line, *line2;
    char *val, *val2;
    char *desc, *desc2;
    char *saveLine, *saveVal;
    int i;
    int foundFirst = 0;
    int wroteFirst = 0;
    int foundValue = 0;

    buff = strdup(str);
    i = strlen(str);
    buff2 = (char*) malloc(i + 2);
    buff2[0] = '\0';
    line2 = (char*) malloc(i + 1);
    line2[0] = '\0';
    val2 = (char*) malloc(i + 1);
    val2[0] = '\0';
    desc2 = (char*) malloc(i + 1);
    desc2[0] = '\0';

    // loop through lines
    line = strtok_r(buff, "\n", &saveLine);
    while (line) {
        foundValue = 0;
        strcpy(line2, line); // save a copy of the line before strtok
        val = strtok_r(line, ":", &saveVal);
        desc = strtok_r(NULL, ":", &saveVal);
        if (val && desc) {
            trimBlanks(val);
            trimBlanks(desc);

            // only use if val is a single word
            if (strpbrk(val, " \t\r\n") == NULL) {
                if (foundFirst) {
                    if (!wroteFirst) {
                        clo_writeXmlStartTag(fout, level, "validValues");
                        wroteFirst = 1;
                    }
                    writeXmlValidValue(fout, level + 1, val2, desc2);
                }
                foundFirst = 1;
                foundValue = 1;
                strcpy(val2, val);
                strcpy(desc2, desc);
            }
        }

        if (foundFirst) {
            if (!foundValue) {
                trimBlanks(line2);
                strcat(desc2, " ");
                strcat(desc2, line2);
            }
        } else {
            trimBlanks(line2);
            strcat(buff2, line2);
            strcat(buff2, " ");
        }
        line = strtok_r(NULL, "\n", &saveLine);
    }

    if (foundFirst) {
        writeXmlValidValue(fout, level + 1, val2, desc2);
        clo_writeXmlEndTag(fout, level, "validValues");
    }

    trimBlanks(buff2);
    clo_writeXmlTag(fout, level, "description", buff2);

    free(buff);
    free(buff2);
    free(line2);
    free(val2);
    free(desc2);
}

/**
 * write the whole option to the XML file
 *
 * @param fout file handle to write to
 * @param option pointer to option structure to write out
 */
static void writeXmlOption(FILE *fout, clo_option_t* option) {
    if (option->dataType == CLO_TYPE_HELP)
        return;

    fprintf(fout, "        <option type=\"%s\">\n",
            clo_dataTypeToString(option->dataType));
    clo_writeXmlTag(fout, 3, "name", option->key);
    if (option->valStr)
        clo_writeXmlTag(fout, 3, "value", option->valStr);
    else if (option->defaultVal)
        clo_writeXmlTag(fout, 3, "value", option->defaultVal);
    else
        clo_writeXmlStartTag(fout, 3, "value/");
    if (option->defaultVal)
        clo_writeXmlTag(fout, 3, "default", option->defaultVal);
    clo_writeXmlTag(fout, 3, "source", option->source);
    int j;
    if (option->numAliases > 0) {
        clo_writeXmlStartTag(fout, 3, "aliases");
        for (j = 0; j < option->numAliases; j++) {
            clo_writeXmlTag(fout, 4, "alias", option->aliases[j]);
        }
        clo_writeXmlEndTag(fout, 3, "aliases");
    }
    writeXmlDescription(fout, 3, option->desc);
    char str[50];
    sprintf(str, "%d", option->position);
    clo_writeXmlTag(fout, 3, "position", str);

    clo_writeXmlEndTag(fout, 2, "option");
}

/**
 * Write an XML file representation for this option list
 *
 * @param list option list to write out
 * @param filename file to write to
 * @return 0 if OK else -1 is returned
 */
int clo_writeXmlFile(clo_optionList_t *list, const char *filename) {
    clo_option_t* option;
    int i;
    FILE *fout;
    char *key;

    assert(list);
    assert(filename);

    // open file
    fout = fopen(filename, "w");
    if (fout == NULL) {
        fprintf(stderr,
                "clo_writeXmlFile - Could not open \"%s\" for writing.\n",
                filename);
        return -1;
    }

    fprintf(fout, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
    fprintf(fout,
            "<paramInfo xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
    fprintf(fout,
            "            xsi:noNamespaceSchemaLocation=\"http://seadas.gsfc.nasa.gov/software/schemas/ParamInfo-1.0.xsd\">\n");

    clo_writeXmlStartTag(fout, 1, "programMetaData");
    clo_writeXmlTag(fout, 2, "hasParFile", "true");
    clo_writeXmlTag(fout, 2, "parFileOptionName", "par");
    for (i = 0; i < clo_programMetadataList.numEntries; i++) {
        clo_writeXmlTag(fout, 2, clo_programMetadataList.entries[i]->tag,
                clo_programMetadataList.entries[i]->value);
    }
    clo_writeXmlEndTag(fout, 1, "programMetaData");
    clo_writeXmlStartTag(fout, 1, "options");
    if (selectOptionKeys) {
        i = 0;
        key = selectOptionKeys[i];
        while (key) {
            option = clo_findOption(list, key);
            if (option) {
                writeXmlOption(fout, option);
            } else {
                fprintf(stderr,
                        "clo_writeXmlFile - Could not find option \"%s\" in option list.\n",
                        key);
            }
            i++;
            key = selectOptionKeys[i];
        }
    } else {
        for (i = 0; i < list->numOptions; i++) {
            option = list->options[i];
            writeXmlOption(fout, option);
        }
    }

    // write out the positional command line params
    for (i = 0; i < list->positionNumOptions; i++) {
        option = list->positionOptions[i];
        writeXmlOption(fout, option);
    }

    clo_writeXmlEndTag(fout, 1, "options");
    clo_writeXmlEndTag(fout, 0, "paramInfo");
    fclose(fout);
    return 0;
}

/* -------------------------------------------------------------- */
//void clo_resetFileRecursion()
