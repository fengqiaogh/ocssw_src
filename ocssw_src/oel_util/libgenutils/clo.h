#ifndef CLO_H
#define CLO_H

#include <stdio.h>
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif


/**\file 
   Command Line Option (CLO) is a library for reading and documenting
   the command line options for a program.  The library can also read
   parameter files.

   First the program should set the version string using
   #clo_setVersionString() function.  Then the option list is created using
   #clo_createList().  The list is filled by defining options and assigning 
   default values and descriptions using the #clo_addOption() function.  

   Next the program should call #clo_readArgs() which first reads the 
   command line arguments into the options list and looks for pre-defined
   options like -help and -version.

   The last occurrence of an option is the one used.  Overwriting the previous
   option value.

   Options on the command line are expressed as "key" or "key=val".  Each 
   option is separated by white space and multiple values are separated by 
   commas, colons, or spaces.

   progName key10 key=val key1=val1,val2,val3 key2="val1 val2 val3"
   
   boolean<br>
   key=0,1,y,n,yes,no,on,off

   int<br>
   key=1,2,548,9<br>
   key=[1,2,3,4]<br>
   key=(1,2,3,4)<br>
   key="[1,2,3,4]"<br>
   key=["1","2",3,4]<br>
   key=1:2<br>

   float<br>
   key=1,8.3,.6, 1.0 , 1<br>

   string<br>
   key=hello,dude<br>
   key="two words","other",what<br>
   key="val1,val2,val3"<br>
   key="val\,with comma, other"<br>


   There are a few key values that are inherently understood by the 
   library.  They are invoked when #clo_readArgs() is called. 

   -version, --version  will print out the version string<br>
   -h, -help, --help  will print the usage<br>
   -dump_options will print out everything help does plus the current 
   value of the option and where the value came from.
   -dump_options_paramfile will write a param file with all of the params and values
   -dump_options_xmlfile will write an xml file with all of the param information

 */

/** how many pointers to allocate space for when adding to an array */
#define CLO_CHUNK_SIZE 32

/** characters that can separate array items */
#define CLO_ARRAY_DELIMITER " \t,:[]()\""

// Define a silly string max length
#define SILLYSTRING 4096

/**
   Data types for the command line arguments.
 */
enum clo_dataType_t {
    CLO_TYPE_BOOL,
    CLO_TYPE_INT,
    CLO_TYPE_INT64,
    CLO_TYPE_FLOAT,
    CLO_TYPE_DOUBLE,
    CLO_TYPE_STRING,
    CLO_TYPE_IFILE,
    CLO_TYPE_OFILE,
    CLO_TYPE_HELP,
    CLO_TYPE_POSITION
};

// forward reference for option structure
struct clo_option_t;

/** 
    Option callback function definition.  This function will get called
    right after this option is read.
    \param option pointer to the option associated with the callback
 */
typedef void (*clo_optionCallback_t)(struct clo_option_t *option);

/**
   Option data structure.
 */
typedef struct clo_option_t {
    char *key; /**< key string  */
    enum clo_dataType_t dataType; /**< option data type (CLO_TYPE_BOOL, CLO_TYPE_INT, CLO_TYPE_FLOAT, CLO_TYPE_DOUBLE, CLO_TYPE_STRING) */
    char *defaultVal; /**< initial value for this option */
    char *desc; /**< description of this option */
    char *valStr; /**< value of option, original string */
    char *source; /**< where did this value come from */
    clo_optionCallback_t cb; /**< function that gets called when this option is read */
    void *cb_data; /**< pointer for user callback data */
    int cb_flag; /**< flag for user callback function*/
    char **strArray; /**< array of strings */
    void *valArray; /**< array of values of the correct data type */
    int count; /**< how many items in the array */
    char **aliases; /**< list of aliases for key */
    int numAliases; /**< number of aliases */
    int aliasStorageSize; /**< number of slots allocated in alias array */
    int position; /**< position of the parameter (CLO_TYPE_POSITION only) */
} clo_option_t;

/** structure for a list of options */
typedef struct clo_optionList_t {
    int storageSize; /**< number of slots allocated in the options array */
    int numOptions; /**< number of options stored in option array */
    clo_option_t **options; /**< array of option pointers */

    int positionStorageSize; /**< number of slots allocated in the options array */
    int positionNumOptions; /**< number of options stored in option array */
    clo_option_t **positionOptions; /**< array of option pointers */
} clo_optionList_t;

/** structure for the extra XML program metadata */
typedef struct clo_programMetadata_t {
    char* tag; /**< XML tag string */
    char* value; /**< XML value string */
} clo_programMetadata_t;

/** structure for a list of XML program metadata */
typedef struct clo_programMetadataList_t {
    int storageSize; /**< number of slots allocated in the options array */
    int numEntries; /**< number of options stored in option array */
    clo_programMetadata_t **entries; /**< array of option pointers */
} clo_programMetadataList_t;

// global library settings
void clo_setIgnoreKeyCase(int val);
int clo_getIgnoreKeyCase();

void clo_setEnableDumpOptions(int val);
int clo_getEnableDumpOptions();

void clo_setEnableExtraOptions(int val);
int clo_getEnableExtraOptions();

void clo_setVersion(const char *str);
void clo_setVersion2(const char* programName, const char* versionStr);
char* clo_getVersion();

void clo_setHelpStr(const char *str);
char* clo_getHelpStr();

void clo_setSelectOptionKeys(char *keys[]);
char **clo_getSelectOptionKeys();

clo_optionList_t* clo_createList();

clo_option_t* clo_createOption(const char *key, enum clo_dataType_t dataType,
        const char *defaultVal, const char *desc);
void clo_addOptionAlias(clo_option_t *option, const char* alias);
void clo_addAlias(clo_optionList_t *list, const char *key, const char* alias);
void clo_insertOption(clo_optionList_t *list, clo_option_t *option);
clo_option_t* clo_addOption(clo_optionList_t *list, const char *key,
        enum clo_dataType_t dataType, const char *defaultVal,
        const char* desc);
clo_option_t* clo_copyOption(clo_option_t *option);
clo_optionList_t* clo_copyList(clo_optionList_t *list);
void clo_deleteOption(clo_option_t *option);
void clo_deleteList(clo_optionList_t *list);
clo_option_t* clo_getOption(clo_optionList_t *list, int i);
clo_option_t* clo_findOption(clo_optionList_t *list, const char *key);
int clo_getNumOptions(clo_optionList_t *list);
char* clo_getOptionRawString(clo_option_t *option);

char* clo_getOptionString(clo_option_t *option);
int clo_getOptionBool(clo_option_t *option);
int clo_getOptionInt(clo_option_t *option);
int64_t clo_getOptionInt64(clo_option_t *option);
float clo_getOptionFloat(clo_option_t *option);
double clo_getOptionDouble(clo_option_t *option);

char** clo_getOptionStrings(clo_option_t *option, int *count);
int* clo_getOptionBools(clo_option_t *option, int *count);
int* clo_getOptionInts(clo_option_t *option, int *count);
float* clo_getOptionFloats(clo_option_t *option, int *count);
double* clo_getOptionDoubles(clo_option_t *option, int *count);

char* clo_getRawString(clo_optionList_t *list, const char *key);
char* clo_getString(clo_optionList_t *list, const char *key);
int clo_getBool(clo_optionList_t *list, const char *key);
int clo_getInt(clo_optionList_t *list, const char *key);
float clo_getFloat(clo_optionList_t *list, const char *key);
double clo_getDouble(clo_optionList_t *list, const char *key);

char** clo_getStrings(clo_optionList_t *list, const char *key, int *count);
int* clo_getBools(clo_optionList_t *list, const char *key, int *count);
int* clo_getInts(clo_optionList_t *list, const char *key, int *count);
float* clo_getFloats(clo_optionList_t *list, const char *key, int *count);
double* clo_getDoubles(clo_optionList_t *list, const char *key, int *count);

int clo_setOptionString(clo_option_t *option, const char *val, const char *source);
int clo_setString(clo_optionList_t *list, const char* key, const char *val, const char *source);

// position option functions
void clo_setEnablePositionOptions(int val);
int clo_getEnablePositionOptions();
int clo_getPositionNumOptions(clo_optionList_t *list);
clo_option_t* clo_getPositionOption(clo_optionList_t *list, int i);
char* clo_getPositionString(clo_optionList_t *list, int pos);

void clo_setEnableParOption(int val);
int clo_getEnableParOption();

void clo_printOptionVal(clo_option_t *option);
void clo_printVals(clo_optionList_t *list);
void clo_printOption(clo_option_t *option);
void clo_printOptions(clo_optionList_t *list);
void clo_dumpOption(clo_option_t *option);
void clo_dumpOptions(clo_optionList_t *list);
void clo_printVersion();
void clo_printHelpString();
void clo_printUsage(clo_optionList_t *list);

void clo_readString(clo_optionList_t *list, const char *str, const char *source);
void clo_readArgs(clo_optionList_t *list, int argc, char *argv[]);
void clo_readArgsPar(clo_optionList_t *list, int argc, char *argv[], int enableFileDescending);
void clo_readOptions(clo_optionList_t *list, clo_optionList_t *readList);
void clo_readFile(clo_optionList_t *list, const char *fileName);

int clo_isOptionSet(clo_option_t *option);
int clo_isSet(clo_optionList_t *list, const char *key);

int clo_writeParameterFile(clo_optionList_t *list, const char *filename);

void clo_addXmlProgramMetadata(const char* tag, const char* value);
void clo_clearXmlProgramMetadata();

void clo_writeXmlStartTag(FILE *fout, int level, const char *tag);
void clo_writeXmlEndTag(FILE *fout, int level, const char *tag);
void clo_writeXmlTag(FILE *fout, int level, const char *tag, const char *value);
int clo_writeXmlFile(clo_optionList_t *list, const char *filename);

void clo_trimDashes(char* str);
char* clo_trimDashesDup(const char* str);

//void clo_resetFileRecursion();


#ifdef __cplusplus
}
#endif

#endif
