#ifndef _PASSTHEBUCK_H_
#define _PASSTHEBUCK_H_

#define LIFE_IS_GOOD  0 /* Successful return status */
#define FILE_ALREADY_EXISTS 1 /* HDF output file already existed */
#define MEMORY_ALLOCATION_ERROR 2 /* malloc() or calloc() failed */
#define HDF_FUNCTION_ERROR 3 /* One of NCSA's functions complained */
#define PROGRAMMER_BOOBOO 4 /* Error in code logic or consistency */
#define CALDATA_NOT_APPENDED 5 /* Calibration data was not appended */
#define MISSING_ENVVAR  6 /* Environment variable not set */
#define STAT_FAILURE  7 /* Call to stat() function failed */
#define FILE_IS_EMPTY  8 /* File size is zero bytes */
#define FOPEN_FAILURE  9 /* Call to fopen() function failed */
#define FREAD_FAILURE  10 /* Call to fread() function failed */

#define PTB(function){   \
  int status = function;  \
  switch(status){   \
    case LIFE_IS_GOOD: break;  \
    default:  return(status); \
  }     \
}

#define DPTB(function){   \
  int status = function;  \
  switch(status){   \
    case LIFE_IS_GOOD: break;  \
    default:  printf("Programming error: file %s, line %d\n",__FILE__,__LINE__); exit(1); \
  }     \
}

#endif
