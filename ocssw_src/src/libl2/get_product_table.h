#ifndef GET_PRODUCT_TABLE
#define GET_PRODUCT_TABLE

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif


/**\file
 * functions used to deal with the product table file.
 */

/**
   Product table data structure.
 */
typedef struct product_table_str {
    char *description; /**< description of the product  */
    char *name; /**< short name of the product */
    char *units; /**< product units */
    char *scaling; /**< scaling type (linear,logarithmic) */
    char *palette; /**< palette file name (without .pal extention) */
    float min; /**< minimum value for scaling */
    float max; /**< max value for scaling */
    char *precision; /**< precision ("B", "I", "F") */
} product_table_t;


/**
 * read the product table file.
 * This routine allocates enough memory to hold the product table and returns
 * a pointer to an array of product_table_t structures to the caller.  It is 
 * the callers responsibility to free the memory.  If there was an error reading
 * the file NULL is returned and *num_entries is set to 0.
 *
 * @param file_name full file name of the product table file
 * @param num_entries retuns the numbers of entries read from the file
 * @return pointer to product table array, NULL if error
 */
product_table_t* get_product_table(char *file_name, int32_t *num_entries);

/**
 * free the product table memory.
 * 
 * @param table pointer to product table array
 * @param num_entries number of entries in the product table
 */
void free_product_table(product_table_t* table, int32_t num_entries);


/** 
 * find the product table entry for name.
 * 
 * @param table pointer to product table array
 * @param num_entries number of entries in the product table
 * @param name product entry to search for
 * @return index for the entry, or -1 if not found
 */
int32_t search_product_table(product_table_t* table, int32_t num_entries, char* name);



#ifdef __cplusplus
}
#endif

#endif
