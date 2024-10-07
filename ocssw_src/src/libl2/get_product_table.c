#include <get_product_table.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>



#define MAX_LINE_LENGTH 1024

product_table_t* get_product_table(char *file_name, int32_t *num_entries) {
    FILE *fp;
    char buf[MAX_LINE_LENGTH];
    product_table_t* table;
    int32_t i;
    char* cptr;

    *num_entries = 0;

    /* open product table file */
    fp = fopen(file_name, "r");
    if (fp == 0x0) {
        fprintf(stderr, "get_product_table - product table \"%s\" not found.\n", file_name);
        return NULL;
    }

    // loop through the file and count the nember of entries
    while (fgets(buf, MAX_LINE_LENGTH, fp) != NULL) {
        if ((buf[0] >= 0x41) && (buf[0] <= 0x5a))
            (*num_entries)++;
    }
    fseek(fp, 0, SEEK_SET);

    if (*num_entries == 0) {
        fprintf(stderr, "get_product_table - no entries found in file \"%s\"\n", file_name);
        return NULL;
    }

    table = (product_table_t*) malloc(sizeof (product_table_t) * *num_entries);
    if (table == NULL) {
        fprintf(stderr, "get_product_table - could not allocate memory.\n");
        return NULL;
    }

    i = 0;
    while (fgets(buf, MAX_LINE_LENGTH, fp) != NULL) {
        if ((buf[0] >= 0x41) && (buf[0] <= 0x5a)) {

            cptr = strtok(buf, ":");
            table[i].description = strdup(cptr);

            cptr = strtok(NULL, ":");
            table[i].name = strdup(cptr);

            cptr = strtok(NULL, ":");
            table[i].units = strdup(cptr);

            cptr = strtok(NULL, ":");
            table[i].scaling = strdup(cptr);

            cptr = strtok(NULL, ":");
            table[i].min = (float) atof(cptr);

            cptr = strtok(NULL, ":");
            table[i].max = (float) atof(cptr);

            cptr = strtok(NULL, ":");
            table[i].precision = strdup(cptr);

            cptr = strtok(NULL, "\n");
            table[i].palette = strdup(cptr);

            i++;
        } // if buff[0] is capital letter
    } // while more lines
    fclose(fp);

    if (i != *num_entries) {
        fprintf(stderr, "get_product_table - number_entries does not equal lines read.\n");
        *num_entries = 0;
        return NULL;
    }


    return table;
}

void free_product_table(product_table_t* table, int32_t num_entries) {
    int32_t i;

    if (table) {
        for (i = 0; i < num_entries; i++) {
            free(table[i].description);
            free(table[i].name);
            free(table[i].units);
            free(table[i].scaling);
            free(table[i].precision);
            free(table[i].palette);
        }
        free(table);
    }
}

int32_t search_product_table(product_table_t* table, int32_t num_entries, char* name) {
    int32_t i;

    if (table == NULL)
        return -1;

    if (name == NULL)
        return -1;

    for (i = 0; i < num_entries; i++) {
        if (table[i].name) {
            if (strcmp(table[i].name, name) == 0) {
                return i;
            }
        }
    }

    return -1;
}


