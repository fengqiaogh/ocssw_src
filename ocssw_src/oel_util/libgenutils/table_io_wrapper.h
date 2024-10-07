#ifndef  _TABLE_IO_WRAPPER_H
#define  _TABLE_IO_WRAPPER_H

int table_column_count(char *input_filename);
int table_row_count(char *input_filename);
int *table_read_i4(char *input_filename, int m, int n);
double *table_read_r8(char *input_filename, int m, int n);
float *table_read_r4(char *input_filename, int m, int n);
void table_free_i4(int *table);
void table_free_r8(double *table);
void table_free_r4(float *table);

#endif
