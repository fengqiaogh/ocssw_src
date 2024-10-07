#include <stdio.h>
#include "hdf.h"
#include "seaproto.h"

int32_t get_beg_ext32(int32_t n_bins_write, int32_t *binnum_data,
        int32_t *basebin, int32_t nrows, int32_t *beg, int32_t *ext) {
    int i, j;
    int32_t last_row;
    int32_t *row;
    int32_t *uniq;

    row = (int32_t *) calloc(n_bins_write, sizeof (int32_t));
    uniq = (int32_t *) calloc(nrows + 1, sizeof (int32_t));

    /* Determine row number of each bin */
    /* -------------------------------- */
    last_row = 0;
    for (j = 0; j < n_bins_write; j++) {
        for (i = last_row; i <= nrows; i++) {
            if (basebin[i] > binnum_data[j]) {
                row[j] = i;
                last_row = i;
                break;
            }
        }
        if (row[j] == 0) {
            printf("Zero row at bin: %d of %d total bins.\n", j, n_bins_write);
            printf("basebin: %ld  binnum_data: %ld\n", (long int) basebin[i],
                    (long int) binnum_data[j]);
            printf("%d %d\n", i, last_row);
        }
    }

    for (i = 0; i < n_bins_write - 1; i++) {
        if (row[i] > row[i + 1]) {
            printf("Improper row number: %d at element: %d\n", row[i], i);
        }
    }

    /* Determine unique row numbers */
    /* ---------------------------- */
    last_row = row[0];
    uniq[0] = 0;
    i = 1;
    for (j = 1; j < n_bins_write; j++) {
        if (row[j] != last_row) {
            uniq[i++] = j;
            last_row = row[j];
        }
    }
    for (j = i; j < nrows; j++)
        uniq[j] = -1;

    /* Determine begining pixel number for each row */
    /* -------------------------------------------- */
    for (j = 0; j < nrows; j++) {
        if (uniq[j] == -1)
            break;
        beg[row[uniq[j]] - 1] = binnum_data[uniq[j]];

        if (row[uniq[j]] - 1 < 0)
            printf("bad row index: %d %d %d\n", j, uniq[j], row[uniq[j]]);
    }

    uniq[i] = n_bins_write;

    /* Determine extent of each row */
    /* ---------------------------- */
    for (j = 0; j < i; j++) {
        if (uniq[j] == -1)
            break;

        ext[row[uniq[j]] - 1] = uniq[j + 1] - uniq[j];
    }

    free(row);
    free(uniq);

    return 0;
}


int64_t get_beg_ext(int32_t n_bins_write, int64_t *binnum_data,
        int64_t *basebin, int32_t nrows, int64_t *beg, int32_t *ext) {
    int i, j;
    int32_t last_row;
    int32_t *row;
    int32_t *uniq;

    row = (int32_t *) calloc(n_bins_write, sizeof (int32_t));
    uniq = (int32_t *) calloc(nrows + 1, sizeof (int32_t));

    /* Determine row number of each bin */
    /* -------------------------------- */
    last_row = 0;
    for (j = 0; j < n_bins_write; j++) {
        for (i = last_row; i <= nrows; i++) {
            if (basebin[i] > binnum_data[j]) {
                row[j] = i;
                last_row = i;
                break;
            }
        }
        if (row[j] == 0) {
            printf("Zero row at bin: %d of %d total bins.\n", j, n_bins_write);
            printf("basebin: %ld  binnum_data: %ld\n", (long int) basebin[i],
                    (long int) binnum_data[j]);
            printf("%d %d\n", i, last_row);
        }
    }

    for (i = 0; i < n_bins_write - 1; i++) {
        if (row[i] > row[i + 1]) {
            printf("Improper row number: %d at element: %d\n", row[i], i);
        }
    }

    /* Determine unique row numbers */
    /* ---------------------------- */
    last_row = row[0];
    uniq[0] = 0;
    i = 1;
    for (j = 1; j < n_bins_write; j++) {
        if (row[j] != last_row) {
            uniq[i++] = j;
            last_row = row[j];
        }
    }
    for (j = i; j < nrows; j++)
        uniq[j] = -1;

    /* Determine begining pixel number for each row */
    /* -------------------------------------------- */
    for (j = 0; j < nrows; j++) {
        if (uniq[j] == -1)
            break;
        beg[row[uniq[j]] - 1] = binnum_data[uniq[j]];

        if (row[uniq[j]] - 1 < 0)
            printf("bad row index: %d %d %d\n", j, uniq[j], row[uniq[j]]);
    }

    uniq[i] = n_bins_write;

    /* Determine extent of each row */
    /* ---------------------------- */
    for (j = 0; j < i; j++) {
        if (uniq[j] == -1)
            break;

        ext[row[uniq[j]] - 1] = uniq[j + 1] - uniq[j];
    }

    free(row);
    free(uniq);

    return 0;
}
