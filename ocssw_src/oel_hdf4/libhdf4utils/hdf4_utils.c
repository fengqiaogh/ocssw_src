#include "hdf4utils.h"

#include <hdf.h>
#include <mfhdf.h>

/*----------------------------------------------------------------------------*/

int load_att(int32_t obj_id, int32_t index, att_struct *att) {
    int32_t status = FAIL;

    /* initialize */
    att->data = NULL;
    att->id = obj_id;
    att->index = index;
    status = SDattrinfo(att->id,
            att->index,
            att->name,
            &att->ntype,
            &att->nvals);
    if (status) return status;

    /* read data */
    if (att->nvals > 0) {
        att->data = malloc(att->nvals * DFKNTsize(att->ntype));
        status = SDreadattr(att->id, att->index, att->data);
    }

    return status;
}

int load_att_byname(int32_t obj_id, const char *attname, att_struct *att) {
    int32_t status;
    status = load_att(obj_id, SDfindattr(obj_id, attname), att);
    return status;
}

void print_att_info(const att_struct att) {
    if (att.data) {
        printf("%8s", "");
        //printf("%2d  ",att.index);
        printf("%-22s  ", att.name);
        printf("%10s [%d]", hdf_typename(att.ntype), att.nvals);
        if (att.nvals == 1)
            printf(" = %s", fmt_hdf_val(att.data, 0, att.ntype));
        //printf("data = %p",att.data);
        printf("\n");
    }
}

void print_att_vals(const att_struct att) {
    int32_t i;
    printf("%-22s = ", att.name);
    if (att.data) {

        /* string */
        if (att.ntype == DFNT_CHAR8)
            printf("\"%s\"", (char *) att.data);

            /* scalar number */
        else
            if (att.nvals == 1)
            printf("%s", fmt_hdf_val(att.data, 0, att.ntype));

            /* array of numbers */
        else {
            printf("[");
            for (i = 0; i < att.nvals; i++) {
                printf("%s", fmt_hdf_val(att.data, i, att.ntype));
                printf("%s", (i < att.nvals - 1) ? ", " : "]");
            }
        }
        printf("\n");
    }
}

int init_attlist(sds_struct *sds) {
    int32_t i;
    int32_t status = FAIL;

    /* allocate array to hold all attribute info */
    sds->atts = malloc(sds->natts * sizeof (att_struct));

    /* load attribute info */
    for (i = 0; i < sds->natts; i++)
        status = load_att(sds->id, i, &sds->atts[i]);

    return status;
}

/*----------------------------------------------------------------------------*/

int init_sds_byname(int32_t fileid, const char *sdsname, sds_struct *sds) {
    int32_t status = FAIL;
    strcpy(sds->name, sdsname);
    sds->data = NULL;
    sds->id = SDselect(fileid, SDnametoindex(fileid, sdsname));
    status = SDgetinfo(sds->id, NULL,
            &sds->ndims, sds->dimlen,
            &sds->ntype, &sds->natts);
    return status;
}

void print_sds_info(const sds_struct s) {
    int32_t i;
    printf("%4s", "");
    //printf("%8d",s.id);
    printf("%-22s", s.name);
    if (s.id != FAIL) {
        printf("%10s [", hdf_typename(s.ntype));
        for (i = 0; i < s.ndims; i++) {
            printf("%d", s.dimlen[i]);
            printf("%s", i < s.ndims - 1 ? "," : "]");
        }
        //printf("%2d attributes.", s.natts);
    }
    printf("\n");
}

int init_sdslist(int32_t fileid, const char *sdsnames[], sds_struct *sdslist) {
    int32_t n_sds, i = 0;

    /* allocate array to hold all SDS info */
    while (sdsnames[i]) i++;
    n_sds = i;
    sdslist = malloc(n_sds * sizeof (sds_struct));

    /* load SDS info */
    for (i = 0; i < n_sds; i++)
        init_sds_byname(fileid, sdsnames[i], &sdslist[i]);

    /* print SDS info */
    for (i = 0; i < n_sds; i++)
        print_sds_info(sdslist[i]);

    return n_sds;
}

int readall_sds(sds_struct *sds) {
    int32_t i;
    int32_t status = FAIL;
    int32_t *start, *edges;
    if (sds->id != FAIL) {
        start = malloc(sds->ndims * sizeof (int32_t));
        edges = malloc(sds->ndims * sizeof (int32_t));
        sds->nvals = 1;
        for (i = 0; i < sds->ndims; i++) {
            start[i] = 0;
            edges[i] = sds->dimlen[i];
            sds->nvals *= sds->dimlen[i];
        }
        sds->data = malloc(sds->nvals * DFKNTsize(sds->ntype));
        status = SDreaddata(sds->id, start, NULL, edges, sds->data);
        free(start);
        free(edges);
    }
    return status;
}

/*----------------------------------------------------------------------------*/

void free_sds(sds_struct *sds) {
    int32_t i;

    /* Free space allocated for SDS data */
    if (sds->data) {
        free(sds->data);
        sds->data = NULL;
    }

    /* Free space allocated for each populated SDS attribute*/
    if (sds->atts) {
        for (i = 0; i < sds->natts; i++)
            if (sds->atts[i].data) {
                free(sds->atts[i].data);
                sds->atts[i].data = NULL;
            }
        free(sds->atts);
        sds->atts = NULL;
    }

    /* Terminate access to the data set */
    SDendaccess(sds->id);
    sds->id = FAIL;

}

/*----------------------------------------------------------------------------*/

const char *hdf_typename(const int32_t ntype) {
    /* Named initializers automatically size the array as needed */
    static const char* typenames[] = {
        [DFNT_UCHAR8] = "uchar8",
        [DFNT_CHAR8] = "char8",
        [DFNT_FLOAT32] = "float32",
        [DFNT_FLOAT64] = "float64",
        [DFNT_INT8] = "int8",
        [DFNT_UINT8] = "uint8",
        [DFNT_INT16] = "int16",
        [DFNT_UINT16] = "uint16",
        [DFNT_INT32] = "int32",
        [DFNT_UINT32] = "uint32",
        [DFNT_INT64] = "int64",
        [DFNT_UINT64] = "uint64"
    };
    static const int32_t maxindex = sizeof (typenames) / sizeof (typenames[0]);
    if ((ntype < 0) || (maxindex <= ntype))
        return "unknown";
    return typenames[ntype];
}

char *fmt_hdf_val(const void *array, const int32_t i, const int32_t ntype) {
    char *val = NULL;
    val = malloc(40 * sizeof (char));
    memset(val, '\0', 40);

    if (array)
        switch (ntype) {
        case DFNT_UCHAR8: /* == 3 */
            sprintf(val, "%d", ((uchar8 *) array)[i]);
            break;
        case DFNT_CHAR8: /* == 4 */
            sprintf(val, "%c", ((char8 *) array)[i]);
            break;
        case DFNT_FLOAT32: /* == 5 */
            sprintf(val, "%f", ((float *) array)[i]);
            break;
        case DFNT_FLOAT64: /* == 6 */
            sprintf(val, "%f", ((double *) array)[i]);
            break;
        case DFNT_INT8: /* == 20 */
            sprintf(val, "%d", ((int8_t *) array)[i]);
            break;
        case DFNT_UINT8: /* == 21 */
            sprintf(val, "%u", ((uint8_t *) array)[i]);
            break;
        case DFNT_INT16: /* == 22 */
            sprintf(val, "%d", ((int16_t *) array)[i]);
            break;
        case DFNT_UINT16: /* == 23 */
            sprintf(val, "%u", ((uint16_t *) array)[i]);
            break;
        case DFNT_INT32: /* == 24 */
            sprintf(val, "%ld", (long) ((int32_t *) array)[i]);
            break;
        case DFNT_UINT32: /* == 25 */
            sprintf(val, "%lu", (unsigned long) ((uint32_t *) array)[i]);
            break;
        case DFNT_INT64: /* == 26 */
            sprintf(val, "%lld", ((long long *) array)[i]);
            break;
        case DFNT_UINT64: /* == 27 */
            sprintf(val, "%llu", ((unsigned long long *) array)[i]);
            break;
        default:
            fprintf(stderr, "Unknown type code: %d\n", ntype);
            break;
        } /* end switch */
    return (val);
}

void fopen_warn(const char *filename, const char *file, const int32_t line) {
    fprintf(stderr, "-W- %s:%d: Cannot open file: %s\n",
            file, line, filename);
}

void fopen_err(const char *filename, const char *file, const int32_t line) {
    fprintf(stderr, "-E- %s:%d: Cannot open file: %s\n",
            file, line, filename);
    exit(1);
}

/*----------------------------------------------------------------------------*/

int parse_odl(const char *odltext, const char *object, char *value) {
    char *p = NULL; /* pointer to current character */
    char *q = NULL; /* pointer to end of list */
    char endstring[4] = ""; /* marks end of list */

    /* advance pointer to first instance of specified string */
    /*   (will trip on ambiguous text)   */
    if ((p = strstr(odltext, object)) == NULL)
        return 0;

    /* advance pointer to next instance of "VALUE" */
    if ((p = strstr(p, "VALUE")) == NULL)
        return 0;

    /* advance pointer past next instance of "=" */
    if ((p = strstr(p, "=")) == NULL)
        return 0;
    p++;

    /* advance pointer to next non-blank character */
    while (*p == ' ') p++;

    /* if it's an open paren, look for close paren to end value list */
    /* - otherwise use newline. */
    strcpy(endstring, (*p == '(') ? ")" : "\n");

    /* find end of value list */
    if ((q = strstr(p, endstring)) == NULL)
        return 0;

    /* advance pointer to next character not in [("] */
    while ((*p == '(') || (*p == '"')) p++;

    /* copy each remaining character not in [")] */
    /* put on a single line */
    while ((*p) && (p < q)) {
        if ((*p != '"') && (*p != '\n')) {
            *value = *p;
            value++;
        }
        if (*p == ' ') {
            while (*p == ' ') p++;
        }/* skip multiple blanks */
        else p++;
    }

    /* terminate the string */
    *value = 0;

    return strlen(value);
}

int get_hdfeos_meta(int32_t sd_id, char *attribute, char *name, char *data) {
    char eosmeta[EOSMETALEN] = "";
    if (SDreadattr(sd_id, SDfindattr(sd_id, attribute), (VOIDP) eosmeta) != 0) {
        printf("Error reading attribute %s\n", attribute);
        return -1;
    }
    return (parse_odl(eosmeta, name, data) == 0) ? -1 : 0;
}
