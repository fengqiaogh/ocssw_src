#include "bin_util.h"
#include "L3Shape.h"

#include <mfhdf.h>

using namespace std;

namespace Hdf {

int create_vdata(int32_t fileid, int32_t vgid, int32_t *vdata_id,
        const char *vdata_name, const char *class_name, int32_t n_flds,
        char const * const fldname[], int32_t type[], int32_t noext, int32_t *aid) {
    int32_t status;
    char buffer[1000];
    static int32_t prod_count = 0;

    // Setup new vdata
    // ---------------
    *vdata_id = VSattach(fileid, -1, "w");
    //    printf("Setting Up: %s %ld\n", vdata_name, *vdata_id);

    for (int32_t i = 0; i < n_flds; i++) {
        status = VSfdefine(*vdata_id, fldname[i], type[i], 1);
        if (status != 0) {
            printf("Error defining \"%s\".\n", fldname[i]);
        }
    }

    /* Set fieldnames */
    /* -------------- */
    strcpy(buffer, fldname[0]);
    for (int32_t i = 1; i < n_flds; i++) {
        strcat(buffer, ",");
        strcat(buffer, fldname[i]);
    }
    status = VSsetfields(*vdata_id, buffer);

    status = Vinsert(vgid, *vdata_id);
    status = VSsetname(*vdata_id, vdata_name);
    status = VSsetclass(*vdata_id, class_name);

    VSsetblocksize(*vdata_id, 4096 * 6);

    // Setup External Data Vdatas
    // --------------------------
    if (strcmp(class_name, "DataSubordinate") == 0 && noext == 0) {

        static FILE *sfile;

        memset(buffer, 0, sizeof (buffer));
        status = VSwrite(*vdata_id, (uint8_t *) buffer, 1, FULL_INTERLACE);
        VSdetach(*vdata_id);
        *vdata_id = VSattach(fileid, VSfind(fileid, vdata_name), "w");

        //  synthesize file name
        char *outname;
        intn access, attach;
        Hfidinquire(fileid, &outname, &access, &attach);
        strcpy(buffer, outname);

        // Remove trailing ".main" if it exists
        char *cptr = strstr(buffer, ".main");
        if (strlen(buffer) + buffer - cptr == 5)
            *cptr = 0;

        sprintf(&buffer[strlen(buffer)], ".x%02d", prod_count);

        //  open file and write name of main file
        //  in first 512 bytes
        sfile = fopen(buffer, "w");

        if (strrchr(outname, '/') != NULL)
            strcpy(buffer, strrchr(outname, '/') + 1);
        else
            strcpy(buffer, outname);

        cptr = strstr(buffer, ".main");
        if (strlen(buffer) + buffer - cptr == 5)
            *cptr = 0;

        status = fwrite(buffer, 512, 1, sfile);
        fclose(sfile);

        //  convert to external element
        strcpy(buffer, outname);

        cptr = strstr(buffer, ".main");
        if (strlen(buffer) + buffer - cptr == 5)
            *cptr = 0;

        sprintf(&buffer[strlen(buffer)], ".x%02d", prod_count);
        *aid = HXcreate(fileid, DFTAG_VS, (uint16_t) VSfind(fileid, vdata_name),
                buffer, 512, 0);

        VSdetach(*vdata_id);
        *vdata_id = VSattach(fileid, VSfind(fileid, vdata_name), "w");
        VSseek(*vdata_id, 0);

        prod_count++;
    }

    return 0;
}

int32_t write_vdata(int vdata_id, int32_t n_recs_to_write, void *data) {
    int32_t n_write = VSwrite(vdata_id, (const uint8_t*) data,
            n_recs_to_write, FULL_INTERLACE);

    if (n_write != n_recs_to_write) {
        cout << "Error in VSwrite (write_vdata): "
                << n_recs_to_write << " " << n_write << endl;
        exit(-1);
    }

    return n_write;
}

int read_binList(int n_elem, int32_t vdata_id_binlist, binListStruct* binList, l3::L3Shape* shape) {
    if (n_elem == 0) {
        printf("The number of elements requested to be read is 0.\n");
    }

    // Offset within BinList memory struct
    int32_t binlist_offset[] = {0, 4, 6, 8, 12, 16, 20};

    // Size of fields in BinList Vdata
    int32_t binlist_size[] = {4, 2, 2, 2, 4, 1, 4};

    uint8_t *bytebuf;
    bytebuf = (uint8_t *) malloc(n_elem * 19);

    int records_read;

    records_read = VSread(vdata_id_binlist, bytebuf, n_elem, FULL_INTERLACE);

    // Loop through bins
    for (int32_t j = 0; j < n_elem; j++) {
        int binlist_idx = j * 19;
        for (size_t i = 0; i < 7; i++) {

            // Point to output location in BinList structure (memory)
            char *ptr = (char *) &binList[j] + binlist_offset[i];
            memcpy(ptr, &bytebuf[binlist_idx], binlist_size[i]);

            binlist_idx += binlist_size[i];

            if (i == 0) {
                int32_t bin;
                memcpy(&bin, &bytebuf[j * 19], 4);
                float lat;
                float lon;
                shape->bin2latlon(bin, lat, lon);
                memcpy((char *) &binList[j] + 24, &lat, 4);
                memcpy((char *) &binList[j] + 28, &lon, 4);
            }
        }
    }
    free(bytebuf);

    return records_read;
} // end read_binList

int write_binList(int n_elem, int32_t vdata_id_binlist, binListStruct* binList) {
    if (n_elem == 0) {
        printf("No elements are requested to be written in \"write_binList\".\n");
    }

    // Offset within BinList memory struct
    int32_t binlist_struct_offset[] = {0, 4, 6, 8, 12, 16, 20};

    // Offset within BinList Vdata
    int32_t binlist_vdata_offset[] = {0, 4, 6, 8, 10, 14, 15};

    // Size of fields in BinList Vdata
    int32_t binlist_size[] = {4, 2, 2, 2, 4, 1, 4};

    uint8_t *bytebuf;
    bytebuf = (uint8_t *) malloc(n_elem * 19);

    for (int32_t j = 0; j < n_elem; j++) {
        // Loop through BinList fields
        for (size_t i = 0; i < 7; i++) {
            char *ptr = (char *) &binList[j] + binlist_struct_offset[i];
            memcpy(&bytebuf[j * 19 + binlist_vdata_offset[i]], ptr,
                    binlist_size[i]);
        }
    }
    int records_written =
            VSwrite(vdata_id_binlist, bytebuf, n_elem, FULL_INTERLACE);

    if (records_written != n_elem) {
        cout << "Error in VSwrite (write_binList): "
                << n_elem << " " << records_written << endl;
        exit(-1);
    }

    free(bytebuf);

    return records_written;
} // end write_binList

int write_prodData(int n_elem, int32_t vdata_id_proddata, float *data,
        binListStruct* binList) {
    float dta[2];
    binListStruct binListRec;
    int records_written;

    for (int32_t j = 0; j < n_elem; j++) {
        memcpy((void *) &binListRec, &binList[j], sizeof (binListStruct));
        //      cout << binListRec.weights << endl;
        dta[0] = data[j] * binListRec.weights;

        // Can't really compute sum_2 for a single data value per bin.
        //      dta[1] = data[j] * data[j] * binListRec.weights;
        dta[1] = 0.0;
        records_written = VSwrite(vdata_id_proddata, (const uint8_t*) dta, 1,
                FULL_INTERLACE);

        if (records_written != 1) {
            cout << "Error in VSwrite (write_prodData): "
                    << 1 << " " << records_written << endl;
            exit(-1);
        }

    }
    return 0;
}

int copy_prodData(int n_elem, int32_t *binsToCopy,
        char const * const fldname3[],
        int32_t in_vdata_id_proddata,
        int32_t out_vdata_id_proddata) {
    float dta[2];
    char buffer[1000];

    strcpy(buffer, fldname3[0]);
    strcat(buffer, ",");
    strcat(buffer, fldname3[1]);
    intn status = VSsetfields(in_vdata_id_proddata, buffer);
    if (status != 0) {
        printf("Error setting fields: %s\n", buffer);
        exit(1);
    }

    for (int32_t j = 0; j < n_elem; j++) {
        VSseek(in_vdata_id_proddata, binsToCopy[j]);
        int records_read = VSread(in_vdata_id_proddata,
                (uint8_t *) & dta[0], 1, FULL_INTERLACE);

        int records_written = VSwrite(out_vdata_id_proddata,
                (const uint8_t*) dta, 1, FULL_INTERLACE);

        if (records_read != records_written) {
            cout << "Error in copy: "
                    << records_read << " " << records_written << endl;
            exit(-1);
        }
    }
    return 0;
}

int create_compound(hid_t group_id, const char *dataset_name, hid_t *dataset_id,
        hid_t *type_id, size_t typesize, int32_t n_flds,
        char const * const fldname[], size_t offset[], hid_t type[],
        hid_t *filespace, hid_t dataspace) {
    herr_t status;

    *type_id = H5Tcreate(H5T_COMPOUND, typesize);
    for (int32_t i = 0; i < n_flds; i++) {
        status = H5Tinsert(*type_id, fldname[i], offset[i], type[i]);
        if (status < 0) {
            cout << "H5Tinsert error for: " << dataset_name
                    << " (field: " << fldname[i] << ")." << endl;
            exit(1);
        }
    }

    hid_t cparms = H5Pcreate(H5P_DATASET_CREATE);

    hsize_t chuck_dims[1] = {1};
    status = H5Pset_chunk(cparms, 1, chuck_dims);
    *dataset_id = H5Dcreate1(group_id, dataset_name, *type_id, dataspace,
            cparms);
    if (*dataset_id < 0) {
        cout << "Cannot create dataset for: " << dataset_name << "." << endl;
        exit(1);
    }

    status = H5Pclose(cparms);

    *filespace = H5Dget_space(*dataset_id);
    if (*filespace < 0) {
        cout << "Cannot get filespace for: " << dataset_name << "." << endl;
        exit(1);
    }

    cout << "\"" << dataset_name << "\"" << " successfully created." << endl;
    cout << "Dataset id:   " << *dataset_id << endl;
    cout << "Datatype id:  " << *type_id << endl;
    cout << "Filespace id: " << *filespace << endl << endl;
    return 0;
}
}
