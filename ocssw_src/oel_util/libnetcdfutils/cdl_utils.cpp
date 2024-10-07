#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libgen.h>

#include <sstream>
#include <iomanip> 
#include "nc4utils.h"
#include "cdl_utils.h"
#include "netcdf.h"
#include "timeutils.h"

ncdfFile::ncdfFile() {
    ncid = -1;
}

ncdfFile::~ncdfFile() {

}

/*----------------------------------------------------------------- */
/* Create a netcdf4 level1 file from cdl file                       */

/* ---------------------------------------------------------------- */
int ncdfFile::cdlCreate(char* l1_filename, char* cdl_filename,
        int32_t numScans) {

    int status, i;
    //   idDS ds_id;
    status = nc_create(l1_filename, NC_NETCDF4, &ncid);
    check_err(status, __LINE__, __FILE__);

    ifstream cdl_data_structure;
    string line;
    string dataStructureFile(cdl_filename);

    expandEnvVar(&dataStructureFile);

    cdl_data_structure.open(dataStructureFile.c_str(), ifstream::in);
    if (cdl_data_structure.fail() == true) {
        cout << "\"" << dataStructureFile.c_str() << "\" not found" << endl;
        exit(1);
    }

    // Find "dimensions" section of CDL file
    while (1) {
        getline(cdl_data_structure, line);

        size_t firstNonBlank = line.find_first_not_of(" ");
        if (firstNonBlank != string::npos)
            if (line.compare(firstNonBlank, 2, "//") == 0) continue;

        size_t pos = line.find("dimensions:");
        if (pos == 0) break;
    }

    // Define dimensions from "dimensions" section of CDL file
    ndims = 0;
    while (1) {
        getline(cdl_data_structure, line);

        size_t firstNonBlank = line.find_first_not_of(" ");
        if (firstNonBlank != string::npos)
            if (line.compare(firstNonBlank, 2, "//") == 0) continue;

        size_t pos = line.find(" = ");
        if (pos == string::npos) break;

        uint32_t dimSize;
        istringstream iss(line.substr(pos + 2, string::npos));
        iss >> dimSize;

        iss.clear();
        iss.str(line);
        iss >> skipws >> line;

        cout << "Dimension Name: " << line.c_str() << " Dimension Size: "
                << dimSize << endl;

        status = nc_def_dim(ncid, line.c_str(), dimSize, &dimid[ndims++]);
        check_err(status, __LINE__, __FILE__);
    } // while loop

    ngrps = 0;
    // Loop through groups
    while (1) {
        getline(cdl_data_structure, line);

        size_t firstNonBlank = line.find_first_not_of(" ");
        if (firstNonBlank != string::npos)
            if (line.compare(firstNonBlank, 2, "//") == 0) continue;

        // Check if end of CDL file
        // If so then close CDL file and return
        if (line.substr(0, 1).compare("}") == 0) {
            cdl_data_structure.close();
            return 0;
        }

        // Check for beginning of new group
        size_t pos = line.find("group:");

        // If found then create new group and variables
        if (pos == 0) {

            // Parse group name
            istringstream iss(line.substr(6, string::npos));
            iss >> skipws >> line;
            cout << "Group: " << line.c_str() << endl;

            // Create NCDF4 group
            status = nc_def_grp(ncid, line.c_str(), &this->gid[ngrps]);
            check_err(status, __LINE__, __FILE__);

            ngrps++;

            int numDims = 0;
            int varDims[NC_MAX_DIMS];
            size_t dimSize[NC_MAX_DIMS];
            char dimName[NC_MAX_NAME + 1];
            string sname;
            string lname;
            string standard_name;
            string units;
            string flag_values;
            string flag_meanings;
            double valid_min = 0.0;
            double valid_max = 0.0;
            double fill_value = 0.0;
            size_t *chunksize;

            int ntype = 0;

            // Loop through datasets in group
            // Skip until "variables:" found
            while (1) {
                getline(cdl_data_structure, line);
                if (line.find("variables:") != string::npos) break;
            }

            while (1) {
                getline(cdl_data_structure, line);

                if (line.length() == 0) continue;
                if (line.substr(0, 1).compare("\r") == 0) continue;
                if (line.substr(0, 1).compare("\n") == 0) continue;

                size_t firstNonBlank = line.find_first_not_of(" ");
                if (firstNonBlank != string::npos)
                    if (line.compare(firstNonBlank, 2, "//") == 0) continue;

                size_t pos = line.find(":");

                // No ":" found, new dataset or empty line or end-of-group
                if (pos == string::npos) {

                    if (numDims > 0) {
                        // Create previous dataset
                        createNCDF(gid[ngrps - 1],
                                sname.c_str(), lname.c_str(),
                                standard_name.c_str(), units.c_str(),
                                (void *) &fill_value,
                                flag_values.c_str(), flag_meanings.c_str(),
                                valid_min, valid_max, ntype, numDims, varDims, chunksize);

                        flag_values.assign("");
                        flag_meanings.assign("");
                        free(chunksize);
                    }

                    valid_min = 0.0;
                    valid_max = 0.0;
                    fill_value = 0.0;

                    if (line.substr(0, 10).compare("} // group") == 0) break;

                    // Parse variable type
                    string varType;
                    istringstream iss(line);
                    iss >> skipws >> varType;

                    // Get corresponding NC variable type
                    if (varType.compare("char") == 0) ntype = NC_CHAR;
                    else if (varType.compare("byte") == 0) ntype = NC_BYTE;
                    else if (varType.compare("short") == 0) ntype = NC_SHORT;
                    else if (varType.compare("int") == 0) ntype = NC_INT;
                    else if (varType.compare("long") == 0) ntype = NC_INT;
                    else if (varType.compare("float") == 0) ntype = NC_FLOAT;
                    else if (varType.compare("real") == 0) ntype = NC_FLOAT;
                    else if (varType.compare("double") == 0) ntype = NC_DOUBLE;
                    else if (varType.compare("ubyte") == 0) ntype = NC_UBYTE;
                    else if (varType.compare("ushort") == 0) ntype = NC_USHORT;
                    else if (varType.compare("uint") == 0) ntype = NC_UINT;
                    else if (varType.compare("int64") == 0) ntype = NC_INT64;
                    else if (varType.compare("uint64") == 0) ntype = NC_UINT64;

                    // Parse short name (sname)
                    pos = line.find("(");
                    size_t posSname = line.substr(0, pos).rfind(" ");
                    sname.assign(line.substr(posSname + 1, pos - posSname - 1));
                    cout << "sname: " << sname.c_str() << endl;

                    // Parse variable dimension info
                    this->parseDims(line.substr(pos + 1, string::npos),
                            &numDims, varDims);
                    for (int i = 0; i < numDims; i++) {
                        nc_inq_dim(ncid, varDims[i], dimName, &dimSize[i]);
                        cout << line.c_str() << " " << i << " " << dimName
                                << " " << dimSize[i] << endl;
                    }

                    chunksize = (size_t *) calloc(numDims, sizeof (size_t));

                } else {
                    // Parse variable attributes
                    size_t posEql = line.find("=");
                    size_t pos1qte = line.find("\"");
                    size_t pos2qte = line.substr(pos1qte + 1, string::npos).find("\"");
                    //    cout << line.substr(pos+1, posEql-pos-2).c_str() << endl;

                    string attrName = line.substr(pos + 1, posEql - pos - 2);

                    // Get long_name
                    if (attrName.compare("long_name") == 0) {
                        lname.assign(line.substr(pos1qte + 1, pos2qte));
                        //             cout << "lname: " << lname.c_str() << endl;
                    }
                        // Get units
                    else if (attrName.compare("units") == 0) {
                        units.assign(line.substr(pos1qte + 1, pos2qte));
                        //             cout << "units: " << units.c_str() << endl;
                    }
                        // Get _FillValue
                    else if (attrName.compare("_FillValue") == 0) {
                        iss.clear();
                        iss.str(line.substr(posEql + 1, string::npos));
                        iss >> fill_value;
                        //             cout << "_FillValue: " << fill_value << endl;
                    }
                        // Get flag_values
                    else if (attrName.compare("flag_values") == 0) {
                        flag_values.assign(line.substr(pos1qte + 1, pos2qte));
                    }
                        // Get flag_meanings
                    else if (attrName.compare("flag_meanings") == 0) {
                        flag_meanings.assign(line.substr(pos1qte + 1, pos2qte));
                    }
                        // Get valid_min
                    else if (attrName.compare("valid_min") == 0) {
                        iss.clear();
                        iss.str(line.substr(posEql + 1, string::npos));
                        iss >> valid_min;
                        //             cout << "valid_min: " << valid_min << endl;
                    }
                        // Get valid_max
                    else if (attrName.compare("valid_max") == 0) {
                        iss.clear();
                        iss.str(line.substr(posEql + 1, string::npos));
                        iss >> valid_max;
                        //             cout << "valid_max: " << valid_max << endl;
                    }                        // Get Chunk sizes
                    else if (attrName.compare("_ChunkSizes") == 0) {
                        size_t posComma = line.find(",");
                        status = posEql + 1;

                        for (i = 0; i < numDims; i++) {
                            iss.clear();
                            if (i == numDims - 1) {
                                iss.str(line.substr(status, string::npos));
                            } else {
                                iss.str(line.substr(status, posComma));
                                status = posComma + 1;
                                posComma += line.substr(posComma + 1, string::npos).find(",") + 1;
                            }
                            iss >> chunksize[i];
                        }
                        //             cout << "_FillValue: " << fill_value << endl;
                    }

                } // if ( pos == string::npos)
            } // datasets in group loop
        } // New Group loop
    } // Main Group loop


    return 0;
}

/*----------------------------------------------------------------- */
/* Create an Generic NETCDF4 level1 file                            */

/* ---------------------------------------------------------------- */
int ncdfFile::cdlCreateDim(char* l1_filename, char* cdl_filename, const char** dim_names, size_t* dim_size, size_t n_dims,
        size_t numScans) {

    // Let user define the dimensions and sizes
    int status, i;
    //  idDS ds_id;
    status = nc_create(l1_filename, NC_NETCDF4, &ncid);
    check_err(status, __LINE__, __FILE__);

    ifstream cdl_data_structure;
    string line;
    string dataStructureFile(cdl_filename);

    expandEnvVar(&dataStructureFile);

    cdl_data_structure.open(dataStructureFile.c_str(), ifstream::in);
    if (cdl_data_structure.fail() == true) {
        cout << "\"" << dataStructureFile.c_str() << "\" not found" << endl;
        exit(1);
    }


    // Define dimensions from "dimensions" section of CDL file
    ndims = 0;

    for (size_t i = 0; i < n_dims; i++) {

        cout << "Dimension Name: " << dim_names[i] << " Dimension Size: "
                << dim_size[i] << endl;

        status = nc_def_dim(ncid, dim_names[i], dim_size[i], &dimid[ndims++]);
        check_err(status, __LINE__, __FILE__);
    } // for loop

    // Find "dimensions" section of CDL file
    while (1) {
        getline(cdl_data_structure, line);

        size_t firstNonBlank = line.find_first_not_of(" ");
        if (firstNonBlank != string::npos)
            if (line.compare(firstNonBlank, 2, "//") == 0) continue;

        size_t pos = line.find("dimensions:");
        if (pos == 0) break;
    }

    // Then skip over the dimensions section
    while (1) {
        getline(cdl_data_structure, line);

        size_t firstNonBlank = line.find_first_not_of(" ");
        if (firstNonBlank != string::npos)
            if (line.compare(firstNonBlank, 2, "//") == 0) continue;

        size_t pos = line.find(" = ");
        if (pos == string::npos) break;

        uint32_t dimSize;
        istringstream iss(line.substr(pos + 2, string::npos));
        iss >> dimSize;

        iss.clear();
        iss.str(line);
        iss >> skipws >> line;

        //     cout << "Dimension Name: " << line.c_str() << " Dimension Size: "
        //          << dimSize << endl;

    } // while loop

    // Loop through global attributes
    string attVal;
    char history[256];

    while (1) {
        getline(cdl_data_structure, line);
        size_t pos = line.find("group:");
        // If found then break from loop and process groups and variables
        if (pos == 0) {
            break;
        }

        if (line.length() == 0) continue;
        if (line.substr(0, 1).compare("\r") == 0) continue;
        if (line.substr(0, 1).compare("\n") == 0) continue;

        size_t firstNonBlank = line.find_first_not_of(" ");
        if (firstNonBlank != string::npos)
            if (line.compare(firstNonBlank, 2, "//") == 0) continue;

        pos = line.find(":");

        //  ":" found, new global attribute to create
        if (pos != string::npos) {
            size_t posEql = line.find("=");
            size_t pos1qte = line.find("\"");
            size_t pos2qte = line.substr(pos1qte + 1, string::npos).find("\"");
            cout << line.substr(pos + 1, posEql - pos - 2).c_str() << endl;

            string attrName = line.substr(pos + 1, posEql - pos - 2);
            // Get long_name
            attVal.assign(line.substr(pos1qte + 1, pos2qte));
            if (attrName.compare("product_name") == 0) {
                //             cout << "lname: " << lname.c_str() << endl;
                status = nc_put_att_text(ncid, NC_GLOBAL, attrName.c_str(), strlen(l1_filename) + 1, l1_filename);
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), attrName.c_str());
                    exit(1);
                }
            } else if (attrName.compare("date_created") == 0) {
                status = nc_put_att_text(ncid, NC_GLOBAL, attrName.c_str(), strlen(ydhmsf(now(), 'G')) + 1, ydhmsf(now(), 'G'));
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), attrName.c_str());
                    exit(1);
                }
            } else if (attrName.compare("history") == 0) {
                sprintf(history, "Generated by avirisbil2oci; cdlfile=%s", cdl_filename);
                status = nc_put_att_text(ncid, NC_GLOBAL, attrName.c_str(), strlen(history) + 1, history);
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), attrName.c_str());
                    exit(1);
                }
            } else if (attrName.compare("orbit_number") == 0) {
                int32_t vr[1];
                vr[0] = 0;
                status = nc_put_att_int(ncid, NC_GLOBAL, attrName.c_str(), NC_INT, 1, &vr[0]);
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), attrName.c_str());
                    exit(1);
                }
            } else if (attrName.compare("license") == 0) {
                printf("length of name=%ld\n", attrName.length());
                status = nc_put_att_text(ncid, NC_GLOBAL, "license",
                        strlen(line.substr(pos1qte + 1, pos2qte).c_str()) + 1,
                        line.substr(pos1qte + 1, pos2qte).c_str());
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), history);
                    exit(1);
                }
            } else if (attrName.compare("creator_name") == 0) {
                printf("length of name=%ld\n", attrName.length());
                status = nc_put_att_text(ncid, NC_GLOBAL, "creator_name",
                        strlen(attVal.c_str()) + 1,
                        attVal.c_str());
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), history);
                    exit(1);
                }

            } else if (attrName.compare("creator_url") == 0) {
                printf("length of name=%ld\n", attrName.length());
                status = nc_put_att_text(ncid, NC_GLOBAL, "creator_url",
                        strlen(attVal.c_str()) + 1,
                        attVal.c_str());
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), history);
                    exit(1);
                }
            } else if (attrName.compare("publisher_name") == 0) {
                printf("length of name=%ld\n", attrName.length());
                status = nc_put_att_text(ncid, NC_GLOBAL, "publisher_name",
                        strlen(attVal.c_str()) + 1,
                        attVal.c_str());
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), history);
                    exit(1);
                }

            } else if (attrName.compare("publisher_url") == 0) {
                printf("length of name=%ld\n", attrName.length());
                status = nc_put_att_text(ncid, NC_GLOBAL, "publisher_url",
                        strlen(attVal.c_str()) + 1,
                        attVal.c_str());
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), history);
                    exit(1);
                }



            } else {
                printf("value = %s\n", line.substr(pos1qte + 1, pos2qte).c_str());
                status = nc_put_att_text(ncid, NC_GLOBAL, attrName.c_str(),
                        strlen(attVal.c_str()) + 1,
                        attVal.c_str());
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for %s\n",
                            __FILE__, __LINE__, nc_strerror(status), attrName.c_str());
                    exit(1);
                }

            }
        }
    }
    ngrps = 0;
    // Loop through groups
    while (1) {

        size_t firstNonBlank = line.find_first_not_of(" ");
        if (firstNonBlank != string::npos)
            if (line.compare(firstNonBlank, 2, "//") == 0) {
                getline(cdl_data_structure, line);
                continue;
            }
        // Check if end of CDL file
        // If so then close CDL file and return
        if (line.substr(0, 1).compare("}") == 0) {
            cdl_data_structure.close();
            return 0;
        }

        // Check for beginning of new group
        size_t pos = line.find("group:");

        // If found then create new group and variables
        if (pos == 0) {

            // Parse group name
            istringstream iss(line.substr(6, string::npos));
            iss >> skipws >> line;
            cout << "Group: " << line.c_str() << endl;

            // Create NCDF4 group
            status = nc_def_grp(ncid, line.c_str(), &this->gid[ngrps]);
            check_err(status, __LINE__, __FILE__);

            ngrps++;

            int numDims = 0;
            int varDims[NC_MAX_DIMS];
            size_t dimSize[NC_MAX_DIMS];
            char dimName[NC_MAX_NAME + 1];
            string sname;
            string lname;
            string standard_name;
            string units;
            string flag_values;
            string flag_meanings;
            double valid_min = 0.0;
            double valid_max = 0.0;
            double fill_value = 0.0;
            size_t *chunksize;

            int ntype = 0;

            // Loop through datasets in group
            // Skip until "variables:" found
            while (1) {
                getline(cdl_data_structure, line);
                if (line.find("variables:") != string::npos) break;
            }

            while (1) {
                getline(cdl_data_structure, line);

                if (line.length() == 0) continue;
                if (line.substr(0, 1).compare("\r") == 0) continue;
                if (line.substr(0, 1).compare("\n") == 0) continue;

                size_t firstNonBlank = line.find_first_not_of(" ");
                if (firstNonBlank != string::npos)
                    if (line.compare(firstNonBlank, 2, "//") == 0) continue;

                size_t pos = line.find(":");

                // No ":" found, new dataset or empty line or end-of-group
                if (pos == string::npos) {

                    if (numDims > 0) {
                        // Create previous dataset
                        createNCDF(gid[ngrps - 1],
                                sname.c_str(), lname.c_str(),
                                standard_name.c_str(), units.c_str(),
                                (void *) &fill_value,
                                flag_values.c_str(), flag_meanings.c_str(),
                                valid_min, valid_max, ntype, numDims, varDims, chunksize);

                        flag_values.assign("");
                        flag_meanings.assign("");
                        free(chunksize);
                    }

                    valid_min = 0.0;
                    valid_max = 0.0;
                    fill_value = 0.0;

                    if (line.substr(0, 10).compare("} // group") == 0) break;

                    // Parse variable type
                    string varType;
                    istringstream iss(line);
                    iss >> skipws >> varType;

                    // Get corresponding NC variable type
                    if (varType.compare("char") == 0) ntype = NC_CHAR;
                    else if (varType.compare("byte") == 0) ntype = NC_BYTE;
                    else if (varType.compare("short") == 0) ntype = NC_SHORT;
                    else if (varType.compare("int") == 0) ntype = NC_INT;
                    else if (varType.compare("long") == 0) ntype = NC_INT;
                    else if (varType.compare("float") == 0) ntype = NC_FLOAT;
                    else if (varType.compare("real") == 0) ntype = NC_FLOAT;
                    else if (varType.compare("double") == 0) ntype = NC_DOUBLE;
                    else if (varType.compare("ubyte") == 0) ntype = NC_UBYTE;
                    else if (varType.compare("ushort") == 0) ntype = NC_USHORT;
                    else if (varType.compare("uint") == 0) ntype = NC_UINT;
                    else if (varType.compare("int64") == 0) ntype = NC_INT64;
                    else if (varType.compare("uint64") == 0) ntype = NC_UINT64;

                    // Parse short name (sname)
                    pos = line.find("(");
                    size_t posSname = line.substr(0, pos).rfind(" ");
                    sname.assign(line.substr(posSname + 1, pos - posSname - 1));
                    cout << "sname: " << sname.c_str() << endl;

                    // Parse variable dimension info
                    this->parseDims(line.substr(pos + 1, string::npos),
                            &numDims, varDims);
                    for (int i = 0; i < numDims; i++) {
                        nc_inq_dim(ncid, varDims[i], dimName, &dimSize[i]);
                        cout << line.c_str() << " " << i << " " << dimName
                                << " " << dimSize[i] << endl;
                    }
                    chunksize = (size_t *) calloc(numDims, sizeof (size_t));

                } else {
                    // Parse variable attributes
                    size_t posEql = line.find("=");
                    size_t pos1qte = line.find("\"");
                    size_t pos2qte = line.substr(pos1qte + 1, string::npos).find("\"");
                    //    cout << line.substr(pos+1, posEql-pos-2).c_str() << endl;

                    string attrName = line.substr(pos + 1, posEql - pos - 2);

                    // Get long_name
                    if (attrName.compare("long_name") == 0) {
                        lname.assign(line.substr(pos1qte + 1, pos2qte));
                        //             cout << "lname: " << lname.c_str() << endl;
                    }
                        // Get units
                    else if (attrName.compare("units") == 0) {
                        units.assign(line.substr(pos1qte + 1, pos2qte));
                        //             cout << "units: " << units.c_str() << endl;
                    }
                        // Get _FillValue
                    else if (attrName.compare("_FillValue") == 0) {
                        iss.clear();
                        iss.str(line.substr(posEql + 1, string::npos));
                        iss >> fill_value;
                        //             cout << "_FillValue: " << fill_value << endl;
                    }
                        // Get flag_values
                    else if (attrName.compare("flag_values") == 0) {
                        flag_values.assign(line.substr(pos1qte + 1, pos2qte));
                    }
                        // Get flag_meanings
                    else if (attrName.compare("flag_meanings") == 0) {
                        flag_meanings.assign(line.substr(pos1qte + 1, pos2qte));
                    }
                        // Get valid_min
                    else if (attrName.compare("valid_min") == 0) {
                        iss.clear();
                        iss.str(line.substr(posEql + 1, string::npos));
                        iss >> valid_min;
                        //             cout << "valid_min: " << valid_min << endl;
                    }
                        // Get valid_max
                    else if (attrName.compare("valid_max") == 0) {
                        iss.clear();
                        iss.str(line.substr(posEql + 1, string::npos));
                        iss >> valid_max;
                        //             cout << "valid_max: " << valid_max << endl;
                    }                        // Get Chunk sizes
                    else if (attrName.compare("_ChunkSizes") == 0) {
                        size_t posComma = line.find(",");
                        status = posEql + 1;

                        for (i = 0; i < numDims; i++) {
                            iss.clear();
                            if (i == numDims - 1) {
                                //            		printf("iss=%s\n",line.substr(status, string::npos).c_str());
                                iss.str(line.substr(status, string::npos));
                            } else {
                                iss.str(line.substr(status, posComma));
                                //            		printf("iss=%s\n",line.substr(status, posComma).c_str());
                                status = posComma + 1;
                                posComma += line.substr(posComma + 1, string::npos).find(",") + 1;
                            }
                            iss >> chunksize[i];
                        }
                        //             cout << "_FillValue: " << fill_value << endl;
                    }

                } // if ( pos == string::npos)
            } // datasets in group loop
        } // New Group loop
        getline(cdl_data_structure, line);
    } // Main Group loop


    return 0;
}

int ncdfFile::parseDims(string dimString, int *numDims, int *varDims) {

    size_t dimSize, curPos = 0;
    char dimName[NC_MAX_NAME + 1];

    *numDims = 0;

    while (1) {
        size_t pos = dimString.find(",", curPos);
        if (pos == string::npos)
            pos = dimString.find(")");

        string varDimName;
        istringstream iss(dimString.substr(curPos, pos - curPos));
        iss >> skipws >> varDimName;

        for (int i = 0; i < ndims; i++) {
            nc_inq_dim(ncid, dimid[i], dimName, &dimSize);
            if (varDimName.compare(dimName) == 0) {
                varDims[(*numDims)++] = dimid[i];
                break;
            }
        }
        if (dimString.substr(pos, 1).compare(")") == 0) break;

        curPos = pos + 1;
    }

    return 0;
}

int ncdfFile::getGid(const char *grpName) {

    int status;
    int grpID;
    status = nc_inq_grp_ncid(ncid, grpName, &grpID);
    check_err(status, __LINE__, __FILE__);

    return grpID;
}

int ncdfFile::getNcid() {

    return ncid;
}

int ncdfFile::close() {

    nc_close(ncid);

    return 0;
}

