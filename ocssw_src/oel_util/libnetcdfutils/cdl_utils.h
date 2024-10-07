#include <stdint.h>
#include <iostream>
#include <fstream>

#define PBUFFER_SIZE 32768

using namespace std;

class ncdfFile {
    int ncid;

    int ngrps;
    int gid[10];

    int ndims;
    int dimid[1000];

public:
    ncdfFile();
    ~ncdfFile();

    int cdlCreate(char* l1_filename, char* cdl_filename, int32_t numScans);

    int cdlCreateDim(char* l1_filename, char* cdl_filename, const char** dim_names, size_t* dim_size, size_t n_dims, size_t numScans);

    int parseDims(string dimString, int *numDims, int *varDims);

    int getGid(const char *grpName);
    int getNcid();
    int close();
};

inline
int expandEnvVar(string *sValue) {
    if ((*sValue).find_first_of("$") == string::npos) return 0;
    string::size_type posEndIdx = (*sValue).find_first_of("/");
    if (posEndIdx == string::npos) return 0;
    char *envVar_str = getenv((*sValue).substr(1, posEndIdx - 1).c_str());
    if (envVar_str == 0x0) {
        printf("Environment variable: %s not defined.\n", envVar_str);
        exit(1);
    }
    *sValue = envVar_str + (*sValue).substr(posEndIdx);

    return 0;
}
