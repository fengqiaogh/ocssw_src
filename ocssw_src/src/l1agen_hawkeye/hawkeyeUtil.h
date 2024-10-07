#include "netcdf.h"
#include "nc4utils.h"
#include "timeutils.h"

const int PHDRLEN = 6; // length of primary header
const int SHDRLEN = 6; // length of secondary header
const int PUSLEN = 3;  // length of PUS header

using namespace std;


int createFile( const char* filename, const char* cdlfile,
                size_t sdim, int *ncid, int *gid);
 
int parseDims( int ncid, int ndims, string dimString, int *numDims,
               int *dimid, int *varDims);

int createNCDF( int ncid, const char *sname, const char *lname, 
                const char *standard_name, const char *units,
                void *fill_value,
                const char *flag_values, const char *flag_meanings,
                double low, double high, 
                float scale_factor, float add_offset,
                int nt, int rank, int *dimids);

int tepoch2yds( double tepoch, int32_t *iyr, int32_t *idy, double *sec);

/*
int get_packet_from_frame( ifstream *framefile, uint8_t **packet,
                           uint32_t &packet_length, int first);
*/

int readFrame ( ifstream *framefile, uint8_t frame[892], int& prevFrameCnt,
                int ierror[2]);

int getSHpacket( ifstream *framefile, uint8_t frame[892], int& framePtr,
                 uint8_t **packet, int& packetLength, int& prevFrameCnt,
                 int& frameDrop, bool isVerbose);

inline
int expandEnvVar( string *sValue) {
  if ( (*sValue).find_first_of( "$" ) == string::npos) return 0;
  string::size_type posEndIdx = (*sValue).find_first_of( "/" );
  if ( posEndIdx == string::npos) return 0;
  char *envVar_str = getenv((*sValue).substr( 1, posEndIdx-1 ).c_str());
  if (envVar_str == 0x0) {
    printf("Environment variable: %s not defined.\n", envVar_str);
    exit(1);
  }
  *sValue = envVar_str + (*sValue).substr( posEndIdx);

  return 0;
}

