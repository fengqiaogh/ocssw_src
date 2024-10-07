
   
using namespace std;

class l1aFile {

  int ngrps;

  int ndims;
  int dimid[1000];
  
  uint32_t SC_records_val;

 public:
  l1aFile();
  ~l1aFile();

  int ncid;
  int gid[10];

  int createl1( char* l1_filename, uint32_t nSC,
                uint32_t imgWidth, uint32_t imgHeight,
                uint32_t fndWidth, uint32_t fndHeight);
 
  int parseDims( string dimString, int *numDims, int *varDims);

  int close();
};

l1aFile::l1aFile() {
  ncid = -1;
  SC_records_val = 0;
}


l1aFile::~l1aFile() {
  
}

typedef struct {
  int16_t iyr;
  int16_t iday;
  double sec;
  uint8_t bus_telemetry[27];
} sensor_t;

typedef struct {
  int16_t iyr;
  int16_t iday;
  double sec;
  float gyro[3];
  float mag1[3];
  float mag2[3];
  float rwheels[3];
  float sunb[3];
} psensor_t;

typedef struct {
  int16_t iyr;
  int16_t iday;
  double sec;
  float quat[4];
} attitude_t;

typedef struct {
  int16_t iyr;
  int16_t iday;
  double sec;
  float pos[3];
  float vel[3];
} propagator_t;


int createNCDF( int ncid, const char *sname, const char *lname, 
                const char *standard_name, const char *units,
                void *fill_value,
                const char *flag_values, const char *flag_meanings,
                double low, double high, int nt,
                int rank, int *dimids);

int unpack_seahawk_adcs( uint8_t *apkt, double startHWKTime, double stopHWKTime,
                         sensor_t *sensor, psensor_t *psensor,
                         attitude_t *attitude, propagator_t *propagator);

#define MAXNPACKETS 3000

#define MAXIMGHEIGHT 6000


#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )

enum {SENSOR, PSENSOR, ATTITUDE, PROPAGATOR};

/*
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
*/
