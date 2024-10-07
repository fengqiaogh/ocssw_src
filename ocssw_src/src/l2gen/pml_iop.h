#ifndef tjsm_pml_iop
#define tjsm_pml_iop

/* Switch on/off iteration flag for scattering */
#define NFLAG 1

/* Number of bands within the LUTs */
#define NB 6

/* Define maximum bands for other processors */
#define MAX_BANDS 16

/* Functions */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define radians(degrees) ((degrees) * M_PI / 180.0)
#define degrees(radians) ((radians) * 180.0 / M_PI)

#endif
