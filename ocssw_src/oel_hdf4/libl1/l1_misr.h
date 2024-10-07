#ifndef L1_MISR_H
#define L1_MISR_H

#define N_CAMERAS 9

#define D2R 0.017453292520

typedef struct misr_struct {
  uint8_t multipleInput;
  
  int32_t startBlock;
  int32_t endBlock;

  int32_t fileID[N_CAMERAS];
  int32_t blockTimeID[N_CAMERAS];

  int32_t ocean_block_numbers[180];
  int8_t isOceanBlock[180];
  int8_t offset[180];
  
  double radScaleFactors[4];
  
  double SolAzimuth[180*8][32];
  double SolZenith[180*8][32];

  double SenAzimuth[N_CAMERAS][180*8][32];
  double SenZenith[N_CAMERAS][180*8][32];
} misr_t;

int openl1_misr(filehandle *l1file);
int readl1_misr(filehandle *l1file, l1str *l1rec);
int closel1_misr(filehandle *l1file);

#endif
