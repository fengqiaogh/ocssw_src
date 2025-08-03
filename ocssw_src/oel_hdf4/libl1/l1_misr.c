#include "l1.h"
#include "l1_misr.h"
#include <hdf4utils.h>

#include <hdf.h>
#include <mfhdf.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_multifit.h>

#include <libgen.h>

static gsl_spline2d *spline_misr;
static gsl_interp_accel *xacc, *yacc;
static double xa[32], ya[2];

static int32_t geoFileId;
static int32_t lat_id;
static int32_t lon_id;

static int32_t sd_id[N_CAMERAS];
static int32_t red_id[N_CAMERAS];
static int32_t grn_id[N_CAMERAS];
static int32_t blu_id[N_CAMERAS];

static int32_t nir_id[N_CAMERAS];

static int icamera;
static int single_ifile=0;

int getRadScaleFactors( char *file, double rad_scl_factors[4]);
int reduce_res(uint16_t rad_data[4][2048]);
int interp_values( double *grid_values, double interpolated_values[16][512]);
int interp_values_dbl( int32_t diff_offset, double *grid_values,
                       double interpolated_values[16][512]);

int openl1_misr(filehandle *file) {
    int32_t i;
    int32_t sds_id;
    int32_t gmp_id;
    int32_t refn;
    int32_t start[3]={0,0,0};
    int32_t count[3]={180,8,32};
    intn status;
    char camera_type[N_CAMERAS][3]=
      {"DF","CF","BF","AF","AN","AA","BA","CA","DA"};
        char camera_type_mc[N_CAMERAS][3]=
      {"Df","Cf","Bf","Af","An","Aa","Ba","Ca","Da"};
    char namebuf[512];

    
    misr_t *private_data = file->private_data;
    //(misr_t *) calloc(1, sizeof(misr_t));

    //file->private_data = private_data;

    for (i=0; i<N_CAMERAS; i++) sd_id[i] = -1;
                                         
    if ( private_data->multipleInput == 1) {
      // Multiple input files

      for (i=0; i<N_CAMERAS; i++) {
        strcpy(namebuf, dirname(file->name));
        strcat(namebuf, "/");
        strncat(namebuf, basename(file->name), 39);
        strcat(namebuf, camera_type[i]);
        strcat(namebuf, basename(file->name)+41);
        sd_id[i] = SDstart(namebuf, DFACC_RDONLY);
        if (sd_id[i] == FAIL) {
          fprintf(stderr, "-E- %s line %d: SDstart(%s, %d) failed.\n",
                  __FILE__, __LINE__, namebuf, DFACC_RDONLY);
          return (HDF_FUNCTION_ERROR);
        }

        // Read radiance scale factors
        if (i == 0)
          getRadScaleFactors( namebuf, private_data->radScaleFactors);
        
        // Get BlockTime IDs
        private_data->fileID[i] = Hopen(namebuf, DFACC_RDONLY, 0);
        Vstart ( private_data->fileID[i]);
        refn = VSfind( private_data->fileID[i], "PerBlockMetadataTime");
        private_data->blockTimeID[i] =
          VSattach( private_data->fileID[i], refn, "r");
      }  // camera loop
      icamera = 0;
      //      private_data->isSingleFile = 0;
    } else {
      // Single input file
      
      for (i=0; i<N_CAMERAS; i++) {
        if ( strncmp(basename(file->name)+39, camera_type[i], 2) == 0) {
          icamera = i;
          single_ifile = 1;
          //          private_data->isSingleFile = 1;
          sd_id[icamera] = SDstart(file->name, DFACC_RDONLY);
          if (sd_id[icamera] == FAIL) {
            fprintf(stderr, "-E- %s line %d: SDstart(%s, %d) failed.\n",
                    __FILE__, __LINE__, file->name, DFACC_RDONLY);
            return (HDF_FUNCTION_ERROR);
          }

          // Read radiance scale factors
          getRadScaleFactors( file->name, private_data->radScaleFactors);
          // Get BlockTime IDs
          private_data->fileID[icamera] = Hopen(file->name, DFACC_RDONLY, 0);
          Vstart ( private_data->fileID[icamera]);
          refn = VSfind( private_data->fileID[icamera], "PerBlockMetadataTime");
          private_data->blockTimeID[icamera] =
            VSattach( private_data->fileID[icamera], refn, "r");

          break;
        }
      }  // camera loop
    }
      
    for (i=0; i<N_CAMERAS; i++) {
      if (sd_id[i] == -1) continue;
      red_id[i] =
        SDselect(sd_id[i],SDnametoindex(sd_id[i], "Red Radiance/RDQI"));
      grn_id[i] =
        SDselect(sd_id[i], SDnametoindex(sd_id[i], "Green Radiance/RDQI"));
      blu_id[i] =
        SDselect(sd_id[i], SDnametoindex(sd_id[i], "Blue Radiance/RDQI"));
      nir_id[i] =
        SDselect(sd_id[i], SDnametoindex(sd_id[i], "NIR Radiance/RDQI"));
    }  // camera loop
    
    // Open the AGP geofile
    //
    // Dimensions:  Block    180
    //              Row      128
    //              Col      512

    // -111 = Fill above data
    // -222 = Fill below data
    // -333 = Fill IPI invalid
    // -444 = Fill to side of data
    // -555 = Fill not processed
    // -999 = Fill IPI error
    
    geoFileId = SDstart(file->geofile, DFACC_RDONLY);
    if (geoFileId == FAIL) {
        fprintf(stderr, "-E- %s line %d: SDstart(%s, %d) failed.\n",
                __FILE__, __LINE__, file->geofile, DFACC_RDONLY);
        return (HDF_FUNCTION_ERROR);
    }

    lat_id = SDselect(geoFileId, SDnametoindex(geoFileId, "GeoLatitude"));
    lon_id = SDselect(geoFileId, SDnametoindex(geoFileId, "GeoLongitude"));


    // Open the GMP file
    //
    // Dimensions:  Block    180
    //              Row        8
    //              Col       32
    //
    // Camera order: DF CF BF AF AN AA BA CA DA
    //

    // Azimuth: 0.0 to 360.0
    // Zenith:  0.0 to 90.0
    
    // Fill Value: -555
    gmp_id = SDstart(file->gmpfile, DFACC_RDONLY);
    if (gmp_id == FAIL) {
        fprintf(stderr, "-E- %s line %d: SDstart(%s, %d) failed.\n",
                __FILE__, __LINE__, file->gmpfile, DFACC_RDONLY);
        return (HDF_FUNCTION_ERROR);
    }

    // Read Solar Azimuth
    sds_id = SDselect(gmp_id, SDnametoindex(gmp_id, "SolarAzimuth"));
    status = SDreaddata(sds_id, start, NULL, count,
                        (VOIDP) &private_data->SolAzimuth);
    if (status != 0) {
      printf("-E- %s line %d:  Error reading SDS %s from %s.\n", __FILE__,
             __LINE__, "SolAzimuth", file->gmpfile);
      exit(1);
    }
    SDendaccess(sds_id);

    // Read Solar Zenith
    sds_id = SDselect(gmp_id, SDnametoindex(gmp_id, "SolarZenith"));
    status = SDreaddata(sds_id, start, NULL, count,
                        (VOIDP) &private_data->SolZenith);
    if (status != 0) {
      printf("-E- %s line %d:  Error reading SDS %s from %s.\n", __FILE__,
             __LINE__, "SolZenith", file->gmpfile);
      exit(1);
    }
    SDendaccess(sds_id);

    // Read Sensor Azimuth/Zenith
    for (i=0; i<N_CAMERAS; i++) {
      if (sd_id[i] == -1) continue;
      
      strcpy(namebuf, camera_type_mc[i]);
      strcat(namebuf, "Azimuth");
      sds_id = SDselect(gmp_id, SDnametoindex(gmp_id, namebuf));
      status = SDreaddata(sds_id, start, NULL, count,
                          (VOIDP) &private_data->SenAzimuth[i]);
      if (status != 0) {
        printf("-E- %s line %d:  Error reading SDS %s from %s.\n", __FILE__,
               __LINE__, namebuf, file->gmpfile);
        exit(1);
      }
      SDendaccess(sds_id);

      strcpy(namebuf, camera_type_mc[i]);
      strcat(namebuf, "Zenith");
      sds_id = SDselect(gmp_id, SDnametoindex(gmp_id, namebuf));
      status = SDreaddata(sds_id, start, NULL, count,
                          (VOIDP) &private_data->SenZenith[i]);
      if (status != 0) {
        printf("-E- %s line %d:  Error reading SDS %s from %s.\n", __FILE__,
               __LINE__, namebuf, file->gmpfile);
        exit(1);
      }
      SDendaccess(sds_id);
    }  // camera loop
    
    SDend(gmp_id);


    // Read start and end block numbers
    SDreadattr(sd_id[icamera],
               SDfindattr(sd_id[icamera], "Start_block"),
               (VOIDP) &private_data->startBlock);
    SDreadattr(sd_id[icamera],
               SDfindattr(sd_id[icamera], "End block"),
               (VOIDP) &private_data->endBlock);
    
    file->nbands = 4;
    file->npix = 512;
    file->nscan = 128*(private_data->endBlock);
    if (single_ifile == 0) file->nscan *= N_CAMERAS;


    // Read ocean block numbers (1-based)
    SDreadattr(sd_id[icamera],
               SDfindattr(sd_id[icamera], "Ocean_blocks.numbers"),
               (VOIDP) private_data->ocean_block_numbers);
    memset(&private_data->isOceanBlock, 0, 180);
    for (i=0; i<180; i++) {
      int j = private_data->ocean_block_numbers[i];
      if (j != 0) private_data->isOceanBlock[j-1] = 1;
    }

    // Compute block offsets
    memset(&private_data->offset, 0, 180);
    for (i=0; i<179; i++) {
      double p1 = private_data->SolZenith[8*i+6][15];
      double p2 = private_data->SolZenith[8*i+7][15];
      double f1 = private_data->SolZenith[8*(i+1)+0][15];
      //double f2 = private_data->SolZenith[8*(i+1)+1][15];
      double f1_l = private_data->SolZenith[8*(i+1)+0][14];
      double f1_r = private_data->SolZenith[8*(i+1)+0][16];
      double diff1a = fabs((p2-p1)-(f1-p2));
      double diff1b = fabs((p2-p1)-(f1_l-p2));
      double diff1c = fabs((p2-p1)-(f1_r-p2));
      if (diff1a < diff1b && diff1a < diff1c)
        private_data->offset[i] =  0;
      if (diff1b < diff1a && diff1b < diff1c)
        private_data->offset[i] = -1;
      if (diff1c < diff1a && diff1c < diff1b)
        private_data->offset[i] = +1;
    }

    
    // Set up interpolation grid

    // 512 / 32 = 16
    // (16-1)/2 = 7.5
    const gsl_interp2d_type *T = gsl_interp2d_bilinear;
    spline_misr = gsl_spline2d_alloc(T, 32, 2);
    xacc = gsl_interp_accel_alloc();
    yacc = gsl_interp_accel_alloc();

    for (size_t i=0; i<32; i++) xa[i] = 7.5 + 16*i;
    ya[0] = -0.5;
    ya[1] = 15.5;

    return 0;
}


int getRadScaleFactors( char *file, double rad_scl_factors[4]) {

  int i,j;
  int32_t file_id, vg_band_ref, refn, vg_id, vd_id;
  int32_t nObj, tag, ref;
  char name[100];
  char *grpTags[4]={"BlueBand", "GreenBand", "RedBand", "NIRBand"};
  
  file_id = Hopen(file, DFACC_RDONLY, 0);
  Vstart (file_id);

  // Loop over bands
  for (i=0; i<4; i++) {

    // Find "Grid Attributes" group for each band
    vg_band_ref = Vfind(file_id, grpTags[i]);
    refn = vg_band_ref;
    while (1) {
      refn = Vgetid(file_id, refn);
      vg_id = Vattach(file_id, refn, "r");
      Vgetname(vg_id, name);
      if (strcmp(name, "Grid Attributes") == 0) break;
      Vdetach(vg_id);
    }
    nObj = Vntagrefs(vg_id);

    // Search within"Grid Attributes" group for "Scale factor" vdata 
    for (j=0; j<nObj; j++) {
      Vgettagref(vg_id, j, &tag, &ref);
      vd_id = VSattach(file_id, ref, "r");
      VSgetname(vd_id, name);
      if (strcmp(name, "Scale factor") == 0) {
        VSread(vd_id, (uint8 *) &rad_scl_factors[i], 1, NO_INTERLACE);
        VSdetach(vd_id);
        break;
      }
      VSdetach(vd_id);
    }
    rad_scl_factors[i] /= 10;
  } // band loop
  Vdetach(vg_id);
  Vend(file_id);
  Hclose(file_id);

  return 0;
}


int readl1_misr(filehandle *l1file, l1str *l1rec) {

  int i,j,k, jcamera;

  int32_t start[3]={0,0,0};
  int32_t count[3]={1,1,512};
  
  //static int32_t last_recnum=-1;
  static int32_t last_recnum_gp=-1;

  //int32_t diff_offset, interp_index;
  int32_t recnum, recnum_red, recnum_gp, block;
  
  ushort rad_data[4][4][2048]; 

  double dbl_data[2*32];
  double dbl_xy[2*32];
  ushort *usptr;
  char timestring[28];
  int32_t year, day, msec;
  int16_t yr16, day16;
  double sec;
 
  static double sla_interp_data[16][512];
  static double slz_interp_data[16][512];
  static double sna_interp_data[N_CAMERAS][16][512];
  static double snz_interp_data[N_CAMERAS][16][512];

  static double xy_interp_data[2][16][512];
  static int32_t gp_index;

  static int first=1;
  intn status;

  misr_t *private_data = l1file->private_data;

  // Initialize to fill value
  if (first == 1) {
    for (i=0; i<16; i++) {
      for (j=0; j<512; j++) {
        sla_interp_data[i][j] = -32767;
        slz_interp_data[i][j] = -32767;
      }
    }
  
    for (k=0; k<N_CAMERAS; k++) {
      for (i=0; i<16; i++) {
        for (j=0; j<512; j++) {
          //          sna_interp_data[k][i][j] = -32767;
          //snz_interp_data[k][i][j] = -32767;

          sna_interp_data[k][i][j] = -555;
          snz_interp_data[k][i][j] = -555;
        }
      }
    }
  }
  first = 0;
  
  recnum = l1rec->iscan;
  if (single_ifile == 0) recnum /= N_CAMERAS;
  
  block = recnum / 128; 
  
  // If no ocean data set time fill and bad nav and bail
  //  if(private_data->isOceanBlock[block] == 0) {
  // l1rec->scantime = BAD_FLT;
  // for (i=0; i<l1rec->npix; i++) l1rec->navfail[i] = 1;
  //return 0;
  //}
  
  // Bail if record already processed
  //if (recnum == last_recnum) last_recnum_gp = -1;

  // Get corresponding camera number for scan
  if (single_ifile) jcamera = icamera; else jcamera = l1rec->iscan % N_CAMERAS;
  // Get scantime (from PerBlockMetadataTime vdata)
  //i = 0;
  //  while (1) {
  //if (block+1 == private_data->ocean_block_numbers[i]) break;
  //i++;
  //}

  VSseek(private_data->blockTimeID[jcamera], block);
  VSread(private_data->blockTimeID[jcamera], (uint8 *) timestring, 1,
         NO_INTERLACE);

  // Convert time string (isodate) to seconds
  isodate2ydmsec(timestring, &year, &day, &msec);
  yr16 = year;
  day16 = day;
  sec = (double) msec / 1000;

  if (yr16 > 1999)
      l1rec->scantime = yds2unix(yr16, day16, sec);
  else
      return 0;
  ///////////////////////
  // Read lon/lat data //
  ///////////////////////
  start[0] = block;
  start[1] = recnum - start[0]*128;

  status = SDreaddata(lat_id, start, NULL, count, (VOIDP) l1rec->lat);
  if (status != 0) {
    printf("-E- %s line %d:  Error reading SDS %s from %s.\n", __FILE__,
           __LINE__, "GeoLatitude", l1file->geofile);
    exit(1);
  }

  status = SDreaddata(lon_id, start, NULL, count, (VOIDP) l1rec->lon);
  if (status != 0) {
    printf("-E- %s line %d:  Error reading SDS %s from %s.\n", __FILE__,
           __LINE__, "GeoLongitude", l1file->geofile);
    exit(1);
  }

  
  //////////////////////////////////////////////////////////////
  // Generate interpolated Azimuth/Zenith values if necessary //
  //////////////////////////////////////////////////////////////
  recnum_gp = (recnum - 8) / 16;
  if (recnum_gp != last_recnum_gp && recnum >= 8) {

    gp_index = 0;
    
//    printf("recnum_gp: %d  block: %d  %d  jcamera: %d\n",
//           recnum_gp, block, recnum_gp % 8, jcamera);

//    interp_index = recnum_gp % 8;
//    if (interp_index == 7) {
//      diff_offset = private_data->offset[block];
//    } else {
//      diff_offset = 0;
//    }

    if (single_ifile || (l1rec->iscan % N_CAMERAS) == 0) {
      //////////////////////////////////////
      // Interpolate Solar Azimuth values //
      //////////////////////////////////////
      for (i=0; i<32; i++) {
        dbl_data[i]    = private_data->SolAzimuth[recnum_gp+0][i];
        if (dbl_data[i] >= 0 && dbl_data[i] < 180)
          dbl_data[i] = dbl_data[i] + 360;
        dbl_data[i+32] = private_data->SolAzimuth[recnum_gp+1][i];
        if (dbl_data[i+32] >= 0 && dbl_data[i+32] < 180)
          dbl_data[i+32] = dbl_data[i+32] + 360;
      }
      interp_values( dbl_data, sla_interp_data);

      /////////////////////////////////////
      // Interpolate Solar Zenith values //
      /////////////////////////////////////
      for (i=0; i<32; i++) {
        dbl_data[i]    = private_data->SolZenith[recnum_gp+0][i];
        dbl_data[i+32] = private_data->SolZenith[recnum_gp+1][i];
      }
      interp_values( dbl_data, slz_interp_data);
    }
    ///////////////////////////////////////
    // Interpolate Sensor Azimuth values //
    ///////////////////////////////////////

    // In order to avoid problems interpolating over 0/360 jump the
    // unit vector xy-values will be used rather than the azimuth angle
    
    for (j=0; j<N_CAMERAS; j++) {
      if (sd_id[j] == -1) continue;
      
      for (i=0; i<32; i++) {
        dbl_xy[i] = -555;
        if ( private_data->SenAzimuth[j][recnum_gp+0][i] > 0)
          dbl_xy[i]    = cos(private_data->SenAzimuth[j][recnum_gp+0][i]*D2R);
        dbl_xy[i+32] = -555;
        if ( private_data->SenAzimuth[j][recnum_gp+1][i] > 0)
          dbl_xy[i+32] = cos(private_data->SenAzimuth[j][recnum_gp+1][i]*D2R);
      }
      interp_values( dbl_xy, xy_interp_data[0]);

      for (i=0; i<32; i++) {
        dbl_xy[i] = -555;
        if ( private_data->SenAzimuth[j][recnum_gp+0][i] > 0)
          dbl_xy[i]    = sin(private_data->SenAzimuth[j][recnum_gp+0][i]*D2R);
        dbl_xy[i+32] = -555;
        if ( private_data->SenAzimuth[j][recnum_gp+1][i] > 0)
          dbl_xy[i+32] = sin(private_data->SenAzimuth[j][recnum_gp+1][i]*D2R);
      }
      interp_values( dbl_xy, xy_interp_data[1]);

      for (i=0; i<16; i++) {
        for (k=0; k<512; k++) {
          if ( xy_interp_data[1][i][k] != -555 &&
               xy_interp_data[0][i][k] != -555)
            sna_interp_data[j][i][k] = atan2(xy_interp_data[1][i][k],
                                             xy_interp_data[0][i][k]) / D2R;
          else
            sna_interp_data[j][i][k] = -555;
        }
      }
    } // camera loop
        
    
    //////////////////////////////////////
    // Interpolate Sensor Zenith values //
    //////////////////////////////////////
    for (j=0; j<N_CAMERAS; j++) {
      if (sd_id[j] == -1) continue;
      for (i=0; i<32; i++) {
        dbl_data[i]    = private_data->SenZenith[j][recnum_gp+0][i];
        dbl_data[i+32] = private_data->SenZenith[j][recnum_gp+1][i];
      }
      interp_values( dbl_data, snz_interp_data[j]);
    }
    
    last_recnum_gp = recnum_gp;

  }  // Generate interpolated GMP values ... if (recnum_gp != last_recnum_gp)

  
  ///////////////////////////
  // Process radiance data //
  ///////////////////////////
  recnum_red = 4*recnum;
  start[0] = recnum_red / 512;
  start[1] = recnum_red - start[0]*512;
  count[1] = 4;
  count[2] = 2048;

  //  printf("Reading radiance record: %d\n", recnum);

  // Red
  status = SDreaddata(red_id[jcamera], start, NULL, count, (VOIDP) rad_data[2]);
  if (status != 0) {
    printf("-E- %s line %d:  Error reading SDS %s from %s.\n", __FILE__,
           __LINE__, "Red Radiance/RDQI", l1file->name);
    exit(1);
  }
  usptr = &rad_data[2][0][0];
  for (size_t i=0; i<count[1]*count[2]; i++) *usptr++ /= 4;
  
  reduce_res(rad_data[2]);

  if ((jcamera+1) != 5) {
    start[0] = recnum / 128;
    start[1] = recnum - start[0]*128;
    count[1] = 1;
    count[2] = 512;
  }

  // Blue
  status = SDreaddata(blu_id[jcamera], start, NULL, count, (VOIDP) rad_data[0]);
  if (status != 0) {
    printf("-E- %s line %d:  Error reading SDS %s from %s.\n", __FILE__,
           __LINE__, "Blue Radiance/RDQI", l1file->name);
    exit(1);
  }
  usptr = &rad_data[0][0][0];
  for (size_t i=0; i<count[1]*count[2]; i++) *usptr++ /= 4;
  if ((jcamera+1) == 5) reduce_res(rad_data[0]);

  
  // Green
  status = SDreaddata(grn_id[jcamera], start, NULL, count, (VOIDP) rad_data[1]);
  if (status != 0) {
    printf("-E- %s line %d:  Error reading SDS %s from %s.\n", __FILE__,
           __LINE__, "Green Radiance/RDQI", l1file->name);
    exit(1);
  }
  usptr = &rad_data[1][0][0];
  for (size_t i=0; i<count[1]*count[2]; i++) *usptr++ /= 4;
  if ((jcamera+1) == 5) reduce_res(rad_data[1]);

  // NIR
  status = SDreaddata(nir_id[jcamera], start, NULL, count, (VOIDP) rad_data[3]);
  if (status != 0) {
    printf("-E- %s line %d:  Error reading SDS %s from %s.\n", __FILE__,
           __LINE__, "NIR Radiance/RDQI", l1file->name);
    exit(1);
  }
  usptr = &rad_data[3][0][0];
  for (size_t i=0; i<count[1]*count[2]; i++) *usptr++ /= 4;
  if ((jcamera+1) == 5) reduce_res(rad_data[3]);


  // Convert scaled radiance values to floating point and copy to Lt
  for (size_t ip=8; ip<504; ++ip) {
    if (rad_data[0][0][ip] < 16378)
      l1rec->Lt[4*ip+0] = rad_data[0][0][ip]*private_data->radScaleFactors[0];
    if (rad_data[1][0][ip] < 16378)
      l1rec->Lt[4*ip+1] = rad_data[1][0][ip]*private_data->radScaleFactors[1];
    if (rad_data[2][0][ip] < 16378)
      l1rec->Lt[4*ip+2] = rad_data[2][0][ip]*private_data->radScaleFactors[2];
    if (rad_data[3][0][ip] < 16378)
      l1rec->Lt[4*ip+3] = rad_data[3][0][ip]*private_data->radScaleFactors[3];
  }

  // Copy Azimuth/Zenith angles to l1rec arrays
  for (size_t ip=8; ip<504; ++ip) {
    if (sla_interp_data[gp_index][ip] != -555) {
      l1rec->sola[ip] = sla_interp_data[gp_index][ip];
      if (l1rec->sola[ip] > 360) l1rec->sola[ip] = l1rec->sola[ip] - 360;
      if (l1rec->sola[ip] > 180) l1rec->sola[ip] = l1rec->sola[ip] - 360;
    }
    
    if (slz_interp_data[gp_index][ip] != -555) {
      l1rec->solz[ip] = slz_interp_data[gp_index][ip];
      l1rec->csolz[ip] = cos(l1rec->solz[ip] / OEL_RADEG);
    }
    
    if (sna_interp_data[jcamera][gp_index][ip] != -555) {
      l1rec->sena[ip] = sna_interp_data[jcamera][gp_index][ip];
      if (l1rec->sena[ip] > 360) l1rec->sena[ip] = l1rec->sena[ip] - 360;
      if (l1rec->sena[ip] > 180) l1rec->sena[ip] = l1rec->sena[ip] - 360;
    }
    
    if (snz_interp_data[jcamera][gp_index][ip] != -555) {
      l1rec->senz[ip] = snz_interp_data[jcamera][gp_index][ip];
      l1rec->csenz[ip] = cos(l1rec->senz[ip] / OEL_RADEG);
    }

    if (sla_interp_data[gp_index][ip] != -555 &&
        sna_interp_data[jcamera][gp_index][ip] != -555) {
      l1rec->delphi[ip] = l1rec->sena[ip] - 180.0 - l1rec->sola[ip];

      if (l1rec->delphi[ip] < -180.)
        l1rec->delphi[ip] += 360.0;
      else if (l1rec->delphi[ip] > 180.0)
        l1rec->delphi[ip] -= 360.0;
    }
  }
  
  if (single_ifile || (l1rec->iscan % N_CAMERAS) == 8) gp_index++;

  //last_recnum = recnum;
  
  return 0;
}


int interp_values( double *grid_values, double interpolated_values[16][512]) {

  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  static gsl_spline2d *spline;
  static gsl_interp_accel *xacc, *yacc;
  static double xa[32], ya[2];
  
  double dbl_data[2*32];

  static int first=1;
 
  // 512 / 32 = 16
  // (16-1)/2 = 7.5

  if (first) {
    // Set up interpolation grid
    spline = gsl_spline2d_alloc(T, 32, 2);
    xacc = gsl_interp_accel_alloc();
    yacc = gsl_interp_accel_alloc();

    for (size_t i=0; i<32; i++) xa[i] = 7.5 + 16*i;
    ya[0] = -0.5;
    ya[1] = 15.5;

    first = 0;
  }
  
  for (size_t i=0; i<32; i++) {
    gsl_spline2d_set(spline, dbl_data, i, 0, grid_values[i]);
    gsl_spline2d_set(spline, dbl_data, i, 1, grid_values[i+32]);
  }
    
  gsl_spline2d_init(spline, xa, ya, dbl_data, 32, 2);

  // -111 = Fill above data
  // -222 = Fill below data
  // -333 = Fill IPI invalid
  // -444 = Fill to side of data
  // -555 = Fill not processed
  // -999 = Fill IPI error
  
  for (size_t j=0; j<16; ++j) {
    double yi = (double) j;
    for (size_t i=8; i<504; ++i) {
      double xi = (double) i;

      if (grid_values[(i-8)/16] == -111 || grid_values[(i-8)/16] == -222 ||
          grid_values[(i-8)/16] == -333 || grid_values[(i-8)/16] == -444 ||
          grid_values[(i-8)/16] == -555 || grid_values[(i-8)/16] == -999 ||

          grid_values[(i-8)/16+1] == -111 || grid_values[(i-8)/16+1] == -222 ||
          grid_values[(i-8)/16+1] == -333 || grid_values[(i-8)/16+1] == -444 ||
          grid_values[(i-8)/16+1] == -555 || grid_values[(i-8)/16+1] == -999 ||

          grid_values[(i-8)/16+32] == -111 || grid_values[(i-8)/16+32] == -222 ||
          grid_values[(i-8)/16+32] == -333 || grid_values[(i-8)/16+32] == -444 ||
          grid_values[(i-8)/16+32] == -555 || grid_values[(i-8)/16+32] == -999 ||

          grid_values[(i-8)/16+32+1] == -111 || grid_values[(i-8)/16+32+1] == -222 ||
          grid_values[(i-8)/16+32+1] == -333 || grid_values[(i-8)/16+32+1] == -444 ||
          grid_values[(i-8)/16+32+1] == -555 || grid_values[(i-8)/16+32+1] == -999)
        
        interpolated_values[j][i] = -555;
      else
        interpolated_values[j][i] =
          gsl_spline2d_eval(spline, xi, yi, xacc, yacc);


    }
  }

  return 0;
}


int interp_values_dbl( int32_t diff_offset, double *grid_values,
                       double interpolated_values[16][512]) {

  double dbl_data[2*32];
  
  if (diff_offset == 0) {
    for (size_t i=0; i<32; i++) {
      gsl_spline2d_set(spline_misr, dbl_data, i, 0, (double) grid_values[i]);
      gsl_spline2d_set(spline_misr, dbl_data, i, 1, (double) grid_values[i+32]);
    }

    gsl_spline2d_init(spline_misr, xa, ya, dbl_data, 32, 2);

    for (size_t j=0; j<16; ++j) {
      double yi = (double) j;
      for (size_t i=8; i<504; ++i) {
        double xi = (double) i;
        //printf("%d %f %d %f\n", j, yi, i, xi);
        interpolated_values[j][i] =
          gsl_spline2d_eval(spline_misr, xi, yi, xacc, yacc);
      }
    }
  } else if (diff_offset == -1) {

    for (size_t i=0; i<32; i++)
      gsl_spline2d_set(spline_misr, dbl_data, i, 0, (double) grid_values[i]);
    for (size_t i=0; i<31; i++)
      gsl_spline2d_set(spline_misr, dbl_data, i+1, 1,
                       (double) grid_values[i+32]);

    for (size_t j=0; j<8; ++j) {
      double yi = (double) j;
      for (size_t i=8; i<504; ++i) {
        double xi = (double) i;
        //printf("%d %f %d %f\n", j, yi, i, xi);
        interpolated_values[j][i] =
          gsl_spline2d_eval(spline_misr, xi, yi, xacc, yacc);
      }
    }

    for (size_t i=0; i<32; i++)
      gsl_spline2d_set(spline_misr, dbl_data, i, 0, (double) grid_values[i]);
    for (size_t i=0; i<32; i++)
      gsl_spline2d_set(spline_misr, dbl_data, i, 1,
                       (double) grid_values[i+32]);
      
  } else {
    for (size_t i=1; i<32; i++)
      gsl_spline2d_set(spline_misr, dbl_data, i-1, 1,
                       (double) grid_values[i]);
    
    for (size_t j=0; j<16; ++j) {
      double yi = (double) j;
      for (size_t i=8; i<504; ++i) {
        double xi = (double) i;
        //printf("%d %f %d %f\n", j, yi, i, xi);
        interpolated_values[j][i] =
          gsl_spline2d_eval(spline_misr, xi, yi, xacc, yacc);
      }
    }

  }

  return 0;
}


int reduce_res(uint16_t rad_data[4][2048]) {

  static int first=1;
  static gsl_matrix *X, *cov;
  static gsl_vector *y, *w, *c;

  static gsl_multifit_linear_workspace *work;

  double chisq;
  
  if (first) {

    X = gsl_matrix_alloc (16, 3);
    y = gsl_vector_alloc (16);
    w = gsl_vector_alloc (16);

    c = gsl_vector_alloc (3);
    cov = gsl_matrix_alloc (3, 3);

    work = gsl_multifit_linear_alloc (16, 3);
    
    for (size_t i=0; i<16; i++) {
      gsl_matrix_set (X, i, 0, 1.0);

      int row, col;
      row = i / 4;
      col = i - row*4;

      double x_val = col - 1.5;
      double y_val = 1.5 - row;
      
      gsl_matrix_set (X, i, 1, x_val);
      gsl_matrix_set (X, i, 2, y_val);

      first = 0;
    }
  }

  for (size_t ip=0; ip<512; ip++) {
    int ngood = 0;
    for (size_t irow=0; irow<4; irow++) {
      for (size_t icol=0; icol<4; icol++) {
        int index = irow*4+icol;
        double d_val = (double) rad_data[irow][4*ip+icol];
        gsl_vector_set (y, index, d_val);
        if (rad_data[irow][4*ip+icol] >= 16378) {
          gsl_vector_set (w, index, 0.0);
        } else {
          gsl_vector_set (w, index, 1.0);
          ngood++;
        }
      }
    }
    if (ngood >= 5) {

      /*
      double *ptr = gsl_matrix_ptr(X,0,0);
      for (size_t i=0; i<16; i++)
        printf("%lf %lf %lf\n", ptr[i*3+0],ptr[i*3+1],ptr[i*3+2]);
      printf("\n");

      ptr = gsl_vector_ptr(w,0);
      for (size_t i=0; i<16; i++) printf("%lf\n", ptr[i]);
      printf("\n");

      ptr = gsl_vector_ptr(y,0);
      for (size_t i=0; i<16; i++) printf("%lf\n", ptr[i]);
      printf("\n");
      */
      gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
      double fit_val = gsl_vector_get(c,0);
      rad_data[0][ip] = (ushort) ( fit_val + 0.5);
    } else {
      rad_data[0][ip] = 16378;
    }
  }  // i loop

  return 0;
}


int closel1_misr(filehandle *file) {

  int i;
  
  SDendaccess(lat_id);
  SDendaccess(lon_id);

  for (i=0; i<N_CAMERAS; i++) {
    if (sd_id[i] == -1) continue;
    SDendaccess(blu_id[i]);
    SDendaccess(grn_id[i]);
    SDendaccess(red_id[i]);
    SDendaccess(nir_id[i]);

    SDend(sd_id[i]);
  }
  //    VSdetach(vd_id);
  //Vend(file_id);
  //Hclose(file_id);
  
  return 0;
}


