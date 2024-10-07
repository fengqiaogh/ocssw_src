/*
    Modification history:
    Programmer       Organization      Date      Description of change
    --------------   ------------    --------    ---------------------
    Joel Gales       Futuretech      058/06/01    Original Development

 */


#include <stdio.h>
#include <math.h>
#include "hdf.h"

unsigned char read_mask(double lon, double lat, int16 close) {

    int16 nrec;
    int16 header[1024];
    int16 ptr_arr[1024 * 64];
    static int16 cmpz[8 * 128];
    int32 offset = -1;

    int32 ilon, ilat, i, j;

    unsigned char mask[128 * 128];
    unsigned char zero = 0;
    unsigned char one = 1;
    unsigned char mask_lonlat = 0;

    static int32 first = 1;
    static unsigned int pow2[16] = {1, 2, 4, 8, 16, 32, 64, 128, 256,
        512, 1024, 2048, 4096, 8192, 16384, 32768};
    static FILE *maskfile;

    if (close == 1) {
        fclose(maskfile);
        return -one;
    }

    if (first == 1) {
        maskfile = fopen("/home/simbiosd/data/common/landmask_modland.dat", "r");
        fread(header, sizeof (int16), 1024, maskfile);
        nrec = header[1];

        fseek(maskfile, sizeof (header), SEEK_SET);
        fread(ptr_arr, sizeof (int16), 1024 * 64, maskfile);
        first = 0;
    }

    ilon = floor(lon) + 180;
    ilat = floor(lat) + 90;

    if (ptr_arr[360 * ilat + ilon] == 1) mask_lonlat = one;

    if (ptr_arr[360 * ilat + ilon] > 65) {
        if (sizeof (int16)*1024 * ptr_arr[360 * ilat + ilon] != offset) {
            offset = sizeof (int16)*1024 * ptr_arr[360 * ilat + ilon];

            fseek(maskfile, offset, SEEK_SET);
            fread(cmpz, sizeof (int16), 8 * 128, maskfile);
        }

        i = (lon - (ilon - 180)) * 128;
        j = (lat - (ilat - 90)) * 128;

        if (i < 0) i = 128 + i;
        if (j < 0) j = 128 + j;

        if ((cmpz[j * 8 + (i / 16)] & pow2[15 - (i % 16)]) != 0) mask_lonlat = one;
        else mask_lonlat = zero;
    }

    return mask_lonlat;
}

#if 1

/* 
cc -o read_mask read_mask.c -g -I/systems04/home/simbiosd/lib/hdf/include -lm 

 */

main() {
    int i, j;
    double inlon, inlat;
    unsigned char mask_lonlat, map[180 * 360];

    FILE *fp;

    /*

      for (i=-90; i<90; i++) {
        for (j=-180; j<180; j++) {

          mask_lonlat=read_mask((double) (j*0.02+10), (double) (i*0.02-2), (int16) 0);
          map[(i+90) * 360 + (j+180)] = mask_lonlat;
        }
      }
     */

loop:

    printf("Input lon lat: ");
    scanf("%lf %lf", &inlon, &inlat);
    if (inlon >= -180) {
        mask_lonlat = read_mask(inlon, inlat, (int16) 0);
        printf("%lf %lf %d\n", inlon, inlat, mask_lonlat);
        goto loop;
    }

    read_mask(inlon, inlat, (int16) 1);

    fp = fopen("map.dat", "wb");
    fwrite(map, 1, 180 * 360, fp);
    close(fp);

    /* map=bytarr(360,180)&openr,1,'map.dat'&readu,1,map&close,1&tv,bytscl(map,min=0,max=1) */

}


#endif

