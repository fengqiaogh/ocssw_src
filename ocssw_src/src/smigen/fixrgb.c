#include <stdio.h>

/*
  This routine zeros any pixels with zero values in ANY of the RGB pnm file.
  It removes the blue/pink stripes that can occur in the SEAWIFS TC images.
 */

int main(int argc, char *argv[]) {
    unsigned char c[3];
    unsigned char zero = 0;
    FILE * fp[3];
    int n, m;

    fp[0] = fopen(argv[1], "rb+");
    n = 0;
    m = 0;
    while (1) {
        fread(&c[0], 1, 1, fp[0]);
        m++;
        if (c[0] == 0x0a) n++;
        if (n == 4) break;
    }

    fp[1] = fopen(argv[2], "rb+");
    fseek(fp[1], m, SEEK_SET);

    fp[2] = fopen(argv[3], "rb+");
    fseek(fp[2], m, SEEK_SET);

    while (1) {
        fread(&c[0], 1, 1, fp[0]);
        if (feof(fp[0]) != 0) break;
        fread(&c[1], 1, 1, fp[1]);
        fread(&c[2], 1, 1, fp[2]);

        if (c[0] == 0 && c[1] == 0 && c[2] == 0) continue;

        if (c[0] == 0 || c[1] == 0 || c[2] == 0) {
            fseek(fp[0], -1, SEEK_CUR);
            fseek(fp[1], -1, SEEK_CUR);
            fseek(fp[2], -1, SEEK_CUR);

            fwrite(&zero, 1, 1, fp[0]);
            fwrite(&zero, 1, 1, fp[1]);
            fwrite(&zero, 1, 1, fp[2]);
        }
    }

    fclose(fp[0]);
    fclose(fp[1]);
    fclose(fp[2]);

    return 0;
}
