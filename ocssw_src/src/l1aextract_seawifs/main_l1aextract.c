#include "regen.h"
#include "regen_proto.h"

int main(int32 argc, char *argv[]) {

    char *infile, *outfile;
    int32 status, spix, epix, sscan, escan, pix_sub, sc_sub;

    /*** check Usage ***/

    if (argc != 9) {
        system("clear");
        printf("\n\n\nUsage: <%s> ", argv[0]);
        printf("infile spix epix sscan escan pix_sub sc_sub outfile"
                "\n   where:"
                "\n\tinfile   - input l1a data HDF file"
                "\n\tspix     - start pixel number"
                "\n\tepix     - end pixel number"
                "\n\tsscan    - start scan line"
                "\n\tescan    - end scan line"
                "\n\tpix_sub  - pixel subsampling rate"
                "\n\tsc_sub   - scan line subsampling rate"
                "\n\toutfile  - output file name"
                "\n\nNOTE: Start and End pixel/scanline will be set to input file's"
                " nsamp\n\t and nrec respectively, if the given values are out"
                " of range.\n\n");
        exit(EXIT_FAILURE);
    }

    /*** load input parameters into local variables */

    infile = argv[1];
    outfile = argv[8];
    spix = atoi(argv[2]);
    epix = atoi(argv[3]);
    sscan = atoi(argv[4]);
    escan = atoi(argv[5]);
    pix_sub = atoi(argv[6]);
    sc_sub = atoi(argv[7]);

    /*** call regen */

    status = regen(infile, &spix, &epix, &sscan, &escan, pix_sub, sc_sub,
            "ALL", outfile);

    if (status != EXIT_SUCCESS) {
        printf("\nError: status returned = %d", status);
        printf("\n%s\n", ERR_MSG);
        exit(EXIT_FAILURE);
    }

    printf("\nThe output file %s has been successfully generated\n", outfile);

    return EXIT_SUCCESS;
}
