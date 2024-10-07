/* ========================================================================
 * MSl1tcpcbox - converts any L1B file to a ppm file
 * 
 * Synopsis:
 *
 *   MSl1tcpcbox  input-l1[a|b]-filename output-ppm-filename
 *
 * Description:
 * 
 * Modification history:
 *
 *     Programmer     Organization      Date       Description of change
 *   --------------   ------------    --------     ---------------------
 *   Bryan A. Franz   GSC             Jul 1998     Original development
 *   Joel Gales	      Futuretech      Jan 2001     Added 8 and 24-bit hdf
 *   Norman Kuring    NASA/GSFC       Feb 2005     New scaling for true color
 *   Sean Bailey      Futuretech      Apr 2005	   Modified MSl1brsgen to create
 *					           MSl1tcpcbox by incorporating
 *						   code from swtcpcbox
 *   Norman Kuring    NASA/GSFC       Mar 2006     New, more accurate estimate
 *						   of pixel boundaries
 *   Norman Kuring    NASA/GSFC       Mar 2006     Changed binsoverlain from
 *						   a static array to a
 *						   dynamically allocated one
 *						   to get rid of that pesky
 *						   MAXBINSPERPIXEL warning.
 *   Norman Kuring    NASA/GSFC       Jul 2006     Removed sawtooth edges from
 *						   output image and fixed bug
 *						   that caused edge pixels to
 *						   wrap from right to left side
 *						   and vice versa.
 *
 *   Sean Bailey      Futuretech      Nov 2007     Made scaling consistent with
 *                          l1brsgen, use a sensor dependent defaults, and
 *                          generally clean up the code
 *   Wayne Robinson   SAIC            Sep 2010     correct delta long 
 *                                               computation for dateline wrap
 * ======================================================================== */
#include "l12_proto.h"
#include "input_struc.h"

#include "miscstruct.h"
#include "miscanfill.h"
#include "mipoly.h"
#include "clo.h"
#include <X11/X.h>
#include <X11/Xlib.h>
#include <png.h>

#include <geotiffio.h>
#include <geo_normalize.h>
#include <geo_tiffp.h>
#include <geo_keyp.h>
#include <xtiffio.h>
#include <geokeys.h>
#include <scene_meta.h>
#include <libnav.h>

//int l1mapgen_usage(char *prog);

#define ROOT2 1.4142135623730950488016887242096981

#define INT32   int32_t 
#define FLOAT32 float
#define BYTE    unsigned char

#define BINBELOWTHRESH  110

#define MALLOC(ptr,typ,num) {      \
  (ptr) = (typ *)malloc((num) * sizeof(typ));    \
  if((ptr) == NULL){       \
    fprintf(stderr,"-E- %s line %d: Memory allocation failure.\n", \
    __FILE__,__LINE__);       \
    exit(EXIT_FAILURE);       \
  }         \
}
#define CALLOC(ptr,typ,num) {                                           \
  (ptr) = (typ *)calloc((num) , sizeof(typ));                           \
  if((ptr) == NULL){                                                    \
    fprintf(stderr,"-E- %s line %d: Memory allocation failure.\n",      \
    __FILE__,__LINE__);                                                 \
    exit(EXIT_FAILURE);                                                 \
  }                                                                     \
}

BYTE logscale(float val, float min, float max);
BYTE linscale(float val, float min, float max);

double rint(double);
int scan_convert(XPoint *);
int collect_bins(int, XPoint *, int *);

static int *binsoverlain = NULL;
static int binsperpixel = 0;
static int numoverlain;
static int width, height;
static double imscale;

/* -------------------------------------------------------------------- *
 *                              main                                    *
 * -------------------------------------------------------------------- */
int main(int argc, char *argv[]) {

    int32_t iscan = 0; /* input scan number                  */
    int32_t spix = 0; /* start pixel for subscene process   */
    int32_t epix = -1; /* end pixel for subscene process     */
    int32_t sscan = 0; /* start scan for subscene process    */
    int32_t escan = -1; /* end scan for subscene process      */
    int32_t dpix = 1; /* pixel increment for sub-sampling   */
    int32_t dscan = 1; /* scan subsampling increment         */

    int32_t ip; /* input pixel number                 */

    int32_t npix = 0; /* Number of output pixels per scan   */
    int32_t nscan = 0; /* Number of output scans             */

    BYTE rgb[3];
    static int32_t r;
    static int32_t g;
    static int32_t b;

    int want_linscale = 0;
    int want_atmocor = 0;
    float min = 0.01;
    float max = 0.9;

    double esdist;
    int32_t sfactor = 1;
    int32_t iw;

    l1str l1rec; /* generic level-1b scan structure      */
    filehandle l1file; /* input file handle                    */
    FILE *outfp = NULL;

    char outfile[FILENAME_MAX] = "\0";

    int32_t bb;

    unsigned char *mean;
    int32_t *count;
    int32_t progress;
    float64 * sum[3];
    uint32_t numbins, bin, numoutpix;
    XPoint corners[8]; /* actually, 4 pixel corners + 4 side midpoints */
    float64 threshold, outpixpercent;
    float32 north, south, west, east;
    float32 nrad, srad, wrad, erad;

    float32 sr_r, sr_g, sr_b;

    int i;

    /* hold all of the command line options */
    clo_optionList_t* list;

    float toa_reflect(l1str * l1rec, int32_t ip, int32_t ib);

    if (argc == 1) {
        l2gen_usage("l1mapgen");
        return 0;
    }

    // see if help on command line
    for (i = 0; i < argc; i++) {
        if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-help") == 0)) {
            l2gen_usage("l1mapgen");
            return 1;
        }
    }

    cdata_();
    filehandle_init(&l1file);
    msl12_input_init();

    // create an option list
    list = clo_createList();

    /* initialize the option list with descriptions and default values */
    l2gen_init_options(list, "l1mapgen");

    // change the default values for this program
    clo_setString(list, "atmocor", "off", "default");
    clo_setString(list, "proc_land", "1", "default");
    clo_setString(list, "sl_pixl", "0", "default");

    /* Parse input parameters */
    if (msl12_option_input(argc, argv, list, "l1mapgen", &l1file) != 0) {
        printf("-E- %s: Error parsing input parameters.\n", argv[0]);
        exit(FATAL_ERROR);
    }

    strcpy(outfile, input->ofile[0]);

    width = input->width;
    threshold = input->threshold;
    want_linscale = input->stype;
    if (input->atmocor == 1)
        want_atmocor = 1;

    min = input->datamin;
    max = input->datamax;

    north = input->north;
    south = input->south;
    east = input->east;
    west = input->west;

    /*                                                                  */
    /* Open input file and get sensor and scan information from handle  */
    /*                                                                  */
    if (openl1(&l1file) != 0) {
        printf("-E- %s: Error opening %s for reading.\n", argv[0],
                l1file.name);
        exit(FATAL_ERROR);
    }

    r = bindex_get(input->rgb[0]);
    g = bindex_get(input->rgb[1]);
    b = bindex_get(input->rgb[2]);
    if (r < 0 || g < 0 || b < 0) {
        printf("-E- Invalid RGB set for this sensor\n");
        if (r < 0)
            printf("   Invalid R: %d\n", input->rgb[0]);
        if (g < 0)
            printf("   Invalid G: %d\n", input->rgb[1]);
        if (b < 0)
            printf("   Invalid B: %d\n", input->rgb[2]);
        exit(FATAL_ERROR);
    }


    /* Allocate memory for L1 scan data */

    if (alloc_l1(&l1file, &l1rec) == 0) {
        printf("-E- %s: Unable to allocate L1 record.\n", argv[0]);
        exit(FATAL_ERROR);
    }

    /* Set the end pixel if it was not set by command argument	        */
    if (l1_input->epixl < 1 || l1_input->epixl > l1file.npix)
        l1_input->epixl = l1file.npix;
    if (l1_input->eline < 1 || l1_input->eline > l1file.nscan)
        l1_input->eline = l1file.nscan;
    if (l1_input->spixl < 1)
        l1_input->spixl = 1;
    if (l1_input->sline < 1)
        l1_input->sline = 1;
    sfactor = 1;
    spix = MAX(l1_input->spixl - 1, 0);
    epix = MIN(l1_input->epixl - 1, l1file.npix - 1);
    sscan = MAX(l1_input->sline - 1, 0);
    escan = MIN(l1_input->eline - 1, l1file.nscan - 1);
    dpix = sfactor;
    dscan = sfactor;

    npix = (epix - spix) / sfactor + 1;
    nscan = (escan - sscan) / sfactor + 1;

    l1file.spix = spix;
    l1file.epix = epix;

    l1_input->sline = sscan + 1;
    l1_input->eline = escan + 1;
    l1_input->dline = dscan;
    l1_input->spixl = spix + 1;
    l1_input->epixl = epix + 1;
    l1_input->dpixl = dpix;

    /* Determin scene boundaries */
    if (north == -999 || south == -999 || east == -999 || south == -999) {
        fprintf(stderr,
                "Not all boundaries defined, grabbing from file...patience please...\n");

        // add a fix for seawifs, can not read seawifs lines out of order.
        if (l1file.sensorID == SEAWIFS) {
            int32_t sd_id = SDstart(l1file.name, DFACC_RDONLY);
            if (sd_id == -1) {
                printf("-E- lonlat2pixline: Error opening %s for reading.\n",
                        l1file.name);
                exit(1);
            }
            float tmpf;
            if (north == -999) {
                READ_GLBL_ATTR_E("Northernmost Latitude", &tmpf);
                north = ceilf(tmpf);
            }
            if (south == -999) {
                READ_GLBL_ATTR_E("Southernmost Latitude", &tmpf);
                south = floorf(tmpf);
            }
            if (east == -999) {
                READ_GLBL_ATTR_E("Easternmost Longitude", &tmpf);
                east = ceilf(tmpf);
            }
            if (west == -999) {
                READ_GLBL_ATTR_E("Westernmost Longitude", &tmpf);
                west = floorf(tmpf);
            }
            SDend(sd_id);

        } else {
            // not a sweawifs L1 file

            scnstr *meta = scene_meta_get();

            for (iscan = sscan; iscan <= escan; iscan += 10) {

                if (readl1_lonlat(&l1file, iscan, &l1rec) != 0) {
                    fprintf(stderr, "-E- %s Line %d: error reading %s at scan %d.\n",
                            __FILE__, __LINE__, l1file.name, iscan);
                    exit(FATAL_ERROR);
                }

                scene_meta_put(&l1rec);
            }

            if (north == -999)
                north = ceilf(meta->northern_lat);
            if (south == -999)
                south = floorf(meta->southern_lat);
            if (east == -999)
                east = ceilf(meta->eastern_lon);
            if (west == -999)
                west = floorf(meta->western_lon);

            if (strcmp(meta->start_node, meta->end_node) != 0) {
                printf("-E- Crossing the pole! Won't continue.\n");
                closel1(&l1file);
                exit(FATAL_ERROR);
            }

        } // not seawifs


    }
    fprintf(stderr, "N: %f S: %f E: %f W:%f\n", north, south, east, west);
    nrad = north * PI / 180;
    srad = south * PI / 180;
    wrad = west * PI / 180;
    erad = east * PI / 180;

    if (nrad <= srad) {
        fprintf(stderr,
                "The northernmost boundary must be greater than the ");
        fprintf(stderr, "southernmost boundary.\n");
        exit(FATAL_ERROR);
    }
    /* Get the size of the output image. */
    if (width < 32) {
        fprintf(stderr,
                "Width (%d) is too small to produce a useful image.\n",
                width);
        l2gen_usage("l1mapgen");
        exit(FATAL_ERROR);
    }
    if (wrad < erad)
        height = rint((nrad - srad) * width / (erad - wrad));
    else
        height = rint((nrad - srad) * width / (erad - wrad + 2 * PI));

    numbins = width * height;
    imscale = height / (nrad - srad);

    for (bb = 0; bb < 3; bb++) {
        CALLOC(sum[bb], float64, numbins);
    }
    CALLOC(count, int32_t, numbins);
    CALLOC(mean, unsigned char, numbins * 3);

    /* Use the following variable to show this program's progress. */
    progress = (int) ceil(((double) nscan / 78));


    printf("Using r,g,b = %d, %d, %d\n", input->rgb[0], input->rgb[1],
            input->rgb[2]);

    /* Read file scan by scan, scale radiances, and write. */
    for (iscan = sscan; iscan <= escan; iscan++) {

        if (iscan % progress == 0)
            fprintf(stderr, ".");
        if (readl1(&l1file, iscan, &l1rec) != 0) {
            fprintf(stderr,
                    "-E- %s Line %d: error reading %s at scan %d.\n",
                    __FILE__, __LINE__, l1file.name, iscan);
            exit(FATAL_ERROR);
        }

        if (want_atmocor) {
            if (loadl1(&l1file, &l1rec) != 0) {
                fprintf(stderr,
                        "-E- %s Line %d: error loading %s at scan %d.\n",
                        __FILE__, __LINE__, l1file.name, iscan);
                exit(FATAL_ERROR);
            }
        } else {
            /* Get correction for Earth-Sun distance and apply to Fo  */
            int16_t year, day;
            double dsec;
            unix2yds(l1rec.scantime, &year, &day, &dsec);
            int32_t msec = (int32_t) (dsec * 1.e3);
            int32_t yr = (int32_t) year;
            int32_t dy = (int32_t) day;
            esdist = esdist_(&yr, &dy, &msec);
            l1rec.fsol = pow(1.0 / esdist, 2);
            for (iw = 0; iw < l1file.nbands; iw++) {
                l1rec.Fo[iw] = l1file.Fobar[iw] * l1rec.fsol;
            }
        }


        /*  Equations used to compute pixel boundaries below are taken from:
                Map Projections -- A Working Manual
                by John P. Snyder
                U.S. Geological Survey Professional Paper 1395
                Fourth Printing 1997
         */

        for (ip = 0; ip < npix; ip++) {

            double plat, plon; /* pixel center */
            double sinplat, cosplat;
            double c[2]; /* edge/corner distance */
            double sinc[2], cosc[2];
            double sinaz[8], cosaz[8];
            int i; /* loop index */
            int out = 0; /* Number of pixel "corners"   */

            /* that fall outside the image */

            plat = l1rec.lat[ip] * PI / 180;
            plon = l1rec.lon[ip] * PI / 180;

            sinplat = sin(plat);
            cosplat = cos(plat);

            if (ip < npix - 1) {

                double nlat, nlon; /* next pixel center */
                double dlat, dlon; /* delta lat and lon */
                double sindlat2, sindlon2, cosnlat;
                double az; /* azimuth to next pixel */
                double phw; /* pixel half width */

                nlat = l1rec.lat[ip + 1] * PI / 180;
                nlon = l1rec.lon[ip + 1] * PI / 180;

                cosnlat = cos(nlat);

                dlat = nlat - plat;
                dlon = nlon - plon;
                /*
                 * WDR if del long is big in radians, assume dateline wrap
                 */
                if (dlon > 5.) dlon = nlon - plon - 2 * PI;
                if (dlon < -5.) dlon = nlon + 2 * PI - plon;

                sindlat2 = sin(dlat / 2);
                sindlon2 = sin(dlon / 2);

                phw = asin(sqrt(sindlat2 * sindlat2 + /* Page 30, */
                        sindlon2 * sindlon2 * /* Equation */
                        cosplat * cosnlat)); /*     5-3a */

                /*
                   Make the pixel's coverage slightly broader to avoid
                   black pin holes in the output image where the estimated
                   pixel boundaries do not quite line up.  Don't increase
                   the size too much, however, because that will cause
                   unwanted image blurring.
                 */
                phw *= 1.6;

                c[0] = phw; /* center-to-edge   distance */
                c[1] = phw * ROOT2; /* center-to-corner distance */

                sinc[0] = sin(c[0]);
                sinc[1] = sin(c[1]);
                cosc[0] = cos(c[0]);
                cosc[1] = cos(c[1]);

                az = /* Page 30, Equation 5-4b */
                        atan2(cosnlat * sin(dlon),
                        cosplat * sin(nlat) -
                        sinplat * cosnlat * cos(dlon)
                        );

                sinaz[0] = sin(az);
                cosaz[0] = cos(az);
                for (i = 1; i < 8; i++) {
                    az += PI / 4;
                    sinaz[i] = sin(az);
                    cosaz[i] = cos(az);
                }
            }

            for (i = 0; i < 8; i++) {
                double lat, lon;
                int ec; /* edge = 0, corner = 1 */

                ec = i & 1;

                /* Page 31, Equation 5-5 */
                lat = asin(sinplat * cosc[ec]
                        + cosplat * sinc[ec] * cosaz[i]);

                lon = plon /* Page 31, Equation 5-6 */
                        + atan2(sinc[ec] * sinaz[i], cosplat * cosc[ec]
                        - sinplat * sinc[ec] * cosaz[i]);

                if (wrad > erad) { /* 180-degree meridian included */
                    if (lon >= wrad && lon > erad) {
                        corners[i].x = (lon - wrad) * imscale;
                    } else if (lon < wrad && lon <= erad) {
                        corners[i].x = (lon - wrad + 2 * PI) * imscale;
                    } else if (lon - erad < wrad - lon) {
                        corners[i].x = (lon - wrad + 2 * PI) * imscale;
                    } else {
                        corners[i].x = (lon - wrad) * imscale;
                    }
                } else { /* 180-degree meridian not included */
                    corners[i].x = (lon - wrad) * imscale;
                }

                corners[i].y = (nrad - lat) * imscale;

                if (corners[i].x < 0 || corners[i].x >= width
                        || corners[i].y < 0 || corners[i].y >= height)
                    out++;
            }

            /*
               If out == 8, then the entire pixel is outside the
               image area, so I skip to the next one.
             */
            if (out == 8)
                continue;

            if (l1rec.Lt[ip * l1rec.l1file->nbands + r] <= 0.01)
                rgb[0] = rgb[1] = rgb[2] = 255;
            else if (want_atmocor) {
                sr_r = l1rec.rhos[ip * l1rec.l1file->nbands + r];
                sr_g = l1rec.rhos[ip * l1rec.l1file->nbands + g];
                sr_b = l1rec.rhos[ip * l1rec.l1file->nbands + b];

                /* Cheat blue band for MODIS */
                /*
                   if ( (b == 2)
                   &&   (l1file.sensorID == 5 || l1file.sensorID == 6) )
                   sr_b *= 0.8;
                 */

                /* Cheat if red band saturates */
                /* if (sr_r > 1.) sr_r = 0.975 * sr_g; */

            } else {
                sr_r = toa_reflect(&l1rec, ip, r);
                sr_g = toa_reflect(&l1rec, ip, g);
                sr_b = toa_reflect(&l1rec, ip, b);
            }

            if (want_linscale == 1) {
                rgb[0] = linscale(sr_r, min, max);
                rgb[1] = linscale(sr_g, min, max);
                rgb[2] = linscale(sr_b, min, max);
            } else {
                rgb[0] = logscale(sr_r, min, max);
                rgb[1] = logscale(sr_g, min, max);
                rgb[2] = logscale(sr_b, min, max);
            }

            if (l1file.sensorID == OCTS && rgb[2] > 234) {
                if (rgb[0] < 230)
                    rgb[0] = 250;
                if (rgb[1] < 230)
                    rgb[1] = 251;
                rgb[2] = 252;
            }

            /*
               Use scan conversion to determine which bins are overlain
               by this pixel.  The global array, binsoverlain, is populated
               by the scan_convert function.
             */
            numoverlain = 0;
            scan_convert(corners);

            for (i = 0; i < numoverlain; i++) {
                bin = binsoverlain[i];
                /* Accumulate the sums for each bin. */
                for (bb = 0; bb < 3; bb++) {
                    sum[bb][bin] += rgb[bb];
                }
                count[bin]++;
            }

        } /* End for pixel */
    } /* End for scan */

    fprintf(stderr, "\n");

    closel1(&l1file);

    /* Calculate the means for each bin. */
    numoutpix = 0;
    for (bin = 0; bin < numbins; bin++) {
        if (count[bin] == 0) {
            mean[3 * bin] = mean[3 * bin + 1] = mean[3 * bin + 2] = 0;
            continue;
        }
        mean[3 * bin] = sum[0][bin] / count[bin];
        mean[3 * bin + 1] = sum[1][bin] / count[bin];
        mean[3 * bin + 2] = sum[2][bin] / count[bin];
        if ((mean[3 * bin] + mean[3 * bin + 1] + mean[3 * bin + 2]) == 0) {
            mean[3 * bin] = mean[3 * bin + 1] = mean[3 * bin + 2] = 1;
        }
        numoutpix++;
    }

    /* Don't generate an image for underwhelming results. */
    outpixpercent = (double) 100 * numoutpix / numbins;

    if (outpixpercent < threshold) {
        fprintf(stderr, "Number of output pixels (%u of a possible %u) ",
                numoutpix, numbins);
        fprintf(stderr,
                "< threshold (%.2lf%%). Output image not generated.\n",
                threshold);
        exit(BINBELOWTHRESH);
    }


    // write out the different file formats
    // PPM file
    if (strcmp(input->oformat, "PPM") == 0) {
        if ((outfp = fopen(outfile, "w")) == NULL) {
            printf("%s: Error: Unable to open %s for writing.\n",
                    argv[0], outfile);
            exit(FATAL_ERROR);
        }
        fprintf(outfp, "P6\n");
        fprintf(outfp, "%d\n", width);
        fprintf(outfp, "%d\n", height);
        fprintf(outfp, "255\n");
        fwrite(mean, 3, numbins, outfp);
        fclose(outfp);

        // PNG file
    } else if (strcmp(input->oformat, "PNG") == 0) {
        if ((outfp = fopen(outfile, "w")) == NULL) {
            printf("%s: Error: Unable to open %s for writing.\n",
                    argv[0], outfile);
            exit(FATAL_ERROR);
        }

        png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                NULL, NULL, NULL);
        if (!png_ptr) {
            fprintf(stderr, "%s: Error: Unable to create PNG write structure.\n", argv[0]);
            exit(FATAL_ERROR);
        }

        png_infop info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr) {
            png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
            fprintf(stderr, "%s: Error: Unable to create PNG info structure.\n", argv[0]);
            exit(FATAL_ERROR);
        }
        if (setjmp(png_jmpbuf(png_ptr))) {
            png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
            fprintf(stderr, "%s: Error: Unable to call PNG setjmp().\n", argv[0]);
            exit(FATAL_ERROR);
        }
        png_init_io(png_ptr, outfp);

        uint8_t * row_pointers[height];
        for (i = 0; i < height; i++)
            row_pointers[i] = mean + (i * width * 3);
        png_set_rows(png_ptr, info_ptr, row_pointers);

        png_set_IHDR(png_ptr, info_ptr, width, height,
                8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

        png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

        /* clean up after the write, and free any memory allocated */
        png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
        fclose(outfp);

        // geoTIFF file
    } else if (strcmp(input->oformat, "TIFF") == 0) {
        TIFF *tiff;
        GTIF *gtif;

        tiff = XTIFFOpen(outfile, "w");
        if (tiff == NULL) {
            fprintf(stderr, "Could not open outgoing image\n");
            exit(FATAL_ERROR);
        }
        gtif = GTIFNew(tiff);
        if (gtif == NULL) {
            fprintf(stderr, "Could not create geoTIFF structure\n");
            exit(FATAL_ERROR);
        }

        // Write the tiff tags to the file
        TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, width);
        TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, height);
        //TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
        TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
        TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);
        TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 3);

        // write geo TIFF  tags
        double tiepoints[6] = {0, 0, 0, 0, 0, 0};
        double pixscale[3] = {0, 0, 0};

        // pixel width
        pixscale[0] = (east - west) / width;

        // pixel height
        pixscale[1] = (north - south) / height;

        // set top left corner pixel lat, lon
        tiepoints[3] = west;
        tiepoints[4] = north;

        TIFFSetField(tiff, GTIFF_PIXELSCALE, 3, pixscale);
        TIFFSetField(tiff, GTIFF_TIEPOINTS, 6, tiepoints);

        // write geo TIFF keys
        GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelGeographic);
        GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);
        GTIFKeySet(gtif, GeographicTypeGeoKey, TYPE_SHORT, 1, GCS_WGS_84);

        GTIFWriteKeys(gtif);

        // Actually write the image
        if (TIFFWriteEncodedStrip(tiff, 0, mean, width * height * 3) == 0) {
            fprintf(stderr, "Could not write TIFF image\n");
            exit(FATAL_ERROR);
        }

        GTIFFree(gtif);
        XTIFFClose(tiff);


    } else {
        fprintf(stderr, "%s Error: oformat=%s not valid.\n", argv[0], input->oformat);
        exit(FATAL_ERROR);
    }

    exit(EXIT_SUCCESS);
}

/*------------------------------------------------------------------------------

   Function: scan_convert

   Returns type: int

   Description:	This function is derived from the scan-conversion
                code used by an X-server for filling in polygons.

   Parameters: (in calling order)
      Type		Name		I/O	Description
      ----		----		---	-----------
      struct _DDXPoint*	ptsIn		In	vertex specifications

   Modification history:
      Programmer	Date		Description of change
      ----------	----		---------------------
      Norman Kuring	15-Sep-1992	Adaptation from X-server code
                                        for my own nefarious purposes.

------------------------------------------------------------------------------*/

/***********************************************************
Copyright 1987 by Digital Equipment Corporation, Maynard, Massachusetts,
and the Massachusetts Institute of Technology, Cambridge, Massachusetts.

                        All Rights Reserved

Permission to use, copy, modify, and distribute this software and its 
documentation for any purpose and without fee is hereby granted, 
provided that the above copyright notice appear in all copies and that
both that copyright notice and this permission notice appear in 
supporting documentation, and that the names of Digital or MIT not be
used in advertising or publicity pertaining to distribution of the
software without specific, written prior permission.  

DIGITAL DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING
ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL
DIGITAL BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR
ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
SOFTWARE.

 ******************************************************************/
/* $XConsortium: mipolygen.c,v 5.0 89/06/09 15:08:35 keith Exp $ */

extern void miloadAET(), micomputeWAET(), miFreeStorage();
extern int miCreateETandAET(), miInsertionSort();

/*
 *     Written by Brian Kelleher;  Oct. 1985
 *
 *     Routine to fill a polygon.  Two fill rules are
 *     supported: frWINDING and frEVENODD.
 *
 *     See fillpoly.h for a complete description of the algorithm.
 *
 *     Modified by Norman Kuring to fill a bin-number array instead of
 *     filling shapes on the server.
 */

int scan_convert(XPoint * ptsIn) {
    int count = 8;
    register EdgeTableEntry *pAET; /* the Active Edge Table   */
    register int y; /* the current scanline    */
    register int nPts = 0; /* number of pts in buffer */
    register EdgeTableEntry *pWETE; /* Winding Edge Table      */
    register ScanLineList *pSLL; /* Current ScanLineList    */
    register XPoint *ptsOut; /* ptr to output buffers   */
    int *wdth;
    XPoint FirstPoint[NUMPTSTOBUFFER]; /* the output buffers */
    int FirstWidth[NUMPTSTOBUFFER];
    EdgeTableEntry *pPrevAET; /* previous AET entry      */
    EdgeTable ET; /* Edge Table header node  */
    EdgeTableEntry AET; /* Active ET header node   */
    EdgeTableEntry *pETEs; /* Edge Table Entries buff */
    ScanLineListBlock SLLBlock; /* header for ScanLineList */
    int fixWAET = 0;

    if (!
            (pETEs =
            (EdgeTableEntry *) malloc(sizeof (EdgeTableEntry) * count))) {
        return (0);
    }
    ptsOut = FirstPoint;
    wdth = FirstWidth;
    if (!miCreateETandAET(count, ptsIn, &ET, &AET, pETEs, &SLLBlock)) {
        free(pETEs);
        return (0);
    }
    pSLL = ET.scanlines.next;

    /*
     *  for each scanline
     */
    for (y = ET.ymin; y < ET.ymax; y++) {
        /*
         *  Add a new edge to the active edge table when we
         *  get to the next edge.
         */
        if (pSLL && y == pSLL->scanline) {
            miloadAET(&AET, pSLL->edgelist);
            micomputeWAET(&AET);
            pSLL = pSLL->next;
        }
        pPrevAET = &AET;
        pAET = AET.next;
        pWETE = pAET;

        /*
         *  for each active edge
         */
        while (pAET) {
            /*
             *  if the next edge in the active edge table is
             *  also the next edge in the winding active edge
             *  table.
             */
            if (pWETE == pAET) {
                ptsOut->x = pAET->bres.minor;
                ptsOut++->y = y;
                *wdth++ = pAET->nextWETE->bres.minor - pAET->bres.minor;
                nPts++;

                /*
                 *  send out the buffer
                 */
                if (nPts == NUMPTSTOBUFFER) {
                    collect_bins(nPts, FirstPoint, FirstWidth);
                    ptsOut = FirstPoint;
                    wdth = FirstWidth;
                    nPts = 0;
                }

                pWETE = pWETE->nextWETE;
                while (pWETE != pAET)
                    EVALUATEEDGEWINDING(pAET, pPrevAET, y, fixWAET);
                pWETE = pWETE->nextWETE;
            }
            EVALUATEEDGEWINDING(pAET, pPrevAET, y, fixWAET);
        }

        /*
         *  reevaluate the Winding active edge table if we
         *  just had to resort it or if we just exited an edge.
         */
        if (miInsertionSort(&AET) || fixWAET) {
            micomputeWAET(&AET);
            fixWAET = 0;
        }
    }

    /*
     *     Get any spans that we missed by buffering
     */
    collect_bins(nPts, FirstPoint, FirstWidth);
    free(pETEs);
    miFreeStorage(SLLBlock.next);
    return (1);
}

int collect_bins(int number_of_initial_points,
        XPoint * initial_point, int *span_width) {

    int i;
    short x, y;
    int end;

    /* The variables, width, height, binsoverlain, and numoverlain, are
     *  static global varibles. */

    for (i = 0; i < number_of_initial_points; i++) {

        y = initial_point[i].y;
        if (y >= height || y < 0)
            continue;

        for (x = initial_point[i].x,
                end = initial_point[i].x + span_width[i]; x < end; x++) {
            if (x >= width || x < 0)
                continue;
            if (numoverlain >= binsperpixel) {
                binsperpixel += 1024;
                binsoverlain =
                        (int *) realloc(binsoverlain,
                        binsperpixel * sizeof (int));
                if (binsoverlain == NULL) {
                    fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
                    fprintf(stderr,
                            "Memory allocation error (numoverlain=%d binsperpixel=%d\n",
                            numoverlain, binsperpixel);
                    exit(EXIT_FAILURE);
                }
                return (0);
            }
            binsoverlain[numoverlain++] = y * width + x;
        }
    }
    return (1);
}
