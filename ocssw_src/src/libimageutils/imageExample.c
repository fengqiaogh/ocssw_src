#include <stdio.h>
#include <stdlib.h>

#include <imageutils.h>
#include <png.h>

void write_png_file(char* filename, int width, int height, uint8_t* image_data,
        int num_colors, uint8_t* palette) {

    FILE* outfp = fopen(filename, "w");
    if (!outfp) {
        fprintf(stderr, "-E- Unable to open file %s.\n", filename);
        exit(EXIT_FAILURE);
    }

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
            NULL, NULL, NULL);
    if (!png_ptr) {
        fprintf(stderr, "-E- Unable to create PNG write structure.\n");
        exit(EXIT_FAILURE);
    }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
        fprintf(stderr, "-E- Unable to create PNG info structure.\n");
        exit(EXIT_FAILURE);
    }
    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
        fprintf(stderr, "-E- Unable to call PNG setjmp().\n");
        exit(EXIT_FAILURE);
    }
    png_init_io(png_ptr, outfp);

    // write the header
    png_set_IHDR(png_ptr, info_ptr, width, height,
            8, PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    // set the color palette
    png_set_PLTE(png_ptr, info_ptr, (png_const_colorp) palette, num_colors);
    png_write_info(png_ptr, info_ptr);

    // write lines
    int y;
    uint8_t* ptr = image_data;
    for (y = 0; y < height; y++) {
        png_write_row(png_ptr, (png_bytep) ptr);
        ptr += width;
    }

    // end writing
    png_write_end(png_ptr, info_ptr);
    png_destroy_write_struct(&png_ptr, &info_ptr);

    fclose(outfp);
}

int main(int c, char *v[]) {
    if (c < 3) {
        fprintf(stderr, "usage: %s ppm_file n_colors\n", v[0]);
        return 0;
    }

    int numColors = atoi(v[2]) ? : 16; /* GCC extension */

    img_rgb_t* im = img_read_ppm(v[1]);
    if (!im) {
        printf("could not read PPM file %s\n", v[1]);
        exit(EXIT_FAILURE);
    }

    // allocate output image and color palette
    uint8_t* palette = (uint8_t*) malloc(numColors * 3);
    uint8_t* out_image = (uint8_t*) malloc(im->w * im->h);

    img_color_palette_quantization(im, numColors, palette, out_image);

    write_png_file("out.png", im->w, im->h, out_image, numColors, palette);

    //    color_quant(im, numColors, 0);
    //    write_ppm(im, "out.ppm");


    free(im);

    return 0;
}