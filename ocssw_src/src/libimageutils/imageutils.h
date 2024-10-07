#ifndef IMAGE_UTILS_H
#define IMAGE_UTILS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**\file 
 * imageutils is a library of functions that work on images.
 */

typedef struct {
    int w, h; // width and height of image
    uint8_t *pix; // image data, first pixel is top left, next pixel is 
    // to the right. each pixel is 3 bytes (red, green, blue)
} img_rgb_t;


img_rgb_t* img_new(int w, int h);
void img_free(img_rgb_t* img);

int img_write_ppm(img_rgb_t* img, char *filename);
img_rgb_t* img_read_ppm(char *filename);

void img_color_quant(img_rgb_t* img, int n_colors, int dither);

void img_color_palette_quantization(img_rgb_t* in_image, int num_colors,
        uint8_t* palette, uint8_t* out_image);



#ifdef __cplusplus
}
#endif

#endif
