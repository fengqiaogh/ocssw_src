/*
 * This was adapted from the code found here.
 * 
 * https://rosettacode.org/wiki/Color_quantization/C
 */


#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include <imageutils.h>

#define ON_INHEAP 1



typedef struct oct_node_t oct_node_t, *oct_node;

struct oct_node_t {
    int64_t r, g, b; /* sum of all child node colors */
    int count, heap_idx;
    uint8_t n_kids, kid_idx, flags, depth;
    oct_node kids[8], parent;
};

typedef struct {
    int alloc, n;
    oct_node* buf;
} node_heap;


static oct_node pool = 0;

img_rgb_t* img_new(int w, int h) {
    img_rgb_t* img = malloc(sizeof (img_rgb_t) + h * w * 3);
    img->w = w;
    img->h = h;
    img->pix = (uint8_t *) (img + 1);
    return img;
}

void img_free(img_rgb_t* img) {
    // this works since the image structure and pixels memory
    // were allocated in one big block
    free(img);
}

int img_write_ppm(img_rgb_t* img, char *filename) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) return 0;
    fprintf(fp, "P6\n%d %d\n255\n", img->w, img->h);
    fwrite(img->pix, 1, 3 * img->w * img->h, fp);
    fclose(fp);
    return 1;
}

int read_num(FILE *f) {
    int n;
    while (!fscanf(f, "%d ", &n)) {
        if ((n = fgetc(f)) == '#') {
            while ((n = fgetc(f)) != '\n')
                if (n == EOF) return 0;
        } else
            return 0;
    }
    return n;
}

img_rgb_t* img_read_ppm(char *filename) {
    FILE *fp = fopen(filename, "rb");
    int w, h, maxval;
    img_rgb_t* img = 0;
    if (!fp) return 0;

    if (fgetc(fp) != 'P' || fgetc(fp) != '6' || !isspace(fgetc(fp)))
        goto bail;

    w = read_num(fp);
    h = read_num(fp);
    maxval = read_num(fp);
    if (!w || !h || !maxval) goto bail;

    img = img_new(w, h);
    fread(img->pix, 1, 3 * w * h, fp);
bail:
    if (fp) fclose(fp);
    return img;
}

/* cmp function that decides the ordering in the heap.  This is how we determine
   which octree node to fold next, the heart of the algorithm. */
int cmp_node(oct_node a, oct_node b) {
    if (a->n_kids < b->n_kids) return -1;
    if (a->n_kids > b->n_kids) return 1;

    int ac = a->count >> a->depth;
    int bc = b->count >> b->depth;
    return ac < bc ? -1 : ac > bc;
}

void down_heap(node_heap *h, oct_node p) {
    int n = p->heap_idx, m;
    while (1) {
        m = n * 2;
        if (m >= h->n) break;
        if (m + 1 < h->n && cmp_node(h->buf[m], h->buf[m + 1]) > 0) m++;

        if (cmp_node(p, h->buf[m]) <= 0) break;

        h->buf[n] = h->buf[m];
        h->buf[n]->heap_idx = n;
        n = m;
    }
    h->buf[n] = p;
    p->heap_idx = n;
}

void up_heap(node_heap *h, oct_node p) {
    int n = p->heap_idx;
    oct_node prev;

    while (n > 1) {
        prev = h->buf[n / 2];
        if (cmp_node(p, prev) >= 0) break;

        h->buf[n] = prev;
        prev->heap_idx = n;
        n /= 2;
    }
    h->buf[n] = p;
    p->heap_idx = n;
}

void heap_add(node_heap *h, oct_node p) {
    if ((p->flags & ON_INHEAP)) {
        down_heap(h, p);
        up_heap(h, p);
        return;
    }

    p->flags |= ON_INHEAP;
    if (!h->n) h->n = 1;
    if (h->n >= h->alloc) {
        while (h->n >= h->alloc) h->alloc += 1024;
        h->buf = realloc(h->buf, sizeof (oct_node) * h->alloc);
    }

    p->heap_idx = h->n;
    h->buf[h->n++] = p;
    up_heap(h, p);
}

oct_node pop_heap(node_heap *h) {
    if (h->n <= 1) return 0;

    oct_node ret = h->buf[1];
    h->buf[1] = h->buf[--h->n];

    h->buf[h->n] = 0;

    h->buf[1]->heap_idx = 1;
    down_heap(h, h->buf[1]);

    return ret;
}

oct_node node_new(uint8_t idx, uint8_t depth, oct_node p) {
    static int len = 0;
    if (len <= 1) {
        oct_node p = calloc(sizeof (oct_node_t), 2048);
        p->parent = pool;
        pool = p;
        len = 2047;
    }

    oct_node x = pool + len--;
    x->kid_idx = idx;
    x->depth = depth;
    x->parent = p;
    if (p) p->n_kids++;
    return x;
}

void node_free() {
    oct_node p;
    while (pool) {
        p = pool->parent;
        free(pool);
        pool = p;
    }
}

/* adding a color triple to octree */
oct_node node_insert(oct_node root, uint8_t *pix) {
    uint8_t i, bit, depth = 0;

    for (bit = 1 << 7; ++depth < 8; bit >>= 1) {
        i = !!(pix[1] & bit) * 4 + !!(pix[0] & bit) * 2 + !!(pix[2] & bit);
        if (!root->kids[i])
            root->kids[i] = node_new(i, depth, root);

        root = root->kids[i];
    }

    root->r += pix[0];
    root->g += pix[1];
    root->b += pix[2];
    root->count++;
    return root;
}

/* remove a node in octree and add its count and colors to parent node. */
oct_node node_fold(oct_node p) {
    if (p->n_kids) abort();
    oct_node q = p->parent;
    q->count += p->count;

    q->r += p->r;
    q->g += p->g;
    q->b += p->b;
    q->n_kids--;
    q->kids[p->kid_idx] = 0;
    return q;
}

/* traverse the octree just like construction, but this time we replace the pixel
   color with color stored in the tree node */
void color_replace(oct_node root, uint8_t *pix) {
    uint8_t i, bit;

    for (bit = 1 << 7; bit; bit >>= 1) {
        i = !!(pix[1] & bit) * 4 + !!(pix[0] & bit) * 2 + !!(pix[2] & bit);
        if (!root->kids[i]) break;
        root = root->kids[i];
    }

    pix[0] = root->r;
    pix[1] = root->g;
    pix[2] = root->b;
}

uint8_t find_color_index(oct_node root, uint8_t* pix) {
    uint8_t i, bit;

    for (bit = 1 << 7; bit; bit >>= 1) {
        i = !!(pix[1] & bit) * 4 + !!(pix[0] & bit) * 2 + !!(pix[2] & bit);
        if (!root->kids[i]) break;
        root = root->kids[i];
    }

    return root->heap_idx - 1;
}

/* Building an octree and keep leaf nodes in a bin heap.  Afterwards remove first node
   in heap and fold it into its parent node (which may now be added to heap), until heap
   contains required number of colors. */
void img_color_quant(img_rgb_t* im, int n_colors, int dither) {
    int i;
    uint8_t *pix = im->pix;
    node_heap heap = {0, 0, 0};

    oct_node root = node_new(0, 0, 0), got;
    for (i = 0; i < im->w * im->h; i++, pix += 3)
        heap_add(&heap, node_insert(root, pix));

    while (heap.n > n_colors + 1)
        heap_add(&heap, node_fold(pop_heap(&heap)));

    double c;
    for (i = 1; i < heap.n; i++) {
        got = heap.buf[i];
        c = got->count;
        got->r = got->r / c + .5;
        got->g = got->g / c + .5;
        got->b = got->b / c + .5;
    }

    for (i = 0, pix = im->pix; i < im->w * im->h; i++, pix += 3)
        color_replace(root, pix);

    node_free();
    free(heap.buf);
}

/**
 * Take the input 24 bit full color RGB image and convert it into a color
 * mapped image and color palette which have already been allocated by the caller.
 * 
 * @param in_image 24 bit RGB input image
 * @param num_colors number of colors in the color palette
 * @param palette color palette (allocated by caller) 3 bytes (r,g,b) for each entry
 * @param out_image output image data (allocated by caller)
 */
void img_color_palette_quantization(img_rgb_t* in_image, int num_colors,
        uint8_t* palette, uint8_t* out_image) {
    int i;
    uint8_t *pix = in_image->pix;
    node_heap heap = {0, 0, 0};

    // load up the oct tree and the heap
    oct_node root = node_new(0, 0, 0), got;
    for (i = 0; i < in_image->w * in_image->h; i++, pix += 3)
        heap_add(&heap, node_insert(root, pix));

    // reduce the colors to the requested number
    while (heap.n > num_colors + 1)
        heap_add(&heap, node_fold(pop_heap(&heap)));

    // write out the color palette
    double c;
    pix = palette;
    for (i = 1; i < heap.n; i++) {
        got = heap.buf[i];
        c = got->count;
        *(pix++) = got->r / c + .5;
        *(pix++) = got->g / c + .5;
        *(pix++) = got->b / c + .5;
    }
    for (; i <= num_colors; i++) { // zero out the unused palette entries
        *(pix++) = 0;
        *(pix++) = 0;
        *(pix++) = 0;
    }

    // write out the image data
    for (i = 0, pix = in_image->pix; i < in_image->w * in_image->h; i++, pix += 3) {
        out_image[i] = find_color_index(root, pix);
    }

    node_free();
    free(heap.buf);
}

