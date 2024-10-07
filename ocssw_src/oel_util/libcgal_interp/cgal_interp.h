#ifndef CGAL_INTERP_H
#define CGAL_INTERP_H

int cgal_nnc( int nxy, float *x, float *y, int nxyq, float *xq, float *yq);
int cgal_interp2( int nxy, float *x, float *y, float *v, int nxyq, float *vq, int nnc_tag);
int cgal_release_tag( int nnc_tag);

#endif
