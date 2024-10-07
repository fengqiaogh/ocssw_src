#include <stdlib.h>
#include <stdio.h>
#include "get_ctht.h"

/*
  file ams_oe_inversion.c, A. Sayer's inversion routine, currently limited to 
    y_true of size 8 (8 rho bands) and xa of 3 (3 prods)
*/
int32_t ams_oe_inversion( double *rhot, gsl_matrix *gm_sy, double *xa, 
  gsl_matrix *gm_sa, double *lut, int32_t *lut_dims, int32_t *min_cost_loc, 
  oe_info_str *oe_info )

 /*
    ams_oe_inversion - Generic optimal estimation (OE) inversion routine

    Returns int32_t  0 if all good

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     double *          rhot             I      sample reflectance to match
     gsl_matrix *      gm_sy            I      Measurement/forward model 
                                               covariance matrix
     double *          xa               I      A priori state vector (product 
                                               we want to get)
     gsl_matrix *      gm_sa            I      A priori state covariance matrix
     double *          lut              I      table of rho for this pixel 
                                               view and sfc p level
     int32_t *         lut_dims         I      all dims of the LUT
     int32_t *         min_cost_loc     I      first guess at x.
     oe_info_str *     oe_info          O      a structure of the inversion 
                                               results and uncertainty

  W. Robinson, SAIC, 18 Aug 2023, converted from A. Sayer's IDL code

There are other controls: nitr_max=nitr_max, djconv=djconv,mqstep=mqstep,
mq0=mq0, verbose=verbose - perhaps some/all will be added

*/
  {
  int32_t nx, ny, ix, iy, idim, nitr, nitr_max, conv_flag, mqstep;
  int32_t tmp_attempt_ctr, new_k_init;
  double j_cost_m, j_cost_a, j_cost, j_old, djconv, mq, mq0;
  double yres[8], xres[3], xpdx[3];
  double j_old_t1, j_old_t2;  //  will need alloc if variable;

  static int32_t f_init = 0;
  static double *y_init, *k_init;
  static gsl_vector *gv_yres, *gv_xres, *gv_inv_sy_yres, *gv_inv_sa_xres;
  static gsl_vector *gv_xpdx, *gv_x_fg, *gv_xa, *gv_r_minus_xa, *gv_jp_t2;
  static gsl_vector *gv_y_init_rho, *gv_jp_t1, *gv_jp_t1_f2;
  static gsl_vector *gv_jp, *gv_dx, *gv_nx, *gv_ny;
  static gsl_matrix *gm_sa_inv, *gm_sy_inv_k_init, *gm_ident;
  static gsl_matrix *gm_jpp_t1, *gm_k_init, *gm_jpp, *gm_scl_ident;
  static gsl_matrix *gm_jpp_ident, *gm_jpp_aux, *gm_jpp_ident_inv;
  static gsl_matrix *gm_sx, *gm_sx_t2_f1, *gm_sx_t2, *gm_sx_uinv;
  static gsl_matrix *gm_gain, *gm_av_kernel, *gm_sy_inv, *gm_g_pt1;
  static gsl_matrix *gm_g_pt2, *gm_g_pt3, *gm_g_pt4;

  gsl_matrix *mat_xfr;  // for doing some transfers

  nx = gm_sa->size1;
  ny = lut_dims[3];

  if( f_init == 0 )
    {
    y_init = (double *)malloc( ny * sizeof(double) );
    k_init = (double  *)malloc( ny * nx * sizeof(double) );
    f_init = 1;

    gm_jpp_t1 = gsl_matrix_alloc( nx, nx );
    gv_x_fg = gsl_vector_alloc(nx);
    gm_ident = gsl_matrix_alloc( nx, nx );
    gv_yres = gsl_vector_alloc(ny);
    gv_xres = gsl_vector_alloc(nx);
    gv_xpdx = gsl_vector_alloc( nx );
    gv_inv_sa_xres = gsl_vector_alloc(nx);
    gv_inv_sy_yres = gsl_vector_alloc(ny);
    gm_sy_inv_k_init = gsl_matrix_alloc( ny, nx );
    gm_k_init = gsl_matrix_alloc( nx, ny );
    gv_jp_t1_f2 = gsl_vector_alloc( ny );
    gm_jpp = gsl_matrix_alloc( nx,nx );
    gv_xa = gsl_vector_alloc( nx );
    gv_r_minus_xa = gsl_vector_alloc( nx );
    gv_jp_t2 = gsl_vector_alloc( nx );
    gv_y_init_rho = gsl_vector_alloc( ny );
    gv_jp_t1 = gsl_vector_alloc( nx );
    gv_jp = gsl_vector_alloc( nx );
    gm_scl_ident = gsl_matrix_alloc( nx, nx );
    gm_jpp_ident = gsl_matrix_alloc( nx, nx );
    gm_jpp_aux = gsl_matrix_alloc( nx, nx );
    gv_dx = gsl_vector_alloc( nx );
    gv_nx = gsl_vector_alloc( nx );
    gv_ny = gsl_vector_alloc( ny );
    gm_sx = gsl_matrix_alloc( nx, nx );
    gm_sx_t2_f1 = gsl_matrix_alloc( ny, nx );
    gm_sx_t2 = gsl_matrix_alloc( nx, nx );
    gm_sx_uinv = gsl_matrix_alloc( nx, nx );
    gm_av_kernel = gsl_matrix_alloc( nx, nx );
    gm_g_pt1 = gsl_matrix_alloc( nx, ny );
    gm_g_pt2 = gsl_matrix_alloc( nx, nx );
    gm_g_pt3 = gsl_matrix_alloc( nx, nx );
    gm_g_pt4 = gsl_matrix_alloc( nx, nx );
    gm_gain = gsl_matrix_alloc( nx, ny );
    //  some output struct init
    oe_info->x_prd_state = (float *) malloc( nx * sizeof(float) );
    }
  //  we'll initialize the optional controls
  j_cost = 1.e20;  //  Cost: initialise to high value
  j_old = j_cost;  //  Previous attempt's cost
  nitr = 0;  //  Number of iterations attempted
  conv_flag = 0;  //  Flag for convergence (0=not converged, 1=converged)
  //  Default parameter values, if keywords not set
  nitr_max = 50;  //  Maximum number of iterations
  djconv = 0.01;  //  Threshold change in cost
  mqstep = 5;  //  Marquardt step parameter
  mq0 = 0.00001;  //  Marquardt parameter initial value, Don't start this 
                 //  too high or can oscillate.
  int32_t snax[3];
  new_k_init = 1;

  for( idim = 0; idim < 3; idim++ )
    snax[idim] = lut_dims[idim];
 /*
  *  Compute initial guess, y, and cost - this can also be set to xa
  */
  for( ix = 0; ix < nx; ix++ )
    gsl_vector_set( gv_x_fg, ix, min_cost_loc[ix] );
  //  set identity matrix
  gsl_matrix_set_identity( gm_ident );
 /*
  *  Get initial values of y from x
  */
  my_lut_funct_3d_oe( gv_x_fg, lut, snax, ny, y_init, k_init );
 /*
  * Calculate cost
  */
  //  set xres, yres and vector versions
  for( idim = 0; idim < ny; idim++ )
    yres[idim] = y_init[idim] - rhot[idim];

  for( idim = 0; idim < ny; idim++ )
    gsl_vector_set( gv_yres, idim, yres[idim] );

 for( idim = 0; idim < nx; idim++ )
    xres[idim] = gv_x_fg->data[idim] - xa[idim];

  for( idim = 0; idim < nx; idim++ )
    gsl_vector_set( gv_xres, idim, xres[idim] );

  // Do line:  j_old=yres#la_invert(Sy)#yres + xres#la_invert(Sa)#xres
  gm_sy_inv = invert_a_matrix( gm_sy );
  gm_sa_inv = invert_a_matrix( gm_sa );

  // term 1: yres#la_invert(Sy)#yres
  gsl_blas_dgemv( CblasNoTrans, 1., gm_sy_inv, gv_yres, 0., gv_inv_sy_yres );
  gsl_blas_ddot( gv_yres, gv_inv_sy_yres, &j_old_t1 );

  // term 2: xres#la_invert(Sa)#xres
  gsl_blas_dgemv( CblasNoTrans, 1., gm_sa_inv, gv_xres, 0., gv_inv_sa_xres );
  gsl_blas_ddot( gv_xres, gv_inv_sa_xres, &j_old_t2 );

  j_old = j_old_t1 + j_old_t2;
  //  Finished making j_old

 /*
  *  For first attempt set mq0 to trace of J''
  *  jpp=ym.k#la_invert(Sy)#transpose(ym.k) + la_invert(Sa)
  *    ym.k -> k_init  Sy -> gm_sy  Sa -> gm_sa
  */
  //  1st, la_invert(Sy)#transpose(ym.k) sy_inv_k_init
  for( iy = 0; iy < ny; iy++ )
    for( ix = 0; ix < nx; ix++ )
      gsl_matrix_set( gm_k_init, ix, iy, *( k_init + ix + nx * iy ) );

  gsl_blas_dgemm( CblasNoTrans, CblasTrans, 1., gm_sy_inv, gm_k_init, 0., gm_sy_inv_k_init );
  //  finally, worked.  on to the other parts and adding it together
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1., gm_k_init, gm_sy_inv_k_init, 0, gm_jpp_t1 );

  //  the 2nd term is just an inverse of Sa so just have to get jpp
  gsl_matrix_memcpy( gm_jpp, gm_sa_inv );
  gsl_matrix_add( gm_jpp, gm_jpp_t1 );

  //  and to get mq0,= trace(jpp/nx) no gsl for that
  double tsum = 0.;
  for( ix = 0; ix < nx; ix++ )
    tsum += gsl_matrix_get( gm_jpp, ix, ix );
  mq0 *= tsum / nx;

 /*
  *  Next, Do the iteration
  *  looks like only 1 new matrix, vector operation set along with 
  *  previous matrix results
  *  So, get the whole loop to end sketched out
  */
  while( ( conv_flag == 0 ) && ( nitr <= nitr_max ) )
    {
    nitr++;
    //  Reset Marquardt parameter to default value
    //   variable 'mq' represents value used in this step
    mq = mq0;

    //  Set max number of attempts per iteration in case of
    //  e.g. numerical instability problems
    tmp_attempt_ctr = 0;

    while( ( j_cost >= j_old ) && ( tmp_attempt_ctr <= 100 ) )
      {
      tmp_attempt_ctr++;
      //  This should not result in an infinite loop
      //  unless gradients are calculated incorrectly
      //  or a parameter is unstable,as it will
      //  tend towards making very small steps, which
      //  must lead to a decreased cost, unless you
      //  are already at a minimum in which case
      //  the cost function should have converged.

      // OK, sep 29 23 update to call my_lut_funct_3d_oe here
      // Call forward model to get prediction and gradients.
      // Necessary as otherwise things would get out of sync if an iteration
      // attempt is unsuccessful
      my_lut_funct_3d_oe( gv_x_fg, lut, snax, ny, y_init, k_init );
      //  don't forget to make gsl version of k_init
      for( iy = 0; iy < ny; iy++ )
        for( ix = 0; ix < nx; ix++ )
          gsl_matrix_set( gm_k_init, ix, iy, *( k_init + ix + nx * iy ) );

      //  Increase Marquardt parameter
      mq *= mqstep;
      //  Get step for next iteration attempt:
      //  Need gradient of cost function, J', denoted jp here
      //  *** NEW jp=ym.k#la_invert(Sy)#(ym.y-y_true) + la_invert(Sa)#(x-xa)
      //             [3, 8]  [8,8]      [8]             [3,3]         3]
      //             siz 3                                   siz 3
      //          siz 3
      //  jp= k_init # gm_sy_inv # ( y_init - rhot ) + gm_sa_inv # ( rhot - xa )

      //  need to make gv_xa from xa here 
      for( ix = 0; ix < nx; ix++ )
        gsl_vector_set( gv_xa, ix, xa[ix] );

      //  create rhot - xa
      gsl_vector_memcpy( gv_r_minus_xa, gv_x_fg );
      gsl_vector_sub( gv_r_minus_xa, gv_xa );
      //  gv_r_minus_xa is (x-xa) above

      // make 2nd term, la_invert(Sa)#(x-xa)
      // or: gm_sa_inv # ( rhot - xa ) (la_invert(Sa)#(x-xa))
      gsl_blas_dgemv( CblasNoTrans, 1., gm_sa_inv, gv_r_minus_xa, 0., 
        gv_jp_t2 );

      //  1st term, k_init # gm_sy_inv # ( y_init - rhot )
      //  do  y_init - rhot part
      for( iy = 0; iy < ny; iy++ )
        gsl_vector_set( gv_y_init_rho, iy, y_init[iy] - rhot[iy] );

      //  intermediate - la_invert(Sy)#(ym.y-y_true)
      gsl_blas_dgemv( CblasNoTrans, 1., gm_sy_inv, gv_y_init_rho, 0., gv_jp_t1_f2);

      //  ok, ym.k#la_invert(Sy)#(ym.y-y_true)
      gsl_blas_dgemv( CblasNoTrans, 1., gm_k_init, gv_jp_t1_f2, 0., gv_jp_t1 );

      //  finally, jp_t1 + jp_t2
      gsl_vector_memcpy( gv_jp, gv_jp_t1 );
      gsl_vector_add( gv_jp, gv_jp_t2 );
      //  OK, jp is made
      
      //  and J'', jpp
      //  Use the jpp computed before  for 1st time through
      // printf( "%s, %d, I: Note that we already have jpp from before the while\n", __FILE__, __LINE__ );
      if( new_k_init == 1 )
        {
        //  WDR NOTE that this is a repeat of code above - should make sub
       /*
        *  For first attempt set mq0 to trace of J''
        *  jpp=ym.k#la_invert(Sy)#transpose(ym.k) + la_invert(Sa)
        *   ym.k -> k_init  Sy -> gm_sy  Sa -> gm_sa
        */
        //  1st, la_invert(Sy)#transpose(ym.k)
        for( iy = 0; iy < ny; iy++ )
          for( ix = 0; ix < nx; ix++ )
            gsl_matrix_set( gm_k_init, ix, iy, *( k_init + ix + nx * iy ) );

        gsl_blas_dgemm( CblasNoTrans, CblasTrans, 1., gm_sy_inv, gm_k_init, 
          0., gm_sy_inv_k_init );
        //  finally, worked.  on to the other parts and adding it together
        gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1., gm_k_init, 
          gm_sy_inv_k_init, 0, gm_jpp_t1 );
      
        //  the 2nd term is just an inverse of Sa so just have to get jpp
        gsl_matrix_memcpy( gm_jpp, gm_sa_inv );
        gsl_matrix_add( gm_jpp, gm_jpp_t1 );
        }
      new_k_init = 1;

      //  This is the step size dx
      //  *** NEW TOO dx = -1*la_invert(jpp + mq*ident_matrix)#jp

      //  mq*ident_matrix: gm_scl_ident
      gsl_matrix_memcpy( gm_scl_ident, gm_ident );
      gsl_matrix_scale( gm_scl_ident, mq );

      //  jpp + mq*ident_matrix: gm_jpp_ident
      gsl_matrix_memcpy( gm_jpp_ident, gm_scl_ident );
      gsl_matrix_memcpy( gm_jpp_aux, gm_jpp );
      gsl_matrix_add( gm_jpp_ident, gm_jpp_aux ); 

      // la_invert(jpp + mq*ident_matrix):  gm_jpp_ident_inv
      gm_jpp_ident_inv = invert_a_matrix( gm_jpp_ident );

      //  -1*la_invert(jpp + mq*ident_matrix)#jp:  gv_dx
      gsl_blas_dgemv( CblasNoTrans, -1., gm_jpp_ident_inv, gv_jp, 0., gv_dx );
      gsl_matrix_free( gm_jpp_ident_inv );
      //  DX MADE
      
      //  curtail the search if the dx has a NaN
      for( ix = 0; ix < nx; ix++ )
        {
        if( isnan( gv_dx->data[ix] ) )
          {
          gsl_vector_set_zero( gv_dx );
          tmp_attempt_ctr=1000;
          break;
          }
        }
      //  Create new x+dx array, make sure remains within bounds
      for( ix = 0; ix < nx; ix++ )
        {
        xpdx[ix] = gv_x_fg->data[ix] + gv_dx->data[ix];  //  go gsl?
        xpdx[ix] = ( xpdx[ix] <= 1.e-5 ) ? 1.e-5 : xpdx[ix];
        xpdx[ix] = ( xpdx[ix] > lut_dims[ix] - 1 ) ? lut_dims[ix] - 1 : xpdx[ix];
        }
      //  Try next step: see if cost lower
      for( ix = 0; ix < nx; ix++ )
        gsl_vector_set( gv_xpdx, ix, xpdx[ix] );
      my_lut_funct_3d_oe( gv_xpdx, lut, snax, ny, y_init, k_init );
      //  get gm_k_init
      for( iy = 0; iy < ny; iy++ )
        for( ix = 0; ix < nx; ix++ )
          gsl_matrix_set( gm_k_init, ix, iy, *( k_init + ix + nx * iy ) );

      //  re-compute the x,yres and vector versions
      for( idim = 0; idim < ny; idim++ )
        {
        yres[idim] = y_init[idim] - rhot[idim];
        gsl_vector_set( gv_yres, idim, yres[idim] );
        }
      for( idim = 0; idim < nx; idim++ )
        {
        xres[idim] = gv_x_fg->data[idim] - xa[idim];
        gsl_vector_set( gv_xres, idim, xres[idim] );
        }

      for( idim = 0; idim < ny; idim++ )
        gsl_vector_set( gv_yres, idim, yres[idim] );
      //  re-compute the Measurement contribution to cost
      //  and A priori contribution to cost,  currently repeat of dbl_int1, 2
      //  derivation higher up
      // term 1
      gsl_blas_dgemv( CblasNoTrans, 1., gm_sy_inv, gv_yres, 0., gv_ny );
      gsl_blas_ddot( gv_yres, gv_ny, &j_cost_m );

      // term 2
      gsl_blas_dgemv( CblasNoTrans, 1., gm_sa_inv, gv_xres, 0., gv_nx );
      gsl_blas_ddot( gv_xres, gv_nx, &j_cost_a );

      //  j_cost_m =  Measurement contribution to cost or
      //     jm=yres#la_invert(Sy)#yres
      //  j_cost_a = A priori contribution to cost or
      //     ja=xres#la_invert(Sa)#xres
      j_cost = j_cost_m + j_cost_a;

      if( isnan( j_cost ) )
        {
        printf( "%s, %d, E: Meaningless cost found\n", __FILE__, __LINE__ );
        exit(1);
        }
//  temp to print out the values during iteration
/*
printf( "\nEnd step loop, nitr: %d, tmp_attempt_ctr: %d\n", nitr, 
  tmp_attempt_ctr );
printf( "dx values:\n" );
gsl_vector_fprintf( stdout, gv_dx, "%f" );
printf( "mq value: %f\n", mq );
printf( "jp values:\n" );
gsl_vector_fprintf( stdout, gv_jp, "%f" );
printf( "jpp:\n" );
print_mat_contents( gm_jpp, nx, nx );
printf( "old, new cost at this point: %f, %f\n", j_old, j_cost );
printf( "END contents\n" );
*/
      }  //  End loop trying to find next step

      //  Update state vector
      gsl_vector_memcpy( gv_x_fg, gv_xpdx );

      //  If cost change small, or measurement cost is very low, we have converged
      if( ( fabs( j_cost - j_old ) < ( djconv * ny ) ) ||
          ( j_cost_m < djconv * ny ) ) conv_flag = 1;
//  temp WDR show items used to decide convergence
// printf( "CONV INFO, j_cost: %f, j_old: %f, djconv: %f, ny: %d, j_cost_m: %f\n", j_cost, j_old, djconv, ny, j_cost_m );

      // If cost decreased but we've not converged, carry on
      // using smaller Marquardt parameter for next step
      j_old = j_cost;
      mq0 = mq / mqstep;
// printf( "\nEnd iteration loop, nitr: %d, tmp_attempt_ctr: %d\n", nitr, 
//  tmp_attempt_ctr );
    }  //  End iteration loop

//;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
//; Calculate information about the solution ;
//;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

//  Calculate uncertainty on retrieved state
//  Assume model parameter error is included in Sy
//  Put Sy as double here as otherwise from time to time float issues can
//  make the expression uninvertable
// Sx=la_invert(la_invert(Sa) + ym.k#la_invert(double(Sy))#transpose(ym.k))
//  make gm_sx_t2_f1 = la_invert(double(Sy))#transpose(ym.k)
gsl_blas_dgemm( CblasNoTrans, CblasTrans, 1., gm_sy_inv, gm_k_init, 0., gm_sx_t2_f1 );
//   make gm_sx_t2 = ym.k#la_invert(double(Sy))#transpose(ym.k)
gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1., gm_k_init, gm_sx_t2_f1, 0., gm_sx_t2 );

gsl_matrix_memcpy( gm_sx_uinv, gm_sa_inv );
gsl_matrix_add( gm_sx_uinv, gm_sx_t2 );  //  may want better descrip of term1
   //  it is prev computed ym.k#la_invert(Sy)#transpose(ym.k))
//  la_invert( gm_sx_uninv )
mat_xfr = invert_a_matrix( gm_sx_uinv );
gsl_matrix_memcpy( gm_sx, mat_xfr );
gsl_matrix_free( mat_xfr );


// Calculate gain matrix: gm_gain
// G = la_invert((la_invert(Sa) + $
//   ym.k # la_invert(double(Sy)) # transpose(ym.k))) # ym.k#la_invert(Sy)
//  this is 
//  gm_g = la_invert( ( gm_sa_inv + gm_k_init # gm_sy_inv # gm_k_init^T ) )
//     # gm_k_init # gm_sy_inv
//
// part 1: gm_k_init # gm_sy_inv -> gm_g_pt1 // appears 2x in above
gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1., gm_k_init, gm_sy_inv, 0.,
  gm_g_pt1 );

// part 2: gm_g_pt1 # gm_k_init^T -> gm_g_pt2
gsl_blas_dgemm( CblasNoTrans, CblasTrans, 1., gm_g_pt1, gm_k_init, 0., gm_g_pt2 );

// part 3: gm_sa_inv + gm_g_pt3
gsl_matrix_memcpy( gm_g_pt3, gm_sa_inv );
gsl_matrix_add( gm_g_pt3, gm_g_pt2 );

// part 4: invert( gm_g_pt3 ) -> gm_g_pt4
mat_xfr = invert_a_matrix( gm_g_pt3 );
gsl_matrix_memcpy( gm_g_pt4, mat_xfr );
gsl_matrix_free( mat_xfr );

//  so we did gm_k_init # gm_sy_inv as gm_g_pt1
//  just get the gain:  gm_g_pt4 # gm_g_pt1
gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1., gm_g_pt4, gm_g_pt1, 0., gm_gain );

//  END Gain matrix comp

//  Now use this to calculate averaging kernel, A=GKt  [nx,nx]
//  gm_gain # gm_k_init^T = gm_av_kernel[nx,nx]
gsl_blas_dgemm( CblasNoTrans, CblasTrans, 1., gm_gain, gm_k_init, 
  0., gm_av_kernel );

// Return output
//  fill a structure with this
for( ix = 0; ix < nx; ix++ )  //  the 'x'
  oe_info->x_prd_state[ix] = gsl_vector_get( gv_x_fg, ix );
oe_info->gm_sx = gm_sx;
oe_info->conv_flag = conv_flag;
oe_info->cost = j_cost;
oe_info->acost = j_cost_a;
oe_info->gm_ak = gm_av_kernel;
oe_info->gm_gain = gm_gain;
oe_info->nitr = nitr;

//  And end?
gsl_matrix_free( gm_sa_inv );
gsl_matrix_free( gm_sy_inv );
return 0;
    
// printf( "%s, %d, I after hairy matrix mult\n", __FILE__, __LINE__ );
  //  and done
  return 0;
  }

  //  Well, we have to move on to my_lut_funct_3d_oe (not sure it is only '3'

int32_t my_lut_funct_3d_oe( gsl_vector *gv_x_fg, double *lut, int32_t *snax, 
  int32_t ny, double *y_init, double *k_init )
/*
    my_lut_funct_3d_oe - Get initial values of y from x

    Returns int32_t, 0 if all good

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     gsl_vector *      gv_x_fg          I      LUT indicies
     double *          lut              I      Specific LUT for this pixel
     int32_t *         snax             I      size of LUT less last, band dim
     int32_t           ny               I      # of bands
     double *          y_init           O      initial y value
     double *          k_init           O      nx by ny matrix

  W. Robinson, SAIC, 18 Aug 2023, converted from A. Sayer's IDL code

*/
  {
  int32_t ibnd, icor, nx, iprd;
  float ar_wt[8], ar_wtu[8], orig_wt[3];
  double yl, yu, dx, dlut;
  double xl[3], xu[3];
  int64_t ar_off[8], ar_offu[8], tmp_off;

  // printf( "%s, %d, I: entering my_lut_funct_3d_oe\n", __FILE__, __LINE__ );
  //  Interpolate to get Y - no need to re-derive points and weights
  iint_3d( snax, gv_x_fg->data, ar_off, ar_wt, orig_wt );
  //  apply for each band
  for( ibnd = 0; ibnd < ny; ibnd++ )
    {
    y_init[ibnd] = 0.;
    for( icor = 0; icor < 8; icor++ )
      {
      tmp_off = snax[0] * snax[1] * snax[2];
      y_init[ibnd] += *( lut + tmp_off * ibnd + ar_off[icor] ) * ar_wt[icor];
      }
    }
  //  Estimate K by finite difference
  dlut = .1;
  nx = 3;
  
  //  no need for zeroing the k_init, x_sizes are [snax,ny]
  //  get lower, upper 1st guess value
  for( iprd = 0; iprd < nx; iprd++ )
    {
    xl[iprd] = gv_x_fg->data[iprd];  //  lower step
    xu[iprd] = gv_x_fg->data[iprd];  //  upper step
    }
  //  I'm not sure of the alg below as xl and xu are used only after 
  //  changing  them in the ibnd location - don't know - WDR
  for( iprd = 0; iprd < nx; iprd++ )
    {
    // Calculate indices to compare for this point in state space
    // Make sure we don't go out of LUT bounds
    xl[iprd] = gv_x_fg->data[iprd] - dlut;
    if( xl[iprd] < 0. ) xl[iprd] = 0.;
    xu[iprd] = gv_x_fg->data[iprd] + dlut;
    if( xu[iprd] > snax[iprd]-1 ) xu[iprd] = snax[iprd]-1;

    dx = xu[iprd] - xl[iprd];
    //  Get weighting function for this parameter for each band
    //  get the offsets and weights for this perturbation
    iint_3d( snax, xl, ar_off, ar_wt, orig_wt );
    iint_3d( snax, xu, ar_offu, ar_wtu, orig_wt );
    //  and apply them to get yl, yu
    for( ibnd = 0; ibnd < ny; ibnd++ )
      {
      yl = 0.;
      yu = 0.;
      for( icor = 0; icor < 8; icor++ )
        {
        tmp_off = snax[0] * snax[1] * snax[2];
        yl += *( lut + tmp_off * ibnd + ar_off[icor] ) * ar_wt[icor];
        yu += *( lut + tmp_off * ibnd + ar_offu[icor] ) * ar_wtu[icor];
        }
      *( k_init + iprd + 3 * ibnd ) = ( yu - yl ) / dx;
      }
    }
  //  and the k_init and y_init are made
  // printf( "%s, %d, I: end of y_init, k_init set\n", __FILE__, __LINE__ );
  return 0;
  }
