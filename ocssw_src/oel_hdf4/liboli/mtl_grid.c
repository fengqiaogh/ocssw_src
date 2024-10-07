/*
 * Name:	mtl_grid.c
 *
 * Purpose:	Routines that create and use the metadata grid file for per-pixel
 *		view and sun angle generation.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mtl_grid.h"
#define	TINY	1.0e-12

int calc_band_grid(
  WRS2 parms,				/* I: WRS2 system parameters */
  SMETA smeta,				/* I: Scene metadata structure */
  double dellon,               	        /* I: Path longitude adjustment (in degrees) */
  VECTOR ecfsun,			/* I: ECF sun vector */
  int band_index,			/* I: Index for current band */
  MTL_GRID_BAND *bgrid )		/* O: Band grid structure */
{
  int nrow;				/* Row index */
  int ncol;				/* Column index */
  int sca;				/* SCA index */

  /* Allocate the grid */
  bgrid->band = smeta.band_smeta[band_index].band;
  bgrid->nsca = smeta.band_smeta[band_index].nsca;
   /* Sudipta new addition for OLI to extract SCA and Det numbers */
  bgrid->num_detectors = smeta.band_smeta[band_index].num_detectors;
  /* End Sudipta new addition for OLI */
  if ( (bgrid->sgrid = (MTL_GRID_SCA *)malloc( bgrid->nsca * sizeof( MTL_GRID_SCA ) )) == NULL )
  {
    printf("Unable to allocate grid for band %d\n", bgrid->band );
    return(ERROR);
  }

  /* Loop through the SCAs */
  for ( sca=0; sca<bgrid->nsca; sca++ )
  {
    bgrid->sgrid[sca].min_line = 1.0e9;
    bgrid->sgrid[sca].max_line = -1.0e9;
    bgrid->sgrid[sca].min_samp = 1.0e9;
    bgrid->sgrid[sca].max_samp = -1.0e9;
    /* Loop through the grid */
    for ( nrow=0; nrow<NUM_GRID_ROW; nrow++ )
      for ( ncol=0; ncol<NUM_GRID_COL; ncol++ )
      {
        bgrid->sgrid[sca].pgrid[nrow][ncol].ndet = GRID_COL_BASE + (double)ncol * GRID_COL_INC;
        bgrid->sgrid[sca].pgrid[nrow][ncol].frow = GRID_ROW_BASE + (double)nrow * GRID_ROW_INC;
        if ( calc_grid_point( parms, smeta, dellon, band_index, sca+1, ecfsun, &bgrid->sgrid[sca].pgrid[nrow][ncol] ) != SUCCESS )
        {
          printf("Error generating grid point:  SCA: %d, Row: %d, Col: %d\n", sca+1, nrow, ncol);
          free( bgrid->sgrid );
          return(ERROR);
        }
        if ( bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_line < bgrid->sgrid[sca].min_line )
             bgrid->sgrid[sca].min_line = bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_line;
        if ( bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_line > bgrid->sgrid[sca].max_line )
             bgrid->sgrid[sca].max_line = bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_line;
        if ( bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_samp < bgrid->sgrid[sca].min_samp )
             bgrid->sgrid[sca].min_samp = bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_samp;
        if ( bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_samp > bgrid->sgrid[sca].max_samp )
             bgrid->sgrid[sca].max_samp = bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_samp;
      }
    /* Calculate line/sample to ndet/frow mapping */
    if ( calc_grid_fit( &bgrid->sgrid[sca] ) != SUCCESS )
    {
        printf("Error calculating grid mapping coefficients for SCA %d\n", sca+1);
        free( bgrid->sgrid );
        return(ERROR);
    }
  }

  return(SUCCESS);
}


// need a forward declaration
int calc_grid_point_ls(
  WRS2 parms,				/* I: WRS2 system parameters */
  SMETA smeta,				/* I: Scene metadata structure */
  double dellon,                        /* I: Path longitude adjustment (in degrees) */
  int band_index,			/* I: Index for current band */
  int sca,				/* I: SCA number (1 to nsca) */
  VECTOR ecfsun,			/* I: ECF sun vector */
  MTL_GRID_PT *gpt );			/* I/O: Grid point with only ndet and frow on input */


/* variant for older LS */
int calc_band_grid_ls(
  WRS2 parms,				/* I: WRS2 system parameters */
  SMETA smeta,				/* I: Scene metadata structure */
  double dellon,               	        /* I: Path longitude adjustment (in degrees) */
  VECTOR ecfsun,			/* I: ECF sun vector */
  int band_index,			/* I: Index for current band */
  MTL_GRID_BAND_LS *bgrid )		/* O: Band grid structure */
{
  int nrow;				/* Row index */
  int ncol;				/* Column index */
  int sca = 0;				/* SCA index */

  /* Allocate the grid */
  bgrid->band = smeta.band_smeta_ls[band_index].band;
  bgrid->nsca = smeta.band_smeta_ls[band_index].nsca;
  if ( (bgrid->sgrid = (MTL_GRID_SCA_LS *)malloc( bgrid->nsca * sizeof( MTL_GRID_SCA_LS ) )) == NULL )
  {
    printf("Unable to allocate grid for band %d\n", bgrid->band );
    return(ERROR);
  }

  /* Loop through the SCAs */
  for ( sca=0; sca<bgrid->nsca; sca++ )
  {
    bgrid->sgrid[sca].min_line = 1.0e9;
    bgrid->sgrid[sca].max_line = -1.0e9;
    bgrid->sgrid[sca].min_samp = 1.0e9;
    bgrid->sgrid[sca].max_samp = -1.0e9;
    /* Loop through the grid */
    for ( nrow=0; nrow<NUM_GRID_ROW_LS; nrow++ )
      for ( ncol=0; ncol<NUM_GRID_COL; ncol++ )
      {
        bgrid->sgrid[sca].pgrid[nrow][ncol].npix = GRID_COL_BASE + (double)ncol * GRID_COL_INC;
        bgrid->sgrid[sca].pgrid[nrow][ncol].frow = GRID_ROW_BASE_LS + (double)nrow * GRID_ROW_INC;
        bgrid->sgrid[sca].pgrid[nrow][ncol].ndet = GRID_COL_BASE + (double)ncol * GRID_COL_INC;
        if ( calc_grid_point_ls( parms, smeta, dellon, band_index, sca+1, ecfsun, &bgrid->sgrid[sca].pgrid[nrow][ncol] ) != SUCCESS )
        {
          printf("Error generating grid point:  SCA: %d, Row: %d, Col: %d\n", sca+1, nrow, ncol);
          free( bgrid->sgrid );
          return(ERROR);
        }
        if ( bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_line < bgrid->sgrid[sca].min_line )
             bgrid->sgrid[sca].min_line = bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_line;
        if ( bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_line > bgrid->sgrid[sca].max_line )
             bgrid->sgrid[sca].max_line = bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_line;
        if ( bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_samp < bgrid->sgrid[sca].min_samp )
             bgrid->sgrid[sca].min_samp = bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_samp;
        if ( bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_samp > bgrid->sgrid[sca].max_samp )
             bgrid->sgrid[sca].max_samp = bgrid->sgrid[sca].pgrid[nrow][ncol].l1t_samp;
      }
    /* Calculate line/sample to ndet/frow mapping */
    if ( calc_grid_fit_ls( &bgrid->sgrid[sca] ) != SUCCESS )
    {
        printf("Error calculating grid mapping coefficients for SCA %d\n", sca+1);
        free( bgrid->sgrid );
        return(ERROR);
    }
  }

  return(SUCCESS);
}

int calc_grid_fit(
  MTL_GRID_SCA *sgrid )		/* I/O: SCA grid array to fit polynomial coefficients for */
{
  int nrow;				/* Row index */
  int ncol;				/* Column index */
  double gline;				/* Normalized grid line coordinate */
  double gsamp;				/* Normalized grid sample coordinate */
  double norm[4][4];			/* Normal equations matrix */
  double rhs1[4];			/* NDet constant vector */
  double rhs2[4];			/* FRow constant vector */

  /* Initialize the normal equations */
  norm[0][0] = 0.0;	norm[0][1] = 0.0;	norm[0][2] = 0.0;	norm[0][3] = 0.0;
  norm[1][0] = 0.0;	norm[1][1] = 0.0;	norm[1][2] = 0.0;	norm[1][3] = 0.0;
  norm[2][0] = 0.0;	norm[2][1] = 0.0;	norm[2][2] = 0.0;	norm[2][3] = 0.0;
  norm[3][0] = 0.0;	norm[3][1] = 0.0;	norm[3][2] = 0.0;	norm[3][3] = 0.0;
  rhs1[0] = 0.0;	rhs1[1] = 0.0;		rhs1[2] = 0.0;		rhs1[3] = 0.0;
  rhs2[0] = 0.0;	rhs2[1] = 0.0;		rhs2[2] = 0.0;		rhs2[3] = 0.0;

  /* Loop through the grid */
  for ( nrow=0; nrow<NUM_GRID_ROW; nrow++ )
    for ( ncol=0; ncol<NUM_GRID_COL; ncol++ )
    {
        gline = (sgrid->pgrid[nrow][ncol].l1t_line - sgrid->min_line) / (sgrid->max_line - sgrid->min_line);
        gsamp = (sgrid->pgrid[nrow][ncol].l1t_samp - sgrid->min_samp) / (sgrid->max_samp - sgrid->min_samp);
        norm[0][0] += 1.0;
        norm[0][1] += gline;
        norm[0][2] += gsamp;
        norm[0][3] += gline*gsamp;
        rhs1[0] += sgrid->pgrid[nrow][ncol].ndet;
        rhs2[0] += sgrid->pgrid[nrow][ncol].frow;
        norm[1][1] += gline*gline;
        norm[1][2] += gline*gsamp;
        norm[1][3] += gline*gline*gsamp;
        rhs1[1] += gline*sgrid->pgrid[nrow][ncol].ndet;
        rhs2[1] += gline*sgrid->pgrid[nrow][ncol].frow;
        norm[2][2] += gsamp*gsamp;
        norm[2][3] += gsamp*gline*gsamp;
        rhs1[2] += gsamp*sgrid->pgrid[nrow][ncol].ndet;
        rhs2[2] += gsamp*sgrid->pgrid[nrow][ncol].frow;
        norm[3][3] += gline*gsamp*gline*gsamp;
        rhs1[3] += gline*gsamp*sgrid->pgrid[nrow][ncol].ndet;
        rhs2[3] += gline*gsamp*sgrid->pgrid[nrow][ncol].frow;
    }

  norm[1][0] = norm[0][1];
  norm[2][0] = norm[0][2];
  norm[3][0] = norm[0][3];
  norm[2][1] = norm[1][2];
  norm[3][1] = norm[1][3];
  norm[3][2] = norm[2][3];

  /* Invert the normal equation matrix */
  if ( simple_inverse( 4, norm ) != SUCCESS )
  {
      printf("Error inverting normal equations.\n");
      return(ERROR);
  } 

  /* Calculate line/samp to ndet and frow coefficients */
  for ( nrow=0; nrow<4; nrow++ )
  {
      sgrid->ls_to_ndet[nrow] = norm[nrow][0]*rhs1[0] + norm[nrow][1]*rhs1[1] + norm[nrow][2]*rhs1[2] + norm[nrow][3]*rhs1[3];
      sgrid->ls_to_frow[nrow] = norm[nrow][0]*rhs2[0] + norm[nrow][1]*rhs2[1] + norm[nrow][2]*rhs2[2] + norm[nrow][3]*rhs2[3];
  }

  return(SUCCESS);
}

/* variant of above for older LS */
int calc_grid_fit_ls(
  MTL_GRID_SCA_LS *sgrid )		/* I/O: SCA grid array to fit polynomial coefficients for */
{
  int nrow;				/* Row index */
  int ncol;				/* Column index */
  double gline;				/* Normalized grid line coordinate */
  double gsamp;				/* Normalized grid sample coordinate */
  double norm[4][4];			/* Normal equations matrix */
  double rhs1[4];			/* NDet constant vector */
  double rhs2[4];			/* FRow constant vector */

  /* Initialize the normal equations */
  norm[0][0] = 0.0;	norm[0][1] = 0.0;	norm[0][2] = 0.0;	norm[0][3] = 0.0;
  norm[1][0] = 0.0;	norm[1][1] = 0.0;	norm[1][2] = 0.0;	norm[1][3] = 0.0;
  norm[2][0] = 0.0;	norm[2][1] = 0.0;	norm[2][2] = 0.0;	norm[2][3] = 0.0;
  norm[3][0] = 0.0;	norm[3][1] = 0.0;	norm[3][2] = 0.0;	norm[3][3] = 0.0;
  rhs1[0] = 0.0;	rhs1[1] = 0.0;		rhs1[2] = 0.0;		rhs1[3] = 0.0;
  rhs2[0] = 0.0;	rhs2[1] = 0.0;		rhs2[2] = 0.0;		rhs2[3] = 0.0;

  /* Loop through the grid */
  for ( nrow=0; nrow<NUM_GRID_ROW_LS; nrow++ )
    for ( ncol=0; ncol<NUM_GRID_COL; ncol++ )
    {
        gline = (sgrid->pgrid[nrow][ncol].l1t_line - sgrid->min_line) / (sgrid->max_line - sgrid->min_line);
        gsamp = (sgrid->pgrid[nrow][ncol].l1t_samp - sgrid->min_samp) / (sgrid->max_samp - sgrid->min_samp);
        norm[0][0] += 1.0;
        norm[0][1] += gline;
        norm[0][2] += gsamp;
        norm[0][3] += gline*gsamp;
        rhs1[0] += sgrid->pgrid[nrow][ncol].ndet;
        rhs2[0] += sgrid->pgrid[nrow][ncol].frow;
        norm[1][1] += gline*gline;
        norm[1][2] += gline*gsamp;
        norm[1][3] += gline*gline*gsamp;
        rhs1[1] += gline*sgrid->pgrid[nrow][ncol].ndet;
        rhs2[1] += gline*sgrid->pgrid[nrow][ncol].frow;
        norm[2][2] += gsamp*gsamp;
        norm[2][3] += gsamp*gline*gsamp;
        rhs1[2] += gsamp*sgrid->pgrid[nrow][ncol].ndet;
        rhs2[2] += gsamp*sgrid->pgrid[nrow][ncol].frow;
        norm[3][3] += gline*gsamp*gline*gsamp;
        rhs1[3] += gline*gsamp*sgrid->pgrid[nrow][ncol].ndet;
        rhs2[3] += gline*gsamp*sgrid->pgrid[nrow][ncol].frow;
    }

  norm[1][0] = norm[0][1];
  norm[2][0] = norm[0][2];
  norm[3][0] = norm[0][3];
  norm[2][1] = norm[1][2];
  norm[3][1] = norm[1][3];
  norm[3][2] = norm[2][3];

  /* Invert the normal equation matrix */
  if ( simple_inverse( 4, norm ) != SUCCESS )
  {
      printf("Error inverting normal equations.\n");
      return(ERROR);
  }

  /* Calculate line/samp to ndet and frow coefficients */
  for ( nrow=0; nrow<4; nrow++ )
  {
      sgrid->ls_to_ndet[nrow] = norm[nrow][0]*rhs1[0] + norm[nrow][1]*rhs1[1] + norm[nrow][2]*rhs1[2] + norm[nrow][3]*rhs1[3];
      sgrid->ls_to_frow[nrow] = norm[nrow][0]*rhs2[0] + norm[nrow][1]*rhs2[1] + norm[nrow][2]*rhs2[2] + norm[nrow][3]*rhs2[3];
  }

  return(SUCCESS);
}

int simple_inverse(
  int n,				/* I: Dimension of matrix to invert */
  double a[4][4] )			/* I/O: Matrix to invert */
{
  int i, j, k;				/* Loop indices */

  for ( k=0; k<n; k++ )
  {
    if ( fabs( a[k][k] ) < TINY ) return(ERROR);
    for ( j=0; j<n; j++ )
      if ( j != k ) a[k][j] /= a[k][k];
    a[k][k] = 1.0 / a[k][k];
    for ( i=0; i<n; i++ )
      if ( i != k )
      {
        for ( j=0; j<n; j++ )
          if ( j != k ) a[i][j] -= a[i][k]*a[k][j];
        a[i][k] *= -a[k][k];
      }
  }

  return(SUCCESS);
}
int calc_grid_point(
  WRS2 parms,				/* I: WRS2 system parameters */
  SMETA smeta,				/* I: Scene metadata structure */
  double dellon,                        /* I: Path longitude adjustment (in degrees) */
  int band_index,			/* I: Index for current band */
  int sca,				/* I: SCA number (1 to nsca) */
  VECTOR ecfsun,			/* I: ECF sun vector */
  MTL_GRID_PT *gpt )			/* I/O: Grid point with only ndet and frow on input */
{
  VECTOR pos;           /* Spacecraft ECEF position vector */
  VECTOR vel;           /* Spacecraft ECEF velocity vector */
  VECTOR slos;          /* Spacecraft line-of-sight vector */
  IAS_VECTOR gpos;	/* Ground point ECEF vector */
  VECTOR ecfview;	/* ECF view vector */
  VECTOR cursun;	/* ECF sun vector rotated based on frow */
  double gp_lat;        /* Ground point latitude */
  double gp_lon;        /* Ground point longitude */
  double proj_x;        /* Ground point projection X coordinate */
  double proj_y;        /* Ground point projection Y coordinate */
  double ecf2lsr[3][3];	/* ECEF to LSR rotation matrix */
  double sunrot;	/* Solar rotation angle due to fractional WRS row */
  double d2r;           /* Conversion from degrees to radians */

  /* Compute conversion constant */
  d2r = atan(1.0)/45.0;
  sunrot = 360.0*d2r*parms.cycle/(double)(parms.numpath*parms.numrow)*gpt->frow;

  /* Rotate the sun vector to the current time */
  cursun.x =  ecfsun.x*cos(sunrot) + ecfsun.y*sin(sunrot);
  cursun.y = -ecfsun.x*sin(sunrot) + ecfsun.y*cos(sunrot);
  cursun.z =  ecfsun.z;

  /* Find the nominal ephemeris state vector */
  pathrow_to_posvel( parms, smeta.wrs_path, dellon, (double)smeta.wrs_row+gpt->frow, &pos, &vel );

  /* Calculate the instrument LOS vector */
  if ( calc_los( smeta.band_smeta[band_index], sca, gpt->ndet, &slos ) != SUCCESS )
  {
      printf("Error generating LOS vector for grid point.\n");
      return(ERROR);
  }

  /* Project to ground latitude/longitude */
  if ( calc_yaw_steering_gp( parms, pos, vel, slos, smeta.roll_angle, &gp_lat, &gp_lon ) < 0 )
  {
      printf("Error projecting nominal scene center.\n");
      return(ERROR);
  }
  gp_lat *= d2r;
  gp_lon *= d2r;

  /* Convert lat/lon to ECEF */
  (void)smeta_geodetic_to_ecef( gp_lat, gp_lon, 0.0, &gpos );

  /* Compute view vector */
  slos.x = pos.x - gpos.x;
  slos.y = pos.y - gpos.y;
  slos.z = pos.z - gpos.z;
  unitvec( slos, &ecfview );

  /* Construct the ECEF to LSR rotation matrix */
  (void)smeta_geodetic_to_ecf2lsr( gp_lat, gp_lon, ecf2lsr );

  /* Convert vectors to LSR */
  rotatevec( ecf2lsr, ecfview, &gpt->V );
  rotatevec( ecf2lsr, cursun, &gpt->S );

  /* Apply scene map projection */
  if ( smeta_geodetic_to_proj( smeta.projection, gp_lat, gp_lon, &proj_x, &proj_y ) != SUCCESS )
  {
      printf("Error converting lat/lon to projection X/Y.\n");
      return(ERROR);
  }

  /* Reduce to L1T line/sample */
  if ( smeta_proj_to_l1t( &smeta, smeta.band_smeta[band_index].band, proj_x, proj_y, &gpt->l1t_line, &gpt->l1t_samp ) != SUCCESS )
  {
      printf("Errro converting X/Y to L1T line/sample.\n");
      return(ERROR);
  }

  return(SUCCESS);
}

/* Variant of above for older LS */
int calc_grid_point_ls(
  WRS2 parms,				/* I: WRS2 system parameters */
  SMETA smeta,				/* I: Scene metadata structure */
  double dellon,                        /* I: Path longitude adjustment (in degrees) */
  int band_index,			/* I: Index for current band */
  int sca,				/* I: SCA number (1 to nsca) */
  VECTOR ecfsun,			/* I: ECF sun vector */
  MTL_GRID_PT *gpt )			/* I/O: Grid point with only ndet and frow on input */
{
  VECTOR pos;           /* Spacecraft ECEF position vector */
  VECTOR vel;           /* Spacecraft ECEF velocity vector */
  VECTOR slos;          /* Spacecraft line-of-sight vector */
  IAS_VECTOR gpos;	/* Ground point ECEF vector */
  VECTOR ecfview;	/* ECF view vector */
  VECTOR cursun;	/* ECF sun vector rotated based on frow */
  double gp_lat;        /* Ground point latitude */
  double gp_lon;        /* Ground point longitude */
  double proj_x;        /* Ground point projection X coordinate */
  double proj_y;        /* Ground point projection Y coordinate */
  double ecf2lsr[3][3];	/* ECEF to LSR rotation matrix */
  double sunrot;	/* Solar rotation angle due to fractional WRS row */
  double d2r;           /* Conversion from degrees to radians */

  /* Compute conversion constant */
  d2r = atan(1.0)/45.0;
  sunrot = 360.0*d2r*parms.cycle/(double)(parms.numpath*parms.numrow)*gpt->frow;

  /* Rotate the sun vector to the current time */
  cursun.x =  ecfsun.x*cos(sunrot) + ecfsun.y*sin(sunrot);
  cursun.y = -ecfsun.x*sin(sunrot) + ecfsun.y*cos(sunrot);
  cursun.z =  ecfsun.z;

  /* Find the nominal ephemeris state vector */
  pathrow_to_posvel( parms, smeta.wrs_path, dellon, (double)smeta.wrs_row+gpt->frow, &pos, &vel );

  /* Calculate the instrument LOS vector */
  if ( calc_los_ls( smeta.band_smeta_ls[band_index], (double)smeta.wrs_row+gpt->frow, gpt->npix,
                 gpt->ndet, &slos ) != SUCCESS )
  {
      printf("Error generating LOS vector for grid point.\n");
      return(ERROR);
  }

  /* Project to ground latitude/longitude */
  if ( calc_yaw_steering_gp_ls( parms, pos, vel, slos, smeta.roll_angle, &gp_lat, &gp_lon ) < 0 )
  {
      printf("Error projecting nominal scene center.\n");
      return(ERROR);
  }
  gp_lat *= d2r;
  gp_lon *= d2r;

  /* Convert lat/lon to ECEF */
  (void)smeta_geodetic_to_ecef( gp_lat, gp_lon, 0.0, &gpos );

  /* Compute view vector */
  slos.x = pos.x - gpos.x;
  slos.y = pos.y - gpos.y;
  slos.z = pos.z - gpos.z;
  unitvec( slos, &ecfview );

  /* Construct the ECEF to LSR rotation matrix */
  (void)smeta_geodetic_to_ecf2lsr( gp_lat, gp_lon, ecf2lsr );

  /* Convert vectors to LSR */
  rotatevec( ecf2lsr, ecfview, &gpt->V );
  rotatevec( ecf2lsr, cursun, &gpt->S );

  /* Apply scene map projection */
  if ( smeta_geodetic_to_proj( smeta.projection, gp_lat, gp_lon, &proj_x, &proj_y ) != SUCCESS )
  {
      printf("Error converting lat/lon to projection X/Y.\n");
      return(ERROR);
  }

  /* Reduce to L1T line/sample */
  if ( smeta_proj_to_l1t( &smeta, smeta.band_smeta_ls[band_index].band, proj_x, proj_y, &gpt->l1t_line, &gpt->l1t_samp ) != SUCCESS )
  {
      printf("Errro converting X/Y to L1T line/sample.\n");
      return(ERROR);
  }

  return(SUCCESS);
}

int smeta_angles(
  double line,			/* I: L1T line coordinate */
  double samp,			/* I: L1T sample coordinate */
  MTL_GRID_BAND bgrid,		/* I: MTL grid structure for current band */
  /* Sudipta new addition for OLI to extract SCA and Det numbers */
  short *sca_num,		/* O: SCA containing this point */
  short *det_num,		/* O: Detector number for this point */
  /* End Sudipta new addition for OLI */
  double *sat_zn,		/* O: Viewing zenith angle scaled to 0.01 degrees */
  double *sat_az,		/* O: Viewing azimuth angle scaled to 0.01 degrees */
  double *sun_zn,		/* O: Solar zenith angle scaled to 0.01 degrees */
  double *sun_az	)		/* O: Solar azimuth angle scaled to 0.01 degrees */
{
  int sca;			/* SCA index */
  int nsca;			/* Number of SCAs that see this point */
  double sline;			/* Scaled line number */
  double ssamp;			/* Scaled sample number */
  double ndet;			/* Normalized detector coordinate */
  double frow;			/* Fractional row coordinate */
  double w0, w1, w2, w3;	/* Interpolation weights */
  double sat_zen;		/* Floating point view zenith angle */
  double sat_azm;		/* Floating point view azimuth angle */
  double sun_zen;		/* Floating point sun zenith angle */
  double sun_azm;		/* Floating point sun azimuth angle */
  double hdist;			/* Horizontal component of vector */
  int cell_row;			/* Index of grid cell row */
  int cell_col;			/* Index of grid cell column */
  VECTOR sat[2];		/* View vector */
  VECTOR sun[2];		/* Solar vector */
  /* Sudipta new addition for OLI to extract SCA and Det numbers */
  short curr_sca[2];		/* Current SCA number(s) */
  short curr_det[2];		/* Current detector number(s) */
  /* End Sudipta new addition for OLI */
  double r2d00;			/* Units conversion factor */

  /* Scaling factor */
  r2d00 = 45.0/atan(1.0);

  /* Initialize SCA count and loop through SCAs */
  nsca = 0;
  for ( sca=0; sca<bgrid.nsca; sca++ )
  {
//      if ( nsca > 1 ) continue;
      if ( line < bgrid.sgrid[sca].min_line ) continue;
      if ( line > bgrid.sgrid[sca].max_line ) continue;
      if ( samp < bgrid.sgrid[sca].min_samp ) continue;
      if ( samp > bgrid.sgrid[sca].max_samp ) continue;

      /* Calculate normalized detector */
      sline = (line - bgrid.sgrid[sca].min_line) / (bgrid.sgrid[sca].max_line - bgrid.sgrid[sca].min_line);
      ssamp = (samp - bgrid.sgrid[sca].min_samp) / (bgrid.sgrid[sca].max_samp - bgrid.sgrid[sca].min_samp);
      ndet = bgrid.sgrid[sca].ls_to_ndet[0]
           + bgrid.sgrid[sca].ls_to_ndet[1]*sline
           + bgrid.sgrid[sca].ls_to_ndet[2]*ssamp
           + bgrid.sgrid[sca].ls_to_ndet[3]*sline*ssamp;
      if ( ndet < -1.0 || ndet > 1.0 ) continue;

      /* Calculate fractional row */
      frow = bgrid.sgrid[sca].ls_to_frow[0]
           + bgrid.sgrid[sca].ls_to_frow[1]*sline
           + bgrid.sgrid[sca].ls_to_frow[2]*ssamp
           + bgrid.sgrid[sca].ls_to_frow[3]*sline*ssamp;
      if ( frow < -0.7 || frow > 0.7 ) continue;

      /* Sudipta new addition for OLI to extract SCA and Det numbers */
       /* Capture the current SCA and detector numbers */
      curr_sca[nsca] = sca + 1;
      curr_det[nsca] = (short)floor((ndet + 1.0)/2.0 * (double)bgrid.num_detectors + 0.5);
      /* End Sudipta new addition for OLI */

      /* Figure out which grid cell we are in */
      cell_col = (int)floor( (ndet - GRID_COL_BASE)/GRID_COL_INC );
      if ( cell_col > (NUM_GRID_COL-2) ) continue;
      cell_row = (int)floor( (frow - GRID_ROW_BASE)/GRID_ROW_INC );
      if ( cell_row > (NUM_GRID_ROW-2) ) continue;

      /* Calculate offsets within the grid cell */
      ndet = (ndet - GRID_COL_BASE)/GRID_COL_INC - (double)cell_col;
      frow = (frow - GRID_ROW_BASE)/GRID_ROW_INC - (double)cell_row;

      /* Interpolate vectors */
      w0 = (1.0-ndet)*(1.0-frow);
      w1 = (1.0-ndet)*frow;
      w2 = ndet*(1.0-frow);
      w3 = ndet*frow;
      sat[nsca].x = bgrid.sgrid[sca].pgrid[cell_row][cell_col].V.x    *w0
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col].V.x  *w1
                  + bgrid.sgrid[sca].pgrid[cell_row][cell_col+1].V.x  *w2
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col+1].V.x*w3;
      sat[nsca].y = bgrid.sgrid[sca].pgrid[cell_row][cell_col].V.y    *w0
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col].V.y  *w1
                  + bgrid.sgrid[sca].pgrid[cell_row][cell_col+1].V.y  *w2
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col+1].V.y*w3;
      sat[nsca].z = sqrt( 1.0 - sat[nsca].x*sat[nsca].x - sat[nsca].y*sat[nsca].y );
      sun[nsca].x = bgrid.sgrid[sca].pgrid[cell_row][cell_col].S.x    *w0
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col].S.x  *w1
                  + bgrid.sgrid[sca].pgrid[cell_row][cell_col+1].S.x  *w2
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col+1].S.x*w3;
      sun[nsca].y = bgrid.sgrid[sca].pgrid[cell_row][cell_col].S.y    *w0
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col].S.y  *w1
                  + bgrid.sgrid[sca].pgrid[cell_row][cell_col+1].S.y  *w2
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col+1].S.y*w3;
      sun[nsca].z = sqrt( 1.0 - sun[nsca].x*sun[nsca].x - sun[nsca].y*sun[nsca].y );
      nsca++;
  }

  /* Set angles to zero if point not in scene */
  if ( nsca < 1 )
  {
      *sat_zn = 0;
      *sat_az = 0;
      *sun_zn = 0;
      *sun_az = 0;
      /* Sudipta new addition for OLI to extract SCA and Det numbers */
      *sca_num = 0;
      *det_num = 0;
      /* End Sudipta new addition for OLI */
      return(SUCCESS);
  }

  /* Calculate the angles */
  sat_zen = 0.0;	sat_azm = 0.0;	sun_zen = 0.0;	sun_azm = 0.0;
  for ( sca=0; sca<nsca; sca++ )
  {
      /* Reduce vectors to angles */
      if ( sat[sca].z > 0.0 ) sat_zen += acos( sat[sca].z );
      hdist = sat[sca].x*sat[sca].x + sat[sca].y*sat[sca].y;
      if ( hdist > 0.0 ) sat_azm += atan2( sat[sca].x, sat[sca].y );
      if ( sun[sca].z > 0.0 ) sun_zen += acos( sun[sca].z );
      hdist = sun[sca].x*sun[sca].x + sun[sca].y*sun[sca].y;
      if ( hdist > 0.0 ) sun_azm += atan2( sun[sca].x, sun[sca].y );
  }
  sat_zen /= (double)nsca;
  sat_azm /= (double)nsca;
  sun_zen /= (double)nsca;
  sun_azm /= (double)nsca;

  *sat_zn = r2d00*sat_zen;
  *sat_az = r2d00*sat_azm;
  *sun_zn = r2d00*sun_zen;
  *sun_az = r2d00*sun_azm;
  /* Sudipta new addition for OLI to extract SCA and Det numbers */
  *sca_num = curr_sca[0];
  *det_num = curr_det[0];

    if ( nsca > 1 )
    {
        if ( MIN( curr_det[0], bgrid.num_detectors-curr_det[0] ) < MIN( curr_det[1], bgrid.num_detectors-curr_det[1] ) )
        {
            *sca_num = curr_sca[1];
            *det_num = curr_det[1];
        }
    }
    /* End Sudipta new addition for OLI */

  return(SUCCESS);
}

/* varinat of smeta_angles for older LS */
int smeta_angles_ls(
  double line,			/* I: L1T line coordinate */
  double samp,			/* I: L1T sample coordinate */
  MTL_GRID_BAND_LS bgrid,		/* I: MTL grid structure for current band */
  short *sat_zn,		/* O: Viewing zenith angle scaled to 0.01 degrees */
  short *sat_az,		/* O: Viewing azimuth angle scaled to 0.01 degrees */
  short *sun_zn,		/* O: Solar zenith angle scaled to 0.01 degrees */
  short *sun_az	)		/* O: Solar azimuth angle scaled to 0.01 degrees */
{
  int sca;			/* SCA index */
  int nsca;			/* Number of SCAs that see this point */
  double sline;			/* Scaled line number */
  double ssamp;			/* Scaled sample number */
  double ndet;			/* Normalized detector coordinate */
  double frow;			/* Fractional row coordinate */
  double w0, w1, w2, w3;	/* Interpolation weights */
  double sat_zen;		/* Floating point view zenith angle */
  double sat_azm;		/* Floating point view azimuth angle */
  double sun_zen;		/* Floating point sun zenith angle */
  double sun_azm;		/* Floating point sun azimuth angle */
  double hdist;			/* Horizontal component of vector */
  int cell_row;			/* Index of grid cell row */
  int cell_col;			/* Index of grid cell column */
  VECTOR sat[2];		/* View vector */
  VECTOR sun[2];		/* Solar vector */
  double r2d00;			/* Units conversion factor */

  /* Scaling factor */
  r2d00 = 4500.0/atan(1.0);

  /* Initialize SCA count and loop through SCAs */
  nsca = 0;
  for ( sca=0; sca<bgrid.nsca; sca++ )
  {
      if ( nsca > 1 ) continue;
      if ( line < bgrid.sgrid[sca].min_line ) continue;
      if ( line > bgrid.sgrid[sca].max_line ) continue;
      if ( samp < bgrid.sgrid[sca].min_samp ) continue;
      if ( samp > bgrid.sgrid[sca].max_samp ) continue;

      /* Calculate normalized detector */
      sline = (line - bgrid.sgrid[sca].min_line) / (bgrid.sgrid[sca].max_line - bgrid.sgrid[sca].min_line);
      ssamp = (samp - bgrid.sgrid[sca].min_samp) / (bgrid.sgrid[sca].max_samp - bgrid.sgrid[sca].min_samp);
      ndet = bgrid.sgrid[sca].ls_to_ndet[0]
           + bgrid.sgrid[sca].ls_to_ndet[1]*sline
           + bgrid.sgrid[sca].ls_to_ndet[2]*ssamp
           + bgrid.sgrid[sca].ls_to_ndet[3]*sline*ssamp;
      if ( ndet < -1.0 || ndet > 1.0 ) continue;

      /* Calculate fractional row */
      frow = bgrid.sgrid[sca].ls_to_frow[0]
           + bgrid.sgrid[sca].ls_to_frow[1]*sline
           + bgrid.sgrid[sca].ls_to_frow[2]*ssamp
           + bgrid.sgrid[sca].ls_to_frow[3]*sline*ssamp;
      if ( frow < -0.6 || frow > 0.6 ) continue;

      /* Figure out which grid cell we are in */
      cell_col = (int)floor( (ndet - GRID_COL_BASE)/GRID_COL_INC );
      if ( cell_col > (NUM_GRID_COL-2) ) continue;
      cell_row = (int)floor( (frow - GRID_ROW_BASE_LS)/GRID_ROW_INC );
      if ( cell_row > (NUM_GRID_ROW_LS-2) ) continue;

      /* Calculate offsets within the grid cell */
      ndet = (ndet - GRID_COL_BASE)/GRID_COL_INC - (double)cell_col;
      frow = (frow - GRID_ROW_BASE_LS)/GRID_ROW_INC - (double)cell_row;

      /* Interpolate vectors */
      w0 = (1.0-ndet)*(1.0-frow);
      w1 = (1.0-ndet)*frow;
      w2 = ndet*(1.0-frow);
      w3 = ndet*frow;
      sat[nsca].x = bgrid.sgrid[sca].pgrid[cell_row][cell_col].V.x    *w0
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col].V.x  *w1
                  + bgrid.sgrid[sca].pgrid[cell_row][cell_col+1].V.x  *w2
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col+1].V.x*w3;
      sat[nsca].y = bgrid.sgrid[sca].pgrid[cell_row][cell_col].V.y    *w0
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col].V.y  *w1
                  + bgrid.sgrid[sca].pgrid[cell_row][cell_col+1].V.y  *w2
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col+1].V.y*w3;
      sat[nsca].z = sqrt( 1.0 - sat[nsca].x*sat[nsca].x - sat[nsca].y*sat[nsca].y );
      sun[nsca].x = bgrid.sgrid[sca].pgrid[cell_row][cell_col].S.x    *w0
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col].S.x  *w1
                  + bgrid.sgrid[sca].pgrid[cell_row][cell_col+1].S.x  *w2
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col+1].S.x*w3;
      sun[nsca].y = bgrid.sgrid[sca].pgrid[cell_row][cell_col].S.y    *w0
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col].S.y  *w1
                  + bgrid.sgrid[sca].pgrid[cell_row][cell_col+1].S.y  *w2
                  + bgrid.sgrid[sca].pgrid[cell_row+1][cell_col+1].S.y*w3;
      sun[nsca].z = sqrt( 1.0 - sun[nsca].x*sun[nsca].x - sun[nsca].y*sun[nsca].y );
      nsca++;
  }

  /* Set angles to zero if point not in scene */
  if ( nsca < 1 )
  {
      *sat_zn = 0;
      *sat_az = 0;
      *sun_zn = 0;
      *sun_az = 0;
      return(SUCCESS);
  }

  /* Calculate the angles */
  sat_zen = 0.0;	sat_azm = 0.0;	sun_zen = 0.0;	sun_azm = 0.0;
  for ( sca=0; sca<nsca; sca++ )
  {
      /* Reduce vectors to angles */
      if ( sat[sca].z > 0.0 ) sat_zen += acos( sat[sca].z );
      hdist = sat[sca].x*sat[sca].x + sat[sca].y*sat[sca].y;
      if ( hdist > 0.0 ) sat_azm += atan2( sat[sca].x, sat[sca].y );
      if ( sun[sca].z > 0.0 ) sun_zen += acos( sun[sca].z );
      hdist = sun[sca].x*sun[sca].x + sun[sca].y*sun[sca].y;
      if ( hdist > 0.0 ) sun_azm += atan2( sun[sca].x, sun[sca].y );
  }
  sat_zen /= (double)nsca;
  sat_azm /= (double)nsca;
  sun_zen /= (double)nsca;
  sun_azm /= (double)nsca;

  *sat_zn = (short)floor( r2d00*sat_zen + 0.5 );
  *sat_az = (short)floor( r2d00*sat_azm + 0.5 );
  *sun_zn = (short)floor( r2d00*sun_zen + 0.5 );
  *sun_az = (short)floor( r2d00*sun_azm + 0.5 );

  return(SUCCESS);
}
