/*
 *	Name:	mtl_geometry.c
 *
 *	Purpose:
 *	Source file containing the functions used to perform geometric operations
 * 	on the data extracted from the MTL file.
 */

#include "mtl_geometry.h"

/*
 *	Name:	calc_los
 *
 *	Purpose:
 *	Source file containing the function used to compute OLI and TIRS line-of-sight
 *	vectors and transform them to the spacecraft coordinate system.
 */

int calc_los(
        SMETA_BAND band_smeta, /* I: Metadata band structure */
        int sca, /* I: SCA number */
        double ndet, /* I: Normalized detector coordinate */
        VECTOR *los) /* O: Line of sight unit vector */ {
    VECTOR ilos; /* Instrument LOS */
    VECTOR slos; /* Spacecraft LOS */

    /* Check SCA number */
    if (sca < 1 || sca > band_smeta.nsca) {
        printf("Invalid SCA %d requested.\n", sca);
        return (ERROR);
    }

    /* Calculate the instrument LOS vector */
    ilos.x = band_smeta.legendre[sca - 1][0][0]
            + band_smeta.legendre[sca - 1][0][1] * ndet
            + band_smeta.legendre[sca - 1][0][2] * (1.5 * ndet * ndet - 0.5)
            + band_smeta.legendre[sca - 1][0][3] * ndet * (2.5 * ndet * ndet - 1.5);
    ilos.y = band_smeta.legendre[sca - 1][1][0]
            + band_smeta.legendre[sca - 1][1][1] * ndet
            + band_smeta.legendre[sca - 1][1][2] * (1.5 * ndet * ndet - 0.5)
            + band_smeta.legendre[sca - 1][1][3] * ndet * (2.5 * ndet * ndet - 1.5);
    ilos.z = 1.0;
    unitvec(ilos, &slos);

    /* Rotate to spacecraft coordinates */
    rotatevec(band_smeta.align, slos, los);

    return (SUCCESS);
}

/* Different varant for older LS sensors */
int calc_los_ls(
        SMETA_BAND_LS band_smeta, /* I: Metadata band structure */
        double frow, /* I: Fractional row */
        double npix, /* I: Sample coordinate */
        double ndet, /* I: Normalized detector coordinate */
        VECTOR *los) /* O: Line of sight unit vector */ {
    VECTOR ilos; /* Instrument LOS */
    VECTOR slos; /* Spacecraft LOS */
    double along_angle; /* Along scan angle */
    double across_angle; /* Across scan angle */
    //double scan_number;			/* Current scan number */

    /* Calculate the instrument LOS vector */
    //scan_number = band_smeta.nscans * (1.0 + ndet)/2.0;

    across_angle = 0.0;
    along_angle = -band_smeta.fascan + (1.0 + npix) / 2.0 * (band_smeta.fascan + band_smeta.rascan) - band_smeta.aoffset;

    ilos.x = sin(across_angle);
    ilos.y = -cos(across_angle) * sin(along_angle);
    ilos.z = cos(across_angle) * cos(along_angle);
    unitvec(ilos, &slos);

    /* Rotate to spacecraft coordinates */
    rotatevec(band_smeta.align, slos, los);

    return (SUCCESS);
}

/******************************************************************************
NAME:           calc_yaw_steering_gp

PURPOSE:
Calculates the ground point intersection latitude and longitude for a specified
input instrument LOS vector, spacecraft ephemeris, and roll angle. Uses the
yaw steering pointing logic to do the calculation.

RETURN VALUE:
None.

NOTES:

ALGORITHM REFERENCES:
 ******************************************************************************
                        Property of the U.S. Government
                            USGS EROS Data Center
 ******************************************************************************/

int calc_yaw_steering_gp(
        WRS2 parms, /* I: WRS2 system parameters */
        VECTOR pos, /* I: Position vector */
        VECTOR vel, /* I: Velocity vector */
        VECTOR i_los, /* I: Instrument line-of-sight unit vector */
        double roll, /* I: Off-nadir roll angle (in degrees) */
        double *gp_lat, /* O: Ground point geodetic latitude (in degrees) */
        double *gp_lon) /* O: Ground point geodetic longitude (in degrees) */ {
    double a[3][3];
    VECTOR omega_los, /* LOS angular rate vector */
            u_los, /* LOS vector */
            o_x_u, /* angular rate cross LOS */
            rp, /* Earth point vector */
            vp, /* Earth point velocity */
            er, /* Earth rotation vector */
            er_x_rp, /* Earth rotation cross position */
            vgt, /* Target velocity */
            u_x_vgt, /* LOS cross target velocity */
            x, /* YSF X axis */
            y; /* YSF Y axis */
    double e_a, /* Earth semimajor axis */
            e_b, /* Earth semiminor axis */
            e_rate, /* Earth rotation rate */
            d2r, /* Degrees to radians conversion */
            scale, /* Orbit to Earth scale factor */
            rsc, /* Orbital radius */
            rho, /* Orbital altitude */
            rho_dot, /* Radial velocity */
            angmo, /* Scalar angular momentum */
            sang, /* Rescaled (for ellipsoid) roll angle */
            cosr, /* Cosine of roll angle */
            sinr; /* Sine of roll angle */

    /* Initialize Earth constants */
    e_a = parms.semimajor;
    e_b = parms.semiminor;
    e_rate = parms.inertialvel;
    d2r = atan(1.0) / 45.0;

    /* Compute roll trig functions */
    cosr = cos(roll * d2r);
    sinr = sin(roll * d2r);

    /* Compute orbital radius and angular momentum */
    rsc = vecnorm(pos);
    crossprod(pos, vel, &omega_los);
    angmo = vecnorm(omega_los);

    /* Construct the LOS */
    u_los.x = -cosr * pos.x / rsc + sinr * omega_los.x / angmo;
    u_los.y = -cosr * pos.y / rsc + sinr * omega_los.y / angmo;
    u_los.z = -cosr * pos.z / rsc + sinr * omega_los.z / angmo;

    /* Compute LOS angular rate */
    omega_los.x = cosr * angmo / (rsc * rsc)*(cosr * omega_los.x / angmo + sinr * pos.x / rsc);
    omega_los.y = cosr * angmo / (rsc * rsc)*(cosr * omega_los.y / angmo + sinr * pos.y / rsc);
    omega_los.z = cosr * angmo / (rsc * rsc)*(cosr * omega_los.z / angmo + sinr * pos.z / rsc);

    /* Compute orbit to Earth scale factor */
    scale = pos.x / e_a * pos.x / e_a + pos.y / e_a * pos.y / e_a + pos.z / e_b * pos.z / e_b;
    sang = u_los.x / e_a * u_los.x / e_a + u_los.y / e_a * u_los.y / e_a + u_los.z / e_b * u_los.z / e_b;
    sang = (-pos.x / e_a * u_los.x / e_a - pos.y / e_a * u_los.y / e_a - pos.z / e_b * u_los.z / e_b) / sqrt(sang * scale);
    sang = MIN(MAX(-1.0, sang), 1.0);
    sang = acos(sang);
    cosr = cos(sang);
    sinr = sin(sang);
    scale = 1.0 / scale - sinr*sinr;
    if (scale > 0.0) {
        scale = sqrt(scale);
    } else {
        printf("Roll angle too large, LOS misses Earth\n");
        return (-1);
    }
    rho = rsc * (cosr - scale);
    rho_dot = dotprod(pos, vel) / rsc * (cosr + sinr * sinr / scale);

    /* LOS rate cross LOS */
    crossprod(omega_los, u_los, &o_x_u);

    /* Construct ground target point */
    rp.x = pos.x + rho * u_los.x;
    rp.y = pos.y + rho * u_los.y;
    rp.z = pos.z + rho * u_los.z;
    vp.x = vel.x + rho_dot * u_los.x + rho * o_x_u.x;
    vp.y = vel.y + rho_dot * u_los.y + rho * o_x_u.y;
    vp.z = vel.z + rho_dot * u_los.z + rho * o_x_u.z;

    /* Construct the Earth rotation vector */
    er.x = 0.0;
    er.y = 0.0;
    er.z = e_rate;

    /* Compute er cross rp and ground target velocity */
    crossprod(er, rp, &er_x_rp);
    vgt.x = vp.x - er_x_rp.x;
    vgt.y = vp.y - er_x_rp.y;
    vgt.z = vp.z - er_x_rp.z;
    crossprod(u_los, vgt, &u_x_vgt);
    unitvec(u_x_vgt, &y);
    crossprod(y, u_los, &x);

    /* Construct ECI2YSF orientation matrix */
    a[0][0] = x.x;
    a[0][1] = x.y;
    a[0][2] = x.z;
    a[1][0] = y.x;
    a[1][1] = y.y;
    a[1][2] = y.z;
    a[2][0] = u_los.x;
    a[2][1] = u_los.y;
    a[2][2] = u_los.z;

    /* Rotate input LOS vector using the orientation matrix */
    u_los.x = a[0][0] * i_los.x + a[1][0] * i_los.y + a[2][0] * i_los.z;
    u_los.y = a[0][1] * i_los.x + a[1][1] * i_los.y + a[2][1] * i_los.z;
    u_los.z = a[0][2] * i_los.x + a[1][2] * i_los.y + a[2][2] * i_los.z;

    /* Compute orbit to Earth scale factor */
    scale = pos.x / e_a * pos.x / e_a + pos.y / e_a * pos.y / e_a + pos.z / e_b * pos.z / e_b;
    sang = u_los.x / e_a * u_los.x / e_a + u_los.y / e_a * u_los.y / e_a + u_los.z / e_b * u_los.z / e_b;
    sang = (-pos.x / e_a * u_los.x / e_a - pos.y / e_a * u_los.y / e_a - pos.z / e_b * u_los.z / e_b) / sqrt(sang * scale);
    sang = MIN(MAX(-1.0, sang), 1.0);
    sang = acos(sang);
    cosr = cos(sang);
    sinr = sin(sang);
    scale = 1.0 / scale - sinr*sinr;
    if (scale > 0.0) {
        scale = sqrt(scale);
    } else {
        printf("Roll angle too large, LOS misses Earth");
        return (-1);
    }
    rho = rsc * (cosr - scale);

    /* Construct ground target point */
    rp.x = pos.x + rho * u_los.x;
    rp.y = pos.y + rho * u_los.y;
    rp.z = pos.z + rho * u_los.z;

    /* Convert ground point to lat/lon */
    *gp_lat = atan(rp.z / sqrt(rp.x * rp.x + rp.y * rp.y) * e_a / e_b * e_a / e_b) / d2r;
    *gp_lon = atan2(rp.y, rp.x) / d2r;

    return (0);
}

/******************************************************************************
NAME:           calc_yaw_steering_gp_ls

PURPOSE:
Calculates the ground point intersection latitude and longitude for a specified
input instrument LOS vector, spacecraft ephemeris, and roll angle. Uses the
yaw steering pointing logic to do the calculation.

RETURN VALUE:
None.

NOTES:

ALGORITHM REFERENCES:
 ******************************************************************************
                        Property of the U.S. Government
                            USGS EROS Data Center
 ******************************************************************************/

int calc_yaw_steering_gp_ls(
        WRS2 parms, /* I: WRS2 system parameters */
        VECTOR pos, /* I: Position vector */
        VECTOR vel, /* I: Velocity vector */
        VECTOR i_los, /* I: Instrument line-of-sight unit vector */
        double roll, /* I: Off-nadir roll angle (in degrees) */
        double *gp_lat, /* O: Ground point geodetic latitude (in degrees) */
        double *gp_lon) /* O: Ground point geodetic longitude (in degrees) */ {
    double a[3][3];
    VECTOR omega_los, /* LOS angular rate vector */
            u_los, /* LOS vector */
            o_x_u, /* angular rate cross LOS */
            rp, /* Earth point vector */
            vp, /* Earth point velocity */
            er, /* Earth rotation vector */
            er_x_rp, /* Earth rotation cross position */
            vgt, /* Target velocity */
            u_x_vgt, /* LOS cross target velocity */
            x, /* YSF X axis */
            y; /* YSF Y axis */
    double e_a, /* Earth semimajor axis */
            e_b, /* Earth semiminor axis */
            e_rate, /* Earth rotation rate */
            d2r, /* Degrees to radians conversion */
            scale, /* Orbit to Earth scale factor */
            rsc, /* Orbital radius */
            rho, /* Orbital altitude */
            rho_dot, /* Radial velocity */
            angmo, /* Scalar angular momentum */
            sang, /* Rescaled (for ellipsoid) roll angle */
            cosr, /* Cosine of roll angle */
            sinr; /* Sine of roll angle */

    /* Initialize Earth constants */
    e_a = parms.semimajor;
    e_b = parms.semiminor;
    e_rate = parms.inertialvel;
    d2r = atan(1.0) / 45.0;

    /* Compute roll trig functions */
    cosr = cos(roll * d2r);
    sinr = sin(roll * d2r);

    /* Compute orbital radius and angular momentum */
    rsc = vecnorm(pos);
    crossprod(pos, vel, &omega_los);
    angmo = vecnorm(omega_los);

    /* Construct the LOS */
    u_los.x = -cosr * pos.x / rsc + sinr * omega_los.x / angmo;
    u_los.y = -cosr * pos.y / rsc + sinr * omega_los.y / angmo;
    u_los.z = -cosr * pos.z / rsc + sinr * omega_los.z / angmo;

    /* Compute LOS angular rate */
    omega_los.x = cosr * angmo / (rsc * rsc)*(cosr * omega_los.x / angmo + sinr * pos.x / rsc);
    omega_los.y = cosr * angmo / (rsc * rsc)*(cosr * omega_los.y / angmo + sinr * pos.y / rsc);
    omega_los.z = cosr * angmo / (rsc * rsc)*(cosr * omega_los.z / angmo + sinr * pos.z / rsc);

    /* Compute orbit to Earth scale factor */
    scale = pos.x / e_a * pos.x / e_a + pos.y / e_a * pos.y / e_a + pos.z / e_b * pos.z / e_b;
    sang = u_los.x / e_a * u_los.x / e_a + u_los.y / e_a * u_los.y / e_a + u_los.z / e_b * u_los.z / e_b;
    sang = (-pos.x / e_a * u_los.x / e_a - pos.y / e_a * u_los.y / e_a - pos.z / e_b * u_los.z / e_b) / sqrt(sang * scale);
    sang = MIN(MAX(-1.0, sang), 1.0);
    sang = acos(sang);
    cosr = cos(sang);
    sinr = sin(sang);
    scale = 1.0 / scale - sinr*sinr;
    if (scale > 0.0) {
        scale = sqrt(scale);
    } else {
        printf("Roll angle too large, LOS misses Earth\n");
        return (-1);
    }
    rho = rsc * (cosr - scale);
    rho_dot = dotprod(pos, vel) / rsc * (cosr + sinr * sinr / scale);

    /* LOS rate cross LOS */
    crossprod(omega_los, u_los, &o_x_u);

    /* Construct ground target point */
    rp.x = pos.x + rho * u_los.x;
    rp.y = pos.y + rho * u_los.y;
    rp.z = pos.z + rho * u_los.z;
    vp.x = vel.x + rho_dot * u_los.x + rho * o_x_u.x;
    vp.y = vel.y + rho_dot * u_los.y + rho * o_x_u.y;
    vp.z = vel.z + rho_dot * u_los.z + rho * o_x_u.z;

    /* Construct the Earth rotation vector */
    er.x = 0.0;
    er.y = 0.0;
    er.z = e_rate;

    /* Compute er cross rp and ground target velocity */
    crossprod(er, rp, &er_x_rp);
    vgt.x = vp.x; // - er_x_rp.x;
    vgt.y = vp.y; // - er_x_rp.y;
    vgt.z = vp.z; // - er_x_rp.z;
    crossprod(u_los, vgt, &u_x_vgt);
    unitvec(u_x_vgt, &y);
    crossprod(y, u_los, &x);

    /* Construct ECI2YSF orientation matrix */
    a[0][0] = x.x;
    a[0][1] = x.y;
    a[0][2] = x.z;
    a[1][0] = y.x;
    a[1][1] = y.y;
    a[1][2] = y.z;
    a[2][0] = u_los.x;
    a[2][1] = u_los.y;
    a[2][2] = u_los.z;

    /* Rotate input LOS vector using the orientation matrix */
    u_los.x = a[0][0] * i_los.x + a[1][0] * i_los.y + a[2][0] * i_los.z;
    u_los.y = a[0][1] * i_los.x + a[1][1] * i_los.y + a[2][1] * i_los.z;
    u_los.z = a[0][2] * i_los.x + a[1][2] * i_los.y + a[2][2] * i_los.z;

    /* Compute orbit to Earth scale factor */
    scale = pos.x / e_a * pos.x / e_a + pos.y / e_a * pos.y / e_a + pos.z / e_b * pos.z / e_b;
    sang = u_los.x / e_a * u_los.x / e_a + u_los.y / e_a * u_los.y / e_a + u_los.z / e_b * u_los.z / e_b;
    sang = (-pos.x / e_a * u_los.x / e_a - pos.y / e_a * u_los.y / e_a - pos.z / e_b * u_los.z / e_b) / sqrt(sang * scale);
    sang = MIN(MAX(-1.0, sang), 1.0);
    sang = acos(sang);
    cosr = cos(sang);
    sinr = sin(sang);
    scale = 1.0 / scale - sinr*sinr;
    if (scale > 0.0) {
        scale = sqrt(scale);
    } else {
        printf("Roll angle too large, LOS misses Earth");
        return (-1);
    }
    rho = rsc * (cosr - scale);

    /* Construct ground target point */
    rp.x = pos.x + rho * u_los.x;
    rp.y = pos.y + rho * u_los.y;
    rp.z = pos.z + rho * u_los.z;

    /* Convert ground point to lat/lon */
    *gp_lat = atan(rp.z / sqrt(rp.x * rp.x + rp.y * rp.y) * e_a / e_b * e_a / e_b) / d2r;
    *gp_lon = atan2(rp.y, rp.x) / d2r;

    return (0);
}

/*
 * Name:	frame_scene.c
 *
 * Purpose:	Calculate the nominal scene frame and adjust the longitude
 *		to fit the observed metadata frame.
 *
 */

int frame_scene(
        WRS2 wrsparms, /* I: Structure of WRS-2 parameters */
        SMETA smeta, /* I: Input scene metadata structure */
        double *dellon, /* O: WRS longitude offset */
        double *delrow) /* O: Scene length in fractional WRS rows */ {
    int band; /* Current band number */
    int band_index; /* Index of current band in metadata structure */
    int sca; /* Current SCA number */
    double ndet; /* Normalized detector coordinate */
    double start_row; /* Fractional row at scene start */
    double stop_row; /* Fractional row at scene end */
    double max_line; /* Maximum product line number */
    double max_samp; /* Maximum product sample number */
    double l1t_line[2]; /* Ground point L1T line coordinate */
    double l1t_samp[2]; /* Ground point L1t sample coordinate */
    double ul_ls[2]; /* Upper left corner line/sample */
    double ur_ls[2]; /* Upper right corner line/sample */
    double lr_ls[2]; /* Lower right corner line/sample */
    double ll_ls[2]; /* Lower left corner line/sample */
    double dr_right; /* Scene size on right side */
    double dr_left; /* Scene size on left side */
    double sc_lat; /* Scene center latitude */
    double sc_lon; /* Scene center longitude */
    double dlon; /* Longitude offset */
    double radius; /* Earth radius in parallel of latitude */
    double e2; /* Earth eccentricity squared */
    double d2r; /* Conversion from degrees to radians */

    /* Compute conversion constant */
    d2r = atan(1.0) / 45.0;

    /* Construct the nominal scene frame */
    /* Find the upper left corner */
    band = 9;
    sca = 1;
    ndet = -1.0;
    /* Sudipta new addition for OLI to extract SCA and Det numbers */
    if ((band_index = smeta_band_number_to_index(&smeta, band)) < 0) {
        band = 10;
        sca = 3;
        ndet = 1.0; /* TIRS only */
    }
    if ((band_index = smeta_band_number_to_index(&smeta, band)) < 0) {
        printf("Scene framing bands 9 and 10 are not available.\n");
        return (ERROR);
    }
    /* End Sudipta new addition for OLI */
    start_row = (double) smeta.wrs_row - 0.7;
    if (project_point(wrsparms, smeta, *dellon, band, sca, start_row, ndet, &l1t_line[0], &l1t_samp[0]) != SUCCESS) {
        printf("Error projecting scene UL corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }
    start_row = (double) smeta.wrs_row - 0.6;
    if (project_point(wrsparms, smeta, *dellon, band, sca, start_row, ndet, &l1t_line[1], &l1t_samp[1]) != SUCCESS) {
        printf("Error projecting scene UL corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }
    start_row = (double) smeta.wrs_row - (0.7 * l1t_line[1] - 0.6 * l1t_line[0]) / (l1t_line[1] - l1t_line[0]);
    if (project_point(wrsparms, smeta, *dellon, band, sca, start_row, ndet, &ul_ls[0], &ul_ls[1]) != SUCCESS) {
        printf("Error projecting scene UL corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }

    /* Find the lower right corner */
    band = 9;
    sca = 14;
    ndet = 1.0;
    band_index = smeta_band_number_to_index(&smeta, band);
    max_line = (double) (smeta.band_smeta[band_index].l1t_lines - 1);
    max_samp = (double) (smeta.band_smeta[band_index].l1t_samps - 1);

    stop_row = (double) smeta.wrs_row + 0.6;
    if (project_point(wrsparms, smeta, *dellon, band, sca, stop_row, ndet, &l1t_line[0], &l1t_samp[0]) != SUCCESS) {
        printf("Error projecting scene LR corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }
    stop_row = (double) smeta.wrs_row + 0.7;
    if (project_point(wrsparms, smeta, *dellon, band, sca, stop_row, ndet, &l1t_line[1], &l1t_samp[1]) != SUCCESS) {
        printf("Error projecting scene LR corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }
    stop_row = (double) smeta.wrs_row + (0.6 * (l1t_line[1] - max_line) + 0.7 * (max_line - l1t_line[0])) / (l1t_line[1] - l1t_line[0]);
    if (project_point(wrsparms, smeta, *dellon, band, sca, stop_row, ndet, &lr_ls[0], &lr_ls[1]) != SUCCESS) {
        printf("Error projecting scene LR corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }

    /* Find provisional upper right corner */
    if (project_point(wrsparms, smeta, *dellon, band, sca, start_row, ndet, &ur_ls[0], &ur_ls[1]) != SUCCESS) {
        printf("Error projecting scene UR corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }

    /* Find the corner of the last odd SCA */
    sca = 13;
    if (project_point(wrsparms, smeta, *dellon, band, sca, start_row, ndet, &l1t_line[0], &l1t_samp[0]) != SUCCESS) {
        printf("Error projecting scene UR corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }

    /* Calculate the size of the scene, in rows */
    dr_right = (start_row - stop_row)
            * ((ur_ls[0] - l1t_line[0])*(ul_ls[1] - l1t_samp[0]) - (ur_ls[1] - l1t_samp[0])*(ul_ls[0] - l1t_line[0]))
            / ((ur_ls[1] - lr_ls[1])*(ul_ls[0] - l1t_line[0]) - (ur_ls[0] - lr_ls[0])*(ul_ls[1] - l1t_samp[0]));
    dr_right = stop_row - start_row - dr_right;

    /* Locate the provisional LL corner */
    sca = 1;
    ndet = -1.0;
    if (project_point(wrsparms, smeta, *dellon, band, sca, stop_row, ndet, &ll_ls[0], &ll_ls[1]) != SUCCESS) {
        printf("Error projecting scene LL corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }

    /* Find the corner of the first even SCA */
    sca = 2;
    if (project_point(wrsparms, smeta, *dellon, band, sca, stop_row, ndet, &l1t_line[0], &l1t_samp[0]) != SUCCESS) {
        printf("Error projecting scene LL corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }

    /* Calculate the size of the scene, in rows */
    dr_left = (start_row - stop_row)
            * ((ll_ls[0] - l1t_line[0])*(lr_ls[1] - l1t_samp[0]) - (ll_ls[1] - l1t_samp[0])*(lr_ls[0] - l1t_line[0]))
            / ((ll_ls[1] - ul_ls[1])*(lr_ls[0] - l1t_line[0]) - (ll_ls[0] - ul_ls[0])*(lr_ls[1] - l1t_samp[0]));
    dr_left = stop_row - start_row - dr_left;

    /* Combine the results */
    *delrow = (dr_left + dr_right) / 2.0;

    /* Locate the UR corner from the LR corner */
    sca = 14;
    ndet = 1.0;
    if (project_point(wrsparms, smeta, *dellon, band, sca, stop_row - (*delrow), ndet, &ur_ls[0], &ur_ls[1]) != SUCCESS) {
        printf("Error projecting scene UR corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }

    /* Locate the LL corner from the UL corner */
    sca = 1;
    ndet = -1.0;
    if (project_point(wrsparms, smeta, *dellon, band, sca, start_row + (*delrow), ndet, &ll_ls[0], &ll_ls[1]) != SUCCESS) {
        printf("Error projecting scene LL corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }

    /* Calculate longitude offset in meters */
    dlon = smeta.band_smeta[band_index].pixsize * (max_samp - MAX(MAX(ul_ls[1], ur_ls[1]), MAX(ll_ls[1], lr_ls[1]))
            - MIN(MIN(ul_ls[1], ur_ls[1]), MIN(ll_ls[1], lr_ls[1]))) / 2.0;

    /* Convert to longitude offset in degrees */
    pathrow_to_latlon(wrsparms, smeta.wrs_path, (double) smeta.wrs_row, &sc_lat, &sc_lon);
    e2 = 1.0 - wrsparms.semiminor / wrsparms.semimajor * wrsparms.semiminor / wrsparms.semimajor;
    radius = cos(d2r * sc_lat) * wrsparms.semimajor / sqrt(1.0 - e2 * sin(d2r * sc_lat) * sin(d2r * sc_lat));
    *dellon = *dellon + dlon / radius / d2r;

    return (SUCCESS);
}

/* another variant for frame_scene for older LS */
int frame_scene_ls(
        WRS2 wrsparms, /* I: Structure of WRS-2 parameters */
        SMETA smeta, /* I: Input scene metadata structure */
        double *dellon, /* O: WRS longitude offset */
        double *delrow) /* O: Scene length in fractional WRS rows */ {
    int band; /* Current band number */
    int band_index; /* Index of current band in metadata structure */
    double nsamp; /* Sample coordinate */
    double nline; /* Line coordinate */
    double start_row; /* Fractional row at scene start */
    double stop_row; /* Fractional row at scene end */
    double max_line; /* Maximum product line number */
    double max_samp; /* Maximum product sample number */
    double l1t_line[2]; /* Ground point L1T line coordinate */
    double l1t_samp[2]; /* Ground point L1t sample coordinate */
    double ul_ls[2]; /* Upper left corner line/sample */
    double ur_ls[2]; /* Upper right corner line/sample */
    double lr_ls[2]; /* Lower right corner line/sample */
    double ll_ls[2]; /* Lower left corner line/sample */
    double sc_lat; /* Scene center latitude */
    double sc_lon; /* Scene center longitude */
    double dlon; /* Longitude offset */
    double radius; /* Earth radius in parallel of latitude */
    double e2; /* Earth eccentricity squared */
    double d2r; /* Conversion from degrees to radians */

    /* Compute conversion constant */
    d2r = atan(1.0) / 45.0;

    /* Construct the nominal scene frame */
    /* Find the upper left corner */
    /* if( smeta.sensor_id == IAS_MSS &&
        (smeta.spacecraft_id == IAS_L1 || smeta.spacecraft_id == IAS_L2 || smeta.spacecraft_id == IAS_L3) )
        band = 4;
    else if( smeta.sensor_id == IAS_MSS &&
        (smeta.spacecraft_id == IAS_L4 || smeta.spacecraft_id == IAS_L5) )
        band = 1;
    else if( (smeta.spacecraft_id == IAS_L7) || (smeta.spacecraft_id == IAS_L5) ) */
    band = 4;
    /* else if( smeta.spacecraft_id == IAS_L4 )
        band = 0;
    else
    {
        printf("Error in sensor id.\n");
        smeta_release_projection();
        return(ERROR);
    } */

    nsamp = -1.0;
    nline = -1.0;
    start_row = (double) smeta.wrs_row - 0.6;
    if (project_point_ls(wrsparms, smeta, *dellon, band, start_row, nsamp, nline, &l1t_line[0], &l1t_samp[0]) != SUCCESS) {
        printf("Error projecting scene UL corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }
    start_row = (double) smeta.wrs_row - 0.5;
    if (project_point_ls(wrsparms, smeta, *dellon, band, start_row, nsamp, nline, &l1t_line[1], &l1t_samp[1]) != SUCCESS) {
        printf("Error projecting scene UL corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }
    start_row = (double) smeta.wrs_row - (0.6 * l1t_line[1] - 0.5 * l1t_line[0]) / (l1t_line[1] - l1t_line[0]);
    // start_row = (double)smeta.wrs_row - (0.8*l1t_line[1] - 0.5*l1t_line[0])/(l1t_line[1] - l1t_line[0]);
    if (project_point_ls(wrsparms, smeta, *dellon, band, start_row, nsamp, nline, &ul_ls[0], &ul_ls[1]) != SUCCESS) {
        printf("Error projecting scene UL corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }

    /* Find the lower right corner */
    nsamp = +1.0;
    nline = +1.0;
    band_index = smeta_band_number_to_index(&smeta, band);
    max_line = (double) (smeta.band_smeta_ls[band_index].l1t_lines - 1);
    max_samp = (double) (smeta.band_smeta_ls[band_index].l1t_samps - 1);

    stop_row = (double) smeta.wrs_row + 0.5;
    if (project_point_ls(wrsparms, smeta, *dellon, band, stop_row, nsamp, nline, &l1t_line[0], &l1t_samp[0]) != SUCCESS) {
        printf("Error projecting scene LR corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }
    stop_row = (double) smeta.wrs_row + 0.6;
    if (project_point_ls(wrsparms, smeta, *dellon, band, stop_row, nsamp, nline, &l1t_line[1], &l1t_samp[1]) != SUCCESS) {
        printf("Error projecting scene LR corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }
    stop_row = (double) smeta.wrs_row + (0.5 * (l1t_line[1] - max_line) + 0.6 * (max_line - l1t_line[0])) / (l1t_line[1] - l1t_line[0]);
    // stop_row = (double)smeta.wrs_row + (0.5*(l1t_line[1] - max_line) + 0.8*(max_line - l1t_line[0]))/(l1t_line[1] - l1t_line[0]);
    if (project_point_ls(wrsparms, smeta, *dellon, band, stop_row, nsamp, nline, &lr_ls[0], &lr_ls[1]) != SUCCESS) {
        printf("Error projecting scene LR corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }

    /* Find provisional upper right corner */
    nsamp = +1.0;
    nline = -1.0;
    if (project_point_ls(wrsparms, smeta, *dellon, band, start_row, nsamp, nline, &ur_ls[0], &ur_ls[1]) != SUCCESS) {
        printf("Error projecting scene UR corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }

    /* Locate the provisional LL corner */
    nsamp = -1.0;
    nline = +1.0;
    if (project_point_ls(wrsparms, smeta, *dellon, band, stop_row, nsamp, nline, &ll_ls[0], &ll_ls[1]) != SUCCESS) {
        printf("Error projecting scene LL corner to L1T line/sample.\n");
        smeta_release_projection();
        return (ERROR);
    }

    /* Combine the results */
    *delrow = stop_row - start_row;

    /* Calculate longitude offset in meters */
    dlon = smeta.band_smeta_ls[band_index].pixsize * (max_samp - MAX(MAX(ul_ls[1], ur_ls[1]), MAX(ll_ls[1], lr_ls[1]))
            - MIN(MIN(ul_ls[1], ur_ls[1]), MIN(ll_ls[1], lr_ls[1]))) / 2.0;

    /* Convert to longitude offset in degrees */
    pathrow_to_latlon(wrsparms, smeta.wrs_path, (double) smeta.wrs_row, &sc_lat, &sc_lon);
    e2 = 1.0 - wrsparms.semiminor / wrsparms.semimajor * wrsparms.semiminor / wrsparms.semimajor;
    radius = cos(d2r * sc_lat) * wrsparms.semimajor / sqrt(1.0 - e2 * sin(d2r * sc_lat) * sin(d2r * sc_lat));
    *dellon = *dellon + dlon / radius / d2r;

    return (SUCCESS);
}

int project_point(
        WRS2 wrsparms, /* I: Structure of WRS-2 parameters */
        SMETA smeta, /* I: Input scene metadata structure */
        double dellon, /* I: WRS longitude offset */
        int band, /* I: Band number (user band) */
        int sca, /* I: SCA number (1 to nsca) */
        double frow, /* I: Fractional WRS row coordinate */
        double ndet, /* I: Normalized detector coordinate */
        double *l1t_line, /* O: Projected line coordinate */
        double *l1t_samp) /* O: Projected sample coordinate */ {
    VECTOR pos; /* Spacecraft ECEF position vector */
    VECTOR vel; /* Spacecraft ECEF velocity vector */
    VECTOR slos; /* Spacecraft line-of-sight vector */
    int band_index; /* Index of current band in metadata structure */
    double gp_lat; /* Ground point latitude */
    double gp_lon; /* Ground point longitude */
    double proj_x; /* Ground point projection X coordinate */
    double proj_y; /* Ground point projection Y coordinate */
    double d2r; /* Conversion from degrees to radians */

    /* Compute conversion constant */
    d2r = atan(1.0) / 45.0;

    /* Find the nominal ephemeris state vector */
    pathrow_to_posvel(wrsparms, smeta.wrs_path, dellon, frow, &pos, &vel);

    /* Find the current band in the metadata */
    band_index = smeta_band_number_to_index(&smeta, band);

    /* Calculate the instrument LOS vector */
    if (calc_los(smeta.band_smeta[band_index], sca, ndet, &slos) != SUCCESS) {
        printf("Error generating LOS vector.\n");
        return (ERROR);
    }

    /* Project to ground latitude/longitude */
    if (calc_yaw_steering_gp(wrsparms, pos, vel, slos, smeta.roll_angle, &gp_lat, &gp_lon) < 0) {
        printf("Error projecting nominal scene center.\n");
        return (ERROR);
    }
    gp_lat *= d2r;
    gp_lon *= d2r;

    /* Apply scene map projection */
    if (smeta_geodetic_to_proj(smeta.projection, gp_lat, gp_lon, &proj_x, &proj_y) != SUCCESS) {
        printf("Error converting scene center lat/lon to projection X/Y.\n");
        return (ERROR);
    }

    /* Reduce to L1T line/sample */
    if (smeta_proj_to_l1t(&smeta, band, proj_x, proj_y, l1t_line, l1t_samp) != SUCCESS) {
        printf("Errro converting scene center X/Y to L1T line/sample.\n");
        return (ERROR);
    }

    return (SUCCESS);
}

/* Another variant for project_point for older LS */
int project_point_ls(
        WRS2 wrsparms, /* I: Structure of WRS-2 parameters */
        SMETA smeta, /* I: Input scene metadata structure */
        double dellon, /* I: WRS longitude offset */
        int band, /* I: Band number (user band) */
        double frow, /* I: Fractional WRS row coordinate */
        double npix, /* I: Normalized detector coordinate */
        double ndet, /* I: Line detector coordinate */
        double *l1t_line, /* O: Projected line coordinate */
        double *l1t_samp) /* O: Projected sample coordinate */ {
    VECTOR pos; /* Spacecraft ECEF position vector */
    VECTOR vel; /* Spacecraft ECEF velocity vector */
    VECTOR slos; /* Spacecraft line-of-sight vector */
    int band_index; /* Index of current band in metadata structure */
    double gp_lat; /* Ground point latitude */
    double gp_lon; /* Ground point longitude */
    double proj_x; /* Ground point projection X coordinate */
    double proj_y; /* Ground point projection Y coordinate */
    double d2r; /* Conversion from degrees to radians */

    /* Compute conversion constant */
    d2r = atan(1.0) / 45.0;

    /* Find the nominal ephemeris state vector */
    pathrow_to_posvel(wrsparms, smeta.wrs_path, dellon, frow, &pos, &vel);

    /* Find the current band in the metadata */
    band_index = smeta_band_number_to_index(&smeta, band);

    /* Calculate the instrument LOS vector */
    if (calc_los_ls(smeta.band_smeta_ls[band_index], frow, npix, ndet, &slos) != SUCCESS) {
        printf("Error generating LOS vector.\n");
        return (ERROR);
    }

    /* Project to ground latitude/longitude */
    if (calc_yaw_steering_gp_ls(wrsparms, pos, vel, slos, smeta.roll_angle, &gp_lat, &gp_lon) < 0) {
        printf("Error projecting nominal scene center.\n");
        return (ERROR);
    }
    gp_lat *= d2r;
    gp_lon *= d2r;

    /* Apply scene map projection */
    if (smeta_geodetic_to_proj(smeta.projection, gp_lat, gp_lon, &proj_x, &proj_y) != SUCCESS) {
        printf("Error converting scene center lat/lon to projection X/Y.\n");
        return (ERROR);
    }

    /* Reduce to L1T line/sample */
    if (smeta_proj_to_l1t(&smeta, band, proj_x, proj_y, l1t_line, l1t_samp) != SUCCESS) {
        printf("Errro converting scene center X/Y to L1T line/sample.\n");
        return (ERROR);
    }

    return (SUCCESS);
}

/*
 * Name:	get_ecfsun.c
 *
 * Purpose:	Convert the metadata sun angles into an ECEF solar vector.
 *		As a first approximation the sun is assumed not to change
 *		position during the scene as the actual change should be less
 *		than 10 arc-minutes.
 *
 */

int get_ecfsun(
        WRS2 wrsparms, /* Structure of WRS-2 parameters */
        SMETA smeta, /* Input scene metadata structure */
        double dellon, /* WRS longitude offset */
        VECTOR *ecfsun) /* Sun vector in ECEF coordinates */ {
    VECTOR pos; /* Spacecraft ECEF position vector */
    VECTOR vel; /* Spacecraft ECEF velocity vector */
    VECTOR los; /* Spacecraft line-of-sight vector */
    double drow; /* Fractional WRS row */
    double gp_lat; /* Ground point latitude */
    double gp_lon; /* Ground point longitude */
    double d2r; /* Conversion from degrees to radians */
    double e2; /* WGS84 eccentricity squared */
    double ecf2lsr[3][3]; /* ECEF to Local Space Rectangular transformation matrix */
    VECTOR lsrsun; /* Sun vector in LSR coordinates */

    /* Compute conversion constant */
    d2r = atan(1.0) / 45.0;
    e2 = 1.0 - wrsparms.semiminor / wrsparms.semimajor * wrsparms.semiminor / wrsparms.semimajor;

    /* Calculate the scene center latitude/longitude */
    drow = (double) smeta.wrs_row;
    pathrow_to_posvel(wrsparms, smeta.wrs_path, dellon, drow, &pos, &vel);
    los.x = 0.0;
    los.y = 0.0;
    los.z = 1.0; /* boresight */
    if (calc_yaw_steering_gp(wrsparms, pos, vel, los, smeta.roll_angle, &gp_lat, &gp_lon) < 0) {
        printf("Error projecting nominal scene center.\n");
        return (ERROR);
    }
    gp_lat *= d2r;
    gp_lon *= d2r;

    /* Convert latitude to geocentric */
    /* This looks wrong but this is how the metadata angles
       are calculated, so we need to back out the geocentric 
       - geodetic effect */
    gp_lat = atan(tan(gp_lat) * (1.0 - e2));

    /* Compute ECEF solar vector */
    lsrsun.x = cos(d2r * smeta.sun_elevation) * sin(d2r * smeta.sun_azimuth);
    lsrsun.y = cos(d2r * smeta.sun_elevation) * cos(d2r * smeta.sun_azimuth);
    lsrsun.z = sin(d2r * smeta.sun_elevation);
    (void) smeta_geodetic_to_ecf2lsr(gp_lat, gp_lon, ecf2lsr);
    ecfsun->x = ecf2lsr[0][0] * lsrsun.x + ecf2lsr[1][0] * lsrsun.y + ecf2lsr[2][0] * lsrsun.z;
    ecfsun->y = ecf2lsr[0][1] * lsrsun.x + ecf2lsr[1][1] * lsrsun.y + ecf2lsr[2][1] * lsrsun.z;
    ecfsun->z = ecf2lsr[0][2] * lsrsun.x + ecf2lsr[1][2] * lsrsun.y + ecf2lsr[2][2] * lsrsun.z;

    return (SUCCESS);
}

int get_ecfsun_ls(
        WRS2 wrsparms, /* Structure of WRS-2 parameters */
        SMETA smeta, /* Input scene metadata structure */
        double dellon, /* WRS longitude offset */
        VECTOR *ecfsun) /* Sun vector in ECEF coordinates */ {
    VECTOR pos; /* Spacecraft ECEF position vector */
    VECTOR vel; /* Spacecraft ECEF velocity vector */
    VECTOR los; /* Spacecraft line-of-sight vector */
    double drow; /* Fractional WRS row */
    double gp_lat; /* Ground point latitude */
    double gp_lon; /* Ground point longitude */
    double d2r; /* Conversion from degrees to radians */
    double e2; /* WGS84 eccentricity squared */
    double ecf2lsr[3][3]; /* ECEF to Local Space Rectangular transformation matrix */
    VECTOR lsrsun; /* Sun vector in LSR coordinates */

    /* Compute conversion constant */
    d2r = atan(1.0) / 45.0;
    e2 = 1.0 - wrsparms.semiminor / wrsparms.semimajor * wrsparms.semiminor / wrsparms.semimajor;

    /* Calculate the scene center latitude/longitude */
    drow = (double) smeta.wrs_row;
    pathrow_to_posvel(wrsparms, smeta.wrs_path, dellon, drow, &pos, &vel);
    los.x = 0.0;
    los.y = 0.0;
    los.z = 1.0; /* boresight */
    if (calc_yaw_steering_gp_ls(wrsparms, pos, vel, los, smeta.roll_angle, &gp_lat, &gp_lon) < 0) {
        printf("Error projecting nominal scene center.\n");
        return (ERROR);
    }
    gp_lat *= d2r;
    gp_lon *= d2r;

    /* Convert latitude to geocentric */
    /* This looks wrong but this is how the metadata angles
       are calculated, so we need to back out the geocentric
       - geodetic effect */
    gp_lat = atan(tan(gp_lat) * (1.0 - e2));

    /* Compute ECEF solar vector */
    lsrsun.x = cos(d2r * smeta.sun_elevation) * sin(d2r * smeta.sun_azimuth);
    lsrsun.y = cos(d2r * smeta.sun_elevation) * cos(d2r * smeta.sun_azimuth);
    lsrsun.z = sin(d2r * smeta.sun_elevation);
    (void) smeta_geodetic_to_ecf2lsr(gp_lat, gp_lon, ecf2lsr);
    ecfsun->x = ecf2lsr[0][0] * lsrsun.x + ecf2lsr[1][0] * lsrsun.y + ecf2lsr[2][0] * lsrsun.z;
    ecfsun->y = ecf2lsr[0][1] * lsrsun.x + ecf2lsr[1][1] * lsrsun.y + ecf2lsr[2][1] * lsrsun.z;
    ecfsun->z = ecf2lsr[0][2] * lsrsun.x + ecf2lsr[1][2] * lsrsun.y + ecf2lsr[2][2] * lsrsun.z;

    return (SUCCESS);
}

/*
 * Name:	pathrow_to_latlon.c
 *
 * Purpose:	Computes the latitude and longitude corresponding to an input WRS-2
 *              path/row including fractional row.
 *
 */

void pathrow_to_latlon(
        WRS2 parms, /* I: WRS2 system parameters */
        int path, /* I: WRS2 path */
        double wrsrow, /* I: WRS2 fractional row */
        double *lat, /* O: Nominal latitude */
        double *lon) /* O: Nominal longitude */ {
    double cta; /* Central travel angle */
    double dnodelon; /* Longitude of the descending node */
    double dlon; /* Longitude offset from descending node */
    double gc_lat; /* Geocentric latitude */
    double d2r; /* Degrees to radians conversion */

    /* Initialize constants */
    d2r = atan(1.0) / 45.0;

    /* Calculate CTA from row */
    cta = 360.0 * (wrsrow - (double) parms.row0lat) / (double) parms.numrow;

    /* Calculate latitude from CTA */
    gc_lat = asin(-sin(d2r * cta) * sin(d2r * parms.orbitincl));
    *lat = atan(tan(gc_lat) * parms.semimajor / parms.semiminor * parms.semimajor / parms.semiminor) / d2r;

    /* Calculate longitude of descending node from path and CTA */
    dnodelon = 360.0 + 180.0 + parms.path1long - (double) (path - 1)*360.0 / (double) parms.numpath;
    dnodelon = fmod(dnodelon, 360.0) - 180.0 - cta * (double) parms.cycle / (double) parms.numpath;

    /* Calculate the longitude offset */
    dlon = atan2(tan(d2r * gc_lat) / tan(parms.orbitincl), cos(d2r * cta) / cos(d2r * gc_lat)) / d2r;
    *lon = dnodelon - dlon;
    if (*lon > 180.0) *lon -= 360.0;
    if (*lon < -180.0) *lon += 360.0;

    return;
}

/*
 * Name:	pathrow_to_posvel.c
 *
 * Purpose:	Computes the nominal satellite position and velocity (ECEF) vectors
 *              that correspond to the specified WRS path/row including fractional row
 *              and path longitude adjustments.
 *
 */

void pathrow_to_posvel(
        WRS2 parms, /* I: WRS2 system parameters */
        int path, /* I: WRS2 path */
        double dellon, /* I: Path longitude adjustment (in degrees) */
        double wrsrow, /* I: WRS2 fractional row */
        VECTOR *pos, /* O: Nominal position vector */
        VECTOR *vel) /* O: Nominal velocity vector */ {
    VECTOR on; /* Orbit normal vector */
    VECTOR dnode; /* Descending node vector */
    VECTOR y; /* In-plane normal vector */
    double cta; /* Central travel angle */
    double dnodelon; /* Longitude of the descending node */
    double orate; /* Orbital angular rate */
    double d2r; /* Degrees to radians conversion */

    /* Initialize constants */
    d2r = atan(1.0) / 45.0;

    /* Calculate CTA from row */
    cta = 360.0 * (wrsrow - (double) parms.row0lat) / (double) parms.numrow;

    /* Calculate longitude of descending node from path and CTA */
    dnodelon = 360.0 + 180.0 + parms.path1long - (double) (path - 1)*360.0 / (double) parms.numpath + dellon;
    dnodelon = fmod(dnodelon, 360.0) - 180.0 - cta * (double) parms.cycle / (double) parms.numpath;

    /* Compute the orbit normal vector */
    on.x = -sin(d2r * parms.orbitincl) * sin(d2r * dnodelon);
    on.y = sin(d2r * parms.orbitincl) * cos(d2r * dnodelon);
    on.z = cos(d2r * parms.orbitincl);

    /* Compute the descending node vector */
    dnode.x = cos(d2r * dnodelon);
    dnode.y = sin(d2r * dnodelon);
    dnode.z = 0.0;
    crossprod(on, dnode, &y);

    /* Construct the position vector */
    pos->x = (dnode.x * cos(d2r * cta) + y.x * sin(d2r * cta)) * parms.orbitrad;
    pos->y = (dnode.y * cos(d2r * cta) + y.y * sin(d2r * cta)) * parms.orbitrad;
    pos->z = (dnode.z * cos(d2r * cta) + y.z * sin(d2r * cta)) * parms.orbitrad;

    /* Construct the velocity vector */
    orate = 360.0 * (double) parms.numpath / (double) parms.cycle / 86400.0 * d2r;
    on.x *= orate;
    on.y *= orate;
    on.z *= orate;
    crossprod(on, *pos, vel);

    return;
}

/*
 * Name:	vector_math.c
 *
 * Purpose:	Code units to calculate vector functions.
 *
 */

double dotprod(/* Returns dot product of two VECTOR structures */
        VECTOR a, /* I: First input vector */
        VECTOR b) /* I: Second input vector */ {
    return ( a.x * b.x + a.y * b.y + a.z * b.z);
}

void crossprod(/* Computes vector cross product */
        VECTOR a, /* I: First input vector */
        VECTOR b, /* I: Second input vector */
        VECTOR *c) /* O: Output vector */ {
    c->x = a.y * b.z - a.z * b.y;
    c->y = a.z * b.x - a.x * b.z;
    c->z = a.x * b.y - a.y * b.x;
    return;
}

double vecnorm(/* Returns vector magnitude */
        VECTOR a) /* I: Input vector */ {
    return ( sqrt(dotprod(a, a)));
}

void unitvec(/* Normalizes input vector */
        VECTOR a, /* I: Input vector */
        VECTOR *b) /* O: Normalized output vector */ {
    double mag; /* Vector magnitude */

    mag = vecnorm(a);
    if (mag > 0.0) {
        b->x = a.x / mag;
        b->y = a.y / mag;
        b->z = a.z / mag;
    } else {
        b->x = a.x;
        b->y = a.y;
        b->z = a.z;
    }
    return;
}

void rotatevec(
        double mat[3][3], /* 3-by-3 rotation matrix */
        VECTOR a, /* Original vector */
        VECTOR *b) /* Rotated vector */ {
    b->x = mat[0][0] * a.x + mat[0][1] * a.y + mat[0][2] * a.z;
    b->y = mat[1][0] * a.x + mat[1][1] * a.y + mat[1][2] * a.z;
    b->z = mat[2][0] * a.x + mat[2][1] * a.y + mat[2][2] * a.z;

    return;
}
