#include "l1.h"
#include <nc4utils.h>
#include <hdf4utils.h>

int get_sdps_cld_mask( l1str *l1rec, int32_t ip, char *cld_category ){
/*
 get_sdps_cld_mask - get the cloud mask catergoy, and the cloud/no cloud 
    from the cloud mask file made for SDPS
 
  Returns: int the cloud value: 0 no cloud, 1 cloud
  l1str *l1rec:  supply the cloud mask file name, controls
  int32_t ip:  pixel #
  char *cld_category:  the cloud category: 0 - confident_cloudy,
      1 - probablyclear, 2 - probablyclear, 3 - probablyclear
  this will return the above category even for the new format, although it returns 
  only 0 for cloud, 3 for clear

  WDR, 29 Mar 2021
  WDR, 9 Nov 2023, also be able to read the new format with SDS cloud_flag
       in group geophysical_data
 */
    static size_t start[] = {0,0}, count[] = {0,0}, fil_np;
    static int firstCall = 1;
    static int ncid, d_id, dim_id;
    static signed char* msk_lcl;
    static int lastScan = -1;
    static int new_cl_msk_fmt;
    int grp_id;

    if (firstCall) {
        char sdsname[HDF4_UTILS_MAX_NAME] = "";
        char file [FILENAME_MAX] = "";

        if (!l1_input->cld_msk_file[0]) {
            printf("-E- %s line %d: no cloud flag file provided\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        strcpy(file, l1_input->cld_msk_file);
        printf("Loading cloud mask flag file %s\n", file);

        /* Open the file */
        if (nc_open(file, 0, &ncid) != NC_NOERR) {
            fprintf(stderr,
                "-E- %s %d: file: %s is not netcdf, not acceptable cloud maskfile\n",
                __FILE__, __LINE__, file);
            exit(1);
        }
        // get the # pixels 
        if (nc_inq_dimid(ncid, "pixels_per_line", &dim_id ) != NC_NOERR) {
            fprintf(stderr, 
                "-E- %s %d: file: %s Could not read the pixels_per_line id\n",
                __FILE__, __LINE__, file);
            exit(1);
        }

        if (nc_inq_dimlen(ncid, dim_id, &fil_np ) != NC_NOERR) {
            fprintf(stderr,
                "-E- %s %d: file: %s Could not read the pixels_per_line\n",
                __FILE__, __LINE__, file);
            exit(1);
        }
        //  try reading the new group/sds name.  if group fails, 
        //  assume old format
        if( nc_inq_grp_ncid( ncid, "geophysical_data", &grp_id ) != NC_NOERR) {
            new_cl_msk_fmt = 0;
            fprintf(stderr,
                "-I- %s %d: file: %s Cloud mask read switch to old format\n",
                __FILE__, __LINE__, file);
            strcpy(sdsname, "CF_CATEGORY");
            if( nc_inq_varid(ncid, sdsname, &d_id ) != NC_NOERR) {
                fprintf(stderr, 
                    "-E- %s %d: file: %s Could not find SDS CF_CATEGORY\n",
                    __FILE__, __LINE__, file);
                exit(1);
            }
        } else {
            //   for new mask format 
            ncid = grp_id;
            new_cl_msk_fmt = 1;
            if(l1_input->cloud_mask_opt == 0) {
                strcpy(sdsname, "cloud_flag");
            } else if(l1_input->cloud_mask_opt == 1) {
                strcpy(sdsname, "cloud_flag_dilated");
            } else {
                fprintf(stderr,
                    "-E- %s:%d - cloud_mask_opt=%d Illegal option\n",
                    __FILE__, __LINE__, l1_input->cloud_mask_opt);
                exit(1);
            }
            if( nc_inq_varid(grp_id, sdsname, &d_id ) != NC_NOERR) {
                fprintf(stderr, 
                    "-E- %s %d: file: %s Could not find SDS CF_CATEGORY\n",
                    __FILE__, __LINE__, file);
                exit(1);
            }
        }
        //  get the range to process and set it
        start[1] = l1_input->spixl - 1;  // spixl, epixl are 1-origin!
        count[1] = l1_input->epixl - l1_input->spixl + 1;
        count[0] = 1;
        // start[0] will be the line to read
        if( ( msk_lcl = (signed char *) 
            malloc( count[1] * sizeof( signed char *) ) ) == NULL ) {
            fprintf(stderr,
                "-E- %s %d: file: %s Could not allocate cloud mask storage\n",
                __FILE__, __LINE__, file);
            exit(1);
        }
        firstCall = 0;
    //  end setup
    }

    if (lastScan != l1rec->iscan) {
        start[0] = l1rec->iscan;
        if( nc_get_vara_schar( ncid, d_id, start, count, msk_lcl ) != NC_NOERR ) {
            fprintf(stderr,
                "-E- %s %d: Could not read the cloud mask file line\n",
                __FILE__, __LINE__);
            exit(1);
        }
        lastScan = l1rec->iscan;
    }

    if( new_cl_msk_fmt == 0 ) {
        *cld_category = msk_lcl[ip];
    } else {
        if( msk_lcl[ip] == BAD_BYTE )
            *cld_category = BAD_BYTE;
        else if( msk_lcl[ip] == 0 )
            *cld_category = 3;
        else *cld_category = 0;
    }

    if (*cld_category == BAD_BYTE )
        return (BAD_BYTE);
    else if (*cld_category < 2 )
        return (1);
    else
        return (0);
}

char get_cloudmask_meris(l1str *l1rec, int32_t ip) {
    // Cloud Masking for MERIS and OLCI

    static int ib443, ib490, ib560, ib620, ib665, ib681, ib709, ib754, ib865, ib885, firstCall = 1;
    int ipb;
    float *rhos = l1rec->rhos, *cloud_albedo = l1rec->cloud_albedo;
    float ftemp, cldtmp;
    char flagcld;
    static productInfo_t* rhos_product_info;

    if (firstCall == 1) {
        ib443 = bindex_get(443);
        ib490 = bindex_get(490);
        ib560 = bindex_get(560);
        ib620 = bindex_get(620);
        ib665 = bindex_get(665);
        ib681 = bindex_get(681);
        ib709 = bindex_get(709);
        ib754 = bindex_get(754);
        ib865 = bindex_get(865);
        ib885 = bindex_get(885);

        if (ib443 < 0 || ib490 < 0 || ib560 < 0 ||  ib620 < 0 || ib665 < 0 || ib681 < 0 || ib709 < 0 || ib754 < 0 || ib865 < 0 || ib885 < 0) {
            printf("get_habs_cldmask: incompatible sensor wavelengths for this algorithm\n");
            exit(EXIT_FAILURE);
        }

        if (rhos_product_info == NULL) {
            rhos_product_info = allocateProductInfo();
            findProductInfo("rhos", l1rec->l1file->sensorID, rhos_product_info);
        }

        firstCall = 0;
    }
    flagcld = 0;
    ipb = l1rec->l1file->nbands*ip;

    if (rhos[ipb + ib443] >= 0.0 &&
        rhos[ipb + ib560] >= 0.0 &&
        rhos[ipb + ib620] >= 0.0 &&
        rhos[ipb + ib665] >= 0.0 &&
        rhos[ipb + ib681] >= 0.0 &&
        rhos[ipb + ib709] >= 0.0 &&
        rhos[ipb + ib754] >= 0.0) {
        // turbidity signal in water
        if ((rhos[ipb + ib443] <= rhos_product_info->validMax) &&
            (rhos[ipb + ib665] <= rhos_product_info->validMax) &&
            (rhos[ipb + ib620] <= rhos_product_info->validMax) &&
            (rhos[ipb + ib681] <= rhos_product_info->validMax) &&
            (rhos[ipb + ib754] <= rhos_product_info->validMax)) {

            ftemp = (rhos[ipb + ib620] + rhos[ipb + ib665] + rhos[ipb + ib681]) - 3 * rhos[ipb + ib443] - \
                    (rhos[ipb + ib754] - rhos[ipb + ib443]) / (754 - 443)*(620 + 665 + 681 - 3 * 443);
        } else {
            ftemp = 0;
        }
        //switch to rhos_865 where cldalb fails
        if (cloud_albedo[ip] >= 0.0) {
            cldtmp = cloud_albedo[ip];
        } else {
            cldtmp = rhos[ipb + ib865];
        }

        //remove turbidity signal from cloud albedo
        if (ftemp > 0) cldtmp = cldtmp - 3 * ftemp;

        if (cldtmp > 0.08) {
            flagcld = 1;
        }

        // to deal with scum look at relative of NIR and blue for lower albedos
        if ((rhos[ipb + ib754] <= rhos_product_info->validMax) &&
            (rhos[ipb + ib709] <= rhos_product_info->validMax) &&
            (rhos[ipb + ib681] <= rhos_product_info->validMax) &&
            (rhos[ipb + ib560] <= rhos_product_info->validMax)) {

            if ((rhos[ipb + ib443] <= rhos_product_info->validMax) &&
                (rhos[ipb + ib490] <= rhos_product_info->validMax)) {
                if ((rhos[ipb + ib754] + rhos[ipb + ib709]) > (rhos[ipb + ib443] + rhos[ipb + ib490]) && cldtmp < 0.1 && (rhos[ipb + ib560] > rhos[ipb + ib681])) flagcld = 0;
            }
            if (rhos[ipb + ib665] <= rhos_product_info->validMax) {
                if (((rhos[ipb + ib754] + rhos[ipb + ib709]) - (rhos[ipb + ib665] + rhos[ipb + ib681])) > 0.01 && cldtmp < 0.15 && (rhos[ipb + ib560] > rhos[ipb + ib681])) flagcld = 0;
                if ((((rhos[ipb + ib754] + rhos[ipb + ib709]) - (rhos[ipb + ib665] + rhos[ipb + ib681])) / cldtmp) > 0.1 && (rhos[ipb + ib560] > rhos[ipb + ib681])) flagcld = 0;
            }
        }
        if (rhos[ipb + ib665] <= rhos_product_info->validMax) {
            if ((rhos[ipb + ib665] > 0.1) && (cldtmp > 0.15)) {
                flagcld = 1;
            }
        }
        if ((rhos[ipb + ib865] >= rhos_product_info->validMin) &&
            (rhos[ipb + ib865] <= rhos_product_info->validMax)) {
            // flag high glint
            if (rhos[ipb + ib865] - cloud_albedo[ip] > 0.25) {
                flagcld = 1;
            }
        } else {
                flagcld = 1;
            }
        }

    if (l1rec->iscan == l1rec->l1file->nscan) {
        freeProductInfo(rhos_product_info);
    }

    return (flagcld);
}

char get_cloudmask_modis(l1str *l1rec, int32_t ip) {
    // Cloud Masking for MODIS
    static int firstCall = 1;
    static int ib469, ib555, ib645, ib667, ib859, ib1240, ib2130;
    int ipb;
    float *rhos = l1rec->rhos, *cloud_albedo = l1rec->cloud_albedo;
    float ftemp, ftemp2, ftemp3;
    float cloudthr = 0.027;
    char flagcld;

    if (firstCall == 1) {
        ib469 = bindex_get(469);
        ib555 = bindex_get(555);
        ib645 = bindex_get(645);
        ib667 = bindex_get(667);
        ib859 = bindex_get(859);
        ib1240 = bindex_get(1240);
        ib2130 = bindex_get(2130);
        if (ib469 < 0 || ib555 < 0 || ib645 < 0 || ib667 < 0 || ib859 < 0 || ib1240 < 0 || ib2130 < 0) {
            printf("get_habs_cldmask: incompatible sensor wavelengths for this algorithm\n");
            exit(EXIT_FAILURE);
        }
        firstCall = 0;
    }

    ipb = l1rec->l1file->nbands*ip;
    flagcld = 0;
    ftemp = 0; //rhos[ipb+ib667];
    //      first correct for turbid water

    if (rhos[ipb + ib667] < 0.0) ftemp = 0.0;
    ftemp2 = cloud_albedo[ip] - ftemp;

    if (ftemp2 > 0.027) flagcld = 1;

    //        non-water check  1240 is bright relative to 859 and the combination is bright
    //        this may hit glint by accident, need to be checked.

    if (rhos[ipb + ib1240] / rhos[ipb + ib859] > 0.5 && (rhos[ipb + ib1240] + rhos[ipb + ib2130]) > 0.10) flagcld = 1;


    //        now try to correct for glint
    //        region check was thrown out {IF (region = "OM") cloudthr = 0.04} rjh 11/2/2015

    ftemp = rhos[ipb + ib645] - rhos[ipb + ib555] + (rhos[ipb + ib555] - rhos[ipb + ib859])*(645.0 - 555.0) / (859.0 - 555.0);
    ftemp2 = cloud_albedo[ip] + ftemp;
    if (ftemp2 < cloudthr) flagcld = 0;
    if (rhos[ipb + ib859] / rhos[ipb + ib1240] > 4.0) flagcld = 0;

    //     scum areas

    if ((rhos[ipb + ib859] - rhos[ipb + ib469]) > 0.01 && cloud_albedo[ip] < 0.30) flagcld = 0;
    if ((rhos[ipb + ib859] - rhos[ipb + ib645]) > 0.01 && cloud_albedo[ip] < 0.15) flagcld = 0;
    if (rhos[ipb + ib1240] < 0.2)
        ftemp2 = ftemp2 - (rhos[ipb + ib859] - rhos[ipb + ib1240]) * fabs(rhos[ipb + ib859] - rhos[ipb + ib1240]) / cloudthr;

    ftemp3 = ftemp2;
    if (ftemp2 < cloudthr * 2) {
        if ((rhos[ipb + ib555] - rhos[ipb + ib1240]) > (rhos[ipb + ib469] - rhos[ipb + ib1240])) {
            ftemp3 = ftemp2 - (rhos[ipb + ib555] - rhos[ipb + ib1240]);
        } else {
            ftemp3 = ftemp2 - (rhos[ipb + ib469] - rhos[ipb + ib1240]);
        }
    }

    if (ftemp3 < cloudthr) flagcld = 0;

    return (flagcld);
}

char get_cldmask(l1str *l1rec, int32_t ip) {
    //function for cloud mask by pixel
    switch (l1rec->l1file->sensorID) {
    case MERIS:
    case OLCIS3A:
    case OLCIS3B:
        return (get_cloudmask_meris(l1rec, ip));
        break;
    case MODISA:
    case MODIST:
        return (get_cloudmask_modis(l1rec, ip));
        break;
    default:
        printf("HABS cldmsk not supported for this sensor (%s).\n",
                sensorId2SensorName(l1rec->l1file->sensorID));
        exit(EXIT_FAILURE);
    }
    return (0);
}
