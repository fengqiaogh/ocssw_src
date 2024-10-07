#include <mfhdf.h>
#include "ancil.h"
#include "ancnrt_proto.h"
#include "o3_toms.h"
#include <stdio.h>
#include <ctype.h>

int rd_toms_ascii(char *file, toms_txt_info_struc *toms_info)
/*******************************************************************

   rd_toms_ascii

   purpose: read the controls and data from a TOMS ascii text file

   Returns type: int - 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      text file from TOMS project
      toms_txt_info_struc *  toms_info   O      structure containing all the 
                                                information and data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 6 Dec 2013      replacement for rdgrid.f and more 
                                        flexible
 
 *******************************************************************/ {
    FILE *fid;
    int32_t doy, year, asc_node_h, asc_node_m, nlon, nlat, toms_typ;
    int32_t iv, ip, il, ilrev, end_ilin, end_olin, cptr;
    float slon, elon, del_lon, slat, elat, del_lat;
    char rstr[100];
    char inst_str[27], asc_node_tod[4];
    int16 *ptr;
    char *upcase();
    /*
     *  open the file and get the header information
     */
    if ((fid = fopen(file, "r")) == NULL) {
        printf("%s, %d E: Could not open TOMS text file: %s\n",
                __FILE__, __LINE__, file);
        return 1;
    }
    if (fgets(rstr, 90, fid) == NULL) {
        printf("%s, %d E: Unexpected end of grid file found, file: %s\n",
                __FILE__, __LINE__, file);
        return 1;
    }
    if ((sscanf(&rstr[6], "%3d", &doy) != 1) ||
            (sscanf(&rstr[18], "%4d", &year) != 1) ||
            (sscanf(&rstr[23], "%s", inst_str) != 1) ||
            (sscanf(&rstr[50], "%s", toms_info->gen_str) != 1) ||
            (sscanf(&rstr[71], "%2d", &asc_node_h) != 1) ||
            (sscanf(&rstr[74], "%2d", &asc_node_m) != 1) ||
            (sscanf(&rstr[77], "%s", asc_node_tod) != 1)) {
        printf("%s, %d E: Trouble reading header 1 of text file: %s\n",
                __FILE__, __LINE__, file);
        return 1;
    }
    toms_typ = TOMS_TYP_UNKNOWN;
    if (strcmp(inst_str, "OMI") == 0) toms_typ = TOMS_TYP_OMITOMS;
    if (strcmp(inst_str, "NIMBUS-7/TOMS") == 0) toms_typ = TOMS_TYP_N7TOMS;
    if (strcmp(inst_str, "EP/TOMS") == 0) toms_typ = TOMS_TYP_EPTOMS;

    toms_info->toms_typ = toms_typ;
    toms_info->year = year;
    toms_info->doy = doy;
    /*
     *  make the ascending node string
     */
    if (strncmp(upcase(asc_node_tod), "PM", 2) == 0)
        asc_node_h += 12;
    sprintf(toms_info->node_time, "%04d%03d%02d%02d00000", year, doy,
            asc_node_h, asc_node_m);
    /*
     *  get 2nd line with longitude info
     */
    if (fgets(rstr, 90, fid) == NULL) {
        printf("%s, %d E: Unexpected end of grid file found, file: %s\n",
                __FILE__, __LINE__, file);
        return 1;
    }
    if ((sscanf(&rstr[13], "%4d", &nlon) != 1) ||
            (sscanf(&rstr[35], "%7f", &slon) != 1) ||
            (sscanf(&rstr[48], "%7f", &elon) != 1) ||
            (sscanf(&rstr[60], "%4f", &del_lon) != 1)) {
        printf("%s, %d E: Trouble reading header 1 of text file: %s\n",
                __FILE__, __LINE__, file);
        return 1;
    }

    slon = -slon;
    toms_info->nlon = nlon;
    toms_info->slon = slon;
    toms_info->elon = elon;
    toms_info->del_lon = del_lon;
    /*
     *  assume W -> E but we can check what's in the file
     *  in various files of this sort, the 'W', 'E' characters may be located
     *  in duifferent columns - this handles both possibilities
     */
    if (((toupper(rstr[42]) != 'W') && (toupper(rstr[43]) != 'W')) ||
            ((toupper(rstr[55]) != 'E') && (toupper(rstr[56]) != 'E'))) {
        printf("%s, %d E: grid data is not W -> E %s\n",
                __FILE__, __LINE__, file);
        return 1;
    }
    /*
     *  get 3rd line with latitude info
     */
    if (fgets(rstr, 90, fid) == NULL) {
        printf("%s, %d E: Unexpected end of grid file found, file: %s\n",
                __FILE__, __LINE__, file);
        return 1;
    }
    if ((sscanf(&rstr[13], "%4d", &nlat) != 1) ||
            (sscanf(&rstr[35], "%7f", &slat) != 1) ||
            (sscanf(&rstr[48], "%7f", &elat) != 1) ||
            (sscanf(&rstr[60], "%4f", &del_lat) != 1)) {
        printf("%s, %d E: Trouble reading header 1 of text file: %s\n",
                __FILE__, __LINE__, file);
        return 1;
    }

    slat = -slat;
    toms_info->nlat = nlat;
    toms_info->slat = slat;
    toms_info->elat = elat;
    toms_info->del_lat = del_lat;
    /*
     *  assume S -> N but we can check what's in the file
     *  and make the start lat, lon be (-)
     */
    if (((toupper(rstr[42]) != 'S') && (toupper(rstr[43]) != 'S')) ||
            ((toupper(rstr[55]) != 'N') && (toupper(rstr[56]) != 'N'))) {
        printf("%s, %d E: grid data is not S -> N %s\n",
                __FILE__, __LINE__, file);
        return 1;
    }
    /*
     *  allocate the space for the ozone and read it in
     *  reverse array in latitude so it goes N - S and W - E
     */
    if ((toms_info->datarr = (int16 *) malloc(nlat * nlon * sizeof (int16)))
            == NULL) {
        printf("%s, %d E: Unable to allocate space for the ozone storage\n",
                __FILE__, __LINE__);
        printf("           File: %s\n", file);
        return 1;
    }
    for (il = 0; il < nlat; il++) {
        ilrev = nlat - il - 1;
        end_olin = 0;
        end_ilin = 0;
        ip = 0;
        while (end_olin == 0) {
            fgets(rstr, 90, fid);
            /* this forces any space chars to be 0 except 1st */
            /* so sscanf stays with right 3 char group */
            for (iv = 1; iv < 89; iv++)
                if (*(rstr + iv) == ' ') *(rstr + iv) = '0';
            cptr = 1;
            end_ilin = 0;
            while (end_ilin == 0) {
                ptr = &toms_info->datarr[ ip + nlon * ilrev ];
                if (sscanf((rstr + cptr), "%3hd", ptr) != 1) {
                    end_ilin = 1;
                } else {
                    cptr += 3;
                    if (ip == (nlon - 1)) {
                        end_ilin = 1;
                        end_olin = 1;
                    }
                    ip++;
                }
            }
        }
    }

    return 0;
}
