
#include "filetype.h"

#include <check.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <sensorDefs.h>
#include <sensorInfo.h>

#if 0
// All types to check, delete from here when implemented:
FT_L1BNCDF,
FT_L1HDF,
FT_L1XCAL,
FT_L2HDF,
FT_L2NCDF,
FT_L3BIN,
FT_L3MAP,

FT_AVIRIS,
FT_CLASSAVHRR,

FT_OCML1B,
FT_OCML1BDB,


FT_OLCI,
FT_OLIL1B,
FT_OCIA,
FT_OSMIL1A,
FT_PRISM,
FT_SGLI, // SGLI files (hdf5)

FT_L5TML1B,
FT_L7ETML1B,

FT_MSIL1C /* Sentinel MSI L1C */


//These weren't found:

FT_MODISL1B, // MODIS L1B (hdf4, ocean-color band subset)
FT_VIIRSGEO, // VIIRS Geolocation (hdf5)

FT_OCTSL1B,

FT_MERISCC,
FT_MERISL2,

#endif

static void ck_assert_file_format(file_format format, file_type type, int32_t sensor_id, int32_t subsensor_id) {
    ck_assert_int_eq(format.type, type);
    ck_assert_int_eq(format.sensor_id, sensor_id);
    ck_assert_int_eq(format.subsensor_id, subsensor_id);
}






START_TEST(octs) {
    ck_assert_file_format(getFormat("../l2gen/octs/O1997079195357.L1A_GAC.subline"), FT_OCTSL1A, OCTS, -1);
    ck_assert_file_format(getFormat("../l2gen/octs/O1997079195357.L2_GAC_IOP.subline.nc"), FT_L2NCDF, OCTS, -1);
}
END_TEST

START_TEST(ocm2) {
    ck_assert_file_format(getFormat("../l2gen/ocm2/O2_01FEB2013_026_013_GAN_L1B_ST_S.hdf"), FT_OCM2L1B, OCM2, -1);
    ck_assert_file_format(getFormat("../l2gen/ocm2/O2_2012032171306.L2_OC.nc"), FT_L2NCDF, OCM2, -1);
}
END_TEST

START_TEST(sundry_l1s) {
    ck_assert_file_format(getFormat("../l2gen/czcs/C1982103032950.L1A_MLAC.subline"), FT_CZCSL1A, CZCS, -1);
    ck_assert_file_format(getFormat("../l2gen/hico/H2014099041255.L1B_ISS"), FT_HICOL1B, HICO, -1);
    ck_assert_file_format(getFormat("../l2gen/mos/M1998059110129.L1B_HWFF.nc"), FT_L1BNCDF, MOS, -1);
    ck_assert_file_format(getFormat("../l2gen/mos/M1998059110129.L1B_HWFF_00023B23"), FT_MOSL1B, MOS, -1);
    ck_assert_file_format(getFormat("../l2gen/seawifs/S1998199173926.L1A_GAC.subline"), FT_SEAWIFSL1A, SEAWIFS, SEAWIFS_GAC);
    ck_assert_file_format(getFormat("../l2gen/seawifs/S2002079071209.L1A_MLAC.subline"), FT_SEAWIFSL1A, SEAWIFS, SEAWIFS_LAC);
}
END_TEST

START_TEST(meris) {
    ck_assert_file_format(getFormat("../l2gen/meris/M2008261123930.L2_LAC_IOP.subpix.nc"), FT_L2NCDF, MERIS, -1);
    ck_assert_file_format(getFormat("../l2gen/meris/M2008261123930.L2_LAC_OC.subline.nc"), FT_L2NCDF, MERIS, -1);
    ck_assert_file_format(getFormat("../l2gen/meris/MER_RR__1PRACR20080917_123930_000026292072_00124_34246_0000.N1.subline"), FT_MERISL1B, MERIS, -1);
    ck_assert_file_format(getFormat("../l2gen/meris/MER_RR__1PRACR20080917_123930_000026292072_00124_34246_0000.N1.subpix"), FT_MERISL1B, MERIS, -1);
}
END_TEST


START_TEST(aqua) {
    ck_assert_file_format(getFormat("../l2gen/aqua/A2008080124500.L1B_LAC.subline"), FT_MODISL1B, MODISA, MODIS_AQUA);
    ck_assert_file_format(getFormat("../l2gen/aqua/A2008080124500.GEO.subline"), FT_MODISGEO, MODISA, MODIS_AQUA);
    ck_assert_file_format(getFormat("../l2gen/aqua/A2008080124500.L2_LAC_SST.subline.nc"), FT_L2NCDF, MODISA, MODIS_AQUA);
    ck_assert_file_format(getFormat("../l2gen/aqua/A2008080195500.L2_LAC_OC.subline.nc"), FT_L2NCDF, MODISA, MODIS_AQUA);
    ck_assert_file_format(getFormat("../l2gen/aqua/A2008080195500.L2_LAC_OC.subline_swir.nc"), FT_L2NCDF, MODISA, MODIS_AQUA);
    ck_assert_file_format(getFormat("../l2gen/aqua/A2008080214500.L2_LAC_IOP.subline.nc"), FT_L2NCDF, MODISA, MODIS_AQUA);
    ck_assert_file_format(getFormat("../l2gen/aqua/A2008080214500.L2_LAC_SST.subline.nc"), FT_L2NCDF, MODISA, MODIS_AQUA);
    ck_assert_file_format(getFormat("../l2gen/aqua/A2008080220000.L2_LAC_OC.subpix_250m.nc"), FT_L2NCDF, MODISA, MODIS_AQUA);
    ck_assert_file_format(getFormat("../l2gen/aqua/A2008080220000.L2_LAC_OC.subpix_500m.nc"), FT_L2NCDF, MODISA, MODIS_AQUA);
    ck_assert_file_format(getFormat("../l2gen/aqua/A2008080220000.L2_LAC_OC.subpix.nc"), FT_L2NCDF, MODISA, MODIS_AQUA);
}
END_TEST

START_TEST(terra) {
    ck_assert_file_format(getFormat("../l2gen/terra/T2008080095000.L1B_LAC.subline"), FT_MODISL1B, MODIST, MODIS_TERRA);
    ck_assert_file_format(getFormat("../l2gen/terra/T2008080122500.GEO.subline"), FT_MODISGEO, MODIST, MODIS_TERRA);
    ck_assert_file_format(getFormat("../l2gen/terra/T2008080095000.L2_LAC_SST4.subline.nc"), FT_L2NCDF, MODIST, MODIS_TERRA);
    ck_assert_file_format(getFormat("../l2gen/terra/T2008080122500.L2_LAC_OC.subline.nc"), FT_L2NCDF, MODIST, MODIS_TERRA);
    ck_assert_file_format(getFormat("../l2gen/terra/T2008080122500.L2_LAC_SST.subline.nc"), FT_L2NCDF, MODIST, MODIS_TERRA);
}
END_TEST

START_TEST(viirs_npp) {
    ck_assert_file_format(getFormat("../l2gen/viirs/SVM01_npp_d20140320_t2221504_e2223145_b12407_obpg_ops.h5"), FT_VIIRSL1B, VIIRSN, VIIRS_NPP);
    ck_assert_file_format(getFormat("../l2gen/viirs/V2014079222148.L2_NPP_OC.nc"), FT_L2NCDF, VIIRSN, VIIRS_NPP);
    ck_assert_file_format(getFormat("../l2gen/viirs/V2017226170000.L1A_SNPP.GEO_extract.nc"), FT_VIIRSGEONC, VIIRSN, VIIRS_NPP);
    ck_assert_file_format(getFormat("../l2gen/viirs/V2017226170000.L1A_SNPP_extract.nc"), FT_VIIRSL1A, VIIRSN, VIIRS_NPP);
    ck_assert_file_format(getFormat("../l2gen/viirs/V2017226170000.L1B_SNPP_extract.nc"), FT_VIIRSL1BNC, VIIRSN, VIIRS_NPP);
}
END_TEST

Suite *stub_suite(void) {
    Suite *s = suite_create("Stub");

    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, aqua);
    tcase_add_test(tc_core, sundry_l1s);
    tcase_add_test(tc_core, terra);
    tcase_add_test(tc_core, viirs_npp);
    tcase_add_test(tc_core, octs);
    tcase_add_test(tc_core, meris);
    tcase_add_test(tc_core, ocm2);
    suite_add_tcase(s, tc_core);

    return s;
}

int main() {
    int number_failed;

    Suite *s = stub_suite();
    SRunner *sr = srunner_create(s);

    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
