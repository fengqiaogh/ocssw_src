
#include <genutils.h>

#include <check.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char* ocvarroot = "/base/ocssw/var";
const char* ocdataroot = "/base/ocssw/share";
const char* ocsswroot = "/base/ocssw";

START_TEST(ocvar){
    char str[50];
    strcpy(str, ocvarroot);
    strcat(str, "/cool/dir");
    char* str2 = replace_ocroots(str);
    ck_assert_str_eq(str2, "$OCVARROOT/cool/dir");
    free(str2);
}
END_TEST

START_TEST(ocdata){
    char str[50];
    strcpy(str, ocdataroot);
    strcat(str, "/cool/dir");
    char* str2 = replace_ocroots(str);
    ck_assert_str_eq(str2, "$OCDATAROOT/cool/dir");
    free(str2);
}
END_TEST

START_TEST(ocssw){
    char str[50];
    strcpy(str, ocsswroot);
    strcat(str, "/cool/dir");
    char* str2 = replace_ocroots(str);
    ck_assert_str_eq(str2, "$OCSSWROOT/cool/dir");
    free(str2);
}
END_TEST

START_TEST(many){
    char str[500];
    
    strcpy(str, ocsswroot);
    strcat(str, "/cool/dir1;");
    strcat(str, ocsswroot);
    strcat(str, "/cool/dir2;");

    strcat(str, ocvarroot);
    strcat(str, "/cool/dir3;");
    strcat(str, ocvarroot);
    strcat(str, "/cool/dir4;");

    strcat(str, ocdataroot);
    strcat(str, "/cool/dir5;");
    strcat(str, ocdataroot);
    strcat(str, "/cool/dir6");

    char* str2 = replace_ocroots(str);
    ck_assert_str_eq(str2, "$OCSSWROOT/cool/dir1;$OCSSWROOT/cool/dir2;$OCVARROOT/cool/dir3;$OCVARROOT/cool/dir4;$OCDATAROOT/cool/dir5;$OCDATAROOT/cool/dir6");
    free(str2);
}
END_TEST

Suite* stub_suite(void){
    Suite *s = suite_create("Stub");

    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, ocvar);
    tcase_add_test(tc_core, ocdata);
    tcase_add_test(tc_core, ocssw);
    tcase_add_test(tc_core, many);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(int argc, char **argv){
    int number_failed;

    // modify the environment variables
    setenv("OCVARROOT", ocvarroot, 1);
    setenv("OCDATAROOT", ocdataroot, 1);
    setenv("OCSSWROOT", ocsswroot, 1);

    Suite *s = stub_suite();
    SRunner *sr = srunner_create(s);

    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

