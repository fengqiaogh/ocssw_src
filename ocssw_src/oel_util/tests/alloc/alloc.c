
#include <allocate2d.h>

#include <check.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

START_TEST(allocate2d){
	size_t w = 5, h = 10;
	// NB the order of the for loops
	char **t1 = allocate2d_char(h, w);
	for (size_t i=0;i<h;i++){
		for (size_t j=0;j<w;j++){
			t1[i][j] = (i * w) + j;
		}
	}
	for (size_t i=0;i<h;i++){
		for (size_t j=0;j<w;j++){
			ck_assert_int_eq(t1[i][j], (i * w) + j);
		}
	}
	free2d_char(t1);
}
END_TEST

Suite* stub_suite(void){
    Suite *s = suite_create("Stub");

    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, allocate2d);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(int argc, char **argv){
    int number_failed;

    Suite *s = stub_suite();
    SRunner *sr = srunner_create(s);

    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

