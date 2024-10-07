
#include <vincenty.h>

#include <check.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>


// https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
bool dbl_equal_eps(double A, double B, double maxRelDiff){
    // Calculate the difference.
	double diff = fabs(A - B);
    A = fabs(A);
    B = fabs(B);
    // Find the largest
    double largest = (B > A) ? B : A;

    if (diff <= largest * maxRelDiff)
        return true;
    return false;
}
bool dbl_equal(double A, double B){
    return dbl_equal_eps(A, B, DBL_EPSILON);
}

START_TEST(edge_cases){
	ck_assert(dbl_equal(vincenty_distance(0.000000, 0.000000, 0.000000, 0.000000), 0.000000));
	ck_assert(dbl_equal_eps(vincenty_distance(0.000000, 0.000000, 0.000000, -2.000000), 222639, 0.0000001));
	ck_assert(dbl_equal_eps(vincenty_distance(0.000000, 0.000000, 0.000000, -1.000000), 111319.5, 0.0000001));
	ck_assert(dbl_equal_eps(vincenty_distance(0.000000, 0.000000, 0.000000, 1.000000), 111319.5, 0.0000001));
	ck_assert(dbl_equal_eps(vincenty_distance(0.000000, 0.000000, 0.000000, 2.000000), 222639, 0.0000001));
}
END_TEST

Suite* stub_suite(void){
    Suite *s = suite_create("Stub");

    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, edge_cases);
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
