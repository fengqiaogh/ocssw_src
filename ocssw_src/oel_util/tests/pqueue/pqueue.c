
#include <pqueue.h>

#include <check.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

void print_double(void *item){
	printf("%.2f", *(double*)item);
}
bool confirm_order(pqueue *pq, int n, ...){
	va_list args;
	va_start(args, n);
	void *item;
	int i;
	for (i=0;i<n;i++){
		item = pqueue_pull(pq);
		if (item == NULL || *(double*)item != va_arg(args, double)){
			va_end(args);
			return false;
		}
	}
	va_end(args);
	return true;
}
void pqueue_push_multiple(pqueue *pq, int n, ...){
	va_list args;
	va_start(args, n);
	int i;
	for (i=0;i<n;i++){
		double priority = va_arg(args, double);
		double *item = va_arg(args, double*);
		pqueue_push(pq, priority, item);
	}
	va_end(args);
}

START_TEST(push_multiple){
	pqueue *pq = pqueue_create(0);
	double a = 9, b = 5, c = 3;
	double one = 1, three = 3, five = 5;
	pqueue_push_multiple(pq, 3, one, &a, five, &b, three, &c);
	ck_assert(confirm_order(pq, 3, a, c, b));
	pqueue_destroy(pq);
}
END_TEST

START_TEST(push_to_start){
	pqueue *pq = pqueue_create(0);
	double a = 9, b = 5, c = 3;
	double one = 1, three = 3, five = 5;
	pqueue_push_multiple(pq, 3, three, &a, five, &b, one, &c);
	ck_assert(confirm_order(pq, 3, c, a, b));
	pqueue_destroy(pq);
}
END_TEST

START_TEST(move_middle){
	pqueue *pq = pqueue_create(0);
	double a = 9, b = 5, c = 3;
	double one = 1, three = 3, five = 5;
	pqueue_push_multiple(pq, 3, one, &a, five, &b, three, &c);
	pqueue_repush(pq, 6, &c);
	ck_assert(confirm_order(pq, 3, a, b, c));
	pqueue_destroy(pq);
}
END_TEST

START_TEST(move_last){
	pqueue *pq = pqueue_create(0);
	double a = 9, b = 5, c = 3;
	double one = 1, three = 3, five = 5;
	pqueue_push_multiple(pq, 3, one, &a, five, &b, three, &c);
	pqueue_repush(pq, 6, &b);
	ck_assert(confirm_order(pq, 3, a, c, b));
	pqueue_destroy(pq);
}
END_TEST

START_TEST(move_first){
	pqueue *pq = pqueue_create(0);
	double a = 9, b = 5, c = 3;
	double one = 1, three = 3, five = 5;
	pqueue_push_multiple(pq, 3, one, &a, five, &b, three, &c);
	pqueue_repush(pq, 6, &a);
	ck_assert(confirm_order(pq, 3, c, b, a));
	pqueue_destroy(pq);
}
END_TEST

Suite* stub_suite(void){
    Suite *s = suite_create("Stub");

    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, push_multiple);
    tcase_add_test(tc_core, push_to_start);
    tcase_add_test(tc_core, move_middle);
    tcase_add_test(tc_core, move_last);
    tcase_add_test(tc_core, move_first);
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
