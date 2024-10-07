
#include <shash.h>

#include <check.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

START_TEST(set_1){
	shash *p1 = shash_create(0);
	char k1[] = "key1";
	char k2[] = "key2";
	char v1_in[] = "value1";
	char v1_2_in[] = "value1.2";
	char v2_in[] = "value2";
	char v2_const_in[] = "value3";

	ck_assert(shash_set(p1, k1, v1_in) == 0);
	const char *v1_out = shash_get(p1, k1);
	ck_assert_str_eq(v1_in, v1_out);

	ck_assert(shash_set(p1, k1, v1_2_in) == 1);
	const char *v1_2_out = shash_get(p1, k1);
	ck_assert_str_eq(v1_2_in, v1_2_out);

	ck_assert(shash_set(p1, k2, v2_in) == 0);
	const char *v2_out = shash_get(p1, k2);
	ck_assert_str_eq(v2_in, v2_out);

	ck_assert(shash_set(p1, "key2", v2_const_in) == 1);
	const char *v2_const_out = shash_get(p1, "key2");
	ck_assert_str_eq(v2_const_in, v2_const_out);
	shash_destroy(p1);
}
END_TEST

START_TEST(get_out_of_scope){
	shash *p1 = shash_create(0);
	const char *v2_const_out = shash_get(p1, "key2");
	ck_assert(v2_const_out == NULL);
	shash_destroy(p1);
}
END_TEST

START_TEST(set_out_of_scope){
	shash *p1 = shash_create(0);
	int ret;
	ret = shash_set(p1, "key2", "value3");
	ck_assert(ret == 0);
	ck_assert_str_eq(shash_get(p1, "key2"), "value3");
	shash_destroy(p1);
}
END_TEST

START_TEST(iterating){
	shash *p1 = shash_create(0);
	int key_count;
	char *keys[] = {"key1", "key2", "key3"};
	char *values[] = {"value1", "value2", "value3"};

	for (key_count = 0; key_count < 3; key_count++){
		ck_assert(shash_set(p1, keys[key_count], values[key_count]) == 0);
	}

	int replace_value = 0;
	for (replace_value = 0; replace_value < 2; replace_value++){
		const char *key;
		const char *value;
		int rewind_count = 0;
		for (rewind_count = 0; rewind_count < 3; rewind_count++){
			key = value = NULL;
			int found[key_count];
			int k;
			for (k=0;k<key_count;k++){
				found[k] = 0;
			}

//			shash_rewind(p1); //auto-rewind
			while (!shash_next(p1, &key, &value)){
				ck_assert(key != NULL && value != NULL);
				for (k=0;k<key_count;k++){
					if (!strcmp(keys[k], key) && !strcmp(values[k], value)){
						found[k]++;
					}
				}
			}

			for (k=0;k<key_count;k++){
				ck_assert_int_eq(found[k], 1);
			}
		}
		values[0] = "value4";
		ck_assert_int_eq(shash_set(p1, keys[0], values[0]), 1);
	}
	shash_destroy(p1);
}
END_TEST

START_TEST(iterating_one_value){
	shash *p1 = shash_create(0);
	int key_count;
	char *keys[] = {"key1"};
	char *values[] = {"value1"};

	for (key_count = 0; key_count < 1; key_count++){
		ck_assert(shash_set(p1, keys[key_count], values[key_count]) == 0);
	}

	int replace_value = 0;
	for (replace_value = 0; replace_value < 2; replace_value++){
		const char *key;
		const char *value;
		int rewind_count = 0;
		for (rewind_count = 0; rewind_count < 3; rewind_count++){
			key = value = NULL;
			int found[key_count];
			int k;
			for (k=0;k<key_count;k++){
				found[k] = 0;
			}

			shash_rewind(p1);
			while (!shash_next(p1, &key, &value)){
				ck_assert(key != NULL && value != NULL);
				for (k=0;k<key_count;k++){
					if (!strcmp(keys[k], key) && !strcmp(values[k], value)){
						found[k]++;
					}
				}
			}

			for (k=0;k<key_count;k++){
				ck_assert_int_eq(found[k], 1);
			}
		}
		values[0] = "value4";
		ck_assert_int_eq(shash_set(p1, keys[0], values[0]), 1);
	}
	shash_destroy(p1);
}
END_TEST

START_TEST(remove_1){
	shash *p1 = shash_create(0);
	int key_count;
	char *keys[] = {"key1", "key2", "key3"};
	char *values[] = {"value1", "value2", "value3"};

	for (key_count = 0; key_count < 3; key_count++){
		ck_assert_int_eq(shash_set(p1, keys[key_count], values[key_count]), 0);
	}

	int remove_value = 0;
	for (remove_value = 0; remove_value < key_count-1; remove_value++){
		const char *key;
		const char *value;
		key = value = NULL;
		int found[key_count];
		int k;
		for (k=0;k<key_count;k++){
			found[k] = 0;
		}

		shash_rewind(p1);
		while (!shash_next(p1, &key, &value)){
			ck_assert(key != NULL && value != NULL);
			for (k=0;k<key_count;k++){
				if (!strcmp(keys[k], key) && !strcmp(values[k], value)){
					found[k]++;
				}
			}
		}
		k = remove_value;
		for (;k<key_count;k++){
			ck_assert_int_eq(found[k], 1);
		}
		ck_assert_int_eq(shash_remove(p1, keys[remove_value]), 0);
	}
	ck_assert_ptr_eq(shash_get(p1, keys[remove_value-1]), NULL);
	ck_assert_ptr_ne(shash_get(p1, keys[remove_value]), NULL);
	shash_destroy(p1);
}
END_TEST

static const char characters[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
static void random_string(char *buffer, int length){
	int i=0;
	for (i=0;i<length;i++){
		buffer[i] = characters[rand() % (sizeof(characters)-1)];
	}
}

START_TEST(stress_test){
	shash *p1 = shash_create(0);
	int i, j, contains_dupes, count = 1000, length = 8;
	char keys[count][length+1];
	char values[count][length+1];

	for (i=0;i<count;i++){
		random_string(keys[i], length);
		random_string(values[i], length);
		keys[i][length] = '\0';
		values[i][length] = '\0';
	}
	do {
		contains_dupes = 0;
		for (i=0;i<count;i++){
			for (j=i+1;j<count;j++){
				if (!strcmp(keys[i], keys[j])){
					contains_dupes = 1;
					random_string(keys[j], length);
				}
			}
		}
	} while (contains_dupes);

	for (i=0;i<count;i++){
		ck_assert_int_eq(shash_set(p1, keys[i], values[i]), 0);
	}

	for (i=0;i<count;i++){
		ck_assert_str_eq(values[i], shash_get(p1, keys[i]));
	}
	shash_destroy(p1);
}
END_TEST

Suite* stub_suite(void){
    Suite *s = suite_create("Stub");

    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, set_1);
    tcase_add_test(tc_core, remove_1);
    tcase_add_test(tc_core, iterating);
    tcase_add_test(tc_core, iterating_one_value);
    tcase_add_test(tc_core, stress_test);
    suite_add_tcase(s, tc_core);

    TCase *tc_segfaults = tcase_create("Segfaults");
    tcase_add_test(tc_segfaults, set_1);
    tcase_add_test(tc_segfaults, get_out_of_scope);
    tcase_add_test(tc_segfaults, set_out_of_scope);
    suite_add_tcase(s, tc_segfaults);

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
