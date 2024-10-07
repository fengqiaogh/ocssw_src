
#include <argpar.h>

#include <check.h>
#include <errno.h>
#include <ftw.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

FILE *argpar_ostream;
const char *argpar_program_name = "argpar-test";
static const char doc[] = "This is an argpar test.";
static const char args_doc[] = "[none]";
static const argpar_option options[] = {
	{ "ifile", 'f', "FILE", 0, "input file" },
	{ "ofile", 'o', "FILE", 0, "output file" },
	{ 0,0,0,0, "This is a header:", 2 },
	{ 0,0,0, OPTION_DOC, "This is some docs." },
	{ 0,0,0,0, "Header #1:", 1 },
	{ 0,0,0, OPTION_DOC, "This is some docs." },
	{ "efile", 'e', "FILE", 0, "error file", 2 },
	{ "long", 'l', 0, OPTION_INT, "long, not casted" },
	{ "int", 'i', 0, OPTION_INT, "casted as int", 1 },
    { "dbl", 'd', 0, OPTION_DBL, "casted as dbl", 1 },
		{ "double", 0, 0, OPTION_ALIAS },
    { "p", 'p', 0, 0, "load parfile via child", 1 },
	{ "enum", 'm', 0, 0, "this option has enum values", 1 },
		{ "val1", 0, 0, OPTION_ENUM, "description of val1", 1 },
		{ "val2", 0, 0, OPTION_ENUM, "description of val2", 1 },
		{ "default", 0, 0, OPTION_ATTR, "val1" },
	{0}
};

char stdout_file[64];

typedef struct arguments {
	const char *argv[32];
	int argc;
	const char *ifile;
	const char *ofile;
	char *ofile_src;
	const char *efile;
	char *efile_src;
	int i; bool i_error;
	double d; bool d_error;
	char *d_alias;
	long l; bool l_error;

	bool ended;
	const char *argv_unknown[32];
	int argc_unknown;

	int unknown;

} arguments;


static void clean_arguments(arguments *a){
    if (a->ofile_src){
        free(a->ofile_src);
    }
    if (a->efile_src){
        free(a->efile_src);
    }
    if (a->d_alias){
        free(a->d_alias);
    }
}

static int parse_options(int key, char *argv, argpar_state *state) {
    arguments *arguments = state->input;
	switch (key){
	case 'f':
		arguments->ifile = argv;
		break;
	case 'o':
		arguments->ofile = argv;
        if (arguments->ofile_src){
            free(arguments->ofile_src);
            arguments->ofile_src = NULL;
        }
        if (state->parfile){
            arguments->ofile_src = strdup(state->parfile);
        }
		break;
	case 'e':
		arguments->efile = state->arg_value;
        if (arguments->efile_src){
            free(arguments->efile_src);
            arguments->efile_src = NULL;
        }
		if (state->parfile){
		    arguments->efile_src = strdup(state->parfile);
		}
		break;
	case 'i':
		arguments->i = state->argv_as_int;
		arguments->i_error = state->argv_as_int_err;
		break;
	case 'd':
		arguments->d = state->argv_as_dbl;
		arguments->d_error = state->argv_as_dbl_err;
        if (arguments->d_alias){
            free(arguments->d_alias);
            arguments->d_alias = NULL;
        }
        if (state->arg_alias){
            arguments->d_alias = strdup(state->arg_alias);
        }
		break;
	case 'l':
		arguments->l = state->argv_as_int;
		arguments->l_error = state->argv_as_int_err;
		break;
    case 'p':
        ck_assert_int_eq(argpar_parse_file(state->argpar, state->arg_value, 0, arguments), 0);
        break;
	case ARGPAR_KEY_ARG:
		arguments->argv[arguments->argc++] = argv;
		break;
	case ARGPAR_KEY_UNKNOWN:
		arguments->unknown++;
		break;
	case ARGPAR_KEY_END:
		arguments->ended = true;
		break;
	}

	return 0;
}

argpar params = { options, parse_options, args_doc, doc };

START_TEST(no_args_usage){
    argpar_ostream = tmpfile();

	arguments arguments = {};
	ck_assert_int_eq(argpar_parse_args(&params, 0, NULL, 0, NULL, &arguments), ARGPAR_ERR_USAGE);

    argpar_clean(&params);
    fclose(argpar_ostream);
}
END_TEST

START_TEST(no_real_args){
    argpar_ostream = tmpfile();

    int argc = 2;
    char *argv[] = {(char*)argpar_program_name, "--"};
    arguments arguments = {
        .argc = 0, .argc_unknown = 0, .i = -1, .d = -1
    };
    ck_assert_int_eq(argpar_parse_args(&params, argc, argv, 0, NULL, &arguments), 0);
    ck_assert_int_eq(arguments.argc, 0);
    ck_assert_int_eq(arguments.argc_unknown, 0);
    ck_assert_int_eq(arguments.i, -1);
    ck_assert(arguments.d == -1);

    clean_arguments(&arguments);

    argpar_clean(&params);
    fclose(argpar_ostream);
}
END_TEST

START_TEST(auto_help){
    argpar_ostream = tmpfile();
    int argc = 2;
    char *argv_const[] = {(char*)argpar_program_name, "help=1"};
	char *argv[argc];
	int i;
	for (i=0;i<argc;i++){
		argv[i] = strdup(argv_const[i]);
	}
    arguments arguments = {.argc=0};

    ck_assert_int_eq(argpar_parse_args(&params, argc, argv, ARGPAR_NO_EXIT, NULL, &arguments), ARGPAR_ERR_ABORT);

	for (i=0;i<argc;i++){
		free(argv[i]);
	}

    clean_arguments(&arguments);

    argpar_clean(&params);
    fclose(argpar_ostream);
}
END_TEST

START_TEST(unknown_option){
    argpar_ostream = tmpfile();

	int argc = 7;
	char *argv_const[] = {(char*)argpar_program_name, "ifile=input", "arg1", "unknown=bad", "--", "ofile=output", "arg3"};
	char *argv[argc];
	int i;
	for (i=0;i<argc;i++){
		argv[i] = strdup(argv_const[i]);
	}
	arguments arguments = {
		.argc = 0, .i = -1, .d = -1,
		.ofile = NULL, .ended = false
	};
	ck_assert_int_eq(argpar_parse_args(&params, argc, argv, 0, NULL, &arguments), ARGPAR_ERR_UNKNOWN);

	for (i=0;i<argc;i++){
		free(argv[i]);
	}

    clean_arguments(&arguments);

    argpar_clean(&params);
    fclose(argpar_ostream);
}
END_TEST

START_TEST(no_key_args){
    argpar_ostream = tmpfile();

	int argc = 7;
	const char *argv_const[] = {
		argpar_program_name, "ifile=_ifile", "ofile=_ofile", "efile=_efile",
		"int=1", "dbl=2", "long=3"
	};
	char *argv[argc];
	int i;
	for (i=0;i<argc;i++){
		argv[i] = strdup(argv_const[i]);
	}
	arguments arguments = {
		.argc = 0, .argc_unknown = 0, .i = -1, .d = -1, .l = -1
	};
	ck_assert_int_eq(argpar_parse_args(&params, argc, argv, 0, NULL, &arguments), 0);
	ck_assert_int_eq(arguments.argc, 0);
	ck_assert_int_eq(arguments.argc_unknown, 0);
	ck_assert_int_eq(arguments.i, 1);
	ck_assert_int_eq(arguments.i_error, 0);
	ck_assert_int_eq(arguments.d, 2);
	ck_assert_int_eq(arguments.d_error, 0);
	ck_assert_int_eq(arguments.l, 3);
	ck_assert_int_eq(arguments.l_error, 0);

	for (i=0;i<argc;i++){
		free(argv[i]);
	}

    clean_arguments(&arguments);

    argpar_clean(&params);
    fclose(argpar_ostream);
}
END_TEST

START_TEST(aliases){
    argpar_ostream = tmpfile();

	int argc = 3;
	const char *argv_const[] = {
		argpar_program_name, "dbl=1", "double=2"
	};
	char *argv[argc];
	int i;
	for (i=0;i<argc;i++){
		argv[i] = strdup(argv_const[i]);
	}
	arguments arguments = {
		.argc = 0, .argc_unknown = 0, .i = -1, .d = -1, .l = -1
	};
	ck_assert_int_eq(argpar_parse_args(&params, argc, argv, 0, NULL, &arguments), 0);
	ck_assert_int_eq(arguments.d, 2);
	ck_assert_int_eq(arguments.d_error, 0);
	ck_assert_str_eq(arguments.d_alias, "double");

	for (i=0;i<argc;i++){
		free(argv[i]);
	}

    clean_arguments(&arguments);

    argpar_clean(&params);
    fclose(argpar_ostream);
}
END_TEST

START_TEST(int_errors){
    argpar_ostream = tmpfile();

	int tests = 8;
	const char *int_strs[] = { "1", "2.0", "456", "-123", "0", "-0", "q", "345a" };
	const int ints[]       = {  1,   2,     456,   -123,   0,    0,   0,   345 };
	const bool errors[]    = {  0,   1,     0,     0,      0,    0,   1,   1 };
	char *argv[2] = {(char*)argpar_program_name, malloc(0)};
	int i;
	for (i=0;i<tests;i++){
		argv[1] = realloc(argv[1], sizeof(char)*(strlen(int_strs[i])+5));
		strcpy(argv[1], "int=");
		strcat(argv[1], int_strs[i]);
		arguments arguments = { .i = -1 };

		ck_assert_int_eq(argpar_parse_args(&params, 2, argv, 0, NULL, &arguments), 0);
		ck_assert_int_eq(arguments.i, ints[i]);
		ck_assert_int_eq(arguments.i_error, errors[i]);

	    clean_arguments(&arguments);
	}
	free(argv[1]);

    argpar_clean(&params);
    fclose(argpar_ostream);
}
END_TEST

START_TEST(dbl_errors){
    argpar_ostream = tmpfile();

	int tests = 9;
	const char *dbl_strs[] = { "1", "2.0", "456", "-123", "0", "-0", "q", "345a", "-2e5" };
	const double dbls[]    = {  1,   2,     456,   -123,   0,    0,   0,   345,   -200000 };
	const bool errors[]    = {  0,   0,     0,     0,      0,    0,   1,   1,     0 };
	char *argv[2] = {(char*)argpar_program_name, malloc(0)};
	int i;
	for (i=0;i<tests;i++){
		argv[1] = realloc(argv[1], sizeof(char)*(strlen(dbl_strs[i])+5));
		strcpy(argv[1], "dbl=");
		strcat(argv[1], dbl_strs[i]);
		arguments arguments = { .i = -1 };

		ck_assert_int_eq(argpar_parse_args(&params, 2, argv, 0, NULL, &arguments), 0);
		ck_assert(arguments.d == dbls[i]);
		ck_assert_int_eq(arguments.d_error, errors[i]);

	    clean_arguments(&arguments);
	}
	free(argv[1]);

    argpar_clean(&params);
    fclose(argpar_ostream);
}
END_TEST


START_TEST(parfile1){
    argpar_ostream = tmpfile();

	int argc = 5;
	const char *argv_const[] = {
			argpar_program_name, "parfile=parfiles/initial.par", "ofile=keyarg_ofile1", "parfile=parfiles/final.par",
			"ofile=keyarg_ofile2"
	};
	char *argv[argc];
	int i;
	for (i=0;i<argc;i++){
		argv[i] = strdup(argv_const[i]);
	}
	arguments arguments = {
		.argc = 0, .argc_unknown = 0, .i = -1, .d = -1, .l = -1
	};
	ck_assert_int_eq(argpar_parse_args(&params, argc, argv, 0, NULL, &arguments), 0);
	ck_assert_int_eq(arguments.argc, 0);
	ck_assert_int_eq(arguments.argc_unknown, 0);
	ck_assert_int_eq(arguments.i, -1);
	ck_assert_int_eq(arguments.i_error, 0);
	ck_assert_int_eq(arguments.d, 2.2);
	ck_assert_int_eq(arguments.d_error, 0);
	ck_assert_int_eq(arguments.l, -11);
	ck_assert_int_eq(arguments.l_error, 0);
	ck_assert_str_eq(arguments.ifile, "initial_input");
	ck_assert_str_eq(arguments.ofile, "keyarg_ofile2");
    ck_assert_ptr_eq(arguments.ofile_src, NULL);
    ck_assert_str_eq(arguments.efile, "defaults_error");
    ck_assert_str_eq(arguments.efile_src, "parfiles/defaults.par");

	for (i=0;i<argc;i++){
		free(argv[i]);
	}

    clean_arguments(&arguments);

    argpar_clean(&params);
    fclose(argpar_ostream);
}
END_TEST

static const char doc_parent[] = "This is an argpar test.";
static const char args_doc_parent[] = "[none]";
static const argpar_option options_parent[] = {
	{ "dbl", 'd', 0, OPTION_DBL, "casted as dbl", 1 },
	{ "parg", 'p', 0, OPTION_INT, "new arg for parent", 1 },
	{0}
};

typedef struct arguments_parent {
	arguments _;
	int parent_arg;
	int argc;
} arguments_parent;

static int parse_options_parent(int key, char *argv, argpar_state *state) {
	arguments_parent *arguments = state->input;
	switch (key){
	case 'd':
		arguments->_.d = state->argv_as_dbl * 2;
		break;
	case 'p':
		arguments->parent_arg = state->argv_as_int;
		break;
	case ARGPAR_KEY_ARG:
		if (arguments->argc){
			return ARGPAR_ERR_UNKNOWN;
		}
		arguments->argc++;
		break;
	case ARGPAR_KEY_INIT:
		state->child_inputs[0] = state->input;
		break;
	}
	return 0;
}
argpar_child children[] = {{&params}, {0}};
argpar params_parent = { options_parent, parse_options_parent, args_doc_parent, doc_parent, children };

START_TEST(children1){
    argpar_ostream = tmpfile();

    int argc = 6;
    char *argv_const[] = {(char*)argpar_program_name, "dbl=2", "parg=5", "int=6", "arg1", "arg2"};
    char *argv[argc];
    int i;
    for (i=0;i<argc;i++){
        argv[i] = strdup(argv_const[i]);
    }
    arguments arguments = {
        .argc = 0, .i = -1, .d = -1,
        .ofile = NULL, .ended = false
    };
    arguments_parent arguments_parent = {arguments, 0, 0};
    ck_assert_int_eq(argpar_parse_args(&params_parent, argc, argv, 0, NULL, &arguments_parent), 0);
    ck_assert(arguments_parent._.d == 4);
    ck_assert_int_eq(arguments_parent._.i, 6);
    ck_assert_int_eq(arguments_parent.parent_arg, 5);
    ck_assert_int_eq(arguments_parent._.argc, 1);
    ck_assert_int_eq(arguments_parent.argc, 1);

    for (i=0;i<argc;i++){
        free(argv[i]);
    }

    clean_arguments(&arguments);

    argpar_clean(&params_parent);
    fclose(argpar_ostream);
}
END_TEST

START_TEST(parse_parfile){
    argpar_ostream = tmpfile();

    int argc = 2;
    char *argv_const[] = {(char*)argpar_program_name, "dbl=2"};
    char *argv[argc];
    int i;
    for (i=0;i<argc;i++){
        argv[i] = strdup(argv_const[i]);
    }
    arguments arguments = {
        .argc = 0, .i = -1, .d = -1, .l = 4,
        .ofile = NULL, .ended = false
    };
    arguments_parent arguments_parent = {arguments, 0, 0};
    ck_assert_int_eq(argpar_parse_args(&params_parent, argc, argv, 0, NULL, &arguments_parent), 0);
    ck_assert(arguments_parent._.d == 4);

    ck_assert_int_eq(argpar_parse_file(&params_parent, "parfiles/initial.par", 0, &arguments_parent), 0);
    ck_assert(arguments_parent._.d == 2.2);
    ck_assert_str_eq(arguments_parent._.efile, "defaults_error");
    ck_assert_str_eq(arguments_parent._.ofile_src, "parfiles/defaults.par");
    ck_assert_str_eq(arguments_parent._.efile_src, "parfiles/defaults.par");


    ck_assert_int_eq(arguments_parent._.l, -11);

    for (i=0;i<argc;i++){
        free(argv[i]);
    }

    clean_arguments(&arguments_parent._);

    argpar_clean(&params_parent);
    fclose(argpar_ostream);
}
END_TEST

START_TEST(parse_parfile_in_child){
    argpar_ostream = tmpfile();

    int argc = 2;
    char *argv_const[] = {(char*)argpar_program_name, "p=parfiles/initial.par"};
    char *argv[argc];
    int i;
    for (i=0;i<argc;i++){
        argv[i] = strdup(argv_const[i]);
    }
    arguments arguments = {
        .argc = 0, .i = -1, .d = -1, .l = 4,
        .ofile = NULL, .ended = false
    };

    arguments_parent arguments_parent = {arguments, 0, 0};
    ck_assert_int_eq(argpar_parse_args(&params_parent, argc, argv, 0, NULL, &arguments_parent), 0);

    ck_assert_int_eq(arguments_parent._.d, 2.2);
    ck_assert_str_eq(arguments_parent._.efile, "defaults_error");
    ck_assert_int_eq(arguments_parent._.l, -11);

    for (i=0;i<argc;i++){
        free(argv[i]);
    }

    clean_arguments(&arguments_parent._);

    argpar_clean(&params_parent);
    fclose(argpar_ostream);
}
END_TEST

START_TEST(skip_parfiles){
    argpar_ostream = tmpfile();

    int argc = 3;
    char *argv_const[] = {(char*)argpar_program_name, "ofile=keyarg_ofile1",  "parfile=parfiles/initial.par"};
    char *argv[argc];
    int i;
    for (i=0;i<argc;i++){
        argv[i] = strdup(argv_const[i]);
    }
    arguments arguments = {
        .argc = 0, .argc_unknown = 0, .i = -1, .d = -1
    };
    ck_assert_int_eq(argpar_parse_args(&params, argc, argv, ARGPAR_SKIP_PARFILES, NULL, &arguments), 0);
    ck_assert_str_eq(arguments.ofile, "keyarg_ofile1");

    for (i=0;i<argc;i++){
        free(argv[i]);
    }

    clean_arguments(&arguments);

    argpar_clean(&params);
    fclose(argpar_ostream);
}
END_TEST

START_TEST(no_keyargs){
    argpar_ostream = tmpfile();

    int argc = 2;
    char *argv_const[] = {(char*)argpar_program_name, "int"};
    char *argv[argc];
    int i;
    for (i=0;i<argc;i++){
        argv[i] = strdup(argv_const[i]);
    }
    arguments arguments = {
        .argc = 0, .argc_unknown = 0, .i = -1, .d = -1
    };
    ck_assert_int_eq(argpar_parse_args(&params, argc, argv, ARGPAR_NO_KEYARGS, NULL, &arguments), 0);
    ck_assert_int_eq(arguments.i, 1);
    ck_assert_int_eq(arguments.unknown, 0);

    for (i=0;i<argc;i++){
        free(argv[i]);
    }

    clean_arguments(&arguments);

    argpar_clean(&params);
    fclose(argpar_ostream);
}
END_TEST

START_TEST(no_keyargs_parfile){
    argpar_ostream = tmpfile();

    int argc = 2;
    char *argv_const[] = {(char*)argpar_program_name, "parfile=parfiles/no_keyargs.par"};
    char *argv[argc];
    int i;
    for (i=0;i<argc;i++){
        argv[i] = strdup(argv_const[i]);
    }
    arguments arguments = {
        .argc = 0, .argc_unknown = 0, .i = -1, .d = -1
    };
    ck_assert_int_eq(argpar_parse_args(&params, argc, argv, ARGPAR_NO_KEYARGS, NULL, &arguments), 0);
    ck_assert_int_eq(arguments.i, 1);

    for (i=0;i<argc;i++){
        free(argv[i]);
    }

    clean_arguments(&arguments);

    argpar_clean(&params);
    fclose(argpar_ostream);
}
END_TEST

START_TEST(split_strings){
    char *str1, *str2, *str3, *str4, *str5, *str6, **ret1;

    str1 = strdup("test1,test2,test3");
    str2 = strdup("test1 test2 test3");
    str3 = strdup("test1, test2, test3");
    str4 = strdup("test1");
    str5 = strdup("test1,,test2");
    str6 = strdup("test1,test2,test3,test4,test5,test6,test7,test8,test9,test10");

    ret1 = argpar_split_str(str1, ",");
    ck_assert_str_eq(ret1[0], "test1");
    ck_assert_str_eq(ret1[1], "test2");
    ck_assert_str_eq(ret1[2], "test3");
    ck_assert_ptr_eq(ret1[3], NULL);
    free(ret1);

    ret1 = argpar_split_str(str2, " ");
    ck_assert_str_eq(ret1[0], "test1");
    ck_assert_str_eq(ret1[1], "test2");
    ck_assert_str_eq(ret1[2], "test3");
    ck_assert_ptr_eq(ret1[3], NULL);
    free(ret1);

    ret1 = argpar_split_str(str3, ", ");
    ck_assert_str_eq(ret1[0], "test1");
    ck_assert_str_eq(ret1[1], "");
    ck_assert_str_eq(ret1[2], "test2");
    ck_assert_str_eq(ret1[3], "");
    ck_assert_str_eq(ret1[4], "test3");
    ck_assert_ptr_eq(ret1[5], NULL);
    free(ret1);

    ret1 = argpar_split_str(str4, ",");
    ck_assert_str_eq(ret1[0], "test1");
    ck_assert_ptr_eq(ret1[1], NULL);
    free(ret1);

    ret1 = argpar_split_str(str5, ",");
    ck_assert_str_eq(ret1[0], "test1");
    ck_assert_str_eq(ret1[1], "");
    ck_assert_str_eq(ret1[2], "test2");
    ck_assert_ptr_eq(ret1[3], NULL);
    free(ret1);

    char test_str[8];
    ret1 = argpar_split_str(str6, ",");
    ck_assert_ptr_ne(ret1, NULL);
    for (int i=0;i<10;i++){
        sprintf(test_str, "test%d", i+1);
        ck_assert_str_eq(ret1[i], test_str);
    }
    ck_assert_ptr_eq(ret1[10], NULL);
    free(ret1);

    free(str1); free(str2); free(str3); free(str4); free(str5); free(str6);
}
END_TEST

START_TEST(split_strings_trim){
    char *str1, *str3, **ret1;

    str1 = strdup("test1,test2,test3");
    str3 = strdup(" test1,test2 , test3 ");

    ret1 = argpar_split_trim(str1, ",");
    ck_assert_str_eq(ret1[0], "test1");
    ck_assert_str_eq(ret1[1], "test2");
    ck_assert_str_eq(ret1[2], "test3");
    ck_assert_ptr_eq(ret1[3], NULL);
    free(ret1);

    ret1 = argpar_split_trim(str3, ",");
    ck_assert_str_eq(ret1[0], "test1");
    ck_assert_str_eq(ret1[1], "test2");
    ck_assert_str_eq(ret1[2], "test3");
    ck_assert_ptr_eq(ret1[3], NULL);
    free(ret1);

    free(str1); free(str3);
}
END_TEST


START_TEST(split_int){
    char *str1, *str2, *str3, *str4;
    int *ret1;

    str1 = strdup("1, 2 , 3");
    str2 = strdup("1");
    str3 = strdup("asdf");
    str4 = strdup("1,,3");

    ret1 = argpar_split_int(str1, ",");
    ck_assert(ret1[0] == 1);
    ck_assert(ret1[1] == 2);
    ck_assert(ret1[2] == 3);
    ck_assert(ret1[3] == NULL_INT);
    free(ret1);

    ret1 = argpar_split_int(str2, ",");
    ck_assert(ret1[0] == 1);
    ck_assert(ret1[1] == NULL_INT);
    free(ret1);

    ret1 = argpar_split_int(str3, ",");
    ck_assert_ptr_eq(ret1, NULL);

    ret1 = argpar_split_int(str4, ",");
    ck_assert(ret1[0] == 1);
    ck_assert(ret1[1] == EMPTY_INT);
    ck_assert(ret1[2] == 3);
    ck_assert(ret1[3] == NULL_INT);
    free(ret1);

    free(str1); free(str2); free(str3); free(str4);
}
END_TEST

START_TEST(split_dbl){
    char *str1, *str2, *str3, *str4;
    double *ret1;

    str1 = strdup("1.1, 2.2 , 3.3");
    str2 = strdup("1.1");
    str3 = strdup("asdf");
    str4 = strdup("1.1,,3.3");

    ret1 = argpar_split_dbl(str1, ",");
    ck_assert(ret1[0] == 1.1);
    ck_assert(ret1[1] == 2.2);
    ck_assert(ret1[2] == 3.3);
    ck_assert(isnan(ret1[3]));
    free(ret1);

    ret1 = argpar_split_dbl(str2, ",");
    ck_assert(ret1[0] == 1.1);
    ck_assert(isnan(ret1[1]));
    free(ret1);

    ret1 = argpar_split_dbl(str3, ",");
    ck_assert_ptr_eq(ret1, NULL);

    ret1 = argpar_split_dbl(str4, ",");
    ck_assert(ret1[0] == 1.1);
    ck_assert(isinf(ret1[1]));
    ck_assert(ret1[2] == 3.3);
    ck_assert(isnan(ret1[3]));
    free(ret1);

    free(str1); free(str2); free(str3); free(str4);
}
END_TEST

START_TEST(accept_any){
    argpar_ostream = tmpfile();

    int argc = 2;
    char *argv_const[] = {(char*)argpar_program_name, "unknown=option"};
    char *argv[argc];
    int i;
    for (i=0;i<argc;i++){
        argv[i] = strdup(argv_const[i]);
    }
    arguments arguments = {
        .argc = 0, .argc_unknown = 0, .i = -1, .d = -1
    };
    ck_assert_int_eq(argpar_parse_args(&params, argc, argv, ARGPAR_ACCEPT_ANY, NULL, &arguments), 0);
    ck_assert_int_eq(arguments.unknown, 1);

    for (i=0;i<argc;i++){
        free(argv[i]);
    }

    clean_arguments(&arguments);

    argpar_clean(&params);
    fclose(argpar_ostream);
}
END_TEST

Suite* stub_suite(void){
    Suite *s = suite_create("Stub");

    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, no_args_usage);
    tcase_add_test(tc_core, no_real_args);
    tcase_add_test(tc_core, auto_help);
    tcase_add_test(tc_core, no_key_args);
    tcase_add_test(tc_core, unknown_option);
    tcase_add_test(tc_core, int_errors);
    tcase_add_test(tc_core, dbl_errors);
    tcase_add_test(tc_core, parfile1);
    tcase_add_test(tc_core, children1);
    tcase_add_test(tc_core, parse_parfile);
    tcase_add_test(tc_core, parse_parfile_in_child);
    tcase_add_test(tc_core, skip_parfiles);
    tcase_add_test(tc_core, no_keyargs);
    tcase_add_test(tc_core, no_keyargs_parfile);
    tcase_add_test(tc_core, split_strings);
    tcase_add_test(tc_core, split_strings_trim);
    tcase_add_test(tc_core, split_int);
    tcase_add_test(tc_core, split_dbl);
    tcase_add_test(tc_core, aliases);
    tcase_add_test(tc_core, accept_any);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(int argc, char **argv){

//    int argc2 = 3;
//    char *argv_const[] = {(char*)argpar_program_name, "help=params", "parfile=initial"};
//    char *argv2[argc];
//    int i;
//    for (i=0;i<argc2;i++){
//        argv2[i] = strdup(argv_const[i]);
//    }
//    arguments arguments = {0};
//
//    printf("Return value: %d\n", argpar_parse_args(&params, argc2, argv2, 0, NULL, &arguments));
//
//    clean_arguments(&arguments);
//
//    argpar_clean(&params);
//
//    for (i=0;i<argc2;i++){
//        free(argv2[i]);
//    }
//
//    return 0;

    int number_failed;

    Suite *s = stub_suite();
    SRunner *sr = srunner_create(s);

    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
