
#include <argpar.h>

#include <stdlib.h>
#include <string.h>


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

static int parse_options(int key, char *argv, argpar_state *state) {
	return 0;
}

argpar params = { options, parse_options, args_doc, doc };

static const char doc_parent[] = "This is an argpar test.";
static const char args_doc_parent[] = "[none]";
static const argpar_option options_parent[] = {
	{ "dbl", 'd', 0, OPTION_DBL, "casted as dbl", 1 },
	{ "parg", 'p', 0, OPTION_INT, "new arg for parent", 1 },
	{0}
};

static int parse_options_parent(int key, char *argv, argpar_state *state) {
	return 0;
}

argpar_child children[] = {{&params}, {0}};
argpar params_parent = { options_parent, parse_options_parent, args_doc_parent, doc_parent, children };

int main(){
	int argc = 6;
	const char *argv_const[] = {
			argpar_program_name, "parfile=parfiles/initial.par", "ofile=keyarg_ofile1", "parfile=parfiles/final.par",
			"ofile=keyarg_ofile2", "help=params"
	};
	char *argv[argc];
	int i;
	for (i=0;i<argc;i++){
		argv[i] = strdup(argv_const[i]);
	}
	int ret = argpar_parse_args(&params, argc, argv, ARGPAR_STORE_PARAMS | ARGPAR_NO_EXIT, NULL, NULL);

	for (i=0;i<argc;i++){
		free(argv[i]);
	}

	argpar_clean(&params);
    return ret;
}
