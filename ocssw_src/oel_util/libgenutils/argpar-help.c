#include "argpar.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define OPT_COL       2		/* column in which long options start */
#define DOC_OPT_COL   2		/* column in which doc options start */
#define OPT_DOC_COL  29		/* column in which option text starts */
#define DOC_COL       2		/* column in which group documents are printed */
#define ENUM_COL      8		/* column in which enum documents are printed */
#define ALIAS_COL     4		/* column in which alias documents are printed */
#define HEADER_COL    1		/* column in which group headers are printed */
#define WRAP_INDENT  12		/* indentation of wrapped lines */
#define RMARGIN      79		/* right margin used for wrapping */

static int print_wrapped(argpar* p, const char *doc, size_t *current_column, int indent, int filter_key) {
    char *text;
    if (filter_key && p->help_filter) {
        text = p->help_filter(filter_key, doc, NULL);
        if (text == doc) {
            text = strdup(doc);
            if (text == NULL) {
                return ENOMEM;
            }
        }
    } else if (doc) {
        text = strdup(doc);
        if (text == NULL) {
            return ENOMEM;
        }
    } else {
        return -1;
    }
    if (text) {
        char *d = text;
        size_t doc_length = strlen(d);
        switch (filter_key) {
        case ARGPAR_KEY_HELP_ARGS_DOC:
            doc_length++;
            fprintf(argpar_ostream, " ");
            break;
        case ARGPAR_KEY_HELP_POST_DOC:
            case ARGPAR_KEY_HELP_EXTRA:
            fprintf(argpar_ostream, "\n");
            if (current_column) {
                *current_column = 0;
            }
            break;
        }
        size_t desired_length = RMARGIN;
        if (current_column) {
            desired_length -= *current_column;
        }
        while (doc_length && doc_length > desired_length) {
            char *doc_end = d;
            size_t length_to_print = 0;

            while (length_to_print <= desired_length && *doc_end != '\n') {
                length_to_print++;
                doc_end++;
            }

            char tmp_char = 0;
            if (!isspace(*doc_end)) {
                while (!isspace(*doc_end) && length_to_print) {
                    doc_end--;
                    length_to_print--;
                }
                if (length_to_print == 0) {
                    length_to_print = desired_length;
                }
            }

            tmp_char = d[length_to_print];
            d[length_to_print] = '\0';

            fprintf(argpar_ostream, "%s\n", d);
            fprintf(argpar_ostream, "%*s", indent, "");

            if (tmp_char == '\n') {
                fprintf(argpar_ostream, "\n");
            }

            if (tmp_char) {
                d[length_to_print] = tmp_char;
                tmp_char = 0;
            }
            d += length_to_print;
            while (isspace(*d)) {
                d++;
            }
            doc_length = strlen(d);
            desired_length = RMARGIN;
        }
        fprintf(argpar_ostream, "%s", d);
        if (current_column) {
            *current_column = strlen(d);
        }
        switch (filter_key) {
        case ARGPAR_KEY_HELP_POST_DOC:
            case ARGPAR_KEY_HELP_EXTRA:
            fprintf(argpar_ostream, "\n");
            if (current_column) {
                *current_column = 0;
            }
            break;
        }
        if (text != doc) {
            free(text);
        }
    }
    return 0;
}

static const char *get_arg_text(const argpar_option *o) {
    const char *arg;
    if (o->arg) {
        arg = o->arg;
    } else if (o->flags & OPTION_DBL) {
        arg = "DBL";
    } else if (o->flags & OPTION_INT) {
        arg = "INT";
    } else {
        arg = "ARG";
    }
    return arg;
}

typedef enum {
    HEADER = 1, OPTION = 2, DOC = 3
} option_type;
typedef enum {
    MERGED = 0, NOT_MERGED = 1
} child_type;
typedef struct item_counts {
    unsigned header;
    unsigned option;
    unsigned doc;
} item_counts;

// all the group/option stuff should probably be pre-processed, like in argp (see clusters), but it works fine
static int print_all_options(const argpar_child *c, bool is_parent); // should prototype all group ones so this looks less weird

static int _print_group(const argpar_child *c, bool is_parent, option_type type, int group, child_type merged, item_counts *counts) {
//	printf("%s %p, group %d, type %d, not merged %d\n", c->header ? "child" : "parent", c->argpar, group, type, merged);
    argpar *p = c->argpar;
    bool printed = false;
    if ((merged == NOT_MERGED && c->header && group == c->group) || (merged == MERGED && !c->header)) {
        printed = true;
        if (c->header) {
//			if (counts->doc || counts->option || counts->header){
//				fprintf(argpar_ostream, "\n");
//			}
            if (strlen(c->header)) {
                fprintf(argpar_ostream, "\n%*s", HEADER_COL, "");
                size_t cur_column = HEADER_COL;
                if (print_wrapped(p, c->header, &cur_column, HEADER_COL, ARGPAR_KEY_HELP_HEADER)) {
                    return 1;
                }
//				fprintf(argpar_ostream, "\n");
                counts->header += 1;
            }
            argpar_child merged_child = { c->argpar, c->flags, NULL, 0 };
            print_all_options(&merged_child, false);
        } else {
            const argpar_option *o = p->options;
            int i = -1, cur_group = 0;
            const argpar_option *base_option = NULL;
            while (o[++i].name || o[i].doc) {
                int is_doc = o[i].flags & OPTION_DOC;
                if (!group && (o[i].group || (!o[i].name && !is_doc))) {
                    break;
                } else if (o[i].group) {
                    cur_group = o[i].group;
                } else if (!o[i].name && !is_doc) {
                    if (cur_group < 0) {
                        cur_group--;
                    } else {
                        cur_group++;
                    }
                }

                if (cur_group == group && !(o[i].flags & OPTION_HIDDEN) && !((o[i].flags & OPTION_PARENT) && !is_parent) && !((o[i].flags & OPTION_CHILD) && is_parent)) {
                    switch (type) {
                    case HEADER:
                        if (!o[i].name && !(o[i].flags & OPTION_DOC)) {
                            fprintf(argpar_ostream, "\n%*s", HEADER_COL, "");
                            size_t cur_column = HEADER_COL;
                            if (print_wrapped(p, o[i].doc, &cur_column, HEADER_COL, ARGPAR_KEY_HELP_HEADER)) {
                                return 1;
                            }
                            fprintf(argpar_ostream, "\n");
                            counts->header += 1;
                        }
                        break;
                    case OPTION:
                        if (o[i].name) {
                            if (!(counts->doc || counts->option || counts->header)) {
                                fprintf(argpar_ostream, "\n");
                            }
                            if (!(o[i].flags & (OPTION_DOC | OPTION_ALIAS))){
                            	base_option = &o[i];
                            }
                            size_t cur_column = 0;
                            if ((o[i].flags & OPTION_ENUM) == OPTION_ENUM) {
                                fprintf(argpar_ostream, "%*s%s", ENUM_COL, "", o[i].name);
                                cur_column = ENUM_COL + strlen(o[i].name);
                            } else if ((o[i].flags & OPTION_ALIAS) == OPTION_ALIAS) {
								fprintf(argpar_ostream, "%*s%s", ALIAS_COL, "", o[i].name);
								cur_column = ALIAS_COL + strlen(o[i].name);
                            } else if (o[i].flags & OPTION_DOC) {
                                fprintf(argpar_ostream, "%*s%s", DOC_OPT_COL, "", o[i].name);
                                cur_column = DOC_OPT_COL + strlen(o[i].name);
                            } else {
                                fprintf(argpar_ostream, "%*s", OPT_COL, "");
                                const char *arg = get_arg_text(&o[i]);
                                fprintf(argpar_ostream, "%s=%s", o[i].name, arg);
                                cur_column = OPT_COL + strlen(o[i].name) + 1 + strlen(arg); // equals sign = 1
                            }

                            if (o[i].doc) {
                                if ((cur_column + 1) > OPT_DOC_COL) { // +1 to ensure a space between opt and doc
                                    fprintf(argpar_ostream, "\n%*s", (int) (cur_column = WRAP_INDENT), "");
                                } else {
                                    fprintf(argpar_ostream, "%*s", (int) (OPT_DOC_COL - cur_column), "");
                                    cur_column = OPT_DOC_COL;
                                }
                                if (print_wrapped(p, o[i].doc, &cur_column, WRAP_INDENT, o[i].key)) {
                                    return 1;
                                }
                            } else if ((o[i].flags & OPTION_ALIAS) == OPTION_ALIAS){
                                if ((cur_column + 1) > OPT_DOC_COL) { // +1 to ensure a space between opt and doc
                                    fprintf(argpar_ostream, "\n%*s", (int) (cur_column = WRAP_INDENT), "");
                                } else {
                                    fprintf(argpar_ostream, "%*s", (int) (OPT_DOC_COL - cur_column), "");
                                    cur_column = OPT_DOC_COL;
                                }
                                char new_doc[12 + strlen(base_option->name)];
                                sprintf(new_doc, "Alias for '%s'", base_option->name);
                                if (print_wrapped(p, new_doc, &cur_column, WRAP_INDENT, o[i].key)) {
                                    return 1;
                                }
                            }
                            fprintf(argpar_ostream, "\n");
                            counts->option += 1;
                        }
                        break;
                    case DOC:
                        if (!o[i].name && (o[i].flags & OPTION_DOC)) {
//							if (counts->doc || counts->option || counts->header){
                            fprintf(argpar_ostream, "\n");
//							}
                            fprintf(argpar_ostream, "%*s", DOC_COL, "");
                            size_t cur_column = DOC_COL;
                            if (print_wrapped(p, o[i].doc, &cur_column, DOC_COL, ARGPAR_KEY_HELP_OPTION_DOC)) {
                                return 1;
                            }
                            if (!(o[i].flags & OPTION_DOC_NO_BREAK)) {
                                fprintf(argpar_ostream, "\n");
                            }
                            counts->doc += 1;
                        }
                        break;
                    }
                }
            }
        }
    }

    if (p->children && (printed || merged == NOT_MERGED)) {
        const argpar_child *c;
        for (c = p->children; c->argpar; c++) {
            if (_print_group(c, false, type, group, merged, counts)) {
                return 1;
            }
        }
    }

    return 0;
}

static int print_group(const argpar_child *c, bool is_parent, int group, child_type merged, item_counts *counts) {
    option_type type;
    if (merged == MERGED) {
        for (type = HEADER; type <= DOC; type++) {
            if (_print_group(c, is_parent, type, group, merged, counts)) {
                return 1;
            }
        }
    } else {
        if (_print_group(c, is_parent, HEADER, group, merged, counts)) {
            return 1;
        }
    }
    return 0;
}

static int print_args(const argpar_child *c, size_t *cur_column) {
    argpar *p = c->argpar;
    if (p) {
        if (p->args_doc) {
            print_wrapped(p, p->args_doc, cur_column, WRAP_INDENT, ARGPAR_KEY_HELP_ARGS_DOC);
        }
//		if (p->options){
//			const argpar_option *o = p->options;
//			int i = -1;
//			while (o[++i].name || o[i].doc){
//				if (!o[i].name || (o[i].flags & OPTION_HIDDEN) || (o[i].flags & OPTION_NO_USAGE) || (o[i].flags & OPTION_DOC)){
//					continue;
//				} else if ((o[i].flags & OPTION_ARG_OPTIONAL) == 0){ // required option
//					const char *arg = get_arg_text(o);
//					int adding_cols = strlen(o[i].name) + strlen(arg) + 1;
//					if ((*cur_column + adding_cols) > RMARGIN){
//						fprintf(argpar_ostream, "\n%*s", WRAP_INDENT, "");
//						*cur_column = WRAP_INDENT;
//					}
//					fprintf(argpar_ostream, "%s=%s ", o[i].name, arg);
//					*cur_column += adding_cols;
//				}
//			}
//		}
        if (p->children) {
            const argpar_child *c;
            for (c = p->children; c->argpar; c++) {
                if (print_args(c, cur_column)) {
                    return 1;
                }
            }
        }
    }
    return 0;
}

typedef struct group_limits {
    int max, min;
} group_limits;

static int find_group_limits(const argpar_child *c, group_limits *l) {
    argpar *p = c->argpar;
    if (p && p->options) {
        if (c->group > l->max) {
            l->max = c->group;
        }
        if (c->group < l->min) {
            l->min = c->group;
        }
        if (!c->header) {
            const argpar_option *o = p->options;
            int i = -1, cur_group = 0;
            while (o[++i].name || o[i].doc) {
                if (o[i].group) {
                    cur_group = o[i].group;
                } else if (!o[i].name && !(o[i].flags & OPTION_DOC)) {
                    if (cur_group < 0) {
                        cur_group--;
                    } else {
                        cur_group++;
                    }
                }
                if (cur_group > l->max) {
                    l->max = cur_group;
                }
                if (cur_group < l->min) {
                    l->min = cur_group;
                }
            }
            if (p->children) {
                const argpar_child *c;
                for (c = p->children; c->argpar; c++) {
                    if (find_group_limits(c, l)) {
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

static int print_all_options(const argpar_child *c, bool is_parent) {
    group_limits g_lims = { .min = INT_MAX, .max = INT_MIN };
    find_group_limits(c, &g_lims);

    int group_i; // really need to get rid of the following copy pasta
    for (group_i = 0; group_i <= g_lims.max; group_i++) {
        item_counts counts = { 0, 0, 0 };
        if (print_group(c, is_parent, group_i, MERGED, &counts)) {
            return 1;
        }
        if (print_group(c, is_parent, group_i, NOT_MERGED, &counts)) {
            return 1;
        }
    }
    for (group_i = g_lims.min; group_i < 0; group_i++) {
        item_counts counts = { 0, 0, 0 };
        if (print_group(c, is_parent, group_i, MERGED, &counts)) {
            return 1;
        }
        if (print_group(c, is_parent, group_i, NOT_MERGED, &counts)) {
            return 1;
        }
    }
    return 0;
}

static int print_ending_doc(const argpar_child *c) {
    argpar *p = c->argpar;
    if (p->doc) {
        char *ending_doc = strchr(p->doc, '\v');
        if (ending_doc) {
            int printed = print_wrapped(p, ending_doc + 1, 0, 0, ARGPAR_KEY_HELP_POST_DOC);
            if (printed > 0) {
                return 1;
            } else if (printed) {
                fprintf(argpar_ostream, "\n");
            }
        }
    }

    if (p->children) {
        const argpar_child *c;
        for (c = p->children; c->argpar; c++) {
            if (print_ending_doc(c)) {
                return 1;
            }
        }
    }
    return 0;
}

#if HAVE_DECL_PROGRAM_INVOCATION_NAME
// stolen directly from argp
static char *nondestructive_basename(char *name) {
    char *short_name = strrchr (name, '/');
    return short_name ? short_name + 1 : name;
}
#endif

// mostly stolen from argp
static char *program_name() {
#if HAVE_DECL_PROGRAM_INVOCATION_SHORT_NAME
    return program_invocation_short_name;
#elif HAVE_DECL_PROGRAM_INVOCATION_NAME
    return nondestructive_basename(program_invocation_name);
#else
    return (argpar_program_name ? (char*) argpar_program_name : "<cmd>");
#endif
}

static int _argpar_print_usage(argpar *p, unsigned flags) {
    const argpar_child base_argpar = { p, flags, NULL, 0 };
    if (!argpar_ostream) {
        argpar_ostream = stderr;
    }
    char *progname = program_name();
    fprintf(argpar_ostream, "Usage: %s", progname);

    size_t cur_column = 8 + strlen(progname); // "Usage: " + space padding after program name

    print_args(&base_argpar, &cur_column);
    fprintf(argpar_ostream, "\n");

    if (p->doc) {
        char *ending_doc = strchr(p->doc, '\v');
        if (ending_doc) {
            size_t new_doc_length = strlen(p->doc) - strlen(ending_doc);
            char new_doc[new_doc_length + 1];
            strncpy(new_doc, p->doc, new_doc_length);
            new_doc[new_doc_length] = '\0';
            if (print_wrapped(p, new_doc, 0, 0, ARGPAR_KEY_HELP_PRE_DOC)) {
                return 1;
            }
        } else {
            if (print_wrapped(p, p->doc, 0, 0, ARGPAR_KEY_HELP_PRE_DOC)) {
                return 1;
            }
        }
        fprintf(argpar_ostream, "\n");
    }
//	fprintf(argpar_ostream, "\n");

    if (print_all_options(&base_argpar, true)) {
        return 1;
    }

    if (p->children) {
        print_ending_doc(&base_argpar);
    }

    print_wrapped(p, NULL, 0, 0, ARGPAR_KEY_HELP_EXTRA);
    if (p->children) {
        const argpar_child *c;
        for (c = p->children; c->argpar; c++) {
            if (print_wrapped(c->argpar, NULL, 0, 0, ARGPAR_KEY_HELP_EXTRA)) {
                return 1;
            }
        }
    }
    return 0;
}

int argpar_usage(argpar_state *state) {
    return _argpar_print_usage(state->argpar, state->flags);
}
int argpar_usage_default(argpar *argpar) {
    return _argpar_print_usage(argpar, 0);
}
int argpar_help(argpar *argpar, FILE *stream, unsigned flags, char *name) {
    argpar_ostream = stream;
    if (name != NULL) {
        argpar_program_name = name;
    }
    return _argpar_print_usage(argpar, flags);
}
