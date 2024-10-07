#include "argpar.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <jansson.h>


static int print_wrapped(json_t *json, argpar* p, const char *key, const char *doc, int filter_key) {
//    printf("%s\n", __PRETTY_FUNCTION__);
    char *text;
    if (filter_key && p->help_filter) {
        text = p->help_filter(filter_key, doc, NULL);
    } else if (doc) {
        text = (char*)doc;
    } else {
        return 0;
    }
    if (text) {
        json_object_set_new(json, key, json_string(doc));

        if (text != doc) {
            free(text);
        }
    }
    return 0;
}

static const char *get_arg_text(const argpar_option *o) {
//    printf("%s\n", __PRETTY_FUNCTION__);
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
static int print_all_options(json_t *json, const argpar_child *c, bool is_parent); // should prototype all group ones so this looks less weird

static int _print_group(json_t *json, const argpar_child *c, bool is_parent, option_type type, int group, child_type merged, item_counts *counts) {
//    printf("%s\n", __PRETTY_FUNCTION__);
//	printf("%s %p, group %d, type %d, not merged %d\n", c->header ? "child" : "parent", c->argpar, group, type, merged);
    argpar *p = c->argpar;
    bool printed = false;
    if ((merged == NOT_MERGED && c->header && group == c->group) || (merged == MERGED && !c->header)) {
        printed = true;
        if (c->header) {
            if (strlen(c->header)) {
                if (print_wrapped(json, p, "header", c->header, ARGPAR_KEY_HELP_HEADER)) {
                    return 1;
                }
                counts->header += 1;
            }
            argpar_child merged_child = { c->argpar, c->flags, NULL, 0 };
            print_all_options(json, &merged_child, false);
        } else {
            const argpar_option *o = p->options;
            int i = -1, cur_group = 0;
            json_t *valid_values = NULL;
            json_t *aliases = NULL;
            json_t *current_option = NULL;
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

                if (cur_group == group && !((o[i].flags & OPTION_PARENT) && !is_parent) && !((o[i].flags & OPTION_CHILD) && is_parent)) {
//				if (cur_group == group && !(o[i].flags & OPTION_HIDDEN) && !((o[i].flags & OPTION_PARENT) && !is_parent) && !((o[i].flags & OPTION_CHILD) && is_parent)) {
                    switch (type) {
                    case HEADER:
                        if (!o[i].name && !(o[i].flags & OPTION_DOC)) {
//                            if (print_wrapped(json, p, "header", o[i].doc, ARGPAR_KEY_HELP_HEADER)) {
//                                return 1;
//                            }
                            counts->header += 1;
                        }
                        break;
                    case OPTION:
                        if (o[i].name) {
                            if ((o[i].flags & OPTION_ATTR) == OPTION_ATTR) {
                                json_object_set_new(current_option, o[i].name, json_string(o[i].doc));
                            } else if ((o[i].flags & OPTION_ENUM) == OPTION_ENUM) {
                            	if (valid_values == NULL){
                            		valid_values = json_array();
                            		json_object_set_new(current_option, "validValues", valid_values);
                            	}

                                json_t *valid_value = json_object();
                            	json_array_append_new(valid_values, valid_value);
                            	json_object_set_new(valid_value, "value", json_string(o[i].name));
                                if (o[i].doc){
                                	json_object_set_new(valid_value, "description", json_string(o[i].doc));
                                }
                            } else if ((o[i].flags & OPTION_ALIAS) == OPTION_ALIAS) {
                            	if (aliases == NULL){
                            		aliases = json_array();
                            		json_object_set_new(current_option, "aliases", aliases);
                            	}
                            	json_array_append_new(aliases, json_string(o[i].name));
                            } else if (o[i].flags & OPTION_DOC) {
                                json_object_set_new(current_option, "doc", json_string(o[i].name));
                            } else {
								current_option = json_object();
	                            json_array_append_new(json, current_option);
                                const char *arg = get_arg_text(&o[i]);
                                json_object_set_new(current_option, "name", json_string(o[i].name));
                                json_object_set_new(current_option, "arg", json_string(arg));

                                if (o[i].doc) {
                                    if (print_wrapped(current_option, p, "description", o[i].doc, o[i].key)) {
                                        return 1;
                                    }
                                }

                                valid_values = NULL;
                            }

                            counts->option += 1;
                        }
                        break;
                    case DOC:
                        if (!o[i].name && (o[i].flags & OPTION_DOC)) {
//                            if (print_wrapped(json, p, "doc", o[i].doc, ARGPAR_KEY_HELP_OPTION_DOC)) {
//                                return 1;
//                            }
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
            if (_print_group(json, c, false, type, group, merged, counts)) {
                return 1;
            }
        }
    }

    return 0;
}

static int print_group(json_t *json, const argpar_child *c, bool is_parent, int group, child_type merged, item_counts *counts) {
//    printf("%s\n", __PRETTY_FUNCTION__);
    option_type type;
    if (merged == MERGED) {
        for (type = HEADER; type <= DOC; type++) {
            if (_print_group(json, c, is_parent, type, group, merged, counts)) {
                return 1;
            }
        }
    } else {
        if (_print_group(json, c, is_parent, HEADER, group, merged, counts)) {
            return 1;
        }
    }
    return 0;
}

static int print_args(json_t *json, const argpar_child *c) {
//    printf("%s\n", __PRETTY_FUNCTION__);
    argpar *p = c->argpar;
    if (p) {
        if (p->args_doc) {
            print_wrapped(json, p, "args", p->args_doc, ARGPAR_KEY_HELP_ARGS_DOC);
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
                if (print_args(json, c)) {
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
//    printf("%s\n", __PRETTY_FUNCTION__);
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

static int print_all_options(json_t *json, const argpar_child *c, bool is_parent) {
//    printf("%s\n", __PRETTY_FUNCTION__);
    group_limits g_lims = { .min = INT_MAX, .max = INT_MIN };
    find_group_limits(c, &g_lims);

    int group_i; // really need to get rid of the following copy pasta
    for (group_i = 0; group_i <= g_lims.max; group_i++) {
        item_counts counts = { 0, 0, 0 };
        if (print_group(json, c, is_parent, group_i, MERGED, &counts)) {
            return 1;
        }
        if (print_group(json, c, is_parent, group_i, NOT_MERGED, &counts)) {
            return 1;
        }
    }
    for (group_i = g_lims.min; group_i < 0; group_i++) {
        item_counts counts = { 0, 0, 0 };
        if (print_group(json, c, is_parent, group_i, MERGED, &counts)) {
            return 1;
        }
        if (print_group(json, c, is_parent, group_i, NOT_MERGED, &counts)) {
            return 1;
        }
    }
    return 0;
}

static int print_ending_doc(json_t *json, const argpar_child *c) {
//    printf("%s\n", __PRETTY_FUNCTION__);
    argpar *p = c->argpar;
    if (p->doc) {
        char *ending_doc = strchr(p->doc, '\v');
        if (ending_doc) {
            int printed = print_wrapped(json, p, "post_doc", NULL, ARGPAR_KEY_HELP_POST_DOC);
            if (printed > 0) {
                return 1;
            }
        }
    }

    if (p->children) {
        const argpar_child *c;
        for (c = p->children; c->argpar; c++) {
            if (print_ending_doc(json, c)) {
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

    json_t *json = json_object();
    json_t *json_program = json_object();
    json_object_set_new(json, "program", json_program);
    json_object_set_new(json_program, "name", json_string(progname));

    print_args(json_program, &base_argpar);

    if (p->doc) {
        char *ending_doc = strchr(p->doc, '\v');
        if (ending_doc) {
            size_t new_doc_length = strlen(p->doc) - strlen(ending_doc);
            char new_doc[new_doc_length + 1];
            strncpy(new_doc, p->doc, new_doc_length);
            new_doc[new_doc_length] = '\0';
            if (print_wrapped(json_program, p, "pre_doc", new_doc, ARGPAR_KEY_HELP_PRE_DOC)) {
                return 1;
            }
        } else {
            if (print_wrapped(json_program, p, "pre_doc", p->doc, ARGPAR_KEY_HELP_PRE_DOC)) {
                return 1;
            }
        }
    }

    json_t *json_options = json_array();
    json_object_set_new(json, "options", json_options);
    if (print_all_options(json_options, &base_argpar, true)) {
        return 1;
    }

    if (p->children) {
        print_ending_doc(json, &base_argpar);
    }

    print_wrapped(json_program, p, "help_extra", NULL, ARGPAR_KEY_HELP_EXTRA);
    if (p->children) {
        const argpar_child *c;
        for (c = p->children; c->argpar; c++) {
            if (print_wrapped(json_program, c->argpar, "help_extra", NULL, ARGPAR_KEY_HELP_EXTRA)) {
                return 1;
            }
        }
    }

#if JANSSON_VERSION_HEX >= 0x020800
    int ret = json_dumpf(json, argpar_ostream, JSON_INDENT(4));
#else
    int ret = json_dumpf(json, argpar_ostream, JSON_INDENT(4) | JSON_PRESERVE_ORDER);
#endif

    ret |= json_object_clear(json);
    json_decref(json);

    return ret;
}


int argpar_usage_json(argpar_state *state) {
    return _argpar_print_usage(state->argpar, state->flags);
}
int argpar_usage_default_json(argpar *argpar) {
    return _argpar_print_usage(argpar, 0);
}
int argpar_help_json(argpar *argpar, FILE *stream, unsigned flags, char *name) {
    argpar_ostream = stream;
    if (name != NULL) {
        argpar_program_name = name;
    }
    return _argpar_print_usage(argpar, flags);
}
