#include <cstdio>
#include <cstdlib>
#include "l1c_processor.h"
#include "bin_l1b.h"
#include "version_l1cgen.h"

int main(int argc, char *argv[]) {
    // make sure log files get updated for each line of output
    setlinebuf(stdout);
    setlinebuf(stderr);
    std::vector<std::string> cli;
    for (int i = 0; i < argc; i++) {
        cli.emplace_back(argv[i]);
    }
    std::string git_sha = GIT_COMMIT_HASH;
    std::string full_version = version  + std::string(" ") + git_sha;
    L1CProcessor l1_c_processor(program_name, full_version, cli);
    if (!l1c_input(argc, argv, l1_c_processor, program_name, full_version)) {
        fprintf(stderr, "-E-: %s:%d l1cgen: failed to parse/load command line arguments\n", __FILE__,
                __LINE__);
        exit(1);
    }
    l1_c_processor.binL1Bgranules();

    return (EXIT_SUCCESS);
}