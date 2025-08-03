#include "scaled_nc_var.hpp"
#include "l2_reader.hpp"
// Simple program to test l2_reader
// IT is not a part of CTEST. It needs to be developed later.
int main(int argc, char** argv) {
    std::string filename;
    std::string varname;
    if (argc > 2) {
        filename = argv[1];
        varname = argv[2];
    }
    else return 1;
    L2_Reader test = L2_Reader(filename);
    // auto var = find_nc_variable_cpp("par", test->nc_);
    ScaledNcVar var = test.get_variable(varname);
    size_t dim1 = var.getDims()[0].getSize();
    size_t dim2 = var.getDims()[1].getSize();
    std::vector<float> data(dim1 * dim2);
    var.getVar(data.data());
    size_t bound1 = 0;
    size_t bound2 = 0;
    size_t slice = 10;
    if (argc > 5) {
        bound1 = atoi(argv[3]);
        bound2 = atoi(argv[4]);
        slice = atoi(argv[5]);
    }
    for (size_t i =bound1; i < bound1 + slice; i++)
        for (size_t j =bound2; j < bound2 +  slice; j++)
            std::cout << data[i * dim2 + j] << " ";
    test.reset_cache();
    return 0;
}