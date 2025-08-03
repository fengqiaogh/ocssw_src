#ifndef SST_ADT_HPP
#define SST_ADT_HPP

#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/ostreamwrapper.h>
#include <boost/variant.hpp>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <chrono>
#include <unordered_map>
#include <memory>
#include <set>
/**
 * @brief namespace to contains the ADT data structures and functions *
 * @param TestParameters - contains the test name, the max and min threshold, treesum increment (positive or negative). Also optional boolean var and flag, would be needed for AVHRR implementation
 * @param VarsAtPixel - value of the quantity to be "thresholded"
 * @param Treenode - ADT structure
 * @param build_sub_tree - building a subtree from rapid json object
 * @param build_tree - building a tree from a json config file
 */
namespace adt {

    struct TestParameters {
        std::string test_name;
        float min_var, max_var;
        bool not_use_mask = true;
        float treesum_increment;
        int16_t flag;
    };
    struct VarsAtPixel {
        std::unordered_map<std::string, float> vars;
        std::unordered_map<std::string, bool> masks;
        size_t i_p;
        size_t i_scan;
    };

    /**
     * @brief The decision tree
     * @param children - decision branches
     * @param traverse_ptr - traversal function ptr
     * @param set_parameters - parameters of the test
     * @param traverse - traverse driver 
     */
    struct Treenode {
        std::vector<std::pair<Treenode *, bool>> children;

        void (*traverse_ptr)(float *, int16_t *, const TestParameters &, const VarsAtPixel &, Treenode *);

        TestParameters set_parameters;

        void traverse(float *inp_treesum, int16_t *inp_flags_sst, const VarsAtPixel &vars) {
            traverse_ptr(inp_treesum, inp_flags_sst, set_parameters, vars, this);
        }
    };

    void deallocate_tree(Treenode *root);

    typedef void (*Traverse)(float *, int16_t *, const TestParameters &, const VarsAtPixel &, Treenode *);

    void tree_traversal(float *inp_treesum, int16_t *inp_flags_sst, const VarsAtPixel &vars, Treenode *node);

    // resetting the tree to default. By default, of all the children are traversable.
    void tree_traversal(Treenode *node);

    void readTreeTest(const rapidjson::Value &member);

    /**
     * @brief Recursively building the tree     *
     * @param member - rapid json subtree
     * @param node - Treenode
     * @param tests_name - store test names found in the tree
     */
    void build_sub_tree(const rapidjson::Value &member, Treenode *node, std::set<std::string> &tests_name);

    /**
     * @brief Building the from a json file     *
     * @param config_path - path to config file
     * @param decision_tree - Treenode object which is being build
     * @param tests_name - test names recorded and later to be parsed
     */
    void build_tree(const std::string &config_path, Treenode *decision_tree, std::set<std::string> &tests_name);

    /**
   * @brief
   * Prints error message for unrecognized keywords in the ADT config file
   * @param file_path - the file path   *
   * @param key_word - key word
   * @param message - message to print
   */
    void print_error_message_for_adt(const std::string &file_path, const std::string &key_word,
                                     const std::string &message = "Undefined keyword");
}
#endif