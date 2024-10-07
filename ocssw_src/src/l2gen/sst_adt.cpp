#include "sst_adt.hpp"

namespace adt {
    std::string config_adt_path;
    const std::unordered_set<std::string> non_test_keywords = {"true", "false", "thresholds", "addTreesum", "setFlags",
                                                               "decision_tree"};
    /**
     * @brief possible traversal functions.
     */
    std::map<std::string, Traverse> traversal_functions = {
            {"single_max",
                    [](float *treesum, int16_t *flags, const TestParameters &parameters, const VarsAtPixel &vars,
                       Treenode *node) {
                        if (vars.vars.at(parameters.test_name) < parameters.max_var &&
                            vars.masks.at(parameters.test_name)) {
                            node->children.at(1).second = false;
                        } else {
                            node->children.at(0).second = false;
                        }
                    }},
            {"single_min",
                    [](float *treesum, int16_t *flags, const TestParameters &parameters, const VarsAtPixel &vars,
                       Treenode *node) {
                        if (vars.vars.at(parameters.test_name) >= parameters.min_var) {
                            node->children.at(1).second = false;
                        } else {
                            node->children.at(0).second = false;
                        }
                    }},
            {"double_condition",
                    [](float *treesum, int16_t *flags, const TestParameters &parameters, const VarsAtPixel &vars,
                       Treenode *node) {
                        if (vars.vars.at(parameters.test_name) >= parameters.min_var &&
                            vars.vars.at(parameters.test_name) <= parameters.max_var) {
                            node->children.at(1).second = false;
                        } else {
                            node->children.at(0).second = false;
                        }
                    }},
            {"addTreesum",
                    [](float *treesum, int16_t *flags, const TestParameters &parameters, const VarsAtPixel &vars,
                       Treenode *node) {
                        *treesum += parameters.treesum_increment;
                    }},
            {"setFlags",
                    [](float *treesum, int16_t *flags, const TestParameters &parameters, const VarsAtPixel &vars,
                       Treenode *node) {
                        *flags |= parameters.flag;
                    }},
            {"default",
                    [](float *treesum, int16_t *flags, const TestParameters &parameters, const VarsAtPixel &vars,
                       Treenode *node) {
                    }},
    };

    auto set_node_parameters = [](Treenode *node, const std::string &test_name,
                                  const std::map<std::string, boost::variant<int16_t, float, bool>> &ranges) {
        auto params = &node->set_parameters;
        if (ranges.count("threshold_max") > 0)
            params->max_var = boost::get<float>(ranges.at("threshold_max"));
        if (ranges.count("threshold_min") > 0)
            params->min_var = boost::get<float>(ranges.at("threshold_min"));
        if (ranges.count("addTreesum"))
            params->treesum_increment = boost::get<float>(ranges.at("addTreesum"));
        if (ranges.count("setFlags"))
            params->flag = boost::get<int16_t>(ranges.at("setFlags"));
        if (ranges.count("incl"))
            params->not_use_mask = !boost::get<bool>(ranges.at("incl"));
        params->test_name = test_name; //
    };

    void deallocate_tree(Treenode *root) {
        if (root == nullptr)
            return;
        for (auto &child: root->children) {
            deallocate_tree(child.first);
        }
        delete root;
    };

    void tree_traversal(float *inp_treesum, int16_t *inp_flags_sst, const VarsAtPixel &vars, Treenode *node) {
        if (node == nullptr)
            return;
        node->traverse(inp_treesum, inp_flags_sst, vars);
        for (auto &child: node->children) {
            if (child.second) {
                tree_traversal(inp_treesum, inp_flags_sst, vars, child.first);
            }
        }
    }

    void tree_traversal(Treenode *node) {
        if (node == nullptr)
            return;
        for (auto &child: node->children) {
            child.second = true;
            tree_traversal(child.first);
        }
    }

    void readTreeTest(const rapidjson::Value &member) {
        for (rapidjson::Value::ConstMemberIterator iter = member.MemberBegin(); iter != member.MemberEnd(); ++iter) {
            if (iter->value.IsObject())
                readTreeTest(iter->value);
            if (iter->value.IsArray()) {
                const rapidjson::Value &array = iter->value;
                for (const rapidjson::Value &el: array.GetArray()) {
                    readTreeTest(el);
                }
            }
        }
    }

    void build_sub_tree(const rapidjson::Value &member, Treenode *node, std::set<std::string> &tests_name) {

        for (rapidjson::Value::ConstMemberIterator iter = member.MemberBegin(); iter != member.MemberEnd(); ++iter) {
            const std::string key = iter->name.GetString();
            {
                assert(iter->value.IsObject() == 1); //
            }
            const rapidjson::Value &nodes = iter->value;
            // disecting what is that
            if (non_test_keywords.count(key) == 0) //  it's  a test
            {
                std::map<std::string, boost::variant<int16_t, float, bool>> thresholds;
                tests_name.insert(key);
                // here it should read the threshold value
                if (nodes.FindMember("thresholds") != nodes.MemberEnd()) {
                    if (!nodes["thresholds"].IsArray()) {
                        print_error_message_for_adt(config_adt_path, key,
                                                    "Thresholds must be an array. Check the test ");
                    }
                    const rapidjson::Value &arr = nodes["thresholds"];
                    for (const auto &thresh_iter: arr.GetArray()) {
                        for (const auto &field: thresh_iter.GetObject()) {
                            if (field.value.IsFloat())
                                thresholds[field.name.GetString()] = field.value.GetFloat();
                            if (field.value.IsBool())
                                thresholds[field.name.GetString()] = field.value.GetBool();
                            if (field.value.IsInt())
                                thresholds[field.name.GetString()] = (int16_t) field.value.GetInt();
                        }
                    }
                    set_node_parameters(node, key, thresholds);
                } // defining the parameters

                if (thresholds.count("threshold_min") > 0 && thresholds.count("threshold_max") > 0)
                    node->traverse_ptr = traversal_functions.at("double_condition");
                else if (thresholds.count("threshold_min") > 0)
                    node->traverse_ptr = traversal_functions.at("single_min");
                else if (thresholds.count("threshold_max") > 0)
                    node->traverse_ptr = traversal_functions.at("single_max");
                else
                    node->traverse_ptr = traversal_functions.at("default");
            } else // this is a branch, now all the branches are binary, "true" or "false"
            {
                if (key != "true" && key != "false") {
                    print_error_message_for_adt(config_adt_path, key, "Undefined keyword for a branch ");
                }
                assert(((key == "true") || (key == "false")) == 1);
                // here it should read the threshold value
                std::map<std::string, boost::variant<int16_t, float, bool>> thresholds;
                if (nodes.FindMember("addTreesum") != nodes.MemberEnd()) {
                    node->traverse_ptr = traversal_functions.at("addTreesum");
                    assert(nodes["addTreesum"].IsFloat());
                    thresholds["addTreesum"] = nodes["addTreesum"].GetFloat();
                    set_node_parameters(node, "addTreesum_" + key, thresholds);
                    // looking for addTreesum or set flags ( decision has been made)}
                } else if (nodes.FindMember("setFlags") != nodes.MemberEnd()) {
                    node->traverse_ptr = traversal_functions.at("setFlags");
                    assert(nodes["setFlags"].IsInt());
                    thresholds["setFlags"] = (int16_t) nodes["setFlags"].GetInt();
                    set_node_parameters(node, "setFlags_" + key, thresholds);
                    // looking for addTreesum or set flags ( decision has been made)}
                } else {
                    node->traverse_ptr = traversal_functions.at("default");
                    set_node_parameters(node, "branch_node_" + key, thresholds);
                }
            }

            if (nodes.FindMember("children") != nodes.MemberEnd()) {
                const rapidjson::Value &children = nodes["children"];
                if (children.IsArray() != 1) {
                    print_error_message_for_adt(config_adt_path, key,
                                                "Children must be an array. Check the test ");
                };
                size_t number_of_children = children.Size();
                for (size_t _ = 0; _ < number_of_children; _++)
                    node->children.push_back(std::make_pair(new Treenode(), true));

                for (size_t i = 0; i < number_of_children; i++) {
                    build_sub_tree(children[i], node->children[i].first, tests_name);
                }
            }
        }
    }

    void build_tree(const std::string &config_path, Treenode *decision_tree, std::set<std::string> &tests_name) {
        config_adt_path = config_path;
        std::ifstream ifs{config_path};
        if (!ifs.is_open()) {
            std::cerr << "Could not open the ADT json file for reading : " << config_path << std::endl;
            exit(EXIT_FAILURE);
            return;
        } else {
            std::cout << "Reading the ADT parameters from the follwing file : " << config_path << "\n";
        }
        rapidjson::IStreamWrapper isw{ifs};
        rapidjson::Document doc{};
        doc.ParseStream(isw);
        const rapidjson::Value &tree_builder = doc["decision_tree"];
        readTreeTest(tree_builder);
        {
            std::cout << "The json file for the ADT has been read and processed succesfully"
                      << "\n";
        }
        if (tree_builder.FindMember("children") != tree_builder.MemberEnd()) {
            const rapidjson::Value &children = tree_builder["children"];
            if (children.IsArray() != 1) {
                print_error_message_for_adt(config_adt_path, "decision_tree",
                                            "Children must be an array. Check the test ");
            };
            size_t number_of_children = children.Size();
            decision_tree->set_parameters.test_name = "decision_tree";
            decision_tree->traverse_ptr = [](float *, int16_t *, const TestParameters &, const VarsAtPixel &,
                                             Treenode *) {};

            // check for additional traversal functions
            if (tree_builder.FindMember("addTreesum") != tree_builder.MemberEnd()) {
                std::map<std::string, boost::variant<int16_t, float, bool>> thresholds;
                decision_tree->traverse_ptr = traversal_functions.at("addTreesum");
                assert(tree_builder["addTreesum"].IsFloat());
                thresholds["addTreesum"] = tree_builder["addTreesum"].GetFloat();
                set_node_parameters(decision_tree, "decision_tree", thresholds);
                // looking for addTreesum or set flags ( decision has been made)}
            }
            for (size_t _ = 0; _ < number_of_children; _++)
                decision_tree->children.push_back(std::make_pair(new Treenode(), true));
            for (size_t i = 0; i < number_of_children; i++) {
                build_sub_tree(children[i], decision_tree->children[i].first, tests_name);
            }
        } else {
            std::cerr << "Error : Empty decision tree. Exiting " << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    void
    print_error_message_for_adt(const std::string &file_path, const std::string &key_word, const std::string &message) {
        std::ifstream file(file_path);
        if (!file.is_open()) {
            std::cerr << "Problem opening the ADT config " << file_path << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string str;
        size_t line_num = 0;
        while (std::getline(file, str)) {
            if (boost::algorithm::contains(str, key_word))
                std::cerr << message << " : " << key_word << " in line " << line_num << " in file " << file_path
                          << std::endl;
            line_num++;
        }
        file.close();
        exit(EXIT_FAILURE);
    }


}