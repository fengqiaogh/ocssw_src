from typing import Dict, List
import json
import os
# Builds ADT json files for l2gen SST from RSMAS txt trees.
# Txt trees are availible (for now ) only for VIIRS J1
if __name__ == "__main__":
    # attributes
    corr_map: Dict[str, str] = {"sst2b": "SST", "sst3b": "SST3", "sst2bmsst3b": "d_SST_SST3", "dm3750": "BT37_MAXMIN",
                                "m3750": "BT37", "sd.sst": "SST_STD",
                                "dm8550": "BT85_MAXMIN", "dm11000": "BT11", "dm12000": "BT12",
                                "dm3750m12000": "d_BT37_BT12",
                                "dm3750m11000": "d_BT37_BT11", "dm11000m12000": "d_BT11_BT12", "m8550": "BT85",
                                "lat": "lat", "lon": "lon", "anc.wv": "wv", "sd.m3750": "BT37_STD",
                                "sd.m11000": "BT11_STD", "sd.m12000": "SST_STD", "Tdeflong": "Tdeflong",
                                "maxm8550": "BT85_MAX", "minm8550": "BT85_MIN", "m.rho1380": "RHOCIRRUS",
                                "dm.rho678": "RHOred_MAXMIN", "m.rho678": "RHOred", " minm.rho678": "RHOred_MIN",
                                "maxm.rho678": "RHOred_MAX", "maxm.rho748": "RHO07_MAX", "minm.rho748": "RHO07_MIN",
                                "minm.rho1610": "RHO16_MIN", "maxm.rho1610": "RHO16_MAX",
                                "minm.rho1380": "RHOCIRRUS_MIN", "maxm.rho1380": "RHOCIRRUS_MAX", "m.rho1610": "RHO16",
                                "dm.rho1610": "RHO16_MAXMIN",
                                "dm.rho1380": "RHOCIRRUS_MAXMIN", "m.rho748": "RHO07", "maxm.rho.1380": "RHOCIRRUS_MAX"}


    class Tree:
        def __init__(self):
            self.children: List[Tree] = []
            self.threshold = None
            self.addTreesum = None
            self.testname = None
            self.node_type = None

        def add_children(self, children):
            self.children = children

        def set_threshold(self, thresholds: Dict[str, float]):
            self.threshold = thresholds

        def set_treesum(self, treesum):
            self.addTreesum = treesum

        def set_testname(self, name):
            self.testname = name

        def traverse(self):
            # print(self.testname, end="; ")
            # print(self.node_type, end="; ")
            # print(self.threshold, end="; ")
            # print(self.addTreesum)
            for child in self.children:
                child.traverse()


    def build_my_tree(input_tree: Tree):
        names = {}
        for i, child in enumerate(input_tree.children):
            if child.testname in names:
                names[child.testname].append(i)
            else:
                names[child.testname] = [i]
        new_children: List[Tree] = []
        for name, indexes in names.items():
            new_children.append(Tree())
            new_children[-1].testname = name
            new_children[-1].children = [input_tree.children[indexes[0]], input_tree.children[indexes[1]]]
            new_children[-1].threshold = input_tree.children[indexes[0]].threshold
            input_tree.children[indexes[0]].testname = "true"
            input_tree.children[indexes[1]].testname = "false"
            input_tree.children[indexes[0]].threshold = None
            input_tree.children[indexes[1]].threshold = None
        input_tree.children = new_children
        for child in new_children:
            for dec in child.children:
                build_my_tree(dec)
        pass


    def count_order(line: str, pattern="|  "):
        pat_len = len(pattern)
        count = 0
        conc: str = ""
        while conc == line[:count * pat_len]:
            conc += pattern
            count += 1
        return count - 1
        pass


    def traverse_list(lines: List[str], i: int, tree: Tree):
        pattern = "|  "
        # print("current node = ", lines[i][:-1], tree.testname)
        if (i > len(lines) - 2):
            # print("No children ", i)
            return i
        current_order = count_order(lines[i], pattern)
        j = i
        while j < len(lines) - 1:
            next_line = lines[j + 1]
            next_order = count_order(next_line, pattern)
            if next_order > current_order:
                tree.children.append(Tree())
                test_name = next_line[next_order * len(pattern):].split()[0]
                test_order = next_line[next_order * len(pattern):].split()[1]
                test_threshold = float(next_line[next_order * len(pattern):].split()[2][:-1])
                test_treesum = float(next_line[next_order * len(pattern):].split()[3])
                # print(
                #     f"rootname = {tree.testname}, childrenlen = {len(tree.children)},"
                #     f" childname = {test_name}, childorder  = {next_order}, test_order = {test_order}, test_threshold = {test_threshold}, test_treesum = {test_treesum} ")
                tree.children[-1].testname = test_name
                tree.children[-1].threshold = test_threshold
                tree.children[-1].addTreesum = test_treesum
                tree.children[-1].node_type = "threshold_max" if test_order == "<" else "threshold_min"
                j = traverse_list(lines, j + 1, tree.children[-1])
            else:
                return j
        return j
        pass


    def trim_word(inp_str: str):
        index = inp_str.find(")")
        return inp_str[index + 1:]


    def create_dict_from_tree(tree: Tree, json):
        arr = []
        test_name = ""
        for _ in range(len(tree.children)):
            arr.append({})
        if tree.testname == "decision_tree":
            test_name = "decision_tree"
            if len(arr) > 0:
                json[test_name] = {"children": arr}
        elif tree.testname != "true" and tree.testname != "false":
            thresholds = [{"threshold_max": tree.threshold}]
            test_name_or = trim_word(tree.testname)
            if test_name_or in corr_map:
                test_name = corr_map[test_name_or]
            else:
                print(f"{test_name_or} NOT FOUND")
                raise SystemExit('Wrong Keyword')
            json[test_name] = {"thresholds": thresholds, "children": arr}
        else:
            test_name = tree.testname
            if len(arr) > 0:
                json[test_name] = {"addTreesum": tree.addTreesum, "children": arr}
            else:
                json[test_name] = {"addTreesum": tree.addTreesum}
        for i, child in enumerate(tree.children):
            create_dict_from_tree(child, json[test_name]["children"][i])
        pass


    def get_the_tree(path_inp_file, out_path):
        # path_inp_file: str = "/Users/avsemeno/Downloads/NOAA20_R2022_SST_delivery/NOAA_20_ADtree_no_glint_R2022.0_results.txt" #/Users/avsemeno/Downloads/NOAA20_R2022_SST_delivery/NOAA_20_ADtree_mod_glint_R2022.0_results.txt" # /Users/avsemeno/Downloads/NOAA20_R2022_SST_delivery/NOAA_20_ADtree_night_R2022.0_results.txt"
        data = []
        append = False
        with open(path_inp_file, "r") as file:
            lines = file.readlines()
        for line in lines:
            if line == ": 0\n":
                append = True
            if "Legend" in line:
                append = False
            if append:
                data.append(line)
        adt = Tree()
        adt.testname = "decision_tree"
        traverse_list(data, 0, adt)
        build_my_tree(adt)
        adt.traverse()
        jsonf = {}
        create_dict_from_tree(adt, jsonf)
        # print(jsonf)

        with open(out_path, "w") as outfile:
            json.dump(jsonf, outfile, indent=4)
        return jsonf


    ocdata = os.environ.get('OCDATAROOT')
    if ocdata == None:
        raise SystemExit('OCDATAROOT not found')
    path_no_glint = os.path.join(ocdata, "viirs", "j1", "cal", "NOAA_20_ADtree_no_glint_R2022.0_results.txt")
    path_glint = os.path.join(ocdata, "viirs", "j1", "cal", "NOAA_20_ADtree_mod_glint_R2022.0_results.txt")
    path_night = os.path.join(ocdata, "viirs", "j1", "cal", "NOAA_20_ADtree_night_R2022.0_results.txt")
    tree_no_glint = get_the_tree(path_no_glint, "no_glint_day.json")
    tree_glint = get_the_tree(path_glint, "glint_day.json")
    adt_night = {"version": "v6.R2022"}
    tree_night = get_the_tree(path_night, "night.json")
    adt_night["decision_tree"] = tree_night["decision_tree"]
    adt_full = {"version": "v6.R2022", "decision_tree": {"children": [
        {"solz": {"thresholds": [{"threshold_min": 85.0}],
                  "children": [{"true": tree_night["decision_tree"]}, {"false": {"children": [
                      {"glintcoef": {"thresholds": [{"threshold_max": 0.005}],
                                     "children": [{"true": tree_no_glint["decision_tree"]},
                                                  {"false": tree_glint["decision_tree"]}]}}]}}]}}]}}

    with open("cloud_mask_sst.json", "w") as outfile:
        json.dump(adt_full, outfile, indent=2)
    with open("cloud_mask_sst3.json", "w") as outfile:
        json.dump(adt_night, outfile, indent=4)
    pass
