#include <productInfo.h>

#include <stdlib.h>
#include <stdio.h>
#include <genutils.h>
#include <string>
#include <ctype.h>
#include <algorithm>

#include <pugixml.hpp>

using namespace std;
using namespace pugi;

extern "C" const char* getGCMDKeywords(const char* suite) {

    static string keyStr;

    xml_document rootNode;
    xml_node GCMDKeywordsNode;

    char *dataRoot;
    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        printf("-E- OCDATAROOT environment variable is not defined.\n");
        exit(1);
    }
    string GCMDKeywordsXMLFileName = (string) dataRoot + "/common/GCMDKeywords.xml";

    xml_parse_result xmlResult = rootNode.load_file(GCMDKeywordsXMLFileName.c_str());
    if (!xmlResult) {

        // if the XML file does not exist then we will not bother with the GCMDKeywords
//        if (xmlResult.status == status_file_not_found)
//            return NULL;

        printf("-E- %s Line %d: Can not load %s.  %s\n", __FILE__, __LINE__,
                GCMDKeywordsXMLFileName.c_str(), xmlResult.description());
        exit(EXIT_FAILURE);
    }

    GCMDKeywordsNode = rootNode.child("GCMDKeywords");
    if (!GCMDKeywordsNode) {
        printf("-E- %s Line %d: could not find GCMDKeywords tag in XML file = %s\n",
                __FILE__, __LINE__, GCMDKeywordsXMLFileName.c_str());
        exit(EXIT_FAILURE);
    }
    
    // search for the requested suite name
    xml_node suiteNode = GCMDKeywordsNode.find_child_by_attribute("suite", "name", suite);
    if (!suiteNode)
        return NULL;

    xml_node keywordNode = suiteNode.child("keyword");
    bool firstOne = true;
    while (keywordNode) {
        if(firstOne) {
            firstOne = false;
            keyStr.clear();
        } else {
            keyStr += "; ";
        }
        keyStr += keywordNode.child_value();
        keywordNode = keywordNode.next_sibling("keyword");
    }
    
    if(firstOne)
        return NULL;
    else
        return keyStr.c_str();
}
