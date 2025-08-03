#include <genutils.h>
#include <fstream> 
#include <algorithm>
#include <vector>
#include <cctype> 
#include <netcdf>
#include <string.h>
using namespace netCDF;

int check_url(const char *file) {
    return (strstr(file, "s3://") != NULL) || (strstr(file, "https://") != NULL) ||
           (strstr(file, "http://") != NULL);
}

std::vector<std::string> readFileList(std::string fileName) { 

    //Declare vector 
    std::vector<std::string> fileList;

    
    //Read file 
    std::ifstream file(fileName);
    if (check_url(fileName.c_str())) {
        fileList.push_back(fileName);
        return fileList;
    }
    if (!file && !check_url(fileName.c_str())) {
        std::cerr << "Error: The file " << fileName << " could not be read." << std::endl;
        exit(EXIT_FAILURE);
    }
    
    char buffer[1024];
    file.read(buffer, sizeof(buffer));

    for (int i=0; i<file.gcount(); i++) {

        //If file is not ASCII
        if (!(isprint(buffer[i]) || isspace(buffer[i]))) {
            fileList.push_back(fileName);
            return fileList; 
        }
    }

    // reset file
    file.close();
    file.open(fileName);
    
    //Check that each file listed in ifile exists 
    if (file.is_open()) {
        std::string line;

        while (std::getline(file, line)) {

            //remove comment after # sign
            size_t pos = line.find("#");
            if (pos != std::string::npos)
                line.erase(pos, std::string::npos);
        
            //remove any whitespace
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace),line.end());

            // ignore empty line
            if (line.empty())
                continue;

            //check if each file in ifile exists
            std::ifstream file_exists(line);
            if (file_exists.good())
                fileList.push_back(line);
            else {
                if (!check_url(line.c_str())) {
                    std::cerr << "Error: The file " << line << " could not be found." << std::endl;
                    exit(EXIT_FAILURE);
                }
		fileList.push_back(line);
            }

        }
        file.close();
    }
    return fileList;
}
