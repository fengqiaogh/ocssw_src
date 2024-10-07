#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include "wavelength_3d.h"

using namespace std;

extern "C" void get_wavelength3d(const filehandle *l1file, instr *input) {

    // Convert user wavelength_3d_str into string
    string wavelength_3d_str = input->wavelength_3d_str;

    //--------------------------------------------------------------------------------------------------------
    //If user input for wavelength_3d_str is empty, populate wavelength_3d_str indexes and size with all wavelengths
    if (wavelength_3d_str.empty()) {
        input->wavelength_3d = l1file->iwave;
        input->nwavelengths_3d = l1file->nbands;
        input->wavelength_3d_index = (int *)malloc(l1file->nbands * sizeof(int *));
        for (int i = 0; i < l1file->nbands; i++) {
            input->wavelength_3d_index[i] = i;
        }
    } 
    else { // If wavelength_3d user input array not empty
        //--------------------------------------------------------------------------------------------------------
        // Check if user wavelength input is in ascending order:
        int size = 0;
        for (int i = 0; input->wavelength_3d_str[i] != 0; i++) {
            size++;
        }
        // Replace : with , for input wavelength_3d array
        char *temp_wave = (char*)malloc(size * sizeof(char*));
        for (int i = 0; i < size; i++) {
            temp_wave[i] = input->wavelength_3d_str[i];
            if (temp_wave[i] == ':') {
                temp_wave[i] = ',';
            }
        }
        temp_wave[size] = '\0';

        string temp_wave_str = temp_wave;
        int size_wave_ptr = 1;
        char **array = (char **)malloc(size * sizeof(char **));
        array[0] = &temp_wave[0];

        // Save comma separated inputs  
        for (int i = 1; temp_wave[i] != 0; i++) {
            if (temp_wave[i] == ',') {
                size_wave_ptr++;
                array = (char **)realloc(array, size_wave_ptr * sizeof(char **));
                array[(size_wave_ptr) - 1] = &temp_wave[i + 1];
                temp_wave[i] = '\0';
            }  
        }

        // if wavelengths are not in ascending order, exit
        for (int i = 0; i < size_wave_ptr - 1; i++) {
            if (atoi(array[i + 1]) <= atoi(array[i])) {
                cout << "Error: Input wavelengths in \"wavelength_3d\" not in ascending order at " << array[i+1] << endl;
                exit(1);
            }
        }
        //--------------------------------------------------------------------------------------------------------
        //  If input wavelengths are in ascending order, check if they are all found in input->iwave: 
        // If not all wavelengths found, then prompt error
        for (int i = 0; i < size_wave_ptr; i++) {
            bool found = false;
            for (int j = 0; j < (int)l1file->nbands; j++) {
                if(atoi(array[i]) == l1file->iwave[j]) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                cout << "Error: Input wavelength " << array[i] << " in \"wavelength_3d\" not found." << endl;
                exit(1);
            }
        }
        //--------------------------------------------------------------------------------------------------------
        // If input wavelengths are in ascending order and all found in input->iwave:
        // Convert l1file->iwave to int vector
        vector<int> iwave(l1file->iwave, l1file->iwave + l1file->nbands);

        // out_iwave array
        vector<int> out_iwave;

        // Position of current divider
        int pos = 0;
        int num1 = 0;
        int num2 = 0;

        // Iterate through wavelength_3d_str string
        while (pos < (int)wavelength_3d_str.length()) {
            // Find next divider
            int div = wavelength_3d_str.find(",", pos);
            // Check if end of string
            if (div == (int)string::npos) {
                div = wavelength_3d_str.length();
            }

            // Check if range (contains ':')
            int range = wavelength_3d_str.find(":", pos);

            if (range != (int)string::npos && range < div) {
                // Split range into start and end
                int start = stoi(wavelength_3d_str.substr(pos, range - pos));
                int end = stoi(wavelength_3d_str.substr(range + 1, div - range - 1));
    
                // Add numbers in range to out_iwave
                for (int i = start; i <= end; i++) {
                    out_iwave.push_back(i);
                    num1++;
                }
            }
        
            else {
                // Is single number, add it to out_iwave
                int num = stoi(wavelength_3d_str.substr(pos, div - pos));
                out_iwave.push_back(num);
                num2++;
            }

            // Update position
            pos = div + 1;
        }

        int n = 0;
        // Allocate indexes array memory
        input->wavelength_3d_index = (int *)malloc(n * sizeof(int *));
        input->nwavelengths_3d = 0;

        // Iterate through iwave
        int size_check = 0;
        for (int i = 0; i < (int)iwave.size(); i++) {
            // check if iwave[i] is in out_iwave
            bool found = false;
            for(int j = 0; j < (int)out_iwave.size(); j++){
                if(iwave[i] == out_iwave[j]){
                    found = true;
                    // Update num of wavelength_3d
                    (input->nwavelengths_3d)++;
                    size_check++;
                    break;
                }
            }
            if(found) {
                // Update wavelength_3d indexes
                input->wavelength_3d_index = (int *)realloc(input->wavelength_3d_index, (n + 1) * sizeof(int*));
                input->wavelength_3d_index[n] = i;
                n++;
            }        
        }
    //--------------------------------------------------------------------------------------------------------
        // Populate expanded user wavelength input into wavelength_3d array
        input->wavelength_3d = (int *)malloc((input->nwavelengths_3d) * sizeof(int *));
        for (int i = 0; i < input->nwavelengths_3d; i ++) {
            input->wavelength_3d[i] = l1file->iwave[input->wavelength_3d_index[i]];
        }
    //--------------------------------------------------------------------------------------------------------
    }
}

