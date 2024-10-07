// procedure to perform quality control on a VIIRS Level-2 file
// includes: 
// 	- metadata format check
//	- start/end time output

//	Arguments
//     
//     	Name    Type 	I/O 	Description
//     	----	---- 	--- 	-----------
//      ifile   string    I     VIIRS packet file name
//		cfile	string	  I		metadata configuration file (template)
//		ofile	string	  O		quality check output file

// Liang Hong, Dec. 8, 2017   V0.2
// Liang Hong, Sep. 21, 2018	V0.21
// Sean Bailey, Nov. 1, 2019  V0.22 - support new filename convention
// Liang Hong, Oct. 31, 2022  V0.33, fixed a bug in checking and displaying missing global att.

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <libgen.h>
#include <string.h>
#include <string>
#include <regex>
#include <math.h>
#include "netcdf.h"
#include <stdbool.h>
#include "l2qc_viirs.h"

#define VERSION "0.33"

#define MAX_ATT_NUM 500
#define MAX_ATT_LEN 256
#define MAX_ERR_LEN 100

using namespace std;

clo_optionList_t* optionList = NULL;


// The main function that checks if two given strings match. The first
// string may contain wildcard characters

bool match(char *first, char * second) {
    // If we reach at the end of both strings, we are done
    if (*first == '\0' && *second == '\0')
        return true;

    // Make sure that the characters after '*' are present in second string.
    // This function assumes that the first string will not contain two
    // consecutive '*' 
    if (*first == '*' && *(first + 1) != '\0' && *second == '\0')
        return false;

    // If the first string contains '?', or current characters of both 
    // strings match
    if (*first == '?' || *first == *second)
        return match(first + 1, second + 1);

    // If there is *, then there are two possibilities
    // a) We consider current character of second string
    // b) We ignore current character of second string.
    if (*first == '*')
        return match(first + 1, second) || match(first, second + 1);
    return false;
}

static void replacestr(char *line, const char *search, const char *replace) {
    char *sp;

    if ((sp = strstr(line, search)) == NULL) {
        return;
    }
    int search_len = strlen(search);
    int replace_len = strlen(replace);
    int tail_len = strlen(sp + search_len);

    memmove(sp + replace_len, sp + search_len, tail_len + 1);
    memcpy(sp, replace, replace_len);

    replacestr(line, search, replace);
}

int ValidBoundary(const char *strMinMax, float *datamin, float *datamax) {
    const char *ptr;
    char strBoundary[100];
    int index;
    int ch = '/';

    // find the min/max value separator "/"
    if ((ptr = strchr(strMinMax, ch)) == NULL) {
        return false;
    }
    index = ptr - strMinMax;

    // find lower boundary
    memset(strBoundary, '\0', sizeof (strBoundary));
    strncpy(strBoundary, strMinMax, index);
    *datamin = atof(strBoundary);


    // find upper boundary
    memset(strBoundary, '\0', sizeof (strBoundary));
    strcpy(strBoundary, ptr + 1);
    *datamax = atof(strBoundary);

    return true;
}

int main(int argc, char *argv[]) {
    printf("l2qc_viirs Version: %s (%s %s)\n", VERSION, __DATE__, __TIME__);

    char *str_L2file = NULL; // VIIRS L2 file
    char *str_rfile = NULL; // report text file
    char *str_cfile = NULL; // configuration text file
    int bReportFile = 0;
    int i, j, iMetaTplt;
    int status; /* error status */
    int ncid; /* netCDF ID */
    int natts; // # of attributes found in Meta data
    int grp0id, grp1id;
    char att_text[MAX_ATT_LEN + 1];
    float att_float;
    double att_double;
    long att_long;
    char attname[NC_MAX_NAME + 1];
    float datamin, datamax;
    nc_type xtypep;
    char strTimeStamp[16] = "[l2timestring]";
    char strTimeString0[15] = "[timestring0]Z";
    char strL2ProdString[12] = "[L2Product]";
    int nMissingMetaData = 0;
    int nWrongMetaData = 0;
    int nErr = 0;
    int nFlagged = 0;
    char strFlagged[1000] = "L2Flags_set=";
    int nattsp;
    int il2_Tplt_global_att = 0;			// Template global attribute index
    int il2_Tplt_processing_control = 0;	// Template processing control index

    char *l2_Tplt_global_att_name[MAX_ATT_NUM];
    char *l2_Tplt_global_att_val[MAX_ATT_NUM];
    char *l2_Tplt_processing_control_g1[MAX_ATT_NUM];		// processing control group 1
    char *l2_Tplt_processing_control_name[MAX_ATT_NUM];		// processing control attribute name
    char *l2_Tplt_processing_control_val[MAX_ATT_NUM];
    

    char *errstr = (char*) malloc(10000 * MAX_ERR_LEN); // string array to store errors
    char strTplt2replace[2] = "|";			// server database won't take '|', using ',' instead
    char strDB2accept[2] = ",";

    optionList = clo_createList();

    l2qcviirs_init_options(optionList, VERSION);
    if (argc == 1) {
        clo_printUsage(optionList);
        exit(1);
    }
    l2qcviirs_read_options(optionList, argc, argv);

    str_L2file = clo_getString(optionList, "ifile");
    str_cfile = clo_getString(optionList, "cfile");
    str_rfile = clo_getString(optionList, "ofile");
    if (str_rfile != NULL) {
        bReportFile = 1;
    }
     
    string str_L2basename = basename(str_L2file);
    char str_L2_product[10] = "All";

    regex OC ("(.*)(OC)(.*)");
    regex IOP ("(.*)(IOP)(.*)");
    regex SST ("(.*)(SST)(.*)");
    regex SST3 ("(.*)(SST3)(.*)");
    if (regex_match(str_L2basename,OC)){
    	strcpy(str_L2_product, "OC");
    } else if (regex_match(str_L2basename,IOP)) {
    	strcpy(str_L2_product, "IOP");
    } else if (regex_match(str_L2basename,SST3)) {
    	strcpy(str_L2_product, "SST3");
    } else if (regex_match(str_L2basename,SST)) {
    	strcpy(str_L2_product, "SST");
    } 
	
    FILE *outfile;

    // read configuration file
    FILE *configfile = fopen(str_cfile, "r");
    if (!configfile) {
        printf("Could not open configuration file. Exiting application.\n");
        return 1;
    }

    // read in the whole configuration text file as a string
    char *file_contents;
    long input_file_size;
    fseek(configfile, 0, SEEK_END);
    input_file_size = ftell(configfile);
    rewind(configfile);
    file_contents = (char*) calloc(input_file_size+50,1);
    fread(file_contents, sizeof (char), input_file_size, configfile);
    fclose(configfile);

    // parse the configuration file content string
    char *line;
    line = strtok(file_contents, "\n");

    int nAttributes_global=0;
    int nAttributes_processing_control=0;

    while (line != NULL) {
        if (match((char*) "#*", line)) {
            // comment line
            line = strtok(NULL, "\n");
            continue;
        } else if (match((char*) "[*]processing_control*", line)) {
        	// processing control group
        	// e.g. [OC,IOP,SST,SST3]processing_control.flag_percentages.NAVFAIL=0/0
        	
        	// if currently opened L2 product is in [], store attribute; otherwise, skip this line
        	char *tmpstrRBracket; 		// find the first right bracket in line
        	char *tmpstrL2product;		// find first L2 product info in the list
        	char *tmpstrEqual;			// find the first equal sign
        	tmpstrRBracket = strstr(line,"]");
        	tmpstrL2product = strstr(line,str_L2_product);
			tmpstrEqual = strstr(line,"=");
        	if (tmpstrL2product!=NULL) {
				if (tmpstrL2product<tmpstrRBracket) {
					line = tmpstrRBracket + 1; 
					replacestr(line, "processing_control.", "");	// line = flag_percentages.NAVFAIL=0/0
					char *strDot; // find first "." in the line
            		strDot = strstr(line, ".");					
					l2_Tplt_processing_control_g1[il2_Tplt_processing_control] = (char*) calloc(MAX_ATT_LEN,1);	
					strncpy(l2_Tplt_processing_control_g1[il2_Tplt_processing_control], line, strDot-tmpstrRBracket-1);
					tmpstrEqual = strstr(line,"=");
					l2_Tplt_processing_control_name[il2_Tplt_processing_control] = (char*) calloc(MAX_ATT_LEN,1);					
					strncpy(l2_Tplt_processing_control_name[il2_Tplt_processing_control], strDot+1, tmpstrEqual-strDot-1);
					l2_Tplt_processing_control_val[il2_Tplt_processing_control] = (char*) calloc(MAX_ATT_LEN,1);					
					strcpy(l2_Tplt_processing_control_val[il2_Tplt_processing_control], tmpstrEqual+1);		
					il2_Tplt_processing_control++;
                                        line = NULL;
				}
        	}

        } else {
        	// global attribute
        	// e.g. [OC,IOP,SST,SST3]orbit_number=1/4294967295
        	
        	char *tmpstrRBracket; 		// find the first right bracket in line
        	char *tmpstrL2product;		// find first L2 product info in the list
        	char *tmpstrEqual;			// find the first equal sign
        	tmpstrRBracket = strstr(line,"]");
        	tmpstrL2product = strstr(line,str_L2_product);
			tmpstrEqual = strstr(line,"=");
        	if (tmpstrL2product!=NULL) {
				if (tmpstrL2product<tmpstrRBracket) {
					l2_Tplt_global_att_name[il2_Tplt_global_att] = (char*) calloc(MAX_ATT_LEN,1);					
					strncpy(l2_Tplt_global_att_name[il2_Tplt_global_att], tmpstrRBracket+1, tmpstrEqual-tmpstrRBracket-1);
					l2_Tplt_global_att_val[il2_Tplt_global_att] = (char*) calloc(MAX_ATT_LEN,1);					
					strcpy(l2_Tplt_global_att_val[il2_Tplt_global_att], tmpstrEqual+1);		
					il2_Tplt_global_att++;
				}
        	}
        }
        if (line){
            line = strtok(NULL, "\n");
        }
    }
    
    // update number of all attributes to check in the tempalte that applies to this L2 Product
	nAttributes_global = il2_Tplt_global_att;
	nAttributes_processing_control = il2_Tplt_processing_control;

    // open L2 netcdf file
    status = nc_open(str_L2file, NC_NOWRITE, &ncid);
    if (status != NC_NOERR) {
        printf("Unable to open input VIIRS L2 file.\n");
        return 1;
    }
	
    // open report file
    if (bReportFile) {
        outfile = fopen(str_rfile, "a+");
        if (outfile == NULL) {
            printf("Unable to open output report file.\n");
            return 1;
        }
        fprintf(outfile, "basename=%s\n", str_L2basename.c_str());
    }

    // Check meta data format
    status = nc_inq_natts(ncid, &natts);

    int meta_global_idx[nAttributes_global];				// global attributes
    for (i = 0; i < nAttributes_global; i++) {
        // for each item in metadata template
        // default = 0: missing metadata in file
        // -1: found matching item in file
        // 1: found item in file, but not matching
        // 2: wrong format or data range in file
        meta_global_idx[i] = 0;
    }
 
//    int meta_processing_control_idx[nAttributes_processing_control];				// processing control group
//    for (i = 0; i < nAttributes_processing_control; i++) {
//        // for each item in metadata template
//        // default = 0: missing metadata in file
//        // -1: found matching item in file
//        // 1: found item in file, but not matching
//        // 2: wrong format or data range in file
//        meta_processing_control_idx[i] = 0;
//    }   
	
	// Check processing control groups
    status = nc_inq_ncid(ncid, "processing_control", &grp0id);
    for (int ipc = 0; ipc < nAttributes_processing_control; ipc++) {
        status = nc_inq_ncid(grp0id, l2_Tplt_processing_control_g1[ipc], &grp1id);
        status = nc_inq_natts(grp1id, &nattsp);
        for (i = 0; i < nattsp; i++) {
            status = nc_inq_attname(grp1id, NC_GLOBAL, i, attname);
            status = nc_get_att_float(grp1id, NC_GLOBAL, attname, &att_float);
            if (strcmp(attname, l2_Tplt_processing_control_name[ipc]) == 0) {
                if (ValidBoundary(l2_Tplt_processing_control_val[ipc], &datamin, &datamax)) {
                    if ((att_float < datamin) || (att_float > datamax)) {
                        if (nFlagged > 0)
                            strcat(strFlagged, ",");
                        strcat(strFlagged, attname);
                        nFlagged++;
                    }
                }
            }
        }
    }
	
    // read global attributes from L2 file
    for (i = 0; i < natts; i++) {
        // read each of the metadata attributes in L2 file and compare with the template
        // to check missing, un-matched and extra items.
        status = nc_inq_attname(ncid, NC_GLOBAL, i, attname);
        status = nc_inq_atttype(ncid, NC_GLOBAL, attname, &xtypep);
        // read metadata template value
        iMetaTplt = -1;
        char strValTplt[MAX_ATT_LEN + 1], strNameTplt[MAX_ATT_LEN + 1];
        
		// find metadata item in the template, output the index
		// otherwise, skip this attribute	
		for (j = 0; j < nAttributes_global; j++) {
			strcpy(strNameTplt, l2_Tplt_global_att_name[j]);
			if (match(strNameTplt, attname)) {
				iMetaTplt = j;
				break;
			}
		}
		if (iMetaTplt < 0) continue;
		
		strcpy(strValTplt, l2_Tplt_global_att_val[iMetaTplt]);
		
        switch (xtypep) {
        case NC_BYTE:
            //printf("NC_BYTE\n");
            break;
        case NC_CHAR:
            //printf("NC_CHAR\n");
            status = nc_get_att_text(ncid, NC_GLOBAL, attname, att_text);
            char strTmpstring[100];
            char *strDot; // find first "." in file name, e.g. V2015363144800_S20151229130521.L2_SNPP_OC.nc
            strDot = strstr(basename(str_L2file), ".");
            strncpy(strTmpstring, basename(str_L2file)+1,strDot-basename(str_L2file)-1);
            replacestr(strValTplt, strTimeStamp, strTmpstring);
            // replace [L2Product] with OC or IOP...
            replacestr(strValTplt, strL2ProdString, str_L2_product);	
            if (match(strValTplt, strTimeString0)) {
                int iyr, mm, dd, ih, mn, sec, msec;
                bool bPassCheck = true;
                sscanf(att_text, "%d-%d-%dT%d:%d:%d.%dZ", &iyr, &mm, &dd, &ih, &mn, &sec, &msec);
                if ((iyr < 2000) || (iyr > 3000)) bPassCheck = false;
                if ((mm < 1) || (mm > 12)) bPassCheck = false;
                if ((dd < 1) || (dd > 31)) bPassCheck = false;
                if ((ih < 0) || (ih > 24)) bPassCheck = false;
                if ((mn < 0) || (mn > 60)) bPassCheck = false;
                if ((sec < 0) || (sec > 60)) bPassCheck = false;
                if ((msec < 0) || (msec > 999)) bPassCheck = false;
                if (!bPassCheck) {
                    meta_global_idx[iMetaTplt] = 2;
                    nWrongMetaData++;
                    char tmpstr[500];
                    sprintf(tmpstr, "var[%d]: %s = %s [Invalid time string]", i, attname, att_text);
                    replacestr(tmpstr, strTplt2replace, strDB2accept);					// replace "|" with "," for database table "file_data" compatibility
                    strncpy(&errstr[nErr * MAX_ERR_LEN], tmpstr, MAX_ERR_LEN);
                    printf("%s\n", tmpstr);
                    nErr++;
                    break;
                }
            } else if (strstr(strValTplt, "|") != NULL) {
                if (strstr(strValTplt, att_text) == NULL) {
                    meta_global_idx[iMetaTplt] = 1;
                    nWrongMetaData++;
                    char tmpstr[500];
                    sprintf(tmpstr, "%s not matching: Template val = %s; att_text = %s", attname, strValTplt, att_text);
            		replacestr(tmpstr, strTplt2replace, strDB2accept);					// replace "|" with "," for database table "file_data" compatibility 
                    strncpy(&errstr[nErr * MAX_ERR_LEN], tmpstr, MAX_ERR_LEN);
                    printf("%s\n", tmpstr);
                    nErr++;
                    break;
                }
            } else if (!match(strValTplt, att_text)) {
                meta_global_idx[iMetaTplt] = 1;
                nWrongMetaData++;
                char tmpstr[500];
                sprintf(tmpstr, "%s not matching: Template val = %s; att_text = %s", attname, strValTplt, att_text);
           		replacestr(tmpstr, strTplt2replace, strDB2accept);					// replace "|" with "," for database table "file_data" compatibility                
                strncpy(&errstr[nErr * MAX_ERR_LEN], tmpstr, MAX_ERR_LEN);
                printf("%s\n", tmpstr);
                nErr++;
                break;
            }

            meta_global_idx[iMetaTplt] = -1; // means passing check
            break;
        case NC_SHORT:
            //printf("NC_SHORT\n");
            break;
        case NC_LONG:
            //printf("NC_LONG\n");
            status = nc_get_att_long(ncid, NC_GLOBAL, attname, &att_long);

            // Check attribute data value to see if it is within boundary
            if (ValidBoundary(strValTplt, &datamin, &datamax)) {
                if ((att_long >= datamin) && (att_long <= datamax)) {
                    meta_global_idx[iMetaTplt] = -1;
                } else {
                    meta_global_idx[iMetaTplt] = 2;
                    nWrongMetaData++;
                    char tmpstr[500];
                    sprintf(tmpstr, "%s: data not in correct range", attname);
            		replacestr(tmpstr, strTplt2replace, strDB2accept);					// replace "|" with "," for database table "file_data" compatibility                    
                    strncpy(&errstr[nErr * MAX_ERR_LEN], tmpstr, MAX_ERR_LEN);
                    printf("%s\n", tmpstr);
                    nErr++;
                    break;
                }
            } else
                meta_global_idx[iMetaTplt] = -1;

            break;
        case NC_FLOAT:
            //printf("NC_FLOAT\n");
            status = nc_get_att_float(ncid, NC_GLOBAL, attname, &att_float);

            // Check attribute data value to see if it is within boundary
            if (ValidBoundary(strValTplt, &datamin, &datamax)) {
                if ((att_float >= datamin) && (att_float <= datamax)) {
                    meta_global_idx[iMetaTplt] = -1;
                } else {
                    meta_global_idx[iMetaTplt] = 2;
                    nWrongMetaData++;
                    char tmpstr[500];
                    sprintf(tmpstr, "%s: data not in correct range", attname);
            		replacestr(tmpstr, strTplt2replace, strDB2accept);					// replace "|" with "," for database table "file_data" compatibility
                    strncpy(&errstr[nErr * MAX_ERR_LEN], tmpstr, MAX_ERR_LEN);
                    printf("%s\n", tmpstr);
                    nErr++;
                    break;
                }
            } else
                meta_global_idx[iMetaTplt] = -1;

            break;
        case NC_DOUBLE:
            //printf("NC_DOUBLE\n");
            status = nc_get_att_double(ncid, NC_GLOBAL, attname, &att_double);

            // Check attribute data value to see if it is within boundary
            if (ValidBoundary(strValTplt, &datamin, &datamax)) {
                if ((att_double >= datamin) && (att_double <= datamax)) {
                    meta_global_idx[iMetaTplt] = -1;
                } else {
                    meta_global_idx[iMetaTplt] = 2;
                    nWrongMetaData++;
                    char tmpstr[500];
                    sprintf(tmpstr, "%s: data not in correct range", attname);
                    replacestr(tmpstr, strTplt2replace, strDB2accept);					// replace "|" with "," for database table "file_data" compatibility
                    strncpy(&errstr[nErr * MAX_ERR_LEN], tmpstr, MAX_ERR_LEN);
                    printf("%s\n", tmpstr);
                    nErr++;
                    break;
                }
            } else
                meta_global_idx[iMetaTplt] = -1;
            break;
        case NC_STRING:
            //printf("NC_STRING\n");
            //status = nc_get_att_text(ncid, NC_GLOBAL, attname, att_text);
            break;
        default:
            char tmpstr[500];
            sprintf(tmpstr, "Attribute [%d] data type not found in template!", i);
            replacestr(tmpstr, strTplt2replace, strDB2accept);					// replace "|" with "," for database table "file_data" compatibility
            strncpy(&errstr[nErr * MAX_ERR_LEN], tmpstr, MAX_ERR_LEN);
            printf("%s\n", tmpstr);
            nErr++;
        }
		
        if (meta_global_idx[iMetaTplt] != -1) printf("not passing check!\n");

        if (strcmp(attname, "time_coverage_start") == 0) {
            printf("start_time=%s\n", att_text);
            if (bReportFile) fprintf(outfile, "start_time=%s\n", att_text);
        }
        if (strcmp(attname, "time_coverage_end") == 0) {
            printf("stop_time=%s\n", att_text);
            if (bReportFile) fprintf(outfile, "stop_time=%s\n", att_text);
        }
        if (strcmp(attname, "processing_version") == 0) {
            printf("processing_version=%s\n", att_text);
            if (bReportFile) fprintf(outfile, "processing_version=%s\n", att_text);
        }
        if (strcmp(attname, "day_night_flag") == 0) {
            printf("day_night_flag=%s\n", att_text);
            if (bReportFile) fprintf(outfile, "day_night_flag=%s\n", att_text);
        }
    }
	
    for (i = 0; i < nAttributes_global; i++) {
        if (meta_global_idx[i] == 0) {
            nMissingMetaData++;
            printf("L2 file missing att[%d] %s\n", i, l2_Tplt_global_att_name[i]);
        }
    }

    if ((nMissingMetaData == 0) && (nWrongMetaData == 0)) printf("No error found in Metadata.\n");
	
	if (nFlagged < 1)
        printf("No flags for attention in Metadata.\n");
    else if (nFlagged == 1)
        printf("1 flag for attention in Metadata.\n");
    else
        printf("%d flags for attention in Metadata.\n", nFlagged);

    // output report file     
    if (bReportFile) {
        fprintf(outfile, "errors=%d\n", nErr);
        for (i = 0; i < nErr; i++) {
            fprintf(outfile, "err_%d=%.*s\n", i + 1, MAX_ERR_LEN, errstr + i * MAX_ERR_LEN);
        }
        if (nFlagged > 0) fprintf(outfile, "%s\n", strFlagged);
        fclose(outfile);
    }

    status = nc_close(ncid);

    free(errstr);
    for (i = 0; i < nAttributes_global; i++) free(l2_Tplt_global_att_name[i]);
    for (i = 0; i < nAttributes_global; i++) free(l2_Tplt_global_att_val[i]);
    for (i = 0; i < nAttributes_processing_control; i++) free(l2_Tplt_processing_control_name[i]);
    for (i = 0; i < nAttributes_processing_control; i++) free(l2_Tplt_processing_control_val[i]);

    clo_deleteList(optionList);
    return nMissingMetaData+nWrongMetaData+nFlagged+nErr;
}



