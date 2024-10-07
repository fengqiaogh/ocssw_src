#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

#include "hawkeyeUtil.h"


//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     FutureTech     05/21/18 0.10  Original development
//  Joel Gales     FutureTech     09/20/18 0.11  Exit if EOF when searching
//                                               for 2052/2053 packet
//  Joel Gales     SAIC           11/01/18       Define ifct as uint8_t
//                                               and compute uint8_t byte_diff
//                                               = frame[2]-ifct to fix prob
//                                               with dropped frames
//  Joel Gales     SAIC           11/16/18       Fixed corrupted code in
//                                               get_packet_from_frame() that
//                                               broke l0gen_hawkeye
//  Liang Hong     SAIC           07/13/21       Update in handling EOI packet
//                                               crossing frames
//  Liang Hong     SAIC           01/12/22 0.20  added scale and offset attributes

using namespace std;

/*----------------------------------------------------------------- */
/* Create an Generic NETCDF4 file                                   */
/* ---------------------------------------------------------------- */
int createFile( const char* filename, const char* cdlfile,
                size_t sdim, int *ncid, int *gid) {
  
   int status;
   
   status = nc_create( filename, NC_NETCDF4, ncid);
   check_err(status,__LINE__,__FILE__);

   ifstream hawkeye_data_structure;
   string line;
   string dataStructureFile;

   dataStructureFile.assign( cdlfile);
   expandEnvVar( &dataStructureFile);

   hawkeye_data_structure.open( dataStructureFile.c_str(), ifstream::in);
   if ( hawkeye_data_structure.fail() == true) {
     cout << "\"" << dataStructureFile.c_str() << "\" not found" << endl;
     exit(1);
   }

   // Find "dimensions" section of CDL file
   while(1) {
     getline( hawkeye_data_structure, line);
     size_t pos = line.find("dimensions:");
     if ( pos == 0) break;
   }

   // Define dimensions from "dimensions" section of CDL file
   int ndims = 0;
   int dimid[1000];
   while(1) {
     getline( hawkeye_data_structure, line);
     size_t pos = line.find(" = ");
     if ( pos == string::npos) break;

     uint32_t dimSize;
     istringstream iss(line.substr(pos+2, string::npos));
     iss >> dimSize;

     iss.clear(); 
     iss.str( line);
     iss >> skipws >> line;

     cout << "Dimension Name: " << line.c_str() << " Dimension Size: "
    	  << dimSize << endl;

     if (line.compare("number_of_scans") == 0) {
       dimSize = sdim;
     }

     status = nc_def_dim( *ncid, line.c_str(), dimSize, &dimid[ndims++]);
     check_err(status,__LINE__,__FILE__);
     
   } // while loop

   // Read global attributes (string attributes only) 
   while(1) {
     getline( hawkeye_data_structure, line);
     size_t pos = line.find("// global attributes");
     if ( pos == 0) break;
   }

   while(1) {
     getline( hawkeye_data_structure, line);
     size_t pos = line.find(" = ");
     if ( pos == string::npos) break;

     //size_t slen = line.length();

     string attValue;

     // Remove leading and trailing quotes
     attValue.assign(line.substr(pos+4));
     size_t posQuote = attValue.find("\"");
     attValue.assign(attValue.substr(0, posQuote));

     istringstream iss(line.substr(pos+2));
     iss.clear(); 
     iss.str( line);
     iss >> skipws >> line;

     // Skip commented out attributes
     if (line.compare("//") == 0) continue;

     string attName;
     attName.assign(line.substr(1).c_str());

     // Skip non-string attributes
     if (attName.compare("orbit_number") == 0) continue;
     if (attName.compare("history") == 0) continue;
     if (attName.compare("format_version") == 0) continue;
     if (attName.compare("instrument_number") == 0) continue;
     if (attName.compare("pixel_offset") == 0) continue;
     if (attName.compare("number_of_filled_scans") == 0) continue;

     // cout << attName.c_str() << " " << attValue.c_str() << endl;

     status = nc_put_att_text( *ncid, NC_GLOBAL, attName.c_str(), 
                              strlen(attValue.c_str()), attValue.c_str());
     check_err(status,__LINE__,__FILE__);

   } // while(1)


   int ngrps = 0;
   // Loop through groups
   while(1) {
     getline( hawkeye_data_structure, line);

     // Check if end of CDL file
     // If so then close CDL file and return
     if (line.substr(0,1).compare("}") == 0) {
       hawkeye_data_structure.close();
       return 0;
     }

     // Check for beginning of new group
     size_t pos = line.find("group:");

     // If found then create new group and variables
     if ( pos == 0) {

       // Parse group name
       istringstream iss(line.substr(6, string::npos));
       iss >> skipws >> line;

       // Create NCDF4 group
       status = nc_def_grp( *ncid, line.c_str(), &gid[ngrps]);
       check_err(status,__LINE__,__FILE__);

       ngrps++;

       int numDims=0;
       int varDims[NC_MAX_DIMS];
       size_t dimSize[NC_MAX_DIMS];
       char dimName[NC_MAX_NAME+1];
       string sname;
       string lname;
       string standard_name;
       string units;
       string flag_values;
       string flag_meanings;
       double valid_min=0.0;
       double valid_max=0.0;
       double fill_value=0.0;
       float scale_factor=1.0;
       float add_offset=0.0;

       int ntype=0;

       // Loop through datasets in group
       getline( hawkeye_data_structure, line); // skip "variables:"
       while(1) {
         getline( hawkeye_data_structure, line);

         if (line.length() == 0) continue;
         if (line.substr(0,1).compare("\r") == 0) continue;
         if (line.substr(0,1).compare("\n") == 0) continue;

         size_t pos = line.find(":");

         // No ":" found, new dataset or empty line or end-of-group
         if ( pos == string::npos) {

           if ( numDims > 0) {
             // Create previous dataset
             createNCDF( gid[ngrps-1], 
                         sname.c_str(), lname.c_str(),
                         standard_name.c_str(), units.c_str(),
                         (void *) &fill_value, 
                         flag_values.c_str(), flag_meanings.c_str(),
                         valid_min, valid_max, scale_factor, add_offset, ntype, numDims, varDims);

             flag_values.assign("");
             flag_meanings.assign("");
             units.assign("");
           }

           valid_min=0.0;
           valid_max=0.0;
           fill_value=0.0;
           scale_factor=1.0;
           add_offset=0.0;

           if (line.substr(0,10).compare("} // group") == 0) break;

           // Parse variable type
           string varType;
           istringstream iss(line);
           iss >> skipws >> varType;

           // Get corresponding NC variable type
           if ( varType.compare("char") == 0) ntype = NC_CHAR;
           else if ( varType.compare("byte") == 0) ntype = NC_BYTE;
           else if ( varType.compare("short") == 0) ntype = NC_SHORT;
           else if ( varType.compare("int") == 0) ntype = NC_INT;
           else if ( varType.compare("long") == 0) ntype = NC_INT;
           else if ( varType.compare("float") == 0) ntype = NC_FLOAT;
           else if ( varType.compare("real") == 0) ntype = NC_FLOAT;
           else if ( varType.compare("double") == 0) ntype = NC_DOUBLE;
           else if ( varType.compare("ubyte") == 0) ntype = NC_UBYTE;
           else if ( varType.compare("ushort") == 0) ntype = NC_USHORT;
           else if ( varType.compare("uint") == 0) ntype = NC_UINT;
           else if ( varType.compare("int64") == 0) ntype = NC_INT64;
           else if ( varType.compare("uint64") == 0) ntype = NC_UINT64;

           // Parse short name (sname)
           pos = line.find("(");
           size_t posSname = line.substr(0, pos).rfind(" ");
           sname.assign(line.substr(posSname+1, pos-posSname-1));
           cout << "sname: " << sname.c_str() << endl;

           // Parse variable dimension info
           parseDims( *ncid, ndims, line.substr(pos+1, string::npos),
                      &numDims, dimid, varDims);
           for (int i=0; i<numDims; i++) { 
             nc_inq_dim( *ncid, varDims[i], dimName, &dimSize[i]);
             cout << line.c_str() << " " << i << " " << dimName
             	  << " " << dimSize[i] << endl;
           }

         } else {
           // Parse variable attributes
           size_t posEql = line.find("=");
           size_t pos1qte = line.find("\"");
           size_t pos2qte = line.substr(pos1qte+1, string::npos).find("\"");
	   // cout << line.substr(pos+1, posEql-pos-2).c_str() << endl;

           string attrName = line.substr(pos+1, posEql-pos-2);

           // Get long_name
           if ( attrName.compare("long_name") == 0) {
             lname.assign(line.substr(pos1qte+1, pos2qte));
             //             cout << "lname: " << lname.c_str() << endl;
           }

           // Get units
           else if ( attrName.compare("units") == 0) {
             units.assign(line.substr(pos1qte+1, pos2qte));
             //             cout << "units: " << units.c_str() << endl;
           }

           // Get _FillValue
           else if ( attrName.compare("_FillValue") == 0) {
             iss.clear(); 
             iss.str( line.substr(posEql+1, string::npos));
             iss >> fill_value;
             //             cout << "_FillValue: " << fill_value << endl;
           }

           // Get flag_values
           else if ( attrName.compare("flag_values") == 0) {
             flag_values.assign(line.substr(pos1qte+1, pos2qte));
           }

           // Get flag_meanings
           else if ( attrName.compare("flag_meanings") == 0) {
             flag_meanings.assign(line.substr(pos1qte+1, pos2qte));
           }

           // Get valid_min
           else if ( attrName.compare("valid_min") == 0) {
             iss.clear(); 
             iss.str( line.substr(posEql+1, string::npos));
             iss >> valid_min;
             //             cout << "valid_min: " << valid_min << endl;
           }

           // Get valid_max
           else if ( attrName.compare("valid_max") == 0) {
             iss.clear(); 
             iss.str( line.substr(posEql+1, string::npos));
             iss >> valid_max;
             //             cout << "valid_max: " << valid_max << endl;
           }
           
           // Get scale_factor
           else if ( attrName.compare("scale_factor") == 0) {
             iss.clear(); 
             iss.str( line.substr(posEql+1, string::npos));
             iss >> scale_factor;
           }
           
           // Get add_offset
           else if ( attrName.compare("add_offset") == 0) {
             iss.clear(); 
             iss.str( line.substr(posEql+1, string::npos));
             iss >> add_offset;
           }

         } // if ( pos == string::npos)
       } // datasets in group loop
     } // New Group loop
   } // Main Group loop
   
   return 0;
}


int parseDims( int ncid, int ndims, string dimString,
               int *numDims, int *dimid, int *varDims) {

  size_t dimSize, curPos=0;
  char dimName[NC_MAX_NAME+1];

  *numDims = 0;

  while(1) {
    size_t pos = dimString.find(",", curPos);
    if ( pos == string::npos) 
      pos = dimString.find(")");

    string varDimName;
    istringstream iss(dimString.substr(curPos, pos-curPos));
    iss >> skipws >> varDimName;

    for (int i=0; i<ndims; i++) {
      int status = nc_inq_dim( ncid, dimid[i], dimName, &dimSize);
      check_err(status,__LINE__,__FILE__);
      if ( varDimName.compare(dimName) == 0) {
        varDims[(*numDims)++] = dimid[i];
        break;
      }
    }
    if ( dimString.substr(pos, 1).compare(")") == 0) break;

    curPos = pos + 1;
  }

  return 0;
}



int createNCDF( int ncid, const char *sname, const char *lname, 
                const char *standard_name, const char *units,
                void *fill_value,
                const char *flag_values, const char *flag_meanings,
                double low, double high, 
                float scale_factor, float add_offset, 
                int nt, int rank, int *dimids) {

  int32_t varid;
  int status;
  size_t dimlength;
  size_t chunksize[3];
     
  /* Create the NCDF dataset */
  status = nc_def_var(ncid, sname, nt, rank, dimids, &varid);
  if( status != NC_NOERR) {
    printf("-E- %s %d: %s for %s\n", 
	   __FILE__, __LINE__, nc_strerror(status), sname);
    exit(1);
  } 

  // Set fill value
  double fill_value_dbl;
  memcpy( &fill_value_dbl, fill_value, sizeof(double));

  int8_t i8;
  uint8_t ui8;
  int16_t i16;
  uint16_t ui16;
  int32_t i32;
  int32_t ui32;
  float f32;

  //  if ( (low < high) && (low != fill_value_dbl)) {
  if ( low != fill_value_dbl) {
    if ( nt == NC_BYTE) {
      i8 = fill_value_dbl;
      status = nc_def_var_fill( ncid, varid, 0, (void *) &i8);
    } else if ( nt == NC_UBYTE) {
      ui8 = fill_value_dbl;
      status = nc_def_var_fill( ncid, varid, 0, (void *) &ui8);
    } else if ( nt == NC_SHORT) {
      i16 = fill_value_dbl;
      status = nc_def_var_fill( ncid, varid, 0, (void *) &i16);
    } else if ( nt == NC_USHORT) {
      ui16 = fill_value_dbl;
      status = nc_def_var_fill( ncid, varid, 0, (void *) &ui16);
    } else if ( nt == NC_INT) {
      i32 = fill_value_dbl;
      status = nc_def_var_fill( ncid, varid, 0, (void *) &i32);
    } else if ( nt == NC_UINT) {
      ui32 = fill_value_dbl;
      status = nc_def_var_fill( ncid, varid, 0, (void *) &ui32);
    } else if ( nt == NC_FLOAT) {
      f32 = fill_value_dbl;
      status = nc_def_var_fill( ncid, varid, 0, (void *) &f32);
    } else {
      status = nc_def_var_fill( ncid, varid, 0, (void *) &fill_value_dbl);
    }
    check_err(status,__LINE__,__FILE__);
  }

  /* vary chunck size based on dimensions */ 
  int do_deflate = 0;
  if ( rank == 3 && (strncmp(sname, "EV_", 3) == 0)) {
    status = nc_inq_dimlen(ncid, dimids[2], &dimlength);
    chunksize[0] = 1;
    chunksize[1] = 16;
    chunksize[2] = dimlength/10;
    do_deflate = 1;
  }

  /* Set compression */
  if ( do_deflate) {
    /* First set chunking */
    status = nc_def_var_chunking(ncid, varid, NC_CHUNKED, chunksize);
    if (status != NC_NOERR) {
      printf("-E- %s %d: %s for %s\n", __FILE__, __LINE__,
             nc_strerror(status), sname);
      exit(1);
    }

    /* Now we can set compression */
    //    status = nc_def_var_deflate(ncid, varid, NC_NOSHUFFLE, 1, 5);
    status = nc_def_var_deflate(ncid, varid, NC_SHUFFLE, 1, 5);
    if (status != NC_NOERR) {
      printf("-E- %s %d: %s for %s\n", __FILE__, __LINE__,
             nc_strerror(status), sname);
      exit(1);
    }
  }


  /* Add a "long_name" attribute */
  status = nc_put_att_text(ncid, varid, "long_name", strlen(lname), lname);
  if( status != NC_NOERR) {
    printf("-E- %s %d: %s for %s\n", 
	   __FILE__, __LINE__, nc_strerror(status), "long_name");
    exit(1);
  } 

  /* Add a "flag_values" attribute if specified*/
  // Parse string and save as signed bytes
  if ( strcmp( flag_values, "") != 0) {

    size_t curPos=0;

    string fv;
    fv.assign( flag_values);
    size_t pos = fv.find("=", curPos);
    fv = fv.substr(pos+1);

    size_t semicln = fv.find(";");
    pos = 0;

    int8_t vec[1024];
    int n = 0;
    while(pos != semicln) {
      pos = fv.find(",", curPos);
      if ( pos == string::npos) 
        pos = semicln;

      string flag_value;
      istringstream iss(fv.substr(curPos, pos-curPos));
      iss >> skipws >> flag_value;
      vec[n++] = atoi( flag_value.c_str());
      curPos = pos + 1;
    }

    status = nc_put_att_schar(ncid, varid, "flag_values", NC_BYTE, n, vec);
    if( status != NC_NOERR) {
      printf("-E- %s %d: %s for %s\n", 
             __FILE__, __LINE__, nc_strerror(status), "flag_values");
      exit(1);
    } 
  }

  /* Add a "flag_meanings" attribute if specified*/
  if ( strcmp( flag_meanings, "") != 0) {
    status = nc_put_att_text(ncid, varid, "flag_meanings", 
                             strlen(flag_meanings), flag_meanings);
    if( status != NC_NOERR) {
      printf("-E- %s %d: %s for %s\n", 
             __FILE__, __LINE__, nc_strerror(status), "flag_meanings");
      exit(1);
    } 
  }
  
  /* Add scale_factor and add_offset*/
  if ((scale_factor!=1) || (add_offset!=0) ) {
  	status = nc_put_att_float(ncid, varid,"scale_factor",NC_FLOAT,1,&scale_factor);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "scale_factor");
	  exit(1);
	} 
	status = nc_put_att_float(ncid, varid,"add_offset",NC_FLOAT,1,&add_offset);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "add_offset");
	  exit(1);
	} 
  }
  
  /* Add "valid_min/max" attributes if specified */
  if (low < high) {
    switch(nt) {              /* Use the appropriate number type */
    case NC_BYTE:
      {
	uint8_t vr[2];
	vr[0] = (uint8_t)low;
	vr[1] = (uint8_t)high;
	status = nc_put_att_uchar(ncid, varid,"valid_min",NC_BYTE,1,&vr[0]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_min");
	  exit(1);
	} 
	status = nc_put_att_uchar(ncid, varid,"valid_max",NC_BYTE,1,&vr[1]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_max");
	  exit(1);
	} 
      }
      break;
    case NC_UBYTE:
      {
	uint8_t vr[2];
	vr[0] = (uint8_t)low;
	vr[1] = (uint8_t)high;
	status = nc_put_att_uchar(ncid, varid,"valid_min",NC_UBYTE,1,&vr[0]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_min");
	  exit(1);
	} 
	status = nc_put_att_uchar(ncid, varid,"valid_max",NC_UBYTE,1,&vr[1]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_max");
	  exit(1);
	} 
      }
      break;
    case NC_SHORT:
      {
	int16_t vr[2];
	vr[0] = (int16_t)low;
	vr[1] = (int16_t)high;
	status = nc_put_att_short(ncid, varid,"valid_min",NC_SHORT,1,&vr[0]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_min");
	  exit(1);
	} 
	status = nc_put_att_short(ncid, varid,"valid_max",NC_SHORT,1,&vr[1]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_max");
	  exit(1);
	} 
      }
      break;
    case NC_USHORT:
      {
	uint16_t vr[2];
	vr[0] = (uint16_t)low;
	vr[1] = (uint16_t)high;
	status = nc_put_att_ushort(ncid, varid,"valid_min",NC_USHORT,1,&vr[0]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_min");
	  exit(1);
	} 
	status = nc_put_att_ushort(ncid, varid,"valid_max",NC_USHORT,1,&vr[1]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_max");
	  exit(1);
	} 
      }
      break;
    case NC_INT:
      {
	int32_t vr[2];
	vr[0] = (int32_t)low;
	vr[1] = (int32_t)high;
	status = nc_put_att_int(ncid, varid,"valid_min",NC_INT,1,&vr[0]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_min");
	  exit(1);
	} 
	status = nc_put_att_int(ncid, varid,"valid_max",NC_INT,1,&vr[1]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_max");
	  exit(1);
	} 
      }
      break;
    case NC_UINT:
      {
	uint32_t vr[2];
	vr[0] = (uint32_t)low;
	vr[1] = (uint32_t)high;
	status = nc_put_att_uint(ncid, varid,"valid_min",NC_UINT,1,&vr[0]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_min");
	  exit(1);
	} 
	status = nc_put_att_uint(ncid, varid,"valid_max",NC_UINT,1,&vr[1]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_max");
	  exit(1);
	} 
      }
      break;
    case NC_FLOAT:
      {
	float vr[2];
	vr[0] = (float)low;
	vr[1] = (float)high;
	status = nc_put_att_float(ncid, varid,"valid_min",NC_FLOAT,1,&vr[0]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_min");
	  exit(1);
	} 
	status = nc_put_att_float(ncid, varid,"valid_max",NC_FLOAT,1,&vr[1]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_max");
	  exit(1);
	} 
      }
      break;
    case NC_DOUBLE:
      {
	double vr[2];
	vr[0] = low;
	vr[1] = high;
	status = nc_put_att_double(ncid, varid,"valid_min",NC_DOUBLE,1,&vr[0]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_min");
	  exit(1);
	} 
	status = nc_put_att_double(ncid, varid,"valid_max",NC_DOUBLE,1,&vr[1]);
	if( status != NC_NOERR) {
	  printf("-E- %s %d: %s for %s\n", 
		 __FILE__, __LINE__, nc_strerror(status), "valid_max");
	  exit(1);
	} 
      }
      break;
    default:
      fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
      fprintf(stderr,"Got unsupported number type (%d) ",nt);
      fprintf(stderr,"while trying to create NCDF variable, \"%s\", ",sname);
      return(1);
    }
  }           
    
  /* Add a "units" attribute if one is specified */
  if(units != NULL && *units != 0) {
    status = nc_put_att_text(ncid, varid, "units", strlen(units), units);
    if( status != NC_NOERR) {
      printf("-E- %s %d: %s for %s\n", 
	     __FILE__, __LINE__, nc_strerror(status), "units");
      exit(1);
    } 
  }

  /* Add a "standard_name" attribute if one is specified */
  if(standard_name != NULL && *standard_name != 0) {
    status = nc_put_att_text(ncid, varid, "standard_name", 
			     strlen(standard_name), standard_name);
    if( status != NC_NOERR) {
      printf("-E- %s %d: %s for %s\n", 
	     __FILE__, __LINE__, nc_strerror(status), "standard_name");
      exit(1);
    } 
  }
  
  return 0;
}

int tepoch2yds( double tepoch, int32_t *iyr, int32_t *idy, double *sec) {
  
  // Program to convert Hawkeye epoch time to year, day, seconds

  // Correct for leap seconds
  int leap = 0;

  if (tepoch >= 15638401) leap = leap + 1; // Leapsecond on 1 July    2015
  if (tepoch >= 63158402) leap = leap + 1; // Leapsecond on 1 January 2016

  double tlocal = tepoch - leap;

  // Convert epoch time to Julian day 
  int32_t jd = (long) (tlocal/86400) + 2457024; // Julian day of 1/1/2015
  jdate( jd, iyr, idy);

  // Extract seconds of day
  *sec = tlocal - (int) (tlocal / 86400) * 86400;

  return 0;
}

/*
int get_packet_from_frame( ifstream *framefile, uint8_t **packet,
                           uint32_t &packet_length, int first)  {

  int LFRM = 892;
  int nbytep, iptr, pptr, nptr, nbytes, apid, len;
  uint8_t phead[6], ifct;
  
  static int fptr;
  static uint8_t frame[892];
    
  // Check if need to read first frame
  if (first) {
    //cout << "In first 1: " << framefile->tellg() << endl;
    framefile->read( (char *) frame, sizeof(frame));
    fptr = (frame[4] % 8)*256 + frame[5] + 6;
    //cout << "fptr1: " << fptr << endl;
    //cout << "In first 2: " << framefile->tellg() << endl;
  } // first


  // Check for non-continuation frame
  while (fptr == 2052 || fptr == 2053) {
    framefile->read( (char *) frame, sizeof(frame));
    if (framefile->eof()) {
      cout << "EOF before 2052/2053 frame found" << endl;
      exit(110);
    }
    fptr = (frame[4] % 8)*256 + frame[5] + 6;
  }

  //cout << "fptr2: " << fptr << endl;

  // Check for packet header split across frames
  // Populate 6-byte packet header
  if (fptr > (LFRM-6)) {
    nbytep = LFRM - fptr;
    memcpy( phead,  &frame[fptr], nbytep);
    framefile->read( (char *) frame, sizeof(frame));
    fptr = (frame[4] % 8)*256 + frame[5] + 6;
    iptr = 12 - nbytep;
    memcpy( &phead[nbytep], &frame[6], 6-nbytep);
  } else {
    memcpy( phead,  &frame[fptr], 6);
    iptr = fptr+ 6;
    //cout << "iptr: " << iptr << endl;
    nbytep = 0;
  }

  // Check for end-of-file (exposure) packet header
  apid = (phead[0] % 8)*256 + phead[1];
  //cout << "apdi: " << apid << endl;
  if (apid == 2047) {
    cout << "End of exposure" << endl;
    return 1;
  }

  // Get length from packet header
  len = phead[4]*256 + phead[5] + 7;
  packet_length = len;
  //cout << "packet length: " << len << endl;

  // Allocate packet and copy header
  *packet = new uint8_t[len];
  memcpy( &(*packet)[0], phead, 6);
  pptr = 6;

  // Check if entire packet is within current frame
  nptr = iptr + len - 6;
  //cout << "nptr: " << nptr << endl;
  if (nptr < LFRM) {
    memcpy( &(*packet)[6], &frame[iptr], nptr-iptr);
    fptr = nptr;
    //    cout << "Entire packet within frame" << endl;
    //cout << "fptr3: " << fptr << endl << endl;
    return 0;
  }

  // Otherwise get frames to complete packet
  nbytes = LFRM - iptr;
  if (nbytes > 0) memcpy( &(*packet)[pptr], &frame[iptr], nbytes);

  while (nptr >= LFRM) {
    pptr = pptr + nbytes;
    nptr = nptr - LFRM + 6;
    ifct = frame[2];
    framefile->read( (char *) frame, sizeof(frame));
    fptr = (frame[4] % 8)*256 + frame[5] + 6;
    nbytes = len - pptr;
    if (nbytes > (LFRM-6)) nbytes = LFRM-6;

    // Check for dropped frames

    // Must define diff as byte, c++ defaults to int  11/01/18
    uint8_t byte_diff = frame[2]-ifct;
    
    if (byte_diff == 1) {
      if (nbytes > 0) memcpy( &(*packet)[pptr], &frame[6], nbytes);
    } else {
      // If continuation frame
      while (fptr == 2053 && !framefile->eof()) {
        //cout << "continuation frame: " << framefile->tellg() << endl;
        framefile->read( (char *) frame, sizeof(frame));
        fptr = (frame[4] % 8)*256 + frame[5] + 6;
      }
      if (framefile->eof())
        return 1;
      //      cout << framefile->eof() << endl;
      nptr = fptr;
    }
  }

  if (fptr > LFRM) {
    cout << "Bad " << fptr << endl;
    exit(110);
  }
  
  return 0;
}
*/

int readFrame ( ifstream *framefile, uint8_t frame[892], int& prevFrameCnt,
                int ierror[2]) {

  ierror[0] = 0;

  if (!framefile->eof()) {
    framefile->read( (char *) frame, 892);

    // get frameCnt
    int frameCnt = (int) frame[2];

    // Determine the change in frame count
    // If not next mod 256 then frame drop
    int deltaCnt = frameCnt - prevFrameCnt;
    //long long int pos = framefile->tellg();   // V0.81

    if (deltaCnt != 1 && deltaCnt != -255) {
       
      // cout << "Frame Drop: " << prevFrameCnt << " " << frameCnt << " " <<
      //   pos << endl; // V0.82

      // Find first good packet
      int framePtr = 2046;
      int first2bytes = 20*256 + 160;
      while ((framePtr == 2046 || framePtr == 2047 || first2bytes != 5280) && \
             !framefile->eof()) {
        framefile->read( (char *) frame, 892);
        framePtr = (frame[4] % 8)*256 + frame[5];
        first2bytes = frame[0]*256 + frame[1];
      }
      framePtr += 6;
            
      frameCnt = (int) frame[2];
      ierror[0] = 1;
      ierror[1] = framePtr; //i % 886 + 6;
    }

    prevFrameCnt = frameCnt;
  }
  
  return 0;
}


int getSHpacket( ifstream *framefile, uint8_t frame[892], int& framePtr,
                 uint8_t **packet, int& packetLength, int& prevFrameCnt,
                 int& frameDrop, bool isVerbose)  {

  // Note: At normal end of this routine the framePtr is at the byte
  // just after the data section.

  // Note: secondHdr[4]*256+secondHdr[5] includes last 2 bytes of PUS header

  int ierror[2];
  uint8_t secondHdr[6];
  uint8_t pusHdr[3];

  int bytesWritten;
  int prevFramePtr;
  int toWrite2nd;
  int toWritePUS;
  int dataLeftToWrite;
  int packetPtr;
  
  int apid;
  int dataLen;

  // If framePtr just past current frame then read new frame
  if (framePtr == 892) {
    // Read new frame
    readFrame( framefile, frame, prevFrameCnt, ierror);

    // Frame drop
    if (ierror[0] == 1) {
      frameDrop = 1;
      framePtr = ierror[1];
      return 0;
    }
    framePtr = PHDRLEN;
  }

  if (framePtr > (892-6)) {

    //////////////////////////////////////////////////
    //////////////// Incomplete 2nd hdr //////////////
    //////////////////////////////////////////////////

    //         S    E    C    O    N    D    P    U    S
    //        887                 891
    //        888            891
    //        889       891
    //        890  891
    //        891

    if (framePtr > 892) {
      long long int pos = framefile->tellg();   // V0.81
      if (isVerbose) cout << "Bad frame pointer: " << framePtr << " " << pos << endl;

      // Force read of next frame
      framePtr = 892;
      frameDrop = -1;
      return 0;
    }

    // Write partial 2nd header
    memcpy( secondHdr, &frame[framePtr], 892-framePtr);

    // Read new frame
    readFrame( framefile, frame, prevFrameCnt, ierror);
    
    // Frame drop
    if (ierror[0] == 1) {
      frameDrop = 1;
      framePtr = ierror[1];
      return 0;
    }

    // Set prevFramePtr at start of 2nd header in previous frame
    // Set framePtr at byte after primary header in new frame
    prevFramePtr = framePtr;
    framePtr = PHDRLEN;

    // 2nd header bytes written & 2nd header bytes left to write
    bytesWritten = 892 - prevFramePtr;
    toWrite2nd = SHDRLEN - bytesWritten;

    // Write remaining part of 2nd header
    memcpy( &secondHdr[bytesWritten], &frame[framePtr], toWrite2nd);

    // Set framePtr at byte after remainder of 2nd header
    framePtr += toWrite2nd;

    // Write PUS header
    memcpy( pusHdr, &frame[framePtr], PUSLEN);
    
    // V0.82
    if (secondHdr[0]==7) {
       // an End of Exposure packet
       cout<<"EOI packet detected type1"<<endl;
      return 1;       
    }

    // Get length of data
    dataLen = secondHdr[4]*256 + secondHdr[5] - 2;
    if (dataLen < 0) {
      long long int pos = framefile->tellg();   // V0.81
      if (isVerbose) cout << "Bad dataLen type1: " << dataLen << " " << pos << endl;

      // Force read of next frame
      framePtr = 892;
      frameDrop = -1;
      return 0;
    }

    // Write 2nd/PUS headers to packet
    *packet = new uint8_t[SHDRLEN+PUSLEN+dataLen];

    memcpy( &(*packet)[0], secondHdr, SHDRLEN);
    memcpy( &(*packet)[SHDRLEN], pusHdr, PUSLEN);
    packetPtr = SHDRLEN+PUSLEN;
    
    apid = ((*packet)[0] % 8) * 256 + (*packet)[1];
    if (((*packet)[7] < 253 || (*packet)[7] > 255) && apid != 2047) {
      long long int pos = framefile->tellg();   // V0.81
      if (isVerbose) cout << "Bad PUS type1: " << (int) (*packet)[7] << " " << pos << endl;

      // Force read of next frame
      framePtr = 892;
      frameDrop = -1;
      return 0;
    }
    
    // Set prevFramePtr at byte after PUS header (start of data)
    prevFramePtr = framePtr + PUSLEN;
        
    // Set framePtr to byte after end of data
    framePtr = prevFramePtr + dataLen;

    // If data is within current frame then
    // fill remainder of packet and return
    if (framePtr <= 892) {
      memcpy( &(*packet)[packetPtr], &frame[prevFramePtr], dataLen);
      packetLength = SHDRLEN + PUSLEN + dataLen;
      return 0;
    }
     
  } else if (framePtr == (892-6)) {

    //////////////////////////////////////////////////
    //////// Complete 2nd hdr, Missing PUS hdr ///////
    //////////////////////////////////////////////////
        
    //         S    E    C    O    N    D    P    U    S         
    //        886                      891

    // Write 2nd header
    memcpy( secondHdr, &frame[framePtr], SHDRLEN);

    // Read new frame
    readFrame( framefile, frame, prevFrameCnt, ierror);
    
    // Frame drop
    if (ierror[0] == 1) {
      frameDrop = 1;
      framePtr = ierror[1];
      return 0;
    }

    // Write PUS header
    memcpy( pusHdr, &frame[PHDRLEN], PUSLEN);

    // Get length of data
    dataLen = secondHdr[4]*256 + secondHdr[5] - 2;
    if (dataLen < 0) {
      long long int pos = framefile->tellg();   // V0.81
      if (isVerbose) cout << "Bad dataLen type2: " << dataLen << " " << pos << endl;

      // Force read of next frame
      framePtr = 892;
      frameDrop = -1;
      return 0;
    }

    // Write 2nd/PUS headers to packet
    *packet = new uint8_t[SHDRLEN+PUSLEN+dataLen];

    memcpy( &(*packet)[0], secondHdr, SHDRLEN);
    memcpy( &(*packet)[SHDRLEN], pusHdr, PUSLEN);
    packetPtr = SHDRLEN + PUSLEN;
    
    apid = ((*packet)[0] % 8) * 256 + (*packet)[1];
    if (((*packet)[7] < 253 || (*packet)[7] > 255) && apid != 2047) {
      long long int pos = framefile->tellg();   // V0.81
      if (isVerbose) cout << "Bad PUS type2: " << (int) (*packet)[7] << " " << pos << endl;

      // Force read of next frame
      framePtr = 892;
      frameDrop = -1;
      return 0;
    }
    
    // Set prevFramePtr at byte after PUS header (start of data)
    prevFramePtr = PHDRLEN + PUSLEN;

    // Set framePtr to byte after end of data
    framePtr = prevFramePtr + dataLen;

    // If data is within current frame then
    // fill remainder of packet and return
    if (framePtr <= 892) {
      memcpy( &(*packet)[packetPtr], &frame[prevFramePtr], dataLen);
      packetLength = SHDRLEN + PUSLEN + dataLen;
      return 0;
    }
    
  } else if (framePtr > (892-9)) {

    //////////////////////////////////////////////////
    ////// Complete 2nd hdr, Incomplete PUS hdr //////
    //////////////////////////////////////////////////
        
    //         S    E    C    O    N    D    P    U    S
    //        884                                891
    //        885                           891

    // Write 2nd header
    memcpy( secondHdr, &frame[framePtr], SHDRLEN);

    // Write partial PUS header
    memcpy( pusHdr, &frame[framePtr+SHDRLEN], 892-(framePtr+SHDRLEN));

    // Read new frame
    readFrame( framefile, frame, prevFrameCnt, ierror);
    
    // Frame drop
    if (ierror[0] == 1) {
      frameDrop = 1;
      framePtr = ierror[1];
      return 0;
    }

    // Set prevFramePtr at start of PUS header in previous frame
    // Set framePtr at byte after primary header in new frame
    prevFramePtr = framePtr + SHDRLEN;
    framePtr = PHDRLEN;

    // PUS header bytes written & PUS header bytes left to write
    bytesWritten = 892 - prevFramePtr;
    toWritePUS = PUSLEN - bytesWritten;

    // Write remaining part of PUS header
    memcpy( &pusHdr[bytesWritten], &frame[framePtr], toWritePUS);

    // Get length of data
    dataLen = secondHdr[4]*256 + secondHdr[5] - 2;
    if (dataLen < 0) {
      long long int pos = framefile->tellg();   // V0.81
      if (isVerbose) cout << "Bad dataLen type3: " << dataLen << " " << pos << endl;

      // Force read of next frame
      framePtr = 892;
      frameDrop = -1;
      return 0;
    }

    // Write 2nd/PUS headers to packet
    *packet = new uint8_t[SHDRLEN+PUSLEN+dataLen];

    memcpy( &(*packet)[0], secondHdr, SHDRLEN);
    memcpy( &(*packet)[SHDRLEN], pusHdr, PUSLEN);
    packetPtr = SHDRLEN+PUSLEN;
    
    apid = ((*packet)[0] % 8) * 256 + (*packet)[1];
    if (((*packet)[7] < 253 || (*packet)[7] > 255) && apid != 2047) {
      long long int pos = framefile->tellg();   // V0.81
      if (isVerbose) cout << "Bad PUS type3: " << (int) (*packet)[7] << " " << pos << endl;

      // Force read of next frame
      framePtr = 892;
      frameDrop = -1;
      return 0;
    }

    // Set prevFramePtr at byte after PUS header (start of data)
    prevFramePtr = framePtr + toWritePUS;
        
    // Set framePtr to byte after end of data
    framePtr = prevFramePtr + dataLen;

    // If data is within current frame then
    // fill remainder of packet and return
    if (framePtr <= 892) {
      memcpy( &(*packet)[packetPtr], &frame[prevFramePtr], dataLen);
      packetLength = SHDRLEN + PUSLEN + dataLen;
      return 0;
    }

  } else {

    //////////////////////////////////////////////////
    //////////// Complete 2nd hdr/PUS hdr ////////////
    //////////////////////////////////////////////////

    //         S    E    C    O    N    D    P    U    S
    //        883                                     891
    //        882  883                                     891
    //        881  882  883                                     891         

    // Write 2nd/PUS headers
    memcpy( secondHdr, &frame[framePtr], SHDRLEN);
    memcpy( pusHdr, &frame[framePtr+SHDRLEN], PUSLEN);
    
    /*
    // V0.84
    if (secondHdr[0]==7) {
       // an End of Exposure packet
       cout<<"EOI packet detected type4"<<endl;
      return 1;       
    }
	*/
	
    // Get length of data
    dataLen = secondHdr[4]*256 + secondHdr[5] - 2;
    if (dataLen < 0) {
      long long int pos = framefile->tellg();   // V0.81
      if (isVerbose) cout << "Bad dataLen type4: " << dataLen << " " << pos << endl;
	  
      // Force read of next frame
      framePtr = 892;
      frameDrop = -1;
      return 0;
    }

    // Write 2nd/PUS headers to packet
    *packet = new uint8_t[SHDRLEN+PUSLEN+dataLen];

    memcpy( &(*packet)[0], secondHdr, SHDRLEN);
    memcpy( &(*packet)[SHDRLEN], pusHdr, PUSLEN);
    packetPtr = SHDRLEN+PUSLEN;
    
    apid = ((*packet)[0] % 8) * 256 + (*packet)[1];
    if (((*packet)[7] < 253 || (*packet)[7] > 255) && apid != 2047) {
      long long int pos = framefile->tellg();   // V0.81
      if (isVerbose) cout << "Bad PUS type4: " << (int) (*packet)[7] << " " << pos << endl;

      // Force read of next frame
      framePtr = 892;
      frameDrop = -1;
      return 0;
    }

    // Set prevFramePtr at byte after PUS header (start of data)
    prevFramePtr = framePtr + SHDRLEN + PUSLEN;

    // Set framePtr to byte after end of data
    framePtr = prevFramePtr + dataLen;
    
    // If data is within current frame then
    // fill remainder of packet and return
    if (framePtr <= 892) {
      memcpy( &(*packet)[packetPtr], &frame[prevFramePtr], dataLen);
      packetLength = SHDRLEN + PUSLEN + dataLen;
      return 0;
    }
  }

  //////////////////////////////////////////////////
  /////////////// Multi-frame packets //////////////
  //////////////////////////////////////////////////    

  // Write data remaining in current frame
  // Note: If complete 2nd/PUS headers at end of frame then no data to write
  if (prevFramePtr < 892) {
    memcpy( &(*packet)[packetPtr], &frame[prevFramePtr], 892-prevFramePtr);
    packetPtr += 892-prevFramePtr;
  }

  // Determine number of bytes left to write
  dataLeftToWrite = dataLen - (892 - prevFramePtr);

  while (dataLeftToWrite > 0) {

    // If data to write then read new frame
    readFrame( framefile, frame, prevFrameCnt, ierror);
    
    // Frame drop
    if (ierror[0] == 1) {
      frameDrop = 1;
      framePtr = ierror[1];
      return 0;
    }
    
    // Set previous frame to start of data (after primary header)
    prevFramePtr = PHDRLEN;

    // Set frame pointer to minimum of frame length or pointer to byte
    // following end of data in current frame
    if (892 < prevFramePtr+dataLeftToWrite)
      framePtr = 892;
    else
      framePtr =  prevFramePtr + dataLeftToWrite;

    int dataToWrite = framePtr - prevFramePtr;
    
    // Write data to packet
    memcpy( &(*packet)[packetPtr], &frame[prevFramePtr], dataToWrite);
    packetPtr += dataToWrite;
    
    // Compute remaining data to write
    dataLeftToWrite -= framePtr - prevFramePtr;
  }

  packetLength = SHDRLEN + PUSLEN + dataLen;

  return 0;
}


