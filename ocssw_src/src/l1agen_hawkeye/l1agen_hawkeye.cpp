#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

#include "hawkeyeUtil.h"
#include "HawkeyeDecode.h"
#include "nc4utils.h"
#include "timeutils.h"

#include "l1agen_hawkeye.h"
#include "netcdf.h"

#define VERSION "1.0.5"

//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     FutureTech     06/16/17 0.10  Original development
//  Joel Gales     FutureTech     05/01/18 0.20  Generate NETCDF4 output
//                                               using CDL file
//  Joel Gales     FutureTech     08/02/18 0.30  Process telemetry files
//                                               directly to extract HEB &
//                                               ADCS data
//
//  Joel Gales     FutureTech     09/20/18 0.31  Pass number of S/C records
//                                               to createl1
//  Joel Gales     FutureTech     09/27/18 0.40  Fix ADCS/Hawkeye time issues,
//                                               Only write ADCS records within
//                                               Hawkeye time period.  Store
//                                               image parameter data as
//                                               attributes.
//  Joel Gales     FutureTech     10/23/18 0.50  Change software_version,
//                                               fpga_version,
//                                               telemetry_counter, overcurrent
//                                               fields from ushort to short,
//                                               delta_time from uint to int
//  Joel Gales     SAIC           11/01/18 0.60  Fix byte diff bug in
//                                               get_packet_from_frame()
//                                               uint8_t byte_diff=
//                                               frame[2]-ifct;
//  Joel Gales     SAIC           11/06/18 0.70  Find and write all bracketing
//                                               ADCS records within telemetry
//                                               stream.
//  Joel Gales     SAIC           03/18/19 0.80  Skip navigation output if no
//                                               good navigation records.
//  Joel Gales     SAIC           03/25/19 0.85  Fix problem with finder image
//                                               output buffer.
//  Joel Gales     SAIC           03/28/19 0.90  Fix code to support 25 finder
//                                               images.
//  Joel Gales     SAIC           04/03/19 0.91  Add finder_time and
//                                               finder_delta_time fields
//  Joel Gales     SAIC           05/28/19 0.92  Modifiy to read new L0 granule
//                                               format.
//  Joel Gales     SAIC           06/19/19 0.93  Initialize instrument
//                                               telemetry to fill values
//  Joel Gales     SAIC           07/31/19 0.94  Modify createl1 function
//                                               and image/finder reader
//                                               to handle non-standard
//                                               size output
//  Joel Gales     SAIC           08/01/19 0.95  Exit 110 if "empty" L1A will
//                                               be generated
//  Joel Gales     SAIC           10/25/19 0.96  Trap negative dataLen in
//                                               getSHpacket() for all cases
//  Joel Gales     SAIC           10/25/19 0.97  Trap zero finder dimensions
//                                               Close file and return
//  Joel Gales     SAIC           03/20/20 0.98  Trap buffer overrun in
//                                               FindHeader() when writing
//                                               finder images
//  Liang Hong     SAIC           10/02/20 0.981 Check image height limit
//  Liang Hong     SAIC           11/09/20 0.982 Exit with band width/height error
//  Liang Hong     SAIC           01/21/21 0.983 Sensor time, bus_telemetry source
//												 data array boundary fix; 
//												 HawkeyeDecodeSpectralImage avoid 
//												 infinite loop; quit if image height
//                                               over MAX; 
//  Liang Hong     SAIC           03/22/21 0.993 Second of day refers to image start
//												 date
//  Liang Hong     SAIC           03/31/21 0.994 Added date_created attribute
//  Liang Hong     SAIC           07/01/21 0.995 Continue processe with excess S/P packets 
//  Liang Hong     SAIC           07/12/21 0.996 Boundary check for sensor time array
//  Liang Hong     SAIC           08/31/21 0.997 Update band height limit when exceeding
//  Liang Hong     SAIC           09/28/21 0.998 CCD temp. unsigned to signed fix
//  Liang Hong     SAIC           10/11/21 0.999 CCD temp. -2 degC set to invalid
//  Liang Hong	   SAIC           10/27/21 1.0.0 re-formatted a few IF statements
//  Liang Hong     SAIC           01/13/22 1.0.1 abnormal finderscope image data handling
//  Liang Hong	   SAIC           02/17/22 1.0.2 flag scan time when delta time is 0
//  Liang Hong	   SAIC           03/25/22 1.0.3 keep scan time for the 1st line vs 1.0.2
//  Liang Hong     SAIC           04/19/22 1.0.4 fixed false alarm on no finder image
//  Liang Hong     SAIC           04/21/23 1.0.5 fixed missing endtime when deltatime = 0

using namespace std;

// Note: To get an HEB dump just supply to L0 file on the command line.

int main (int argc, char* argv[])
{

  cout << "l1agen_hawkeye " << VERSION << " (" 
       <<  __DATE__ << " " << __TIME__ << ")" << endl;

  if ( argc == 1) {
    cout << endl << 
      "l1agen_hawkeye " <<
      "input_l0_filename input_l0_textfilename output_l1a_filename" << endl;

    return 0;
  }

  // Open L0 telemetry file for reading
  ifstream tlmfile;
  tlmfile.open( argv[1], ios::binary | ifstream::in);

  
  // Open L0 text file
  istringstream istr;
  ifstream l0Text;
  char txtbuf[512];
  double startHWKTime;
  double stopHWKTime;
  bool bImageCrossDay = false;
  int32_t startyr, startdy, stopyr, stopdy;
  double startsec, stopsec;
  uint32_t nSC;
  if (argc > 2) {
    l0Text.open( argv[2], ifstream::in);
    if ( l0Text.fail() != 0) {
      cout << endl << "L0 Text File: " << argv[2] << " not found." << endl;
      exit(1);
    }

    // Read StartHWKTime/StopHWKTime
    while( !l0Text.eof()) {
      size_t found;
      l0Text.getline( txtbuf, 512);
      string sValue = txtbuf;
    
      found = sValue.find("StartHWKTime");
      if ( found != string::npos) {
        found = sValue.find("=");
        istr.str(sValue.substr(found+1));
        istr >> startHWKTime;
        istr.clear();
      }

      found = sValue.find( "StopHWKTime");
      if ( found != string::npos) {
        found = sValue.find("=");
        istr.str(sValue.substr(found+1));
        istr >> stopHWKTime;
        istr.clear();
      }
    }

    nSC = (uint32_t) (stopHWKTime - startHWKTime) + 2;
    
    tepoch2yds(startHWKTime, &startyr, &startdy, &startsec);
    tepoch2yds(stopHWKTime, &stopyr, &stopdy, &stopsec);
    if ((stopyr*1000+stopdy)>(startyr*1000+startdy)) bImageCrossDay = true;
  }
  

  // Allocate ADCS data buffer pointer arrays  
  sensor_t *sensor[MAXNPACKETS];
  psensor_t *psensor[MAXNPACKETS];
  attitude_t *attitude[MAXNPACKETS];
  propagator_t *propagator[MAXNPACKETS];

  uint32_t n[4]={0,0,0,0};
  uint32_t packet_length;

  // Image data buffer
  vector<uint8_t> heb_buf;

  uint8_t **packet = new uint8_t*;
  uint32_t nImages=0, nFinder=0;
    
  //  int first = 1;
  while  (1) {
    int tlmfile_pos = tlmfile.tellg();
    
    //    int status = get_packet_from_frame( &tlmfile, packet, packet_length,
    //                                    first);

    //    cout << tlmfile_pos << endl;
    
    int status;
    uint8_t header[6];
    tlmfile.read( (char *) header, 6);

    if( tlmfile.eof()) break;
    
    int dataLen = header[4]*256 + header[5] - 2;
    packet_length = 6 + 3 + dataLen;
    tlmfile.seekg( -6, ios::cur);

    (*packet) = new uint8_t[packet_length];
    tlmfile.read( (char *) (*packet), packet_length);
    
    //    if (status == 1) {
    //first = 1;
    // break;
    //}
    
    // Each packet contains a CCSDS 6-byte header followed by a SeaHawk PUS
    // subheader. The second byte of the subheader indicates the data type:

    //    253     Hawkeye data
    //    254     GPS data
    //    255     ADCS data

    // The Hawkeye packets are slightly modified versions of the .heb file
    // blocks.  The block start sequence is replaced by the CCSDS header and
    // the first two bytes of the PUS subheader.  The block type is contained
    // in the third byte of the subheader.  Also, the block length field is
    // now the packet length, which is 7 bytes longer than the CCSDS header
    // packet length.
    
    if ((*packet)[7] == 255) {
      short psub = (*packet)[8];
      
      if (psub == 1) {
        sensor[n[0]] = new sensor_t;
      } else if (psub == 2) {
        psensor[n[1]] = new psensor_t;        
      } else if (psub == 3) {
        attitude[n[2]] = new attitude_t;
      } else if (psub == 4) {
        propagator[n[3]] = new propagator_t;
      }
	
	  // Ver 0.995
	  if (n[psub-1] >= MAXNPACKETS-1) {
		cout << "Insufficient pointer allocation for psub: " << psub << endl;
        cout << "Increase value of MAXNPACKETS" << endl;
	  	continue;
	  }
	  
      status = unpack_seahawk_adcs( *packet, startHWKTime, stopHWKTime,
                                    sensor[n[0]], psensor[n[1]],
                                    attitude[n[2]], propagator[n[3]]);

      if (status == 0) n[psub-1]++;
      
      //if (n[psub-1] == MAXNPACKETS) {
      //  cout << "Insufficient pointer allocation for psub: " << psub << endl;
      //  cout << "Increase value of MAXNPACKETS" << endl;
      //  exit(1);
      //}

    } else if ((*packet)[7] == 253) {
      // HEB data
      uint16_t row;

      // Compressed Image/Finder Data
      if ( (*packet)[8] == 2) {
        uint8_t infoByte = (*packet)[15];
        uint8_t sb = infoByte & 0x1f;

        // sd_f: 0 - finder, 1 - spectral (image)
        uint8_t sd_f = (infoByte & 0x40) >> 6;
        uint8_t dark = (infoByte & 0x80) >> 7;
        memcpy(&row, &(*packet)[20], 2);
        row = SWAP_2(row);

        if (sd_f == 0 && sb > nFinder) nFinder = sb;
        if (sd_f == 1 && sb > nImages) nImages = sb;
        
        if ( argc == 2) {
          cout << "SB or Finder #: " << (int) sb <<
            "  SpecData or Finder: " << (int) sd_f <<
            " Dark: " << (int) dark << "  Row #: " << row <<
            "  packet length: " << packet_length <<
            "  tlmfile_pos: " << tlmfile_pos << endl;
        }
      }

      // Image Parameters Data Block
      if ( (*packet)[8] == 3) {
        //cout << "Image Parameter" << endl;
        uint16_t hgt;
        memcpy(&hgt, &(*packet)[48], 2);
        //        uint8_t nFinder = (*packet)[55];
        //cout << "SB hgt: " << SWAP_2(hgt) <<
        //  "  # find img: " << (int) nFinder << endl;
      }
      //if ( (*packet)[8] == 4) cout << "Telemetry" << endl;
      //      if ( (*packet)[8] == 6) cout << "Mission Log" << endl;

      
      // Block start sequence
      uint8_t bss[6] = {0x00,0x00,0xff,0xff,0xff,0xff};

      // Compute new checksum
      uint16_t checksum = 4*255;
      for (size_t i=8; i<packet_length; i++) {
        checksum += (*packet)[i];
      }
      uint16_t swapcheck = SWAP_2(checksum);

      // Write HEB data
      heb_buf.insert(heb_buf.end(), bss, bss+6);
      heb_buf.insert(heb_buf.end(), &(*packet)[8], &(*packet)[packet_length]);
      heb_buf.push_back(swapcheck & 0x00ff);
      heb_buf.push_back((swapcheck & 0xff00) >> 8);
    } else {
      //      cout << (int) (*packet)[7] << " " << (int) (*packet)[8] << endl;

    }
    delete[] (*packet);
  }
  uint32_t hebBufLength = heb_buf.size();
  
  cout << "Number of SENSOR packets:     " << n[SENSOR] << endl;
  cout << "Number of PSENSOR packets:    " << n[PSENSOR] << endl;
  cout << "Number of ATTITUDE packets:   " << n[ATTITUDE] << endl;
  cout << "Number of PROPAGATOR packets: " << n[PROPAGATOR] << endl << endl;

  cout << "Number of Spectral Images:    " << nImages << endl;
  cout << "Number of Finder Images:      " << nFinder << endl << endl;
  
  if ( argc == 2) {
    ofstream hebFile ("heb_buf.bin", ios::out | ios::binary);
    hebFile.write((char *) heb_buf.data(), hebBufLength);
    hebFile.close();
    
    return 0;
  }

  HawkeyeStreamInfo streamInfo;
  string l1a_name;
  l1a_name.assign( argv[3]);

  // Image height is determined by counting the number of "good" rows
  // Image width is determined in ScanRow() (for each row which are identical)
  HawkeyeScanStream( (uint8_t*) heb_buf.data(), hebBufLength, &streamInfo);

  // Bail if empty
  if ( streamInfo.spectralInfo[0].width == 0) {
    cout << "Empty output file: " << l1a_name.c_str() << endl;
    exit(110);
  }
  
   // check image height not to exceed max value in case of corrupted bytes in raw data reading
   // LH, 1/21/2021
   // updated band height checking, LH, 8/31/2021
   uint16_t band_max_height = 0;
   for (size_t i=0; i<8; i++) {
   		if (streamInfo.spectralInfo[i].height>band_max_height && streamInfo.spectralInfo[i].height<=MAXIMGHEIGHT) {
   			band_max_height = streamInfo.spectralInfo[i].height;
   		}
   }
   if (band_max_height==0) {
		cout << "Image height in all bands exceeds MAXIMGHEIGHT " << endl;
		//cout << "Empty output file: " << l1a_name.c_str() << endl;
		//exit(110);
		band_max_height = MAXIMGHEIGHT;
   }
   for (size_t i=0; i<8; i++) {
		if (streamInfo.spectralInfo[i].height>MAXIMGHEIGHT) {
			cout << "Image height " <<streamInfo.spectralInfo[i].height<<" in band " <<(i+1)<< " exceeds MAXIMGHEIGHT " << endl;
			streamInfo.spectralInfo[i].height = band_max_height;
		}
   }
  
  // Full size: 18 byte offset / 2 byte time code and 16 bytes dark pixels
  // Half size: 12 byte offset / 2 byte time code and 10 bytes dark pixels
  uint32_t rowOffset;
  if ( streamInfo.spectralInfo[0].width == 1816) rowOffset = 18;
  if ( streamInfo.spectralInfo[0].width ==  908) rowOffset = 10;

  static l1aFile outfile;
  outfile.createl1( (char *) l1a_name.c_str(), nSC,
                    streamInfo.spectralInfo[0].width+2-rowOffset,
                    streamInfo.spectralInfo[0].height,
                    streamInfo.finderscopeInfo[0].width,
                    streamInfo.finderscopeInfo[0].height);

  // scan_line_attributes
  // parameters_telemetry_data
  // navigation_data
  // earth_view_data

  char buf[32];
  strcpy( buf, unix2isodate(now(), 'G'));

  int status;
  int varid;
  string varname;

  
  ////////////////////// Write Telemetry //////////////////////
  int dimid;
  size_t nTlmBlocks, sizeTelBlk, nCCDtemps, nFPGAvolts, nCurrents, nCCDvolts;

  nc_inq_dimid( outfile.ncid, "number_of_tlm_blocks", &dimid);
  nc_inq_dimlen( outfile.ncid, dimid, &nTlmBlocks);

  nc_inq_dimid( outfile.ncid, "telemetry_block", &dimid);
  nc_inq_dimlen( outfile.ncid, dimid, &sizeTelBlk);

  nc_inq_dimid( outfile.ncid, "ccd_temps", &dimid);
  nc_inq_dimlen( outfile.ncid, dimid, &nCCDtemps);

  nc_inq_dimid( outfile.ncid, "ccd_volts", &dimid);
  nc_inq_dimlen( outfile.ncid, dimid, &nCCDvolts);

  nc_inq_dimid( outfile.ncid, "fpga_volts", &dimid);
  nc_inq_dimlen( outfile.ncid, dimid, &nFPGAvolts);

  nc_inq_dimid( outfile.ncid, "currents", &dimid);
  nc_inq_dimlen( outfile.ncid, dimid, &nCurrents);

  
  uint16_t *telemetry = new uint16_t[nTlmBlocks*sizeTelBlk];
  double *tlm_time = new double[nTlmBlocks];
  uint32_t *tlm_time_stamp = new uint32_t[nTlmBlocks];
  int16_t *software_version = new int16_t[nTlmBlocks];
  int16_t *fpga_version = new int16_t[nTlmBlocks];
  int16_t *telemetry_counter = new int16_t[nTlmBlocks];
  
  float *ccd_temperatures = new float[nTlmBlocks*nCCDtemps];
  for (size_t i=0; i<nTlmBlocks*nCCDtemps; i++) ccd_temperatures[i] = -999.0;
  
  float *fpga_temperature = new float[nTlmBlocks];
  for (size_t i=0; i<nTlmBlocks; i++) fpga_temperature[i] = -999.0;
  
  float *fpga_voltages = new float [nTlmBlocks*nFPGAvolts];
  for (size_t i=0; i<nTlmBlocks*nFPGAvolts; i++) fpga_voltages[i] = -999.0;
 
  int16_t *overcurrent = new int16_t[nTlmBlocks];

  float *current_monitors = new float[nTlmBlocks*nCurrents];
  for (size_t i=0; i<nTlmBlocks*nCurrents; i++) current_monitors[i] = -999.0;
  
  float *ccd_voltages = new float[nTlmBlocks*nCCDvolts];
  for (size_t i=0; i<nTlmBlocks*nCCDvolts; i++) ccd_voltages[i] = -999.0;
  
  float *solenoid_voltage = new float[nTlmBlocks];
  for (size_t i=0; i<nTlmBlocks; i++) solenoid_voltage[i] = -999.0;
  
  int32_t iyr, idy;
  double sec;
  
  for (size_t i=0; i<nTlmBlocks; i++) {
    //    cout << "i: " << i << endl;
    tlm_time_stamp[i] = streamInfo.telemetryInfo[i].timeStamp;

    tlm_time[i] = (streamInfo.imageInfo.epochT0 + tlm_time_stamp[i]) * 0.001;

    // JD 2015/01/01 = 2457023.5
    // JD 1980/01/06 = 2444244.5
    // delta sec = (2457023.5-2444244.5) * 24 * 3600 = 1104105600
    // Add 16 leapsecs between 1980/01/06 and 2015/01/01
    // Convert to milliseconds: 1104105616000
    // Subtract 2^40 = 1099511627776
    // delta = 4593988224
    // Convert to days: 4593988224 / (24*3600*1000) = 53.17116
    // 0.17116 * 24 * 3600 = 14788.224; 

    tlm_time[i] -= (53*24*3600 + 14788.224);
    tepoch2yds( tlm_time[i], &iyr, &idy, &sec);

    tlm_time[i] = sec;
    // if image spans two days
    if (bImageCrossDay && (sec<80000)) tlm_time[i] = sec+86400.0;

    memcpy( &telemetry[i*sizeTelBlk], &streamInfo.telemetryInfo[i].telemetry,
            sizeTelBlk*sizeof(short));
            
    uint32_t nChannels = streamInfo.telemetryInfo[i].noChannels;
    for (size_t j=0; j<nChannels; j++) {
      uint16_t value =
        streamInfo.telemetryInfo[i].telemetry.channel[j].channelValue;

      uint16_t interp =
        streamInfo.telemetryInfo[i].telemetry.channel[j].channelInterp;

      switch (j) {
          
      case TC_SOFTWARE_VERSION:
        if ( interp == 1)
          software_version[i] = value;
        else
          software_version[i] = -999;
        break;

      case TC_FPGA_VERSION:
        if ( interp == 1)
          fpga_version[i] = value;
        else
          fpga_version[i] = -999;
        break;

      case TC_INDEX:
        if ( interp == 1)
          telemetry_counter[i] = value;
        else
          telemetry_counter[i] = -999;
        break;

      case TC_CCD1_TEMP:
      case TC_CCD2_TEMP:
      case TC_CCD3_TEMP:
      case TC_CCD4_TEMP:
        if ( interp == 1)  {
        // ver 0.998: CCD temperatures are signed values but read as unsigned.  LH
          ccd_temperatures[i*nCCDtemps+(j-TC_CCD1_TEMP)] = (value<32768)? (value / 32.0):(value / 32.0-2048);
        }
		if (i*nCCDtemps+(j-TC_CCD1_TEMP)>0) { // ver 0.999: -2 degC readings following an invalid value are flagged as invalid
			if ((ccd_temperatures[i*nCCDtemps+(j-TC_CCD1_TEMP)-1]==-999.0) && (ccd_temperatures[i*nCCDtemps+(j-TC_CCD1_TEMP)]==-2.0)) {
				ccd_temperatures[i*nCCDtemps+(j-TC_CCD1_TEMP)] = -999.0;
			}
		}
        else
          ccd_temperatures[i*nCCDtemps+(j-TC_CCD1_TEMP)] = -999.0;
        break;

      case TC_FPGA_TEMP:
        if ( interp == 1)
          fpga_temperature[i] = (float) value;
        else
          fpga_temperature[i] = -999.0;
        break;

      case TC_FPGA_VAUX:
      case TC_FPGA_VINT:
      case TC_FPGA_VNVP:
        if ( interp == 1)
          fpga_voltages[i*nFPGAvolts+(j-TC_FPGA_VAUX)] = value / 1000.0;
        else
          fpga_voltages[i*nFPGAvolts+(j-TC_FPGA_VAUX)] = -999.0;
        break;

      case TC_CCD_VDD_OC:
        if ( interp == 1)
          overcurrent[i] = value;
        else
          overcurrent[i] = -999;
        break;
        
      case TC_AD7490_CH01:
      case TC_AD7490_CH02:
      case TC_AD7490_CH03:
      case TC_AD7490_CH04:
        if ( interp == 1)
          current_monitors[i*nCurrents+(j-TC_AD7490_CH01)] = value / 100.0;
        else
          current_monitors[i*nCurrents+(j-TC_AD7490_CH01)] = -999.0;
        break;

      case TC_AD7490_CH05:
      case TC_AD7490_CH06:
      case TC_AD7490_CH07:
      case TC_AD7490_CH08:
      case TC_AD7490_CH09:
      case TC_AD7490_CH10:
      case TC_AD7490_CH11:
      case TC_AD7490_CH12: 
      case TC_AD7490_CH13:
      case TC_AD7490_CH14:
      case TC_AD7490_CH15:
        if ( interp == 1)
          ccd_voltages[i*nCCDvolts+(j-TC_AD7490_CH05)] = value / 1000.0;
        else
          ccd_voltages[i*nCCDvolts+(j-TC_AD7490_CH05)] = -999.0;
        break;

      case TC_AD7490_CH16:
        if ( interp == 1)
          solenoid_voltage[i] = value / 1000.0;
        else
          solenoid_voltage[i] = -999.0;
      }
    } // channel loop

  } // tlm block loop


  varname.assign( "telemetry");
  status = nc_inq_varid( outfile.gid[1], varname.c_str(), &varid);
  status = nc_put_var_ushort( outfile.gid[1], varid, telemetry);
  check_err(status,__LINE__,__FILE__);
      
  varname.assign( "tlm_time_stamp");
  status = nc_inq_varid( outfile.gid[1], varname.c_str(), &varid);
  status = nc_put_var_uint( outfile.gid[1], varid, tlm_time_stamp);
  check_err(status,__LINE__,__FILE__);
  
  varname.assign( "tlm_time");
  status = nc_inq_varid( outfile.gid[1], varname.c_str(), &varid);
  status = nc_put_var_double( outfile.gid[1], varid, tlm_time);
  check_err(status,__LINE__,__FILE__);
  
  varname.assign( "software_version");
  status = nc_inq_varid( outfile.gid[1], varname.c_str(), &varid);
  status = nc_put_var_short( outfile.gid[1], varid, software_version);
  check_err(status,__LINE__,__FILE__);
  
  varname.assign( "FPGA_version");
  status = nc_inq_varid( outfile.gid[1], varname.c_str(), &varid);
  status = nc_put_var_short( outfile.gid[1], varid, fpga_version);
  check_err(status,__LINE__,__FILE__);
  
  varname.assign( "telemetry_counter");
  status = nc_inq_varid( outfile.gid[1], varname.c_str(), &varid);
  status = nc_put_var_short( outfile.gid[1], varid, telemetry_counter);
  check_err(status,__LINE__,__FILE__);
  
  varname.assign( "CCD_temperatures");
  status = nc_inq_varid( outfile.gid[1], varname.c_str(), &varid);
  status = nc_put_var_float( outfile.gid[1], varid, ccd_temperatures);
  check_err(status,__LINE__,__FILE__);
  
  varname.assign( "FPGA_temperature");
  status = nc_inq_varid( outfile.gid[1], varname.c_str(), &varid);
  status = nc_put_var_float( outfile.gid[1], varid, fpga_temperature);
  check_err(status,__LINE__,__FILE__);
  
  varname.assign( "FPGA_voltages");
  status = nc_inq_varid( outfile.gid[1], varname.c_str(), &varid);
  status = nc_put_var_float( outfile.gid[1], varid, fpga_voltages);
  check_err(status,__LINE__,__FILE__);
  
  varname.assign( "overcurrent");
  status = nc_inq_varid( outfile.gid[1], varname.c_str(), &varid);
  status = nc_put_var_short( outfile.gid[1], varid, overcurrent);
  check_err(status,__LINE__,__FILE__);
  
  varname.assign( "current_monitors");
  status = nc_inq_varid( outfile.gid[1], varname.c_str(), &varid);
  status = nc_put_var_float( outfile.gid[1], varid, current_monitors);
  check_err(status,__LINE__,__FILE__);
  
  varname.assign( "CCD_voltages");
  status = nc_inq_varid( outfile.gid[1], varname.c_str(), &varid);
  status = nc_put_var_float( outfile.gid[1], varid, ccd_voltages);
  check_err(status,__LINE__,__FILE__);
  
  varname.assign( "solenoid_voltage");
  status = nc_inq_varid( outfile.gid[1], varname.c_str(), &varid);
  status = nc_put_var_float( outfile.gid[1], varid, solenoid_voltage);
  check_err(status,__LINE__,__FILE__);

  // Write image paramter data as attributes
  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "errorCode",
                              NC_USHORT, 1, &streamInfo.imageInfo.errorCode);
  check_err(status,__LINE__,__FILE__);
  
  status = nc_put_att_uint( outfile.gid[1], NC_GLOBAL, "exposureID",
                            NC_UINT, 1, &streamInfo.imageInfo.exposureID);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_uint( outfile.gid[1], NC_GLOBAL, "imageID",
                            NC_UINT, 1, &streamInfo.imageInfo.imageID);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ulonglong( outfile.gid[1], NC_GLOBAL, "epochT0",
                                 NC_UINT64, 1, (const long long unsigned int*)
                                 &streamInfo.imageInfo.epochT0);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_uint( outfile.gid[1], NC_GLOBAL, "hostDeltaEpoch",
                            NC_UINT, 1, &streamInfo.imageInfo.hostDeltaEpoch);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_uint( outfile.gid[1], NC_GLOBAL, "hawkeyeDeltaEpoch",
                            NC_UINT, 1,
                            &streamInfo.imageInfo.hawkeyeDeltaEpoch);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "spectralBinning",
                              NC_USHORT, 1,
                              &streamInfo.imageInfo.spectralBinning);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "finderscopeBinning",
                              NC_USHORT, 1,
                              &streamInfo.imageInfo.finderscopeBinning);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "channelBitfield",
                              NC_USHORT, 1,
                              &streamInfo.imageInfo.channelBitfield);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "ccd1Exposure",
                              NC_USHORT, 1, &streamInfo.imageInfo.ccd1Exposure);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "ccd2Exposure",
                              NC_USHORT, 1, &streamInfo.imageInfo.ccd2Exposure);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "ccd3Exposure",
                              NC_USHORT, 1, &streamInfo.imageInfo.ccd3Exposure);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "ccd4Exposure",
                              NC_USHORT, 1, &streamInfo.imageInfo.ccd4Exposure);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "height",
                              NC_USHORT, 1, &streamInfo.imageInfo.height);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "darkHeight",
                              NC_USHORT, 1, &streamInfo.imageInfo.darkHeight);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "interval",
                              NC_USHORT, 1, &streamInfo.imageInfo.interval);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "oversampling",
                              NC_USHORT, 1, &streamInfo.imageInfo.oversampling);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "finderscopeExposure",
                              NC_USHORT, 1,
                              &streamInfo.imageInfo.finderscopeExposure);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "noFinderscopeImages",
                              NC_USHORT, 1,
                              &streamInfo.imageInfo.noFinderscopeImages);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ubyte( outfile.gid[1], NC_GLOBAL, "compression",
                             NC_UBYTE, 1,
                             (uint8_t *) &streamInfo.imageInfo.compression);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "darkSubtracted",
                              NC_USHORT, 1,
                              &streamInfo.imageInfo.darkSubtracted);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "shutterSolenoid",
                              NC_USHORT, 1,
                              &streamInfo.imageInfo.shutterSolenoid);
  check_err(status,__LINE__,__FILE__);

  status = nc_put_att_ushort( outfile.gid[1], NC_GLOBAL, "readoutOrder",
                              NC_USHORT, 1, &streamInfo.imageInfo.readoutOrder);
  check_err(status,__LINE__,__FILE__);
  
  delete[] ( telemetry);
  delete[] ( tlm_time);
  delete[] ( tlm_time_stamp);
  delete[] ( software_version);
  delete[] ( fpga_version); 
  delete[] ( telemetry_counter); 
  delete[] ( ccd_temperatures);
  delete[] ( fpga_temperature);
  delete[] ( fpga_voltages);
  delete[] ( overcurrent); 
  delete[] ( current_monitors); 
  delete[] ( ccd_voltages);
  delete[] ( solenoid_voltage); 

  // Date created
  strcpy( buf, unix2isodate(now(), 'G'));
  status = nc_put_att_text(outfile.ncid, NC_GLOBAL,
                                     "date_created", strlen(buf), buf);
  
  ////////////////////// Write EV bands //////////////////////
  double scan_time[MAXIMGHEIGHT];
  int32_t delta_time[MAXIMGHEIGHT];
  for (size_t i=0; i<MAXIMGHEIGHT; i++) {
    scan_time[i] = -999;
    delta_time[i] = -999;
  }
  
  int varidDark;
  size_t startDark[3]={0, 0, 0};
  size_t countDark[3]={1, 1, (size_t) rowOffset-2};
  status = nc_inq_varid( outfile.gid[3], "band_dark_pixels", &varidDark);

  for (size_t i=0; i<8; i++) {
    startDark[0] = i;

    int hasAverageDark = streamInfo.spectralInfo[i].averageDarkRowPresent;
    
    uint16_t imgHeight = streamInfo.spectralInfo[i].height + hasAverageDark;

    uint16_t imgWidth = 2 + streamInfo.spectralInfo[i].width;
    uint64_t epochT0 = streamInfo.imageInfo.epochT0;
    uint32_t pixelLen = imgHeight * imgWidth;

    //    cout << epochT0 << endl;
    
    //    uint32_t exposure;
    //if (i / 2 == 0)
    //exposure = streamInfo.imageInfo.ccd1Exposure * BAND_US_PER_EXPOSURE_COUNT;
    //else if (i / 2 == 1)
    //exposure = streamInfo.imageInfo.ccd2Exposure * BAND_US_PER_EXPOSURE_COUNT;
    //else if (i / 2 == 2)
    //exposure = streamInfo.imageInfo.ccd3Exposure * BAND_US_PER_EXPOSURE_COUNT;
    //else
    //exposure = streamInfo.imageInfo.ccd4Exposure * BAND_US_PER_EXPOSURE_COUNT;

    uint16_t *pixels = new uint16_t[pixelLen];
    for ( size_t j=0; j<pixelLen; j++) pixels[j] = 0;
    
    HawkeyeDecodeSpectralImage( i, pixels, pixelLen-imgWidth,
                                NULL, imgWidth);

    size_t start[2]={0, 0};
    size_t count[2]={1, (size_t) imgWidth-rowOffset};

    stringstream varstr;
    varstr << "band_" << setw(1) << i+1;
    status = nc_inq_varid( outfile.gid[3], varstr.str().c_str(), &varid);

    cout << "Writing Band: " << (i+1);
    cout << "  Image size: " << imgWidth-rowOffset << " by "
         << streamInfo.spectralInfo[i].height << endl;

    // skip band with invalid width read from corrupted raw data, LH, 11/02/2020
    if (imgWidth>2000) { // normal image width = 1818
        cout << "Image Width Error in Band " << (i+1)  << endl;
        exit(101);
    }
    
    // Note: Dark row is not written to L1A file
    for ( size_t j=0; j < (size_t) (imgHeight-1); j++) {
      if (i == 0) {
        delta_time[j] = pixels[j*imgWidth]*256*256 + pixels[j*imgWidth+1];
        if ((delta_time[j]==0) && (j>0) ) continue;  // Flagging scan_time if delta_time is 0, LH, 02/17/2022 ; avoid skipping first record, to get time_coverage_start value, LH 3/25/2022
        scan_time[j] = epochT0 + delta_time[j];
        scan_time[j] *= 0.001;
        scan_time[j] -= (53*24*3600 + 14788.224);
        tepoch2yds( scan_time[j], &iyr, &idy, &sec);
        scan_time[j] = sec;
        if (bImageCrossDay && (sec<80000)) scan_time[j] = sec + 86400.0;
        if (j == 0) {
          int16_t year, doy;
          year = iyr;
          doy = idy;
          double myUnixTime= yds2unix(year, doy, sec);
          strcpy( buf, unix2isodate(myUnixTime, 'G'));
          status = nc_put_att_text(outfile.ncid, NC_GLOBAL,
                                     "time_coverage_start", strlen(buf), buf);
        }
      }

      start[0] = j;
      status = nc_put_vara_ushort( outfile.gid[3], varid, start, count, 
                                   &pixels[rowOffset+j*imgWidth]);
      check_err(status,__LINE__,__FILE__);

      startDark[1] = j;
      status = nc_put_vara_ushort( outfile.gid[3], varidDark,
                                   startDark, countDark, 
                                   &pixels[2+j*imgWidth]);
      check_err(status,__LINE__,__FILE__);
    }

    int16_t year, doy;
    year = iyr;
    doy = idy;
    double myUnixTime= yds2unix(year, doy, sec);
    strcpy( buf, unix2isodate(myUnixTime, 'G'));
    status = nc_put_att_text(outfile.ncid, NC_GLOBAL,
                                     "time_coverage_end", strlen(buf), buf);

    delete[] pixels;
  } // Band loop

  varname.assign( "scan_time");
  status = nc_inq_varid( outfile.gid[0], varname.c_str(), &varid);
  status = nc_put_var_double( outfile.gid[0], varid, scan_time);
  check_err(status,__LINE__,__FILE__);

  varname.assign( "delta_time");
  status = nc_inq_varid( outfile.gid[0], varname.c_str(), &varid);
  status = nc_put_var_int( outfile.gid[0], varid, delta_time);
  check_err(status,__LINE__,__FILE__);


  ////////////////////// Write Finder images //////////////////////
  //uint16_t finderHeight = FINDERSCOPE_LIGHT_HEIGHT;
  //uint16_t finderWidth =  FINDERSCOPE_LIGHT_WIDTH;

  uint16_t finderHeight = 0; // streamInfo.finderscopeInfo[0].height;  // ver 1.0.1
  uint16_t finderWidth =  0; // streamInfo.finderscopeInfo[0].width;
  // find max height and width of finderscope images [MAX_FINDERSCOPE_IMAGES]
  for (size_t i=0; i<FINDERSCOPE_MAX_IMAGES; i++) {
  	if (streamInfo.finderscopeInfo[i].height > 0) finderHeight = 480;
  	if (streamInfo.finderscopeInfo[i].width > 0) finderWidth = 756;
  }
  uint32_t finderPixelLen = finderHeight*finderWidth;

  if ( finderPixelLen == 0) {
    cout << "Zero Finder Pixel Length" << endl;
    outfile.close();
    return 0;
  }
  
  uint16_t *finderPixels = new uint16_t[finderPixelLen];
  for (size_t i=0; i<finderPixelLen; i++) finderPixels[i] = 0;
    
  size_t start[3]={0, 0, 0};
  size_t count[3]={1, finderHeight, finderWidth};

  double finder_time[FINDERSCOPE_MAX_IMAGES];
  int32_t finder_delta_time[FINDERSCOPE_MAX_IMAGES];
  uint64_t epochT0 = streamInfo.imageInfo.epochT0;
  
  status = nc_inq_varid( outfile.gid[3], "finder", &varid);

  cout << endl << "Writing Finder Images";
  cout << "  Finder size: " << finderWidth << " by " << finderHeight << endl;
  int finderScopeStatus;
  // for (size_t i=0; i< (size_t) streamInfo.noFinderscopeImages; i++) {
  // There are abnormal cases when finder scope images are stored in sequence as,
  // [5,6,7,8,9,10,11,15,16,17,18,19,2]
  for (size_t i=0; i< FINDERSCOPE_MAX_IMAGES; i++) {
    finder_time[i] = -999;
    if (streamInfo.finderscopeInfo[i].width==0) continue;
    
    finderScopeStatus =
      HawkeyeDecodeFinderscopeImage(i, finderPixels, finderPixelLen);

    // if ( finderScopeStatus == HSE_NO_HEADER_FOUND) break;
    if ( finderScopeStatus == HSE_NO_HEADER_FOUND) continue;

    // We consolidate the finderPixels buffer by removing the first
    // n dark pixels from each row and then zeroing out the remaining rows.

    //    for (size_t j=0; j<(752*480/756); j++) {
    // memmove(&finderPixels[j*FINDERSCOPE_LIGHT_WIDTH],
    //        &finderPixels[j*(FINDERSCOPE_LIGHT_WIDTH+4)+4],
    //        FINDERSCOPE_LIGHT_WIDTH*2);
    //}
    //memset(&finderPixels[(752*480/756)*FINDERSCOPE_LIGHT_WIDTH],
    //     0, (finderPixelLen-(752*480/756)*FINDERSCOPE_LIGHT_WIDTH)*2);
    
    start[0] = i;
    status = nc_put_vara_ushort( outfile.gid[3], varid, start, count, 
                                 finderPixels);
    check_err(status,__LINE__,__FILE__);

    finder_delta_time[i] =  streamInfo.finderscopeInfo[i].timeStamp;
    finder_time[i] = epochT0 + finder_delta_time[i];
    finder_time[i] *= 0.001;
    finder_time[i] -= (53*24*3600 + 14788.224);
    tepoch2yds( finder_time[i], &iyr, &idy, &sec);
    finder_time[i] = sec;
    if (bImageCrossDay && (sec<80000)) finder_time[i] = sec + 86400.0;

  }
  delete[] finderPixels;

  if ( finderScopeStatus == HSE_NO_ERROR) {
    varname.assign( "finder_time");
    status = nc_inq_varid( outfile.gid[0], varname.c_str(), &varid);
    status = nc_put_var_double( outfile.gid[0], varid, finder_time);
    check_err(status,__LINE__,__FILE__);

    varname.assign( "finder_delta_time");
    status = nc_inq_varid( outfile.gid[0], varname.c_str(), &varid);
    status = nc_put_var_int( outfile.gid[0], varid, finder_delta_time);
    check_err(status,__LINE__,__FILE__);
  }

  ////////////////////// Write Navigation //////////////////////

  size_t startNav[2]={0, 0};
  size_t countNav[2]={1, 1};

  int offadcs;
  double min;
  size_t ngood;

  // ATTITUDE
  offadcs = 0;
  min = attitude[0]->sec;
  for (size_t i=1; i<n[ATTITUDE]; i++) {
    //    cout << i << " " << attitude[i]->sec << endl;
    if (attitude[i]->sec < min) {
      offadcs = i;
      min = attitude[i]->sec;
    }
  }

  ngood = nSC;
  for (size_t i=offadcs+1; i<n[ATTITUDE]; i++) {
    if ((attitude[i]->sec-attitude[i-1]->sec) < 0) {
      ngood = i - offadcs;
      break;
    }
  }
  if (ngood > nSC) ngood = nSC;
  if (n[ATTITUDE] == 0) {
    ngood = 0;  // Set to 0 if no attitude records
    cout << "No attitude records written" << endl;
  }
  
  varname.assign( "att_time");
  status = nc_inq_varid( outfile.gid[2], varname.c_str(), &varid);
  for (size_t i=0; i<ngood; i++) {
    startNav[0] = i;
    status = nc_put_vara_double( outfile.gid[2], varid, startNav, countNav,
                                 &attitude[i+offadcs]->sec);
    check_err(status,__LINE__,__FILE__);
  }

  varname.assign( "att_quat");
  status = nc_inq_varid( outfile.gid[2], varname.c_str(), &varid);
  for (size_t i=0; i<ngood; i++) {
    startNav[0] = i;
    countNav[1] = 4;
    status = nc_put_vara_float( outfile.gid[2], varid, startNav, countNav,
                                &attitude[i+offadcs]->quat[0]);
    check_err(status,__LINE__,__FILE__);
  }

  // PROPAGATOR
  offadcs = 0;
  min = propagator[0]->sec;
  for (size_t i=1; i<n[PROPAGATOR]; i++) {
    //    cout << i << " " << propagator[i]->sec << endl;
    if (propagator[i]->sec < min) {
      offadcs = i;
      min = propagator[i]->sec;
    }
  }
  
  ngood = nSC;
  for (size_t i=offadcs+1; i<n[PROPAGATOR]; i++) {
    if ((propagator[i]->sec-propagator[i-1]->sec) < 0) {
      ngood = i - offadcs;
      break;
    }
  }
  if (ngood > nSC) ngood = nSC;
  if (n[PROPAGATOR] == 0) {
    ngood = 0;  // Set to 0 if no propagator records
    cout << "No progagator records written" << endl;
  }
  
  varname.assign( "orb_time");
  status = nc_inq_varid( outfile.gid[2], varname.c_str(), &varid);
  for (size_t i=0; i<ngood; i++) {
    startNav[0] = i;
    status = nc_put_vara_double( outfile.gid[2], varid, startNav, countNav,
                                 &propagator[i+offadcs]->sec);
    check_err(status,__LINE__,__FILE__);
  }
  
  varname.assign( "orb_pos");
  status = nc_inq_varid( outfile.gid[2], varname.c_str(), &varid);
  for (size_t i=0; i<ngood; i++) {
    startNav[0] = i;
    countNav[1] = 3;
    status = nc_put_vara_float( outfile.gid[2], varid, startNav, countNav,
                                &propagator[i+offadcs]->pos[0]);
    check_err(status,__LINE__,__FILE__);
  }

  varname.assign( "orb_vel");
  status = nc_inq_varid( outfile.gid[2], varname.c_str(), &varid);
  for (size_t i=0; i<ngood; i++) {
    startNav[0] = i;
    countNav[1] = 3;
    status = nc_put_vara_float( outfile.gid[2], varid, startNav, countNav,
                                &propagator[i+offadcs]->vel[0]);
    check_err(status,__LINE__,__FILE__);
  }


  // SENSOR
  offadcs = 0;
  min = sensor[0]->sec;
  for (size_t i=1; i<n[SENSOR]; i++) {
    //    cout << i << " " << sensor[i]->sec << endl;
    if (sensor[i]->sec < min) {
      offadcs = i;
      min = sensor[i]->sec;
    }
  }
  
  ngood = nSC;
  for (size_t i=offadcs+1; i<n[SENSOR]; i++) {
    if ((sensor[i]->sec-sensor[i-1]->sec) < 0) {
      ngood = i - offadcs;
      break;
    }
  }
  if (ngood > nSC) ngood = nSC; // LH, 7/12/2021
  if (ngood > n[SENSOR]) ngood = n[SENSOR];  // LH, 01/15/2021
  if (ngood == 0) {
    cout << "No sensor records written" << endl;
  }

  varname.assign( "sensor_time");
  status = nc_inq_varid( outfile.gid[2], varname.c_str(), &varid);
  for (size_t i=0; i<ngood; i++) {
    startNav[0] = i;
    status = nc_put_vara_double( outfile.gid[2], varid, startNav, countNav,
                                 &sensor[i+offadcs]->sec);
    check_err(status,__LINE__,__FILE__);
  }

  varname.assign( "sensor_bus_telemetry");
  status = nc_inq_varid( outfile.gid[2], varname.c_str(), &varid);
  for (size_t i=0; i<ngood; i++) {
    startNav[0] = i;
    countNav[1] = 27;
    status = nc_put_vara_ubyte( outfile.gid[2], varid, startNav, countNav,
                                &sensor[i+offadcs]->bus_telemetry[0]);
    check_err(status,__LINE__,__FILE__);
  }


  // PSENSOR
  offadcs = 0;
  min = psensor[0]->sec;
  for (size_t i=1; i<n[PSENSOR]; i++) {
    //    cout << i << " " << psensor[i]->sec << endl;
    if (psensor[i]->sec < min) {
      offadcs = i;
      min = psensor[i]->sec;
    }
  }
  
  ngood = nSC;
  for (size_t i=offadcs+1; i<n[PSENSOR]; i++) {
    if ((psensor[i]->sec-psensor[i-1]->sec) < 0) {
      ngood = i - offadcs;
      break;
    }
  }
  if (ngood > nSC) ngood = nSC;
  if (n[PSENSOR] == 0) {
    ngood = 0;  // Set to 0 if no psensor records
    cout << "No psensor records written" << endl;
  }
  
  varname.assign( "processed_sensor_time");
  status = nc_inq_varid( outfile.gid[2], varname.c_str(), &varid);
  for (size_t i=0; i<ngood; i++) {
    startNav[0] = i;
    status = nc_put_vara_double( outfile.gid[2], varid, startNav, countNav,
                                 &psensor[i+offadcs]->sec);
    check_err(status,__LINE__,__FILE__);
  }

  varname.assign( "gyro_rates");
  status = nc_inq_varid( outfile.gid[2], varname.c_str(), &varid);
  for (size_t i=0; i<ngood; i++) {
    startNav[0] = i;
    countNav[1] = 3;
    status = nc_put_vara_float( outfile.gid[2], varid, startNav, countNav,
                                &psensor[i+offadcs]->gyro[0]);
    check_err(status,__LINE__,__FILE__);
  }

  varname.assign( "mag1");
  status = nc_inq_varid( outfile.gid[2], varname.c_str(), &varid);
  for (size_t i=0; i<ngood; i++) {
    startNav[0] = i;
    countNav[1] = 3;
    status = nc_put_vara_float( outfile.gid[2], varid, startNav, countNav,
                                &psensor[i+offadcs]->mag1[0]);
    check_err(status,__LINE__,__FILE__);
  }

  varname.assign( "mag2");
  status = nc_inq_varid( outfile.gid[2], varname.c_str(), &varid);
  for (size_t i=0; i<ngood; i++) {
    startNav[0] = i;
    countNav[1] = 3;
    status = nc_put_vara_float( outfile.gid[2], varid, startNav, countNav,
                                &psensor[i+offadcs]->mag2[0]);
    check_err(status,__LINE__,__FILE__);
  }

  varname.assign( "rwheels");
  status = nc_inq_varid( outfile.gid[2], varname.c_str(), &varid);
  for (size_t i=0; i<ngood; i++) {
    startNav[0] = i;
    countNav[1] = 3;
    status = nc_put_vara_float( outfile.gid[2], varid, startNav, countNav,
                                &psensor[i+offadcs]->rwheels[0]);
    check_err(status,__LINE__,__FILE__);
  }

  varname.assign( "sun_vector");
  status = nc_inq_varid( outfile.gid[2], varname.c_str(), &varid);
  for (size_t i=0; i<ngood; i++) {
    startNav[0] = i;
    countNav[1] = 3;
    status = nc_put_vara_float( outfile.gid[2], varid, startNav, countNav,
                                &psensor[i+offadcs]->sunb[0]);
    check_err(status,__LINE__,__FILE__);
  }

  outfile.close();

  return 0;
}


/*----------------------------------------------------------------- */
/* Create an Generic NETCDF4 level1 file                            */
/* ---------------------------------------------------------------- */
int l1aFile::createl1( char* l1_filename, uint32_t nSC,
                       uint32_t imgWidth, uint32_t imgHeight,
                       uint32_t fndWidth, uint32_t fndHeight) {
  
   int status;
   
   status = nc_create( l1_filename, NC_NETCDF4, &ncid);
   check_err(status,__LINE__,__FILE__);

   ifstream hawkeye_l1a_data_structure;
   string line;
   string dataStructureFile;

   dataStructureFile.assign("$OCDATAROOT/hawkeye/Hawkeye_Level-1A_Data_Structure.cdl");
   expandEnvVar( &dataStructureFile);

   hawkeye_l1a_data_structure.open( dataStructureFile.c_str(), ifstream::in);
   if ( hawkeye_l1a_data_structure.fail() == true) {
     cout << "\"" << dataStructureFile.c_str() << "\" not found" << endl;
     exit(1);
   }

   // Find "dimensions" section of CDL file
   while(1) {
     getline( hawkeye_l1a_data_structure, line);
     size_t pos = line.find("dimensions:");
     if ( pos == 0) break;
   }

   // Define dimensions from "dimensions" section of CDL file
   ndims = 0;
   //int32_t numScans=0;
   while(1) {
     getline( hawkeye_l1a_data_structure, line);
     size_t pos = line.find(" = ");
     if ( pos == string::npos) break;

     uint32_t dimSize;
     istringstream iss(line.substr(pos+2, string::npos));
     iss >> dimSize;

     iss.clear(); 
     iss.str( line);
     iss >> skipws >> line;

     //     cout << "Dimension Name: " << line.c_str() << " Dimension Size: "
     //	  << dimSize << endl;

     //if (line.compare("number_of_scans") == 0 && numScans > 0) {
     //  dimSize = numScans;
     //}

     if (line.compare("number_of_SC_records") == 0) {
       dimSize = nSC;
     }

     if (line.compare("number_of_scans") == 0) {
       dimSize = imgHeight;
     }

     if (line.compare("number_of_pixels") == 0) {
       dimSize = imgWidth;
     }

     if (line.compare("finder_lines") == 0) {
       dimSize = fndHeight;
     }

     if (line.compare("finder_pixels") == 0) {
       dimSize = fndWidth;
     }

     status = nc_def_dim( ncid, line.c_str(), dimSize, &dimid[ndims++]);
     check_err(status,__LINE__,__FILE__);
     
   } // while loop

   // Read global attributes (string attributes only) 
   while(1) {
     getline( hawkeye_l1a_data_structure, line);
     size_t pos = line.find("// global attributes");
     if ( pos == 0) break;
   }

   while(1) {
     getline( hawkeye_l1a_data_structure, line);
     size_t pos = line.find(" = ");
     if ( pos == string::npos) break;

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

     status = nc_put_att_text(ncid, NC_GLOBAL, attName.c_str(), 
                              strlen(attValue.c_str()), attValue.c_str());
     check_err(status,__LINE__,__FILE__);

   } // while(1)


   ngrps = 0;
   // Loop through groups
   while(1) {
     getline( hawkeye_l1a_data_structure, line);

     // Check if end of CDL file
     // If so then close CDL file and return
     if (line.substr(0,1).compare("}") == 0) {
       hawkeye_l1a_data_structure.close();
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
       status = nc_def_grp( ncid, line.c_str(), &gid[ngrps]);
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
       getline( hawkeye_l1a_data_structure, line); // skip "variables:"
       while(1) {
         getline( hawkeye_l1a_data_structure, line);

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
           //cout << "sname: " << sname.c_str() << endl;

           // Parse variable dimension info
           this->parseDims( line.substr(pos+1, string::npos), 
                      &numDims, varDims);
           for (int i=0; i<numDims; i++) { 
             nc_inq_dim( ncid, varDims[i], dimName, &dimSize[i]);
             //cout << line.c_str() << " " << i << " " << dimName
             //	  << " " << dimSize[i] << endl;
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

         } // if ( pos == string::npos)
       } // datasets in group loop
     } // New Group loop
   } // Main Group loop
   
   return 0;
}


int l1aFile::parseDims( string dimString, int *numDims, int *varDims) {

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


int l1aFile::close() {

  int status = nc_close(ncid);

  if (status != NC_NOERR)
  {
     printf("-E - nc_close failed for ncid: %i\n",ncid);
  }
  return 0;
}


int unpack_seahawk_adcs( uint8_t *apkt, double startHWKTime, double stopHWKTime,
                         sensor_t *sensor, psensor_t *psensor,
                         attitude_t *attitude, propagator_t *propagator) {

  int16_t i16;
  int32_t i32;
  int32_t iyr, idy;
  int32_t startdy, stopdy;
  double sec;

  short psub = apkt[8];

  // Get time

  // JD 2015/01/01 = 2457023.5
  // JD 1980/01/06 = 2444244.5
  // delta sec = (2457023.5-2444244.5) * 24 * 3600 = 1104105600
  // Add 16 leapsecs between 1980/01/06 and 2015/01/01 
  
  memcpy(&i32, &apkt[9], 4);
  double tepoch = (double) __builtin_bswap32(i32) - (1104105600+16);

  //  cout << tepoch - startHWKTime << endl;
  if ( (startHWKTime-tepoch) > 2 || (tepoch - stopHWKTime) > 1) {
    return -1;
  }
  
  tepoch2yds( startHWKTime, &iyr, &idy, &sec);
  startdy = iyr*1000 + idy;
  tepoch2yds( stopHWKTime, &iyr, &idy, &sec);
  stopdy = iyr*1000 + idy;
    
  tepoch2yds( tepoch, &iyr, &idy, &sec);
  memcpy(&i32, &apkt[13], 4);
  sec += ((double) __builtin_bswap32(i32)) / 4294967296.0;
  memcpy(&i16, &apkt[17], 2);
  sec -= ((double) SWAP_2(i16)) / 1000;
  // if the image spans over 00:00, the sec of day for earlier part of image is from previous day, 
  // ending with a value close to 86400
  if ((stopdy>startdy) && (sec<80000)) sec += 86400.0; 
  
  switch (psub) {
   case 1:
     sensor->iyr = iyr;
     sensor->iday = idy;
     sensor->sec = sec;
     memcpy(sensor->bus_telemetry, &apkt[19], sizeof(sensor->bus_telemetry));
     break;

  case 2:
    psensor->iyr = iyr;
    psensor->iday = idy;
    psensor->sec = sec;
    for (size_t i=0; i<3; i++) {
      memcpy(&i32, &apkt[19+4*i], 4);
      i32 = __builtin_bswap32(i32);
      memcpy(&psensor->gyro[i], &i32, 4);

      memcpy(&i32, &apkt[31+4*i], 4);
      i32 = __builtin_bswap32(i32);
      memcpy(&psensor->mag1[i], &i32, 4);
      
      memcpy(&i32, &apkt[43+4*i], 4);
      i32 = __builtin_bswap32(i32);
      memcpy(&psensor->mag2[i], &i32, 4);

      memcpy(&i32, &apkt[55+4*i], 4);
      i32 = __builtin_bswap32(i32);
      memcpy(&psensor->rwheels[i], &i32, 4);

      memcpy(&i32, &apkt[67+4*i], 4);
      i32 = __builtin_bswap32(i32);
      memcpy(&psensor->sunb[i], &i32, 4);
    }
    break;

  case 3:
    attitude->iyr = iyr;
    attitude->iday = idy;
    attitude->sec = sec;
    for (size_t i=0; i<4; i++) {
      memcpy(&i32, &apkt[19+4*i], 4);
      i32 = __builtin_bswap32(i32);
      memcpy(&attitude->quat[i], &i32, 4);
    }
    break;

  case 4:
    propagator->iyr = iyr;
    propagator->iday = idy;
    propagator->sec = sec;
    for (size_t i=0; i<3; i++) {
      memcpy(&i32, &apkt[19+4*i], 4);
      i32 = __builtin_bswap32(i32);
      memcpy(&propagator->pos[i], &i32, 4);

      memcpy(&i32, &apkt[31+4*i], 4);
      i32 = __builtin_bswap32(i32);
      memcpy(&propagator->vel[i], &i32, 4);
    }
    break;
  }
  
  return 0; 
}


