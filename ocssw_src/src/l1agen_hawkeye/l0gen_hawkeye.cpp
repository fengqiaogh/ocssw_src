#include <stdint.h>
#include <libgen.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <getopt.h>

#include "hawkeyeUtil.h"

#define VERSION "1.0.0_2023-05-10"

//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     FutureTech     08/08/18 0.10  Original development
//  Joel Gales     FutureTech     08/30/18 0.11  Output L0 filenames generated
//                                               from image start time, Add
//                                               support for text output files
//
//  Joel Gales     FutureTech     09/20/18 0.20  Add duration to textfile,
//                                               generate output name from
//                                               starttime
//  Joel Gales     FutureTech     09/27/18 0.30  Add Start/Stop HWK time to
//                                               textfile (Delete duration).
//                                               Pass through StartTime code
//                                               only once.
//  Joel Gales     SAIC           10/26/18 0.40  Code cleanup, Prepend tlm
//                                               filename, exposure ID, and
//                                               image counter to HWK L0 name
//  Joel Gales     SAIC           11/16/18 0.50  Add support for GPS telemetry
//                                               output telemetry file
//  Joel Gales     SAIC           03/07/19 0.55  Fix error in last L0 outfile
//                                               name
//  Joel Gales     SAIC           05/23/19 0.60  Overhaul to handle frame drops
//  Joel Gales     SAIC           06/19/19 0.70  Force new frame read if bad
//                                               packet header
//  Joel Gales     SAIC           07/09/19 0.71  Trap negative dataLen
//                                               filename
//  Joel Gales     SAIC           07/29/19 0.80  Trap bad framePtr value in
//                                               getSHpacket.  Check for
//                                               output filesize
//  Liang Hong     SAIC           07/13/21 0.82  Handle EOI packet crossing frames
//  Liang Hong     SAIC           08/06/21 0.83  Handle missing EOI packet images
//  Liang Hong     SAIC           08/12/21 0.84  Updated frame pointer initialization,
//                                               packet reading, and EOI detections
//  Liang Hong     SAIC           09/01/21 0.85  Updated End of Exposure detection
//  Liang Hong     SAIC           05/10/23 1.00  Add verbose option

using namespace std;

#define DATAFRAMEBYTES 892

int getZuluTime_L0Name( double tlm_time,
                        string *zuluTime, string *packetOutName);

#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )

int main (int argc, char* argv[])
{
  cout << "l0gen_hawkeye " << VERSION << " (" 
       <<  __DATE__ << " " << __TIME__ << ")" << endl;

  if ( argc == 1) {
    cout << endl << 
      "l0gen_hawkeye_new input_telemetry_filename" << endl;

    return 0;
  }
  
  ifstream framefile;
  framefile.open( argv[1], ios::binary );
  //  int lastpos = framefile.tellg();
  
  bool isVerbose = false;
  int c;
  while (1) {
    static struct option long_options[] = {
      {"Verbose", no_argument,0,'v'}
    };

    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    c = getopt_long( argc, argv, "v", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
      {
      case 'v':
        isVerbose = true; 
        break;
        
      default:
        abort ();
      }
  }
  
  ofstream packetOut;
  ofstream gpsOut;
  ofstream textOut;
  
  string packetOutname="TBD";
  string zuluStartTime, zuluStopTime;
  string prePend;

  double startHWKTime=0.0;
  double stopHWKTime=0.0;

  uint32_t exposureID;
  //uint32_t deltaTime;
  //uint32_t imageRow;
  //uint32_t finderNumber;
  
  uint32_t prevExposureID=999999999;
  //uint32_t prevRow;
  //uint32_t prevFinder;
  //int prevDeltaTime;
  
  int packetLength;
  
  // get length of input framefile:
  // framefile.seekg (0, ios::end);
  // uint32_t granuleLength = framefile.tellg();
  // framefile.seekg (0, ios::beg);

  uint8_t **packet = new uint8_t*;
  
  //  int write=0;
  int firstGPS=1;
  
  //  int imageCounter=0;

  uint8_t frame[DATAFRAMEBYTES];

  // Skip fill records and find first good packet
  int framePtr = 2046;
  while (framePtr == 2046 || framePtr == 2047) {
    framefile.read( (char *) frame, DATAFRAMEBYTES);
    framePtr = (frame[4] % 8)*256 + frame[5];
  }
  framePtr += 6;
  
  int frameDrop = 0;
  
  int prevFrameCnt = frame[2];

  //  int framePtr = 6;

  // packet[7]
  // ---------
  // 253 -- Image data
  // 254 -- GPS data
  // 255 -- ADCS data

  // Image Data (packet[8])
  // ----------
  // 1 -- Uncompressed Image Data
  // 2 -- Compressed Image Data
  // 3 -- Image parameters
  // 4 -- Telemetry
  // 5 -- End of File
  // 6 -- Mission Log

  packetOut.open("tempL0_0", ios::out | ios::trunc | ios::binary);
  int tempCnt = 0;
  int exposureCount = 0;
  while (1) {
  	frameDrop = 0; // V0.84

    int result_getSHpacket = getSHpacket( &framefile, frame, framePtr, packet, packetLength,
                 prevFrameCnt, frameDrop, isVerbose);
    
	// V0.82, End of exposure packet detected at frame boundary
	if (result_getSHpacket==1){	
		cout<<"exposureID="<<exposureID<<endl;
		cout << "end of exposure" << endl;
		packetOut.close();

		tempCnt++;
		string tempName;
		tempName.assign("tempL0_");
		tempName.append(to_string(tempCnt));
		packetOut.open(tempName.c_str(), ios::out | ios::trunc | ios::binary);

		//prevRow = 0;
		//prevFinder = 0;

		framefile.read( (char *) frame, DATAFRAMEBYTES);

		// int priPacketLen = (frame[4] % 8) * 256 + frame[5];
		framePtr = (frame[4] % 8) * 256 + frame[5];
		while (framePtr == 2046 || framePtr == 2047) {
			framefile.read( (char *) frame, DATAFRAMEBYTES);
			framePtr = (frame[4] % 8) * 256 + frame[5];
		}
		prevFrameCnt = frame[2];
		framePtr += 6;	// V0.84
		prevExposureID = 999999999;
		continue;
	} else {
		if (frameDrop == -1)  continue; 

		long long int pos = framefile.tellg();
		//    cout << "pos: " << pos << endl;
		if ( framefile.eof() || pos == -1) break;

		// frame drop
		if (frameDrop == 1) {
		  frameDrop = 2;
		  continue;
		}

		if ((frame[4] % 8)*256 + frame[5] == 2046) {
		  continue;
		}
	}

    //////////////////////////////////////////////
    ////////// End of Exposure Packet ////////////
    //////////////////////////////////////////////
    //int apid = ((*packet)[0] % 8) * 256 + (*packet)[1];  // V0.85
    //if (apid == 2047) {									// V0.85
    if ((*packet)[0]==7 && (*packet)[1]==255) {				// updated End of Exposure detection, V0.85
	  cout<<"exposureID="<<exposureID<<endl; // V0.82
      cout << "end of exposure" << endl;
      packetOut.close();

      tempCnt++;
      string tempName;
      tempName.assign("tempL0_");
      tempName.append(to_string(tempCnt));
      packetOut.open(tempName.c_str(), ios::out | ios::trunc | ios::binary);

      //prevRow = 0;
      //prevFinder = 0;
      
      framefile.read( (char *) frame, DATAFRAMEBYTES);

      // int priPacketLen = (frame[4] % 8) * 256 + frame[5];
      framePtr = (frame[4] % 8) * 256 + frame[5];
      while (framePtr == 2046 || framePtr == 2047) {
        framefile.read( (char *) frame, DATAFRAMEBYTES);
        framePtr = (frame[4] % 8) * 256 + frame[5];
      }
      prevFrameCnt = frame[2];
      framePtr += 6;  // V0.84
      prevExposureID = 999999999; 
      continue;
    }


    //////////////////////////////////////////////
    ///////////// Parameter Packet ///////////////
    //////////////////////////////////////////////
    if ((*packet)[7] == 253 && (*packet)[8] == 3) {

      uint64_t ll;
      uint32_t ul;
      
      // Decode Exposure ID
      ul =  (uint32_t) (*packet)[15] << 24;
      ul += (uint32_t) (*packet)[16] << 16;
      ul += (uint32_t) (*packet)[17] << 8;
      ul += (*packet)[18];
      exposureID = ul;
      cout<<"read exposureID="<<exposureID<<endl; 
        
      // Decode Epoch Time
      // Code taken from DecodeImageParams() in HawkeyeDecode.c
      ll =  (uint64_t) (*packet)[26] << 40;
      ll += (uint64_t) (*packet)[27] << 32;
      ll += (uint64_t) (*packet)[28] << 24;
      ll += (uint64_t) (*packet)[29] << 16;
      ll += (uint64_t) (*packet)[30] << 8;
      ll += (*packet)[31];
      ll = (ll & 0x0FFFFFFFFFF0) >> 4;

      //      double tlm_time;

      //tlm_time = ll * 0.001;
      //tlm_time -= (53*24*3600 + 14788.224);
      //cout << "imag tepoch: " << setprecision(11) << setw(15) << tlm_time << endl;

	  if (exposureID<310000000 && exposureID>180000000) {	// V0.85, avoid bogus exposure id
      //if (frameDrop == 2) {
           if (prevExposureID != 999999999 && prevExposureID != exposureID)   frameDrop = 3;
      //}
       		prevExposureID = exposureID;
       }
    }
      

    //////////////////////////////////////////////
    /////////////// Image Packet /////////////////
    //////////////////////////////////////////////
    if ((*packet)[7] == 253 && (*packet)[8] == 2) {
     uint16_t row;
     memcpy(&row, &(*packet)[20], 2);
     row = SWAP_2(row);

     //uint8_t infoByte = (*packet)[15];
     //uint8_t sb = infoByte & 0x1f;

     // sd_f: 0 - finder, 1 - spectral (image)
     //uint8_t sd_f = (infoByte & 0x40) >> 6;
     //     uint8_t dark = (infoByte & 0x80) >> 7;

     //uint32_t ul;
      
     //Decode Delta-Hawkeye time
     //ul =  (uint32_t) (*packet)[16] << 24;
     //ul += (uint32_t) (*packet)[17] << 16;
     //ul += (uint32_t) (*packet)[18] << 8;
     //ul += (*packet)[19];
     //deltaTime = (ul >> 4) & 0xFFFFFF;

/*    V0.84
     if (frameDrop == 2) {
       if (sd_f == 1 && row < prevRow)
         frameDrop = 3;
       if (sd_f == 0 && sb < prevFinder)
         frameDrop = 3;
       //       if (sd_f == 1 && (((int) deltaTime) - prevDeltaTime) > 100*1000)
       // frameDrop = 3;
     }

     if (sd_f == 1) {
       //       cout << "band, row: " << (int) sb << " " << (int) row << endl;
       prevRow = row;
     } else {
       //cout << "find, row: " << (int) sb << " " << (int) row << endl;
       prevFinder = sb;
     }
     //if (sd_f == 1) prevDeltaTime = deltaTime;
*/  

    }


    //////////////////////////////////////////////
    //////////////// GPS Packet //////////////////
    //////////////////////////////////////////////
    if ((*packet)[7] == 254) {

      //    cout << "GPS Packet" << (int) (*packet)[8] << " " << packetLength << endl;
      
      if (firstGPS) {
        string gpsName;
        gpsName.assign(basename(argv[1]));
        size_t strpos = gpsName.find(".tlm");
        gpsName = gpsName.substr(0, strpos);
        gpsName.append(".gps"); 

        gpsOut.open(gpsName.c_str(), ios::out | ios::trunc | ios::binary);
        firstGPS = 0;
      }
      
      gpsOut.write( (char *) &(*packet)[0], 161);

      
      //////////////////////////////////////////////
      /////////////// ADCS Packet //////////////////
      //////////////////////////////////////////////
    } else if ((*packet)[7] == 255) {
      // ADCS packet
      //      cout << "ADCS packet: " << (int) (*packet)[3] << " " <<
      // (int) (*packet)[8] << " " << packetLength << endl;


      int16_t i16;
      int32_t i32;
      int32_t iyr, idy;
      double sec;

      // Get time

      // JD 2015/01/01 = 2457023.5
      // JD 1980/01/06 = 2444244.5
      // delta sec = (2457023.5-2444244.5) * 24 * 3600 = 1104105600
      // Add 16 leapsecs between 1980/01/06 and 2015/01/01 
  
      memcpy(&i32, &(*packet)[9], 4);
      double tepoch = (double) __builtin_bswap32(i32) - (1104105600+16);
      //      cout << "adcs tepoch: " << setprecision(11) << setw(20) << tepoch << endl;
      tepoch2yds( tepoch, &iyr, &idy, &sec);
      memcpy(&i32, &(*packet)[13], 4);
      sec += ((double) __builtin_bswap32(i32)) / 4294967296.0;
      memcpy(&i16, &(*packet)[17], 2);
      sec -= ((double) SWAP_2(i16)) / 1000;
      //cout << "adcs sec: " << sec << endl;
    }


    if ( frameDrop == 3) {
      packetOut.close();
	  cout<<"exposureID="<<exposureID<<endl; // V0.82
      cout << "Close image file" << endl;     // V0.82
	  
      tempCnt++;
      string tempName;
      tempName.assign("tempL0_");
      tempName.append(to_string(tempCnt));
      
      packetOut.open(tempName.c_str(), ios::out | ios::trunc | ios::binary);

      //prevRow = 0;
      //prevFinder = 0;
      
      frameDrop = 0;
      // prevExposureID = 999999999; V0.85
    }

        
    if ((*packet)[7] != 0) {
      packetOut.write( (char *) (*packet), packetLength);
      delete[] (*packet);
    }
  } // end main loop

  delete packet;

  tempCnt++;
  for (int i=0; i<tempCnt; i++) {
    string tempName;
    tempName.assign("tempL0_");
    tempName.append(to_string(i));

    ifstream tempFile;
    uint8_t sndpus[9];
    uint8_t *data;
    
    tempFile.open( tempName.c_str(), ios::binary );

    // Check file size -- Skip if 4503 bytes
    int fsize = tempFile.tellg();
    tempFile.seekg( 0, std::ios::end);
    fsize = (int) tempFile.tellg() - fsize;
    //    cout << fsize << endl;
    if ( fsize == 4503 || fsize <= 10000000) {
      remove( tempName.c_str());
      continue;
    }
    
    tempFile.seekg( 0, std::ios::beg);

    while (!tempFile.eof()) {
      tempFile.read( (char *) sndpus, 9);
      int datalen = sndpus[4]*256 + sndpus[5] - 2;
      data = new uint8_t[datalen];
      tempFile.read( (char *) data, datalen);

      if (sndpus[7] == 253 && sndpus[8] == 3) {
        // Parameter Packet

        uint64_t ll;
        uint32_t ul;
      
        // Decode Exposure ID
        ul =  (uint32_t) data[15-9] << 24;
        ul += (uint32_t) data[16-9] << 16;
        ul += (uint32_t) data[17-9] << 8;
        ul += data[18-9];
        exposureID = ul;
      
        // Decode Epoch Time
        // Code taken from DecodeImageParams() in HawkeyeDecode.c
        ll =  (uint64_t) data[26-9] << 40;
        ll += (uint64_t) data[27-9] << 32;
        ll += (uint64_t) data[28-9] << 24;
        ll += (uint64_t) data[29-9] << 16;
        ll += (uint64_t) data[30-9] << 8;
        ll += data[31-9];
        ll = (ll & 0x0FFFFFFFFFF0) >> 4;

        //Decode Delta-Hawkeye time
        ul =  (uint32_t) data[36-9] << 24;
        ul += (uint32_t) data[37-9] << 16;
        ul += (uint32_t) data[38-9] << 8;
        ul += data[39-9];
        ul = (ul >> 4) & 0xFFFFFF;

        double tlm_time;

        tlm_time = ll * 0.001;
        tlm_time -= (53*24*3600 + 14788.224);
        getZuluTime_L0Name( tlm_time, &zuluStartTime, &packetOutname);
        startHWKTime = tlm_time;

        tlm_time = (ll + ul) * 0.001;
        tlm_time -= (53*24*3600 + 14788.224);
        getZuluTime_L0Name( tlm_time, &zuluStopTime, NULL);
        stopHWKTime = tlm_time;

        break;
        // Note HWK time is seconds from 2015/01/01
      }
    }
    tempFile.close();

    // Generate prePend
    prePend.assign(basename(argv[1]));
    size_t strpos = prePend.find(".dec");
    prePend = prePend.substr(0, strpos);
    prePend.append("_");
    
    std::ostringstream ss;
    exposureCount++;
    ss << setw(4) << setfill('0') << exposureCount;
    prePend.append(ss.str().c_str());
    prePend.append("_");
    
    prePend.append(to_string(exposureID));
    prePend.append("_");

    string outName;
    outName.assign(prePend);
    outName.append(packetOutname);

    rename( tempName.c_str(), outName.c_str());

    outName.append(".txt");

    textOut.open(outName.c_str(), ios::out | ios::trunc);
    textOut << "TlmFile=" << basename( argv[1]) << endl;
    textOut << "StartTime=" << zuluStartTime.c_str() << endl;
    textOut << "StopTime =" << zuluStopTime.c_str()  << endl;
    textOut << "ExposureID=" << to_string(exposureID).c_str()  << endl;
    textOut << "StartHWKTime=" << to_string(startHWKTime).c_str()  << endl;
    textOut << "StopHWKTime =" << to_string(stopHWKTime).c_str()  << endl;
    textOut.close();
  } // for (int i=0; i<tempCnt; i++)
    
  framefile.close();
  gpsOut.close();

  return 0;
}


int getZuluTime_L0Name( double tlm_time,
                        string *zuluTime, string *packetOutName) {

  int32_t iyr, idy;
  int16_t mon, idm, hr, min;
  double sec;

  tepoch2yds( tlm_time, &iyr, &idy, &sec);
  yd2md((int16_t) iyr, (int16_t) idy, &mon, &idm);

  //  cout << iyr << " " << idy << " " << sec << endl;

  stringstream ss;

  // yyyy-mn-dyThr:mn:ss.sss
  ss = stringstream();
  ss << setw(4) << to_string(iyr) << "-";
  ss << setw(2) << setfill('0') << mon << "-";
  ss << setw(2) << setfill('0') << idm << "T";
  hr =  (int) sec/3600;
  min = (int) ((sec - hr*3600) / 60);
  sec = sec - hr*3600 - min*60; 
  ss << setw(2) << setfill('0') << hr << ":";
  ss << setw(2) << setfill('0') << min << ":";
  ss << fixed << setw(6) << setprecision(3) << setfill('0') << sec;

  zuluTime->assign( ss.str());

  //  HWKyyyydoyhrmnsc.L0_SE1
  if (packetOutName != NULL) {
    ss = stringstream();
    ss << "HWK" << setw(4) << to_string(iyr);
    ss << setw(3) << setfill('0') << to_string(idy);
    ss << setw(2) << setfill('0') << hr;
    ss << setw(2) << setfill('0') << min;
    ss << setw(2) << setfill('0') << (int) sec;
    ss << ".L0_SE1";
    packetOutName->assign( ss.str());
  }

  return 0;
}
