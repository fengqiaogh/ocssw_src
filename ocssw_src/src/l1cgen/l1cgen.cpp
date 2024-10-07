
//**********************************************************
// l1cgen firstly defined as mainL1C.cpp
//  Originally created by Martin Montes on 3/27/2021 (first version)

//**********************************************************
#include <allocate2d.h>
#include <allocate3d.h>
#include "allocate4d.h"
#include <iostream>
#include "l1c_input.h"
#include "l1c_filehandle.h"
#include "l1c_str.h"
#include "l1c.h"
#include "hawkeye_methods.h"
#include <chrono>
#include <sys/stat.h>
#include "l2_str.h"

#include <filehandle.h>
#include <l1.h>
#include <genutils.h>

#include "l1c_latlongrid.h"
#include <netcdf>

#include <global_attrs.h>

using namespace std;
using namespace l1c;
using namespace std::chrono;
using namespace netCDF;
using namespace netCDF::exceptions;

#define PROGRAM_NAME "l1cgen"
#define VERSION "5.63 7/2/2024"


int main(int argc, char** argv) {
    int status;
    std::string str;
    double etime;
    string ifile_str;
    char* ifile_char;
    auto start = high_resolution_clock::now();   
    vector<pair<float, int>> vp;

    L1C_input l1cinput;
    l1c_filehandle l1cfile;
    l1c_str l1cstr;
    l2_str l2str;
    L1C* ptl1c = new L1C();
    filehandle l1file;
    l1str l1rec;
    //-------------------------------------------------------------------
    if (argc == 1) {
        l1cinput.l1c_usage(PROGRAM_NAME, VERSION);
        return 1;
    }

    for (int i = 0; i < argc; i++) {
        if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-help") == 0)) {
            l1cinput.l1c_usage(PROGRAM_NAME, VERSION);
        }
    }

    cout << PROGRAM_NAME << " " << VERSION << endl;
    l1cfile.version = VERSION;
    l1cfile.progname=PROGRAM_NAME;
    // grab CLI info
    l1cinput.l1c_inputmain(argc, argv, &l1cinput, &l1cfile, PROGRAM_NAME, VERSION);

    ifile_str = l1cinput.files[0];
    ifile_char = (char*)ifile_str.c_str();
    file_format format = getFormat(ifile_char);
    l1cfile.format = format.type;

    if (format.type != FT_INVALID) 
    {
        if (ptl1c->load_l1c_filehandle4(&l1cfile, &l1cinput) != 0) {
            printf("-E- %s: Error loading %sl1c filehandle.\n", argv[0], l1cfile.l1b_name.c_str());
            exit(1);
        }
           //----- READING HKT FILES -----------------------------------------------------------------
        //**********************************************************************************************************-
        if (l1cfile.l1c_pflag == 5) {
            cout << "Reading telemetry from L1A files..--> SOCEA -- L1C grid.." << endl;

            string history = "", stri;
            for (int i = 0; i < argc; i++) {
                stri = argv[i];
                history += stri + " ";
            }

            strcpy(l1cinput.history, history.c_str());

            if (status = ptl1c->open_l1atol1c3(&l1cinput, &l1cfile) > 1) {
                printf("-E- %d: status=2  error opening..L1A file for sensor...............\n", format.type);
                exit(1);
            }

            if (ptl1c != nullptr)
                delete ptl1c;

            auto stop1 = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop1 - start);
            etime = duration.count();
            etime /= (1000000 * 60);
            cout << "done processing options #5..(creating SOCEA L1C common grid from telemetry " << etime
                 << "minutes..." << endl;
        }

         else if (l1cfile.l1c_pflag == 7) {
            string senstr;
            string l1c_str = l1cinput.l1c_grid;
            cout << "producing L1C grid from L1C granules  and at CTH "
                    "......................................................"
                 << endl;
         
            if (l1cinput.sensor == 34) {
                senstr = "SPEXone";
                l1cfile.nbinx = 25;
            } else if (l1cinput.sensor == 30) {
                senstr = "OCI";
                l1cfile.nbinx = 519;
            } else if (l1cinput.sensor == 35) {
                senstr = "HARP2";
                l1cfile.nbinx = 457;
            } else if (l1cinput.sensor == 28) {
                senstr = "MISR";
                l1cfile.nbinx = 81;
            } else {
                if (l1cinput.verbose)
                    cout << "sensor by default is OCI option 2....." << endl;
                senstr = "OCI";
                l1cfile.nbinx = 519;
            }

            if (l1cinput.verbose) {
                cout << "number of L1B files to be processed..." << 1 << endl;
                cout << "number of L1C granules for screening scantime..." << l1cinput.files_l1c.size()
                     << endl;
            }

            ptl1c->l1_cloud_correct(&l1cinput, &l1cfile);

            if (ptl1c != nullptr)
                delete ptl1c;

            auto stop1 = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop1 - start);
            etime = duration.count();
            etime /= (1000000 * 60);
            cout << "done processing options #7..(cloud height parallax correction for L1B and L1C) in "
                 << etime << "minutes..." << endl;
        }

        //**********************************************************************************************************************
        else if (l1cfile.l1c_pflag == 8) 
        {
            //---------------------------------------------------------
            // DISCRETE BINNING -------------------------
            //---------------------------------------------------------
            if (l1cinput.verbose)
                cout << "init filehandle and add options.." << endl;
            filehandle_init(&l1file);
            l1_input_init();
            clo_optionList_t* optionList = clo_createList();
            l1_add_options(optionList);
            // add extra options from CLI and coming in msl12_defaults.par
            clo_addOption(optionList, "suite", CLO_TYPE_STRING, "OC", "suite of OC products");
            clo_addOption(optionList, "aermodfile", CLO_TYPE_STRING, "Unspecified",
                          "path dir for aerosol files anc data");
            clo_addOption(optionList, "aer_wave_short", CLO_TYPE_INT, "765",
                          "default shortest wavelength used for aerosol correction, epsilon");
            clo_addOption(optionList, "aer_wave_long", CLO_TYPE_INT, "865",
                          "default longest wavelength used for aerosol correction, epsilon");
            clo_addOption(optionList, "mumm_alpha", CLO_TYPE_FLOAT, "1.72",
                          "mumm alpha for AC Rudick turbid waters");
            clo_addOption(optionList, "mumm_gamma", CLO_TYPE_FLOAT, "1.0",
                          "mumm gamma for AC Rudick turbid waters");
            clo_addOption(optionList, "mumm_epsilon", CLO_TYPE_FLOAT, "1.0",
                          "mumm epsilon for AC Rudick turbid waters");
            clo_addOption(optionList, "chloc2_wave", CLO_TYPE_INT, "[-1,-1]",
                          "sensor wavelengths for OC2 chlorophyll\n        algorithm");
            clo_addOption(optionList, "chloc2_coef", CLO_TYPE_FLOAT, "[0.0,0.0,0.0,0.0,0.0]",
                          "coefficients for OC2\n        chlorophyll algorithm");
            clo_addOption(optionList, "chloc3_wave", CLO_TYPE_INT, "[-1,-1]",
                          "sensor wavelengths for OC3 chlorophyll\n        algorithm");
            clo_addOption(optionList, "chloc3_coef", CLO_TYPE_FLOAT, "[0.0,0.0,0.0,0.0,0.0]",
                          "coefficients for OC3\n        chlorophyll algorithm");
            clo_addOption(optionList, "chloc4_wave", CLO_TYPE_INT, "[-1,-1]",
                          "sensor wavelengths for OC4 chlorophyll\n        algorithm");
            clo_addOption(optionList, "chloc4_coef", CLO_TYPE_FLOAT, "[0.0,0.0,0.0,0.0,0.0]",
                          "coefficients for OC4\n        chlorophyll algorithm");
            clo_addOption(optionList, "kd2_wave", CLO_TYPE_INT, "[-1,-1]",
                          "sensor wavelengths for polynomial Kd(490)\n        algorithm");
            clo_addOption(optionList, "kd2_coef", CLO_TYPE_FLOAT, "[0.0,0.0,0.0,0.0,0.0,0.0]",
                          "sensor wavelengths\n        for polynomial Kd(490) algorithm");
            clo_addOption(optionList, "coccolith", CLO_TYPE_FLOAT, "[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]",
                          "coefficients for coccolith  algorithm");
            clo_addOption(optionList, "qaa_wave", CLO_TYPE_INT, NULL, "sensor wavelengths for QAA algorithm");
            clo_addOption(optionList, "giop_wave", CLO_TYPE_FLOAT, "-1",
                          "optimization comma-separated list, default is all visible bands (400-700nm)");
            clo_addOption(optionList, "aer_wave_base", CLO_TYPE_INT, "865",
                          "base sensor wavelength for aerosol \n         extrapolation");
            clo_addOption(optionList, "aer_opt", CLO_TYPE_INT, "99", "aerosol mode option");
            clo_addOption(optionList, "gas_opt", CLO_TYPE_INT, "1", "gaseous transmittance bitmask selector");
            clo_addOption(optionList, "mbac_wave", CLO_TYPE_INT, NULL,
                          "bands used for mbac atmospheric correction");
            clo_addOption(optionList, "wavelength_3d", CLO_TYPE_STRING, NULL,
                          "wavelength_3d input, written in ascending order\n"
                          "        with format 'wavelength_3d=nnn,nnn,nnn' where nnn is a sensor wavelength\n"
                          "        or a range of wavelengths as follows 'nnn:nnn'");
            clo_addOption(optionList, "uncertaintyfile", CLO_TYPE_IFILE, NULL, "uncertainty LUT");
            clo_addOption(optionList, "aermodels", CLO_TYPE_STRING, NULL, "aerosol models");
            clo_addOption(optionList, "oxaband_opt", CLO_TYPE_INT, "0", "oxygen A-band correction");
            clo_addOption(optionList, "filter_opt", CLO_TYPE_BOOL, NULL, "filtering input data option");
            clo_addOption(optionList, "filter_file", CLO_TYPE_IFILE, "$OCDATAROOT/sensor/sensor_filter.dat", "\n        data file for input filtering");
            clo_addOption(optionList, "absaer_opt", CLO_TYPE_INT, "0", "absorbing aerosol flagging option\n"
                            "        0: disabled\n"
                            "        1: use rhow_412 aerosol index test\n"
                            "        2: GMAO ancillary aerosol test");
            l1_read_default_files(optionList, &l1file, ifile_char);
            l1_load_options(optionList, &l1file);
            clo_deleteList(optionList);

            string l1c_str = l1cinput.l1c_grid;
            string l1c_str2 = l1cinput.l1c_anc;

            const char* l1c_grid = l1c_str.c_str();
            const char* l1c_anc = l1c_str2.c_str();

            bin_str binl1c;
            string history = "", stri;
            for (int i = 0; i < argc; i++) {
                stri = argv[i];
                history += stri + " ";
            }

            binl1c.history = history;
            binl1c.version = l1cfile.version;
            strcpy(binl1c.l1c_anc, l1cinput.l1c_anc);
            binl1c.verbose = l1cinput.verbose;
            binl1c.cloudem_flag = l1cinput.cloud_correct;
            binl1c.cloud_type = l1cinput.cloud_type;
            binl1c.dem_flag = l1cinput.demcloud_flag;

            int nfiles = l1cfile.ifiles.size();
            int timebin_ok[nfiles];
            int32_t npix_tot = 0;
            int timeflag = -1;
            const char* char_file;

            // checking time limits -------------
            for (int f = 0; f < nfiles; f++) {
                string str_file = l1cinput.files[f];
                char_file = str_file.c_str();
                if (l1cinput.verbose)
                    cout << "screening time for  l1b.." << str_file << "granule # " << f + 1 << "of.."
                         << l1cinput.files.size() << endl;
                strcpy(l1file.name, char_file);

                timeflag = check_l1c_time(char_file, l1c_grid, &binl1c);

                if (l1cinput.verbose)
                    cout << "time flag.." << timeflag << endl;
                if (timeflag == 0) {  // inside the grid
                    timebin_ok[f] = 0;
                } else
                    timebin_ok[f] = 1;
                if (l1cinput.verbose)
                    cout << "timebinok.." << timebin_ok[f] << endl;
            }

            if (l1cinput.verbose)
                cout << "checking time limits of ANC L1C file" << endl;
            if (binl1c.l1c_anc[0] != '\0' && strcmp(binl1c.l1c_anc, char_file) == 0) {
                timeflag = check_l1c_time(l1c_anc, l1c_grid, &binl1c);
            }

            //-------- creating output L1C file ------------------------
            string str_file;
            if (nfiles >= 2) {
                str_file = l1cinput.files[1];  // index 1 assuming 3 files, this is the middle one
            } else if (nfiles == 1) {
                str_file = l1cinput.files[0];
            } else {
                cout << "ERROR- no L1B files provided!!" << endl;
                exit(1);
            }

            char_file = str_file.c_str();
            strcpy(l1file.name, char_file);
            openl1(&l1file);
            if (alloc_l1(&l1file, &l1rec) == 0) {
                cout << "Unable to allocate L1 record....." << endl;
                exit(1);
            }
            int32_t spix = MAX(l1_input->spixl - 1, 0);
            int32_t epix = MAX(l1_input->epixl - 1, l1file.npix - 1);
            int32_t dpix = MAX(l1_input->dpixl, 1);
            int32_t sscan = MAX(l1_input->sline - 1, 0);
            int32_t escan = MAX(l1_input->eline - 1, l1file.nscan - 1);
            int32_t dscan = MAX(l1_input->dline - 1, 1);

            l1file.spix = spix;
            l1file.epix = epix;

            int32_t npix = (epix - spix) / dpix + 1;
            int32_t nscan = (escan - sscan) / dscan + 1;

            l1file.npix = npix;
            l1file.nscan = nscan;

            short** gdindex = nullptr;

            if (l1cinput.verbose)
                cout << "number of 'working' bands..." << l1file.nbands << endl;

            ifile_char = l1file.name;
            file_format format = getFormat(ifile_char);

            NcFile* nc_output;
            string l1c_full_str, ofile, senstr;
            if (format.type == FT_OCIL1B)
                senstr = "OCI";

            // check if ofile is empty
            if (strcmp(l1cinput.ofile, "") == 0) {
                if (l1cinput.verbose)
                    cout << "ofile not provided - use DEFAULT based on L1C grid.." << l1c_str << endl;
                string substr = ".";
                size_t ixstr1 = l1c_str.find(substr);
                if (ixstr1 == string::npos) {
                    cout << "ERROR building L1C ofile from L1C grid" << endl;
                    exit(1);
                }

                char arr[] = ".";
                size_t ixstr2 = l1c_str.find(arr, ixstr1 + 1);
                if (ixstr2 == string::npos) {
                    cout << "ERROR building L1C ofile from L1C grid" << endl;
                    exit(1);
                }

                ofile = "PACE_" + senstr + l1c_str.substr(ixstr1, ixstr2) + ".nc";
                if (l1cinput.verbose)
                    cout << "ofile default..." << ofile << endl;
                l1c_full_str = ofile;
                strcpy(l1cinput.ofile, ofile.c_str());
            } else {
                l1c_full_str = string(l1cinput.ofile);
                if (l1cinput.verbose)
                    cout << "ofile.." << l1cinput.ofile << endl;
            }

            binl1c.full_l1cgrid = l1c_full_str;

            cout << "creating FULL L1C file---------------" << l1c_full_str << endl;

            // empty file
            try {
                nc_output = new NcFile(l1cinput.ofile, NcFile::replace);
            } catch (NcException& e) {
                cerr << e.what() << "l1cgen l1c_pflag= 8:: Failure to replace FULL L1C grid: " + l1c_full_str
                     << endl;
                exit(1);
            }

            history = get_history(nc_output);
            history.append(call_sequence(argc, argv));
            set_global_attrs(nc_output, history);

            nc_output->close();

            if (l1cinput.verbose)
                cout << "adding more metadata to the FULL L1C file---------------" << l1c_full_str << endl;

            try {
                nc_output = new NcFile(l1cinput.ofile, NcFile::write);

            } catch (NcException& e) {
                cerr << e.what() << "l1cgen l1c_pflag= 8:: Failure to write FULL L1C grid: " + l1c_full_str
                     << endl;
                exit(1);
            }

            strcpy(binl1c.outlist, l1cinput.outlist);
            if(l1cinput.pversion[0])
                binl1c.pversion = l1cinput.pversion;
            if(l1cinput.doi[0])
                binl1c.doi = l1cinput.doi;

            meta_l1c_full(&l1file, &binl1c, l1cinput.l1c_grid,
                          nc_output);  // sensor view bands, geolocation, calls global attributes

            readl1(&l1file, 1, &l1rec);
            closel1(&l1file);

            gdindex = allocate2d_short(npix, 2);
            int32_t skipped[nfiles], nscan_tot = 0;

            for (int f = 0; f < nfiles; f++) {
                skipped[f] = 0;
                int firstcall=0;
                if (timebin_ok[f] == 0) {
                    l1cfile.l1b_name = l1cfile.ifiles[f];
                    const char* char_file = l1cfile.l1b_name.c_str();

                    //----binning --------------------------
                    strcpy(l1file.name, char_file);
                    openl1(&l1file);
                    if (l1cinput.verbose)
                        cout << "binning filename :" << l1file.name << endl;

                    if (alloc_l1(&l1file, &l1rec) == 0) {
                        cout << "Unable to allocate L1 record....." << endl;
                        exit(1);
                    }
                    int32_t spix = MAX(l1_input->spixl - 1, 0);
                    int32_t epix = MAX(l1_input->epixl - 1, l1file.npix - 1);
                    int32_t dpix = MAX(l1_input->dpixl, 1);
                    int32_t sscan = MAX(l1_input->sline - 1, 0);
                    int32_t escan = MAX(l1_input->eline - 1, l1file.nscan - 1);
                    int32_t dscan = MAX(l1_input->dline - 1, 1);

                    nscan_tot += l1file.nscan;

                    l1file.spix = spix;
                    l1file.epix = epix;

                    int32_t npix = (epix - spix) / dpix + 1;
                    int32_t nscan = (escan - sscan) / dscan + 1;

                    l1file.npix = npix;
                    l1file.nscan = nscan;
                    npix_tot += npix * nscan;

                    double timetemp = 0.0, time_offset = 0.;
                    for (int sline = 0; sline < l1file.nscan; sline++) {
                        readl1(&l1file, sline, &l1rec);
                    
                        // set timetemp to the beginning of the day
                        // only need to calculate for the first line
                        if (sline == 0) {
                            int16_t syear, smon, sday;
                            double secs;
                            unix2ymds(l1rec.scantime, &syear, &smon, &sday, &secs);
                            time_offset = ymds2unix(syear, smon, sday, 0.0);
                        }
                        timetemp = l1rec.scantime - time_offset;

                        int result = search_l1c(&l1file, &l1rec, &binl1c, gdindex);
                        if (result == 110) {
                            if (l1cinput.verbose) {
                                cout << "WARNING---- some lat/lon out of boundaries or FILLVALUES  at line # "
                                     << sline + 1 << endl;
                                cout << "skipping line # :" << sline + 1 << endl;
                            }
                            skipped[f] += 1;
                            binl1c.outpix += npix;
                        } else {
                            if (l1cinput.cloud_correct > 0) {
                                result =
                                    parallax(&l1file, l1c_anc, l1c_grid, &l1rec, &binl1c, gdindex, nc_output,sline,firstcall);
                                if (result == 110) {
                                    if (l1cinput.verbose) {
                                        cout << "WARNING---- some lat/lon out of boundaries or FILLVALUES  "
                                                "at line # "
                                             << sline + 1 << endl;
                                        cout << "skipping line # :" << sline + 1 << endl;
                                    }
                                    skipped[f] += 1;
                                }
                            }
                           
                            bintime_l1c(&l1file, &l1rec, &binl1c, gdindex, timetemp,
                                        nc_output);                                 // TIME VARS
                            bin_l1c(&l1file, &l1rec, &binl1c, gdindex, nc_output);  // NO TIME VARS
                            rmse_l1c_alt(&l1file, &binl1c, &l1rec, gdindex);        // diff2 height
                        }
                        if (sline % 100 == 0)
                            cout << "scanline " << sline + 1
                                 << ", pixels binned " << binl1c.inpix
                                 << ", pixels outside the grid " << binl1c.outpix << endl;
                    }

                    closel1(&l1file);
                } else {
                    if (l1cinput.verbose)
                        cout << "L1B granule # " << f + 1
                             << "is SKIPPED due to time beyond L1C granule boundaries" << endl;
                }
            }  // end files

            meta_l1c_bin(&l1file, &binl1c, nc_output);  // time offset/ mean geometry and I
            meta_l1c_altvar(&binl1c, nc_output);        // RMS or stdev  height

            free(gdindex);
            nc_output->close();
            ptl1c->add_proc_group_l1c(&l1cinput,&l1cfile,l1cinput.ofile);
    

            int32_t sk_tot = 0;
            for (int f = 0; f < nfiles; f++) {
                cout << "file # :" << f + 1 << " skipped lines # :" << skipped[f] << endl;
                if (nscan_tot > 0) {
                    cout << "percentage of missing scans due to geolocation # :"
                         << 100.0 * skipped[f] / nscan_tot << endl;
                } else {
                    cout << "WARNING: L1B FILES HAVE NO DATA!! " << endl;
                    break;
                }
                sk_tot += skipped[f];
            }

            cout << "#badgeo " << binl1c.badgeo << "#inpix: " << binl1c.inpix << "#outpix: " << binl1c.outpix
                 << "tot pixels: " << npix_tot << "missing scans %.." << 100.0 * sk_tot / nscan_tot << endl;

            binl1c.close_bin(&binl1c);
            delete (l1_input);

            if (ptl1c != nullptr)
                delete ptl1c;

            auto stop1 = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop1 - start);
            etime = duration.count();
            etime /= (1000000 * 60);

            cout << "done processing options #8..(using L1C functions from libl1 and l1 readers) in " << etime
                 << "minutes..." << endl;

            if (binl1c.inpix == 0 && binl1c.outpix == 0 && binl1c.badgeo > 0) {
                cout << "WARNING: L1B files have EMPTY OR /BAD geolocation data" << endl;
                return (110);
            } else if (binl1c.inpix == 0 && binl1c.outpix != 0 && binl1c.badgeo == 0) {
                cout << "WARNING: L1B files have geolocation data but outside the grid" << endl;
                return (130);
            } else
                return (0);
        }  
    }     

    else {
        printf("-E- %s: Error opening %s for reading unknown sensor...................\n", argv[0],
               l1cfile.l1b_name.c_str());
        return (1);
    }

    return (0);
}
