#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 15:56:47 EST 2019

@author: dshea
"""

import sys
import argparse
import subprocess
import multiprocessing as mp

machineInfo = [
        # RHEL 8
        {
           "name" : "linux_64", 
           "login" : "-p 10022 gs616-container103",
           "initStr" : "source .bash_profile",
           "exitStr" : ""
        },
        # Ubuntu 20.04
        {
            "name" : "odps", 
            "login" : "analysis701",
            "initStr" : "source .profile",
            "exitStr" : ""
        },
        # macOS intel
        {
            "name" : "macosx_intel", 
            "login" : "seadas4",
            "initStr" : "source .bash_profile",
            "exitStr" : " && fix_mac_rpath.py"
        },
        # macOS arm64
        {
            "name" : "macosx_arm64", 
            "login" : "seadas5",
            "initStr" : "source .zshrc",
            "exitStr" : ""
        },
        # Poseidon HPC, RHEL 8
        {
           "name" : "linux_hpc", 
           "login" : "poseidon",
           "initStr" : "source .bash_profile",
           "exitStr" : ""
        },
]

firstRunStr = " && cd \$OCSSWROOT/.. && if [ ! -d ocssw ]; then mkdir ocssw; fi && if [ ! -d ocssw/share ]; then mkdir ocssw/share; fi && if [ ! -d ocssw/testdata ]; then mkdir ocssw/testdata; fi && if [ ! -d ocssw/var ]; then mkdir ocssw/var; fi && if [ ! -d ocssw/opt ]; then mkdir ocssw/opt; fi   && if [ ! -d polarimetry ]; then mkdir polarimetry; fi"
saveOcsswStr = " && rm -rf saveOcssw && mkdir saveOcssw && cd ocssw && mv share testdata var ../saveOcssw"
restoreOcsswStr = " && rm -rf share/modis && mv ../saveOcssw/share ../saveOcssw/testdata ../saveOcssw/var ."
saveOptStr = " && mv opt ../saveOcssw"
restoreOptStr = " && mv ../saveOcssw/opt ."
getOcsswStr = " && cd \$OCSSWROOT/.. && rm -rf ocssw && git clone https://oceandata.sci.gsfc.nasa.gov/rcs/obpg/ocssw.git && cd ocssw && source OCSSW_bash.env"
getSubmodulesStr = " && git submodule init && git submodule update"
buildOptStr = " && ./get_lib3_src.sh && cd opt/src && ./BuildIt.py && cd ../.."
buildOcsswStr = " && mkdir build && cd build && cmake .. -DBUILD_ALL=1 && make -j 20 install"
getViirsStr = " && cd \$OCSSWROOT/.. && rm -rf viirs_l1 && git clone https://oceandata.sci.gsfc.nasa.gov/rcs/viirs/viirs_l1.git && cd viirs_l1"
buildViirsStr = " && mkdir build && cd build && cmake .. && make -j 20 install"
getFocsStr = " && cd \$OCSSWROOT/.. && rm -rf focs && git clone https://oceandata.sci.gsfc.nasa.gov/rcs/obpg/focs.git && cd focs"
buildFocsStr = " && mkdir build && cd build && cmake .. && make -j 20 install"
getDtdbStr = " && cd \$OCSSWROOT/.. && rm -rf dtdb && git clone https://oceandata.sci.gsfc.nasa.gov/rcs/aerosol/dtdb.git && cd dtdb"
buildDtdbStr = " && mkdir build && cd build && cmake .. && make -j 20 install"

getHarp2Str = " && cd \$OCSSWROOT/../polarimetry && rm -rf harp && git clone https://oceandata.sci.gsfc.nasa.gov/rcs/polarimetry/harp.git && cd harp"
buildHarp2HippStr = " && cd hipp && mkdir build && cd build && cmake .. && make install && cd .."

getL1bcgen_spexoneStr = " && cd \$OCSSWROOT/../polarimetry && rm -rf spex && git clone https://oceandata.sci.gsfc.nasa.gov/rcs/polarimetry/spex.git && cd spex"
buildL1bcgen_spexoneStr = " && mkdir build && cd build && cmake .. -DCMAKE_EXE_LINKER_FLAGS=\\\"-Wl,-rpath,\\\\\$ORIGIN/../opt/lib -L\$LIB3_DIR/lib -lnetcdf_c++4 -lnetcdf -lhdf5 -llapack -lblas\\\" -DCMAKE_PREFIX_PATH=\$LIB3_DIR -DL1BC_ONLY=1 -DCMAKE_BUILD_TYPE=Release && make -j 20 && cp spexone \$OCSSWROOT/bin/l1bcgen_spexone"

getFastmapolStr = " && cd \$OCSSWROOT/../polarimetry && rm -rf fastmapol && git clone https://oceandata.sci.gsfc.nasa.gov/rcs/polarimetry/fastmapol.git && cd fastmapol"
buildFastmapolStr = " && mkdir build && cd build && cmake .. && make install"

getUaaStr = " && cd \$OCSSWROOT/../sat && rm -rf unified_dtdb_aerosol && git clone https://oceandata.sci.gsfc.nasa.gov/rcs/sat/unified_dtdb_aerosol.git && cd unified_dtdb_aerosol"
buildUaaStr = " && mkdir build && cd build && cmake .. && make install"


def doIt(cmd, logFilename):
    logFile = open(logFilename, "w")
    result = subprocess.run(cmd, shell=True, stderr=subprocess.STDOUT, stdout=logFile)
    logFile.close()
    sys.exit(result.returncode)


def run():
    processes = []

    parser = argparse.ArgumentParser(description="Build OCSSW on all of our architectures")
    parser.add_argument("-t", "--tag", default=None,
                        help="git tag or branch that you want to build")
    parser.add_argument("-a", "--arch", default=None,
                        help="comma separated list of architectures that you want to build (linux_64,linux_hpc,odps,macosx_intel,macosx_arm64)")
    parser.add_argument("--build_opt", default=False, action="store_true", 
                        help="build opt (lib3) first")


    args = parser.parse_args()

    # make sure the arch list is valid
    if args.arch:
        machineList = []
        for info in machineInfo:
            machineList.append(info["name"])
        archList = args.arch.split(',')
        for arch in archList:
            if arch not in machineList:
                print("Architecture", arch, "is not in the list of supported architectures")
                sys.exit(1)

    runList = []

    archFound = False
    for info in machineInfo:
        if args.arch:
            if info["name"] not in args.arch:
                continue
        archFound = True
        
        runList.append(info["name"])
        print("\n---making", info["name"])
        
        commandLine = "ssh " + info["login"] + ' "set -x; ' + info["initStr"] + firstRunStr + saveOcsswStr
        if not args.build_opt:
            commandLine += saveOptStr
        commandLine += getOcsswStr
        if args.tag:
            commandLine += " && git checkout " + args.tag
        commandLine += getSubmodulesStr
        if args.build_opt:
            commandLine += buildOptStr
        else:
            commandLine += restoreOptStr
        commandLine += restoreOcsswStr
        commandLine += " && rm -rf ../saveOcssw"

        # build VIIRS
        commandLine += getViirsStr
        if args.tag:
            commandLine += " && git checkout " + args.tag
        commandLine += getSubmodulesStr
        commandLine += buildViirsStr

        # build FOCS
        # temporary exclusion for the tensor flow on macs
        if (info["name"] != "macosx_intel") and (info["name"] != "macosx_arm64"):
            commandLine += getFocsStr
            if args.tag:
                commandLine += " && git checkout " + args.tag
            commandLine += getSubmodulesStr
            commandLine += buildFocsStr

        # build dtdb
        commandLine += getDtdbStr
        if args.tag:
            commandLine += " && git checkout " + args.tag
        commandLine += getSubmodulesStr
        commandLine += buildDtdbStr

        # build OCSSW
        commandLine += " && cd \$OCSSWROOT"
        commandLine += buildOcsswStr

        # for now only build for ODPS
        if info["name"] == "odps":
            # get harp2 and build hipp
            commandLine += getHarp2Str
            if args.tag:
                commandLine += " && git checkout " + args.tag
            commandLine += buildHarp2HippStr

            # get fasstmapol and build
            commandLine += getFastmapolStr
            if args.tag:
                commandLine += " && git checkout " + args.tag
            commandLine += buildFastmapolStr

            # get UAA and build
            commandLine += getUaaStr
            if args.tag:
                commandLine += " && git checkout " + args.tag
            commandLine += buildUaaStr

        #
        # for now only build on poseidon 
        #
        if info["name"] == "linux_hpc":
            # build l1bcgen_spexone
            commandLine += getL1bcgen_spexoneStr
            if args.tag:
                commandLine += " && git checkout " + args.tag
            commandLine += buildL1bcgen_spexoneStr


        commandLine += info["exitStr"]
        commandLine += '"'

        logFilename = "build." + info["name"] + ".log"

        processes.append(mp.Process(target=doIt, name=info["name"], args=(commandLine, logFilename)))

    if not archFound:
        print("Architecture", args.arch, "is not valid")
        sys.exit(1)

    print("Starting Processes")
    
    # Run processes
    for p in processes:
        p.start()

    print("Waiting for Processes")

    # Exit the completed processes
    for p in processes:
        p.join()

    print("Checking results")
   
    # check the return result
    result = True
    for p in processes:
        if p.exitcode != 0:
            result = False
            print(p.name, "return code =", p.exitcode)
        else:
            print(p.name, "Success")
   
    if result:
        print()
        print("Everyting completed successfully.")
        return 0
    
    return 1


if __name__ == "__main__":
    sys.exit(run())
