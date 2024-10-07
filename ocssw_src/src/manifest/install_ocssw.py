#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 11:08:27 2018

@author: dshea
"""

from operator import indexOf
import sys
import re
import argparse
import os
import shutil
import subprocess
import urllib.parse
import hashlib
import json
import datetime
import tempfile
import tarfile
import smtplib 
from email.message import EmailMessage

#sys.path.insert(0, os.path.dirname(__file__))
import manifest as mf

# global variables
MANIFEST_BASENAME = "manifest.json"
BUNDLELIST_BASENAME = "bundleList.json"
INSTALL_OCSSW_BASENAME = "install_ocssw.json"
versionString = "8.0"
baseUrl = "https://oceandata.sci.gsfc.nasa.gov/manifest/tags"
manifestCommand = os.path.dirname(__file__) + "/manifest.py"
currentThing = 1
totalNumThings = 0

##########################################################################
#  WARNING - The initialBundleList should match the latest bundleList.json
##########################################################################
# mapping of bundle names to directories and extra options
bundleList = []
initialBundleList = [
    {"name":"root", "dir":".", "help":"random files in the root dir", "extra":"--exclude . --include bundleList.json --include OCSSW_bash.env --include OCSSW.env", "commandLine":True},
    {"name":"bin_linux_64", "dir":"bin_linux_64", "help":"executables for Linux", "commandLine":False}, 
    {"name":"bin_linux_hpc", "dir":"bin_linux_hpc", "help":"executables for Linux HPC", "commandLine":False}, 
    {"name":"bin_macosx_intel", "dir":"bin_macosx_intel", "help":"executables for Mac", "commandLine":False},
    {"name":"bin_macosx_arm64", "dir":"bin_macosx_arm64", "help":"executables for Mac ARM64", "commandLine":False},
    {"name":"bin_odps", "dir":"bin_odps", "help":"executables for ODPS", "commandLine":False},
    {"name":"lib_linux_64", "dir":"lib_linux_64", "help":"shared libraries for Linux", "commandLine":False}, 
    {"name":"lib_linux_hpc", "dir":"lib_linux_hpc", "help":"shared libraries for Linux HPC", "commandLine":False}, 
    {"name":"lib_macosx_intel", "dir":"lib_macosx_intel", "help":"shared libraries for Mac", "commandLine":False},
    {"name":"lib_macosx_arm64", "dir":"lib_macosx_arm64", "help":"shared libraries for Mac ARM64", "commandLine":False},
    {"name":"lib_odps", "dir":"lib_odps", "help":"shared libraries for ODPS", "commandLine":False},
    {"name":"opt_linux_64", "dir":"opt_linux_64", "help":"3rd party library for Linux", "extra": "--exclude src", "commandLine":False},
    {"name":"opt_linux_hpc", "dir":"opt_linux_hpc", "help":"3rd party library for Linux HPC", "extra": "--exclude src", "commandLine":False},
    {"name":"opt_macosx_intel", "dir":"opt_macosx_intel", "help":"3rd party library for Mac", "extra": "--exclude src", "commandLine":False},
    {"name":"opt_macosx_arm64", "dir":"opt_macosx_arm64", "help":"3rd party library for Mac ARM64", "extra": "--exclude src", "commandLine":False},
    {"name":"opt_odps", "dir":"opt_odps", "help":"3rd party library for ODPS", "extra": "--exclude src", "commandLine":False},

    {"name":"ocssw_src", "dir":"ocssw_src", "help":"OCSSW source code", "commandLine":False},
    {"name":"opt_src", "dir":"opt/src", "help":"3rd party library sources", "extra":"--exclude .buildit.db", "commandLine":True},
    
    {"name":"afrt", "dir":"share/afrt", "help":"Ahmad-Fraser RT data", "commandLine":True},
    {"name":"aquaverse", "dir":"share/aquaverse", "help":"Algorithm based on Mixture Density Networks", "commandLine":True},
    {"name":"avhrr", "dir":"share/avhrr", "help":"AVHRR", "commandLine":True},
    {"name":"aviris", "dir":"share/aviris", "help":"AVIRIS", "commandLine":True},
    {"name":"common", "dir":"share/common", "help":"common", "commandLine":True},
    {"name":"czcs", "dir":"share/czcs", "help":"CZCS", "commandLine":True},
    {"name":"eval", "dir":"share/eval", "help":"evaluation", "commandLine":True},
    {"name":"goci", "dir":"share/goci", "help":"GOCI", "commandLine":True},
    {"name":"harp2", "dir":"share/harp2", "help":"HARP 2 Polarimeter", "commandLine":True},
    {"name":"hawkeye", "dir":"share/hawkeye", "help":"Hawkeye", "commandLine":True},
    {"name":"hico", "dir":"share/hico", "help":"HICO", "commandLine":True},
    {"name":"l5tm", "dir":"share/l5tm", "help":"l5tm", "commandLine":True},
    {"name":"l7etmp", "dir":"share/l7etmp", "help":"l7etmp", "commandLine":True},
    {"name":"meris", "dir":"share/meris", "help":"MERIS", "commandLine":True},
    {"name":"misr", "dir":"share/misr", "help":"MISR", "commandLine":True},
    {"name":"modis", "dir":"share/modis", "help":"MODIS common", "extra":"--exclude aqua --exclude terra", "commandLine":False},
    {"name":"modisa", "dir":"share/modis/aqua", "help":"MODIS AQUA", "commandLine":True},
    {"name":"modist", "dir":"share/modis/terra", "help":"MODIS TERRA", "commandLine":True},
    {"name":"msi", "dir":"share/msi", "help":"MSI Sentinel 2 common", "extra":"--exclude s2a --exclude s2b", "commandLine":False},
    {"name":"msis2a", "dir":"share/msi/s2a", "help":"MSI Sentinel 2A", "commandLine":True},
    {"name":"msis2b", "dir":"share/msi/s2b", "help":"MSI Sentinel 2B", "commandLine":True},
    {"name":"oci", "dir":"share/oci", "help":"PACE OCI", "commandLine":True},
    {"name":"ocm1", "dir":"share/ocm1", "help":"OCM1", "commandLine":True},
    {"name":"ocm2", "dir":"share/ocm2", "help":"OCM2", "commandLine":True},
    {"name":"ocrvc", "dir":"share/ocrvc", "help":"OC Virtual Constellation", "commandLine":True},
    {"name":"octs", "dir":"share/octs", "help":"OCTS", "commandLine":True},
    {"name":"olci", "dir":"share/olci", "help":"OLCI Sentinel 3 common", "extra":"--exclude s3a --exclude s3b", "commandLine":False},
    {"name":"olcis3a", "dir":"share/olci/s3a", "help":"OLCI Sentinel 3A", "commandLine":True},
    {"name":"olcis3b", "dir":"share/olci/s3b", "help":"OLCI Sentinel 3B", "commandLine":True},
    {"name":"oli", "dir":"share/oli", "help":"OLI Landsat", "extra":"--exclude l8 --exclude l9", "commandLine":False},
    {"name":"olil8", "dir":"share/oli/l8", "help":"OLI Landsat 8", "commandLine":True},
    {"name":"olil9", "dir":"share/oli/l9", "help":"OLI Landsat 9", "commandLine":True},
    {"name":"prism", "dir":"share/prism", "help":"PRISM", "commandLine":True},
    {"name":"sabiamar", "dir":"share/sabiamar", "help":"Sabiamar", "commandLine":True},
    {"name":"seawifs", "dir":"share/seawifs", "help":"SeaWiFS", "commandLine":True},
    {"name":"sgli", "dir":"share/sgli", "help":"SGLI", "commandLine":True},
    {"name":"spexone", "dir":"share/spexone", "help":"SPEX One Polarimeter", "commandLine":True},
    {"name":"viirs", "dir":"share/viirs", "extra":"--exclude dem --exclude j1 --exclude j2 --exclude npp", "help":"VIIRS common", "commandLine":False},
    {"name":"viirsdem", "dir":"share/viirs/dem", "help":"VIIRS Digital Elevation", "commandLine":True},
    {"name":"viirsj1", "dir":"share/viirs/j1", "help":"VIIRS JPSS1", "commandLine":True},
    {"name":"viirsj2", "dir":"share/viirs/j2", "help":"VIIRS JPSS2", "commandLine":True},
    {"name":"viirsn", "dir":"share/viirs/npp", "help":"VIIRS NPP", "commandLine":True},
    {"name":"wv3", "dir":"share/wv3", "help":"WV3", "commandLine":True},
    {"name":"aerosol", "dir":"share/aerosol", "help":"aerosol processing with dtdb", "commandLine":True},
    {"name":"cloud", "dir":"share/cloud", "help":"cloud properties processing", "commandLine":True},
    {"name":"telemetery", "dir":"share/telemetry", "help":"telemetry packet descriptions", "commandLine":True},
    {"name":"benchmark", "dir":"benchmark", "help":"benchmark MOSIS Aqua, level0 -> level3 Mapped", "commandLine":True},
    {"name":"viirs_l1_benchmark", "dir":"viirs_l1_benchmark", "help":"VIIRS benchmark data", "commandLine":True},
    {"name":"viirs_l1_bin_macosx_intel", "dir":"viirs_l1_bin_macosx_intel", "help":"Subset of binary files for VIIRS", "commandLine":False},
    {"name":"viirs_l1_bin_macosx_arm64", "dir":"viirs_l1_bin_macosx_arm64", "help":"Subset of binary files for VIIRS", "commandLine":False},
    {"name":"viirs_l1_bin_linux_64", "dir":"viirs_l1_bin_linux_64", "help":"Subset of binary files for VIIRS", "commandLine":False},
    {"name":"viirs_l1_bin_odps", "dir":"viirs_l1_bin_odps", "help":"Subset of binary files for VIIRS", "commandLine":False}
]

# list of bundles that have luts
lutBundles = ["seawifs", "modisa", "modist", "viirsn", "viirsj1", "viirsj2", "oci"]


def findBundleInfo(bundleName):
    for bundleInfo in bundleList:
        if bundleInfo["name"] == bundleName:
            return bundleInfo
    return None

def getArch():
    """
    Return the system arch string.
    """
    (sysname, _, _, _, machine) = os.uname()
    if sysname == 'Darwin':
        if machine == 'x86_64' or machine == 'i386':
            return 'macosx_intel'
        if machine == 'arm64':
            return 'macosx_arm64'
        print("unsupported Mac machine =", machine)
        exit(1)
    if sysname == 'Linux':
        if machine == 'x86_64':
            return 'linux_64'
        print("Error: can only install OCSSW software on 64bit Linux")
        exit(1)
    if sysname == 'Windows':
        print("Error: can not install OCSSW software on Windows")
        exit(1)
    print('***** unrecognized system =', sysname, ', machine =', machine)
    print('***** defaulting to linux_64')
    return 'linux_64'

def runCommand(command, pipe=False, shellVal=False):
    if pipe:
        proc = subprocess.run(command, stdout=subprocess.PIPE, shell=shellVal)
    else:
        proc = subprocess.run(command, shell=shellVal)
    if proc.returncode != 0:
        print("Error: trying to run = ", command)
        sys.exit(1)
    if pipe:
        return proc.stdout

def listTags(options):
    manifest_options = mf.create_default_options()
    manifest_options.wget = options.wget
    mf.list_tags(manifest_options, None)

def checkTag(options):
    manifest_options = mf.create_default_options()
    manifest_options.tag = options.tag
    manifest_options.wget = options.wget
    if not mf.check_tag(manifest_options, None):
        print("tag=%s," % (options.tag), "Does not exist.")
        sys.exit(1)

def recordInstalledTag(options):
    installInfo = {}
    installInfo["tag"] = options.tag
    installInfo["version"] = versionString
    installInfo["install_date"] = datetime.datetime.now().isoformat(timespec="seconds")

    # save as a JSON file
    if not os.path.isdir(options.install_dir):
        os.makedirs(options.install_dir)

    fileName = "%s/%s" % (options.install_dir, INSTALL_OCSSW_BASENAME)
    with open(fileName, 'w') as outfile:
        json.dump(installInfo, outfile, indent=4, sort_keys=True)

def getInstalledTag(options):
    installOcsswFilename = "%s/%s" % (options.install_dir, INSTALL_OCSSW_BASENAME)
    try:
        # try reading the install_ocssw.json file
        with open(installOcsswFilename, 'rb') as installOcsswFile:
            installInfo = json.load(installOcsswFile)
            return installInfo["tag"]
    except:
        # if that does not work, grab the last tag out of the manifest.json file
        try:
            manifestFilename = "%s/%s" % (options.install_dir, MANIFEST_BASENAME)
            with open(manifestFilename, 'rb') as manifestFile:
                manifest = json.load(manifestFile)
                return manifest["tags"][-1]
        except:
            # if nothing works, do this
            return "Unknown"

def installedTag(options):
    print("installedTag =", getInstalledTag(options))

def installBundle(options, bundleInfo):
    global currentThing
    if options.verbose:
        print()
    print("Installing (" + str(currentThing), "of", str(totalNumThings) + ") -", bundleInfo["name"], flush=True)
    currentThing += 1

    manifest_options = mf.create_default_options()
    manifest_options.verbose = options.verbose
    manifest_options.name = bundleInfo["name"]
    manifest_options.tag = options.tag
    manifest_options.dest_dir = "%s/%s" % (options.install_dir, bundleInfo["dir"])
    manifest_options.save_dir = options.save_dir
    manifest_options.local_dir = options.local_dir
    manifest_options.wget = options.wget

    mf.download(manifest_options, None)

    if options.clean:
        command = "%s clean %s/%s" % (manifestCommand, options.install_dir, bundleInfo["dir"])
        if "extra" in bundleInfo:
            command += " " + bundleInfo["extra"]
        runCommand(command, shellVal=True)

def updateLuts(options, lut):
    runner = options.install_dir + "/bin/ocssw_runner"
    if not os.path.isfile(runner):
        print("Error - bin directory needs to be installed.")
        exit(1)
    if options.verbose:
        print()
    print("Installing lut -", lut)
    command = "%s --ocsswroot %s update_luts %s" % (runner, options.install_dir, lut)
    runCommand(command, shellVal=True)

def getBundleListTag(manifestFilename):
    try:
        with open(manifestFilename, 'rb') as manifestFile:
            manifest = json.load(manifestFile)
            return manifest['files'][BUNDLELIST_BASENAME]['tag']
    except json.JSONDecodeError:
        print(manifestFilename, "is not a manifest file")
        sys.exit(1)
    except FileNotFoundError:
        print(manifestFilename, "Not found")
        sys.exit(1)
    except KeyError:
        print(manifestFilename, "is corrupt")
        sys.exit(1)
    print("could not find bundeList tag in", manifestFilename)
    sys.exit(1)

def downloadBundleList(options):
    global bundleList

    if not options.tag:
        print("\nWARNING: --tag is required to get the proper bundle list.\n")
        bundleList = initialBundleList
        return
    else:
        checkTag(options)

    tempDir = tempfile.TemporaryDirectory(prefix="install_ocssw-")
    tempBundleFilePath = dest = "%s/%s" % (tempDir.name, BUNDLELIST_BASENAME)
    if options.local_dir:
        manifestFilename = "%s/%s/root/%s" % (options.local_dir, options.tag, MANIFEST_BASENAME)
        bundleFilename = "%s/%s/root/%s" % (options.local_dir, getBundleListTag(manifestFilename),BUNDLELIST_BASENAME)
        if not os.path.isfile(bundleFilename):
            print(bundleFilename, "file does not exist")
            sys.exit(1)
        shutil.copy(bundleFilename, tempBundleFilePath)
    elif options.wget:
        command = "cd %s; wget -q %s/%s/root/%s" % (tempDir.name, options.base_url, options.tag, MANIFEST_BASENAME)
        runCommand(command, shellVal=True)
        manifestFilename = "%s/%s" % (tempDir.name, MANIFEST_BASENAME)
        bundleListUrl = "%s/%s/root/%s" % (options.base_url, getBundleListTag(manifestFilename), BUNDLELIST_BASENAME)
        command = "cd %s; wget -q %s" % (tempDir.name, bundleListUrl)
        runCommand(command, shellVal=True)
    else:
        manifestUrl = "%s/%s/root/%s" % (options.base_url, options.tag, MANIFEST_BASENAME)
        parts = urllib.parse.urlparse(manifestUrl)
        host = parts.netloc
        request = parts.path
        status = mf.httpdl(host, request, localpath=tempDir.name, outputfilename=MANIFEST_BASENAME, force_download=True)
        if status != 0:
            print("Error downloading", manifestUrl, ": return code =", status)
            sys.exit(1)
        manifestFilename = "%s/%s" % (tempDir.name, MANIFEST_BASENAME)
        bundleListUrl = "%s/%s/root/%s" % (options.base_url, getBundleListTag(manifestFilename), BUNDLELIST_BASENAME)
        parts = urllib.parse.urlparse(bundleListUrl)
        host = parts.netloc
        request = parts.path
        status = mf.httpdl(host, request, localpath=tempDir.name, outputfilename=BUNDLELIST_BASENAME, force_download=True)
        if status != 0:
            print("Error downloading", bundleListUrl, ": return code =", status)
            sys.exit(1)
    with open(tempBundleFilePath, 'rb') as bundleListFile:
        bundleList = json.load(bundleListFile)


def downloadFile(options, tag, bundleInfo, fileName, destDir):
    dest = "%s/%s" % (destDir, fileName)
    dest_dir = os.path.dirname(dest)
    if not os.path.isdir(dest_dir):
        os.makedirs(dest_dir)

    if options.local_dir:
        src = "%s/%s/%s/%s" % (options.local_dir, options.tag, bundleInfo["name"], fileName)
        if options.verbose:
            print("Copying %s from %s" % (fileName, src))
        shutil.copy(src, dest)
        return True
    
    url = "%s/%s/%s/%s" % (baseUrl, tag, bundleInfo["name"], fileName)
    if options.verbose:
        print("Downloading %s from %s" % (fileName, url))
    if options.wget:
        if os.path.isfile(dest):
            os.remove(dest)
        command = "cd %s; wget -q %s" % (dest_dir, url)
        runCommand(command, shellVal=True)
    else:
        parts = urllib.parse.urlparse(url)
        host = parts.netloc
        request = parts.path
        status = mf.httpdl(host, request, localpath=dest_dir, 
                        outputfilename=os.path.basename(dest),
                        verbose=options.verbose,
                        force_download=True)
        if status != 0:
            print("Error downloading", dest, ": return code =", status)
            return False
    return True

def bundleStatus(options, bundleInfo, fileLst_chg, fileLst_relPath_chg, fileLst_del):
    returnStatus = True
    statusTag = "tempStatusTag"
    #os.chdir(options.install_dir)

    currentManifest = {}
    tempDir = tempfile.TemporaryDirectory(prefix="install_ocssw-")
    currentManifestFilename = "%s/%s" % (tempDir.name, MANIFEST_BASENAME)
    if downloadFile(options, options.tag, bundleInfo, MANIFEST_BASENAME, tempDir.name):
        with open(currentManifestFilename, 'rb') as manifestFile:
            currentManifest = json.load(manifestFile)

    command = "cd %s; " % (options.install_dir) 
    command += manifestCommand + " generate"
    command += " -n " + bundleInfo["name"]
    command += " -t " + statusTag
    command += " -b " + currentManifestFilename
    command += " " + bundleInfo["dir"]
    if "extra" in bundleInfo:
        command += " " + bundleInfo["extra"]

    # command = [manifestCommand, "generate", 
    #            "-n", bundleInfo["name"], 
    #            "-t", statusTag,
    #            "-b", currentManifestFilename,
    #            bundleInfo["dir"]]
    # if "extra" in bundleInfo:
    #     command += bundleInfo["extra"].split()
    statusManifest = json.loads(runCommand(command, pipe=True, shellVal=True))

    if "tags" in currentManifest:
        tagList = currentManifest["tags"]
        versionStr = "("
        if len(tagList) > 1:
            versionStr += tagList[0] + ".."
        versionStr += tagList[-1] + ")"
        print(bundleInfo["name"], versionStr)
    else:
        print(bundleInfo["name"])

    
        
    # list new files
    fileList = []
    for f, info in statusManifest["files"].items():
        if ("files" not in currentManifest) or (f not in currentManifest["files"]):
            fileList.append(f)
            if options.create_change:
                fileLst_chg.append("%s/%s/%s" % (options.install_dir, bundleInfo["dir"], f))
                fileLst_relPath_chg.append("%s/%s" % (bundleInfo["dir"], f))
    if fileList:
        returnStatus = False
        print("  New Files:")
        for f in fileList:
            print("    %s/%s" % (bundleInfo["dir"], f))
        
    # list modified files
    fileList = []
    for f, info in statusManifest["files"].items():
        if ("files" in currentManifest) and (f in currentManifest["files"]):
            if "symlink" in currentManifest["files"][f]:
                 # both symlink
                if "symlink" in info:                                              
                    if info["symlink"] != currentManifest["files"][f]["symlink"]:
                        fileList.append(f)  
                        if options.create_change:
                            fileLst_chg.append("%s/%s/%s" % (options.install_dir, bundleInfo["dir"], f))
                            fileLst_relPath_chg.append("%s/%s" % (bundleInfo["dir"], f))
                # remote symlink 
                else:
                    fileList.append(f)
                    if options.create_change:
                        fileLst_chg.append("%s/%s/%s" % (options.install_dir, bundleInfo["dir"], f))
                        fileLst_relPath_chg.append("%s/%s" % (bundleInfo["dir"], f))
            # local symlink                  
            elif "symlink" in info: 
                fileList.append(f)
                if options.create_change:
                    fileLst_chg.append("%s/%s/%s" % (options.install_dir, bundleInfo["dir"], f))
                    fileLst_relPath_chg.append("%s/%s" % (bundleInfo["dir"], f))
            # neither symlink              
            elif info["tag"] != currentManifest["files"][f]["tag"]:
                localFilename = "%s/%s/%s" % (options.install_dir, bundleInfo["dir"], f)
                with open(localFilename, "rb") as localFile:
                    bytes = localFile.read()
                    localFileChecksum = hashlib.sha256(bytes).hexdigest() 
                    if "checksum" in currentManifest["files"][f]:    
                        if localFileChecksum != currentManifest["files"][f]["checksum"]:
                            fileList.append(f)   
                            if options.create_change:
                                fileLst_chg.append("%s/%s/%s" % (options.install_dir, bundleInfo["dir"], f))
                                fileLst_relPath_chg.append("%s/%s" % (bundleInfo["dir"], f))
    if fileList:                
        returnStatus = False
        print("  Modified Files:")
        for f in fileList:
            print("    %s/%s" % (bundleInfo["dir"], f))
            if options.diff:
                if "symlink" in currentManifest["files"][f]:
                    # both symlink
                    if  "symlink" in "symlink" in statusManifest["files"][f]:
                        print("diff %s remote..local" % (f))
                        print("<      " + currentManifest["files"][f]["tag"] + "     = " + currentManifest["files"][f]["symlink"])
                        print("---")
                        print(">      current dir = " + statusManifest["files"][f]["symlink"])
                    # remote symlink
                    else:
                        print("diff %s remote..local" % (f))
                        print("<      " + currentManifest["files"][f]["tag"] + "     = " + currentManifest["files"][f]["symlink"])
                        print("---")
                        print(">      current dir = Not a symlink")
                # local symlink
                elif "symlink" in "symlink" in statusManifest["files"][f]:
                        print("diff %s remote..local" % (f))
                        print("<      " + currentManifest["files"][f]["tag"] + "     = Not a symlink")
                        print("---")
                        print(">      current dir = " + statusManifest["files"][f]["symlink"])
                #neither symlink
                else:
                    try:
                        if re.search('root', bundleInfo["name"]):
                            localFilename = "%s/%s" % (options.install_dir, f)
                        else:
                            localFilename = "%s/%s/%s" % (options.install_dir, bundleInfo["dir"], f)
                        with open(localFilename) as file_1:
                            file_1_text = file_1.readlines()
                        remoteFilename = "%s/%s" % (tempDir.name, f)
                        file_tag = currentManifest["files"][f]["tag"]
                        if downloadFile(options, file_tag, bundleInfo, f, tempDir.name):
                            if options.difftool is not None:
                                command = "%s %s %s" % (options.difftool, remoteFilename, localFilename)
                                subprocess.run(command, shell=True)
                            else:
                                print("diff %s remote..local" % (f))
                                command = "diff %s %s" % (remoteFilename, localFilename)
                                subprocess.run(command, shell=True)
                            
                    except:
                        continue

    # deleted file     
    fileList = []
    if "files" in currentManifest:
        for f, info in currentManifest["files"].items():
            if f not in statusManifest["files"]:
                fileList.append(f)
                if options.create_change:
                    # fileLst_del.append("%s/%s/%s" % (options.install_dir, bundleInfo["dir"], f))
                    fileLst_del.append("%s/%s" % (bundleInfo["dir"], f))
    if fileList:
        returnStatus = False
        print("  Deleted Files:")
        for f in fileList:
            print("    %s/%s" % (bundleInfo["dir"], f))
    return returnStatus

def status(options):
    global bundleList

    # theInstalledTag = getInstalledTag(options)
    # print("installedTag =", theInstalledTag)

    # if not options.tag:
    #     options.tag = theInstalledTag

    theStatus = True
    bundle_found = False
    theInstalledTag = ''
    fileList_change = []
    #relative path of the fileList_change    
    fileList_relPath_change = []              
    fileList_delete = []
     

    #loop through provided bundles e.g. --hico
    for bundleInfo in bundleList:
        if hasattr(options, bundleInfo["name"]) and getattr(options, bundleInfo["name"]):
            # if not re.search("root", bundleInfo["name"]):
            if not theInstalledTag:
                theInstalledTag = getInstalledTag(options)
                print("installedTag =", theInstalledTag)

                if not options.tag:
                    options.tag = theInstalledTag
            if os.path.exists("%s/%s/%s" % (options.install_dir, bundleInfo["dir"], MANIFEST_BASENAME)):
                bundle_found = True
                if not bundleStatus(options, bundleInfo, fileList_change, fileList_relPath_change, fileList_delete):
                    theStatus = False

    if not bundle_found:
        startingDir = os.path.abspath(os.getcwd())
        for bundleInfo in bundleList:
            bundleDir = os.path.abspath(os.path.join(options.install_dir, bundleInfo["dir"]))
            if startingDir not in bundleDir:
                continue
            if not re.search("root", bundleInfo["name"]):
                if not theInstalledTag:
                    theInstalledTag = getInstalledTag(options)
                    print("installedTag =", theInstalledTag)

                    if not options.tag:
                        options.tag = theInstalledTag

                if os.path.exists("%s/%s/%s" % (options.install_dir, bundleInfo["dir"], MANIFEST_BASENAME)):
                    bundle_found = True
                    if not bundleStatus(options, bundleInfo, fileList_change, fileList_relPath_change, fileList_delete):
                        theStatus = False

    #loop through all the bundles
    if not bundle_found:
        theInstalledTag = getInstalledTag(options)
        print("installedTag =", theInstalledTag)

        if not options.tag:
            options.tag = theInstalledTag

        for bundleInfo in bundleList:    
                if os.path.exists("%s/%s/%s" % (options.install_dir, bundleInfo["dir"], MANIFEST_BASENAME)):
                    bundle_found = True
                    if not bundleStatus(options, bundleInfo, fileList_change, fileList_relPath_change, fileList_delete):
                        theStatus = False

    if options.create_change:
        uid = runCommand("whoami", True)
        date = datetime.date.today().strftime("%Y-%m-%d")
        if len(fileList_change) != 0:
            file_name = options.tag + "_" + date + "_" + str(uid, 'UTF-8').strip()
            tar_filename = file_name + ".tar"
            tar = tarfile.open(file_name + ".tar", "w")
            print("Creating " + tar_filename)
            for index, name in enumerate(fileList_change):
                tar.add(name, fileList_relPath_change[index])
            tar.close()

        if len(fileList_delete) != 0:
            remove_filename = file_name + ".remove"
            print("Creating " + remove_filename)
            with open(file_name + ".remove", "w") as fp:
                for name in fileList_delete:
                    # write each item on a new line
                    fp.write("%s\n" % name)
        
        if options.deliver_change or options.deliver_email_from is not None:
            glusteruser_dir = "/glusteruser/analyst/share_delivery"

            if(os.path.exists(glusteruser_dir)):
                if len(fileList_change) != 0:
                    print("Copying .tar to /glusteruser/analyst/share_delivery")
                    shutil.copyfile(tar_filename, os.path.join(glusteruser_dir, tar_filename))
                if len(fileList_delete) != 0:
                    print("Copying .remove files to /glusteruser/analyst/share_delivery")
                    shutil.copyfile(remove_filename, os.path.join(glusteruser_dir, remove_filename))

                if len(fileList_change) != 0 or len(fileList_delete) != 0:
                    message = EmailMessage()
                    if len(fileList_delete) == 0:
                        print("Sending email to swdevels@oceancolor.gsfc.nasa.gov notifying the delivery of the .tar files")
                        message.set_content(tar_filename)
                        message['Subject'] = '.tar file copied to /glusteruser/analyst/share_delivery'
                    elif len(fileList_change) == 0:
                        print("Sending email to swdevels@oceancolor.gsfc.nasa.gov notifying the delivery of the .remove files")
                        message.set_content(remove_filename)
                        message['Subject'] = '.remove file copied to /glusteruser/analyst/share_delivery'
                    else: 
                        print("Sending email to swdevels@oceancolor.gsfc.nasa.gov notifying the delivery of the .tar and .remove files")
                        message.set_content(tar_filename + "\n" + remove_filename)
                        message['Subject'] = '.tar and.remove file copied to /glusteruser/analyst/share_delivery'

                    if options.deliver_email_from is not None:
                        message['From'] = options.deliver_email_from
                    elif options.deliver_change:
                        install_ocssw_config_filename = os.path.join(os.path.expanduser('~'), ".install_ocssw")
                        with open(install_ocssw_config_filename, "r") as file_email:
                            message['From'] = file_email.readlines()[0].strip().split('deliver_email_from=')[1]
                    message['To'] = 'swdevels@oceancolor.gsfc.nasa.gov'

                    smtp_server = smtplib.SMTP('localhost')
                    smtp_server.send_message(message)
                    smtp_server.quit()
            else:
                print("You do not have access to /glusteruser/analyst/share_dir")
      

    return theStatus

def update(options):
    global bundleList
    global totalNumThings

    theInstalledTag = getInstalledTag(options)
    print("installedTag =", theInstalledTag)

    # count the number of bundles to install
    totalNumThings = 0
    for bundleInfo in bundleList:
        if os.path.exists("%s/%s/%s" % (options.install_dir, bundleInfo["dir"], MANIFEST_BASENAME)):
            totalNumThings += 1

    for bundleInfo in bundleList:
        if os.path.exists("%s/%s/%s" % (options.install_dir, bundleInfo["dir"], MANIFEST_BASENAME)):
            installBundle(options, bundleInfo)

    recordInstalledTag(options)

def run():
    global totalNumThings
    global bundleList

    # first make a parser to download the bundleInfo file
    parser = argparse.ArgumentParser(description="Install OCSSW bundles", add_help=False)
    parser.add_argument("-t", "--tag", default=None,
                        help="tag that you want to install")
    parser.add_argument("-b", "--base_url", default=baseUrl, 
                        help="remote url for the bundle server")
    parser.add_argument("-l", "--local_dir", default=None, 
                        help="local directory to use for bundle source instead of the bundle server")
    parser.add_argument("--wget", default=False, action="store_true", 
                        help="use wget for file download")
    parser.add_argument("--list_tags", default=False, action="store_true", 
                        help="list the tags that exist on the server")
    parser.add_argument("--version", default=False, action="store_true", 
                        help="print this program's version")

    options1 = parser.parse_known_args()
    if not options1[0].list_tags and not options1[0].version:
        downloadBundleList(options1[0])

    epilog = """
        The content of the config file ~/.install_ocssw --
            deliver_email_from=a valid email address that is on the list of swdevels@oceancolor.gsfc.nasa.gov 
    """
    # now make the real parser
    parser = argparse.ArgumentParser(description="Install OCSSW bundles", epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter,)
    parser.add_argument("--version", default=False, action="store_true", 
                        help="print this program's version")
    parser.add_argument("--list_tags", default=False, action="store_true", 
                        help="list the tags that exist on the server")
    parser.add_argument("--installed_tag", default=False, action="store_true", 
                        help="list the currently installed tag")
    parser.add_argument("--status", default=False, action="store_true", 
                        help="compare the main tag manifest to the files in the bundle directories")
    parser.add_argument("--update", default=False, action="store_true", 
                        help="update all installed bundles to the tag given")
    parser.add_argument("--diff", default=False, action="store_true", 
                        help="show the difference between the old and new text files")
    parser.add_argument("--difftool", default=None, 
                        help="show the difference between the old and new text files using diff tool such as meld")
    parser.add_argument("--create_change", default=False, action="store_true", 
                        help="create a tar file of new and modified files and a text file with a list of files to be deleted")
    parser.add_argument("--deliver_change", default=False, action="store_true", 
                        help="create and copy to /glusteruser/analyst/share_delivery/ a tar file of new and modified files and a text file with a list of files to be deleted, and email to swdevels@oceancolor.gsfc.nasa.gov from the email address stored in .install_ocssw")
    parser.add_argument("--deliver_email_from", default=None, 
                        help="create and copy to /glusteruser/analyst/share_delivery/ a tar file of new and modified files and a text file with a list of files to be deleted, and email to swdevels@oceancolor.gsfc.nasa.gov from the email address you provide")
    parser.add_argument("-t", "--tag", default=None,
                        help="tag that you want to install")
    parser.add_argument("-i", "--install_dir", default=os.environ.get("OCSSWROOT", None), 
                        help="root directory for bundle installation (default=$OCSSWROOT)")
    parser.add_argument("-b", "--base_url", default=baseUrl, 
                        help="remote url for the bundle server")
    parser.add_argument("-l", "--local_dir", default=None, 
                        help="local directory to use for bundle source instead of the bundle server")
    parser.add_argument("-s", "--save_dir", default=None, 
                        help="local directory to save a copy of the downloaded bundles")
    parser.add_argument("-c", "--clean", default=False, action="store_true", 
                        help="delete extra files in the destination directory")
    parser.add_argument("--wget", default=False, action="store_true", 
                        help="use wget for file download")
    parser.add_argument("-a", "--arch", default=None, 
                        help="use this architecture instead of guessing the local machine (linux_64,linux_hpc,macosx_intel,macosx_arm64,odps)")
    parser.add_argument("-v", "--verbose", action="count", default=0, help="increase output verbosity")

    # add weird bundle switches
    parser.add_argument("--bin", default=False, action="store_true", 
                        help="install binary executables")
    parser.add_argument("--opt", default=False, action="store_true", 
                        help="install 3rd party programs and libs")
    parser.add_argument("--src", default=False, action="store_true", 
                        help="install source files")
    parser.add_argument("--luts", default=False, action="store_true", 
                        help="install LUT files")
    parser.add_argument("--viirs_l1_bin", default=False, action="store_true", 
                        help="install VIIRS binary executables subset")

    # add bundles from the bundle list
    for bundleInfo in bundleList:
        if bundleInfo["commandLine"]:
            parser.add_argument("--" + bundleInfo["name"], default=False, action="store_true", 
                                help="install " + bundleInfo["help"] + " files")

    # add argument
    parser.add_argument("--direct_broadcast", default=False, action="store_true", 
                        help="toggle on bundles needed for MODIS direct broadcast")
    parser.add_argument("--seadas", default=False, action="store_true", 
                        help="toggle on the base set of bundles for SeaDAS")
    parser.add_argument("--odps", default=False, action="store_true", 
                        help="toggle on the base set of bundles for ODPS systems")
    parser.add_argument("--viirs_l1", default=False, action="store_true", 
                        help="install everything to run and test the VIIRS executables")
    parser.add_argument("--all", default=False, action="store_true", 
                        help="toggle on all satellite bundles")

    options = parser.parse_args()

    # print version
    if options.version:
        print(os.path.basename(sys.argv[0]), versionString)
        sys.exit(0)
    
    if options.list_tags:
        listTags(options)
        sys.exit(0)

    if options.installed_tag:
        installedTag(options)
        sys.exit(0)

    # make sure arch is set
    if not options.arch:
        if options.odps:
            options.arch = "odps"
        else:
            options.arch = getArch()

    if options.status:
        options.diff = False
        options.create_change = False
        if status(options):
            sys.exit(0)
        else:
            sys.exit(1)

    if options.difftool is not None:
        options.diff = True
    
    if options.diff:
        if status(options):
            sys.exit(0)
        else:
            sys.exit(1)
        
    if options.deliver_change or options.deliver_email_from is not None:
        options.create_change = True

    if options.create_change:
        if status(options):
            sys.exit(0)
        else:
            sys.exit(1)

    if not options.tag:
        print("--tag is required")
        sys.exit(1)

    if not options.install_dir:
        print("--install_dir is required if OCSSWROOT enviroment variable is not set")
        sys.exit(1)

    if options.update:
        update(options)
        sys.exit(0)

    # add convience arguments
    if options.benchmark:
        options.seadas = True
        options.modisa = True

    if options.viirs_l1:
        options.root = True
        options.viirs_l1_bin = True
        options.viirs_l1_benchmark = True
        options.opt = True
        options.luts = True
        options.common = True
        options.ocrvc = True
        options.viirsn = True
        options.viirsj1 = True
        options.viirsdem = True

    if options.direct_broadcast:
        options.seadas = True
        options.modisa = True
        options.modist = True
        
    if options.seadas:
        options.root = True
        options.bin =True
        options.opt = True
        options.luts = True
        options.common = True
        options.ocrvc = True
        
    if options.odps:
        options.root = True
        options.bin =True
        options.opt = True
        options.common = True
        options.ocrvc = True
        options.all = True
        
    if options.all:
        for bundleInfo in bundleList:
            if bundleInfo["commandLine"]:
                if "share/" in bundleInfo["dir"]:
                    setattr(options, bundleInfo["name"], True)

    # unset silly sensors for ODPS
    if options.odps:
        options.opt_src = False
        if hasattr(options, "benchmark"):
            options.benchmark = False
        if hasattr(options, "afrt"):
            options.afrt = False
        if hasattr(options, "aquaverse"):
            options.aquaverse = False
        options.avhrr = False
        options.aviris = False
        options.l5tm = False
        options.l7etmp = False
        if hasattr(options, "misr"):
            options.misr = False
        options.mos = False
        options.msi = False
        if hasattr(options, "msis2a"):
            options.msis2a = False
            options.msis2b = False
        options.ocm1 = False
        options.ocia = False
        options.ocip = False
        options.ocm2 = False
        options.oli = False
        if hasattr(options, "olil8"):
            options.olil8 = False
            options.olil9 = False
        options.osmi = False
        options.prism = False
        options.sabiamar = False
        options.sgli = False
        options.wv3 = False
        options.benchmark = False
        if hasattr(options, "viirs_l1_benchmark"):
            options.viirs_l1_benchmark = False

    # take care of the sub-sensors
    if options.modisa or options.modist:
        options.modis = True
    if hasattr(options, "msis2a"):
        if options.msis2a or options.msis2b:
            options.msi = True
    if options.viirsn or options.viirsj1:
        options.viirs = True
    if hasattr(options, "viirsj2"):
        if options.viirsj2:
            options.viirs = True
    if hasattr(options, "olcis3a"):
        if options.olcis3a or options.olcis3b:
            options.olci = True
    if hasattr(options, "olil8"):
        if options.olil8 or options.olil9:
            options.oli = True
            options.olil8 = True   # This is needed since L9 is a symbolic link to L8 LUTS

    # make sure bin and viirs_l1_bin are not both set
    if options.bin and options.viirs_l1_bin:
        print("Error: Can not install --bin and --viirs_l1_bin")
        sys.exit(1)

    # turn off binary options if installdir is a git repo
    if os.path.isdir(options.install_dir + "/.git"):
        if options.root or options.bin or options.opt or options.viirs_l1_bin:
            print("Turning off root, bin, opt and viirs_l1_bin since you are installing into a git repo")
        options.root = False
        options.bin = False
        options.opt = False
        options.viirs_l1_bin = False

    # count the things we are going to install
    totalNumThings = 0
    for bundleInfo in bundleList:
        if hasattr(options, bundleInfo["name"]) and getattr(options, bundleInfo["name"]):
            totalNumThings += 1
        
    # count bin and lib bundles
    if options.bin:
        totalNumThings += 2

    # count viirs_l1_bin and lib bundles
    if options.viirs_l1_bin:
        totalNumThings += 2

    # count the opt bundle
    if options.opt:
        totalNumThings += 1

    # count source bundles (ocssw-src and opt/src)
    if options.src:
        totalNumThings += 2

    # record the installed tag if any bundles are getting installed
    if totalNumThings > 0:
        recordInstalledTag(options)

    # now install the bundles

    # do the weird bundles
    if options.bin:
        bundleName = "bin_" + options.arch
        bundleInfo = findBundleInfo(bundleName)
        bundleInfo["dir"] = "bin"   # fix the install directory
        installBundle(options, bundleInfo)

    if options.viirs_l1_bin:
        bundleName = "viirs_l1_bin_" + options.arch
        bundleInfo = findBundleInfo(bundleName)
        if not bundleInfo:
            print("Error: tag does not contain the viirs_l1_bin bundles")
            sys.exit(1)
        bundleInfo["dir"] = "bin"   # fix the install directory
        installBundle(options, bundleInfo)

    if options.bin or options.viirs_l1_bin:
        bundleName = "lib_" + options.arch
        bundleInfo = findBundleInfo(bundleName)
        bundleInfo["dir"] = "lib"   # fix the install directory
        installBundle(options, bundleInfo)

    if options.opt:
        bundleName = "opt_" + options.arch
        bundleInfo = findBundleInfo(bundleName)
        bundleInfo["dir"] = "opt"   # fix the install directory
        installBundle(options, bundleInfo)

    if options.src:
        bundleInfo = findBundleInfo("opt_src")
        installBundle(options, bundleInfo)
        bundleInfo = findBundleInfo("ocssw_src")
        installBundle(options, bundleInfo)
        
    # do the normal bundles
    for bundleInfo in bundleList:
        if hasattr(options, bundleInfo["name"]) and getattr(options, bundleInfo["name"]):
            installBundle(options, bundleInfo)
    
    # update luts
    if options.luts:
        commonUpdated = False
        for bundleInfo in bundleList:
            if hasattr(options, bundleInfo["name"]) and getattr(options, bundleInfo["name"]):
                if bundleInfo["name"] in lutBundles:
                    updateLuts(options, bundleInfo["name"])
                    commonUpdated = True
        # make sure var/common gets updated
        if not commonUpdated:
            updateLuts(options, "common")

    print("Done\n")

if __name__ == "__main__":
    sys.exit(run())
