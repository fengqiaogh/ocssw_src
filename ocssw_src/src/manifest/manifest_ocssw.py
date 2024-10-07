#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 15:11:27 2018

@author: dshea
"""

import sys
import argparse
import os
import hashlib
import json
import subprocess
import shutil


MANIFEST_BASENAME = "manifest.json"
BUNDLELIST_BASENAME = "bundleList.json"

# configInfo contains global info about the ocssw_manifest
#   baseDir - base directory holding the manifest bundles
#   tag     - current installed tag
configInfo = {}
configInfoFilename = ".manifest_ocssw.json"
startingDir = ""

manifestCommand = os.path.dirname(__file__) + "/manifest.py"
manifestFilename = "manifest.json"

remoteHosts = ["ocssw@gs616-container106"]
remoteRepoDir = "/data1/manifest/tags"

##########################################################################
#  WARNING - The initialBundleList should match the latest bundleList.json
##########################################################################
# this is written to the bundle list json file when initializing a new manifest_ocssw dir
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

def add_subparser_init_dir(subparsers):
    newParser = subparsers.add_parser('init', help="Initialize the directory to be the root of the OCSSW manifest bundles")
    newParser.add_argument("dir", help="directory to initialize", nargs='?', default='.')
    newParser.set_defaults(func=init_dir)

def add_subparser_diff(subparsers):
    newParser = subparsers.add_parser('diff', help="Show the difference between the old and new text files")
    newParser.add_argument("tag", nargs='?', help="tag name to use")
    newParser.set_defaults(func=diff) 
    newParser.set_defaults(difftool=None)

def add_subparser_difftool(subparsers):
    newParser = subparsers.add_parser('difftool', help="Show the difference between the old and new text files using diff tool such as meld")
    newParser.add_argument("difftool", nargs='?', help="tool name to use for diff")
    newParser.add_argument("tag", nargs='?', help="tag name to use")
    newParser.set_defaults(func=diff)   

def add_subparser_status(subparsers):
    newParser = subparsers.add_parser('status', help="List the files that differ from the manifest")
    newParser.add_argument("tag", nargs='?', help="tag name to use")
    newParser.set_defaults(func=status)

def add_subparser_push(subparsers):
    newParser = subparsers.add_parser('push', help="Create a new tag and push to the server")
    newParser.add_argument("tag", help="tag name to use")
    newParser.add_argument("-o", "--overwrite", default=False, action="store_true", help="overwrite the tag if it already exists")
    newParser.set_defaults(func=push)

def add_subparser_pull(subparsers):
    newParser = subparsers.add_parser('pull', help="Download a tag from the server")
    newParser.add_argument("tag", help="tag name to use")
    newParser.add_argument("-c", "--clean", default=False, action="store_true", help="clean the directory after downloading")
    newParser.set_defaults(func=pull)

def add_subparser_list(subparsers):
    newParser = subparsers.add_parser('list', help="List the tags on the remote system")
    newParser.set_defaults(func=list_tags)


def init_dir(options):
    global configInfo

    if os.path.exists(options.dir):
        if not os.path.isdir(options.dir):
            print("Error:", options.dir, "exists and is not a directory")
    else:
        os.mkdir(options.dir)

    fileName = os.path.join(os.path.abspath(options.dir), configInfoFilename)
    if os.path.exists(fileName):
        print("Error: manifest_ocssw config file", fileName, "already exists")
        sys.exit(1)

    # fill up an initial config structure
    configInfo["baseDir"] = os.path.abspath(options.dir)
    configInfo["tag"] = "Unknown"

    # save as a JSON file
    with open(fileName, 'w') as outfile:
        json.dump(configInfo, outfile, indent=4, sort_keys=True)

    fileName = os.path.join(os.path.abspath(options.dir), "bundleList.json")
    with open(fileName, 'w') as outfile:
        json.dump(initialBundleList, outfile, indent=4, sort_keys=True)


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


def tagExists(tag):
    command = ["ssh", "-q", remoteHosts[0], "[ -d %s/%s ]" % (remoteRepoDir, tag) ]
    proc = subprocess.run(command)
    if proc.returncode == 0:
        return True
    else:
        return False
    
    
def bundleExists(tag, name):
    command = ["ssh", "-q", remoteHosts[0], "[ -d %s/%s/%s ]" % (remoteRepoDir, tag, name) ]
    proc = subprocess.run(command)
    if proc.returncode == 0:
        return True
    else:
        return False

def findBundleInList(name, bundleList):
    if bundleList:
        for info in bundleList:
            if info["name"] == name:
                return info
    return None

def diff(options, bundleInfo):
    global startingDir

    # only do a status if user is sitting in the bundel directory
    bundleDir = os.path.abspath(os.path.join(configInfo["baseDir"], bundleInfo["dir"]))
    if startingDir not in bundleDir:
        return

    print(bundleInfo["name"], end=' ', flush=True)

    statusTag = "tempStatusTag"
    os.chdir(configInfo["baseDir"])

    currentManifest = {}
    if hasattr(options, "tag") and options.tag:
        if tagExists(options.tag):
            if bundleExists(options.tag, bundleInfo["name"]):
                command = ["scp", "-q", "%s:%s/%s/%s/%s" % (remoteHosts[0], remoteRepoDir, options.tag, bundleInfo["name"], manifestFilename), "/tmp"]
                runCommand(command)
                with open("/tmp/" + manifestFilename, 'rb') as manifest:
                    currentManifest = json.load(manifest)
                os.remove("/tmp/" + manifestFilename)    
            else:
                print("Warning: bundle %s does not exist in tag %s on server."  % (bundleInfo["name"], options.tag))
        else:
            print("Error: tag:", options.tag, "does not exist on server")
            sys.exit(1)
    else:
        # open current manifest
        if bundleExists(configInfo["tag"], bundleInfo["name"]):
            command = ["scp", "-q", "%s:%s/%s/%s/%s" % (remoteHosts[0], remoteRepoDir, configInfo["tag"], bundleInfo["name"], manifestFilename), "/tmp"]
            runCommand(command)
            with open("/tmp/" + manifestFilename, 'rb') as manifest:
                currentManifest = json.load(manifest)
            os.remove("/tmp/" + manifestFilename)    

    if not os.path.isdir(bundleInfo["dir"]):
        print("Warning: bundle directory %s does not exist."  % (bundleInfo["dir"]))
        return
        
    command = [manifestCommand, "generate", 
               "-n", bundleInfo["name"], 
               "-t", statusTag,
               bundleInfo["dir"]]
    if "extra" in bundleInfo:
        command += bundleInfo["extra"].split()
    statusManifest = json.loads(runCommand(command, pipe=True))

    if "tags" in currentManifest:
        tagList = currentManifest["tags"]
        versionStr = "("
        if len(tagList) > 1:
            versionStr += tagList[0] + ".."
        versionStr += tagList[-1] + ")"
        print(versionStr)
    else:
        print()
        
    # list new files
    fileList = []
    for f, info in statusManifest["files"].items():
        if ("files" not in currentManifest) or (f not in currentManifest["files"]):
            fileList.append(f)
    if fileList:
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
                # remote symlink 
                else: 
                    fileList.append(f)
            # local symlink  
            elif "symlink" in info:
                    fileList.append(f)                   
            elif info["tag"] != currentManifest["files"][f]["tag"]:
                localFilename = "%s/%s" % (bundleInfo["dir"], f)
                with open(localFilename, "rb") as localFile:
                    bytes = localFile.read()
                    localFileChecksum = hashlib.sha256(bytes).hexdigest() 
                    if "checksum" in currentManifest["files"][f]:    
                        if localFileChecksum != currentManifest["files"][f]["checksum"]:
                                fileList.append(f)   

    if fileList:
        print("  Modified Files:")
        for f in fileList:
            print("    %s/%s" % (bundleInfo["dir"], f))
            if "symlink" in currentManifest["files"][f]:
                # both symlink
                if "symlink" in statusManifest["files"][f]:     
                    print("diff %s remote..local" % (f))
                    print("<      " + currentManifest["files"][f]["tag"] + "     = " + currentManifest["files"][f]["symlink"])
                    print("---")
                    print(">      current dir = " + statusManifest["files"][f]["symlink"])
                # remot symlink
                else:
                    print("diff %s remote..local" % (f))
                    print("<      " + currentManifest["files"][f]["tag"] + "     = " + currentManifest["files"][f]["symlink"])
                    print("---")
                    print(">      current dir = Not a symlink")
            # local symlink
            elif "symlink" in statusManifest["files"][f]:
                print("diff %s remote..local" % (f))
                print("<      " + currentManifest["files"][f]["tag"] + "     = Not a symlink")
                print("---")
                print(">      current dir = " + statusManifest["files"][f]["symlink"])
            # neither symlink
            else:
                try:
                    localFilename = "%s/%s" % (bundleInfo["dir"], f)
                    subproc = subprocess.run(["file", '-L', localFilename], capture_output=True, text=True, shell=False)
                    if "ASCII" in subproc.stdout:
                        remoteFilename = "%s/%s" % ("/tmp/", os.path.basename(f))
                        file_tag = currentManifest["files"][f]["tag"]
                        command = ["scp", "-q", "%s:%s/%s/%s/%s" % (remoteHosts[0], remoteRepoDir, file_tag, bundleInfo["name"], f), "/tmp"]             
                        status = subprocess.run(command, shell=False).returncode
                        if status == 0:
                            if options.difftool is not None:
                                command = [options.difftool, remoteFilename, localFilename]
                                subprocess.run(command, shell=False)
                            else:
                                print("diff %s remote..local" % (f))
                                command = ["diff",remoteFilename, localFilename]
                                subprocess.run(command, shell=False)
                        os.remove(remoteFilename)
                except:
                    continue



    # deleted file     
    fileList = []
    if "files" in currentManifest:
        for f, info in currentManifest["files"].items():
            if f not in statusManifest["files"]:
                fileList.append(f)
    if fileList:
        print("  Deleted Files:")
        for f in fileList:
            print("    %s/%s" % (bundleInfo["dir"], f))

def status(options, bundleInfo):
    global startingDir

    # only do a status if user is sitting in the bundel directory
    bundleDir = os.path.abspath(os.path.join(configInfo["baseDir"], bundleInfo["dir"]))
    if startingDir not in bundleDir:
        return

    print(bundleInfo["name"], end=' ', flush=True)

    statusTag = "tempStatusTag"
    os.chdir(configInfo["baseDir"])

    currentManifest = {}
    if hasattr(options, "tag") and options.tag:
        if tagExists(options.tag):
            if bundleExists(options.tag, bundleInfo["name"]):
                command = ["scp", "-q", "%s:%s/%s/%s/%s" % (remoteHosts[0], remoteRepoDir, options.tag, bundleInfo["name"], manifestFilename), "/tmp"]
                runCommand(command)
                with open("/tmp/" + manifestFilename, 'rb') as manifest:
                    currentManifest = json.load(manifest)
                os.remove("/tmp/" + manifestFilename)    
            else:
                print("Warning: bundle %s does not exist in tag %s on server."  % (bundleInfo["name"], options.tag))
        else:
            print("Error: tag:", options.tag, "does not exist on server")
            sys.exit(1)
    else:
        # open current manifest
        if bundleExists(configInfo["tag"], bundleInfo["name"]):
            command = ["scp", "-q", "%s:%s/%s/%s/%s" % (remoteHosts[0], remoteRepoDir, configInfo["tag"], bundleInfo["name"], manifestFilename), "/tmp"]
            runCommand(command)
            with open("/tmp/" + manifestFilename, 'rb') as manifest:
                currentManifest = json.load(manifest)
            os.remove("/tmp/" + manifestFilename)    

    if not os.path.isdir(bundleInfo["dir"]):
        print("Warning: bundle directory %s does not exist."  % (bundleInfo["dir"]))
        return
        
    command = [manifestCommand, "generate", 
               "-n", bundleInfo["name"], 
               "-t", statusTag,
               bundleInfo["dir"]]
    if "extra" in bundleInfo:
        command += bundleInfo["extra"].split()
    statusManifest = json.loads(runCommand(command, pipe=True))

    if "tags" in currentManifest:
        tagList = currentManifest["tags"]
        versionStr = "("
        if len(tagList) > 1:
            versionStr += tagList[0] + ".."
        versionStr += tagList[-1] + ")"
        print(versionStr)
    else:
        print()
        
    # list new files
    fileList = []
    for f, info in statusManifest["files"].items():
        if ("files" not in currentManifest) or (f not in currentManifest["files"]):
            fileList.append(f)
    if fileList:
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
                # remote symlink 
                else: 
                    fileList.append(f)
            # local symlink  
            elif "symlink" in info:
                    fileList.append(f) 
            # neither symlink                  
            elif info["tag"] != currentManifest["files"][f]["tag"]:
                localFilename = "%s/%s" % (bundleInfo["dir"], f)
                with open(localFilename, "rb") as localFile:
                    bytes = localFile.read()
                    localFileChecksum = hashlib.sha256(bytes).hexdigest() 
                    if "checksum" in currentManifest["files"][f]:    
                        if localFileChecksum != currentManifest["files"][f]["checksum"]:
                                fileList.append(f)   
    if fileList:
        print("  Modified Files:")
        for f in fileList:
            print("    %s/%s" % (bundleInfo["dir"], f))

    # deleted file     
    fileList = []
    if "files" in currentManifest:
        for f, info in currentManifest["files"].items():
            if f not in statusManifest["files"]:
                fileList.append(f)
    if fileList:
        print("  Deleted Files:")
        for f in fileList:
            print("    %s/%s" % (bundleInfo["dir"], f))

def push(options, bundleInfo):
    os.chdir(configInfo["baseDir"])
    command = [manifestCommand, "generate", 
               "-n", bundleInfo["name"], 
               "-t", options.tag,
               bundleInfo["dir"]]
    if "extra" in bundleInfo:
        command += bundleInfo["extra"].split()
    newManifest = json.loads(runCommand(command, pipe=True))

    print(bundleInfo["name"])
    os.chdir(bundleInfo["dir"])
 
    # list files to upload
    fileList = []
    for f, info in newManifest["files"].items():
        if info["tag"] == options.tag:
            fileList.append(f)

    # check if files were deleted
    filesWereDeleted = False
    if bundleExists(configInfo["tag"], bundleInfo["name"]):
        command = ["scp", "-q", "%s:%s/%s/%s/%s" % (remoteHosts[0], remoteRepoDir, configInfo["tag"], bundleInfo["name"], manifestFilename), "/tmp"]
        runCommand(command)
        with open("/tmp/" + manifestFilename, 'rb') as manifest:
            currentManifest = json.load(manifest)
    
            for f, info in currentManifest["files"].items():
                if f not in newManifest["files"]:
                    filesWereDeleted = True
                    break
        os.remove("/tmp/" + manifestFilename)

    if fileList or filesWereDeleted:
	# save manifest file
        with open(manifestFilename, 'w') as outfile:
            json.dump(newManifest, outfile, indent=4, sort_keys=True)

        # make tmp file with list of files to upload        
        tmpFilename = "/tmp/changeList.txt"
        with open(tmpFilename, 'w') as changeFile:
            changeFile.write(manifestFilename)
            changeFile.write("\n")
            for f in fileList:
                if not os.path.islink(f):
                    changeFile.write(f)
                    changeFile.write("\n")

        # upload files
        for host in remoteHosts:
            command = ["ssh", "-q",  host, "mkdir -p %s/%s/%s" % (remoteRepoDir, options.tag, bundleInfo["name"]) ]
            runCommand(command)

            command = ["rsync", "-hav", "--files-from=" + tmpFilename, ".", "%s:%s/%s/%s/" % (host, remoteRepoDir, options.tag, bundleInfo["name"]) ]
            runCommand(command)
        os.remove(tmpFilename)
                
    else:
        print("  No files Changed, creating symbolic link")
        
        # add tag to manifest
        command = [manifestCommand, "add-tag", options.tag]
        newManifest = json.loads(runCommand(command, pipe=True))
        with open(manifestFilename, 'w') as outfile:
            json.dump(newManifest, outfile, indent=4, sort_keys=True)

        firstTag = newManifest['tags'][0]
        
        # update remote manifest and make symbolic link
        for host in remoteHosts:
            command = ["ssh", "-q",  host, "mkdir -p %s/%s" % (remoteRepoDir, options.tag)]
            runCommand(command)
            
            command = [ "scp", "-q", manifestFilename, "%s:%s/%s/%s/%s" % (host, remoteRepoDir, firstTag, bundleInfo["name"], manifestFilename) ]
            runCommand(command)

            command = ["ssh", "-q", host, "ln -s ../%s/%s %s/%s/%s" %
                       (firstTag,    bundleInfo["name"],
                        remoteRepoDir, options.tag, bundleInfo["name"]) ]
            runCommand(command)

def pull(options, bundleInfo):
    print(bundleInfo["name"])
    if bundleExists(options.tag, bundleInfo["name"]):
        os.chdir(configInfo["baseDir"])
        command = [manifestCommand, "download",
                   "-n", bundleInfo["name"], 
                   "-t", options.tag,
                   "-d", bundleInfo["dir"]]
        runCommand(command)
        if options.clean:
            command = [manifestCommand, "clean", bundleInfo["dir"]]
            if "extra" in bundleInfo:
                command += bundleInfo["extra"].split()
            runCommand(command)
    else:
        print("Warning: Bundle %s does not exist in tag %s on the server" %(bundleInfo["name"], options.tag))
        if options.clean:
            if os.path.isdir(bundleInfo["dir"]):
                print("Cleaning bundle directory", bundleInfo["dir"])
                shutil.rmtree(bundleInfo["dir"])
    
def list_tags(options):
    command = ["ssh", "-q", remoteHosts[0], "ls " + remoteRepoDir]
    runCommand(command)


def run():
    global configInfo
    global startingDir
    
    parser = argparse.ArgumentParser(description="Work on the minifest directories for OCSSW")
    parser.set_defaults(func=status)
    subparsers = parser.add_subparsers()

    add_subparser_init_dir(subparsers)
    add_subparser_diff(subparsers)
    add_subparser_difftool(subparsers)
    add_subparser_status(subparsers)
    add_subparser_push(subparsers)
    add_subparser_pull(subparsers)
    add_subparser_list(subparsers)

    options = parser.parse_args()
    
    # check if we are initializing a new OCSSW manifest
    if options.func == init_dir:
        options.func(options)
        return 0
    
    # check if we are listing the tags
    if options.func == list_tags:
        options.func(options)
        return 0
    
    # if making a new tag, make sure tag does not exist or the --overwite flag is used
    if options.func == push:
        if tagExists(options.tag):
            if options.overwrite:
                for host in remoteHosts:
                    command = ["ssh", "-q", host, "rm -rf " + remoteRepoDir + "/" + options.tag]
                    runCommand(command)
            else:
                print("Error: tag", options.tag, "exists.")
                print("Use --overwrite (-o) to continue (will delete old version)")
                return 1
        
    # find the config file
    baseDir = os.path.abspath(os.getcwd())
    startingDir = baseDir
    while(baseDir != "/"):
        fileName = os.path.join(baseDir, configInfoFilename)
        if os.path.exists(fileName):
            with open(fileName, 'rb') as configFile:
                configInfo = json.load(configFile)
                if configInfo["baseDir"] != baseDir:
                    print("Error: config file has been moved")
                    print("    real baseDir =      ", baseDir)
                    print("    configFile baseDir =", configInfo["baseDir"])
                    return 1
                else:
                    break
        baseDir = os.path.dirname(baseDir)
    if baseDir == "/":
        print("Error: manifest_ocssw config file not found in current or parent directories")
        print("       Run: manifest_ocssw init")
        return 1
            
    print("current tag =", configInfo["tag"])

    #load local bundleList json file
    localBundleList = None
    localBundleListFile = "%s/%s" % (baseDir, BUNDLELIST_BASENAME)
    if not os.path.isfile(localBundleListFile):
        print("Must have a %s file in the root directory" % (BUNDLELIST_BASENAME))
        sys.exit(1)
    with open(localBundleListFile, 'rb') as jsonFile:
        localBundleList = json.load(jsonFile)
    if not localBundleList:
        print("Must have a valid %s file in the root directory" % (BUNDLELIST_BASENAME))
        sys.exit(1)

    # download remote bundleList json files
    remoteBundleList = None
    remoteBundleListFile = "/tmp/%s" % (BUNDLELIST_BASENAME)
    remoteManifest = None
    remoteManifestFile = "/tmp/%s" % (MANIFEST_BASENAME)
    tag = configInfo["tag"]
    if options.func == pull:
        tag = options.tag
    elif options.func == status:
        if hasattr(options, "tag") and options.tag:
            tag = options.tag

    if tagExists(tag):
        # download manifest file first
        command = ["scp", "-q", "%s:%s/%s/root/%s" % 
            (remoteHosts[0], remoteRepoDir, tag, MANIFEST_BASENAME), 
            remoteManifestFile]
        runCommand(command)
        with open(remoteManifestFile, 'rb') as jsonFile:
            remoteManifest = json.load(jsonFile)
        os.remove(remoteManifestFile)

        # download bundleList
        command = ["scp", "-q", "%s:%s/%s/root/%s" % 
                   (remoteHosts[0], remoteRepoDir, remoteManifest['files'][BUNDLELIST_BASENAME]['tag'], 
                    BUNDLELIST_BASENAME), remoteBundleListFile]
        runCommand(command)
        with open(remoteBundleListFile, 'rb') as jsonFile:
            remoteBundleList = json.load(jsonFile)
        os.remove(remoteBundleListFile)
        
    # check for new local bundles
    for bundleInfo in localBundleList:
        if not findBundleInList(bundleInfo["name"], remoteBundleList):
            print("New local bundle -", bundleInfo["name"])

    # check for deleted bundles
    if remoteBundleList:
        for bundleInfo in remoteBundleList:
            if not findBundleInList(bundleInfo["name"], localBundleList):
                print("Deleted local bundle -", bundleInfo["name"])

    bundleList = localBundleList
    if options.func == pull:
        bundleList = remoteBundleList

    # actually run the function
    #for bundleInfo in bundleList:
    for bundleInfo in bundleList:
        options.func(options, bundleInfo)

    # update the tag in the configInfo
    if options.func == push or options.func == pull:
        configInfo["tag"] = options.tag
        
        # save as a JSON file
        fileName = os.path.join(baseDir, configInfoFilename)
        with open(fileName, 'w') as outfile:
            json.dump(configInfo, outfile, indent=4, sort_keys=True)
        
    return 0


if __name__ == "__main__":
    sys.exit(run())
    
