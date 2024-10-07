#! /usr/bin/env python3

import os
import sys
import argparse
from pathlib import Path

version = "2.0"
seperatorStr = "=============================================="
local = False

def read_next_test(fp):
    line = ""

    # init the output dictionary
    out = {
        "valid": False,     # did we fill up the record
        "failed": False,    # did the test Pass of Fail
        "name": "",         # name of the test
        "log": [],          # all of the lines of the test output
        "command": "",      # test command line
        "dir": "",          # dir the test was run from
        "cmpDir": "",       # dir nccmp was run from
        "cdataDir": "",     # dir to run cdata from
        "newFile": [],      # new files just created
        "oldFile": []       # old test files use as comparison
    }

    # find first line of test
    while True:
        line = fp.readline()
        if not line:
            return False
        if " Testing: " in line:
            parts = line.split()
            if len(parts) == 3:
                out["name"] = parts[2]
                break

    # add log header
    out["log"].append(seperatorStr)
    out["log"].append(seperatorStr)
    out["log"].append(seperatorStr)
    
    # search for command
    while True:
        line = fp.readline()
        if not line:
            return False
        out["log"].append(line)
        if line.startswith("Command: "):
            parts = line.split(' ', 1)
            out["command"] = parts[1]
            break;

    # search for dir
    while True:
        line = fp.readline()
        if not line:
            return False
        out["log"].append(line)
        if line.startswith("Directory: "):
            parts = line.split()
            out["dir"] = parts[1]
            break;

    # find new and old filenames
    found = False
    out["cmpDir"] = out["dir"]
    parts = out["command"].split('&&')
    for cmd in parts:
        cmd = cmd.strip()
        if cmd.startswith("cd "):
            cmdParts = cmd.split()
            newDir = cmdParts[1]
            if newDir[0] != '/':
                newPath = Path(out["cmpDir"]) / newDir
                newDir = str(newPath.resolve())
            out["cmpDir"] = newDir
        elif cmd.startswith("nccmp ") or cmd.startswith("compare_hdf4 ") or cmd.startswith("diff ") or cmd.startswith("tiffcmp ") :
            cmdParts = cmd.split()
            out["newFile"].append(out["cmpDir"] + "/" + cmdParts[-1].strip().rstrip('"'))
            out["oldFile"].append(out["cmpDir"] + "/" + cmdParts[-2].strip())
            found = True
        elif cmd.startswith("compare "):
            cmdParts = cmd.split()
            out["newFile"].append(out["cmpDir"] + "/" + cmdParts[-2].strip())
            out["oldFile"].append(out["cmpDir"] + "/" + cmdParts[-3].strip())
            found = True
    if not found:
        return out

    # figure out the cdata dir
    line = out["cmpDir"]
    parts = line.split("/")
    if parts[-2] == "testdata":
        out["cdataDir"] = line
    elif parts[-3] == "testdata":
        parts.pop(-1)
        out["cdataDir"] = "/".join(parts)
    else:
        print("ERROR: can't find cdataDir")
        print("    name   =", out["name"])
        print("    cmpDir =", out["cmpDir"])
        print()
        return out

    # delete cdataDir from old and new filenames
    s = out["cdataDir"] + "/"
    for i in range(len(out["newFile"])):
        out["newFile"][i] = out["newFile"][i].replace(s, "")
        out["oldFile"][i] = out["oldFile"][i].replace(s, "")

    # search for passed or failed
    lastLine = ""
    while True:
        line = fp.readline()
        if not line:
            return False
        out["log"].append(line)
        if lastLine.startswith("-----------------------------------------"):
            if line == "Test Passed.\n":
                out["failed"] = False
                break;
            if line == "Test Failed.\n":
                out["failed"] = True
                break;
        lastLine = line

    out["valid"] = True
    return out


def runCmd(cmd):
    print("+" + cmd)
    os.system(cmd)


def update(fp):
    cdataSet = set()
    while(True):
        # read the next test record
        rec = read_next_test(fp)
        if not rec:
            break
        
        if rec["valid"] and rec["failed"]:

            cmdList = []
            cmdList.append("cd " + rec["cdataDir"])
            for i in range(len(rec["newFile"])):
                cmdList.append("cp " + rec["newFile"][i] + " " + rec["oldFile"][i])
                if not local:
                    cmdList.append("cdata add " + rec["oldFile"][i])

            for line in rec["log"]:
                print(line.rstrip())
            for line in cmdList:
                print("+" + line)
            print(seperatorStr)

            val = input("update(u)/ignore(i) (defaut=u): ")
            if val == '' or val == 'u' or val == 'U':
                print("updating")
                runCmd("; ".join(cmdList))
                cdataSet.add(rec["cdataDir"])
            else:
                print("ignoring")

    print()
    if not local:
        for line in list(cdataSet):
            runCmd("cd " + line + "; cdata push")
    return 0


def listFailed(fp):
    while(True):
        # read the next test record
        rec = read_next_test(fp)
        if not rec:
            break
        
        if rec["valid"] and rec["failed"]:
            for line in rec["log"]:
                print(line.rstrip())
    return 0


def main():
    global local

    parser = argparse.ArgumentParser(prog="cdata-log")
    parser.add_argument('--version', action='version', version='%(prog)s ' + version)
    parser.add_argument('--list-failed', action='store_true', default=False,
                        help="list log of failed tests")
    parser.add_argument('--local', action='store_true', default=False,
                        help="only copy updated test files to local dir")

    args = parser.parse_args()

    filename = "Testing/Temporary/LastTest.log"

    if not os.path.exists(filename):
        print("ERROR:", filename, "does not exist")
        print("ERROR: cdata-log must be run from the same directory as ctest")
        return 1

    if args.local:
        local = True

    with open(filename) as fp:
        if args.list_failed:
            return listFailed(fp)
        else:
            return update(fp)

if __name__ == "__main__":
    sys.exit(main())
