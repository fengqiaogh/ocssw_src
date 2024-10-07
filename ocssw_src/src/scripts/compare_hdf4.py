#!/usr/bin/env python3

"""
Created on Thu Sep 29 12:19:06 2016

@author: dshea
"""

import argparse
import subprocess
import sys


#
# make an item class
#
class DiffItem:
    
    def __init__(self, kind, name):
        self.kind = kind
        self.name = name
        self.lines = []
        
    def addLine(self, line):
        self.lines.append(line)

    def printItem(self):
        for line in self.lines:
            print(line)


startItems = {}
dataItems = {}
globalAttributeItems = {}


def addItem(item):
    if item.lines:
        if item.kind == "start":
            startItems[item.name] = item
        if item.kind == "data":
            dataItems[item.name] = item
        elif item.kind == "globalAttribute":
            globalAttributeItems[item.name] = item

#
# start of main program
#
parser = argparse.ArgumentParser(description='Compare two HDF4 files.')

parser.add_argument('-G', '--globalExclude', help='Global Attrbutes to ignore')
parser.add_argument('-p', '--relativeLimit', help='relative limit, criteria is |(b-a)/a| > relativeLimit')
parser.add_argument('file1', help='first HDF4 file to compare')
parser.add_argument('file2', help='second HDF4 file to compare')

args = parser.parse_args()

if args.relativeLimit == None:
    proc = subprocess.Popen(["hdiff", "-e", "10",  args.file1, args.file2], stdout=subprocess.PIPE)
else:
    proc = subprocess.Popen(["hdiff", "-e", "10", "-p", args.relativeLimit, args.file1, args.file2], stdout=subprocess.PIPE)

out = proc.communicate()[0]

# if hdiff returns 0 we are done
if(proc.returncode == 0):
    sys.exit(0)

if(proc.returncode != 1):
    print(out.decode())
    print("hdiff error return = ", proc.returncode)
    sys.exit(1)
    
# init the first item
item = DiffItem('start', 'first item')

# run through all the lines and create DiffItems
lines = out.split(b'\n')
for line in lines:
    line = line.decode()  # change from binary string into normal ascii string
    line = line.strip()
        
    if len(line) > 0:
        if line.startswith("---------------------------"):
            continue
        elif line.startswith("position"):
            addItem(item)
            item = DiffItem('data', line.split()[1])
            item.addLine(line)
        elif line.startswith("Attr Name:"):
            addItem(item)
            item = DiffItem('globalAttribute', line.split(None,2)[2])
            item.addLine(line)
        else:
            item.addLine(line)

# add the last item
addItem(item)
                    
# run through all of the exclude items and delete thier corresponding DiffItem
if args.globalExclude != None:
    excludes = args.globalExclude.split(",")
    for exclude in excludes:
        exclude = exclude.strip()
        if exclude in globalAttributeItems:
            del globalAttributeItems[exclude]

# set initial exit code
exitCode = 0

# print out everything
if len(startItems) > 0:
    for key,value in startItems.items():
        print('---------------------------')
        value.printItem() 
    
if len(dataItems) > 0:
    exitCode = 1
    for key,value in dataItems.items():
        print('---------------------------')
        value.printItem() 
        
if len(globalAttributeItems) > 0:
    exitCode = 1
    for key,value in globalAttributeItems.items():
        print("---------------------------")
        value.printItem() 
    

sys.exit(exitCode)

