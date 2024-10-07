#!/usr/bin/env python3
#
import sys
import os
import time
import datetime
import netCDF4 as nc
import numpy as np

if len(sys.argv) != 3:
    print('Usage:')
    callseq = sys.argv[0] + ' file ' + 'attr-file'
    print(callseq)
    print('\nThis script modifies attributes in a netCDF file')
    print('\tfile:\t\tsource file to modify')
    print('\tattr-file:\tfile containig attribue=value pairs to modify\n')
    print('\tTo delete an attribute, set the value to a blank string\n')
    print('\tAttributes listed in the attr-file but not existing in the netCDF file\n\twill be added')
    print('\tFor group attributes, prepend the attribute name with the group,\n\tdelineated with a "/", e.g.:')
    print('\t\tnavigation_data/mycoolattribute=awesomeness')
    sys.exit(1)

inFile = sys.argv[1]
attrFile = sys.argv[2]
rootgrp = nc.Dataset(inFile, 'a')

# generate string to update the history
processtime = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%dT%H:%M:%S')
if 'history' in rootgrp.ncattrs():
    history = nc.Dataset.getncattr(rootgrp,'history').rstrip() + '\n'+'['+processtime+'] '+os.path.basename(sys.argv[0])+ ' ' + sys.argv[1] + ' ' + sys.argv[2]
else:
    history = '['+processtime+'] ' + os.path.basename(sys.argv[0]) + ' ' + sys.argv[1] + ' ' + sys.argv[2]

attrs = open(attrFile,'r')

for attr in attrs:
    attrName,attrValue = attr.split('=')
    attrValue = attrValue.rstrip()
    group = rootgrp
    try:
        parts = attrName.split('/')
        if len(parts) > 1:
            attrName = parts.pop()
            attrGroup = '/'.join(parts).strip()
            group = rootgrp.createGroup(attrGroup)
    except:
        pass

    # Check to see if the attribute exists,
    # if it does, delete it
    try:
        attrValueExisting = nc.Dataset.getncattr(group,attrName)
        nc.Dataset.delncattr(group,attrName)
    except:
        pass

    # if starts with a quote it is a string
    try:
        if len(attrValue):
            if attrValue.startswith('"'):
                attrValue = attrValue.strip('"').encode('ascii')
                nc.Dataset.setncattr(group,attrName,attrValue)
            else:
                # if has a comma it is an array
                if "," in attrValue:
                    if attrValue.startswith('['):
                        attrValue = attrValue.lstrip('[').rstrip(']')
                    partsStr = attrValue.split(',')
                    if '.' in attrValue:
                        parts = np.array(partsStr).astype(np.dtype('f4'))
                    else:
                        parts = np.array(partsStr).astype(np.dtype('i4'))
                    nc.Dataset.setncattr(group,attrName,parts)

                else:
                    if '.' in attrValue:
                        nc.Dataset.setncattr(group,attrName,np.array(attrValue).astype(np.dtype('f4')))
                    else:
                        nc.Dataset.setncattr(group,attrName,np.array(attrValue).astype(np.dtype('i4')))
    except Exception as e:
        sys.exit("Whoops! Something went seriously wrong:\n"+ str(e))
        
# update the history
nc.Dataset.setncattr(rootgrp,'history', history)

        
rootgrp.close()
