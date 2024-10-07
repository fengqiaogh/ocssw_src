#!/usr/bin/env python3
#
import sys
from netCDF4 import Dataset

if len(sys.argv) != 2:
    print('usage:', sys.argv[0], '<L3b file>')
    print('    convert an L3 bin file with the following products to')
    print('    an OCRVC sensor L3 bin file:')
    print('        Rrs_vc_412')
    print('        Rrs_vc_443')
    print('        Rrs_vc_490')
    print('        Rrs_vc_510')
    print('        Rrs_vc_555')
    print('        Rrs_vc_670')
    sys.exit(1)

inFile = sys.argv[1]
rootgrp = Dataset(inFile, 'a')
datagrp = rootgrp.groups['level-3_binned_data']

# make sure it is a L3b file
if 'Level-3 Binned Data' not in rootgrp.title:
    print(inFile, 'is not a L3 bin file')
    sys.exit(1)

# change global attrbutes
rootgrp.title = 'OCRVC Level-3 Binned Data'
rootgrp.instrument = 'OCRVC'
rootgrp.platform = 'OCRVC'
#rootgrp.history += ', Convert to OCRVC sensor'

datagrp.renameVariable('Rrs_vc_412','Rrs_412')
datagrp.renameVariable('Rrs_vc_443','Rrs_443')
datagrp.renameVariable('Rrs_vc_490','Rrs_490')
datagrp.renameVariable('Rrs_vc_510','Rrs_510')
datagrp.renameVariable('Rrs_vc_531','Rrs_531')
datagrp.renameVariable('Rrs_vc_555','Rrs_555')
datagrp.renameVariable('Rrs_vc_670','Rrs_670')

rootgrp.close()

print('Converted', inFile, 'to an OCRVC bin file.')

sys.exit(0)
