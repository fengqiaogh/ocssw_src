#! /usr/bin/env python3

import numpy as np
from scipy import ndimage
from pyhdf.SD import *

from datetime import datetime
from shutil import copyfile
import re

'''
G. Fireman May 2017; adapted from Wayne Robinson's met_reanal_full.pro
'''

NY_OUT = 181  # 'Number of Rows'
NX_OUT = 360  # 'Number of Columns'

PARAM_RANGES = {
   'z_wind' : [-50,   50],
   'm_wind' : [-50,   50],
   'press'  : [850, 1084],
   'rel_hum': [  0,  100],
   'p_water': [  0,  200]
}
params = PARAM_RANGES.keys()


def enforce_range(values, param):
   try:
      limits = PARAM_RANGES[param]
      np.clip(values, limits[0], limits[1], out=values)
   except Exception as e:
      print(e)
      print('Unknown parameter referenced in enforce_range: ', param)


def smooth_latlon(old_arr, width=3):
   # y axis = latitude: truncate at poles
   new_arr = ndimage.uniform_filter1d(old_arr, width, axis=0, mode='nearest')
   # x axis = longitude: wrap around 180
   new_arr = ndimage.uniform_filter1d(new_arr, width, axis=1, mode='wrap')
   return new_arr


def fix_bias(old_dat, new_dat, width=3):

   # resize new data to same resolution as old
   zoom = np.divide(old_dat.shape, new_dat.shape)
   new_dat2 = ndimage.interpolation.zoom(new_dat, zoom, order=1, mode='nearest')

   # difference smoothed fields to find bias
   bias = smooth_latlon(old_dat, width) - smooth_latlon(new_dat2, width)

   # return original data corrected with bias
   return old_dat - bias


def merge_met_reanalysis(file_upd, file_nrt, file_tpl, file_out=None):

   try:
      h4_upd = SD(file_upd, SDC.READ)
      h4_nrt = SD(file_nrt, SDC.READ)

      # get metadata from updated file
      try:
         start_time  = getattr(h4_upd, 'Start Time')[:9]  # yyyyjjjhh
         data_source = getattr(h4_upd, 'Data Source')
         p = re.compile('NCEP Reanalysis (\d)')
         reanalysis_type = p.search(data_source).group(1)
      except:
         print('Reanalysis file has non-standard Data Source: "'
               + data_source + '"')
         reanalysis_type = 'X'

      # define output filename
      if file_out is None:
         file_out ="N" + start_time + "_MET_NCEPR" + reanalysis_type + "_6h.hdf"

      # start with realtime file if it's full-size
      nx = getattr(h4_nrt, 'Number of Columns')
      ny = getattr(h4_nrt, 'Number of Rows')

      if (nx == NX_OUT) and (ny == NY_OUT):
         print('Calling apply_reanalysis_bias to make file: ', file_out)
         copyfile(file_nrt, file_out)
         h4_out = SD(file_out, SDC.WRITE)

         apply_reanalysis_bias(h4_nrt, h4_upd, h4_out)
         infiles = file_nrt+' '+file_upd
         copy_atts = ('Title', 'Data Source', 'Data Source Desc')

      # otherwise, start with file size appropriate to observation date
      else:
         print('Calling expand_reanalysis_file to make file: ', file_out)
         is_modern = (int(start_time[:4]) > 1995) # expand post-CZCS era only
         if is_modern:
            copyfile(file_tpl, file_out) # start with large dummy file
         else:
            copyfile(file_nrt, file_out) # start with small realtime file
         h4_out = SD(file_out, SDC.WRITE)
         expand_reanalysis_file(h4_nrt, h4_upd, h4_out, expand_array=is_modern)
         infiles = file_upd
         copy_atts = ('Title', 'Data Source', 'Data Source Desc',
                      'Start Time', 'Start Year', 'Start Day', 'Start Millisec',
                      'End Time'  , 'End Year'  , 'End Day'  , 'End Millisec'  )

      # update output file attributes
      setattr(h4_out, 'Product Name', file_out)
      setattr(h4_out, 'Data Center', 'NASA/GSFC Ocean Biology Processing Group')
      setattr(h4_out, 'Mission', 'General')
      setattr(h4_out, 'Satellite Platform', ' ')
      setattr(h4_out, 'Software ID', 'merge_met_reanalysis')
      setattr(h4_out, 'Processing Time',
              datetime.utcnow().strftime('%Y%j%H%M%S%f'))
      setattr(h4_out, 'Input Files', infiles)
      setattr(h4_out, 'Processing Control', ' ')
      for att in copy_atts:   # copy others from reanalysis file
         setattr(h4_out, att, getattr(h4_upd, att))

      print('Success!')

   except Exception as e:
      print(e)
      exit(1)

   finally:
      try:
         h4_upd.end()
         h4_nrt.end()
         h4_out.end()
      except:
         pass


def apply_reanalysis_bias(h4_nrt, h4_upd, h4_out, width=11):

   try:
      for sds in params: # {

         # for wind params, use realtime data
         if (sds == 'z_wind') or (sds == 'm_wind'):
            data_out = h4_nrt.select(sds).get()

         # otherwise, adjust according to reanalysis data
         else:
            data_upd = h4_upd.select(sds).get()
            data_nrt = h4_nrt.select(sds).get()
            enforce_range(data_upd, sds)
            enforce_range(data_nrt, sds)

            # derive bias and apply to realtime data
            data_out = fix_bias(data_nrt, data_upd, width=width)

         # write data to existing file
         enforce_range(data_out, sds)
         h4_out.select(sds).set(data_out)
         # }

   except Exception as e:
      print(e)
      raise


def expand_reanalysis_file(h4_nrt, h4_upd, h4_out, expand_array=False):

   try:
      for sds in params: # {

         # for wind params, use realtime data
         if (sds == 'z_wind') or (sds == 'm_wind'):
            data_upd = h4_nrt.select(sds).get()

         # otherwise, use reanalysis data
         else:
            data_upd = h4_upd.select(sds).get()

         # resize reanalysis data to operational resolution
         if expand_array:
            zoom = (NY_OUT/data_upd.shape[0], NX_OUT/data_upd.shape[1])
            data_out = ndimage.interpolation.zoom(
               data_upd, zoom, order=1, mode='nearest')
         else:
            data_out = data_upd

         # write data to existing file
         enforce_range(data_out, sds)
         h4_out.select(sds).set(data_out)
         # }

   except Exception as e:
      print(e)
      raise


if __name__ == "__main__":
   import argparse
   parser = argparse.ArgumentParser(description='Adjust predicted met files with bias derived from NCEP Reanalysis.')
   parser.add_argument('file_upd', metavar='prer2file', help='predicted file')
   parser.add_argument('file_nrt', metavar='metfile', help='updated file')
   parser.add_argument('file_tpl', metavar='r2_template', help='file template', nargs='?')
   parser.add_argument('file_out', metavar='ofile', help='output file', nargs='?')

   args = parser.parse_args()
   dict_args=vars(args)
   #print(dict_args)

   merge_met_reanalysis(dict_args['file_upd'],
                        dict_args['file_nrt'],
                        dict_args['file_tpl'],
                        dict_args['file_out'])
   exit(0)
