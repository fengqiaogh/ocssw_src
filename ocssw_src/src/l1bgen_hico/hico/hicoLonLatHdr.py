'''
Created on Dec 8, 2015

@author: rhealy
'''

code_version = 1.0
code_author  = 'R. Healy (richard.healy@nasa.gov) SAIC'
code_name    = 'proc_hico.py'
code_date    = '12/9/2015'

class hicoLonLatHdr():
    '''
    This is a hico Lon Lat Header class for constructing an
    object to be passed to the writeHicoLonLatHdr function.
    '''


    def __init__(self,**kwargs):
        '''
        Constructor: Create the header class consisting of history
        and various parameters written to the LonLat Hdr file.

        '''
        
        self.variables = kwargs
        
    def set_variable(self,k,v):
        self.variables[k]= v
        
    def get_variable(self,k):
        return self.variables.get(k,None)

    def createHistory(self,qinfo):
    
        import getpass
        from datetime import date
        from datetime import datetime
        import socket
        timenow = datetime.now()

        
        history = {}
        
        history['geometry_code_version'] = code_version
        history['geometry_code_date']    = code_date
        history['geometry_code_author']  = code_author
        history['geometry_code_name']    = code_name
        history['geometry_code_SPosVelQuatCsv'] = qinfo['This file name']
        history['geometry_source_csv_filename'] = qinfo['Source CSV file']
        history['geometry_subset_csv_filename'] = qinfo['Subset CSV file']
        history['geometry_run_date']            = date.today()
        history['geometry_run_time']            = timenow.strftime("%H:%M:%S")
        history['geometry_run_by']              = getpass.getuser()
        history['geometry_run_on']              = socket.gethostname()
        history['pvq_code_name']                = qinfo['Code name']
        history['pvq_code_version']             = qinfo['Code version']
        history['pvq_code_date']                = qinfo['Code date']
        history['pvq_code_author']              = qinfo['Code author']
        history['pvq_code_computer']            = qinfo['Code executed on computer']
        history['pvq_code_username']            = qinfo['Code executed by username']
#         history['pvq_code_IDL_family']          = qinfo['Code run under IDL osfamily']
#         history['pvq_code_IDL_os']              = qinfo['Code run under IDL os']
#         history['pvq_code_IDL_osname']          = qinfo['Code run under IDL osname']
#         history['pvq_code_IDL_version']         = qinfo['Code run under IDL version']
#         history['pvq_code_IDL_start_time']      = qinfo['Code start time']
#         history['pvq_code_IDL_end_time']        = qinfo['Code end time']
        history['geometry_hico_angle_in_degrees_from_stowed'] = qinfo['Theta (degrees from stowed position)']
        history['geometry_sensor_orientation']  = qinfo['ISS orientation']
        history['geometry_ISS_Z_direction']     = self.get_variable('geometry_ISS_Z_direction')
        history['hico_exposure_interval_s']     = qinfo['Exposure interval (frame time)']
        history['hico_trigger_pulse_width_s']   = qinfo['Trigger pulse width (s)']
        history['iss_ODRC_broadcast_time-gps_time_s']         = qinfo['ODRC broadcast time - gps time (s)']
        history['geometry_deltaUT1_s']          = self.get_variable('geometry_deltaUT1_s')
        history['geometry_deltaT_s']            = self.get_variable('geometry_deltaT_s')
        history['geometry_LOD']                 = self.get_variable('geometry_LOD')
        history['geometry_xp_rad']              = self.get_variable('geometry_xp_rad')
        history['geometry_yp_rad']              = self.get_variable('geometry_yp_rad')
        history['geometry_bore_sight_params']   = self.get_variable('geometry_bore_sight_params')
        
        self.history = history
 
    def writeHicoENVIHeaderFile(self):
        
        fhdr = open(self.get_variable('filename'),"w")
        
        fhdr.write('ENVI\n')
        fhdr.write('description = { ' + self.get_variable('description') + ' }\n')
        for v in sorted(self.variables):
            if v not in 'description':
                fhdr.write('{} = {}\n'.format(v,self.variables[v]))        
            
        fhdr.write('\nhistory = begins\n')
        
        for key in sorted(self.history.keys()):
            fhdr.write('{} = {}\n'.format(key,self.history[key]))

        fhdr.write('\nhistory = ends\n')
        
        fhdr.close()
