'''
Created on Feb 1, 2016

@author: rhealy
'''

code_version = 1.0
code_author  = 'R. Healy (richard.healy@nasa.gov) SAIC'
code_name    = 'interpFields2HicoTimes.py'
code_date    = '2/3/2016'

def qtpow(q,pwr):
    import sys
    import numpy as np

    # NAME:
    #   qtpow
    #
    # Converted to python by R. Healy (2/9/2016) richard.healy@nasa.gov
    #
    # AUTHOR:
    #   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
    #   craigm@lheamail.gsfc.nasa.gov
    #   UPDATED VERSIONs can be found on my WEB PAGE:
    #      http://cow.physics.wisc.edu/~craigm/idl/idl.html
    #
    # PURPOSE:
    #   Raise quaternion Q to the "power" POW
    #
    # MAJOR TOPICS:
    #   Geometry
    #
    # CALLING SEQUENCE:
    #   QNEW = QTPOW(Q, POW)
    #
    # DESCRIPTION:
    #
    #   The function QTPOW raises a quaterion Q to the power P.  The
    #   operation
    #
    #      QNEW = QTPOW(Q, POW)
    #
    #   is equivalent to
    #
    #      QNEW = QTEXP( POW * QTLOG(Q))
    #
    #   which is the same as the definition of raising a real number to
    #   any power (however, QTPOW is faster than using QTLOG and QTEXP).
    #
    #   For integer values of POW, this form of exponentiation is also
    #   directly equivalent to the multiplication of that many Q's
    #   together.

    nq = q.size/4
    npw = pwr.size

    if nq < 1 or npw < 1:
        sys.exit('qtpow: Invalid array size')

    v = q[:,0:3]
    sinth = np.sqrt(np.sum(np.power(v,2),axis=1))
    th = np.arctan2(sinth, q[:,3])
    rat = th*0
    wh = (sinth != 0)
    if wh.size > 0:
        rat[wh] = np.sin(np.multiply(pwr[wh],th[wh]))/sinth[wh]
    q1 = np.zeros(q.shape)
    q1[:,3]   = np.cos(np.multiply(th,pwr))

    q1[:,0:3] = np.multiply(np.tile(rat,(3,1)).T,v)
    return q1

def qtmult(aqt, bqt, inverse1=0, inverse2=0):
    import numpy as np
    import sys

    #This will only work with hico

    sz1 = aqt.shape
    sz2 = bqt.shape
    print('sz1,sz2=',sz1,sz2)
    if sz1[0] < 1 or sz2[0] < 1:
        sys.exit('ERROR A: Q1 and Q2 must be quaternions')
    if sz1[1] != 4 or sz2[1] != 4:
        sys.exit('ERROR B: Q1 and Q2 must be quaternions')

    n1 = aqt.size/4
    n2 = bqt.size/4

    if n1 != n2 and n1 != 1 and n2 != 1:
        sys.exit( 'ERROR: Q1 and Q2 must both have the same number of quaternions')

    #nq = np.max(n1,n2)

    cqt = np.zeros(aqt.shape)

    aqt0 = np.squeeze(aqt[:,0])
    aqt1 = np.squeeze(aqt[:,1])
    aqt2 = np.squeeze(aqt[:,2])
    aqt3 = np.squeeze(aqt[:,3])


    bqt0 = np.squeeze(bqt[:,0])
    bqt1 = np.squeeze(bqt[:,1])
    bqt2 = np.squeeze(bqt[:,2])
    bqt3 = np.squeeze(bqt[:,3])

    if inverse1 > 0:
        aqt0 = -aqt0
        aqt1 = -aqt1
        aqt2 = -aqt2

    if inverse2 > 0:
        bqt0 = -bqt0
        bqt1 = -bqt1
        bqt2 = -bqt2
#     print('aqt1=',aqt1.shape,'aqt2=',aqt2.shape,'aqt3=',aqt3.shape,'aqt0=',aqt0.shape)
#     print('bqt1=',bqt1.shape,'bqt2=',bqt2.shape,'bqt3=',bqt3.shape,'bqt0=',bqt0.shape)
#     print('mult=',np.multiply(aqt0,bqt3).shape)
    cqt[:,0] = np.squeeze([ np.multiply(aqt0,bqt3) + np.multiply(aqt1,bqt2) - np.multiply(aqt2,bqt1) + np.multiply(aqt3,bqt0)] )
    cqt[:,1] = np.squeeze([-np.multiply(aqt0,bqt2) + np.multiply(aqt1,bqt3) + np.multiply(aqt2,bqt0) + np.multiply(aqt3,bqt1)] )
    cqt[:,2] = np.squeeze([ np.multiply(aqt0,bqt1) - np.multiply(aqt1,bqt0) + np.multiply(aqt2,bqt3) + np.multiply(aqt3,bqt2)] )
    cqt[:,3] = np.squeeze([-np.multiply(aqt0,bqt0) - np.multiply(aqt1,bqt1) - np.multiply(aqt2,bqt2) + np.multiply(aqt3,bqt3)] )

    return cqt

def qterp_slerp(t0, q0, t1):
    from value_locate import value_locate
    import numpy as np
    import pandas as pd
    # NAME:
    #   QTERP
    #   Converted to python by R. Healy (richard.healy@nasa.gov) 2/10/16
    #
    # AUTHOR:
    #    Translated to python by Rick Healy (richard.healy@nasa.gov)
    #
    #   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
    #   craigm@lheamail.gsfc.nasa.gov
    #   UPDATED VERSIONs can be found on my WEB PAGE:
    #      http://cow.physics.wisc.edu/~craigm/idl/idl.html
    #
    # PURPOSE:
    #   Smoothly interpolate from a grid of quaternions (spline or slerp)
    #
    # MAJOR TOPICS:
    #   Geometry
    #
    # CALLING SEQUENCE:
    #   QNEW = QTERP(TGRID, QGRID, TNEW, [/SLERP], QDIFF=, [/RESET])
    #
    # DESCRIPTION:
    #
    #  The function QTERP is used to interplate from a set of known unit
    #  quaternions specified on a grid of independent values, to a new set
    #  of independent values.  For example, given a set of quaternions at
    #  specified key times, QTERP can interpolate at any points between
    #  those times.  This has applications for computer animation and
    #  spacecraft attitude control.
    #
    nq = np.int(q0.size/4)

    if nq == 1:
        return np.tile(q0.T,(t1.size,1))

    qdiff = qtmult(q0[0:nq-1,:], q0[1:,:],inverse1=1)
    wh = qdiff[:,3] <0
    qn = qdiff[wh,3]
    if qn.size > 0:
        qdiff[wh,3] = -qdiff[wh,3]

    ii = value_locate(t0, t1)
    hh = np.zeros(len(t1))
    # Maybe there's a better way to do this, but for now this works.
    # The way IDL does this is totally mysterious to me since it
    # subtracts t0(ii) from t1 without reference to ii or the number of ii's
    # and passing hh the same way
    for i in np.arange(0,len(ii)):
        if ii[i]>0 and ii[i]<nq-1:
            hh[i] = (t1[i]-t0[ii[i]])/(t0[ii[i]+1]-t0[ii[i]])
    ii2 = (ii[ii>0])
    ii3 = (ii2[ii2<nq-1])
    return qtmult(q0[ii3,:],qtpow(qdiff[ii3,:],hh[ii3]) )
#
#   I've only implemented the slerp component of qterp
#   since that's all we use. Here's the IDL code for the spline component
#   for later implementation, if desired
#   -rjh 2/10/16
#
#       q1 = (q0(*,0) # t1) * 0
#       for i = 0, 3 do $
#         q1(i,*) = spl_interp(t0, q0(i,*), spl_init(t0, q0(i,*)), t1)
#       tot = sqrt(total(q1^2,1))
#       for i = 0, 3 do $
#         q1(i,*) = q1(i,*) / tot
#       return, q1


class Hico(object):

    def __init__(self,fileName,csvFileName,delta_odrcBTmGPS,iss_orientation,n_pixels=512,delta_texp=0,delta_ticugps=0,delta_tisspvq=0):
        from HicoHeader import HicoHeader
        from astreduc import JDay
        import numpy as np
        import pandas as pd
        from scipy.interpolate import interp1d,UnivariateSpline
        import math
        import sys
        from datetime import datetime

        self.begin_time = datetime.now()
        self.csvFileName = csvFileName
        infile_root = '.'.join(fileName.split('.')[0:-2])
        print('infile_root={}'.format(infile_root))
        self.outname = '{}_pos_vel_quat.csv'.format(infile_root)
        self.anglename = '{}_LonLatViewAngles.bil'.format(infile_root)
        self.n_pixels = n_pixels
        self.delta_odrcBTmGPS = delta_odrcBTmGPS
        # GPS Week 0 began at 00:00:00 UTC 1980 Jan 06
        # sec in day*( (days in 4yr) * 7 sets +
        # yr2008 - first five days of 1980)
        #          + leap seconds since 1980 Jan 6
        self.gps_seconds_2009Jan01_00_00_00 = np.dtype('int32')
        self.gps_seconds_2009Jan01_00_00_00 = 86400*((3 * 365 + 366)*7 +
                                             366 - 5 ) + 15

        self.orientation = iss_orientation
        # width of pulse. PPS records arrival, image not triggered 'til end.
        self.trigger_pulse_width = 0.860e-3


        hdr = HicoHeader(fileName)
        self.header = hdr.header
        self.L0 = hdr.L0
        self.end_time = hdr.L0['end_time'][2]
        self.nls=hdr.L0['n_image']
        self.ffpps=hdr.header['FFpps']
        self.ffpps_sub=hdr.header['FFppsSub']
        self.lfpps=hdr.header['LFpps']
        self.lfpps_sub=hdr.header['LFppsSub']
        self.thetas=hdr.L0['ScenePointingAngle']

        print('estimated start_date=',hdr.L0["start_date"],hdr.L0["start_time"])
        print('estimated end_date=',hdr.L0["end_date"],hdr.L0["end_time"])
        start_t= hdr.L0["start_time"]
        start_t[2]-=15
        self.start_date, self.start_time = self.timeIncrement(hdr.L0["start_date"],start_t)
        end_t= hdr.L0["end_time"]
        end_t[2]+=15
        self.end_date, self.end_time = self.timeIncrement(hdr.L0["end_date"],end_t)

        print('start_date=',self.start_date,self.start_time)
        print('end_date=',self.end_date,self.end_time)
        self.start_struct = self.getTimeStruct(self.start_date,self.start_time)
        self.end_struct = self.getTimeStruct(self.end_date,self.end_time)

        print('Year for jday=',self.start_struct['year'])
        self.Jday_start = JDay(self.start_struct['century']*100+self.start_struct['year'],self.start_struct['month'],self.start_struct['day'],
                          self.start_struct['hour'],self.start_struct['minute'],self.start_struct['second'])
        self.Jday_end = JDay(self.end_struct['century']*100+self.end_struct['year'],self.end_struct['month'],self.end_struct['day'],
                          self.end_struct['hour'],self.end_struct['minute'],self.end_struct['second'])

        if delta_texp <= 0:
            self.exptimes = hdr.header['TriggerCount']* 1.117460459e-6
            self.cexp = '_expDEF'
        else:
            self.exptimes = hdr.header['TriggerCount']* 1.0
            self.cexp = '_exp{}'.format(np.int(delta_texp*1.0e7))

        self.d_ticugps = delta_ticugps
        if delta_ticugps == 0:
            self.cdticu='_dticu0'
        else:
            self.cdticu='_dticu{}'.format(np.int(delta_ticugps*1.0e4))


        self.d_tisspvq = delta_tisspvq
        if delta_tisspvq == 0:
            self.cdtpvq='_dpvq0'
        else:
            self.cdtpvq='_dpvq{}'.format(np.int(delta_tisspvq*1.0e4))

        self.hsdatf = pd.read_csv(csvFileName)

        self.hsdat = self.gethsdatrange(self.hsdatf,self.Jday_start,self.Jday_end)

        #self.hsdat=self.hsdatf.copy()
        self.hsdat=self.hsdatf[self.jtB:self.jtE].set_index(np.arange(0,self.jtE-self.jtB)).copy()
        print('hsdat shape=',self.hsdat.shape)
        self.fixed_hwreg = 0
        if self.hsdat.loc[0,'ICUTIMEHWREGSECOND'] > self.hsdat.loc[len(self.hsdat.loc[:,'ICUTIMEHWREGSECOND'])-1,'ICUTIMEHWREGSECOND']:
            idx = self.hsdat['ICUTIMEHWREGSECOND'] < self.hsdat.loc[0,'ICUTIMEHWREGSECOND']
            self.hsdat.loc[idx,'ICUTIMEHWREGSECOND']+=65536
            self.fixed_hwreg = 1

        self.ffpps_all=self.ffpps+self.ffpps_sub*16.72e-6
        self.lfpps_all=self.lfpps+self.lfpps_sub*16.72e-6

        if self.lfpps_all < self.ffpps_all:
            self.ffpps_all+=65536

        self.hwreg=self.hsdat.loc[:,'ICUTIMEHWREGSECOND'] + self.hsdat.loc[:,'ICUTIMEHWREGSUBSECOND'] * 16.72e-6

        if self.fixed_hwreg == 1:
        # at this point both ffpps_all and lfpps_all are either < 65536 or >= 65536
        # since fixed_hwreg==1, that means we check one and add
            if self.lfpps_all < self.hwreg[0]:
                self.lfpps_all+=65536
            if self.ffpps_all < self.hwreg[0]:
                self.ffpps_all+=65536
        print('lfpps:',self.lfpps_all,self.hwreg[len(self.hwreg)-1])
        if self.lfpps_all < self.hwreg[0]:
            sys.exit('\n================================================\n \
                       interpolate_fields_to_hico_times error: Error    \n \
                       All hsm csv data is AFTER L0 file - check inputs \n \
                       L0 = {}\n CSV = {}                               \n \
                       ================================================ \n \
                      \n'.format(fileName,csvFileName))

        if self.ffpps_all > self.hwreg[len(self.hwreg)-1]:
            sys.exit('\n================================================\n \
                       interpolate_fields_to_hico_times error: Error    \n \
                       All hsm csv data is BEFORE L0 file - check inputs\n \
                       L0 = {}\n CSV = {}                               \n \
                       ================================================ \n \
                      \n'.format(fileName,csvFileName))

        if self.ffpps_all < self.hwreg[0]:
            sys.exit('\n================================================      \n \
                       interpolate_fields_to_hico_times error: Error          \n \
                       hsm csv data starts AFTER L0 starts file - check inputs\n \
                       L0 = {}\n CSV = {}                                     \n \
                       ================================================       \n \
                      \n'.format(fileName,csvFileName))

        if self.lfpps_all > self.hwreg[len(self.hwreg)-1]:
            sys.exit('\n================================================  \n \
                       interpolate_fields_to_hico_times error: Error      \n \
                       hsm csv data ends BEFORE L0 file ends- check inputs\n \
                       L0 = {}\n CSV = {}                                 \n \
                       ================================================   \n \
                      \n'.format(fileName,csvFileName))

        # the last of these is the tick just before the camera starts
        # This is used to get the starting time in the HICO file.
        #self.hwreg_locs_lt=self.hwreg[self.hwreg <= self.ffpps_all]
        #self.just_before=len(self.hwreg_locs_lt)-2
        self.locs_lt = self.hwreg <= self.ffpps_all
        self.idx_lt = self.locs_lt.index[self.locs_lt == True]
        # the first of these is the tick just after the camera stops
        self.hwreg_locs_gt=self.hwreg[self.hwreg >= self.lfpps_all]
        self.just_after=0

        # "t" are the times associated with the ISSPOSITIONS, and, perhaps more
        # basically, with the PPS signals.
        # For this, I need to convert the CC-YY-MM-DD-HH:MM:SS.ssss into "seconds
        # since the start of the file" - NEEDS TO BE WRITTEN - I was going to make
        # another subroutine (well, function of program to IDL) to do this.  "t"
        # has one "time" per lines in the file.
        #
        # really, this is t_icu, the broadcast time that the PPS arrived

        self.t_icu  = self.time_to_seconds(self.hsdat['ICUTIMEISSHOUR'].values,
                                           self.hsdat['ICUTIMEISSMINUTE'].values,
                                           self.hsdat['ICUTIMEISSSECOND'].values,
                                           self.hsdat['ICUTIMEISSSUBSECOND'].values)
        self.t_icugps = self.hsdat.loc[:,'ICUTIMEGPSSECONDS'].values + 1.e-6 * self.hsdat.loc[:,'ICUTIMEISSSUBSECOND'].values


        # need to add subseconds? used for Star tracker
        self.t_issgps = self.hsdat.loc[:,'_ISSGPSTIME'].values + 1.e-6 * self.hsdat.loc[:,'ISSTIMESUBSECOND'].values
        self.t_iss = self.time_to_seconds(self.hsdat['ISSTIMEHOUR'].values,
                                          self.hsdat['ISSTIMEMINUTE'].values,
                                          self.hsdat['ISSTIMESECOND'].values,
                                          self.hsdat['ISSTIMESUBSECOND'].values)
        # Need to take of the Broadcast time - GPS time difference. Even though the
        # BAD times are marked as GPS, they are not held strictly to GPS. So we
        # need to take care of the difference. One of the (now required) input
        # arguments is delta_odrcBTmGPS.
        # Now Bill references the _ISSGPSSECONDS in his email. This seems to be
        # exactly the same as the ICUTIMEGPSSECONDS. Now the notes don't say what
        # the source of the latter is, except when you take the fields that they
        # claim are from the broadcast ancillary data, they add up to the same
        # number of seconds.
        # 2012/12/03 test 0 time offset by commenting

        self.t_issgps-=delta_odrcBTmGPS # This is the sense Bill uses.
        self.t_icugps-=delta_odrcBTmGPS # This is the sense Bill uses.
        self.t_icugps+=self.d_ticugps # 2012/01/17 for testing

        # Now the time of the first line is the time at the location in the CSV
        # file + the difference after since the last update?
        # To get the time start, linearly interpolate ffpps_all into the regime
        # between the PPSCounters (hwreg above) at the known times (t). Voila! You
        # have the start time in seconds since the first time in this (small)
        # HSM_CSV file.

        # hwstart is fraction of interval
        #self.hwstart= (self.ffpps_all-self.hwreg_locs_lt[self.just_before])/(self.hwreg_locs_lt[self.just_before + 1] -
        #    self.hwreg_locs_lt[self.just_before])
        self.hwstart= (self.ffpps_all-self.hwreg[self.idx_lt[-2]])/(self.hwreg[self.idx_lt[-1]] -
            self.hwreg[self.idx_lt[-2]])

        # time_start is the start time of the first line wrt to the start of the
        # "time" file
        self.time_start = self.hwstart*(self.t_icugps[self.idx_lt[-1]] -
                            self.t_icugps[self.idx_lt[-2]]) + self.t_icugps[self.idx_lt[-2]]

        # These are the times we want. The start time is the; all wrt start time in
        # the time file.
        # Addition of 0.860 ms : an email on Aug. 27 by Dan Korwan explains
        # this. Basically, the trigger happens before the image. The recorded time
        # (register) is at the start of the pulse, but the image is not obtained
        # until the end of the pulse, which is 0.860 ms later.
        #       hico_times = time_start + exptimes[i]*indgen(nls[i]) + 0.860d0
        #       hico_times = time_start + exptimes[i]*indgen(nls[i]) + 0.860d-3 ; 2010 Sep 2
        #       hico_times = time_start + exptimes[i]*indgen(nls[i]) + trigger_pulse_width ; 2010 Sep 3
        # Added 1/2 frame time to better center
        self.scanrange = np.arange(self.nls,dtype=float)
        self.testrange = np.multiply(self.scanrange,self.exptimes)
        self.hico_times = np.add(np.add(self.time_start, self.testrange),(self.trigger_pulse_width + 0.5*self.exptimes))
        #self.hico_times = self.time_start + self.exptimes* + self.trigger_pulse_width + 0.5*self.exptimes # 2013Jan17
        self.locs_goodt=(self.hsdat['USGNC_PS_Pointing_Coarse_Time_Tag'] >= (self.hico_times[0]-10)) \
                        & (self.hsdat['USGNC_PS_Pointing_Coarse_Time_Tag'] <= (self.hico_times[-1]+10))
        self.idx_goodt = self.locs_goodt[self.locs_goodt == True].index

        if len(self.idx_goodt) == 0:
            sys.exit('No valid USGNC Coarse times in range of hico_times\n')

        #self.u_usgnc = pd.unique(self.hsdat.loc[self.idx_goodt,'USGNC_PS_Pointing_Coarse_Time_Tag']) #.values.ravel())
        # now only get indexes from unique times
        self.u_usgnc = self.hsdat.loc[self.idx_goodt,'USGNC_PS_Pointing_Coarse_Time_Tag'].duplicated() #.values.ravel())
        self.udx_usgnc = self.u_usgnc.index[self.u_usgnc == False]
        if len(self.udx_usgnc) < 4:
            sys.exit('Less than 4 unique USGNC times during HICO collect\n')

        u_usgnc_coarse=self.hsdat.loc[self.udx_usgnc,'USGNC_PS_Pointing_Coarse_Time_Tag'].values
        u_usgnc_fine=self.hsdat.loc[self.udx_usgnc,'USGNC_PS_PD_Fine_Pointing_Fine_Time_Tag'].values
        u_usgnc_rx=self.hsdat.loc[self.udx_usgnc,'USGNC_PS_Pointing_Inert_Posn_VectorX'].values
        u_usgnc_ry=self.hsdat.loc[self.udx_usgnc,'USGNC_PS_Pointing_Inert_Posn_VectorY'].values
        u_usgnc_rz=self.hsdat.loc[self.udx_usgnc,'USGNC_PS_Pointing_Inert_Posn_VectorZ'].values
        u_usgnc_vx=self.hsdat.loc[self.udx_usgnc,'USGNC_PS_Pointing_Inert_Vel_VectorX'].values
        u_usgnc_vy=self.hsdat.loc[self.udx_usgnc,'USGNC_PS_Pointing_Inert_Vel_VectorY'].values
        u_usgnc_vz=self.hsdat.loc[self.udx_usgnc,'USGNC_PS_Pointing_Inert_Vel_VectorZ'].values
        u_usgnc_q0=np.array(self.hsdat.loc[self.udx_usgnc,'USGNC_PS_Pointing_Current_Inert_Att_Quatrn_0'].values)
        self.u_usgnc_q1=np.array(self.hsdat.loc[self.udx_usgnc,'USGNC_PS_Pointing_Current_Inert_Att_Quatrn_1'].values)
        u_usgnc_q2=np.array(self.hsdat.loc[self.udx_usgnc,'USGNC_PS_Pointing_Current_Inert_Att_Quatrn_2'].values)
        u_usgnc_q3=np.array(self.hsdat.loc[self.udx_usgnc,'USGNC_PS_Pointing_Current_Inert_Att_Quatrn_3'].values)

        self.t_issposvelquat=u_usgnc_coarse + u_usgnc_fine/256.
        self.t_issposvelquat+=self.d_tisspvq
        f = UnivariateSpline(self.t_issposvelquat, u_usgnc_rx)
        self.ISSPOSITIONX =  f(self.hico_times)
        f = UnivariateSpline(self.t_issposvelquat, u_usgnc_ry)
        self.ISSPOSITIONY =  f(self.hico_times)
        f = UnivariateSpline(self.t_issposvelquat, u_usgnc_rz)
        self.ISSPOSITIONZ =  f(self.hico_times)
        f = UnivariateSpline(self.t_issposvelquat, u_usgnc_vx)
        self.ISSVELOCITYX = f(self.hico_times)
        f = UnivariateSpline(self.t_issposvelquat, u_usgnc_vy)
        self.ISSVELOCITYY = f(self.hico_times)
        f = UnivariateSpline(self.t_issposvelquat, u_usgnc_vz)
        self.ISSVELOCITYZ = f(self.hico_times)

        self.issqt = np.matrix([self.u_usgnc_q1.T,u_usgnc_q2.T,u_usgnc_q3.T,u_usgnc_q0.T]).T
        self.issqt_hicotimes = qterp_slerp(self.t_issposvelquat,self.issqt,self.hico_times)

        for k in np.arange(0,self.issqt_hicotimes.shape[1]):
            if any(np.isnan(self.issqt_hicotimes[:,k])):
                sys.exit('returning: NANs detected in interpolated quaternions')

        print('HSTCLOCKTIME0=',self.hsdat['HSTCLOCKTIME0'])
        self.hstclocktime = self.hsdat['HSTCLOCKTIME0']*0.050 + self.hsdat['HSTCLOCKTIME1'] * 62.5e-9
        if all(t==0 for t in self.hstclocktime) or all(t==0 for t in self.hsdat['HSTCLOCKTIME0']):
            # no star tracker telemetry?
            nht=size(self.hico_times)
            hstqt_status=np.array([2]*nht,dtype=int)
            hstqt_hicotimes=np.zeros((4,nht),dtype='f8')
        else:
            # t is the elapsed times wrt to the start time in the
            # file, at the PPS signal.

            self.hstattitudetime=self.hsdat['HSTATTITUDETIME0']*0.050 + self.hsdat['HSTATTITUDETIME1'] * 62.5e-9

            self.hstqt=np.array([self.hsdat['HSTATTITUDEQUATX'],self.hsdat['HSTATTITUDEQUATY'],
                        self.hsdat['HSTATTITUDEQUATZ'],self.hsdat['HSTATTITUDEQUATS']]).T
            print('shapes=',self.hstattitudetime.shape,self.hstqt.shape,self.hstclocktime.shape)

            # The interp below takes the "attitude" times to the "clock times",
            # i.e. interps the quaternions at their measured times to the times of the
            # PPS.

            self.hstqt_pps = qterp_slerp(self.hstattitudetime,self.hstqt,self.hstclocktime)
            if any(math.isnan(t) for t in self.hstqt_pps[:,3]):
                bad=math.isnan(self.hstqt_pps)
                self.hstqt_pps[bad]=0.0
                print('WARNING: NANs detected in unused star tracker quaternions, part a')
                print('reference_file = ',csvFileName)

            # the appropriate times are the ISS times according to Dan
            # Quaternion interpolation from the PPS times to the hico_times.

            self.hstqt_hicotimes=qterp_slerp(self.t_issgps,self.hstqt_pps,self.hico_times)
            f = UnivariateSpline(self.hstattitudetime,self.hsdat['HSTATTITUDESTATUSMODE'])
            ytemp = f(self.hstclocktime)
            f = UnivariateSpline(self.t_issgps,ytemp)
            self.hstqt_status = (f(self.hico_times)).astype(int)

            print('self.hstqt_status=',self.hstqt_status.shape)
            print('hstqt_hicotimes shape=',self.hstqt_hicotimes.shape)
            print('hstqt_hicotimes =',self.hstqt_hicotimes)
            #self.hstqt_hicotimes[0,0] = float('NaN')
            for k in np.arange(0,self.hstqt_hicotimes.shape[1]):
                if any(np.isnan(self.hstqt_hicotimes[:,k])):
                        bad=np.isnan(self.hstqt_hicotimes[:,k])
                        self.hstqt_hicotimes[bad,k]=0.0
                        self.hstqt_status[bad]=2 # "test" wherever there is a NAN
                        print('WARNING: NANs detected in unused star tracker quaternions, part b')
                        print('reference_file = ',csvFileName)

            self.finish_time = datetime.now()

    def gethsdatrange(self,hsdatf,Jday_start,Jday_end):
        from astreduc import JDay
        import numpy as np
        import sys

#         year = 100*hsdatf.loc[0,'ISSTIMECENTURY']+hsdatf.loc[0,'ISSTIMEYEAR']
#         month = hsdatf.loc[0,'ISSTIMEMONTH']
#         day = hsdatf.loc[0,'ISSTIMEDAY']
#         hour = hsdatf.loc[0,'ISSTIMEHOUR']
#         minute = hsdatf.loc[0,'ISSTIMEMINUTE']
#         second = hsdatf.loc[0,'ISSTIMESECOND']

        jt0 = 0
        Jday_now = JDay(100*hsdatf.loc[0,'ISSTIMECENTURY']+hsdatf.loc[0,'ISSTIMEYEAR'],
                        hsdatf.loc[0,'ISSTIMEMONTH'],hsdatf.loc[0,'ISSTIMEDAY'],
                        hsdatf.loc[0,'ISSTIMEHOUR'],hsdatf.loc[0,'ISSTIMEMINUTE'],
                        hsdatf.loc[0,'ISSTIMESECOND'])
        jtend = len(hsdatf)
        while (jt0 < jtend) and (Jday_now < Jday_start):
            jt0+=1
            Jday_now = JDay(100*hsdatf.loc[jt0,'ISSTIMECENTURY']+hsdatf.loc[jt0,'ISSTIMEYEAR'],
                        hsdatf.loc[jt0,'ISSTIMEMONTH'],hsdatf.loc[jt0,'ISSTIMEDAY'],
                        hsdatf.loc[jt0,'ISSTIMEHOUR'],hsdatf.loc[jt0,'ISSTIMEMINUTE'],
                        hsdatf.loc[jt0,'ISSTIMESECOND'])
        if jt0 >= jtend:
            sys.exit('CSV file out of time range')

        jtE = jt0
        while (jtE < jtend) and (Jday_now <= Jday_end):
            jtE+=1
            Jday_now = JDay(100*hsdatf.loc[jtE,'ISSTIMECENTURY']+hsdatf.loc[jtE,'ISSTIMEYEAR'],
                        hsdatf.loc[jtE,'ISSTIMEMONTH'],hsdatf.loc[jtE,'ISSTIMEDAY'],
                        hsdatf.loc[jtE,'ISSTIMEHOUR'],hsdatf.loc[jtE,'ISSTIMEMINUTE'],
                        hsdatf.loc[jtE,'ISSTIMESECOND'])
            print('Hr/Min/Sec=',hsdatf.loc[jtE,'ISSTIMEHOUR'],hsdatf.loc[jtE,'ISSTIMEMINUTE'],
                        hsdatf.loc[jtE,'ISSTIMESECOND'])
        if jtE > jtend:
            sys.exit('CSV file out of time range')


        self.jtB = jt0
        self.jtE = jtE+2

        return hsdatf[self.jtB:self.jtE].set_index(np.arange(0,self.jtE-self.jtB)).copy()

    def createCSVInfo(self):

        import getpass
        import socket,os,platform


        csvinfo = {}
        csvinfo['This file name']=os.path.basename(self.outname)
        csvinfo['Source CSV file']=self.csvFileName
        csvinfo['Subset CSV file']='None'
        csvinfo['Expected Lon/Lat/View angle filename']=os.path.basename(self.anglename)
        csvinfo['Epoch']='2009 Jan 01 00:00:00 UTC'
        csvinfo['Requested number of pixels']='{0:4d}'.format(self.n_pixels)
        csvinfo['Distance Unit']='Feet'
        csvinfo['Central Body']='Earth'
        csvinfo['CoordinateSystem']='J2000'
        csvinfo['Theta (degrees from stowed position)']='{}'.format(self.thetas)

        csvinfo['Code name'] = code_name
        csvinfo['Code version'] = code_version
        csvinfo['Code date']    = code_date
        csvinfo['Code author']  = code_author
        csvinfo['Code executed on computer']    = socket.gethostname()
        csvinfo['Code executed by username']    = getpass.getuser()
        csvinfo['Code run under Python osfamily']    = os.name
        csvinfo['Code run under Python os']    = platform.release()
        csvinfo['Code start time']            = self.begin_time.strftime("%H:%M:%S")
        csvinfo['Code end time']            = self.finish_time.strftime("%H:%M:%S")
        csvinfo['Exposure interval (frame time)'] = '{}'.format(self.exptimes)
        csvinfo['ISS orientation']='{}'.format(self.orientation)
        csvinfo['Trigger pulse width (s)']='{}'.format(self.trigger_pulse_width)
        csvinfo['ODRC broadcast time - gps time (s)']='{}'.format(self.delta_odrcBTmGPS)
        csvinfo['delta_ticugps (s)']='{}'.format(self.d_ticugps)
        csvinfo['delta_tisspvq (s)']='{}'.format(self.d_tisspvq)

        self.csvinfo = csvinfo

        return

    def writeCSVfile(self):
        import numpy as np

        fcsv = open(self.outname,"w")

        fcsv.write('\n')

        for key in sorted(self.csvinfo.keys()):
            fcsv.write('{}, {}\n'.format(key,self.csvinfo[key]))

        fcsv.write('\n')
        fcsv.write('SecondsSinceEpoch, ISSPOSX, ISSPOSY, ISSPOSZ, ISSVELX, ISSVELY, ISSVELZ, ISSQX, ISSQY, ISSQZ, ISSQS, STQX, STQY, STQZ, STQS, HST_ATTITUDE_STATUS\n')

        ofmt='{:20.6f},'*7 +'{:20.8f},'*8 +'{:03d}\n'
        for iqq in np.arange(0,len(self.hico_times)):
            fcsv.write(ofmt.format(self.hico_times[iqq]-self.gps_seconds_2009Jan01_00_00_00,
                                   self.ISSPOSITIONX[iqq],
                                   self.ISSPOSITIONY[iqq],
                                   self.ISSPOSITIONZ[iqq],
                                   self.ISSVELOCITYX[iqq],
                                   self.ISSVELOCITYY[iqq],
                                   self.ISSVELOCITYZ[iqq],
                                   self.issqt_hicotimes[iqq][0],
                                   self.issqt_hicotimes[iqq][1],
                                   self.issqt_hicotimes[iqq][2],
                                   self.issqt_hicotimes[iqq][3],
                                   self.hstqt_hicotimes[iqq][0],
                                   self.hstqt_hicotimes[iqq][1],
                                   self.hstqt_hicotimes[iqq][2],
                                   self.hstqt_hicotimes[iqq][3],
                                   self.hstqt_status[iqq]))


        fcsv.close()


    @staticmethod
    def time_to_seconds(hh,nn,ss,subsec):


        # This is for HICO, CC will be constant.
        # This unwritten subroutine is called to get times in seconds for the
        # HSM_CSV file. At most there is one day turnover - that is the ONLY
        # important turnover for us. This is where HH goes from 23 to 00 since the
        # observations are only ~30 seconds long. Don't care about cc,yy,mm,dd

        # hh,nn,ss,subsec are input arrays of the same length.
        # elapsed is the output array of the same length.


        hhturn = hh < hh[0]
        count = len(hhturn)
        if count > 0:
            hh[hhturn]+=24

        secs_in_day = subsec*1.e-6 + ss + 60.*(nn + hh*60.)
        return secs_in_day-secs_in_day[0]


    @staticmethod
    def time_to_seconds(hh,nn,ss,subsec):


        # This is for HICO, CC will be constant.
        # This unwritten subroutine is called to get times in seconds for the
        # HSM_CSV file. At most there is one day turnover - that is the ONLY
        # important turnover for us. This is where HH goes from 23 to 00 since the
        # observations are only ~30 seconds long. Don't care about cc,yy,mm,dd

        # hh,nn,ss,subsec are input arrays of the same length.
        # elapsed is the output array of the same length.


        hhturn = hh < hh[0]
        count = len(hhturn)
        if count > 0:
            hh[hhturn]+=24

        secs_in_day = subsec*1.e-6 + ss + 60.*(nn + hh*60.)
        return secs_in_day-secs_in_day[0]


    @staticmethod
    def getTimeStruct(date,time):
        import numpy as np

        return {"century":np.int(date[0]/100),"year":(date[0] % 100),"month":date[1],
                       "day":date[2],"hour":np.int(time[0]),
                       "minute":np.int(time[1]),"second":np.int(time[2])}

    @staticmethod
    def timeIncrement(date,time):
        from astreduc import isLeapYear
        import numpy as np

        year  = date[0]
        month = date[1]
        day   = date[2]
        hour  = time[0]
        minute= time[1]
        second= time[2]

        if second >= 60:
            mmm = np.floor(second/60)
            second = np.mod(second,60)
            minute += mmm

        while second < 0:
            second += 60
            minute -= 1

        if minute >= 60:
            hhh = np.floor(minute/60)
            minute = np.mod(minute,60)
            hour+=hhh

        while minute < 0:
            minute += 60
            hour -= 1

        if hour >= 24:
            ddd = np.floor(minute/24)
            hour = np.mod(minute,24)
            day+=ddd

        while hour < 0:
            hour += 24
            day -= 1

        while month > 12:
            month -= 12
            year += 1

        while month < 1:
            month += 12
            year -= 1

        last=np.array([31,31,28+isLeapYear(year),31,30,31,30,31,31,30,31,30,31,31])

        while day > last[month-1]:
            day-=last[month-1]
            month+=1
            if month > 12:
                month-=12
                year+=1
                last=np.array([31,31,28+isLeapYear(year),31,30,31,30,31,31,30,31,30,31,31])

        while day < 1:
            day+=last[month-1]
            month+=1
            if month < 1:
                month+=12
                year-=1
                last=np.array([31,31,28+isLeapYear(year),31,30,31,30,31,31,30,31,30,31,31])

        while month > 12:
            month-=12
            year+=1

        while month < 1:
            month+=12
            year-=1

        return np.array([year,month,day]),np.array([hour,minute,second])


def writeHicoTimes(hico):
    hico.createCSVInfo()
    print(newH.issqt_hicotimes.shape,newH.hstqt_hicotimes.shape)
#    print('a={:20.8f},b={:20.8f},c={:20.8f},d={:20.8f},e={:02d}'.format(hico.hstqt_hicotimes[0][0],hico.hstqt_hicotimes[0][1],hico.hstqt_hicotimes[0][2],hico.hstqt_hicotimes[0][3],hico.hstqt_status[0]))
    hico.writeCSVfile()

    return

def get_odrc_time_offset(filename):
    import sys
    import numpy as np
    for theLine in open(filename,'r') :
        fields = theLine.split('=')
        if len(fields) > 1 and  'odrc_time_offset' in fields[0]:
            return np.float(fields[1])

    sys.exit('odrc_time_offset not in file {}'.format(filename))

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''\
      Generate a HICO geometry file.

            ''', add_help=True)
    parser.add_argument('-ifile', nargs=1, type=str, help=' iss*hico.bil input file (Must be a BIL file) ')
    parser.add_argument('-csvfile', nargs=1, type=str, help=' cvs file ')
    parser.add_argument('-hdr', nargs=1, type=str, help=' header file ')
    parser.add_argument('-orient', default=('-XVV'),nargs=1, type=str, help=' iss orientation file ')
    args = parser.parse_args('-ifile /home/rhealy/src/python/hico/l02l1b/iss.2013067.0308.063527.L0.12933.20130308205743.hico.bil -csvfile /home/rhealy/src/python/hico/l02l1b/iss.2013067.0308.063527.L0.12933.20130308205743.hico.csv -hdr /home/rhealy/src/python/hico/l02l1b/iss.2013067.0308.063527.L0.12933.20130308205743.hico.hdr'.split())

    ifile = ''.join(args.ifile)
    csvfile = ''.join(args.csvfile)
    hdrfile     = ''.join(args.hdr)
    orient = ''.join(args.orient)

    odrc_time_offset = get_odrc_time_offset(hdrfile)

    newH = Hico(ifile,csvfile,odrc_time_offset,orient)

    writeHicoTimes(newH)

#     print(newH.issqt_hicotimes[:,2])
#     print(newH.start_date)
#     print(newH.start_time)
#     print(newH.start_struct)
#     print(newH.end_struct)
#     print(newH.L0["start_time"])
#
#     print('Julday=',newH.Jday_start)
#     print(newH.cexp)
#     print(newH.ffpps_all)
#     print(newH.t_icu)
#     print(newH.hsdat['ICUTIMEISSSECOND'].values)
#     print(newH.ffpps_all)
#     print(newH.hwstart)
#     idx = newH.locs_lt.index[newH.locs_lt == True]
#     print(idx[-2])
#     print(newH.time_start)
#     print(newH.scanrange)
#     print(newH.testrange)
#     print(newH.hico_times)
#     #print(newH.idx_goodt.values)
#     print(newH.udx_usgnc.values)
#     print(newH.ISSPOSITIONX)
#     print(newH.ISSPOSITIONX)
#     print('u_usgnc_q1=',newH.u_usgnc_q1)
#     print('u_usgnc_q1=',len(newH.u_usgnc_q1))
#     print("shape issqt=",newH.issqt.shape)
#     print(newH.issqt[0,1])
#     shp=newH.issqt.shape
#     aout=np.zeros((2,4,47,5,188))
#     print("shape=",shp,aout.shape)
#     print(aout[0].shape)
#     print(len(newH.issqt))
#     print('hico times=',newH.issqt_hicotimes)
#     print('hstclocktime=',newH.hstclocktime)
#     s=newH.hsdat[:10].copy()
#     print(s)
#     print(newH.jtB,newH.jtE)
