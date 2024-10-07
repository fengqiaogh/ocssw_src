import numpy as np
import sys
import pandas as pd
from datetime import datetime as DT


def qtpow(q, pwr):
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
    #   The function QTPOW raises a quaterion Q to the power P.  The operation
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

    v = q[:, 0:3]
    sinth = np.sqrt(np.sum(np.power(v, 2), axis=1))
    th = np.arctan2(sinth, q[:, 3])
    rat = th*0
    wh = (sinth != 0)
    if wh.size > 0:
        rat[wh] = np.sin(np.multiply(pwr[wh], th[wh]))/sinth[wh]
    q1 = np.zeros(q.shape)
    q1[:, 3] = np.cos(np.multiply(th, pwr))
    q1[:, 0:3] = np.multiply(np.tile(rat, (3, 1)).T, v)
    return q1


def qtmult(aqt, bqt, inverse1=0, inverse2=0):
    '''This will only work with hico'''

    sz1 = aqt.shape
    sz2 = bqt.shape
    print('sz1,sz2=', sz1, sz2)
    if sz1[0] < 1 or sz2[0] < 1:
        sys.exit('ERROR A: Q1 and Q2 must be quaternions')
    if sz1[1] != 4 or sz2[1] != 4:
        sys.exit('ERROR B: Q1 and Q2 must be quaternions')

    n1 = aqt.size/4
    n2 = bqt.size/4

    if n1 != n2 and n1 != 1 and n2 != 1:
        sys.exit('ERROR: Q1 and Q2 must both have the same number of quaternions')

    # nq = np.max(n1,n2)
    cqt = np.zeros(aqt.shape)
    aqt0 = np.squeeze(aqt[:, 0])
    aqt1 = np.squeeze(aqt[:, 1])
    aqt2 = np.squeeze(aqt[:, 2])
    aqt3 = np.squeeze(aqt[:, 3])

    bqt0 = np.squeeze(bqt[:, 0])
    bqt1 = np.squeeze(bqt[:, 1])
    bqt2 = np.squeeze(bqt[:, 2])
    bqt3 = np.squeeze(bqt[:, 3])

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
    cqt[:, 0] = np.squeeze([np.multiply(aqt0, bqt3) + np.multiply(aqt1, bqt2) -
                            np.multiply(aqt2, bqt1) + np.multiply(aqt3, bqt0)])
    cqt[:, 1] = np.squeeze([-np.multiply(aqt0, bqt2) + np.multiply(aqt1, bqt3) +
                            np.multiply(aqt2, bqt0) + np.multiply(aqt3, bqt1)])
    cqt[:, 2] = np.squeeze([np.multiply(aqt0, bqt1) - np.multiply(aqt1, bqt0) +
                            np.multiply(aqt2, bqt3) + np.multiply(aqt3, bqt2)])
    cqt[:, 3] = np.squeeze([-np.multiply(aqt0, bqt0) - np.multiply(aqt1, bqt1) -
                            np.multiply(aqt2, bqt2) + np.multiply(aqt3, bqt3)])
    return cqt


def qterp_slerp(t0, q0, t1):
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
    nq = np.int(q0.size / 4)
    if nq == 1:
        return np.tile(q0.T, (t1.size, 1))

    qdiff = qtmult(q0[0:nq-1, :], q0[1:, :], inverse1=1)
    wh = qdiff[:, 3] < 0
    qn = qdiff[wh, 3]
    if qn.size > 0:
        qdiff[wh, 3] = -qdiff[wh, 3]
    # ii = value_locate(t0, t1)
    ii = ValueLocate(t0, t1)
    hh = np.zeros(len(t1))
    # Maybe there's a better way to do this, but for now this works.
    # The way IDL does this is totally mysterious to me since it
    # subtracts t0(ii) from t1 without reference to ii or the number of ii's
    # and passing hh the same way
    for i in np.arange(0, len(ii)):
        if ii[i] > 0 and ii[i] < nq-1:
            hh[i] = (t1[i] - t0[ii[i]]) / (t0[ii[i] + 1] - t0[ii[i]])
    ii2 = (ii[ii > 0])
    ii3 = (ii2[ii2 < nq-1])
    return qtmult(q0[ii3, :], qtpow(qdiff[ii3, :], hh[ii3]))


def ValueLocate(vec, vals):
    '''
    Equivalent to IDL's value_locate routine. Excerpt from exelis follows
    "The VALUE_LOCATE function finds the intervals within a given monotonic
    vector that brackets a given set of one or more search values."
    Assumes vec and val are 1D arrays
    '''
    return np.amax([np.where(v >= vec, np.arange(len(vec)), -1)
                    for v in vals], axis=1)


def ConvLst2DT(date, time):
    """
    Takes a date and time strings and returns a datetime object"""
    ds = ','.join([str(int(x)) for x in date + time])
    try:
        return DT.strptime(ds, '%Y,%m,%d,%H,%M,%S')
    except ValueError:
        print('Conversion error. Check validity of date/time entries:')
        keys = ['yr', 'mo', 'day', 'hr', 'min', 'sec']
        for key, entry in zip(keys, ds.split(',')):
            print('%s: %s' %(key, entry))



def GetOdrcTimeOffset(refHdrFile):
    with open(refHdrFile) as f:
        lines = f.read()
        fline = lines.find('odrc_time_offset')
    if fline >= 0:
        return float(lines[fline:].split('\n')[0].split('=')[1].strip())
    return None


class CNamespace():
    '''
    Class to replace command line argument parser for IPython calls.
    Example Usage: args=Namespace(ifile='',opath='',prsil='',prnoi='')
    **kwargs can be anything needed
    '''

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        return None
