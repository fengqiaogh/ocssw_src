import unittest
from hico.makePosVelQuatCsv import MakePosVelQuatCSV
from hico.HicoL0toL1B import HicoL0toL1b
from hico.auxiliary import CNamespace
import os


class TestMPVQC(unittest.TestCase):
    """
    This test depends on a particular input, bil, file
    coupled with a particular output, csv, file.
    Note that this test will rely on some test files located in the test folder.
    """

    def setup(self):
        testdatadir = os.path.join(os.path.realpath('.'), 'testdata')
        dataBaseName = os.path.commonprefix(os.listdir(testdatadir))
        try:
            pArgs = CNamespace(l0file='%s.bil' % dataBaseName,
                               hdr='%s.rhdr' % dataBaseName, csvfile='%s.csv' % dataBaseName)
            hcPtr = HicoL0toL1b(pArgs)
        except Exception as e:
            print(e)
        finally:
            self.mpvqc = MakePosVelQuatCSV(hcPtr)

    def test_pvqHeader(self):
        pass

    def test_pvqData(self):
        pass


class TestAuxiliary(unittest.TestCase):
    """
    This tests functions in hico/auxiliary.py that MakePosVeluatCSV depends on.
    """

    def test_ConvLst2DT_is_DT(self):
        # start_date_time = ConvLst2DT
        pass

    def test_ConvLst2DT_is_correct(self):
        pass
