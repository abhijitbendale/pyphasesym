"""This module contains unittests for pyphasesym package"""

import unittest
import Image

import numpy as np

from pyphasesym import *

class phasesym_arr_dimension_tests(unittest.TestCase):
    """class for sanity check for array sizes from various functions in
    pyphasesym"""

    imagename = 'data/cameraman.tif'
    image = Image.open(imagename)
    img_array = np.asarray(image)

    def test_phasesym_from_array(self):
        """test for sanity check for dimensions of phaseSym and orientation
        arrays returned by phasesym function"""

        phase_sym, orientation = phasesym_from_array(self.img_array,
                                                    DEFAULT_NSCALE,
                                                    DEFAULT_ORIENT,
                                                    DEFAULT_MWAVELENGTH,
                                                    DEFAULT_MULT,
                                                    DEFAULT_SIGMAONF,
                                                    DEFAULT_DTHETASEGMA,
                                                    DEFAULT_NSTD,
                                                    DEFAULT_POLARITY)
        self.assertEqual(self.img_array.shape, phase_sym.shape)
        self.assertEqual(self.img_array.shape, orientation.shape)

    def test_get_low_pass_filter(self):
        """test to check dimensions returned from filter_initialization function
        of phasesym package"""

        rows, cols = self.img_array.shape
        sintheta, costheta, radius = filter_intialization(rows, cols)
        self.assertEqual(self.img_array.shape, sintheta.shape)
        self.assertEqual(self.img_array.shape, costheta.shape)
        self.assertEqual(self.img_array.shape, radius.shape)

    def test_lp_filter(self):
        """test to check low pass filter dimension from pyphasesym 
        package"""

        rows, cols = self.img_array.shape
        lp_filter = get_low_pass_filter(rows, cols, 0.4, 10.)
        self.assertEqual(self.img_array.shape, lp_filter.shape)

def suite():
    """Create test suites for tests for pyphasesym"""

    test_suite = unittest.makeSuite(phasesym_arr_dimension_tests, 'test')
    return test_suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
