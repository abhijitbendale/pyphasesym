########################################################################
# This module is regression test module for python version of phasesym #
# This module compares between matlab outputs and python outputs for   #
# same inputs as a part of regression tests                            #
#                                                                      #
# This module is part of pyphasesym package                            #
########################################################################
""" This is regression test module for python phasesym"""

import glob

import scipy as sp
import scipy.io as so

from numpy.testing import assert_array_almost_equal

from pyphasesym import phasesym

# Get all the mat files for regression test
MAT_FILES = glob.glob('matfiles/*.mat') 

#-------------------------------------------------------------------------------
def compare_py_mat(mat_file):

    """Regression test function

    Compare output of matlab and python for same inputs
    Inputs obtained from stored mat files 
    Output for matlab obtained from stored mat files
    Output for python obtained by invoking phasesym function from 
    mainPhasesym.py module"""

    # Here you could have **kwargs as dicts so that you don't have to type
    # all the variables. You can just reference the variables

    mat_vars = so.loadmat(mat_file)
    p_phasesym, p_orientation = phasesym(mat_vars['image'], 
                                         mat_vars['scale'], 
                                         mat_vars['orient'],
                                         mat_vars['minWaveLength'], 
                                         mat_vars['mult'], 
                                         mat_vars['sigmaOnf'], 
                                         mat_vars['dThetaOnSigma'], 
                                         mat_vars['k'], 
                                         mat_vars['polarity'])

    # have kwargs
    assert_array_almost_equal(p_phasesym, mat_vars['phaseSym'])
    assert_array_almost_equal(p_orientation, mat_vars['orientation'])

#-------------------------------------------------------------------------------
def test_generator():

    """ Generate tests as long mat files exists files 

    mat files stored in MAT_FILES/ directory of the package"""
    
    for mat_file in MAT_FILES:
        yield compare_py_mat, mat_file 
    

        


