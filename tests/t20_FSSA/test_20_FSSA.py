# This tests the free surface stabilization algorithm, using the same setup as in the
# original paper of Kaus et al. (2010)

# Note that the tests below are only for a few timesteps & lower resolution. 
# The setup, as committed, 

import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re


def test_1():

  # Test simple shear BC's with various orientations
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t20_FSSA/RTI_FSSA.dat -nstep_max 20 -nel_x 50 -nel_z 100'   # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't20_FSSA/RTI_FSSA_1-p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)
    
    key = re.escape("cloned")
    unittest.compareInteger(key,10)
    
  # Create unit test object
  ex1 = pth.pthUnitTest('t20_FSSA_1_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
