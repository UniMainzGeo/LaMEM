
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re


def test_a():

  # Test visco-elasto-plastic localization case on 4 cores, using optimized LaMEM
  ranks = 4
  launch = '../bin/opt/LaMEM -ParamFile ./t5_Perm/Permea.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't5_Perm/Permeability_direct_opt-p4.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t5_Permeability_Direct_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
