
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launch as launch
import re


def test_a():
  
  # Test a falling block case on 1 core, using optimized LaMEM
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t4_Loc/localization_mat_free.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't4_Loc/Loc1_4cores.expected'
  
  def comparefunc(unittest):
 
    key = 'Div_min'
    unittest.compareFloatingPoint(key,1e-7)
  
    key = 'Div_max'
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('unit_Loc1_a',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')
  
  return(ex1)
