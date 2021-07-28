import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re

def test_2fields_dike():

  # Test dike feature using optimized LaMEM
    
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t27_T-dep_Conductivity/t27_TDep_NuK_Conductivity.dat' # This must be a relative path with respect to runLaMEM_Tests.p
  expected_file = 't27_T-dep_Conductivity/t27_TDep_NuK_Conductivity.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t27_TdepCond',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
