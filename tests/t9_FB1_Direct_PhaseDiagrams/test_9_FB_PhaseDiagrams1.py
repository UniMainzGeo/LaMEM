
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re


def test_a():

  # Test a falling block with phase diagram with build-in direct solver on 2 core, using optimized
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t9_PhaseDiagrams/test_9_FallingBlock_PhaseDiagrams.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't9_PhaseDiagrams/test_9_FallingBlock_PhaseDiagrams.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

  # Create unit test object
  ex1 = pth.pthUnitTest('t9_PhaseDiagrams',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
