
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re


def test_a():

  # Test inflow/outflow conditions using optimized LaMEM
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t16_PhaseTransitions/Plume_PhaseTransitions.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't16_PhaseTransitions/PhaseTransitions-p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t16_PhaseTransitions_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
