
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re


def test_a():

  # Test phase transitions with an OPEN top boundary
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
  ex1 = pth.pthUnitTest('t16_PhaseTransitions_OpenTop_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


def test_b():

  # Test phase transitions with an free slip top boundary & pressure shift
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t16_PhaseTransitions/Plume_PhaseTransitions.dat -open_top_bound 0 -act_press_shift 1' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't16_PhaseTransitions/PhaseTransitions-FreeSlip_p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t16_PhaseTransitions_FreeSlip_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
