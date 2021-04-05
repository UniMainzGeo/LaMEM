
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


def test_c():

  # Case with partial melting and a phase transition that depends on the amoujnt of partial melting 
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t16_PhaseTransitions/Plume_PhaseTransitions_Melting.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't16_PhaseTransitions/PhaseTransitions-Melting_p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

    key = re.escape("|eRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

  # Create unit test object
  ex1 = pth.pthUnitTest('t16_PhaseTransitions_Melting_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


def test_d():

  # Tests phase transitions wiuth X/Z and Box coordinates
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t16_PhaseTransitions/Plume_PhaseTransitions_Box_XZ.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't16_PhaseTransitions/PhaseTransitions-XBox-p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)


  # Create unit test object
  ex1 = pth.pthUnitTest('t16_PhaseTransitions_XBox_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


def test_e():

  # Tests phase transition triggered by time
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t16_PhaseTransitions/TimeTransition.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't16_PhaseTransitions/TimeTransition-p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)


  # Create unit test object
  ex1 = pth.pthUnitTest('t16_PhaseTransitions_Time_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
