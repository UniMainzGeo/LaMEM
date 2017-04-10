
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launch as launch
import re


def test_a():

  # Test visco-elasto-plastic localization case on 4 cores, using optimized LaMEM
  ranks = 4
  launch = '../bin/opt/LaMEM -ParamFile ./t4_Loc/localization.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't4_Loc/Loc1_a_MUMPS_VEP_opt-p4.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('Loc1_a_MUMPS_VEP_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

def test_b():

  # Test visco-plastic localization case on 4 cores, using optimized LaMEM
  ranks = 4
  launch = '../bin/opt/LaMEM -ParamFile ./t4_Loc/localization_viscoplastic.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't4_Loc/Loc1_b_MUMPS_VP_opt-p4.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('Loc1_b_MUMPS_VP_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


def test_c():

  # Test visco-elasto-plastic localization case on 1 core, using optimized LaMEM
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t4_Loc/localization.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't4_Loc/Loc1_c_Direct_VEP_opt-p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('Loc1_c_Direct_VEP_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
