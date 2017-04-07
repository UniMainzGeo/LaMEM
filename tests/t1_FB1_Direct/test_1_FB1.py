
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launch as launch
import re


def test_a():

  # Test a falling block case with build-in direct solver on 1 core, using optimized
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t1_FB1_Direct/FallingBlock_mono_PenaltyDirect.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't1_FB1_Direct/FB1_Direct.expected'

  def comparefunc(unittest):

    key = '|Div|_inf'
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('unit_FB1_a_Direct_Core1_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

def test_b():

  # Test a falling block case with build-in direct solver on 1 core, using the debug
  ranks = 1
  launch = '../bin/deb/LaMEM -ParamFile ./t1_FB1_Direct/FallingBlock_mono_PenaltyDirect.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't1_FB1_Direct/FB1_Direct.expected'

  def comparefunc(unittest):

    key = '|Div|_inf'
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('unit_FB1_b_Direct_Core1_deb',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

def test_c():

  # Test a falling block case with direct solver on 2 core, using the optimized version and MUMPS
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t1_FB1_Direct/FallingBlock_mono_PenaltyDirect.dat -jp_pc_factor_mat_solver_package mumps' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't1_FB1_Direct/FB1_Direct_2cores.expected'

  def comparefunc(unittest):

    key = '|Div|_inf'
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('unit_FB1_c_Direct_Core2_MUMPS',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

def test_d():

  # Test a falling block case with direct solver on 4 cores, using the optimized version and PaStiX
  ranks = 4
  launch = '../bin/opt/LaMEM -ParamFile ./t1_FB1_Direct/FallingBlock_mono_PenaltyDirect.dat -jp_pc_factor_mat_solver_package pastix' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't1_FB1_Direct/FB1_Direct_4cores_PaStiX.expected'

  def comparefunc(unittest):

    key = '|Div|_inf'
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('unit_FB1_d_Direct_Core2_PaStiX',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
