
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re


def test_a():

  # Test inflow/outflow conditions using optimized LaMEM
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t19_CompensatedInflow/CompensatedInflow_test_2D.dat -nstep_max 10' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't19_CompensatedInflow/CompensatedInflow-p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t19_CompensatedInflow',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

def test_b():

  # Test inflow/outflow conditions using optimized LaMEM
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t19_CompensatedInflow/CompensatedInflow_test_3D.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't19_CompensatedInflow/CompensatedInflow3D-p2.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t19_CompensatedInflow3D',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
