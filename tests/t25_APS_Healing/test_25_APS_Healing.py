
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re


def test_2D():

  # Test APS healing in 2D using optimized LaMEM
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t25_APS_Healing/APS_Healing2D.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't25_APS_Healing/APS_Healing2D.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t25_APS_Healing2D_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

def test_2cores():    

  # Test APS healing in 2D using optimized LaMEM, fails with 2 cores, works with 1 core
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t25_APS_Healing/APS_Healing2cores.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't25_APS_Healing/APS_Healing2cores.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t25_APS_Healing2cores_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
