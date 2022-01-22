import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re

def test_recharge1():

  # Test dike feature using optimized LaMEM
    
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t28_HeatRecharge/FallingBlockHeatReacharge1.dat' # This must be a relative path with respect to runLaMEM_Tests.p
  expected_file = 't28_HeatRecharge/t28_HeatRecharge1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,5e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,2e-5)

  # Create unit test object
  ex1 = pth.pthUnitTest('t28_recharge1',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

def test_recharge2():

  # Test dike feature using optimized LaMEM
    
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t28_HeatRecharge/FallingBlockHeatReacharge2.dat' # This must be a relative path with respect to runLaMEM_Tests.p
  expected_file = 't28_HeatRecharge/t28_HeatRecharge2.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,3e-6)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,3e-5)

  # Create unit test object
  ex1 = pth.pthUnitTest('t28_recharge2',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
