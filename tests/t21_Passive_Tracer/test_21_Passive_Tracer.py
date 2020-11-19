
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re


def test_a():

  # Test a falling block case with build-in direct solver on 1 core, using optimized
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t21_Passive_Tracer/Passive_tracer_ex2D.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't21_Passive_Tracer/Passive_tracer-2D_p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t21_Passive_Tracer_Always',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

def test_b():

  # Test a falling block case with build-in direct solver on 1 core, using the debug
  ranks = 1

  launch = '../bin/opt/LaMEM -ParamFile ./t21_Passive_Tracer/Passive_tracer_ex2D_Condition.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't21_Passive_Tracer/Passive_tracer-2D_Condition_p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)
    
  #  key = re.escape("Currently active tracers")
   # unittest.compareFloatingPoint(key,10.0)

  # Create unit test object
  ex1 = pth.pthUnitTest('t21_Passive_Tracer_Condition',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

