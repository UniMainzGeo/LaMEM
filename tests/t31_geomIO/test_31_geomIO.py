import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re

def test_geomIO_1():

  # Test if geomIO polygons are read in correctly
    
  ranks = 4
  launch = '../bin/opt/LaMEM -ParamFile ./t31_geomIO/t31_geomIO_1.dat' # This must be a relative path with respect to runLaMEM_Tests.p
  expected_file = 't31_geomIO/t31_geomIO_1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,2e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,5e-7)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,5e-7)

  # Create unit test object
  ex1 = pth.pthUnitTest('t31_geomIO_1',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

def test_geomIO_2():

  # Test if hollow geomIO polygons are read in correctly
    
  ranks = 4
  launch = '../bin/opt/LaMEM -ParamFile ./t31_geomIO/t31_geomIO_2.dat' # This must be a relative path with respect to runLaMEM_Tests.p
  expected_file = 't31_geomIO/t31_geomIO_2.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,2e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,5e-7)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,5e-7)

  # Create unit test object
  ex1 = pth.pthUnitTest('t31_geomIO_2',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
