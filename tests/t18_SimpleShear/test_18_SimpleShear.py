
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re


def test_xz():

  # Test simple shear BC's with various orientations
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t18_SimpleShear/SS.dat'   # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't18_SimpleShear/SimpleShear_xz-p2.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t18_SimpleShear_xz_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

def test_yz():

  # Test simple shear BC's with various orientations
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t18_SimpleShear/SS.dat  -exz_strain_rates 0 -eyz_strain_rates 1e-15 -eyz_num_periods 1'   # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't18_SimpleShear/SimpleShear_yz-p2.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t18_SimpleShear_yz_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

def test_xy():

  # Test simple shear BC's with various orientations
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t18_SimpleShear/SS.dat  -exz_strain_rates 0 -exy_strain_rates 1e-15 -exy_num_periods 1'   # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't18_SimpleShear/SimpleShear_xy-p2.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t18_SimpleShear_xy_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


def test_xz_yz():

  # Test simple shear BC's with various orientations
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t18_SimpleShear/SS.dat  -exz_strain_rates 1e-15 -eyz_strain_rates 1e-15 -eyz_num_periods 1'   # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't18_SimpleShear/SimpleShear_xz_yz-p2.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t18_SimpleShear_xz_yz_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
