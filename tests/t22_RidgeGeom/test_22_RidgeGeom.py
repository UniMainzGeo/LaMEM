
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re


def test_2D():

  # Test ridge geometry 2D using optimized LaMEM
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t22_RidgeGeom/ridge_geom_2D.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't22_RidgeGeom/RidgeGeom2D.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t22_RidgeGeom2D_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)



def test_3D():

  # Test ridge geometry conditions in 3D using optimized LaMEM
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t22_RidgeGeom/ridge_geom_3D_2cores.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't22_RidgeGeom/RidgeGeom3D_2cores.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t22_RidgeGeom3D_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


def test_oblique():

  # Test oblique ridge geometry conditions in 3D using optimized LaMEM   
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t22_RidgeGeom/ridge_geom_oblique_2cores.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't22_RidgeGeom/RidgeGeom_oblique_2cores.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t22_RidgeGeom_oblique_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
