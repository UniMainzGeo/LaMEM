
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re

def test_M1_2D():

  # Test dike feature using optimized LaMEM                                                                                                                                  
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t26_Dike/dike_M1_2D.dat' # This must be a relative path with respect to runLaMEM_Tests.py                                      
  expected_file = 't26_Dike/dike_M1_2D.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object                                                                                                                                                  
  ex1 = pth.pthUnitTest('t26_Dike_opt1',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


def test_M075_2D_2cores():

  # Test dike feature using optimized LaMEM with 2 cores
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t26_Dike/dike_M075_2D_2cores.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't26_Dike/dike_M075_2D_2cores.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t26_Dike_opt2',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


def test_M05_2D():

  # Test dike in 2D  using debugging version of LaMEM
  ranks = 1
  launch = '../bin/deb/LaMEM -ParamFile ./t26_Dike/dike_M05_2D.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't26_Dike/dike_M05_2D.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t26_Dike_deb',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


def test_2cores_2dikes():

# Test dike feature using optimized LaMEM
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t26_Dike/move2Dikes.dat' # This must be a relative path with respect to runLaMEM_Tests.py 
  expected_file = 't26_Dike/move2Dikes.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object     
  ex1 = pth.pthUnitTest('t26_moveDike',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


def heat_kfac():

  ranks = 2
  launch = '../bin/deb/LaMEM -ParamFile ./t26_Dike/dike_heating_kfac.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't26_Dike/dike_heating_kfac.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object                                                                                                                                                
  ex1 = pth.pthUnitTest('t26_dikeHeat_kfac',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


def heat_rhoA():

  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t26_Dike/dike_heating_rhoA.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't26_Dike/dike_heating_rhoA.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
    
  ex1 = pth.pthUnitTest('t26_dikeHeat_rhoA',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
