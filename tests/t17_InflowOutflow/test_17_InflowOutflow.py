
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re


def test_2D():

  # Test inflow/outflow conditions in 2D using optimized LaMEM
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t17_InflowOutflow/PlumeLithos_Interaction_2D.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't17_InflowOutflow/InflowOutflow-2D_p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

    key = re.escape("|eRes|_2")
    unittest.compareFloatingPoint(key,1e-7)
    
  # Create unit test object
  ex1 = pth.pthUnitTest('t17_InflowOutflow2D_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)



def test_3D():

  # Test inflow/outflow conditions in 2D using optimized LaMEM
  ranks = 4
  launch = '../bin/opt/LaMEM -ParamFile ./t17_InflowOutflow/PlumeLithos_Interaction_3D.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't17_InflowOutflow/InflowOutflow-3D_p4.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t17_InflowOutflow3D_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

def test_2D_Pres():

  # Test inflow/outflow conditions in 2D using optimized LaMEM
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t17_InflowOutflow/PlumeLithos_Interaction_2D_Perm.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't17_InflowOutflow/InflowOutflow-2D_Perm_p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

    key = re.escape("|eRes|_2")
    unittest.compareFloatingPoint(key,1e-7)
    
  # Create unit test object
  ex1 = pth.pthUnitTest('t17_InflowOutflow2D_Pres_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)



def test_3D_Pres():

  # Test inflow/outflow conditions in 2D using optimized LaMEM
  ranks = 4
  launch = '../bin/opt/LaMEM -ParamFile ./t17_InflowOutflow/PlumeLithos_Interaction_3D_Perm.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't17_InflowOutflow/InflowOutflow-3D_Perm_p4.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t17_InflowOutflow3D_Pres_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
