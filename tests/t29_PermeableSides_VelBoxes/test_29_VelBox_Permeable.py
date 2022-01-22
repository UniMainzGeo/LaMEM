import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re

def test_permeableSides_VelBoxes():

  # Test dike feature using optimized LaMEM
    
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t29_PermeableSides_VelBoxes/VelBoxes_Permeable_sides.dat' # This must be a relative path with respect to runLaMEM_Tests.p
  expected_file = 't29_PermeableSides_VelBoxes/t29_PermeableSides_VelBoxes.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,5e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,2e-5)

  # Create unit test object
  ex1 = pth.pthUnitTest('t29_permeablesides_velboxes',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
