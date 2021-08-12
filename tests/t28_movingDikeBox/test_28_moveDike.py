import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re

def test_2cores_2dikes():

  # Test dike feature using optimized LaMEM
    
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t28_movingDikeBox/t28_move2Dikes.dat' # This must be a relative path with respect to runLaMEM_Tests.p
  expected_file = 't28_movingDikeBox/t28_Resultmove2Dikes.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t28_moveDike',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
