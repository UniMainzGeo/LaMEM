
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launch as launch
import re


def test_a():

  # Test a falling block case on 1 core, using optimized LaMEM
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t2_FB2_MG/FallingBlock_mono_CoupledMG_RedundantCoarse.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't2_FB2_MG/FB2_MG_1core.expected'

  def comparefunc(unittest):

    key = '|Div|_inf'
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('unit_FB2_a_CoupledMG_Core1',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
