import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re

# This file performs tests for the examples in thsi directory

def makeLocalPathAbsolute(localRelPath) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  return(os.path.join(thisDir,localRelPath))

#------------------------------------------------------------------------------------------------
def test_1():

  # first test with no regularization
  ranks  = 2
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ../Localization/localization.dat -nstep_max 25 -dt_out 0 -nstep_ini 0 -nel_x 64 -nel_z 16') 
  expected_file = makeLocalPathAbsolute('Localization_noRegularization-p2.expected')

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = 'CONVERGED_RTOL iterations'
    unittest.compareInteger(key,0)

  # Create unit test object
  test = pth.pthUnitTest('Localization_noReg',ranks,launch,expected_file)
  test.setVerifyMethod(comparefunc)
  test.appendKeywords('@')
  test.setUseSandbox()      # put test output in seperate directory

  return(test)

#------------------------------------------------------------------------------------------------
def test_2():

  # second test WITH regularization
  ranks  = 2
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ../Localization/localization_eta_min_reg.dat -nstep_max 25 -dt_out 0 -nstep_ini 0 -nel_x 64 -nel_z 16') 
  expected_file = makeLocalPathAbsolute('Localization_Regularization_etaMin-p2.expected')

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = 'CONVERGED_RTOL iterations'
    unittest.compareInteger(key,0)

  # Create unit test object
  test = pth.pthUnitTest('Localization_Reg_EtaMin',ranks,launch,expected_file)
  test.setVerifyMethod(comparefunc)
  test.appendKeywords('@')
  test.setUseSandbox()      # put test output in seperate directory

  return(test)
#------------------------------------------------------------------------------------------------
def test_3():

  # 3rd, with phase-wise regularization
  ranks  = 2
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ../Localization/localization_eta_st_reg.dat -nstep_max 25 -dt_out 0 -nstep_ini 0 -nel_x 64 -nel_z 16') 
  expected_file = makeLocalPathAbsolute('Localization_Regularization_etaSt-p2.expected')

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = 'CONVERGED_RTOL iterations'
    unittest.compareInteger(key,0)

  # Create unit test object
  test = pth.pthUnitTest('Localization_Reg_EtaSt',ranks,launch,expected_file)
  test.setVerifyMethod(comparefunc)
  test.appendKeywords('@')
  test.setUseSandbox()      # put test output in seperate directory

  return(test)