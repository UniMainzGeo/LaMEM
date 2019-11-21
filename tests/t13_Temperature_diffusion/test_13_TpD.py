
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re


def test_a():

  # 1) Create a partitioning file and do not show any output of this
  os.system('../bin/opt/LaMEM -ParamFile ./t13_Temperature_diffusion/t13_Temperature_diffusion.dat -mode save_grid > /dev/null');

  # 2) Run MATLAB to create the Particles input (matlab should be in path)
  os.system('$MATLAB -nojvm -r "cd t13_Temperature_diffusion; Create_Marker; exit" > /dev/null')
  # Test a falling block case with build-in direct solver on 1 core, using optimized
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t13_Temperature_diffusion/t13_Temperature_diffusion.dat -mark_load_file ./markers_pT1/mdb' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't13_Temperature_diffusion/TpD_a.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('TpD_a',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

