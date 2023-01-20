
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re
import subprocess

# This tests a subduction setup, in which we create the setup with a MATLAB script.
# It requires MATLAB to be present on the command-line & a way to ensure that this tests is only executed when this is the case

# Tests are performed on 8 cores, using opt/deb

def test_a():

  # 1) Create a partitioning file and do not show any output of this
  os.system('mpiexec -n 8 ../bin/opt/LaMEM -ParamFile ./t24_Erosion_Sedimentation/Erosion_Sedimentation_2D.dat -mode save_grid > /dev/null');

  # 2) Run MATLAB to create the Particles input  (requires the environmental variable $MATLAB to be defined!)
  os.system('$MATLAB -nojvm -r "cd t24_Erosion_Sedimentation; CreateMarkers_2D; exit" > /dev/null')

  # Run the input script wth matlab-generated particles
  ranks = 8
  launch = '../bin/opt/LaMEM -ParamFile ./t24_Erosion_Sedimentation/Erosion_Sedimentation_2D.dat'
  expected_file = 't24_Erosion_Sedimentation/Erosion_Sedimentation_2D_opt-p8.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,2.5e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t24_Erosion_Sedimentation_2D_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

def test_b():

  # 1) Create a partitioning file and do not show any output of this
  os.system('mpiexec -n 8 ../bin/deb/LaMEM -ParamFile ./t24_Erosion_Sedimentation/Erosion_Sedimentation_2D.dat -mode save_grid > /dev/null');

  # 2) Run MATLAB to create the Particles input  (requires the environmental variable $MATLAB to be defined!)
  os.system('$MATLAB -nojvm -r "cd t24_Erosion_Sedimentation; CreateMarkers_2D; exit" > /dev/null')

  # Run the input script wth matlab-generated particles
  ranks = 8
  launch = '../bin/deb/LaMEM -ParamFile ./t24_Erosion_Sedimentation/Erosion_Sedimentation_2D.dat'
  expected_file = 't24_Erosion_Sedimentation/Erosion_Sedimentation_2D_deb-p8.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,2.5e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('t24_Erosion_Sedimentation_2D_deb',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
