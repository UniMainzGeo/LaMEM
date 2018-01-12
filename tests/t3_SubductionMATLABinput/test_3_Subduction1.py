
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re
import subprocess

# This tests a subduction setup, in which we create the setup with a MATLAB script.
# It requires MATLAB to be present on the command-line & a way to ensure that this tests is only executed when this is the case

# Tests are performed on 1 and 4 cores, using opt/deb

def test_a():

  # 1) Create a partitioning file and do not show any output of this
  os.system('../bin/opt/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_MATLAB_Particles.dat -mode save_grid > /dev/null');

  # 2) Run MATLAB to create the Particles input (matlab should be in path)
  os.system('$MATLAB -nojvm -r "cd t3_SubductionMATLABinput; CreateMarkers_Subduction(1); exit" > /dev/null')

  # Run the input script wth matlab-generated particles
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_MATLAB_Particles.dat -mark_load_file ./markers_p1/mdb'
  expected_file = 't3_SubductionMATLABinput/Sub1_MATLAB_a_Direct_opt-p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('Sub1_MATLAB_a_Direct_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

def test_b():

  # 1) Create a partitioning file and do not show any output of this
  os.system('mpiexec -n 4 ../bin/opt/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_MATLAB_Particles.dat -mode save_grid > /dev/null');

  # 2) Run MATLAB to create the Particles input  (requires the environmental variable $MATLAB to be defined!)
  os.system('$MATLAB -nojvm -r "cd t3_SubductionMATLABinput; CreateMarkers_Subduction(4); exit" > /dev/null')

  # Run the input script wth matlab-generated particles
  ranks = 4
  launch = '../bin/opt/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_MATLAB_Particles.dat -mark_load_file ./markers_p4/mdb -jp_pc_factor_mat_solver_package superlu_dist'
  expected_file = 't3_SubductionMATLABinput/Sub1_MATLAB_b_SUPERLUDIST_opt-p4.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('Sub1_MATLAB_b_SUPERLUDIST_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

def test_c():

  # 1) Create a partitioning file and do not show any output of this
  os.system('mpiexec -n 4 ../bin/deb/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_MATLAB_Particles.dat -mode save_grid > /dev/null');

  # 2) Run MATLAB to create the Particles input  (requires the environmental variable $MATLAB to be defined!)
  os.system('$MATLAB -nojvm -r "cd t3_SubductionMATLABinput; CreateMarkers_Subduction(4); exit" > /dev/null')

  # Run the input script wth matlab-generated particles
  ranks = 4
  launch = '../bin/deb/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_MATLAB_Particles.dat -mark_load_file ./markers_p4/mdb -jp_pc_factor_mat_solver_package superlu_dist'
  expected_file = 't3_SubductionMATLABinput/Sub1_MATLAB_c_SUPERLUDIST_deb-p4.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('Sub1_MATLAB_c_SUPERLUDIST_deb',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


def test_d():

  # 1) Create a partitioning file and do not show any output of this
  os.system('mpiexec -n 8 ../bin/opt/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_VEP.dat -mode save_grid > /dev/null');

  # 2) Run MATLAB to create the Particles input  (requires the environmental variable $MATLAB to be defined!)
  os.system('$MATLAB -nojvm -r "cd t3_SubductionMATLABinput; CreateMarkers_SubductionVEP_parallel; exit" > /dev/null')

  # Run the input script wth matlab-generated particles
  ranks = 8
  launch = '../bin/opt/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_VEP.dat'
  expected_file = 't3_SubductionMATLABinput/Sub1_MATLAB_d_MUMPS_MG_VEP_opt-p8.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('Sub1_MATLAB_d_MUMPS_MG_VEP_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

