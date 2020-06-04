
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re
import subprocess

# This tests the inversion possibilities 

# Subduction setup using adjoint gradients using internal Gradient descent algorithm
def SubductionInversion_GD():

  # Note that we run this at a low resolution to speed up testing & that we only compare the 
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile t7_AdjointGradientInversion/t7_Subduction2D_FreeSlip_Inversion.dat -nel_z 16 -nel_x 64 | grep "| "'
  expected_file = 't7_AdjointGradientInversion/t7_AdjointGradientInversion_1.expected'

  def comparefunc(unittest):

    key = re.escape("F / FINI =")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("| 1. Diff parameter value =")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("| 2. Diff parameter value =")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("| 1. Parameter value =")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("| 2. Parameter value =")
    unittest.compareFloatingPoint(key,1e-5)

  # Create unit test object
  ex1 = pth.pthUnitTest('t7_AdjointGradientInversion_1',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

# Subduction setup using adjoint gradients using TAO
def SubductionInversion_TAO():

  # Note that we run this at a low resolution to speed up testing & that we only compare the 
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile t7_AdjointGradientInversion/t7_Subduction2D_FreeSlip_Inversion.dat -tao_fmin 1e-6 -nel_z 16 -nel_x 64 -Inversion_EmployTAO 1 | grep "| "'
  expected_file = 't7_AdjointGradientInversion/t7_AdjointGradientInversion_2.expected'

  def comparefunc(unittest):

    key = re.escape("| misfit           =")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("| misfit / misfit0 =")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|   1. eta[0] =")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|   2. rho[1] =")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|  adjoint     1:   log10  eta[ 0]")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|  adjoint     2:          rho[ 1]")
    unittest.compareFloatingPoint(key,1e-5)

  # Create unit test object
  ex1 = pth.pthUnitTest('t7_AdjointGradientInversion_2',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


# Subduction setup using FD gradients
def SubductionInversion_FD_TAO():

  # Note that we run this at a low resolution to speed up testing & that we only compare the 
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile t7_AdjointGradientInversion/t7_Subduction2D_FreeSlip_Inversion_FD.dat -tao_fmin 1e-6 -nel_z 16 -nel_x 32 | grep "| "'
  expected_file = 't7_AdjointGradientInversion/t7_AdjointGradientInversion_3.expected'

  def comparefunc(unittest):

    key = re.escape("| misfit           =")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("| misfit / misfit0 =")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|   1. eta[0] =")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|   2. rho[1] =")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|       FD     1:   log10  eta[ 0]")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|       FD     2:          rho[ 1]")
    unittest.compareFloatingPoint(key,1e-5)

  # Create unit test object
  ex1 = pth.pthUnitTest('t7_AdjointGradientInversion_3',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)



# falling sphere ND test