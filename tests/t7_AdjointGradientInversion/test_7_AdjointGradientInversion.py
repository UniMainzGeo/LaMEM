
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
    unittest.compareFloatingPoint(key,1e-3)

    key = re.escape("| misfit / misfit0 =")
    unittest.compareFloatingPoint(key,1e-3)

    key = re.escape("|   1. eta[0] =")
    unittest.compareFloatingPoint(key,1e-3)

    key = re.escape("|   2. rho[1] =")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|       FD     1:   log10  eta[ 0]")
    unittest.compareFloatingPoint(key,1e-3)

    key = re.escape("|       FD     2:          rho[ 1]")
    unittest.compareFloatingPoint(key,1e-5)

  # Create unit test object
  ex1 = pth.pthUnitTest('t7_AdjointGradientInversion_3',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)



# PSD paper inversion for nonlinear materials:
def PSD_Paper_GD_Nonlinear():

  # Note that we run this at a low resolution to speed up testing & that we only compare the 
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile t7_AdjointGradientInversion/t7_PSDInversionPaper.dat -Inversion_rtol 4.6e-2 | grep "| "'
  expected_file = 't7_AdjointGradientInversion/t7_PSDInversionPaper_1.expected'

  def comparefunc(unittest):

    key = re.escape("| LS factor for 1.Parameter = ")
    unittest.compareFloatingPoint(key,1e-1)
    
    key = re.escape("|    F =")
    unittest.compareFloatingPoint(key,1e-3)

    key = re.escape("| 1. Parameter value = ")
    unittest.compareFloatingPoint(key,1e-5)

  # Create unit test object
  ex1 = pth.pthUnitTest('t7_PSDInversion_1',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


# PSD paper inversion for linear materials:
def PSD_Paper_GD_Linear():

  # Note that we run this at a low resolution to speed up testing & that we only compare the 
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile t7_AdjointGradientInversion/t7_PSDInversionPaper.dat  -nel_x 8 -nel_y 8 -nel_z 8  -n[0] 1 -n[1] 1 -n[2] 1  -Value[0] 135 | grep "| "'
  expected_file = 't7_AdjointGradientInversion/t7_PSDInversionPaper_2.expected'

  def comparefunc(unittest):

    key = re.escape("| LS factor for 1.Parameter = ")
    unittest.compareFloatingPoint(key,1e-3)
    
    key = re.escape("|    F =")
    unittest.compareFloatingPoint(key,1e-3)

    key = re.escape("| 1. Parameter value = ")
    unittest.compareFloatingPoint(key,1e-5)

  # Create unit test object
  ex1 = pth.pthUnitTest('t7_PSDInversion_2',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
