
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re
import subprocess

# Compute geodynamic sensitivity kernel for density, using a falling sphere setup
def rho_SensitivityKernel():

  # Run the input script wth matlab-generated particles
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t8_AdjointGradients/t8_AdjointGradients_SensitivityKernel.dat'
  expected_file = 't8_AdjointGradients/t8_Adjoint_rho_SensitivityKernel_p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)


  # Create unit test object
  ex1 = pth.pthUnitTest('t8_Adjoint_rho_SensitivityKernel',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

# Compute gradients for a falling sphere setup using ND units & all parameters
def FallingSphere_ND_all():

  # Run the input script wth matlab-generated particles
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t8_AdjointGradients/t8_AdjointGradients.dat'
  expected_file = 't8_AdjointGradients/t8_AdjointGradients_Sphere_ND_all.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|           delta(rho)[  1]")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|                  eta[  0]")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|   Velocity check            :")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|  adjoint     2:          eta[ 0]")
    unittest.compareFloatingPoint(key,1e-6)

  # Create unit test object
  ex1 = pth.pthUnitTest('t8_AdjointGradients_Sphere_ND_all',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


# Compute gradients for a falling sphere setup using ND units and both adjoint & FD
def FallingSphere_ND_CompareGradients_1():

  # Run the input script wth matlab-generated particles
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t8_AdjointGradients/t8_AdjointGradients_CompareGradients.dat'
  expected_file = 't8_AdjointGradients/t8_AdjointGradients_CompareGradients_1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|       FD     1:          eta[ 1]")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|  adjoint     2:          eta[ 1]")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|       FD     3:          eta[ 0]")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|  adjoint     4:          eta[ 0]")
    unittest.compareFloatingPoint(key,1e-8)

  # Create unit test object
  ex1 = pth.pthUnitTest('t8_AdjointGradients_CompareGradients_1',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

# Compute gradients for a falling sphere setup with nonlinear rheology using ND units and both adjoint & FD
def FallingSphere_ND_CompareGradients_2():

  # Run the input script wth matlab-generated particles
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t8_AdjointGradients/t8_AdjointGradients_CompareGradients_2.dat'
  expected_file = 't8_AdjointGradients/t8_AdjointGradients_CompareGradients_2.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|       FD     1:            n[ 0]")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|  adjoint     2:            n[ 0]")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|   Prefactor A               : ")
    unittest.compareFloatingPoint(key,1e-8)

  # Create unit test object
  ex1 = pth.pthUnitTest('t8_AdjointGradients_CompareGradients_2',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


# Compute gradients and scaling laws for a free-slip subduction setup with plasticity 
def SubductionSetup_Dimensional():

  # Run the input script wth matlab-generated particles
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t8_AdjointGradients/t8_Subduction2D_FreeSlip_DirectSolver.dat'
  expected_file = 't8_AdjointGradients/t8_Subduction2D_FreeSlip_DirectSolver_p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|   Prefactor A               : ")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|                  eta[  0]")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|           delta(rho)[  1]")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|           delta(rho)[  2]")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|                  eta[  2]")
    unittest.compareFloatingPoint(key,1e-8)
    
    key = re.escape(" |                  eta[  1]")
    unittest.compareFloatingPoint(key,1e-8)
    

  # Create unit test object
  ex1 = pth.pthUnitTest('t8_Adjoint_Subduction2D_FreeSlip',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)