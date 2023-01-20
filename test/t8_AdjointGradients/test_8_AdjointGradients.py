
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


# Compute gradients for a falling sphere setup using ND units and both adjoint & FD
def FallingSphere_GEO_CompareGradients_1():

  # Run the input script wth matlab-generated particles
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t8_AdjointGradients/t8_AdjointGradients_CompareGradients_geo.dat | grep "| "'
  expected_file = 't8_AdjointGradients/t8_AdjointGradients_CompareGradients_geo.expected'

  def comparefunc(unittest):


    key = re.escape("|       FD     1:          eta[ 1]")
    unittest.compareFloatingPoint(key,1e-30)

    key = re.escape("|  adjoint     2:          eta[ 1]")
    unittest.compareFloatingPoint(key,1e-30)

    key = re.escape("|       FD     3:          eta[ 0]")
    unittest.compareFloatingPoint(key,1e-28)

    key = re.escape("|  adjoint     4:          eta[ 0]")
    unittest.compareFloatingPoint(key,1e-28)

    key = re.escape("|           delta(rho)[  1]")
    unittest.compareFloatingPoint(key,1e-3)


  # Create unit test object
  ex1 = pth.pthUnitTest('t8_AdjointGradients_CompareGradients_geo',ranks,launch,expected_file)
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
    unittest.compareFloatingPoint(key,2e-6)

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
  launch = '../bin/opt/LaMEM -ParamFile ./t8_AdjointGradients/t8_Subduction2D_FreeSlip_DirectSolver.dat -nel_y 1 '
  expected_file = 't8_AdjointGradients/t8_Subduction2D_FreeSlip_DirectSolver_p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|mRes|_2")
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
    
    key = re.escape(" |      log10       eta[  0]")
    unittest.compareFloatingPoint(key,2e-3)
    
    key = re.escape(" |       FD     7:           fr[ 2]")
    unittest.compareFloatingPoint(key,1e-8)
    
  
  # Create unit test object
  ex1 = pth.pthUnitTest('t8_Adjoint_Subduction2D_FreeSlip',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)



# Compute FD and adjoint gradients for principal stress directions
def PSD_ND():

  # Run the input script wth matlab-generated particles
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t8_AdjointGradients/t8_FB_PSDTest.dat -nel_x 8 -nel_y 8 -nel_z 8  | grep "| "'
  expected_file = 't8_AdjointGradients/t8_FB_PSDTest_p1.expected'

  def comparefunc(unittest):

    key = re.escape("|  adjoint     1:          rho[ 2]")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|       FD     2:          rho[ 2]")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("| Current Cost function = ")
    unittest.compareFloatingPoint(key,1e-6)


  # Create unit test object
  ex1 = pth.pthUnitTest('t8_Adjoint_PSD',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)





# Compute geodynamic sensitivity kernel for density & PSD
def rho_SensitivityKernel_PSD():

  # Run the input script wth matlab-generated particles
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t8_AdjointGradients/t8_AdjointGradients_SensitivityKernel_PSD.dat'
  expected_file = 't8_AdjointGradients/t8_Adjoint_rho_SensitivityKernel_PSD_p2.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

    key = re.escape("| Current Cost function = ")
    unittest.compareFloatingPoint(key,1e-6)

  # Create unit test object
  ex1 = pth.pthUnitTest('t8_Adjoint_rho_SensitivityKernel_PSD',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)




# Compute geodynamic sensitivity kernel for n & PSD from Georg's paper (at lower res)
def n_SensitivityKernelPaper_PSD():

  # Run the input script wth matlab-generated particles
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t8_AdjointGradients/t8_PSDKernelPaper.dat  -nel_x 8  -nel_y 8 -nel_z 8 | grep "| "'
  expected_file = 't8_AdjointGradients/t8_Adjoint_n_SensitivityKernelPaper_PSD.expected'

  def comparefunc(unittest):

    key = re.escape("|   Norm of field gradient vector :")
    unittest.compareFloatingPoint(key,5e-1)

  # Create unit test object
  ex1 = pth.pthUnitTest('t8_Adjoint_n_SensitivityKernel_PSD',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

# Compute geodynamic sensitivity kernel for eta0 & PSD from Georg's paper (at lower res)
def eta0_SensitivityKernelPaper_PSD():

  # Run the input script wth matlab-generated particles
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t8_AdjointGradients/t8_PSDKernelPaper.dat  -nel_x 8  -nel_y 8 -nel_z 8 -Type[0] eta0 | grep "| "'
  expected_file = 't8_AdjointGradients/t8_Adjoint_eta0_SensitivityKernelPaper_PSD.expected'

  def comparefunc(unittest):

    key = re.escape("|   Norm of field gradient vector :")
    unittest.compareFloatingPoint(key,5e-1)

  # Create unit test object
  ex1 = pth.pthUnitTest('t8_Adjoint_eta0_SensitivityKernel_PSD',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

