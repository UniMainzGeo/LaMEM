
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re
import subprocess

# This tests a simple Rayleigh Taylor instability as discussed in the Reuber et al. scaling law paper

# Test is performed on 2 cores, using opt and is valid for the Hard-film RTI
def test_RTI_1():

  # Run the input script wth matlab-generated particles
  ranks = 2
  launch = '../bin/opt/LaMEM -ParamFile ./t6_AdjointGradientScaling/t6_RTI_ScalingLaw.dat'
  expected_file = 't6_AdjointGradientScaling/t6_AdjointGradientScaling_p2.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

    key = re.escape("|   Prefactor A               :")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|   Velocity check            :")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|  adjoint     1:          rho[ 0]")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|  adjoint     2:          eta[ 0]")
    unittest.compareFloatingPoint(key,1e-8)
    
    key = re.escape("|  adjoint     4:          eta[ 1]")
    unittest.compareFloatingPoint(key,1e-8)
    
  # Create unit test object
  ex1 = pth.pthUnitTest('t6_AdjointGradientScalingLaws_p2',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


# Test is performed on 1 cores, using opt and is valid for the Soft-film RTI
def test_RTI_2():

  # Run the input script wth matlab-generated particles
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t6_AdjointGradientScaling/t6_RTI_ScalingLaw.dat -surf_level 0.1 -eta[0] 10 -eta[1] 1 -coord_x -0.4,0.4 -FreeSurf_Wavelength 0.8'
  expected_file = 't6_AdjointGradientScaling/t6_AdjointGradientScaling_SoftFilm_p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

    key = re.escape("|   Prefactor A               :")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|   Velocity check            :")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|  adjoint     1:          rho[ 0]")
    unittest.compareFloatingPoint(key,1e-8)

    key = re.escape("|  adjoint     2:          eta[ 0]")
    unittest.compareFloatingPoint(key,1e-8)
    
    key = re.escape("|  adjoint     4:          eta[ 1]")
    unittest.compareFloatingPoint(key,1e-8)
    
  # Create unit test object
  ex1 = pth.pthUnitTest('t6_AdjointGradientScalingLaws_SoftFilm',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)