
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

