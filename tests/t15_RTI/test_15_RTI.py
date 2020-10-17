import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re
import sys
import subprocess
from t15_RTI.CompareNumericsWithAnalytics import *


def RTI_isovisous_NoSlip():
  # No Slip

  #==============================================
  # Run the input script wth matlab-generated particles
  ranks = 1

  # This runs several executables for different strain rates & renames the directories (such that we can read them later with Python & create a plot)
  launch = ['../bin/opt/LaMEM -ParamFile ./t15_RTI/t15_RTI.dat -coord_x -0.125,0.125 -FreeSurf_Wavelength 0.25','rm -rf RTwav_0; mv Timestep_00000001_1.10000000e+00 RTwav_0', 
            '../bin/opt/LaMEM -ParamFile ./t15_RTI/t15_RTI.dat -coord_x -0.250,0.250 -FreeSurf_Wavelength 0.50','rm -rf RTwav_1; mv Timestep_00000001_1.10000000e+00 RTwav_1',
            '../bin/opt/LaMEM -ParamFile ./t15_RTI/t15_RTI.dat -coord_x -0.500,0.500 -FreeSurf_Wavelength 1.00','rm -rf RTwav_2; mv Timestep_00000001_1.10000000e+00 RTwav_2',
            '../bin/opt/LaMEM -ParamFile ./t15_RTI/t15_RTI.dat -coord_x -0.625,0.625 -FreeSurf_Wavelength 1.25','rm -rf RTwav_3; mv Timestep_00000001_1.10000000e+00 RTwav_3',
            '../bin/opt/LaMEM -ParamFile ./t15_RTI/t15_RTI.dat -coord_x -0.750,0.750 -FreeSurf_Wavelength 1.50','rm -rf RTwav_4; mv Timestep_00000001_1.10000000e+00 RTwav_4',
            '../bin/opt/LaMEM -ParamFile ./t15_RTI/t15_RTI.dat -coord_x -1.000,1.000 -FreeSurf_Wavelength 2.00','rm -rf RTwav_5; mv Timestep_00000001_1.10000000e+00 RTwav_5',
            '../bin/opt/LaMEM -ParamFile ./t15_RTI/t15_RTI.dat -coord_x -2.000,2.000 -FreeSurf_Wavelength 4.00','rm -rf RTwav_6; mv Timestep_00000001_1.10000000e+00 RTwav_6']
  

            
  # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't15_RTI/t15_RTI_IsoViscous-p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)
    #----------------------------  

    try: 
   
      # Load the data using the VTK toolbox; compute analytical solution & create plot
      data = LoadRTI_Data('RTI_test_surf_p00000000');
      AnalyticalSolution_FreeSlip(data);
      PlotRT_Data(data,'t15_RTI_isovisous_NoSlip');

      print('Created output figure ./t15_RTI/t15_RTI_isovisous_NoSlip.png comparing analytics vs. numerics')
    except:
      print('VTK/MatPlotLib/NumPy toolboxes are not installed; will not create plots')
    #----------------------------

  # Create unit test object
  ex1 = pth.pthUnitTest('t15_RTI_isovisous_NoSlip',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  
  return(ex1)