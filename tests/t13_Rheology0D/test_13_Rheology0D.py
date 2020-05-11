
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re
import sys
import subprocess
from t13_Rheology0D.CompareNumericsWithAnalytics import *


# This tests whether Maxwell viscoelasticity works
def ViscoElastic():
  # Visco-elastic rheology
  
  #==============================================
  # Run the input script wth matlab-generated particles
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_VE_0D.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't13_Rheology0D/Rheology_VE_0D-p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)
    #----------------------------  


    try: 
      data        = LoadTimeDependentData('Rheolog0D_VE');                                  # Load the data using the VTK toolbox
      YieldStress = 1e10;                                                     # large value, to deactivate it
      data        = AnalyticalSolution_VE(data,YieldStress);                  # Compute analytical solution & compute error
      PlotData(data,'./t13_Rheology0D/t13_ViscoElastic_output.png');          # Create Plot
  
      print('Created output figure ./t13_Rheology0D/t13_ViscoElastic_output.png comparing analytics vs. numerics')
    except:
      print('VTK/MatPlotLib/NumPy toolboxes are not installed; will not create plots')
    #----------------------------

  # Create unit test object
  ex1 = pth.pthUnitTest('t13_ViscoElastic',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


# This tests whether Maxwell viscoelastoplasticity works (with von Mises plasticity)
def ViscoElastoPlastic():
  # Visco-elasto-plastic rheology
  
  #==============================================
  # Run the input script wth matlab-generated particles
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_VEP_0D.dat' # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't13_Rheology0D/Rheology_VEP_0D-p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)
    #----------------------------  


    try: 
      data        = LoadTimeDependentData('Rheolog0D_VEP');                                        # Load the data using the VTK toolbox
      YieldStress = 1e7;                                                            # large value, to deactivate it
      data        = AnalyticalSolution_VE(data,YieldStress);                        # Compute analytical solution & compute error
      PlotData(data,'./t13_Rheology0D/t13_ViscoElastoPlasticMises_output.png');     # Create Plot
  
      print('Created output figure ./t13_Rheology0D/t13_ViscoElastoPlasticMises_output.png comparing analytics vs. numerics')
    except:
      print('VTK/MatPlotLib/NumPy toolboxes are not installed; will not create plots')
    #----------------------------

  # Create unit test object
  ex1 = pth.pthUnitTest('t13_ViscoElastoPlastic',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
