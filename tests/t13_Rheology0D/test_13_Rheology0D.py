
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
  launch = ['rm -r Timestep*',
            'rm -rf ./t13_Rheology0D/Out_VE_0D; mkdir ./t13_Rheology0D/Out_VE_0D',
            '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_VE_0D.dat',
            'mv Timestep* ./t13_Rheology0D/Out_VE_0D'] # This must be a relative path with respect to runLaMEM_Tests.py
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
      data        = LoadTimeDependentData('./t13_Rheology0D/Out_VE_0D','Rheolog0D_VE')                                  # Load the data using the VTK toolbox
      YieldStress = 1e10                                                     # large value, to deactivate it
      data        = AnalyticalSolution_VE(data,YieldStress)                  # Compute analytical solution & compute error
      PlotTimeDependentData(data,'./t13_Rheology0D/t13_ViscoElastic_output.png')          # Create Plot
  
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
  # Visco-elasto-plastic rheology with linear viscosity
  
  #==============================================
  # Run the input script wth matlab-generated particles
  ranks = 1
  launch = ['rm -rf ./t13_Rheology0D/Out_VEP_0D; mkdir ./t13_Rheology0D/Out_VEP_0D',
            '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_VEP_0D.dat',
             'mv Timestep* ./t13_Rheology0D/Out_VEP_0D'] # This must be a relative path with respect to runLaMEM_Tests.py
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
      data        = LoadTimeDependentData('./t13_Rheology0D/Out_VEP_0D','Rheolog0D_VEP')                                        # Load the data using the VTK toolbox
      YieldStress = 1e7                                                            # large value, to deactivate it
      data        = AnalyticalSolution_VE(data,YieldStress)                        # Compute analytical solution & compute error
      PlotTimeDependentData(data,'./t13_Rheology0D/t13_ViscoElastoPlasticMises_output.png')     # Create Plot
  
      print('Created output figure ./t13_Rheology0D/t13_ViscoElastoPlasticMises_output.png comparing analytics vs. numerics')
    except:
      print('VTK/MatPlotLib/NumPy toolboxes are not installed; will not create plots')
    #----------------------------

  # Create unit test object
  ex1 = pth.pthUnitTest('t13_ViscoElastoPlastic',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


# This tests whether Maxwell viscoelastoplasticity works (with von Mises plasticity)
def ViscoElastoPlastic_DislocationCreep():
  # Visco-elasto-plastic rheology with nonlinear dislocation creep viscosity
  
  #==============================================
  # Run the input script wth matlab-generated particles
  ranks = 1
  launch = ['rm -rf ./t13_Rheology0D/Out_DisCreep_VEP_0D; mkdir ./t13_Rheology0D/Out_DisCreep_VEP_0D',
            '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_DislocationCreep_VEP_0D.dat',
            'mv Timestep* ./t13_Rheology0D/Out_DisCreep_VEP_0D'] # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't13_Rheology0D/Rheology_DislocationCreep_VEP_0D-p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)
    #----------------------------  
    
  
    try: 
      data        = LoadTimeDependentData('./t13_Rheology0D/Out_DisCreep_VEP_0D','Rheolog0D_DislocationCreep_VEP')        # Load the data using the VTK toolbox
      YieldStress = 15e6                                                            # large value, to deactivate it
      data        = AnalyticalSolution_DislocationCreep_VEP(data,YieldStress)       # Compute analytical solution & compute error

      PlotTimeDependentData(data,'./t13_Rheology0D/t13_DislocationCreep_ViscoElastoPlasticMises_output.png')     # Create Plot
  
      print('Created output figure ./t13_Rheology0D/t13_DislocationCreep_ViscoElastoPlasticMises_output.png comparing analytics vs. numerics')
    except:
      print('VTK/MatPlotLib/NumPy toolboxes are not installed; will not create plots')
    #----------------------------

  # Create unit test object
  ex1 = pth.pthUnitTest('t13_ViscoElastoPlastic_DislocationCreep',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)

# This tests whether Maxwell viscoelasticity works (with nonlinear creep laws, requiring local iterations)
def ViscoElastic_DislocationCreep():
  # Visco-elastic rheology with nonlinear dislocation creep viscosity
  
  #==============================================
  # Run the input script wth matlab-generated particles
  ranks = 1
  launch = ['rm -rf ./t13_Rheology0D/Out_DisCreep_VE_0D; mkdir ./t13_Rheology0D/Out_DisCreep_VE_0D',
            '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_DislocationCreep_VE_0D.dat',
            'mv Timestep* ./t13_Rheology0D/Out_DisCreep_VE_0D'] # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't13_Rheology0D/Rheology_DislocationCreep_VE_0D-p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)
    #----------------------------  
    
  
    try: 
      data        = LoadTimeDependentData('./t13_Rheology0D/Out_DisCreep_VE_0D','Rheolog0D_DislocationCreep_VE')                           # Load the data using the VTK toolbox
      YieldStress = 1e10                                                                             # large value, to deactivate it
      data        = AnalyticalSolution_DislocationCreep_VEP(data,YieldStress)                        # Compute analytical solution & compute error

      PlotTimeDependentData(data,'./t13_Rheology0D/t13_DislocationCreep_ViscoElastic_output.png')    # Create Plot
  
      print('Created output figure ./t13_Rheology0D/t13_DislocationCreep_ViscoElastic_output.png comparing analytics vs. numerics')
    except:
      print('VTK/MatPlotLib/NumPy toolboxes are not installed; will not create plots')
    #----------------------------

  # Create unit test object
  ex1 = pth.pthUnitTest('t13_ViscoElastic_DislocationCreep',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


def LinearViscous():
  # linear viscous rheology

  # This computes a solution for different applied background strainrate values
  # The analytical solution is simply T2nd = 2*eta*E2nd, where E2nd are the applied strain rate values
  
  #==============================================
  # Run the input script wth matlab-generated particles
  ranks = 1

  # This runs several executables for different strain rates & renames the directories (such that we can read them later with Python & create a plot)
  launch = ['rm -rf ./t13_Rheology0D/Out_LV; mkdir ./t13_Rheology0D/Out_LV',
            '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_linearViscous_0D.dat -exx_strain_rates  -1e-13','mv Timestep_00000001_2.00000000e-03 ./t13_Rheology0D/Out_LV/Strainrate_0_13', 
            '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_linearViscous_0D.dat -exx_strain_rates  -1e-14','mv Timestep_00000001_2.00000000e-03 ./t13_Rheology0D/Out_LV/Strainrate_1_14', 
            '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_linearViscous_0D.dat -exx_strain_rates  -1e-15','mv Timestep_00000001_2.00000000e-03 ./t13_Rheology0D/Out_LV/Strainrate_2_15', 
            '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_linearViscous_0D.dat -exx_strain_rates  -1e-16','mv Timestep_00000001_2.00000000e-03 ./t13_Rheology0D/Out_LV/Strainrate_3_16', 
            '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_linearViscous_0D.dat -exx_strain_rates  -1e-17','mv Timestep_00000001_2.00000000e-03 ./t13_Rheology0D/Out_LV/Strainrate_4_17']
  
  
  # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't13_Rheology0D/Rheology_linearViscous_0D-p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)
    #----------------------------  

  
    try: 
      # Load the data using the VTK toolbox; compute analtical 
      data        = LoadStrainrateData('./t13_Rheology0D/Out_LV/','Rheolog0D_linearViscous')
      data        = AnalyticalSolution_linearViscous(data)
      PlotStrainrateData(data,'./t13_Rheology0D/t13_linearViscous_output.png')     # Create Plot
   
      print('Created output figure ./t13_Rheology0D/t13_linearViscous_output.png comparing analytics vs. numerics')
    except:
      print('VTK/MatPlotLib/NumPy toolboxes are not installed; will not create plots')
    #----------------------------

  # Create unit test object
  ex1 = pth.pthUnitTest('t13_linearViscous',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  
  return(ex1)


def DislocationCreeplaw():
  # Dislocation Creeplaw (same as in the Gerya textbook)

  # This computes a solution for different applied background strainrate values
  # The analytical solution is simply T2nd = 2*eta*E2nd, where E2nd are the applied strain rate values
  
  #==============================================
  # Run the input script wth matlab-generated particles
  ranks = 1

  # This runs several executables for different strain rates & renames the directories (such that we can read them later with Python & create a plot)
  launch = ['rm -rf ./t13_Rheology0D/Out_DisCr; mkdir ./t13_Rheology0D/Out_DisCr',
            '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_PowerlawCreep_DryOlivine_0D.dat -exx_strain_rates  -1e-13','mv Timestep_00000001_2.00000000e-03 ./t13_Rheology0D/Out_DisCr/Strainrate_0_13', 
            '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_PowerlawCreep_DryOlivine_0D.dat -exx_strain_rates  -1e-14','mv Timestep_00000001_2.00000000e-03 ./t13_Rheology0D/Out_DisCr/Strainrate_1_14', 
            '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_PowerlawCreep_DryOlivine_0D.dat -exx_strain_rates  -1e-15','mv Timestep_00000001_2.00000000e-03 ./t13_Rheology0D/Out_DisCr/Strainrate_2_15', 
            '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_PowerlawCreep_DryOlivine_0D.dat -exx_strain_rates  -1e-16','mv Timestep_00000001_2.00000000e-03 ./t13_Rheology0D/Out_DisCr/Strainrate_3_16', 
            '../bin/opt/LaMEM -ParamFile ./t13_Rheology0D/Rheology_PowerlawCreep_DryOlivine_0D.dat -exx_strain_rates  -1e-17','mv Timestep_00000001_2.00000000e-03 ./t13_Rheology0D/Out_DisCr/Strainrate_4_17']
  
  
  # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't13_Rheology0D/Rheology_DislocationCreeplaw_0D-p1.expected'

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
      data        = LoadStrainrateData('./t13_Rheology0D/Out_DisCr/','Rheolog0D_DryOlivine')
      data        = AnalyticalSolution_DislocationCreep(data, 'DryOlivine')
      PlotStrainrateData(data,'./t13_Rheology0D/t13_DislocationCreeplaw_DryOlivine_output.png')     # Create Plot

      print('Created output figure ./t13_Rheology0D/t13_DislocationCreeplaw_DryOlivine_output.png comparing analytics vs. numerics')
    except:
      print('VTK/MatPlotLib/NumPy toolboxes are not installed; will not create plots')
    #----------------------------

  # Create unit test object
  ex1 = pth.pthUnitTest('t13_DislocationCreeplaw',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  
  return(ex1)