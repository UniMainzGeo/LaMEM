
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launch as launch
import re
import subprocess

# This tests a subduction setup, in which we create the setup with a MATLAB script.
# It requires MATLAB to be present on the command-line & a way to ensure that this tests is only executed when this is the case

# Tests are performed on 1 and 4 cores, using opt/deb

def test_a():
  
  # 1) Create a partitioning file and do not show any output of this
  os.system('../bin/opt/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_MATLAB_Particles_p1.dat -SavePartitioning 1 > /dev/null');
  
  # 2) Move partitioning file to directory we need it
  os.system('mv ProcessorPartitioning_1cpu_1.1.1.bin ./t3_SubductionMATLABinput/')
  
  # 3) Run MATLAB to create the Particles input  (requires the environmental variable $MATLAB to be defined!)
  os.system('$MATLAB -nojvm -r "cd t3_SubductionMATLABinput; CreateMarkers_Subduction(1); exit" > /dev/null')

  # Run the input script wth matlab-generated particles
  ranks = 1
  launch = '../bin/opt/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_MATLAB_Particles_p1.dat' 
  expected_file = 't3_SubductionMATLABinput/Sub1_MG_1core.expected'
  
  def comparefunc(unittest):
 
    key = 'Div_min'
    unittest.compareFloatingPoint(key,1e-7)
  
    key = 'Div_max'
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('unit_Sub1_a_MATLAB_Core1',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')
  
  return(ex1)

def test_b():
  
  # 1) Create a partitioning file and do not show any output of this
  os.system('mpiexec -n 4 ../bin/opt/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_MATLAB_Particles_p4.dat -SavePartitioning 1 > /dev/null');
  
  # 2) Move partitioning file to directory we need it
  os.system('mv ProcessorPartitioning_4cpu_4.1.1.bin ./t3_SubductionMATLABinput/')
  
  # 3) Run MATLAB to create the Particles input  (requires the environmental variable $MATLAB to be defined!)
  os.system('$MATLAB -nojvm -r "cd t3_SubductionMATLABinput; CreateMarkers_Subduction(4); exit" > /dev/null')

  # Run the input script wth matlab-generated particles
  ranks = 4
  launch = '../bin/opt/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_MATLAB_Particles_p4.dat -jp_pc_factor_mat_solver_package superlu_dist' 
  expected_file = 't3_SubductionMATLABinput/Sub1_MG_4core.expected'
  
  def comparefunc(unittest):
 
    key = 'Div_min'
    unittest.compareFloatingPoint(key,1e-7)
  
    key = 'Div_max'
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('unit_Sub1_b_MATLAB_Core4_SUPERLUDIST',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')
  
  return(ex1)

def test_c():
  
  # 1) Create a partitioning file and do not show any output of this
  os.system('mpiexec -n 4 ../bin/deb/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_MATLAB_Particles_p4.dat -SavePartitioning 1 > /dev/null');
  
  # 2) Move partitioning file to directory we need it
  os.system('mv ProcessorPartitioning_4cpu_4.1.1.bin ./t3_SubductionMATLABinput/')
  
  # 3) Run MATLAB to create the Particles input  (requires the environmental variable $MATLAB to be defined!)
  os.system('$MATLAB -nojvm -r "cd t3_SubductionMATLABinput; CreateMarkers_Subduction(4); exit" > /dev/null')

  # Run the input script wth matlab-generated particles
  ranks = 4
  launch = '../bin/deb/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_MATLAB_Particles_p4.dat -jp_pc_factor_mat_solver_package superlu_dist' 
  expected_file = 't3_SubductionMATLABinput/Sub1_MG_4core.expected'
  
  def comparefunc(unittest):
 
    key = 'Div_min'
    unittest.compareFloatingPoint(key,1e-7)
  
    key = 'Div_max'
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('unit_Sub1_c_MATLAB_Core4_SUPERLUDIST',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')
  
  return(ex1)

def test_d():
    # this tests a subduction setup with VEP rheology and pushing block
    
    
    # 1) Create a partitioning file and do not show any output of this
    os.system('mpiexec -n 8 ../bin/opt/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_VEP.dat -SavePartitioning 1 > /dev/null');
        
    # 2) Move partitioning file to directory we need it
    os.system('mv ProcessorPartitioning_8cpu_2.1.4.bin ./t3_SubductionMATLABinput/')
            
    # 3) Run MATLAB to create the Particles input  (requires the environmental variable $MATLAB to be defined!)
    os.system('$MATLAB -nojvm -r "cd t3_SubductionMATLABinput; CreateMarkers_SubductionVEP_parallel; exit" > /dev/null')
                
    # Run the input script wth matlab-generated particles
    ranks = 8
    launch = '../bin/opt/LaMEM -ParamFile ./t3_SubductionMATLABinput/Subduction_VEP.dat'
    expected_file = 't3_SubductionMATLABinput/Sub_VEP_MG_8core.expected'
                            
    def comparefunc(unittest):
                                
        key = 'Div_min'
        unittest.compareFloatingPoint(key,1e-7)
                                        
        key = 'Div_max'
        unittest.compareFloatingPoint(key,1e-7)
                                                
        key = re.escape("|Div|_2")
        unittest.compareFloatingPoint(key,1e-5)
                                                        
        key = re.escape("|mRes|_2")
        unittest.compareFloatingPoint(key,1e-4)
                                                                
    # Create unit test object
    ex1 = pth.pthUnitTest('unit_Sub1_d',ranks,launch,expected_file)
    ex1.setVerifyMethod(comparefunc)
    ex1.appendKeywords('@')
                                                                            
    return(ex1)
