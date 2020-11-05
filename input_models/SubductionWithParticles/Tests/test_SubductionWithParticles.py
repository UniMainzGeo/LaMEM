
import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re
import subprocess

def makeLocalPathAbsolute(localRelPath) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  return(os.path.join(thisDir,localRelPath))

# This tests a subduction setup, in which we create the setup with a MATLAB script.
# It requires MATLAB to be present on the command-line & a way to ensure that this tests is only executed when this is the case
def test_1():
  # linear viscous setup of free subduction with free slip and a plastic crust

  # 1) Run MATLAB to create the Particles input on one processor (matlab should be in path)
  os.system('$MATLAB -nojvm -r "cd ./SubductionWithParticles; CreateMarkers_Subduction_Linear_FreeSlip_parallel; exit" > /dev/null')
  
  # Copy directory and give it a new name
  os.system('mv -f ./SubductionWithParticles/markers ./SubductionWithParticles/markers_linear')
 
  os.system('cp -rf ./SubductionWithParticles/markers_linear .')

  os.system('cd ./TestOutput')

  # Run the input script wth matlab-generated particles
  ranks = 1
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ./SubductionWithParticles/Subduction2D_FreeSlip_Particles_Linear_DirectSolver.dat -mark_load_file ./markers_linear/mdb -nstep_max 2 -dt_out 0') 
  expected_file = makeLocalPathAbsolute('Subduction2D_FreeSlip_Particles_Linear-MUMPS-p1.expected')
  

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('Sub1_Particles_Direct_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)




def test_2():
  # nonlinear setup of free subduction with a free surface and T-dependent viscoelastoplastic rheology

  # 1) Run MATLAB to create the Particles input on one processor (matlab should be in path)
  os.system('$MATLAB -nojvm -r "cd ./SubductionWithParticles; CreateMarkers_Subduction_Tdependent_FreeSurface_parallel; exit" > /dev/null')
 
   # Copy directory and give it a new name
  os.system('mv -f ./SubductionWithParticles/markers ./SubductionWithParticles/markers_nonlinear')

  os.system('rm -rf ./SubductionWithParticles/markers_nonlinear .')

  os.system('cp -r ./SubductionWithParticles/markers_nonlinear .')

  os.system('cd ./TestOutput')

  # Run the input script wth matlab-generated particles
  ranks = 1
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ./SubductionWithParticles/Subduction2D_FreeSurface_Particles_Nonlinear_DirectSolver.dat -mark_load_file ./markers_nonlinear/mdb -nstep_max 2 -dt_out 0') 
  expected_file = makeLocalPathAbsolute('Subduction2D_FreeSurface_Particles_Nonlinear-MUMPS-p1.expected')
  

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,3e-3)

  # Create unit test object
  ex1 = pth.pthUnitTest('Sub2_Particles_Direct_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)



def test_3():
  # This is a 3D example with multigrid solver (as in test1, bt 3D)

  # Copy directory and give it a new name
  os.system('cd ./SubductionWithParticles; mpiexec -n 4 ../../bin/opt/LaMEM -ParamFile ./Subduction3D_FreeSlip_MATLABParticles_Linear_Multigrid.dat -nel_x 64 -nel_y 32 -nel_z 16 -mode save_grid > /dev/null')


  # 2) Run MATLAB to create the Particles input on one processor (matlab should be in path)
  os.system('$MATLAB -nojvm -r "cd ./SubductionWithParticles; CreateMarkers_Subduction_Linear_FreeSlip_parallel_3D; exit" > /dev/null')
  
  # Copy directory and give it a new name
  os.system('mv -f ./SubductionWithParticles/markers ./SubductionWithParticles/markers_3D_linear')
 
  os.system('cp -rf ./SubductionWithParticles/markers_3D_linear .')

  os.system('cd ./TestOutput')

  # Run the input script wth matlab-generated particles
  ranks = 4
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ./SubductionWithParticles/Subduction3D_FreeSlip_MATLABParticles_Linear_Multigrid.dat -mark_load_file ./markers_3D_linear/mdb -nstep_max 2 -nel_x 64 -nel_y 32 -nel_z 16 -dt_out 0 -mode normal') 
  expected_file = makeLocalPathAbsolute('Subduction3D_FreeSlip_Particles_Linear-MUMPS-p4.expected')
  

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  ex1 = pth.pthUnitTest('Sub3D_Particles_Direct_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)
