import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re

# This file performs tests for the examples in thsi directory

def makeLocalPathAbsolute(localRelPath) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  return(os.path.join(thisDir,localRelPath))

#------------------------------------------------------------------------------------------------
def test_1():

  # Test a falling block case with build-in direct solver on 1 core, using optimized LaMEM
  ranks  = 2
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ../BuildInSetups/FallingBlock_DirectSolver.dat -nstep_max 3 -dt_out 0 -nstep_ini 0') 
  expected_file = makeLocalPathAbsolute('FallingBlock_DirectSolver-MUMPS-p2.expected')

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = 'CONVERGED_RTOL iterations'
    unittest.compareInteger(key,0)

  # Create unit test object
  test = pth.pthUnitTest('FallingBlock_Direct_MUMPS_serial',ranks,launch,expected_file)
  test.setVerifyMethod(comparefunc)
  test.appendKeywords('@')
  test.setUseSandbox()      # put test output in seperate directory

  return(test)

#------------------------------------------------------------------------------------------------
def test_2():

  # Test a falling block case with SUPERLU_DIST direct solver
  ranks = 1
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ../BuildInSetups/FallingBlock_DirectSolver.dat -nstep_max 3 -dt_out 0 -nstep_ini 0 -jp_pc_factor_mat_solver_package superlu_dist') 
  expected_file = makeLocalPathAbsolute('FallingBlock_DirectSolver-SUPERLU_DIST-p1.expected')

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = 'CONVERGED_RTOL iterations'
    unittest.compareInteger(key,0)

  # Create unit test object
  test = pth.pthUnitTest('FallingBlock_Direct_SUPERLU_DIST',ranks,launch,expected_file)
  test.setVerifyMethod(comparefunc)
  test.appendKeywords('@')
  test.setUseSandbox()      # put test output in seperate directory

  return(test)

#------------------------------------------------------------------------------------------------
def test_3():

  # Test a falling block case with build-in iterative solver on 1 core
  ranks = 1
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ../BuildInSetups/FallingBlock_IterativeSolver.dat -nstep_max 3') 
  expected_file = makeLocalPathAbsolute('FallingBlock_Iterative-p1.expected')

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = 'CONVERGED_RTOL iterations'
    unittest.compareInteger(key,0)

  # Create unit test object
  test = pth.pthUnitTest('FallingBlock_Iterative',ranks,launch,expected_file)
  test.setVerifyMethod(comparefunc)
  test.appendKeywords('@')
  test.setUseSandbox()      # put test output in seperate directory

  return(test)

#------------------------------------------------------------------------------------------------
def test_4():

  # Test a falling block case with build-in iterative solver on 2 cores
  ranks = 2
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ../BuildInSetups/FallingBlock_IterativeSolver.dat -nstep_max 3') 
  expected_file = makeLocalPathAbsolute('FallingBlock_Iterative_p2.expected')

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = 'CONVERGED_RTOL iterations'
    unittest.compareInteger(key,0)

  # Create unit test object
  test = pth.pthUnitTest('FallingBlock_Iterative_parallel',ranks,launch,expected_file)
  test.setVerifyMethod(comparefunc)
  test.appendKeywords('@')
  test.setUseSandbox()      # put test output in seperate directory

  return(test)

#------------------------------------------------------------------------------------------------
def test_5():
  # Falling block setup with multigrid solver
  ranks = 2
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ../BuildInSetups/FallingBlock_Multigrid.dat -nstep_max 3 -dt_out 0 -nstep_ini 0') 
  expected_file = makeLocalPathAbsolute('FallingBlock_Multigrid-p2.expected')

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = 'CONVERGED_RTOL iterations'
    unittest.compareInteger(key,0)

    key = 'KSP Residual norm'
    unittest.compareFloatingPoint(key,1e-6)
    
  # Create unit test object
  test = pth.pthUnitTest('FallingBlock_Multigrid',ranks,launch,expected_file)
  test.setVerifyMethod(comparefunc)
  test.appendKeywords('@')
  test.setUseSandbox()      # put test output in seperate directory

  return(test)


#------------------------------------------------------------------------------------------------
def test_6():
  # Falling spheres setup with multigrid solver
  ranks = 2
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ../BuildInSetups/FallingSpheres_Multigrid.dat -dt_out 0 -nstep_ini 0') 
  expected_file = makeLocalPathAbsolute('FallingSpheres_Multigrid-p2.expected')

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = 'CONVERGED_RTOL iterations'
    unittest.compareInteger(key,0)

    key = 'KSP Residual norm'
    unittest.compareFloatingPoint(key,1e-6)
    
  # Create unit test object
  test = pth.pthUnitTest('FallingSpheres_Multigrid',ranks,launch,expected_file)
  test.setVerifyMethod(comparefunc)
  test.appendKeywords('@')
  test.setUseSandbox()      # put test output in seperate directory

  return(test)

#------------------------------------------------------------------------------------------------
def test_7():
  # Falling spheres setup with multigrid solver but jacobi smoother
  ranks = 2
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ../BuildInSetups/FallingSpheres_Multigrid.dat -gmg_mg_levels_ksp_type richardson -gmg_mg_levels_pc_type jacobi -gmg_mg_levels_ksp_richardson_scale 0.5') 
  expected_file = makeLocalPathAbsolute('FallingSpheres_Multigrid-Jacobi-p2.expected')

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = 'CONVERGED_RTOL iterations'
    unittest.compareInteger(key,0)

    key = 'KSP Residual norm'
    unittest.compareFloatingPoint(key,1e-6)
    
  # Create unit test object
  test = pth.pthUnitTest('FallingSpheres_Multigrid_Jacobi',ranks,launch,expected_file)
  test.setVerifyMethod(comparefunc)
  test.appendKeywords('@')
  test.setUseSandbox()      # put test output in seperate directory

  return(test)


#------------------------------------------------------------------------------------------------
def test_8():
  # 2D viscoplastic free slip subduction setup with direct solvers
  ranks = 2
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ../BuildInSetups/Subduction2D_FreeSlip_DirectSolver.dat -dt_out 0 -nstep_ini 0 -nel_x 64 -nel_z 16 -nstep_max 3 -rand_noise 0') 
  expected_file = makeLocalPathAbsolute('Subduction2D_FreeSlip_DirectSolver-p2.expected')

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = 'CONVERGED_RTOL iterations'
    unittest.compareInteger(key,0)

    key = 'KSP Residual norm'
    unittest.compareFloatingPoint(key,1e-6)
    
  # Create unit test object
  test = pth.pthUnitTest('Subduction2D_FreeSlip_DirectSolver',ranks,launch,expected_file)
  test.setVerifyMethod(comparefunc)
  test.appendKeywords('@')
  test.setUseSandbox()      # put test output in seperate directory

  return(test)


#------------------------------------------------------------------------------------------------
def test_9():
  # 2D viscous free surface subduction setup with direct solvers
  ranks = 1
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ../BuildInSetups/Subduction2D_FreeSurface_DirectSolver.dat -dt_out 0 -nstep_ini 0 -nel_x 128 -nel_z 32 -nstep_max 3 -rand_noise 0') 
  expected_file = makeLocalPathAbsolute('Subduction2D_FreeSurface_DirectSolver-p1.expected')

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = 'CONVERGED_RTOL iterations'
    unittest.compareInteger(key,0)

    key = 'KSP Residual norm'
    unittest.compareFloatingPoint(key,1e-6)
    
  # Create unit test object
  test = pth.pthUnitTest('Subduction2D_FreeSurface_DirectSolver',ranks,launch,expected_file)
  test.setVerifyMethod(comparefunc)
  test.appendKeywords('@')
  test.setUseSandbox()      # put test output in seperate directory

  return(test)


#------------------------------------------------------------------------------------------------
def test_10():
  # 3D viscous double subduction setup with multigrid solver
  ranks = 2
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ../BuildInSetups/Subduction3D_DoubleSubduction_FreeSlip_Multigrid.dat -dt_out 0 -nstep_ini 0 -nel_x 64 -nel_y 32 -nel_z 16 -nstep_max 3 -rand_noise 0') 
  expected_file = makeLocalPathAbsolute('Subduction3D_FreeSlip_MultigridSolver-p2.expected')

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = 'CONVERGED_RTOL iterations'
    unittest.compareInteger(key,0)

    key = 'KSP Residual norm'
    unittest.compareFloatingPoint(key,1e-6)
    
  # Create unit test object
  test = pth.pthUnitTest('Subduction3D_FreeSlip_MultigridSolver',ranks,launch,expected_file)
  test.setVerifyMethod(comparefunc)
  test.appendKeywords('@')
  test.setUseSandbox()      # put test output in seperate directory

  return(test)


#------------------------------------------------------------------------------------------------
def test_11():
  # 2D asymmetric rifting setup
  ranks = 2
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ../BuildInSetups/Rifting2D_MultigridSolver.dat -dt_out 0 -nstep_ini 0 -nel_x 64 -nel_z 32 -nstep_max 25 -rand_noise 0') 
  expected_file = makeLocalPathAbsolute('Rifting2D_VEP_MultigridSolver-p2.expected')

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = 'CONVERGED_RTOL iterations'
    unittest.compareInteger(key,0)

    key = 'KSP Residual norm'
    unittest.compareFloatingPoint(key,1e-6)
    
    key = 'SNES Function norm'
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  test = pth.pthUnitTest('Rifting2D_VEP_MultigridSolver',ranks,launch,expected_file)
  test.setVerifyMethod(comparefunc)
  test.appendKeywords('@')
  test.setUseSandbox()      # put test output in seperate directory

  return(test)

  
  #------------------------------------------------------------------------------------------------
def test_12():
  # 2D plume-lithosphere interaction
  ranks = 2
  launch =  makeLocalPathAbsolute('../../../bin/opt/LaMEM -ParamFile ../BuildInSetups/PlumeLithosphereInteraction.dat -dt_out 0 -nstep_ini 0 -nel_x 32 -nel_z 32 -nstep_max 3 -rand_noise 0') 
  expected_file = makeLocalPathAbsolute('PlumeLithosphereInteraction-p2.expected')

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-6)

    key = 'CONVERGED_RTOL iterations'
    unittest.compareInteger(key,0)

    key = 'KSP Residual norm'
    unittest.compareFloatingPoint(key,1e-6)
    
    key = 'SNES Function norm'
    unittest.compareFloatingPoint(key,1e-4)

  # Create unit test object
  test = pth.pthUnitTest('PlumeLithosphereInteraction_Direct',ranks,launch,expected_file)
  test.setVerifyMethod(comparefunc)
  test.appendKeywords('@')
  test.setUseSandbox()      # put test output in seperate directory

  return(test)