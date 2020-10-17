import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re
import subprocess



#----------------------------------------------------
# Read data from disk
def LoadData(path):
  # load required modules for reading data
  try: 
    from vtk import vtkXMLPRectilinearGridReader
    from vtk.util import numpy_support as VN
    import numpy as np
    import glob
    import shutil
  except:
    print('VTK toolboxes are not installed; cannot load data')

  # get list of timesteps
  dirlist  = glob.glob(os.path.join(path,'Timestep*'))
  # sort Tiemsteps (last one first)
  dirlist.sort(reverse=True)
  # get last timestep
  step = dirlist[0]

  file_in =  'output.pvtr'
  filesep = '/'
  filename = step+filesep+file_in
  reader = vtkXMLPRectilinearGridReader()
  reader.SetFileName(filename)
  reader.Update()
  data = reader.GetOutput()

  # extract coordinates
  data.x = VN.vtk_to_numpy(data.GetXCoordinates())
  data.y = VN.vtk_to_numpy(data.GetYCoordinates())
  data.z = VN.vtk_to_numpy(data.GetZCoordinates())

  # dimensions
  nx = np.size(data.x); ny = np.size(data.y); nz = np.size(data.z)

  # load data from VTK file & reshape them to 3D arrays
  P             =   VN.vtk_to_numpy(data.GetPointData().GetArray('pressure [MPa]'));      data.P      = np.reshape(P  ,(nz,ny,nx));
  rho           =   VN.vtk_to_numpy(data.GetPointData().GetArray('density [kg/m^3]'));    data.rho    = np.reshape(rho,(nz,ny,nx));
  Pf            =   VN.vtk_to_numpy(data.GetPointData().GetArray('pore_press [MPa]'));    data.Pf     = np.reshape(Pf ,(nz,ny,nx));
  Pl            =   VN.vtk_to_numpy(data.GetPointData().GetArray('litho_press [MPa]'));   data.Pl     = np.reshape(Pl ,(nz,ny,nx));
  T2nd_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('j2_dev_stress [MPa]')); data.T2nd   = np.reshape(T2nd_v ,(nz,ny,nx));
  Phase_v       =   VN.vtk_to_numpy(data.GetPointData().GetArray('phase [ ]'));           data.phase  = np.reshape(Phase_v ,(nz,ny,nx)); 
  SHmaxv        =   VN.vtk_to_numpy(data.GetPointData().GetArray('SHmax [ ]'));           
  data.SHmax_x  =   np.reshape(SHmaxv[:,0],(nz,ny,nx)); data.SHmax_y  =   np.reshape(SHmaxv[:,1],(nz,ny,nx)); data.SHmax_z  =   np.reshape(SHmaxv[:,2],(nz,ny,nx));
  T             =   VN.vtk_to_numpy(data.GetPointData().GetArray('dev_stress [MPa]'));    
  data.Txx      =   np.reshape(T[:,0],(nz,ny,nx));     data.Txy = np.reshape(T[:,1],(nz,ny,nx)); data.Txz = np.reshape(T[:,2],(nz,ny,nx));
  data.Tyy      =   np.reshape(T[:,4],(nz,ny,nx));     data.Tyz = np.reshape(T[:,5],(nz,ny,nx)); data.Tzz = np.reshape(T[:,8],(nz,ny,nx));
  
  # Compute required stresses, etc.
  data.Sxx  = -(-(data.P) + data.Txx)        # change sign, to be consistent with analytics later
  data.Syy  = -(-(data.P) + data.Tyy)       
  data.Szz  = -(-(data.P) + data.Tzz)
  
  # clean up
  for dir in dirlist:
    shutil.rmtree(dir)

  return(data)

#----------------------------------------------------

#----------------------------------------------------
# Compute the analytical solution for this problem
def AnalyticalSolution(data):
   # load required modules for reading data
  try: 
    import numpy as np
  except:
    print('numpy toolbox not installed; cannot manipulate data')


  # Define material properties for each of the phases 
  #  This should obviously be the same as in the input file
  biot    = 1.0
  rho_w   = 1000                                    # density of water
  p_fac   = np.array([0, 0.1, 0.1, 0.3])            # pore fluid factor 
  rho     = data.rho[:,1,1]                         # Density
  v_vec   = np.array([0.4999, 0.27, 0.4999, 0.27])  # Poissons Ratio

  phase_vec     = data.phase[:,1,1]
  nz            = np.size(phase_vec)
  rp            = np.empty_like(phase_vec)
  v             = np.empty_like(phase_vec)
  
  Sv_anal       = np.empty_like(phase_vec)  # vertical stress
  P_hydro       = np.empty_like(phase_vec)  # 
  Pf_anal       = np.empty_like(phase_vec)  # pore fluid

  # Create arrays with material constants
  for i in range(nz-1,-1,-1):
    phase     = int( phase_vec[i] )
    v[i]      =   v_vec[int( phase_vec[i] )]
    rp[i]     =   p_fac[int( phase_vec[i] )]
    P_hydro[i]=  -rho_w*10*data.z[i]*1e3    

  P_hydro[P_hydro<0] = 0

  # Compute vertical stress 
  Sv_anal[nz-1] = 0
  for i in range(nz-2,-1,-1):
    rho_mean = (rho[i+1] + rho[i])/2.0
    dz       = (data.z[i+1] - data.z[i])*1000    # in m
    Sv_anal[i] = Sv_anal[i+1] + dz*10*rho_mean

  Pf_anal = P_hydro + rp*(Sv_anal-P_hydro) 
  
  Sh_anal = (v/(1-v))*(Sv_anal-biot*Pf_anal) + biot*Pf_anal

  # store data
  data.Sv_anal    =   Sv_anal/1e6
  data.Pf_anal    =   Pf_anal/1e6
  data.P_hydro    =   P_hydro/1e6
  data.Sh_anal    =   Sh_anal/1e6
  
  return(data)
#----------------------------------------------------

#----------------------------------------------------
# Plot data from disk
def PlotData(data,plotName):
  # load required modules for reading data
  try: 
    import matplotlib.pyplot as plt
    matplotlib.use('Agg')
    import datetime 
    
  except:
    print('matplotlib toolbox not installed; cannot plot data')

  # Time and date stamp on title
  title_str = datetime.datetime.now().strftime("Test performed %H:%M on %B %d, %Y")

  fig1 = plt.figure()


  ax1 = fig1.add_subplot(121)
  ax1.plot(data.Sv_anal*10, data.z,'r-')   # Vertical stress analytically
  ax1.plot(data.Pf_anal*10, data.z,'b-')   # Pore fluid pressure analytically
  ax1.plot(data.P_hydro*10, data.z,'b--')  # Hydrostatic Pore fluid pressure
  ax1.plot(data.Sh_anal*10, data.z,'g-')   # Horizontal stress analytically
  
  ax1.plot(data.Pf[:,1,1]*10, data.z,'b.')    # pore pressure LaMEM
  ax1.plot(data.Szz[:,1,1]*10, data.z,'r+')   # Vertical stress LaMEM
  ax1.plot(data.Sxx[:,1,1]*10, data.z,'g+')   # Horizonal stress LaMEM

  ax1.legend(['$S_v$ analytical','$P_{fluid}$ analytical','$P_{hydrostatic}$','$S_h$ analytical',
              '$S_v$ LaMEM','$P_{fluid}$ LaMEM','$S_h$ LaMEM'            ],fontsize='x-small')  
  ax1.set_ylabel('Depth [km]')
  ax1.set_xlabel('Pressure & Stress [bar]')
  ax1.axis([0, 1250, -5, 1.0])

  title_str = 'Elastic compressible setup : LaMEM vs. analytics \n ' + title_str
  plt.title(title_str,fontsize='xx-small')

  ax2 = fig1.add_subplot(122)
  ax2.plot(data.T2nd[:,1,1], data.z,'k+-')     # 2nd invariant of deviatoric stress tensor, LaMEM
  ax2.axis([0, 20, -5, 1.0])

  ax2.set_xlabel('T2nd [MPa]')

  ax3 = ax2.twiny() 

  #ax3 = fig1.add_subplot(133)
  ax3.plot(data.rho[:,1,1], data.z,'r-+')     # density (should be depth-depend)
  ax3.axis([0, 2500, -5, 1.0])

  ax3.set_xlabel('Density [kg/m3]',color='r')
  ax3.tick_params(axis='x', labelcolor='r')


  fig1.savefig(plotName,dpi=360)    # save figure to disk

  return()

#----------------------------------------------------



# This tests whether elastic compressibility works, and creates plots if possible
def test_a():
# Tests are performed on 1 and 2 cores, using opt/deb

  
  #==============================================
  # Run the input script wth matlab-generated particles
  ranks = 1
  launch = ['rm -r Timestep*',
            'rm -rf ./t10_Compressibility/Out1Core; mkdir ./t10_Compressibility/Out1Core',
            '../bin/opt/LaMEM -ParamFile ./t10_Compressibility/Compressible1D_withSaltandBasement.dat',
            'mv Timestep* ./t10_Compressibility/Out1Core'] # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't10_Compressibility/test_10_Compressibility_opt-p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

    #----------------------------  
    try: 
      data = LoadData('./t10_Compressibility/Out1Core')      # Load the data using the VTK toolbox
      data = AnalyticalSolution(data)                        # Compute analytical solution
      PlotData(data,'./t10_Compressibility/Compressible1D_output.png')            # Create Plot

      print('Created output figure ./t10_Compressibility/Compressible1D_output.png comparing analytics vs. numerics')
    except:
      print('VTK/MatPlotLib/NumPy toolboxes are not installed; will not create plots')
    #----------------------------

  # Create unit test object
  ex1 = pth.pthUnitTest('t10_Compressibility_Direct_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)


def test_b():
  # Test visco-elasto-plastic localization case on 2 cores, using debug LaMEM
  ranks = 2
  launch = ['rm -r Timestep* 2>/dev/null',
            'rm -rf ./t10_Compressibility/Out2Core; mkdir ./t10_Compressibility/Out2Core',
            '../bin/deb/LaMEM -ParamFile ./t10_Compressibility/Compressible1D_withSaltandBasement.dat',
            'mv Timestep* ./t10_Compressibility/Out2Core 2>/dev/null'] # This must be a relative path with respect to runLaMEM_Tests.py
  expected_file = 't10_Compressibility/Compressibility_Direct_deb-p2.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

    #----------------------------  
    try: 
      data = LoadData('./t10_Compressibility/Out2Core')      # Load the data using the VTK toolbox
      data = AnalyticalSolution(data)                        # Compute analytical solution
      PlotData(data,'./t10_Compressibility/Compressible1D_2Cores_output.png')             # Create Plot
    
      print('Created output figure ./t10_Compressibility/Compressible1D_2Cores_output.png comparing analytics vs. numerics')
    except:
      print('VTK/MatPlotLib/NumPy toolboxes are not installed; will not create plots')
    #----------------------------

  # Create unit test object
  ex1 = pth.pthUnitTest('t10_Compressibility_Direct_deb',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')

  return(ex1)



