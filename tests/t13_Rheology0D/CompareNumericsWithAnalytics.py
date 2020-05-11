
import os
import re
import subprocess


#----------------------------------------------------
# Read data from disk
def LoadTimeDependentData(fname):

  # load required modules for reading data
  try: 
    from vtk import vtkXMLPRectilinearGridReader
    from vtk.util import numpy_support as VN
    import numpy as np
    from pathlib import Path
  except:
    print('VTK toolboxes are not installed; cannot load data')


  # extract directory name and times
  time_vec  = np.zeros(100)
  dir_vec   = np.empty(21,dtype='object')
  num       = 0
  fu        = [f.path for f in os.scandir('.') if f.is_dir()];
  
  for dirname in fu:
    if (dirname.startswith('./Timestep') ):
      number            = int(dirname[11:19]);
      time_vec[number]  = float(dirname[20:]);    # Time 
      dir_vec[number]   = dirname
      num               = num+1;                  # number
      #print(dirname)
      #print(number)
      #print(float(dirname[18:]))

  
  #Initialize values
  T2nd = np.empty(num);
  P    = np.empty(num);
  T    = np.empty(num);
  E2nd = np.empty(num);
  Time = np.empty(num);

  # Loop over all timestep directories
  for i in range(num):
    dir = dir_vec[i]
    
    file_in =  fname+'.pvtr'
    filesep = '/'
    filename = dir+filesep+file_in
    #print ('Loading LaMEM output from directory ', filename)

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
    T2nd_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('j2_dev_stress [MPa]')); 
    E2nd_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('j2_strain_rate [1/s]')); 
    Pres_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('pressure [MPa]')); 
    Temp_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('temperature [C]')); 

    # Compute average value of T2nd 
    T2nd[i] = np.average(T2nd_v);
    E2nd[i] = np.average(E2nd_v);
    P[i]    = np.average(Pres_v);
    T[i]    = np.average(Temp_v);
    Time[i] = time_vec[i];


  # save data
  data.T2nd = T2nd;
  data.P    = P;
  data.T    = T;
  data.E2nd = E2nd;
  data.Time = Time;

  return(data);

#----------------------------------------------------

#----------------------------------------------------
# Read data from disk
def LoadStrainrateData(fname):

  # load required modules for reading data
  try: 
    from vtk import vtkXMLPRectilinearGridReader
    from vtk.util import numpy_support as VN
    import numpy as np
    from pathlib import Path
  except:
    print('VTK toolboxes are not installed; cannot load data')


  # extract directory name and times
  E2nd_vec  = np.zeros(100)
  dir_vec   = np.empty(5,dtype='object')
  num       = 0
  fu        = [f.path for f in os.scandir('.') if f.is_dir()];
  

  for dirname in fu:
    if (dirname.startswith('./Strainrate_') ):
      number            = int(dirname[13]);
      E2nd_vec[number]  = 10**(-float(dirname[15:]))
      dir_vec[number]   = dirname
      num               = num+1;                  # number
     

  #Initialize values
  T2nd    = np.empty(num);
  P       = np.empty(num);
  T       = np.empty(num);
  E2nd    = np.empty(num);
  E2nd_bg = np.empty(num);

  # Loop over all timestep directories
  for i in range(num):
    dir     =   dir_vec[i]
    file_in =   fname+'.pvtr'
    filesep =   '/'

    filename = dir+filesep+file_in
    #print ('Loading LaMEM output from directory ', filename)

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
    T2nd_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('j2_dev_stress [MPa]')); 
    E2nd_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('j2_strain_rate [1/s]')); 
    Pres_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('pressure [MPa]')); 
    Temp_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('temperature [C]')); 

    # Compute average value of T2nd 
    T2nd[i]       = np.average(T2nd_v);
    E2nd[i]       = np.average(E2nd_v);
    P[i]          = np.average(Pres_v);
    T[i]          = np.average(Temp_v);
    E2nd_bg[i]    = E2nd_vec[i];


  # save data for the different strainrate components
  data.T2nd     = T2nd;
  data.P        = P;
  data.T        = T;
  data.E2nd     = E2nd;
  data.E2nd_bg  = E2nd_bg;

  return(data);

#----------------------------------------------------


#----------------------------------------------------
# Compute the analytical solution for this problem
def AnalyticalSolution_VE(data, YieldStress):
   # load required modules for reading data
  try: 
    import numpy as np
  except:
    print('numpy toolbox not installed; cannot manipulate data')


  # Define material properties for each of the phases 
  #  This should obviously be the same as in the input file
  G             = 5e10;             # elastic shear module
  eta           = 1e22;             # viscosity
  str           = 1e-15;            # BG strainrate 
  SecYear       = 3600*24*365.25;   # sec/year
  time_end      = 0.05*1e6*SecYear; 
  num           = 500;
  dt            = time_end/num;

  tau_maxwell   = eta/G;
  
  Time_anal     = np.arange(0, time_end, dt)
  T2nd_anal     = np.arange(0, time_end, dt)

  # Create arrays with material constants
  num=0;
  for time in Time_anal:
    T2nd_anal[num]  = 2*eta*(1-np.exp(-time/tau_maxwell))*str;
    if (T2nd_anal[num]>YieldStress):
      T2nd_anal[num]=YieldStress;
    num             = num+1;

  # Compute analytical solution @ same points as LaMEM solution
  T2nd_LaMEM_anal  = np.zeros(len(data.Time))
  RMS_err_Numerics = np.zeros(len(data.Time))
  num=0;
  for time in data.Time:
    T2nd_LaMEM_anal[num]  = 2*eta*(1-np.exp(-time*SecYear*1e6/tau_maxwell))*str;    # in MPa
    if (T2nd_LaMEM_anal[num]>YieldStress):
      T2nd_LaMEM_anal[num]=YieldStress;

    RMS_err_Numerics[num] = (T2nd_LaMEM_anal[num]/1e6-data.T2nd[num]) #error in MPa 
    num                   = num+1;

  Error_L2 = np.linalg.norm(RMS_err_Numerics);
  print('Test: L2 error-norm [MPa] = ',Error_L2)

  # store data
  data.T2nd_LaMEM_anal    =   T2nd_LaMEM_anal/1e6;
  data.T2nd_anal          =   T2nd_anal/1e6;
  data.Time_anal          =   Time_anal/SecYear/1e6;
  
  return(data);
#----------------------------------------------------


#----------------------------------------------------
# Compute the analytical solution for this problem
def AnalyticalSolution_linearViscous(data):
   # load required modules for reading data
  try: 
    import numpy as np
  except:
    print('numpy toolbox not installed; cannot manipulate data')


  # Define material properties  
  eta           = 1e21;             # viscosity

  E2nd_anal     = 10**np.arange(np.min(np.log10(data.E2nd_bg)), np.max(np.log10(data.E2nd_bg)), .1)
  T2nd_anal     = 10**np.arange(np.min(np.log10(data.E2nd_bg)), np.max(np.log10(data.E2nd_bg)), .1)

  # Create arrays with material constants
  num=0;
  for E2nd in E2nd_anal:
    T2nd_anal[num]  = 2*eta*E2nd;
    num             = num+1;

  # Compute analytical solution @ same points as LaMEM solution
  T2nd_LaMEM_anal  = np.zeros(len(data.E2nd_bg))
  RMS_err_Numerics = np.zeros(len(data.E2nd_bg))
  num=0;
  for E2nd in data.E2nd_bg:
    T2nd_LaMEM_anal[num]  = 2*eta*E2nd;    # in MPa
    RMS_err_Numerics[num] = (T2nd_LaMEM_anal[num]/1e6-data.T2nd[num]) #error in MPa 
    num                   = num+1;

  Error_L2 = np.linalg.norm(RMS_err_Numerics);
  print('Test: L2 error-norm [MPa] = ',Error_L2)

  # store data
  data.T2nd_LaMEM_anal    =   T2nd_LaMEM_anal/1e6;
  data.T2nd_anal          =   T2nd_anal/1e6;
  data.E2nd_anal          =   E2nd_anal;
  
  return(data);
#----------------------------------------------------

#----------------------------------------------------
# Compute the analytical solution for this problem
def AnalyticalSolution_DislocationCreep(data, FlowLaw):
   # load required modules for reading data
  try: 
    import numpy as np
  except:
    print('numpy toolbox not installed; cannot manipulate data')


  # Define material properties  
  if FlowLaw.lower() in ['dryolivine']:
      AD    = 2.5e-17;   # 1/Pa^n/s, 
      n     = 3.5;       # dimensionless
      Ea    = 532000;    # J/mol 
      R     = 8.314;     # J/mol/K

      # Define correction coefficient F2 
      # for a strain rate based viscosity formulation
      F2    = 1/2**((n-1)/n)/3**((n+1)/2/n);


  T       = data.T[0];      # average temperature in domain (in Celcius)

  eterm   = np.exp(Ea/n/R/(T+273));


  E2nd_anal     = 10**np.arange(np.min(np.log10(data.E2nd_bg)), np.max(np.log10(data.E2nd_bg)), .1)
  T2nd_anal     = 10**np.arange(np.min(np.log10(data.E2nd_bg)), np.max(np.log10(data.E2nd_bg)), .1)

  # Create arrays with material constants
  num=0;
  for E2nd in E2nd_anal:
    
    eta_eff         = F2/AD**(1/n)/E2nd**((n-1)/n)*eterm;     # effective viscosity

    T2nd_anal[num]  = 2*eta_eff*E2nd;
    num             = num+1;

  # Compute analytical solution @ same points as LaMEM solution
  T2nd_LaMEM_anal  = np.zeros(len(data.E2nd_bg))
  RMS_err_Numerics = np.zeros(len(data.E2nd_bg))
  num=0;
  for E2nd in data.E2nd_bg:
    eta_eff               = F2/AD**(1/n)/E2nd**((n-1)/n)*eterm;    # effectic viscosity

    T2nd_LaMEM_anal[num]  = 2*eta_eff*E2nd;    # in MPa
    RMS_err_Numerics[num] = (T2nd_LaMEM_anal[num]/1e6-data.T2nd[num]) #error in MPa 
    num                   = num+1;

  Error_L2 = np.linalg.norm(RMS_err_Numerics);
  print('Test: L2 error-norm [MPa] = ',Error_L2)

  # store data
  data.T2nd_LaMEM_anal    =   T2nd_LaMEM_anal/1e6;
  data.T2nd_anal          =   T2nd_anal/1e6;
  data.E2nd_anal          =   E2nd_anal;
  
  return(data);
#----------------------------------------------------

#----------------------------------------------------
# Plot time-dependent data of LaMEM
def PlotTimeDependentData(data,plotName):
  # load required modules for reading data
  try: 
    import matplotlib.pyplot as plt
    import datetime 
  except:
    print('matplotlib toolbox not installed; cannot plot data')

  # Time and date stamp on title
  title_str = datetime.datetime.now().strftime("Test performed %H:%M on %B %d, %Y");

  fig1 = plt.figure()
  ax1 = fig1.add_subplot(111)
  ax1.plot(data.Time_anal,  data.T2nd_anal,'b-');   # Stress analytical
  ax1.plot(data.Time,       data.T2nd,'r+');        # LaMEM
  
  ax1.legend(['analytical','numerical'],fontsize='x-small')  
  ax1.set_ylabel('T2nd [MPa]');
  ax1.set_xlabel('Time [Myrs]');
  ax1.axis([0, 0.05, 0, 21])

  title_str = 'LaMEM vs. analytics \n ' + title_str;
  plt.title(title_str,fontsize='xx-small')

  fig1.savefig(plotName,dpi=360)    # save figure to disk

  return();

#----------------------------------------------------


#----------------------------------------------------
# Plot data from disk
def PlotStrainrateData(data,plotName):
  # load required modules for reading data
  try: 
    import matplotlib.pyplot as plt
    import datetime 
  except:
    print('matplotlib toolbox not installed; cannot plot data')

  # Time and date stamp on title
  title_str = datetime.datetime.now().strftime("Test performed %H:%M on %B %d, %Y");

  fig1 = plt.figure()
  ax1 = fig1.add_subplot(111)
  ax1.loglog(data.E2nd_anal,  data.T2nd_anal,'b-');   # Stress analytical
  ax1.loglog(data.E2nd_bg,    data.T2nd,'r+');        # LaMEM
  
  ax1.legend(['analytical','numerical'],fontsize='x-small')  
  ax1.set_ylabel('T2nd [MPa]');
  ax1.set_xlabel('E2nd [1/s]');

  title_str = 'LaMEM vs. analytics \n ' + title_str;
  plt.title(title_str,fontsize='xx-small')

  fig1.savefig(plotName,dpi=360)    # save figure to disk

  return();

#----------------------------------------------------




#data = LoadData();    # Load the data using the VTK toolbox
#data = AnalyticalSolution_VE(data);                        # Compute analytical solution
#PlotData(data,'t13_ViscoElastic_output.png');             # Create Plot
