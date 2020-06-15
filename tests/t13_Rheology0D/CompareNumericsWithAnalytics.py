import matplotlib
matplotlib.use('Agg')

import os
import re
import subprocess

#----------------------------------------------------
# Read data from disk
def LoadTimeDependentData(path,fname):
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
  # sort Tiemsteps
  dirlist.sort()
  numSteps = len(dirlist)
  
  #Initialize values
  T2nd = np.empty(numSteps)
  P    = np.empty(numSteps)
  T    = np.empty(numSteps)
  E2nd = np.empty(numSteps)
  Time = np.empty(numSteps)

  # Loop over all timestep directories
  for i in range(numSteps):
    dir = dirlist[i]
    
    file_in =  fname+'.pvtr'
    filesep = '/'
    filename = dir+filesep+file_in

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

    # Read Time
    stepName      = dir.split('/')
    stepName      = stepName[3]
    time          = stepName.split('_')
    t             = float(time[2])
    Time[i]       = t

    # load data from VTK file & reshape them to 3D arrays
    T2nd_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('j2_dev_stress [MPa]'))
    E2nd_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('j2_strain_rate [1/s]'))
    Pres_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('pressure [MPa]'))
    Temp_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('temperature [C]'))

    # Compute average value of T2nd 
    T2nd[i]       = np.average(T2nd_v)
    E2nd[i]       = np.average(E2nd_v)
    P[i]          = np.average(Pres_v)
    T[i]          = np.average(Temp_v)


  # save data
  data.T2nd = T2nd
  data.P    = P
  data.T    = T
  data.E2nd = E2nd
  data.Time = Time

  # remove output
  for dir in dirlist:
    shutil.rmtree(dir)

  return(data)

#----------------------------------------------------

#----------------------------------------------------
# Read data from disk
def LoadStrainrateData(path,fname):
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
  dirlist  = glob.glob(os.path.join(path,'Strainrate*'))
  # sort Timesteps
  dirlist.sort()
  numSteps = len(dirlist)   

  #Initialize values
  T2nd    = np.empty(numSteps)
  P       = np.empty(numSteps)
  T       = np.empty(numSteps)
  E2nd    = np.empty(numSteps)
  E2nd_bg = np.empty(numSteps)
  Eta_eff = np.empty(numSteps)
  

  # Loop over all timestep directories
  for i in range(numSteps):
    dir     =   dirlist[i]
    file_in =   fname+'.pvtr'
    filesep =   '/'

    filename = dir+filesep+file_in

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

    # Extract exponent
    stepName      = dir.split('/')
    stepName      = stepName[3]
    expo          = stepName.split('_')
    expo          = float(expo[2])
    E2nd_bg[i]    = 10**(-float(expo))

    # load data from VTK file & reshape them to 3D arrays
    T2nd_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('j2_dev_stress [MPa]'))
    E2nd_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('j2_strain_rate [1/s]'))
    Pres_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('pressure [MPa]'))
    Temp_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('temperature [C]'))

    # Compute average value of T2nd 
    T2nd[i]       = np.average(T2nd_v)
    E2nd[i]       = np.average(E2nd_v)
    P[i]          = np.average(Pres_v)
    T[i]          = np.average(Temp_v)
    Eta_eff[i]    = T2nd[i]*1e6/2/(E2nd[i])

  # save data for the different strainrate components
  data.T2nd       = T2nd
  data.P          = P
  data.T          = T
  data.E2nd       = E2nd
  data.E2nd_bg    = E2nd_bg
  data.Eta_eff    = Eta_eff

  # remove output
  for dir in dirlist:
    shutil.rmtree(dir)

  return(data)

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
  Eta_eff_anal  = np.arange(0, time_end, dt)

  # Create arrays with material constants
  num=0;
  for time in Time_anal:
    T2nd_anal[num]  = 2*eta*(1-np.exp(-time/tau_maxwell))*str;
    if (T2nd_anal[num]>YieldStress):
      T2nd_anal[num]=YieldStress;
    Eta_eff_anal[num] = T2nd_anal[num]/2/str;

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
  data.Eta_eff_anal       =   Eta_eff_anal;

  return(data);
#----------------------------------------------------

#----------------------------------------------------
# Compute the analytical solution for this problem
def AnalyticalSolution_DislocationCreep_VEP(data, YieldStress):
  # Note that this requires local iterations, as only the viscous strain rate should be employed to compute the creeplaw
  # Below we create a comparison of a case w/out and one with local iterations 
   
  # load required modules for reading data
  try: 
    import numpy as np
  except:
    print('numpy toolbox not installed; cannot manipulate data')


  # Define material properties for each of the phases 
  #  This should obviously be the same as in the input file
  G             = 5e10;             # elastic shear module
  
  # these are the parameters for Dry Olivine, to be consistent with Gerya's book; can be changed ofcourse
  AD    = 2.5e-17;   # 1/Pa^n/s, 
  n     = 3.5;       # dimensionless
  Ea    = 532000;    # J/mol 
  R     = 8.3144621;     # J/mol/K

  # Define correction coefficient F2 
  # for a strain rate based viscosity formulation
  F2    = 1/2**((n-1)/n)/3**((n+1)/2/n);

  
  str           = 1e-14;            # BG strainrate 
  SecYear       = 3600*24*365.25;   # sec/year
  time_end      = 0.005*1e6*SecYear; 
  num           = 500;

  T             = data.T[0];      # average temperature in domain (in Celcius)
  eterm         = np.exp(Ea/n/R/(T+273.15));
    
  eta_eff       = F2/AD**(1/n)/str**((n-1)/n)*eterm;     # effective viscosity

  dt            = time_end/num;

  tau_maxwell   = eta_eff/G;
  
  Time_anal         = np.arange(0, time_end, dt)
  T2nd_anal         = np.arange(0, time_end, dt)
  T2nd_anal_noLocal = np.arange(0, time_end, dt)
  Eta_eff_anal      = np.arange(0, time_end, dt)

  # Create arrays with material constants
  num=0;
  for time in Time_anal:
    T2nd_anal_noLocal[num]  = 2*eta_eff*(1-np.exp(-time/tau_maxwell))*str;    # this does NOT employ local iterations
    if (T2nd_anal_noLocal[num]>YieldStress):
      T2nd_anal_noLocal[num]=YieldStress;
    
    num             = num+1;

  # Perform local nonlinear iterations (correct way)
  T2nd_anal[0] = 0;
  e_vis        = str/2;
  for num in range(1, len(Time_anal)):
    Tau     =   T2nd_anal[num-1];
    Tau_new =   Tau;
    dTau    =   1000;
    it      =   1;
    for it in range(50):
   
      Tau_new1    =   Tau + 2*G*dt*(str-e_vis);               # update stress
      eta         =   F2/AD**(1/n)/e_vis**((n-1)/n)*eterm;    # local iterations; employ viscous strain rate here (the wrong approach is to use the full strainrate)
      e_vis       =   Tau_new1/2/eta;                          # viscous strainrate

      if (eta>1e28):
          eta=1e28;
      
      if (Tau_new1>YieldStress):
        Tau_new1 = YieldStress;

      dTau        = Tau_new1-Tau_new;
      Tau_new     = Tau_new1;
      it          = it+1;
    
    T2nd_anal[num]  =   Tau_new;
    
    Eta_eff_anal[num] = T2nd_anal[num]/2/str;
    num             = num+1;

  # Interpolate analytical solution to LaMEM datapoints (linear interpolation)
  T2nd_LaMEM_anal   = np.interp(data.Time, Time_anal/SecYear/1e6, T2nd_anal);
  RMS_err_Numerics  = T2nd_LaMEM_anal/1e6-data.T2nd;                            # error

  # Print error (if required)
  num=0;
  for time in data.Time:
  #  print('T2nd: [LaMEM;Anal;Ratio] = [',data.T2nd[num],';',T2nd_LaMEM_anal[num]/1e6,';',data.T2nd[num]/(T2nd_LaMEM_anal[num]/1e6),']')
    num                   = num+1;



  Error_L2 = np.linalg.norm(RMS_err_Numerics);
  print('Test: L2 error-norm [MPa] = ',Error_L2)

  # store data
  data.T2nd_LaMEM_anal    =   T2nd_LaMEM_anal/1e6;
  data.T2nd_anal_noLocal  =   T2nd_anal_noLocal/1e6;
  data.T2nd_anal          =   T2nd_anal/1e6;
  data.Eta_eff_anal       =   Eta_eff_anal;
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
  Eta_eff_anal  = 10**np.arange(np.min(np.log10(data.E2nd_bg)), np.max(np.log10(data.E2nd_bg)), .1)

  # Create arrays with material constants
  num=0;
  for E2nd in E2nd_anal:
    T2nd_anal[num]  = 2*eta*E2nd;
    Eta_eff_anal[num] = T2nd_anal[num]/2/E2nd;

    num             = num+1;

  # Compute analytical solution @ same points as LaMEM solution
  T2nd_LaMEM_anal  = np.zeros(len(data.E2nd_bg))
  RMS_err_Numerics = np.zeros(len(data.E2nd_bg))
  num=0;
  for E2nd in data.E2nd_bg:
    T2nd_LaMEM_anal[num]  = 2*eta*E2nd;    # in MPa
    RMS_err_Numerics[num] = (T2nd_LaMEM_anal[num]/1e6-data.T2nd[num]) #error in MPa 

    #print('T2nd: [LaMEM;Anal;Ratio] = [',data.T2nd[num],';',T2nd_LaMEM_anal[num]/1e6,';',data.T2nd[num]/(T2nd_LaMEM_anal[num]/1e6),']')

    num                   = num+1;

  Error_L2 = np.linalg.norm(RMS_err_Numerics);
  print('Test: L2 error-norm [MPa] = ',Error_L2)

  # store data
  data.T2nd_LaMEM_anal    =   T2nd_LaMEM_anal/1e6;
  data.T2nd_anal          =   T2nd_anal/1e6;
  data.E2nd_anal          =   E2nd_anal;
  data.Eta_eff_anal       =   Eta_eff_anal;
  
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
      R     = 8.3144621;     # J/mol/K

      # Define correction coefficient F2 
      # for a strain rate based viscosity formulation
      F2    = 1/2**((n-1)/n)/3**((n+1)/2/n);


  T       = data.T[0];      # average temperature in domain (in Celcius)

  eterm   = np.exp(Ea/n/R/(T+273.15));


  E2nd_anal     = 10**np.arange(np.min(np.log10(data.E2nd_bg)), np.max(np.log10(data.E2nd_bg)), .1)
  T2nd_anal     = 10**np.arange(np.min(np.log10(data.E2nd_bg)), np.max(np.log10(data.E2nd_bg)), .1)
  Eta_eff_anal  = 10**np.arange(np.min(np.log10(data.E2nd_bg)), np.max(np.log10(data.E2nd_bg)), .1)

  # Create arrays with material constants
  num=0;
  for E2nd in E2nd_anal:
    eta_eff           = F2/AD**(1/n)/E2nd**((n-1)/n)*eterm;     # effective viscosity

    T2nd_anal[num]    = 2*eta_eff*E2nd;
    Eta_eff_anal[num] = eta_eff;
    num               = num+1;

  # Compute analytical solution @ same points as LaMEM solution
  T2nd_LaMEM_anal  = np.zeros(len(data.E2nd_bg))
  RMS_err_Numerics = np.zeros(len(data.E2nd_bg))
  num=0;
  for E2nd in data.E2nd_bg:
    eta_eff               = F2/AD**(1/n)/E2nd**((n-1)/n)*eterm;    # effectic viscosity

    T2nd_LaMEM_anal[num]  = 2*eta_eff*E2nd;    # in MPa
    RMS_err_Numerics[num] = (T2nd_LaMEM_anal[num]/1e6-data.T2nd[num]) #error in MPa 
    #print('T2nd: [LaMEM;Anal;Ratio] = [',data.T2nd[num],';',T2nd_LaMEM_anal[num]/1e6,';',data.T2nd[num]/(T2nd_LaMEM_anal[num]/1e6),']')
    num                   = num+1;

  Error_L2 = np.linalg.norm(RMS_err_Numerics);
  print('Test: L2 error-norm [MPa] = ',Error_L2)

  # store data
  data.T2nd_LaMEM_anal    =   T2nd_LaMEM_anal/1e6;
  data.T2nd_anal          =   T2nd_anal/1e6;
  data.E2nd_anal          =   E2nd_anal;
  data.Eta_eff_anal       =   Eta_eff_anal;

  return(data);
#----------------------------------------------------

#----------------------------------------------------
# Plot time-dependent data of LaMEM
def PlotTimeDependentData(data,plotName):
  # load required modules for reading data
  try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import datetime 
  except:
    print('matplotlib toolbox not installed; cannot plot data')

  # Time and date stamp on title
  title_str = datetime.datetime.now().strftime("Test performed %H:%M on %B %d, %Y");

  fig1 = plt.figure()
  ax1 = fig1.add_subplot(111)
  ax1.plot(data.Time_anal,  data.T2nd_anal,'b-',linewidth=1);   # Stress analytical
  ax1.plot(data.Time,       data.T2nd,'ro',markersize=1);        # LaMEM
  
  ax1.legend(['analytical','numerical LaMEM'],fontsize='x-small')  
  ax1.set_ylabel('T2nd [MPa]');
  ax1.set_xlabel('Time [Myrs]');
  #ax1.axis([0, 0.0005, 0, 21])

  if hasattr(data, 'T2nd_anal_noLocal'):
    # in case with local iterations, also show the solution w/out local iterations
    ax1.plot(data.Time_anal,  data.T2nd_anal_noLocal,'k--',linewidth=0.5);   # Stress analytical w/out local iterations
    ax1.legend(['analytical (local iterations for e_vis)','numerical LaMEM','analytical (no local iterations, using e_total for eta_eff)'],fontsize='x-small')  


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
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import datetime 
  except:
    print('matplotlib toolbox not installed; cannot plot data')

  # Time and date stamp on title
  title_str = datetime.datetime.now().strftime("Test performed %H:%M on %B %d, %Y");

  fig1 = plt.figure()
  ax1 = fig1.add_subplot(211)
  ax1.loglog(data.E2nd_anal,  data.T2nd_anal,'b-');   # Stress analytical
  ax1.loglog(data.E2nd_bg,    data.T2nd,'r+');        # LaMEM
  
  ax1.legend(['analytical','numerical'],fontsize='x-small')  
  ax1.set_ylabel('T2nd [MPa]');
  #ax1.set_xlabel('E2nd [1/s]');

  title_str = 'LaMEM vs. analytics \n ' + title_str;
  plt.title(title_str,fontsize='xx-small')

  ax2 = fig1.add_subplot(212)
  ax2.loglog(data.E2nd_anal,  data.Eta_eff_anal,'b-');   # Stress analytical
  ax2.loglog(data.E2nd_bg,    data.Eta_eff,'r+');        # LaMEM
  ax2.legend(['analytical','numerical'],fontsize='x-small')  
  ax2.set_xlabel('E2nd [1/s]');
  ax2.set_ylabel('Eta_eff [Pas]');

  fig1.savefig(plotName,dpi=360)    # save figure to disk



  fig2 = plt.figure()
  ax1 = fig2.add_subplot(111)
  ax1.semilogx(data.E2nd,  data.T2nd_LaMEM_anal/data.T2nd,'b-');   # Stress analytical
  
  ax1.set_ylabel('T2nd_LaMEM/T2nd_analytical [ ]');
  ax1.set_xlabel('E2nd [1/s]');
  
  plotName_stress = plotName[:-4]+'_stress.png'
  fig2.savefig(plotName_stress,dpi=360)    # save figure to disk



  return();

#----------------------------------------------------




#data = LoadData();    # Load the data using the VTK toolbox
#data = AnalyticalSolution_VE(data);                        # Compute analytical solution
#PlotData(data,'t13_ViscoElastic_output.png');             # Create Plot
