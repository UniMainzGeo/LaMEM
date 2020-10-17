import os
import re
import subprocess


#----------------------------------------------------
# Read data from disk
def LoadRTI_Data(fname):

  # load required modules for reading data
  try: 
    from vtk import vtkXMLStructuredGridReader
    from vtk.util import numpy_support as VN
    import numpy as np
    from pathlib import Path
  except:
    print('VTK toolboxes are not installed; cannot load data')


  # extract directory name and times
  E2nd_vec  = np.zeros(100)
  dir_vec   = np.empty(7,dtype='object')
  num       = 0
  fu        = [f.path for f in os.scandir('.') if f.is_dir()];

  for dirname in fu:
    if (dirname.startswith('./RTwav_') ):
      print('directory name ', dirname)
      number            = int(dirname[8:]);
      dir_vec[number]   = dirname
      num               = num+1;                  # number
     

  #Initialize values
  Vz_max    = np.empty(num);
  Ampl_max  = np.empty(num);
  lam       = np.empty(num);
  

  # Loop over all timestep directories
  for i in range(num):
    dir     =   dir_vec[i]
    file_in =   fname+'.vts'
    filesep =   '/'

    filename = dir+filesep+file_in
    print ('Loading LaMEM output from directory ', filename)

    reader = vtkXMLStructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()

    # extract coordinates
    #data.x = VN.vtk_to_numpy(data.GetXCoordinates())
    #data.y = VN.vtk_to_numpy(data.GetYCoordinates())
    #data.z = VN.vtk_to_numpy(data.GetZCoordinates())

    # dimensions
    #nx = np.size(data.x); ny = np.size(data.y); nz = np.size(data.z)

    # load data from VTK file & reshape them to 3D arrays
    Vel_v        =   VN.vtk_to_numpy(data.GetPointData().GetVectors('velocity [ ]')); 
    Ampl_v       =   VN.vtk_to_numpy(data.GetPointData().GetArray('amplitude [ ]')); 
    bounds       =   data.GetBounds();
    wavelength   =   bounds[1]-bounds[0];
    #print(Vel_v)
    
    # Compute max velocity
    Vz_max[i]     =  max(Vel_v[:,2]);
    Ampl_max[i]   =  max(Ampl_v);
    lam[i]        =  wavelength;
    print('lamda/Vz_max, Ampl, q',wavelength, Vz_max[i],Ampl_max[i], Vz_max[i]/Ampl_max[i])

  # save data for the different strainrate components
  data.Vz_max     = Vz_max;
  data.Ampl_max   = Ampl_max;
  data.lam        = lam;
  data.q          = Vz_max/Ampl_max;

  return(data);

#----------------------------------------------------

#----------------------------------------------------
# Analytical solution RT
def AnalyticalSolution_FreeSlip(data):
  import numpy as np
  import math

  lam_anal     = np.arange(1e-9, 5, 1/100)

  # Create arrays with material constants
  q_anal  = np.zeros(len(lam_anal))
  num=0;
  for lam in lam_anal:
    k            = 2*math.pi/lam;
    q_anal[num]  = ((k**2+2)*math.exp(-k) - math.exp(-2*k) - 1.0)/(4*k*(-2*k*math.exp(-k) + math.exp(-2*k) - 1));
    
    num             = num+1;

  data.lam_anal = lam_anal;
  data.q_anal   = q_anal;

  return(data);

#----------------------------------------------------
# Plot data from disk
def PlotRT_Data(data,plotName):
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
  ax1.plot(data.lam,       data.q,'ro');      # speed numerical
  ax1.plot(data.lam_anal, data.q_anal,'-');   # analytical solution

  ax1.legend(['analytical','numerical'],fontsize='x-small')  
  ax1.set_ylabel('Growthrate');
  ax1.set_xlabel('wavelength lambda/H');

  title_str = 'LaMEM vs. analytics \n ' + title_str;
  plt.title(title_str,fontsize='xx-small')

  plotName = './t15_RTI/'+plotName;
  fig1.savefig(plotName,dpi=360)    # save figure to disk

  return();

#----------------------------------------------------



#data = LoadRTI_Data('RTI_test_surf_p00000000');

#data = LoadData();    # Load the data using the VTK toolbox
#data = AnalyticalSolution_VE(data);                        # Compute analytical solution
#PlotData(data,'t13_ViscoElastic_output.png');             # Create Plot
