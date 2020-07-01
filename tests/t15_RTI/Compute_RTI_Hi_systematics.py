# This will create a growthrate diagram for a RTI setup & plot it
#  
import os
import numpy as np
from vtk import vtkXMLStructuredGridReader
from vtk import vtkXMLPRectilinearGridReader
from vtk.util import numpy_support as VN
import matplotlib.pyplot as plt 
import math 

# Input parameters:
lam          = 1;   # wavelength
eta_low      = 1;
eta_up       = 1;
numPoints    = 25;

# Basic setup command:
BasicCommand = " ../../bin/opt/LaMEM -ParamFile t15_RTI.dat ";


# Specify wavelength
LamCommand   = "-FreeSurf_Wavelength "+str(lam)+" -coord_x "+str(-lam/2)+","+str(lam/2)

# Specific options for the current setup
SetupOptions = " -noslip 0,0,0,0,1,0 -open_top_bound 1 -nel_z 64 -nel_x 64 -FreeSurf_AmplCos 1e-4 -eta[1] "+str(eta_low)+" -eta[0] "+str(eta_up)+" ";

# Output directory name
OutputDir = "Timestep_00000001_1.10000000e+00/"

# Output file 
OutFile    = "RTI_test_surf_p00000000.vts"
OutFile1   = "RTI_test.pvtr"


log10Hi_vec     = np.arange(-2, -0.01, 2/numPoints)

# Create arrays with material constants
q_num    = np.zeros(len(log10Hi_vec))
Hi_vec   = np.zeros(len(log10Hi_vec))

num = 0;
for log10Hi in log10Hi_vec:

  Hi = 10**log10Hi;
 
  HiCommand    = " -surf_level " + str(Hi);    
  LaunchCommand = BasicCommand+SetupOptions+LamCommand+HiCommand;
  print(LaunchCommand)

  # run LaMEM
  os.system(LaunchCommand)

  # read output
  # surface
  filename = OutputDir+OutFile;
  print(filename)
  reader = vtkXMLStructuredGridReader()
  reader.SetFileName(filename)
  reader.Update()
  data = reader.GetOutput()

  # 3D velocity
  filename1 = OutputDir+OutFile1;
  reader1 = vtkXMLPRectilinearGridReader()
  reader1.SetFileName(filename1)
  reader1.Update()
  data1 = reader1.GetOutput()

  # load data from VTK file & reshape them to 3D arrays
  Vel_v1        =   VN.vtk_to_numpy(data.GetPointData().GetArray('velocity [ ]')); 
  max_Vz        =   max(Vel_v1[:,2]);

  # load data from VTK file & extract growthrate data
  Vel_v        =   VN.vtk_to_numpy(data.GetPointData().GetVectors('velocity [ ]')); 
  Ampl_v       =   VN.vtk_to_numpy(data.GetPointData().GetArray('amplitude [ ]')); 
  bounds       =   data.GetBounds();
  #q_num[num]   =   max(Vel_v[:,2])/max(Ampl_v);
  q_num[num]   =   max_Vz/max(Ampl_v);
  
  Hi_vec[num]  =   Hi;
  num          =   num+1;



# output
print("Hi=",Hi," q=",q_num)



# Create plot -----------------------------------------------------------------------------
fig1 = plt.figure()
plt.loglog(Hi_vec,     q_num,'r+-');      # speed numerical

plt.ylabel('q [ ]');
plt.xlabel('$Hi$ [ ]');

title_str = "LaMEM growthrate $\lambda$="+str(lam)+" $\eta_{lowerLayer}$="+str(eta_low)+" $\eta_{upperLayer}$="+str(eta_up)
plt.title(title_str,fontsize='xx-small')

  
fig1.savefig("Hi_systematics.png",dpi=360)    # save figure to disk
# -----------------------------------------------------------------------------------------



# save output to file
np.savetxt('Hi_systematics.out', (Hi_vec,q_num))   