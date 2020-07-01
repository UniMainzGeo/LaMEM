# This will create a growthrate diagram for a RTI setup & plot it
#  
import os
import numpy as np
from vtk import vtkXMLStructuredGridReader
from vtk.util import numpy_support as VN
import matplotlib.pyplot as plt 
import math 

# Input parameters:
lam          = 1;   # wavelength
Hi           = 0.5; # height interface
eta_low      = 1;
numPoints    = 100;

# Basic setup command:
BasicCommand = " ../../bin/opt/LaMEM -ParamFile t15_RTI.dat ";


# Specify wavelength
LamCommand   = "-FreeSurf_Wavelength "+str(lam)+" -coord_x "+str(-lam/2)+","+str(lam/2)

# Specific options for the current setup
SetupOptions = " -noslip 0,0,0,0,1,0 -open_top_bound 1 -FreeSurf_AmplCos 1e-3 -eta[1]"+str(eta_low)+" -surf_level "+str(Hi)+" ";

# Output directory name
OutputDir = "Timestep_00000001_1.10000000e+00/"

# Output file 
OutFile   = "RTI_test_surf_p00000000.vts"


log10_eta_vec     = np.arange(-3, 3, 6/numPoints)

# Create arrays with material constants
q_num    = np.zeros(len(log10_eta_vec))
eta_vec  = np.zeros(len(log10_eta_vec))

num = 0;
for log10_eta in log10_eta_vec:

  eta = 10**log10_eta;
  EtaCommand    = " -eta[0] " + str(eta);    
  LaunchCommand = BasicCommand+SetupOptions+LamCommand+EtaCommand;
  print(LaunchCommand)

  # run LaMEM
  os.system(LaunchCommand)

  # read output
  filename = OutputDir+OutFile;
  print(filename)
  reader = vtkXMLStructuredGridReader()
  reader.SetFileName(filename)
  reader.Update()
  data = reader.GetOutput()

  # load data from VTK file & extract growthrate data
  Vel_v        =   VN.vtk_to_numpy(data.GetPointData().GetVectors('velocity [ ]')); 
  Ampl_v       =   VN.vtk_to_numpy(data.GetPointData().GetArray('amplitude [ ]')); 
  bounds       =   data.GetBounds();
  q_num[num]   =   max(Vel_v[:,2])/max(Ampl_v);
  eta_vec[num] =   eta;
  num          =   num+1;



# output
print("log10(Eta)=",log10_eta," q=",q_num)



# Create plot -----------------------------------------------------------------------------
fig1 = plt.figure()
plt.loglog(eta_vec,     q_num,'r+-');      # speed numerical

plt.ylabel('q [ ]');
plt.xlabel('$\eta_{upperLayer}$ [ ]');

title_str = "LaMEM growthrate $\lambda$="+str(lam)+" H$_i$="+str(Hi)+" $\eta_{lowerLayer}$="+str(eta_low)
plt.title(title_str,fontsize='xx-small')

  
fig1.savefig("Eta_systematics.png",dpi=360)    # save figure to disk
# -----------------------------------------------------------------------------------------



# save output to file
np.savetxt('Eta_systematics.out', (eta_vec,q_num))   