# This will compute the systematics of growthrate vs eta for a single wavelength & gen
#  
import os
import numpy as np
from vtk import vtkXMLStructuredGridReader
from vtk.util import numpy_support as VN
import matplotlib.pyplot as plt  



# Basic setup command:
BasicCommand = " ../../bin/opt/LaMEM -ParamFile t15_RTI.dat ";

# Specific options for the current setup
SetupOptions = " -noslip 0,0,0,0,1,0 -open_top_bound 1 -FreeSurf_AmplCos 1e-3 -eta[1] 1 -eta[0] 1e-3  -surf_level 0.5 ";

# Output directory name
OutputDir = "Timestep_00000001_1.10000000e+00/"

# Output file 
OutFile   = "RTI_test_surf_p00000000.vts"


lam_anal     = np.arange(0.05, 5, 5/50)

# Create arrays with material constants
q_num  = np.zeros(len(lam_anal))
wav    = np.zeros(len(lam_anal))

num = 0;
for lam in lam_anal:


  LamCommand = "-FreeSurf_Wavelength "+str(lam)+" -coord_x "+str(-lam/2)+","+str(lam/2)

  LaunchCommand = BasicCommand+SetupOptions+LamCommand;
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
  wav[num]     =   bounds[1]-bounds[0];
  q_num[num]   =   max(Vel_v[:,2])/max(Ampl_v);
  num          =   num+1;



# output
print("Wavelength=",wav," q=",q_num)



# Create plot -----------------------------------------------------------------------------
fig1 = plt.figure()
#  ax1 = fig1.add_subplot(111)
plt.plot(wav,      q_num,'r+-');      # speed numerical

plt.ylabel('Growthrate q [ ]');
plt.xlabel('$\lambda/H$ [ ]');

title_str = 'LaMEM growthrate \n ';
plt.title(title_str,fontsize='xx-small')

  
fig1.savefig("Growthrate.png",dpi=360)    # save figure to disk
# -----------------------------------------------------------------------------------------




# save output to file
np.savetxt('growthrate.out', (wav,q_num))   