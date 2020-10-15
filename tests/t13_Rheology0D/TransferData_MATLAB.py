from vtk import vtkXMLPRectilinearGridReader
from vtk.util import numpy_support as VN
import numpy as np
import glob
import shutil
import csv



#data        = LoadStrainrateData('./Out_DisCr/','Rheolog0D_linearViscous')

filename = 'Timestep_00000001_2.00000000e-03/Rheolog0D_CombinedCreeplaws.pvtr';
   
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
T2nd_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('j2_dev_stress [MPa]'))
E2nd_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('j2_strain_rate [1/s]'))
Pres_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('total_pressure [MPa]'))
Temp_v        =   VN.vtk_to_numpy(data.GetPointData().GetArray('temperature [C]'))
Eta_v         =   VN.vtk_to_numpy(data.GetPointData().GetArray('visc_creep [Pa*s]'))


# Compute average value of T2nd 
T2nd       = np.average(T2nd_v)
E2nd       = np.average(E2nd_v)
P          = np.average(Pres_v)
T          = np.average(Temp_v)
Eta        = np.average(Eta_v)

print('Average values: ')
print(' T2nd: ', T2nd)
print(' E2nd: ', E2nd)
print(' P   : ', P)
print(' T   : ', T)

# store data into file
data_out    = np.empty(5)
data_out[0] = T2nd;
data_out[1] = E2nd;
data_out[2] = T;
data_out[3] = P;
data_out[4] = Eta;


# save as textfile (to be read in matlab)
np.savetxt('output.txt', data_out, delimiter=" ", fmt="%e") 
