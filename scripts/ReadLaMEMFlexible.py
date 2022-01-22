#! /usr/bin/env python
#----------------------------------------------------
# flexible script to read LaMEM output
# arguments:
# 1. -n or --name=        : out_file_name of the run
# 2. -s or --steps=       : which timesteps [all or last]
# 3. -f or --format=      : which format to save them in [.mat or .csv]
# 4. -o or --outputs=     : list all LaMEM variables that should be read, can also use 'all'
# 5. -d or --destination= : destination of output files
#----------------------------------------------------
def main():
  # load required modules for reading data
  try: 
    from vtk import vtkXMLGenericDataObjectReader
    #from vtk import vtkXMLPRectilinearGridReader
    from vtk.util import numpy_support as VN
    import sys
    import getopt
    import numpy as np
    import os
    import csv
    import glob
    import scipy.io
  except:
    print 'Import failed. Most likely vtk. \n'

  # setup getopt
  shortOpt = "n:s:f:o:d:h"
  longOpt  = ["name=", "steps=", "format=", "outputs=", "destination=", "help"]

  # get all arguments
  argList = sys.argv[1:]

  # check options given
  try:
    ArgVal, empty = getopt.getopt(argList, shortOpt, longOpt)
  except getopt.error as err:
    # Output error, and return with an error code
    print (str(err))
    sys.exit(2)
  
  # read options
  for arg, val in ArgVal:
    if arg in ("-n", "--name"):
      name       = val
    elif arg in ("-s", "--steps"):
      steps      = val
    elif arg in ("-f", "--format"):
      fileFormat = val
    elif arg in ("-o", "--outputs"):
      outputs    = val.split(' ')
    elif arg in ("-d", "--destination"):
      dest       = val
    elif arg in ("-h", "--help"):
      displayHelp()
      exit()

  # check arguments
  if not('name' in locals()):
    sys.exit('Please specify the name of the run (-n, --name=)')
  if (not('steps' in locals()) or not(steps == 'all' or steps == 'last')):
    print 'Option steps not set or understood. Using: all'
    steps = 'all'
  if (not('fileFormat' in locals()) or not(fileFormat == 'csv' or fileFormat == 'mat')):
    print 'Option format not set or understood. Using: mat'
    fileFormat = 'mat'
  if not('outputs' in locals()):
    sys.exit('No outputs specified')
  else:
    if outputs[0] == 'all':
      allOutFlag = True
      outputs    = []
    else:
      allOutFlag = False
  if ('dest' in locals()):
    if not os.path.exists(dest):
      os.makedirs(dest)
      print 'Created output directory: '+dest
  else:
    print 'Option destination not set or understood. Using Timestep directories'
  
  # check for surface and outputs
  surfaceFlag = False
  for out in outputs:
    if out.startswith('surf_'):
      surfaceFlag = True
      break

  # check for passive tracers
  tracerFlag = False
  for out in outputs:
    if out.startswith('tracer_'):
      tracerFlag = True
      break
  
  # set up names
  namePVD       = name+'.pvtr'
  nameSurfPVD   = name+'_surf.pvts'
  nameTracerPVD = name+'_passive_tracers.pvtu'

  # get list of directories
  dirList = glob.glob('./Timestep_*')
  dirList.sort(reverse=False)
  if steps == 'last':
    dirList = findSteps(dirList,name,'last',surfaceFlag)
  else:
    dirList = findSteps(dirList,name,'all',surfaceFlag)

  # loop through list
  for dir in dirList:
    # initialize reader for .pvtr
    filename = os.path.join(dir,namePVD)
    reader   = vtkXMLGenericDataObjectReader()
    #reader   = vtkXMLPRectilinearGridReader()
    reader.SetFileName(filename)
    reader.Update()
    data     = reader.GetOutput()

    # extract coordinates
    x = VN.vtk_to_numpy(data.GetXCoordinates())
    y = VN.vtk_to_numpy(data.GetYCoordinates())
    z = VN.vtk_to_numpy(data.GetZCoordinates())

    # dimensions
    nx = np.size(x)
    ny = np.size(y)
    nz = np.size(z)

    # get info from filename
    time    = float(dir.split('_')[2])
    stepNum = 'Timestep_'+dir.split('_')[1]

    # set output path
    if ('dest' in locals()):
      outPath = os.path.join(dest,stepNum)
      if not os.path.exists(outPath):
        os.mkdir(outPath)
    else:
      outPath = dir

    # write coordinates
    if fileFormat == 'mat':
      scipy.io.savemat(os.path.join(outPath,'Coords.mat'), {'x_vec':x,'y_vec':y,'z_vec':z})
    else:
      with open(os.path.join(outPath,'Coords.csv'),'wb') as file:
        writer = csv.writer(file)
        writer.writerow(x)
        writer.writerow(y)
        writer.writerow(z)
    
    # write time
    if fileFormat == 'mat':
      scipy.io.savemat(os.path.join(outPath,'Time.mat'), {'time':time})
    else:
      with open(os.path.join(outPath,'Time.csv'),'wb') as file:
        writer = csv.writer(file)
        writer.writerow([time])

    # in case of all outputs
    if allOutFlag:
      numFields = data.GetPointData().GetNumberOfArrays()
      for iField in range(numFields):
        outputs.append(data.GetPointData().GetArrayName(iField).split()[0])


    # loop through outputs
    for out in outputs:
      accessSaveOutput(outPath,data,fileFormat,out,nx,ny,nz)

    if (surfaceFlag or allOutFlag):  
      # initialize reader for .pvts
      filename = os.path.join(dir,nameSurfPVD)
      reader   = vtkXMLGenericDataObjectReader()
      #reader   = vtkXMLPRectilinearGridReader()
      reader.SetFileName(filename)
      reader.Update()
      data     = reader.GetOutput()

      # in case of all outputs
      if allOutFlag:
        numFields = data.GetPointData().GetNumberOfArrays()
        for iField in range(numFields):
          outputs.append('surf_'+data.GetPointData().GetArrayName(iField).split()[0])
      
      # loop through outputs
      for out in outputs:
        accessSurfSaveOutput(outPath,data,fileFormat,out,nx,ny)

    if (tracerFlag or allOutFlag):
      # initialize reader
      filename = os.path.join(dir,nameTracerPVD)
      reader   = vtkXMLGenericDataObjectReader()
      #reader   = vtkXMLPRectilinearGridReader()
      reader.SetFileName(filename)
      reader.Update()
      data     = reader.GetOutput()

      # in case of all outputs
      if allOutFlag:
        numFields = data.GetPointData().GetNumberOfArrays()
        for iField in range(numFields):
          outputs.append('tracer_'+data.GetPointData().GetArrayName(iField).split()[0])

      # loop through outputs
      for out in outputs:
        accessTracerSaveOutput(outPath,data,fileFormat,out,nx,ny)
      
    # print progress
    print 'Finished '+dir



def findSteps(dirList,name,steps,surfaceFlag):
  import os
  # sort reverse so we can start searching from the end
  dirList.sort(reverse=True)
  # make outList
  outList = []
  iter = 0
  print "Searching for "+name
  # loop through all timestep dirs
  for dir in dirList:
    # get list of files
    list  = os.listdir(dir)
    found = lookForRun(list,name,surfaceFlag)
    if (found and steps == 'last'):
      print "Found in directory: "+dir
      outList.append(dir)
      return outList
    elif found:
      outList.append(dir)
      iter += 1
  print "Found in "+str(iter)+" directories"
  return outList

def lookForRun(list,name,surfaceFlag):
  if surfaceFlag:
    if (name+'.pvtr' in list and name+'_surf.pvts' in list): 
      return True
    else:
      return False
  else:
    if name+'.pvtr' in list:
      return True
    else:
      return False  

def accessSaveOutput(outPath, data, fileFormat, out, nx, ny, nz):
  from vtk.util import numpy_support as VN
  import numpy as np
  import os
  import csv
  import scipy.io

  isVec = False
  isTen = False
  # scalars
  if out == 'cont_res':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('cont_res [1/s]'))
    outName  = 'ContRes'
  elif out == 'density':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('density [kg/m^3]'))
    outName  = 'Density'
  elif out == 'j2_dev_stress':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('j2_dev_stress [MPa]'))
    outName  = 'J2_DevStress'
  elif out == 'j2_strain_rate':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('j2_strain_rate [1/s]'))
    outName  = 'J2_StrainRate'
  elif out == 'litho_press':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('litho_press [MPa]'))
    outName  = 'Lith_Pressure'
  elif out == 'over_press':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('over_press [MPa]'))
    outName  = 'OverPressure'
  elif out == 'phase':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('phase [ ]'))
    outName  = 'Phase'
  elif out == 'plast_dissip':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('plast_dissip [W/m^3]'))
    outName  = 'PlastDissip'
  elif out == 'plast_strain':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('plast_strain [ ]'))
    outName  = 'PlastStrain'
  elif out == 'tot_strain':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('tot_strain [ ]'))
    outName  = 'TotStrain'
  elif out == 'pore_press':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('pore_press [MPa]'))
    outName  = 'PorePressure'
  elif out == 'pressure':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('pressure [MPa]'))
    outName  = 'Pressure'
  elif out == 'eff_press':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('eff_press [MPa]'))
    outName  = 'EffPressure'
  elif out == 'rel_dif_rate':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('rel_dif_rate [ ]'))
    outName  = 'RelDifRate'
  elif out == 'rel_dis_rate':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('rel_dis_rate [ ]'))
    outName  = 'RelDisRate'
  elif out == 'rel_prl_rate':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('rel_prl_rate [ ]'))
    outName  = 'RelPrlRate'
  elif out == 'temperature':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('temperature [C]'))
    outName  = 'Temp'
  elif out == 'total_pressure':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('total_pressure [MPa]'))
    outName  = 'TotPressure'
  elif out == 'visc_creep':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('visc_creep [Pa*s]'))
    outName  = 'LogViscCreep'
  elif out == 'visc_total':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('visc_total [Pa*s]'))
    outName  = 'LogViscTotal'
  elif out == 'yield':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('yield [MPa]'))
    outName  = 'Yield'
  # vectors
  elif out == 'EHmax':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('EHmax [ ]'))
    outName  = 'EHmax'
    isVec    = True
  elif out == 'SHmax':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('SHmax [ ]'))
    outName  = 'SHmax'
    isVec    = True
  elif out == 'moment_res':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('moment_res [N/m^3]'))
    outName  = 'MomentRes'
    isVec    = True
  elif out == 'tot_displ':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('tot_displ [km]'))
    outName  = 'TotDisplacement'
    isVec    = True
  elif out == 'velocity':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('velocity [cm/yr]'))
    outName  = 'Velocity'
    isVec    = True
  # tensors
  elif out == 'dev_stress':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('dev_stress [MPa]'))
    outName  = 'DevStress'
    isTen    = True
  elif out == 'strain_rate':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('strain_rate [1/s]'))
    outName  = 'StrainRate'
    isTen    = True
  # surface
  elif (out == 'surf_amplitude' or out == 'surf_topography' or out == 'surf_velocity'):
    # ignore
    return
  # passive tracers
  elif (out == 'tracer_ID' or out == 'tracer_Active' or out == 'tracer_Mf' or out == 'tracer_Mf_Grid'
        or out == 'tracer_Phase' or out == 'tracer_Pressure' or out == 'tracer_Temperature'):
    # ignore
    return
  else:
    print "Unknown output: "+out
    return

  # reshape and swap axes
  if (not isVec and not isTen): 
    Out = np.reshape(Out,(nz,ny,nx))
    Out =	np.swapaxes(np.swapaxes(Out,0,2),1,0)
  elif isVec:
    Out = np.reshape(Out,(nz,ny,nx,3))
    Out = np.swapaxes(np.swapaxes(Out,0,2),1,0)
  elif isTen:
    Out = np.reshape(Out,(nz,ny,nx,9))
    Out = np.swapaxes(np.swapaxes(Out,0,2),1,0)
  
  # save output
  if isTen:
    if fileFormat == 'csv':
      print ".csv files don't work for Tensors, switching to "+outName+".mat"
    scipy.io.savemat(os.path.join(outPath,outName+'.mat'), {outName:Out})
  elif isVec:
    if fileFormat == 'csv':
      print ".csv files don't work for Vectors, switching to "+outName+".mat"
    scipy.io.savemat(os.path.join(outPath,outName+'.mat'), {outName:Out})
  else:
    if (fileFormat == 'csv' and ny == 2):
      Out = np.squeeze(Out[0,:,:])
      with open(os.path.join(outPath,outName+'.csv'),'wb') as file:
        writer = csv.writer(file)
        writer.writerows(Out)
    else:
      if fileFormat == 'csv':
        print ".csv files don't work for Vectors, switching to "+outName+".mat"
      scipy.io.savemat(os.path.join(outPath,outName+'.mat'), {outName:Out})
  
  return

def accessSurfSaveOutput(outPath, data, fileFormat, out, nx, ny):
  from vtk.util import numpy_support as VN
  import numpy as np
  import os
  import csv
  import scipy.io

  isVec = False
  # scalars
  if out == 'surf_amplitude':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('amplitude [km]'))
    Out = np.reshape(Out,(ny,nx))
    outName  = 'SurfAmplitude'
  elif out == 'surf_topography':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('topography [km]'))
    Out = np.reshape(Out,(ny,nx))
    outName  = 'SurfTopography'
  elif out == 'surf_velocity':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('velocity [cm/yr]'))
    Out = np.reshape(Out,(ny,nx,3))
    outName  = 'surfVelocity'
    isVec    = True
  else:
    # ignore
    return
  
  # save output
  if isVec:
    if fileFormat == 'csv':
      OutX = np.squeeze(Out[:,:,0])
      OutY = np.squeeze(Out[:,:,1])
      OutZ = np.squeeze(Out[:,:,2])
      with open(os.path.join(outPath,outName+'_X.csv'),'wb') as file:
        writer = csv.writer(file)
        writer.writerows(OutX)
      with open(os.path.join(outPath,outName+'_Y.csv'),'wb') as file:
        writer = csv.writer(file)
        writer.writerows(OutY)
      with open(os.path.join(outPath,outName+'_Z.csv'),'wb') as file:
        writer = csv.writer(file)
        writer.writerows(OutZ)
    else:
      scipy.io.savemat(os.path.join(outPath,outName+'.mat'), {outName:Out})
  else:
    if fileFormat == 'csv':
      with open(os.path.join(outPath,outName+'.csv'),'wb') as file:
        writer = csv.writer(file)
        writer.writerows(Out)
    else:
      scipy.io.savemat(os.path.join(outPath,outName+'.mat'), {outName:Out})
  
  return

def accessTracerSaveOutput(outPath, data, fileFormat, out, nx, ny):
  from vtk.util import numpy_support as VN
  import os
  import csv
  import scipy.io
  import numpy as np

  # save coordinates
  coords      = VN.vtk_to_numpy(data.GetPoints().GetData())
  with open(os.path.join(outPath,'TracerCoords.csv'),'wb') as file:
    writer = csv.writer(file)
    writer.writerows(coords)

  # save ID
  ID          = VN.vtk_to_numpy(data.GetPointData().GetArray('ID'))
  np.savetxt(os.path.join(outPath,'TracerID.csv'), ID, delimiter="'")

  # scalars
  if out == 'tracer_Active':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('Active'))
    outName  = 'TracerActive'
  elif out == 'tracer_Mf':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('Mf [ ]'))
    outName  = 'TracerMf'
  elif out == 'tracer_Mf_Grid':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('Mf_Grid [ ]'))
    outName  = 'TracerMf_Grid'
  elif out == 'tracer_Phase':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('Phase'))
    outName  = 'TracerPhase'
  elif out == 'tracer_Pressure':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('Pressure [MPa]'))
    outName  = 'TracerPressure'
  elif out == 'tracer_Temperature':
    Out = VN.vtk_to_numpy(data.GetPointData().GetArray('Temperature [C]'))
    outName  = 'TracerTemperature'
  else:
    # ignore
    return
  
  # save output
  if fileFormat == 'csv':
    np.savetxt(os.path.join(outPath,outName+'.csv'), Out, delimiter="'")
  else:
    scipy.io.savemat(os.path.join(outPath,outName+'.mat'), {outName:Out})
  
  return

def displayHelp():
  print "Use with argument-value-pairs"
  print "-n or --name=       : out_file_name of the run"
  print "-s or --steps=      : which timesteps to read [all or last]"
  print "-f or --format=     : which format to save them in [mat or csv]"
  print "-o or --outputs=    : list all LaMEM variables that should be read, can also use all"
  print "             example: -o'density j2_dev_stress surf_velovity'"
  print "-d or --destination=: name of output destination (default is the Timestep directories)"
  print "-h or --help        : display help"
  print "Version             : September 21, 2021"

#---------------------------------------------------
main()
print "Done!"
#---------------------------------------------------