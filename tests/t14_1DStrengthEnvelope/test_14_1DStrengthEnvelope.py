## test 14 1D Strength Envelope

import os
import pyTestHarness.unittest as pth
import pyTestHarness.launcher as launch
import re
import pickle

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
  # sort Tiemsteps
  dirlist.sort()
  numSteps = len(dirlist)

  # preallocate
  nn    = 65
  Time  = np.zeros((1,numSteps))
  Z     = np.zeros((nn,numSteps))
  Temp  = np.zeros((nn,numSteps))
  Phase = np.zeros((nn,numSteps))
  Tau   = np.zeros((nn,numSteps))
  TauY  = np.zeros((nn,numSteps))
  Pres  = np.zeros((nn,numSteps))

  # loop over timesteps
  for dir,iStep in zip(dirlist,range(0,numSteps)):
    file_in  =  'output.pvtr'
    filesep  = '/'
    filename = dir+filesep+file_in
    reader   = vtkXMLPRectilinearGridReader()
    reader.SetFileName(filename)
    reader.Update()
    data     = reader.GetOutput()
  
    # extract coordinates
    x              = VN.vtk_to_numpy(data.GetXCoordinates())
    y              = VN.vtk_to_numpy(data.GetYCoordinates())
    z              = VN.vtk_to_numpy(data.GetZCoordinates())
    Z[:,iStep]     = z
  
    # dimensions
    nx             = np.size(x)
    ny             = np.size(y)
    nz             = np.size(z)
  
    # Read Time
    # split directory name to extract time
    stepName       = dir.split('/')
    stepName       = stepName[3]
    time           = stepName.split('_')
    t              = float(time[2])
    Time[0,iStep]  = t
  
    # load data from VTK file & extratct central column
    T              = VN.vtk_to_numpy(data.GetPointData().GetArray('temperature [C]'))
    T              = np.reshape(T    ,(nz,ny,nx))
    Temp[:,iStep]  = T[:,1,1]
    Ph             = VN.vtk_to_numpy(data.GetPointData().GetArray('phase [ ]'))
    Ph             = np.reshape(Ph   ,(nz,ny,nx))
    Phase[:,iStep] = Ph[:,1,1]
    T2nd           = VN.vtk_to_numpy(data.GetPointData().GetArray('j2_dev_stress [MPa]'))
    T2nd           = np.reshape(T2nd ,(nz,ny,nx))
    Tau[:,iStep]   = T2nd[:,1,1]
    Yi             = VN.vtk_to_numpy(data.GetPointData().GetArray('yield [MPa]'))
    Yi             = np.reshape(Yi   ,(nz,ny,nx)) 
    TauY[:,iStep]  = Yi[:,1,1]
    P              = VN.vtk_to_numpy(data.GetPointData().GetArray('pressure [MPa]'))
    P              = np.reshape(P    ,(nz,ny,nx))
    Pres[:,iStep]  = P[:,1,1]
    
  # combine into one structure
  output = {
    "Time"  : Time,
    "Z"     : Z,
    "Temp"  : Temp,
    "Phase" : Phase,
    "Tau"   : Tau,
    "TauY"  : TauY,
    "Pres"  : Pres
  }

  # remove output
  for dir in dirlist:
    shutil.rmtree(dir)
    
  return(output)
#----------------------------------------------------

#----------------------------------------------------
# compute analytical solution
def Analytical(data):

  import numpy as np

  # list of parameters for analytical solution
  # Phases: 0: Air, 1: WetQuartzite, 2: Granite, 3: DryOlivine
  G         = [2e10, 3e10, 5e10, 7e10]
  Bn        = [5e-19, 1.55371e-17, 1.67675e-25, 1.48058e-16]
  En        = [0, 154e3, 186.5e3, 532e3]
  n         = [1, 2.3, 3.3, 3.5]
  Vn        = [0,0,0,1.7e-5]
  maxPhase  = len(G) - 1
  
  # gas constant
  R         = 8.314463
  # strainrate
  eps_tot   = 1e-15

  # get number of steps
  numSteps  = np.size(data["Time"])

  # copy most things from input
  output = {
    "Time"  : data["Time"],
    "Z"     : data["Z"],
    "Temp"  : data["Temp"],
    "Phase" : data["Phase"],
    "TauY"  : data["TauY"],
    "Pres"  : data["Pres"]
  }

  # preallocate
  nn        = 65
  Tau       = np.zeros((nn,numSteps))

  # compute
  for iStep in range(0,numSteps):
    for iEl in range(0,nn):
      T     = np.squeeze(output["Temp"][iEl,iStep]) + 273.15
      Yi    = np.squeeze(output["TauY"][iEl,iStep])
      P     = np.squeeze(output["Pres"][iEl,iStep]) *1e6
      ph    = np.squeeze(output["Phase"][iEl,iStep])
      ph1   = int(np.floor(ph))
      if ph1 < maxPhase:
        ph2 = ph1 + 1
      else:
        ph2 = ph1
      phfrac= ph - ph1
      Tau1  = 1e-6 * Bn[ph1]**(-1/n[ph1]) * eps_tot**(1/n[ph1]) * np.exp((En[ph1]+Vn[ph1]*P)/(n[ph1]*R*T))
      Tau2  = 1e-6 * Bn[ph2]**(-1/n[ph2]) * eps_tot**(1/n[ph2]) * np.exp((En[ph2]+Vn[ph2]*P)/(n[ph2]*R*T))
      TauAna= Tau1 * (1-phfrac) + Tau2 * phfrac
      if TauAna > Yi:
        TauAna = Yi
      Tau[iEl,iStep] = TauAna
    
  # add Tau to output
  output["Tau"] = Tau

  return(output)
#----------------------------------------------------  

#----------------------------------------------------
# interpolate different runs on the same timesteps (the ones from the first run)
def interpRuns(*args):

  import numpy as np

  numRuns = len(args)

  # set timesteps
  Time     = args[0]["Time"]
  Time     = Time[Time < 0.45]
  numSteps = len(Time)
  
  # preallocate
  nn       = 65
  Z        = np.zeros((nn,numRuns,numSteps))
  Tau      = np.zeros((nn,numRuns,numSteps))
  Temp     = np.zeros((nn,numRuns,numSteps))
  TauF     = np.zeros(( 1,numRuns,numSteps))

  #ipdb.set_trace()

  # interpolate
  for Run,iRun in zip(args,range(0,numRuns)):
    for iEl in range(0,nn):
      Z[iEl,iRun,:]    = np.interp(Time,np.squeeze(Run["Time"]),np.squeeze(Run["Z"][iEl,:]))
      Tau[iEl,iRun,:]  = np.interp(Time,np.squeeze(Run["Time"]),np.squeeze(Run["Tau"][iEl,:]))
      Temp[iEl,iRun,:] = np.interp(Time,np.squeeze(Run["Time"]),np.squeeze(Run["Temp"][iEl,:]))
    TauF[0,iRun,:]     = np.sum(np.squeeze(Tau[:,iRun,:]),0)

  # combine in dict
  Runs = {
    "Time": Time,
    "Z": Z,
    "Tau": Tau,
    "TauF": TauF,
    "Temp": Temp
    }

  return(Runs)
#----------------------------------------------------

#----------------------------------------------------
# make plots
def makePlot(Runs):
  try: 
    import matplotlib
    matplotlib.use('Agg')
 
    import numpy as np
    import matplotlib.pyplot as plt
    import datetime
  except:
    print('matplotlib toolbox not installed; cannot plot data')


  Time = Runs["Time"]
  Z    = Runs["Z"]
  Tau  = Runs["Tau"]
  TauF = Runs["TauF"]
  Temp = Runs["Temp"]

  numSteps = len(Time)

  # figure 1
  fig1, axs = plt.subplots(1,2, gridspec_kw={'width_ratios': [1,1]})
  
  p0, = axs[0].plot(np.squeeze(Tau[:,0,numSteps-1]),np.squeeze(Z[:,0,numSteps-1]),'m',label='VP')
  axs[1].plot(np.squeeze(Temp[:,0,numSteps-1]),np.squeeze(Z[:,0,numSteps-1]),'m')
  
  p4, = axs[0].plot(np.squeeze(Tau[:,4,numSteps-1]),np.squeeze(Z[:,4,numSteps-1]),'k--',label='Analytical Sol.')
  axs[1].plot(np.squeeze(Temp[:,4,numSteps-1]),np.squeeze(Z[:,4,numSteps-1]),'k--')
  
  p1, = axs[0].plot(np.squeeze(Tau[:,1,numSteps-1]),np.squeeze(Z[:,1,numSteps-1]),'gs',label='VEP_5ka')
  axs[1].plot(np.squeeze(Temp[:,1,numSteps-1]),np.squeeze(Z[:,1,numSteps-1]),'gs')
  
  p2, = axs[0].plot(np.squeeze(Tau[:,2,numSteps-1]),np.squeeze(Z[:,2,numSteps-1]),'r*',label='VEP_10ka')
  axs[1].plot(np.squeeze(Temp[:,2,numSteps-1]),np.squeeze(Z[:,2,numSteps-1]),'r*')
  
  p3, = axs[0].plot(np.squeeze(Tau[:,3,numSteps-1]),np.squeeze(Z[:,3,numSteps-1]),'b+',label='VEP_50ka')
  axs[1].plot(np.squeeze(Temp[:,3,numSteps-1]),np.squeeze(Z[:,3,numSteps-1]),'b+')
  
  # add legend
  axs[1].legend(handles=[p0,p1,p2,p3,p4])
  
  # axis labels
  axs[0].set_xlabel(r'$\tau_{II}$ [MPa]', fontsize=14)
  axs[0].set_ylabel('Depth [km]', fontsize=14)
  axs[0].tick_params(axis='both', which='major', labelsize=12)
  axs[1].set_xlabel(r'Temp $[^{\circ}C]$', fontsize=14)
  axs[1].set_yticks([],[])
  
  # time and date stamp on title
  title_str = datetime.datetime.now().strftime("Test performed %H:%M on %B %d, %Y")
  title_str = '1D Strength profile: \n ' + title_str
  axs[0].set_title(title_str,fontsize=10)

  # save figure
  fig1 = plt.gcf()
  fig1.savefig('./t14_1DStrengthEnvelope/1D_StrengthEnvelope.png',dpi=360)
  print('Created figure: ./t14_1DStrengthEnvelope/1D_StrengthEnvelope.png \n')

  
  # figure 2
  fig2 = plt.figure(figsize=(5,3))
  axs2 = fig2.add_subplot(111)
  p0,  = axs2.plot(Time,np.squeeze(TauF[0,0,:]),'m',label='VP')
  p4,  = axs2.plot(Time,np.squeeze(TauF[0,4,:]),'k--',label='Analytical Sol.')
  p1,  = axs2.plot(Time,np.squeeze(TauF[0,1,:]),'gs',label='VEP_5ka')
  p2,  = axs2.plot(Time,np.squeeze(TauF[0,2,:]),'r*',label='VEP_10ka')
  p3,  = axs2.plot(Time,np.squeeze(TauF[0,3,:]),'b+',label='VEP_50ka')

  # add legend
  axs2.legend(handles=[p0,p1,p2,p3,p4])

  # axis labels
  axs2.set_xlabel('time [Ma]', fontsize=14)
  axs2.set_ylabel(r'$\sum \tau_{II}$ [MPa]', fontsize=14)
  axs2.tick_params(axis='both', which='major', labelsize=12)

  # save figure
  fig2 = plt.gcf()
  fig2.tight_layout(pad=1.0)
  fig2.savefig('./t14_1DStrengthEnvelope/1D_StressBuild.png',dpi=360)
  print('Created figure: ./t14_1DStrengthEnvelope/1D_StressBuild.png \n')

#----------------------------------------------------

# first test runs visco-plastic setup with dt = 10 ka
def test_a():
  ranks         = 1
  launch        = ['rm -r Timestep*',
                   'rm -rf ./t14_1DStrengthEnvelope/OutVP; mkdir ./t14_1DStrengthEnvelope/OutVP',
                   '../bin/opt/LaMEM -ParamFile ./t14_1DStrengthEnvelope/1D_VP.dat',
                   'mv Timestep* ./t14_1DStrengthEnvelope/OutVP/']
  expected_file = 't14_1DStrengthEnvelope/t14_1D_VP_Direct_opt-p1.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)
 
  # Create unit test object
  ex1 = pth.pthUnitTest('t14_1D_VP_Direct_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')
 
  return(ex1)

# 2nd test runs visco-elasto-plastic setup with dt = 5 ka
def test_b():
  ranks         = 1
  launch        = ['rm -rf ./t14_1DStrengthEnvelope/OutVEP5; mkdir ./t14_1DStrengthEnvelope/OutVEP5',
                   '../bin/opt/LaMEM -ParamFile ./t14_1DStrengthEnvelope/1D_VEP5.dat',
                   'mv Timestep* t14_1DStrengthEnvelope/OutVEP5/']
  expected_file = 't14_1DStrengthEnvelope/t14_1D_VEP5_Direct_opt-p2.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)
 
  # Create unit test object
  ex1 = pth.pthUnitTest('t14_1D_VEP5_Direct_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')
 
  return(ex1)

# 3rd test runs visco-elasto-plastic setup with dt = 10 ka
def test_c():
  ranks         = 1
  launch        = ['rm -rf ./t14_1DStrengthEnvelope/OutVEP10; mkdir ./t14_1DStrengthEnvelope/OutVEP10',
                   '../bin/opt/LaMEM -ParamFile ./t14_1DStrengthEnvelope/1D_VEP10.dat',
                   'mv Timestep* t14_1DStrengthEnvelope/OutVEP10/']
  expected_file = 't14_1DStrengthEnvelope/t14_1D_VEP10_Direct_opt-p3.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)
 
  # Create unit test object
  ex1 = pth.pthUnitTest('t14_1D_VEP10_Direct_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')
 
  return(ex1)

# 4th test runs visco-plastic setup with dt = 50 ka
def test_d():
  ranks         = 2
  launch        = ['rm -rf ./t14_1DStrengthEnvelope/OutVEP50; mkdir ./t14_1DStrengthEnvelope/OutVEP50',
                   '../bin/opt/LaMEM -ParamFile ./t14_1DStrengthEnvelope/1D_VEP50.dat',
                   'mv Timestep* t14_1DStrengthEnvelope/OutVEP50/ 2>/dev/null']
  expected_file = 't14_1DStrengthEnvelope/t14_1D_VEP50_Direct_opt-p4.expected'

  def comparefunc(unittest):

    key = re.escape("|Div|_inf")
    unittest.compareFloatingPoint(key,1e-7)

    key = re.escape("|Div|_2")
    unittest.compareFloatingPoint(key,1e-5)

    key = re.escape("|mRes|_2")
    unittest.compareFloatingPoint(key,1e-4)

    # read timesteps
    try:
      VP    = LoadData('./t14_1DStrengthEnvelope/OutVP')
      VEP5  = LoadData('./t14_1DStrengthEnvelope/OutVEP5')
      VEP10 = LoadData('./t14_1DStrengthEnvelope/OutVEP10')
      VEP50 = LoadData('./t14_1DStrengthEnvelope/OutVEP50')
      
      # compute analytical solution from VP run
      AnSo  = Analytical(VP)
    
      # interpolate all runs onto the saem timesteps
      Runs  = interpRuns(VP,VEP5,VEP10,VEP50,AnSo)
    
      # plot
      makePlot(Runs)
 
    except:
      print('VTK/MatPlotLib/NumPy toolboxes are not installed; will not create plots')

    
   
  # Create unit test object
  ex1 = pth.pthUnitTest('t14_1D_VEP50_Direct_opt',ranks,launch,expected_file)
  ex1.setVerifyMethod(comparefunc)
  ex1.appendKeywords('@')
 
  return(ex1)