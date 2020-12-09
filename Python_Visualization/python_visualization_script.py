# -*- coding: utf-8 -*-
import sys,os,fnmatch
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import shutil
from time import perf_counter 
import glob
import re
from numpy.ma import masked_array
from  classes_visualization import * 
from vtk.numpy_interface import dataset_adapter as dsa
import getopt
import argparse

"""
Import useful function 
"""
from AUX_FUN import *
from AUX_FUN import _parse_line
from AUX_FUN import _file_list
from AUX_FUN import _figure_maker
from AUX_FUN import _figure_maker_phase
from AUX_FUN import _parse_grid


from initial_input import * 

"""
# setup getopt
shortOpt = "F:pts:tn"
longOpt  = ["Folder=", "ptsave=", "Test_Name="]

  # get all arguments
argList = sys.argv[1:]

  # check options given

  # read options
for arg, val in ArgVal:
    if arg in ("-F", "--Folder"):
      Folder       = val
    elif arg in ("-pts", "--ptsave"):
      ptsave      = val
    elif arg in ("-tn", "--Test_Name"):
      Test_Name = val
"""
parser = argparse.ArgumentParser()
parser.add_argument("Folder", help="Folder of the Tests",type=str)
parser.add_argument("ptsave",help="Folder where to create the output",type=str)
parser.add_argument("Test_Name",help=" Name of the Test",type=str)
args = parser.parse_args()
Folder = args.Folder
ptsave = args.ptsave
Test_Name = args.Test_Name
print("Processing the",Test_Name," in the Folder",Folder," and saving the output into",ptsave)


####################################
if not os.path.isdir(ptsave):
  os.mkdir(ptsave)
ptsave=os.path.join(ptsave,Test_Name)
if not os.path.isdir(ptsave):
  os.mkdir(ptsave)

if not os.path.isdir(ptsave):
  os.mkdir(ptsave)


reader=vtk.vtkXMLGenericDataObjectReader()
VTK_SET=reader.SetFileName
VTK_UPDATE=reader.Update
VTK_RELEASE=reader.UnRegister
VTK_POINTS=vtk.vtkPoints
VTK_OUT=reader.GetOutput
AllPoints=VTK_POINTS()
GetData=AllPoints.GetData
nodes_vtk_array= GetData()

VTK_SET_dyn=reader.SetFileName
VTK_UPDATE_dyn=reader.Update
VTK_RELEASE_dyn=reader.UnRegister
VTK_POINTS_dyn=vtk.vtkPoints
VTK_OUT_dyn=reader.GetOutput
AllPoints_dyn=VTK_POINTS()
GetData_dyn=AllPoints_dyn.GetData
nodes_vtk_array_dyn= GetData_dyn()

reader_ptt = vtk.vtkXMLGenericDataObjectReader()
VTK_SET_ptt=reader_ptt.SetFileName
VTK_UPDATE_ptt=reader_ptt.Update
VTK_RELEASE_ptt=reader_ptt.UnRegister
VTK_POINTS_ptt=vtk.vtkPoints
VTK_OUT_ptt=reader_ptt.GetOutput
AllPoints_ptt=VTK_POINTS()
GetData_ptt=AllPoints_ptt.GetData
nodes_vtk_array_ptt= GetData_ptt()    



                    
fname2 = fname+".pvd"
fname=os.path.join(Folder,Test_Name,fname2)

time, Flist, n_tim_step =_file_list(fname)    # Retrive the file list and the time step information



ipic=0 # ipic counter
tot_time_perf = 0 
for istp in Flist:
    
    t1_start = perf_counter()  
    fn=istp[1:istp.find('/')]
    t_cur=time[ipic]
    tstep=n_tim_step[ipic]
    Filename_dyn=os.path.join(Folder,Test_Name,fn,dyn)
    VTK_SET_dyn(Filename_dyn)
    VTK_UPDATE_dyn()
    VTK_OUT_dyn().GetPoints(AllPoints_dyn)        
    p=0
    nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array_dyn)
    xd,yd,zd= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1] , nodes_nummpy_array[:,2]
    nx, ny, nz=_parse_grid(Filename_dyn,p)
    xd = xd.reshape(nz,ny,nx)
    yd = yd.reshape(nz,ny,nx)
    zd = zd.reshape(nz,ny,nx)
    xd = xd[:,0,:]
    yd = yd[:,0,:]
    zd = zd[:,0,:]
    #Contour:
    
    #Quiver 
    if PTT:
        Filename_ptt=os.path.join(Folder,Test_Name,fn,passive_tracers)
        VTK_SET_ptt(Filename_ptt)
        VTK_UPDATE_ptt()
        VTK_OUT_ptt().GetPoints()
        number_marker = VTK_OUT_ptt().GetNumberOfPoints()
        PTT_Hist = PTT_max_hist(number_marker)
        for i in range(number_marker):
            
            
            PTT_Hist.x[i],PTT_Hist.y[i],PTT_Hist.z[i] = VTK_OUT_ptt().GetPoint(i)
        
        ######################################################################
        # Retrieve ptt vector
        ######################################################################
        for name in PTT_vis.field_ptt:
            buf_ptt=vtk_to_numpy(VTK_OUT_ptt().GetPointData().GetArray(name))
            if name == 'Temperature [C]':
                PTT_Hist.T=buf_ptt
            if name == 'Pressure [MPa]' :
                PTT_Hist.p = buf_ptt
            if name == 'Mf_Grid [ ]'          :
                PTT_Hist.Mf = buf_ptt
            if name == 'Active'         :
                PTT_Hist.Active = buf_ptt
            if name == 'ID'             : 
                PTT_Hist.ID =buf_ptt
            if name == 'Phase'          :
                PTT_Hist.phase = buf_ptt
                
         
    
    
    if DYN:  
        t1_start_A = perf_counter()  

        for iprop in Prop:
            Pic=eval(iprop)
            
            if(Pic.PTT_show.yes == 'yes'):
                Pic.PTT_show.PTT_structure = PTT_Hist 

                
            if Pic.contour.yes == 'yes':
                buf2=reader.GetOutput().GetPointData().GetArray(contour.field)
                buf2=vtk_to_numpy(buf2)
                buf2 = buf2.reshape(nz,ny,nx)
                Pic.contour.Scalar=buf2[:,0,:]
                
            if Pic.quiver.yes == 'yes' or Pic.stream_line.yes == 'yes':
                
                if(Pic.quiver.yes == 'yes'):
                    bufQ=reader.GetOutput().GetPointData().GetArray(quiver.field)
                else:
                    bufQ=reader.GetOutput().GetPointData().GetArray(stream_line.field)

                
                bufQ=vtk_to_numpy(bufQ)
                bufx,bufy,bufz= bufQ[:,0] , bufQ[:,1] , bufQ[:,2]
                bufx = bufx.reshape(nz,ny,nx)
                bufy = bufy.reshape(nz,ny,nx)
                bufz = bufz.reshape(nz,ny,nx)
                if Pic.quiver.yes == 'yes':
                    Pic.quiver.cx=bufx[:,0,:]
                    Pic.quiver.cy=bufy[:,0,:]
                    Pic.quiver.cz=bufz[:,0,:]
                else: 
                    Pic.stream_line.cx=bufx[:,0,:]
                    Pic.stream_line.cy=bufy[:,0,:]
                    Pic.stream_line.cz=bufz[:,0,:]

            if(Pic.n_component > 1 & Pic.n_component <=3):
                bufQ=reader.GetOutput().GetPointData().GetArray(Pic.name)
                bufQ=vtk_to_numpy(bufQ)
                bufx,bufy,bufz= bufQ[:,0] , bufQ[:,1] , bufQ[:,2]
                bufx = bufx.reshape(nz,ny,nx)
                bufy = bufy.reshape(nz,ny,nx)
                bufz = bufz.reshape(nz,ny,nx)
                n_c = 0
                for component in Pic.component_show:
                    if(component == 'mag'):
                        buf = np.sqrt(bufx**2+bufz**2)
                    if(component == 'x'):
                        buf = bufx 
                    if (component == 'z'):
                        buf = bufz
                    label = 0
                    _figure_maker(ptsave,buf,xd,zd,iprop,t_cur,tstep, Pic, label,component)
                    n_c = n_c+1

            else:
                buf=reader.GetOutput().GetPointData().GetArray(Pic.name)
                buf=vtk_to_numpy(buf)
                buf=buf.reshape((nz,ny,nx))
                label=0
                component ='none'
                if(Pic.log =='yes'):
                    buf = np.log10(buf)
                _figure_maker(ptsave,buf,xd,zd,iprop,t_cur,tstep, Pic, label,component)
            
            del Pic, buf
            
        del nodes_nummpy_array
        t2_start_A = perf_counter()  
        print("Dyn properties took ", '{:0.2f}'.format(t2_start_A-t1_start_A),"seconds")
        
        
    if SURF: 
        Filename=os.path.join(Folder,Test_Name,fn,surf)
    if PHASE: 
        t1_start_A = perf_counter()  

        Filename_ph=os.path.join(Folder,Test_Name,fn,phase)
        VTK_SET(Filename_ph)
        VTK_UPDATE()
        VTK_OUT().GetPoints(AllPoints)
        nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
        xp,yp,zp= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1] , nodes_nummpy_array[:,2]
        p=1
        nx1, ny1, nz1=_parse_grid(Filename_ph,p)
        xp=xp.reshape(nz1+1,ny1+1,nx1+1)
        zp=zp.reshape(nz1+1,ny1+1,nx1+1)
        xp =xp[:-1,1,:-1]
        zp =zp[:-1,1,:-1]
        Pic=Phase
        P_vtk = reader.GetOutput().GetCellData().GetArray(Pic.name)
        ph = vtk_to_numpy(P_vtk)
        ph=ph.reshape(nz1,ny1,nx1)
        if(Pic.PTT_show.yes == 'yes'):
            Pic.PTT_show.PTT_structure = PTT_Hist 
            
            
        if Pic.quiver.yes == 'yes' or Pic.stream_line.yes == 'yes':
                VTK_SET(Filename_dyn)
                VTK_UPDATE()    
            
            
                if(Pic.quiver.yes == 'yes'):
                    bufQ=reader.GetOutput().GetPointData().GetArray(quiver.field)
                else:
                    bufQ=reader.GetOutput().GetPointData().GetArray(stream_line.field)

                
                bufQ=vtk_to_numpy(bufQ)
                bufx,bufy,bufz= bufQ[:,0] , bufQ[:,1] , bufQ[:,2]
                bufx = bufx.reshape(nz,ny,nx)
                bufy = bufy.reshape(nz,ny,nx)
                bufz = bufz.reshape(nz,ny,nx)
                if Pic.quiver.yes == 'yes':
                    Pic.quiver.cx=bufx[:,0,:]
                    Pic.quiver.cy=bufy[:,0,:]
                    Pic.quiver.cz=bufz[:,0,:]
                else: 
                    Pic.stream_line.cx=bufx[:,0,:]
                    Pic.stream_line.cy=bufy[:,0,:]
                    Pic.stream_line.cz=bufz[:,0,:]

        if Pic.contour.yes == 'yes':
            Pic.contour.Scalar=buf2[:,1,:]
            
        
        _figure_maker_phase(ptsave,ph,xp,zp,t_cur,tstep,Pic,xd,zd)
        t2_start_A = perf_counter()  

        print("Phase properties took ", '{:0.2f}'.format(t2_start_A-t1_start_A),"seconds")

    
                    
        
                
                
                
    
    ipic+=1
    t1_end = perf_counter()
    tstep_time_pr = t1_end-t1_start
    minutes, seconds = divmod(t1_end-t1_start, 60)
        
       
    print("Timestep, ",tstep ,"took ","{:02}".format(int(minutes)),"minutes and","{:05.2f}".format(seconds),' seconds' )
    tot_time_perf +=tstep_time_pr

hours, rem = divmod(tot_time_perf, 3600)
minutes, seconds = divmod(rem, 60)
print("The visualization of the test took (h:m:s)","{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))

#for itime in Flist:
    

